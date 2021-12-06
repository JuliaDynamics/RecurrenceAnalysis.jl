#=
In this file the core computations for creating a recurrence matrix
are defined (via multiple dispatch).

The low level interface is contained in the function
`recurrence_matrix`, and this is where any specialization should happen.
=#
################################################################################
# AbstractRecurrenceMatrix type definitions and documentation strings
################################################################################
const FAN = NeighborNumber
export FAN

abstract type AbstractRecurrenceMatrix{T} end
const ARM = AbstractRecurrenceMatrix

struct RecurrenceMatrix{T} <: AbstractRecurrenceMatrix{T}
    data::SparseMatrixCSC{Bool,Int}
end
struct CrossRecurrenceMatrix{T} <: AbstractRecurrenceMatrix{T}
    data::SparseMatrixCSC{Bool,Int}
end
struct JointRecurrenceMatrix{T} <: AbstractRecurrenceMatrix{T}
    data::SparseMatrixCSC{Bool,Int}
end

function Base.summary(R::AbstractRecurrenceMatrix)
    N = nnz(R.data)
    return "$(nameof(typeof(R))) of size $(size(R.data)) with $N entries"
end
Base.show(io::IO, R::AbstractRecurrenceMatrix) = println(io, summary(R))

# Propagate used functions:
begin
    extentions = [
        (:Base, (:getindex, :size, :length, :view, :iterate)),
        (:LinearAlgebra, (:diag, :triu, :tril, :issymmetric)),
        (:SparseArrays, (:nnz, :rowvals, :nzrange, :nonzeros))
    ]
    for (M, fs) in extentions
        for f in fs
            @eval $M.$(f)(x::ARM, args...) = $(f)(x.data, args...)
        end
    end
end

for operator in [:(==), :(!=)]
    @eval Base.$operator(x::ARM, y::ARM) = $operator(x.data, y.data)
end

LinearAlgebra.issymmetric(::RecurrenceMatrix{WithinRange}) = true
LinearAlgebra.issymmetric(::JointRecurrenceMatrix{WithinRange}) = true
# column values in sparse matrix (parallel to rowvals)
function colvals(x::SparseMatrixCSC)
    cv = zeros(Int,nnz(x))
    @inbounds for c=1:size(x,2)
        cv[nzrange(x,c)] .= c
    end
    cv
end
colvals(x::ARM) = colvals(x.data)

# Convert to matrices
Base.Array{T}(R::ARM) where T = Matrix{T}(R.data)
Base.Matrix{T}(R::ARM) where T = Matrix{T}(R.data)
Base.Array(R::ARM) = Matrix(R.data)
Base.Matrix(R::ARM) = Matrix(R.data)
SparseArrays.SparseMatrixCSC{T}(R::ARM) where T = SparseMatrixCSC{T}(R.data)
SparseArrays.SparseMatrixCSC(R::ARM) = SparseMatrixCSC(R.data)

"""
    RecurrenceMatrix(x, ε; kwargs...)
    RecurrenceMatrix{FAN}(x, ε; kwargs...)

Create a recurrence matrix from trajectory `x` (either a `Dataset` or a `Vector`).
Objects of type `<:AbstractRecurrenceMatrix` are displayed as a [`recurrenceplot`](@ref).

## Description

The recurrence matrix is a numeric representation of a "recurrence plot" [1, 2],
in the form of a sparse square matrix of Boolean values.

`x` must be a `Vector` or an `AbstractDataset`
(possibly representing an embedded phase space; see [`embed`](@ref)).
If `d(x[i], x[j]) ≤ ε` (with `d` the distance function),
then the cell `(i, j)` of the matrix will have a `true`
value. The criteria to evaluate distances between data points are defined
by the following keyword arguments:

* `scale=1` : a function of the distance matrix (see [`distancematrix`](@ref)),
  or a fixed number, used to scale the value of `ε`. Typical choices are
  `maximum` or `mean`, such that the threshold `ε` is defined as a ratio of the
  maximum or the mean distance between data points, respectively (using
  `mean` or `maximum` calls specialized versions that are faster than the naive
  approach).  Use `1` to keep the distances unscaled (default).
* `fixedrate::Bool=false` : a flag that indicates if `ε` should be
  taken as a target fixed recurrence rate (see [`recurrencerate`](@ref)).
  If `fixedrate` is set to `true`, `ε` must be a value between 0 and 1,
  and `scale` is ignored.
* `metric="euclidean"` : metric of the distances, either `Metric` or a string,
   as in [`distancematrix`](@ref).
* `parallel::Bool=false` : whether to parallelize the computation of the recurrence
   matrix.  This will split the computation of the matrix across multiple threads.

The parametrized constructor `RecurrenceMatrix{NeighborNumber}` creates the recurrence matrix
with a fixed number of neighbors for each point in the phase space, i.e. the number
of recurrences is the same for all columns of the recurrence matrix.
In such case, `ε` is taken as the target fixed local recurrence rate,
defined as a value between 0 and 1, and `scale` and `fixedrate` are ignored.
This is often referred to in the literature as the method of "Fixed Amount of Nearest Neighbors"
(or FAN for short); `RecurrenceMatrix{FAN}` can be used as a convenient alias
for `RecurrenceMatrix{NeighborNumber}`.

If no parameter is specified, `RecurrenceMatrix` returns a
`RecurrenceMatrix{WithinRange}` object, meaning that recurrences will be taken
for pairs of data points whose distance is within the range determined by
the input arguments. Note that while recurrence matrices
with neighbors defined within a given range are always symmetric, those defined
by a fixed amount of neighbors can be non-symmetric.

See also: [`CrossRecurrenceMatrix`](@ref), [`JointRecurrenceMatrix`](@ref) and
use [`recurrenceplot`](@ref) to turn the result of these functions into a plottable format.

## References
[1] : N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems",
*Phys. Reports 438*(5-6), 237-329 (2007).
[DOI:10.1016/j.physrep.2006.11.001](https://doi.org/10.1016/j.physrep.2006.11.001)

[2] : N. Marwan & C.L. Webber, "Mathematical and computational foundations of
recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence
Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).
"""
function RecurrenceMatrix{WithinRange}(x, ε; metric = DEFAULT_METRIC,
    parallel::Bool = length(x) > 500 && Threads.nthreads() > 1,
    kwargs...)
    m = getmetric(metric)
    s = resolve_scale(x, m, ε; kwargs...)
    m = recurrence_matrix(x, m, s, Val(parallel))
    return RecurrenceMatrix{WithinRange}(m)
end

RecurrenceMatrix(args...; kwargs...) = RecurrenceMatrix{WithinRange}(args...; kwargs...)

function RecurrenceMatrix{NeighborNumber}(x, ε; metric = DEFAULT_METRIC,
    parallel::Bool = length(x) > 500 && Threads.nthreads() > 1,
    kwargs...)
    m = getmetric(metric)
    s = get_fan_threshold(x, m, ε)
    m = recurrence_matrix(x, x, m, s, Val(parallel)) # to allow for non-symmetry
    return RecurrenceMatrix{NeighborNumber}(m)
end

"""
    CrossRecurrenceMatrix(x, y, ε; kwargs...)
    CrossRecurrenceMatrix{FAN}(x, y, ε; kwargs...)

Create a cross recurrence matrix from trajectories `x` and `y`.

The cross recurrence matrix is a bivariate extension of the recurrence matrix.
For the time series `x`, `y`, of length `n` and `m`, respectively, it is a
sparse `n×m` matrix of Boolean values, such that if `d(x[i], y[j]) ≤ ε`,
then the cell `(i, j)` of the matrix will have a `true` value.

Note that, unlike univariate recurrence matrices, cross recurrence matrices
are not generally symmetric, regardless of the method used to make them.

See [`RecurrenceMatrix`](@ref) for details, references and keywords.
See also: [`JointRecurrenceMatrix`](@ref).
"""
function CrossRecurrenceMatrix{WithinRange}(x, y, ε; metric = DEFAULT_METRIC,
    parallel::Bool = length(x) > 500 && Threads.nthreads() > 1,
    kwargs...)
    m = getmetric(metric)
    s = resolve_scale(x, y, m, ε; kwargs...)
    m = recurrence_matrix(x, y, m, s, Val(parallel))
    return CrossRecurrenceMatrix{WithinRange}(m)
end

CrossRecurrenceMatrix(args...; kwargs...) = CrossRecurrenceMatrix{WithinRange}(args...; kwargs...)

function CrossRecurrenceMatrix{NeighborNumber}(x, y, ε; metric = DEFAULT_METRIC,
    parallel::Bool = length(x) > 500 && Threads.nthreads() > 1,
    kwargs...)
    m = getmetric(metric)
    s = get_fan_threshold(x, y, m, ε)
    m = recurrence_matrix(x, y, m, s, Val(parallel))
    return CrossRecurrenceMatrix{NeighborNumber}(m)
end

"""
    JointRecurrenceMatrix(x, y, ε; kwargs...)
    JointRecurrenceMatrix{FAN}(x, y, ε; kwargs...)

Create a joint recurrence matrix from `x` and `y`.

The joint recurrence matrix considers the recurrences of the trajectories
of `x` and `y` separately, and looks for points where both recur
simultaneously. It is calculated by the element-wise multiplication
of the recurrence matrices of `x` and `y`. If `x` and `y` are of different
length, the recurrences are only calculated until the length of the shortest one.

See [`RecurrenceMatrix`](@ref) for details, references and keywords.
See also: [`CrossRecurrenceMatrix`](@ref).
"""
function JointRecurrenceMatrix{T}(x, y, ε; kwargs...) where T
    n = min(size(x,1), size(y,1))
    if n == size(x,1) && n == size(y,1)
        rm1 = RecurrenceMatrix{T}(x, ε; kwargs...)
        rm2 = RecurrenceMatrix{T}(y, ε; kwargs...)
    else
        rm1 = RecurrenceMatrix{T}(x[1:n,:], ε; kwargs...)
        rm2 = RecurrenceMatrix{T}(y[1:n,:], ε; kwargs...)
    end
    return JointRecurrenceMatrix{T}(rm1.data .* rm2.data)
end

JointRecurrenceMatrix(args...; kwargs...) = JointRecurrenceMatrix{WithinRange}(args...; kwargs...)

"""
    JointRecurrenceMatrix(R1, R2; kwargs...)

Create a joint recurrence matrix from given recurrence matrices `R1, R2`.
"""
function JointRecurrenceMatrix(
    R1::AbstractRecurrenceMatrix{T}, R2::AbstractRecurrenceMatrix{T}; kwargs...
    ) where T
    R3 = R1.data .* R2.data
    return JointRecurrenceMatrix{T}(R3)
end

################################################################################
# Scaling / fixed rate / fixed amount of neighbors
################################################################################
# here args... is (x, y, metric, ε) or just (x, metric, ε)
function resolve_scale(args...; scale=1, fixedrate=false)
    ε = args[end]
    # Check fixed recurrence rate - ε must be within (0, 1)
    if fixedrate
        sfun = (m) -> quantile(vec(m), ε)
        return resolve_scale(Base.front(args)..., 1.0; scale=sfun, fixedrate=false)
    else
        scale_value = _computescale(scale, Base.front(args)...)
        return ε*scale_value
    end
end

# If `scale` is a function, compute the numeric value of the scale based on the
# distance matrix; otherwise return the value of `scale` itself
_computescale(scale::Real, args...) = scale
_computescale(scale::Function, x, metric::Metric) =
_computescale(scale, x, x, metric)

# generic method that uses `distancematrix`
function _computescale(scale::Function, x, y, metric)
    if x===y
        distances = zeros(Int(length(x)*(length(x)-1)/2), 1)
        c = 0
        @inbounds for i in 1:length(x)-1, j=(i+1):length(x)
            distances[c+=1] = evaluate(metric, x[i], y[j])
        end
    else
        distances = distancematrix(x, y, metric)
    end
    return scale(distances)
end
# specific methods to avoid `distancematrix`
function _computescale(scale::typeof(maximum), x, y, metric::Metric)
    maxvalue = zero(eltype(x))
    if x===y
        @inbounds for i in 1:length(x)-1, j=(i+1):length(x)
            newvalue = evaluate(metric, x[i], y[j])
            (newvalue > maxvalue) && (maxvalue = newvalue)
        end
    else
        @inbounds for xi in x, yj in y
            newvalue = evaluate(metric, xi, yj)
            (newvalue > maxvalue) && (maxvalue = newvalue)
        end
    end
    return maxvalue
end
function _computescale(scale::typeof(mean), x, y, metric::Metric)
    meanvalue = 0.0
    if x===y
        @inbounds for i in 1:length(x)-1, j=(i+1):length(x)
            meanvalue += evaluate(metric, x[i], y[j])
        end
        denominator = length(x)*(length(x)-1)/2
    else
        @inbounds for xi in x, yj in y
            meanvalue += evaluate(metric, xi, yj)
        end
        denominator = Float64(length(x)*length(y))
    end
    return meanvalue/denominator
end

"""
    get_fan_threshold(args) → ε_fan
Compute the adaptive Fixed Amount of Neibours (FAN) `ε_fan`
Here args... is (x, y, metric, ε) or just (x, metric, ε)
"""
function get_fan_threshold(x, y, metric, ε)
    @assert 0 < ε < 1 "Global recurrence rate must be ∈ (0, 1)"
    fan_threshold = zeros(length(y))
    d = distancematrix(x, y, metric)
    if x === y
        ε += 1/length(x)
    end
    for i in 1:size(d, 2)
        fan_threshold[i] = quantile(view(d, : ,i), ε)
    end
    return fan_threshold
end

get_fan_threshold(x, metric, ε) = get_fan_threshold(x, x, metric, ε)


################################################################################
# recurrence_matrix - Low level interface
################################################################################
# TODO: increase the efficiency here by not computing everything!

# Core function

# First, we define the non-parallel versions.

# For two datasets

"""
    recurrence_matrix(x, y, metric::Metric, ε, parallel::Val)

Return a sparse matrix which encodes recurrence points.

Note that `parallel` may be either `Val(true)` or `Val(false)`.
"""
function recurrence_matrix(xx::AbstractDataset, yy::AbstractDataset, metric::Metric, ε, ::Val{false})
    x = xx.data
    y = yy.data
    @assert ε isa Real || length(ε) == length(y)
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
    for j in 1:length(y)
        nzcol = 0
        for i in 1:length(x)
            @inbounds if evaluate(metric, x[i], y[j]) ≤ ( (ε isa Real) ? ε : ε[j] )
                push!(rowvals, i)
                nzcol += 1
            end
        end
        append!(colvals, fill(j, (nzcol,)))
    end
    nzvals = fill(true, (length(rowvals),))
    return sparse(rowvals, colvals, nzvals, length(x), length(y))
end

# Vector version can be more specialized (and metric is irrelevant)
function recurrence_matrix(x::AbstractVector, y::AbstractVector, metric::Metric, ε, ::Val{false})
    @assert ε isa Real || length(ε) == length(y)
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
    for j in 1:length(y)
        nzcol = 0
        for i in 1:length(x)
            if @inbounds abs(x[i] - y[j]) ≤ ( (ε isa Real) ? ε : ε[j] )
                push!(rowvals, i)
                nzcol += 1
            end
        end
        append!(colvals, fill(j, (nzcol,)))
    end
    nzvals = fill(true, (length(rowvals),))
    return sparse(rowvals, colvals, nzvals, length(x), length(y))
end

# For one dataset
function recurrence_matrix(x::AbstractVector, metric::Metric, ε, ::Val{false})
    @assert ε isa Real || length(ε) == length(y)
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
    for j in 1:length(x)
        nzcol = 0
        for i in 1:j
            if @inbounds abs(x[i] - x[j]) ≤ ( (ε isa Real) ? ε : ε[j] )
                push!(rowvals, i)
                nzcol += 1
            end
        end
        append!(colvals, fill(j, (nzcol,)))
    end
    nzvals = fill(true, (length(rowvals),))
    return Symmetric(sparse(rowvals, colvals, nzvals, length(x), length(x)), :U)
end

function recurrence_matrix(xx::AbstractDataset, metric::Metric, ε, ::Val{false})
    x = xx.data
    @assert ε isa Real || length(ε) == length(y)
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
    for j in 1:length(x)
        nzcol = 0
        for i in 1:j
            @inbounds if evaluate(metric, x[i], x[j]) ≤ ( (ε isa Real) ? ε : ε[j] )
                push!(rowvals, i)
                nzcol += 1
            end
        end
        append!(colvals, fill(j, (nzcol,)))
    end
    nzvals = fill(true, (length(rowvals),))
    return Symmetric(sparse(rowvals, colvals, nzvals, length(x), length(x)), :U)
end


# Now, we define the parallel versions of these functions.

# Core function
function recurrence_matrix(xx::AbstractDataset, yy::AbstractDataset, metric::Metric, ε, ::Val{true})
    x = xx.data
    y = yy.data
    @assert ε isa Real || length(ε) == length(y)
    # We create an `Array` of `Array`s, for each thread to have its
    # own array to push to.  This avoids race conditions with
    # multiple threads pushing to the same `Array` (`Array`s are not atomic).
    rowvals = [Vector{Int}() for _ in 1:Threads.nthreads()]
    colvals = [Vector{Int}() for _ in 1:Threads.nthreads()]

    # This is the same logic as the serial function, but parallelized.
    Threads.@threads for j in 1:length(y)
        threadn = Threads.threadid()
        nzcol = 0
        for i in 1:length(x)
            @inbounds if evaluate(metric, x[i], y[j]) ≤ ( (ε isa Real) ? ε : ε[j] )
                push!(rowvals[threadn], i) # push to the thread-specific row array
                nzcol += 1
            end
        end
        append!(colvals[threadn], fill(j, (nzcol,)))
    end
    finalrows = vcat(rowvals...) # merge into one array
    finalcols = vcat(colvals...) # merge into one array
    nzvals = fill(true, (length(finalrows),))
    return sparse(finalrows, finalcols, nzvals, length(x), length(y))
end

# Vector version can be more specialized (and metric is irrelevant)
function recurrence_matrix(x::AbstractVector, y::AbstractVector, metric::Metric, ε, ::Val{true})
    ε isa Real || length(ε) == length(x)
    # We create an `Array` of `Array`s, for each thread to have its
    # own array to push to.  This avoids race conditions with
    # multiple threads pushing to the same `Array` (`Array`s are not atomic).
    rowvals = [Vector{Int}() for _ in 1:Threads.nthreads()]
    colvals = [Vector{Int}() for _ in 1:Threads.nthreads()]

    # This is the same logic as the serial function, but parallelized.
    Threads.@threads for j in 1:length(y)
        threadn = Threads.threadid()
        nzcol = 0
        for i in 1:length(x)
            @inbounds if abs(x[i] - y[j]) ≤ ( (ε isa Real) ? ε : ε[j] )
                push!(rowvals[threadn], i) # push to the thread-specific row array
                nzcol += 1
            end
        end
        append!(colvals[threadn], fill(j, (nzcol,)))
    end
    finalrows = vcat(rowvals...) # merge into one array
    finalcols = vcat(colvals...) # merge into one array
    nzvals = fill(true, (length(finalrows),))
    return sparse(finalrows, finalcols, nzvals, length(x), length(y))
end


function recurrence_matrix(x::AbstractVector, metric::Metric, ε, ::Val{true})
    @assert ε isa Real || length(ε) == length(y)
    # We create an `Array` of `Array`s, for each thread to have its
    # own array to push to.  This avoids race conditions with
    # multiple threads pushing to the same `Array` (`Array`s are not atomic).
    rowvals = [Vector{Int}() for _ in 1:Threads.nthreads()]
    colvals = [Vector{Int}() for _ in 1:Threads.nthreads()]

    # This is the same logic as the serial function, but parallelized.
    Threads.@threads for k in partition_indices(length(x))
        threadn = Threads.threadid()
        for j in k
            nzcol = 0
            for i in 1:j
                @inbounds if abs(x[i] - x[j]) ≤ ( (ε isa Real) ? ε : ε[j] )
                    push!(rowvals[threadn], i) # push to the thread-specific row array
                    nzcol += 1
                end
            end
            append!(colvals[threadn], fill(j, (nzcol,)))
        end
    end
    finalrows = vcat(rowvals...) # merge into one array
    finalcols = vcat(colvals...) # merge into one array
    nzvals = fill(true, (length(finalrows),))
    return Symmetric(sparse(finalrows, finalcols, nzvals, length(x), length(x)), :U)
end

function recurrence_matrix(xx::AbstractDataset, metric::Metric, ε, ::Val{true})
    x = xx.data
    @assert ε isa Real || length(ε) == length(y)
    # We create an `Array` of `Array`s, for each thread to have its
    # own array to push to.  This avoids race conditions with
    # multiple threads pushing to the same `Array` (`Array`s are not atomic).
    rowvals = [Vector{Int}() for _ in 1:Threads.nthreads()]
    colvals = [Vector{Int}() for _ in 1:Threads.nthreads()]

    # This is the same logic as the serial function, but parallelized.
    Threads.@threads for k in partition_indices(length(x))
        threadn = Threads.threadid()
        for j in k
            nzcol = 0
            for i in 1:j
                @inbounds if evaluate(metric, x[i], x[j]) ≤ ( (ε isa Real) ? ε : ε[j] )
                    push!(rowvals[threadn], i) # push to the thread-specific row array
                    nzcol += 1
                end
            end
            append!(colvals[threadn], fill(j, (nzcol,)))
        end
    end
    finalrows = vcat(rowvals...) # merge into one array
    finalcols = vcat(colvals...) # merge into one array
    nzvals = fill(true, (length(finalrows),))
    return Symmetric(sparse(finalrows, finalcols, nzvals, length(x), length(x)), :U)
end

# Transforms the standard RP into a close returns map
function convert_recurrence_matrix(R::SparseMatrixCSC)

    N = size(R)
    # init new matrix
    Y = spzeros(Bool, 2*N[1]+1,N[1])
    # fill rows of Y with the diagonals of R
    # upper triangle
    for i = 0:N[1]-1
       Y[N[1]+i+1,(1:(N[1]-i))] = diag(R,i)
    end
    # lower triangle
    for i = 0:N[1]-1
       Y[N[1]-i+1,(1:(N[1]-i)).+i] = diag(R,-i)
    end
    return Y
end

# Transforms the reverted RP (close returns map) into a normal RP
function revert_close_returns_map(R::SparseMatrixCSC)

    N = size(R)
    # init new matrix
    Y = spzeros(Bool, N[2],N[2])
    # make R to a square matrix, fill the new part with zeros
    Z = [R spzeros(N[1],N[1]+1)]
    Z = Z[end:-1:1,:]

    # fill columns of Y with the diagonals of Z (but only the first N points)
    for i = 1:N[2]
        di = diag(Z,-i)
        Y[:,N[2]-i+1] = di[1:N[2]]
    end
    return Y
end

# from the 3-by-N Matrix, which stores all horizontal lines in the reverted RP,
# (close returns map) [obtained by calling first `revert_recurrence_matrix` on a
# recurrence Matrix R, and consecutevely `horizontalhisto(R)`] construct and return
# a "dummy" close returns map, where all lines are encoded according to their length.
function create_close_returns_map(lines::AbstractMatrix{Int}, size_of_cl_ret_RP::Tuple{Int, Int})
    X = zeros(Int, size_of_cl_ret_RP)
    for i = 1:size(lines,2)
        line_ind = lines[2,i]
        column_ind = lines[3,i]
        for j = 0:lines[1,i]-1
            X[line_ind,column_ind+j] = lines[1,i]
        end
    end
    return sparse(X)
end

# deletes a line, specified in 'l_vec' (line vector, with first line being
# the total line length, the second line the line-index of the starting point
# and the third line the column-index of the starting pint) from the close
# returns "RP" 'R'.
function delete_line_from_cl_ret_RP!(R::SparseMatrixCSC,line::Vector{Int})
    @assert length(line) == 3
    R[line[2], line[3] .+ (1:line[1]).-1] .= 0
end


"""
    skeletonize(R) → R_skel

Skeletonizes the recurrence matrix `R` (see [`recurrence_matrix`](@ref)) by using
the algorithm proposed by Kraemer & Marwan [^Kraemer2019]. This function returns
`R_skel`, a recurrence matrix, which only consists of diagonal lines of
"thickness" one.

## References
[^Kraemer2019]: Kraemer, K.H., Marwan, N. (2019). [Border effect corrections for diagonal line based recurrence quantification analysis measures. Physics Letters A 383(34)](https://doi.org/10.1016/j.physleta.2019.125977).
"""
function skeletonize(X::Union{ARM,SparseMatrixCSC})
    if issymmetric(X)
        symm = true
        # convert lower triangle into a close returns map
        X_cl1 = convert_recurrence_matrix(tril(X))
        # get "horizontal" line distribution with position indices
        lines1 = horizontalhisto(X_cl1)
        lines_copy1 = deepcopy(lines1)

        Nlines1 = size(lines1,2) # number of found lines

        # create a close returns map with horizontal lines represented by
        # numbers, equal to their lengths
        X_hori1 = create_close_returns_map(lines1, size(X_cl1))
    else
        symm = false
        # convert upper and lower triangles into close returns maps
        X_cl1 = convert_recurrence_matrix(tril(X))
        X_cl2 = convert_recurrence_matrix(sparse(triu(X)'))

        # get "horizontal" line distribution with position indices
        lines1 = horizontalhisto(X_cl1)
        lines2 = horizontalhisto(X_cl2)
        lines_copy1 = deepcopy(lines1)
        lines_copy2 = deepcopy(lines2)

        Nlines1 = size(lines1,2) # number of found lines in lower triangle
        Nlines2 = size(lines2,2) # number of found lines in upper triangle

        # create close returns maps with horizontal lines represented by
        # numbers, equal to their lengths
        X_hori1 = create_close_returns_map(lines1, size(X_cl1))
        X_hori2 = create_close_returns_map(lines2, size(X_cl2))
    end

    # scan the lines, start with the longest one and discard all adjacent lines

    # initialize final line matrix
    line_matrix_final = zeros(Int, 3, 1)
    # go through all lines stored in the sorted line matrix
    N, M = size(X_hori1)
    for l_ind = 1:Nlines1
        # check if line is still in the rendered line matrix
        ismember = [lines1[:,l_ind]==lines_copy1[:,i] for i = 1:size(lines_copy1,2)]
        ~any(ismember) ? continue : nothing
        # get index pair for start of the line
        linei, columni = lines1[2,l_ind], lines1[3,l_ind]
        # copy this line in the final line matrix
        line_matrix_final = hcat(line_matrix_final, lines1[:,l_ind])
        # delete this line from the RP
        delete_line_from_cl_ret_RP!(X_hori1, lines1[:,l_ind])
        # go along each point of the line and check for neighbours
        l_max = lines1[1,l_ind]
        for l = 1:l_max
            # scan each line twice - above and underneth
            for index = -1:2:1
                # make sure not to exceed RP-boundaries
                (linei+index > N || linei+index == 0) ? break : nothing
                # if there is a neighbouring point, call recursive scan-function
                (X_hori1[linei+index, columni+l-1] != 0) ? lines_copy1 = scan_lines!(X_hori1, lines_copy1, linei+index, columni+l-1) : nothing
            end
        end
    end

    # if not symmetric input RP, than compute for the upper triangle as well
    if ~symm
        # initialize final line matrix
        line_matrix_final2 = zeros(Int, 3, 1)
        for l_ind = 1:Nlines2
            # check if line is still in the rendered line matrix
            ismember = [lines2[:,l_ind]==lines_copy2[:,i] for i = 1:size(lines_copy2,2)]
            ~any(ismember) ? continue : nothing
            # get index pair for start of the line
            linei, columni = lines2[2,l_ind], lines2[3,l_ind]
            # copy this line in the final line matrix
            line_matrix_final2 = hcat(line_matrix_final2, lines2[:,l_ind])
            # delete this line from the RP
            delete_line_from_cl_ret_RP!(X_hori2, lines2[:,l_ind])
            # go along each point of the line and check for neighbours
            l_max = lines2[1,l_ind]
            for l = 1:l_max
                # scan each line twice - above and underneth
                for scan = 1:2
                    if scan == 1
                        index = 1
                        (linei+index > N) ? break : nothing # make sure not to exceed RP-boundaries
                    else
                        index = -1
                        (linei+index == 0) ? break : nothing # make sure not to exceed RP-boundaries
                    end
                    # if there is a neighbouring point, call recursive scan-function
                    (X_hori2[linei+index,columni+l-1] != 0) ? lines_copy2 = scan_lines!(X_hori2, lines_copy2, linei+index, columni+l-1) : nothing
                end
            end
        end
    end

    # build RP based on the histogramm of the reduced lines
    X_cl_new = spzeros(Bool, N, M)
    ~symm ? X_cl2_new = spzeros(Bool, N, M) : nothing
    # fill up close returns map with lines stored in the new line matrix
    for i = 1:size(line_matrix_final, 2)
        l_max = line_matrix_final[1,i]
        linei = line_matrix_final[2,i]
        columni = line_matrix_final[3,i]
        for j = 1:l_max
            X_cl_new[linei, columni+j-1] = true
        end
    end
    if symm
        XX = revert_close_returns_map(X_cl_new) # revert this close returns map into a legal RP
        X_new = XX .+ XX'
        X_new[diagind(X_new)] .= true # LOI
    else
        # fill up close returns map with lines stored in the new line matrix
        for i = 1:size(line_matrix_final2,2)
            l_max = line_matrix_final2[1,i]
            linei = line_matrix_final2[2,i]
            columni = line_matrix_final2[3,i]
            for j = 1:l_max
                X_cl2_new[linei,columni+j-1] = true
            end
        end
        XX = revert_close_returns_map(X_cl_new) # revert this close returns map into a legal RP
        XXX= revert_close_returns_map(X_cl2_new)
        X_new = XX .+ XXX'
        X_new[diagind(X_new)] .= true # LOI
    end
    return X_new
end
