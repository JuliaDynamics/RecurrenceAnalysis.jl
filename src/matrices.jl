#=
In this file the core computations for creating a recurrence matrix
are defined (via multiple dispatch).

The low level interface is contained in the function
`recurrence_matrix`, and this is where any specialization should happen.
=#

################################################################################
# Distance matrix
################################################################################
"""
    distancematrix(x [, y = x], metric = "euclidean")

Create a matrix with the distances between each pair of points of the
time series `x` and `y` using `metric`.

The time series `x` and `y` can be `Dataset`s or vectors or matrices with data points
in rows.
The data point dimensions (or number of columns) must be the same for `x` and `y`.
The returned value is a `n×m` matrix, with `n` being the length (or number of rows)
of `x`, and `m` the length of `y`.

The metric can be identified by a string, or any of the `Metric`s defined in
the [`Distances` package](https://github.com/JuliaStats/Distances.jl).
The list of strings available to define the metric are:

* `"max"` or `"inf"` for the maximum or L∞ norm
  (`Chebyshev()` in the `Distances` package).
* `"euclidean"` for the L2 or Euclidean norm, used by default
  (`Euclidean()` in `Distances`).
* `"manhattan"`, `"cityblock"`, `"taxicab"` or `"min"` for the Manhattan or L1 norm
  (`Cityblock()` in `Distances`).
"""
distancematrix(x, metric::Union{Metric,String}=DEFAULT_METRIC, parallel = size(x)[1] > 500) = _distancematrix(x, getmetric(metric), Val(parallel))

# For 1-dimensional arrays (vectors), the distance does not depend on the metric
distancematrix(x::Vector, y::Vector, metric=DEFAULT_METRIC, parallel = size(x)[1] > 500) = abs.(x .- y')

# If the metric is supplied as a string, get the corresponding Metric from Distances
distancematrix(x, y, metric::String, parallel = size(x)[1] > 500) = distancematrix(x, y, getmetric(metric), parallel)

const MAXDIM = 9
function distancematrix(x::Tx, y::Ty, metric::Metric=DEFAULT_METRIC, parallel = size(x)[1] > 500) where
         {Tx<:Union{AbstractMatrix, Dataset}} where {Ty<:Union{AbstractMatrix, Dataset}}
    sx, sy = size(x), size(y)
    @assert sx[2] == sy[2] """
        The dimensions of the data points in `x` and `y` must be equal!
        Found dim(x)=$(sx[2]), dim(y)=$(sy[2]).
        """

    if sx[2] ≥ MAXDIM && typeof(metric) == Euclidean # Blas optimization
        return _distancematrix(Matrix(x), Matrix(y), metric, Val(parallel))
    elseif Tx <: Matrix && Ty <: Matrix && metric == Chebyshev()
        return _distancematrix(x, y, metric, Val(parallel))
    else
        return _distancematrix(Dataset(x), Dataset(y), metric, Val(parallel))
    end
end

# Core function for Matrices (wrapper of `pairwise` from the Distances package)
_distancematrix(x::AbstractMatrix, y::AbstractMatrix, metric::Metric, ::Val{false}) = pairwise(metric, x', y', dims=2)

# First we define the serial versions of the functions.
# Core function for Datasets
function _distancematrix(x::Dataset{S,Tx}, y::Dataset{S,Ty},
    metric::Metric, ::Val{false}) where {S, Tx, Ty}

    x = x.data
    y = y.data
    d = zeros(promote_type(Tx,Ty), length(x), length(y))
    for j in 1:length(y)
        for i in 1:length(x)
            @inbounds d[i,j] = evaluate(metric, x[i], y[j])
        end
    end
    return d
end

# Now, we define the parallel versions.

function _distancematrix(x::Dataset{S,Tx}, y::Dataset{S,Ty},
    metric::Metric, ::Val{true}) where {S, Tx, Ty}

    x = x.data
    y = y.data
    d = zeros(promote_type(Tx,Ty), length(x), length(y))
    Threads.@threads for j in 1:length(y)
        for i in 1:length(x)
            @inbounds d[i,j] = evaluate(metric, x[i], y[j])
        end
    end
    return d
end

function _distancematrix(x::Matrix{Tx}, y::Matrix{Ty},
    metric::Metric, ::Val{true}) where {Tx, Ty}

    x = x.data
    y = y.data
    d = zeros(promote_type(Tx,Ty), length(x), length(y))
    Threads.@threads for j in 1:length(y)
        for i in 1:length(x)
            @inbounds d[i,j] = evaluate(metric, x[i, :], y[j, :])
        end
    end
    return d
end


# Now, we define methods for the single trajectory.
# We can speed this up by only calculating the lower triangle,
# which halves the number of computations needed.

# Again, we'll define the serial version first:
function _distancematrix(x::Vector{T}, metric::Metric, ::Val{false}) where T
    d = zeros(T, length(x), length(x))

    for j in 1:length(x)
        for i in 1:j
            @inbounds d[i, j] = abs(x[i] - x[j])
        end
    end

    return Symmetric(d, :U)

end

function _distancematrix(x::Dataset{S, T}, metric::Metric, ::Val{false}) where T where S
    d = zeros(T, length(x), length(x))

    for j in 1:length(x)
        for i in 1:j
            @inbounds d[i, j] = evaluate(metric, x[i], x[j])
        end
    end

    return Symmetric(d, :U)

end

# Now, we define the parallel version.  There's a twist, though.

# The single trajectory case is a little tricky.  If you try it naively,
# using the method we used for the serial computation above, the points are
# partioned out unequally between threads.  This means that performance actually
# **decreases** compared to the full parallel version.  To mitigate this, we need
# to partition the computation equally among all threads.
# The function I've written below does essentially this - given a length,
# it calculates the appropriate partitioning, and returns a list of indices,
# by utilizing the fact that the area is proportional to the square of the height.
# It partitions the "triangle" which needs to be computed into "trapezoids",
# which all have an equal area.
function partition_indices(len)
    indices = Vector{UnitRange{Int}}(undef, Threads.nthreads())
    length = len
    offset = 1

    # Here, we have to iterate in reverse, to get an equal area every time
    for n in Threads.nthreads():-1:1
        partition = round(Int, length / sqrt(n)) # length varies as square root of area
        indices[n] = offset:(partition .+ offset .- 1)
        length -= partition
        offset += partition
    end

    return indices
end

# And now for the actual implementation:
function _distancematrix(x::Vector{T}, metric::Metric, ::Val{true}) where T
    d = zeros(T, length(x), length(x))

    Threads.@threads for k in partition_indices(length(x))
        for j in k
            for i in 1:j
                @inbounds d[i, j] = abs(x[i] - x[j])
            end
        end
    end

    return Symmetric(d, :U)

end

function _distancematrix(x::Dataset{S, T}, metric::Metric, ::Val{true}) where T where S
    d = zeros(T, length(x), length(x))

    Threads.@threads for k in partition_indices(length(x))
        for j in k
            for i in 1:j
                @inbounds d[i, j] = evaluate(metric, x[i], x[j])
            end
        end
    end

    return Symmetric(d, :U)

end

################################################################################
# AbstractRecurrenceMatrix type definitions and documentation strings
################################################################################
abstract type AbstractRecurrenceMatrix end
const ARM = AbstractRecurrenceMatrix

struct RecurrenceMatrix <: AbstractRecurrenceMatrix
    data::SparseMatrixCSC{Bool,Int}
end
struct CrossRecurrenceMatrix <: AbstractRecurrenceMatrix
    data::SparseMatrixCSC{Bool,Int}
end
struct JointRecurrenceMatrix <: AbstractRecurrenceMatrix
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

LinearAlgebra.issymmetric(::RecurrenceMatrix) = true
# column values in sparse matrix (parallel to rowvals)
function colvals(x::SparseMatrixCSC)
    cv = zeros(Int,nnz(x))
    @inbounds for c=1:size(x,2)
        cv[nzrange(x,c)] .= c
    end
    cv
end
colvals(x::ARM) = colvals(x.data)


"""
    RecurrenceMatrix(x, ε; kwargs...)

Create a recurrence matrix from trajectory `x`.
Objects of type `<:AbstractRecurrenceMatrix` are displayed as a [`recurrenceplot`](@ref).

## Description

The recurrence matrix is a numeric representation of a "recurrence plot" [^1, ^2],
in the form of a sparse square matrix of Boolean values.

`x` must be a `Vector` or a `Dataset` or a `Matrix` with data points in rows
(possibly representing and embedded phase space; see [`embed`](@ref)).
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
* `parallel=false` : whether to parallelize the computation of the recurrence
   matrix.  This will split the computation of the matrix across multiple threads.
   This **may not work on Julia versions before v1.3**, so be warned!

See also: [`CrossRecurrenceMatrix`](@ref), [`JointRecurrenceMatrix`](@ref) and
use [`recurrenceplot`](@ref) to turn the result of these functions into a plottable format.

## References
[^1] : N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems",
*Phys. Reports 438*(5-6), 237-329 (2007).

[^2] : N. Marwan & C.L. Webber, "Mathematical and computational foundations of
recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence
Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).
"""
function RecurrenceMatrix(x, ε; metric = DEFAULT_METRIC,
    parallel::Bool = length(x) > 500 && Threads.nthreads() > 1, kwargs...)
    m = getmetric(metric)
    s = resolve_scale(x, m, ε; kwargs...)
    m = recurrence_matrix(x, m, s, Val(parallel))
    return RecurrenceMatrix(m)
end

"""
    CrossRecurrenceMatrix(x, y, ε; kwargs...)

Create a cross recurrence matrix from trajectories `x` and `y`.

The cross recurrence matrix is a bivariate extension of the recurrence matrix.
For the time series `x`, `y`, of length `n` and `m`, respectively, it is a
sparse `n×m` matrix of Boolean values, such that if `d(x[i], y[j]) ≤ ε`,
then the cell `(i, j)` of the matrix will have a `true` value.

See [`RecurrenceMatrix`](@ref) for details, references and keywords.
See also: [`JointRecurrenceMatrix`](@ref).
"""
function CrossRecurrenceMatrix(x, y, ε; metric = DEFAULT_METRIC,
    parallel::Bool = length(x) > 500 && Threads.nthreads() > 1, kwargs...)
    m = getmetric(metric)
    s = resolve_scale(x, y, m, ε; kwargs...)
    m = recurrence_matrix(x, y, m, s, Val(parallel))
    return CrossRecurrenceMatrix(m)
end


"""
    JointRecurrenceMatrix(x, y, ε; kwargs...)

Create a joint recurrence matrix from `x` and `y`.

The joint recurrence matrix considers the recurrences of the trajectories
of `x` and `y` separately, and looks for points where both recur
simultaneously. It is calculated by the element-wise multiplication
of the recurrence matrices of `x` and `y`. If `x` and `y` are of different
length, the recurrences are only calculated until the length of the shortest one.

See [`RecurrenceMatrix`](@ref) for details, references and keywords.
See also: [`CrossRecurrenceMatrix`](@ref).
"""
function JointRecurrenceMatrix(x, y, ε; kwargs...)
    n = min(size(x,1), size(y,1))
    if n == size(x,1) && n == size(y,1)
        rm1 = RecurrenceMatrix(x, ε; kwargs...)
        rm2 = RecurrenceMatrix(y, ε; kwargs...)
    else
        rm1 = RecurrenceMatrix(x[1:n,:], ε; kwargs...)
        rm2 = RecurrenceMatrix(y[1:n,:], ε; kwargs...)
    end
    return JointRecurrenceMatrix(rm1.data .* rm2.data)
end


################################################################################
# Scaling / fixed rate
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
        distances = zeros(Int(length(x)*(length(x)-1)/2))
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


################################################################################
# recurrence_matrix - Low level interface
################################################################################
# TODO: increase the efficiency here by not computing everything!

"""
    recurrence_matrix(x, y, metric::Metric, ε::Real, parallel::Val)

Return a sparse matrix which encodes recurrence points.

Note that `parallel` may be either `Val(true)` or `Val(false)`.
"""
function recurrence_matrix(x::AbstractMatrix, y::AbstractMatrix,
                           metric::Metric, ε, parallel::Val)
    # Convert Matrices to Datasets (better performance in all cases)
    return recurrence_matrix(Dataset(x), Dataset(y), metric, ε, parallel)
end

# First, we define the non-parallel versions.

# Core function
function recurrence_matrix(xx::Dataset, yy::Dataset, metric::Metric, ε, parallel::Val{false})
    x = xx.data
    y = yy.data
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
    for j in 1:length(y)
        nzcol = 0
        for i in 1:length(x)
            @inbounds if evaluate(metric, x[i], y[j]) ≤ ε
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
function recurrence_matrix(x::AbstractVector, y::AbstractVector, metric, ε, parallel::Val{false})
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
    for j in 1:length(y)
        nzcol = 0
        for i in 1:length(x)
            if @inbounds abs(x[i] - y[j]) ≤ ε
                push!(rowvals, i)
                nzcol += 1
            end
        end
        append!(colvals, fill(j, (nzcol,)))
    end
    nzvals = fill(true, (length(rowvals),))
    return sparse(rowvals, colvals, nzvals, length(x), length(y))
end

function recurrence_matrix(x::AbstractVector, metric::Metric, ε::Real, parallel::Val{false})
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
    for j in 1:length(x)
        nzcol = 0
        for i in 1:j
            if @inbounds abs(x[i] - x[j]) ≤ ε
                push!(rowvals, i)
                nzcol += 1
            end
        end
        append!(colvals, fill(j, (nzcol,)))
    end
    nzvals = fill(true, (length(rowvals),))
    return Symmetric(sparse(rowvals, colvals, nzvals, length(x), length(x)), :U)
end

function recurrence_matrix(xx::Dataset, metric::Metric, ε::Real, parallel::Val{false})
    x = xx.data
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
    for j in 1:length(x)
        nzcol = 0
        for i in 1:j
            @inbounds if evaluate(metric, x[i], x[j]) ≤ ε
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
function recurrence_matrix(xx::Dataset, yy::Dataset, metric::Metric, ε, parallel::Val{true})
    x = xx.data
    y = yy.data
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
            @inbounds if evaluate(metric, x[i], y[j]) ≤ ε
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
function recurrence_matrix(x::AbstractVector, y::AbstractVector, metric, ε, parallel::Val{true})
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
            @inbounds if abs(x[i] - y[j]) ≤ ε
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



function recurrence_matrix(xx::AbstractVector, metric::Metric, ε::Real, parallel::Val{true})
    x = xx
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
                @inbounds if abs(x[i] - x[j]) ≤ ε
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

function recurrence_matrix(xx::Dataset, metric::Metric, ε::Real, parallel::Val{true})
    x = xx.data
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
                @inbounds if evaluate(metric, x[i], x[j]) ≤ ε
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
