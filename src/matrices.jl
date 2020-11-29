#=
In this file the core computations for creating a recurrence matrix
are defined (via multiple dispatch).

The low level interface is contained in the function
`recurrence_matrix`, and this is where any specialization should happen.
=#
################################################################################
# AbstractRecurrenceMatrix type definitions and documentation strings
################################################################################
struct FixedRange end
struct FixedAmount end
const FAN = FixedAmount

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

LinearAlgebra.issymmetric(::RecurrenceMatrix{FixedRange}) = true
LinearAlgebra.issymmetric(::JointRecurrenceMatrix{FixedRange}) = true
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
    RecurrenceMatrix{FAN}(...)

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
* `parallel::Bool=false` : whether to parallelize the computation of the recurrence
   matrix.  This will split the computation of the matrix across multiple threads.

The parametrized constructor `RecurrenceMatrix{FixedAmount}` (or in shorter form,
`RecurrenceMatrix{FAN}`) creates the recurrence matrix with a
Fixed Amount of Neighbors for each point in the phase space, i.e. the number
of recurrences is the same for all columns of the recurrence matrix.
In such case, `ε` is taken as the target fixed local recurrence rate,
defined as a value between 0 and 1, and `scale` and `fixedrate` are ignored.

If no parameter is specified, `RecurrenceMatrix` returns a
`RecurrenceMatrix{FixedRange}` object. Note that while recurrence matrices
with neighbors defined by a fixed range are always symmetric, those defined
by a fixed amount of neighbors can be non-symmetric.

See also: [`CrossRecurrenceMatrix`](@ref), [`JointRecurrenceMatrix`](@ref) and
use [`recurrenceplot`](@ref) to turn the result of these functions into a plottable format.

## References
[^1] : N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems",
*Phys. Reports 438*(5-6), 237-329 (2007).

[^2] : N. Marwan & C.L. Webber, "Mathematical and computational foundations of
recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence
Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).
"""
function RecurrenceMatrix{FixedRange}(x, ε; metric = DEFAULT_METRIC,
    parallel::Bool = length(x) > 500 && Threads.nthreads() > 1,
    kwargs...)
    m = getmetric(metric)
    s = resolve_scale(x, m, ε; kwargs...)
    m = recurrence_matrix(x, m, s, Val(parallel))
    return RecurrenceMatrix{FixedRange}(m)
end

RecurrenceMatrix(args...; kwargs...) = RecurrenceMatrix{FixedRange}(args...; kwargs...)

function RecurrenceMatrix{FixedAmount}(x, ε; metric = DEFAULT_METRIC,
    parallel::Bool = length(x) > 500 && Threads.nthreads() > 1,
    kwargs...)
    m = getmetric(metric)
    s = get_fan_threshold(x, m, ε)
    m = recurrence_matrix(x, x, m, s, Val(parallel)) # to allow for non-symmetry
    return RecurrenceMatrix{FixedAmount}(m)
end

"""
    CrossRecurrenceMatrix(x, y, ε; kwargs...)
    CrossRecurrenceMatrix{FAN}(...)

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
function CrossRecurrenceMatrix{FixedRange}(x, y, ε; metric = DEFAULT_METRIC,
    parallel::Bool = length(x) > 500 && Threads.nthreads() > 1,
    kwargs...)
    m = getmetric(metric)
    s = resolve_scale(x, y, m, ε; kwargs...)
    m = recurrence_matrix(x, y, m, s, Val(parallel))
    return CrossRecurrenceMatrix{FixedRange}(m)
end

CrossRecurrenceMatrix(args...; kwargs...) = CrossRecurrenceMatrix{FixedRange}(args...; kwargs...)

function CrossRecurrenceMatrix{FixedAmount}(x, y, ε; metric = DEFAULT_METRIC,
    parallel::Bool = length(x) > 500 && Threads.nthreads() > 1,
    kwargs...)
    m = getmetric(metric)
    s = get_fan_threshold(x, y, m, ε)
    m = recurrence_matrix(x, y, m, s, Val(parallel))
    return CrossRecurrenceMatrix{FixedAmount}(m)
end

"""
    JointRecurrenceMatrix(x, y, ε; kwargs...)
    JointRecurrenceMatrix{FAN}(...)

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

JointRecurrenceMatrix(args...; kwargs...) = JointRecurrenceMatrix{FixedRange}(args...; kwargs...)

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

"""
    recurrence_matrix(x, y, metric::Metric, ε, parallel::Val)

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
# For two datasets
function recurrence_matrix(xx::Dataset, yy::Dataset, metric::Metric, ε, ::Val{false})
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

function recurrence_matrix(xx::Dataset, metric::Metric, ε, ::Val{false})
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
function recurrence_matrix(xx::Dataset, yy::Dataset, metric::Metric, ε, ::Val{true})
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

function recurrence_matrix(xx::Dataset, metric::Metric, ε, ::Val{true})
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
