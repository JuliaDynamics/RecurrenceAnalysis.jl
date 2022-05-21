#=
In this file the core computations for creating a recurrence matrix
are defined (via multiple dispatch).

The low level interface is contained in the function
`recurrence_matrix`, and this is where any specialization should happen.
=#
################################################################################
# AbstractRecurrenceMatrix type hierarchy and extensions of methods
################################################################################
const FAN = NeighborNumber

abstract type AbstractRecurrenceMatrix{RT} end
const ARM = AbstractRecurrenceMatrix

struct RecurrenceMatrix{RT} <: AbstractRecurrenceMatrix{RT}
    data::SparseMatrixCSC{Bool,Int}
    recurrence_type::RT
end
struct CrossRecurrenceMatrix{RT} <: AbstractRecurrenceMatrix{RT}
    data::SparseMatrixCSC{Bool,Int}
    recurrence_type::RT
end
struct JointRecurrenceMatrix{RT} <: AbstractRecurrenceMatrix{RT}
    data::SparseMatrixCSC{Bool,Int}
    recurrence_type::RT
end

function Base.summary(R::AbstractRecurrenceMatrix)
    N = nnz(R.data)
    return "$(nameof(typeof(R))) of size $(size(R.data)) with $N entries"
end
Base.show(io::IO, R::AbstractRecurrenceMatrix) = println(io, summary(R))

# Propagate used functions:
begin
    extentions = [
        (:Base, (:getindex, :size, :length, :view, :iterate,
            :eachindex, :axes, :CartesianIndices)),
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
LinearAlgebra.issymmetric(::ARM) = false

# column values in sparse matrix (parallel to rowvals)
function colvals(x::SparseMatrixCSC)
    cv = zeros(Int, nnz(x))
    @inbounds for c in axes(x,2)
        cv[nzrange(x,c)] .= c
    end
    cv
end
colvals(x::ARM) = colvals(x.data)

# Convert to matrices
Base.Array(R::ARM) = Matrix(R.data)
Base.Matrix(R::ARM) = Matrix(R.data)
SparseArrays.SparseMatrixCSC(R::ARM) = SparseMatrixCSC(R.data)

################################################################################
# Definition of ways to find a recurrence
################################################################################
abstract type AbstractRecurrenceType end

"""
    RecurrenceThreshold(ε::Real)
Recurrences are defined as any point with distance `< ε` from the referrence point.
See [`RecurrenceMatrix`](@ref) for more.
"""
struct RecurrenceThreshold{T<:Real} <: AbstractRecurrenceType
    ε::T
end

"""
    RecurrenceThreshold(r::Real, scale)
Recurrences are defined as any point with distance `< d` from the referrence point,
where `d` is a scaled ratio (specified by `scale`) of the distance matrix.
See [`RecurrenceMatrix`](@ref) for more.
"""
struct RecurrenceThresholdScaled{T<:Real, S} <: AbstractRecurrenceType
    ε::T
    scale::S
end

"""
    GlobalRecurrenceRate(r::Real)
Recurrences are defined as a constant global recurrence rate `r`.
See [`RecurrenceMatrix`](@ref) for more.
"""
struct GlobalRecurrenceRate{T<:Real} <: AbstractRecurrenceType
    r::T
end


"""
    LocalRecurrenceRate(r::Real)
Recurrences are defined as a constant local recurrence rate `r`.
See [`RecurrenceMatrix`](@ref) for more.
"""
struct LocalRecurrenceRate{T<:Real} <: AbstractRecurrenceType
    r::T
end


################################################################################
# Concrete Implementations & Documentation
################################################################################
"""
    RecurrenceMatrix(x, ε::Real; kwargs...)

Create a recurrence matrix from trajectory `x` (either a `Dataset` or a `Vector`)
and with distance threshold `ε`.
Objects of type `<:AbstractRecurrenceMatrix` are displayed as a [`recurrenceplot`](@ref).

See also: [`CrossRecurrenceMatrix`](@ref), [`JointRecurrenceMatrix`](@ref) and
use [`recurrenceplot`](@ref) to turn the result of these functions into a plottable format.

## Keyword Arguments
* `metric = "euclidean"` : metric of the distances, either `Metric` or a string,
   as in [`distancematrix`](@ref).
* `scale = 1` : a function of the distance matrix (see [`distancematrix`](@ref)),
  or a fixed number, used to scale the value of `ε`. Typical choices are
  `maximum` or `mean`, such that the threshold `ε` is defined as a ratio of the
  maximum or the mean distance between data points, respectively (using
  `mean` or `maximum` calls specialized versions that are faster than the naive
  approach).  Use `1` to keep the distances unscaled (default).
* `fixedrate::Bool = false` : a flag that indicates if `ε` should be
  taken as a target fixed recurrence rate (see [`recurrencerate`](@ref)).
  If `fixedrate` is set to `true`, `ε` must be a value between 0 and 1,
  and `scale` is ignored.
* `parallel::Bool` : whether to parallelize the computation of the recurrence
   matrix using threads. Default depends on available threads and length of `x`.

## Description

The recurrence matrix is a numeric representation of a recurrence plot,
described in detail in [^Marwan2007] and [^Marwan2015]. It represents a square
a sparse square matrix of Boolean values that quantifies recurrences in the trajectory,
i.e., points where the trajectory returns close to itself.
If `d(x[i], x[j]) ≤ ε` (with `d` the distance function),
then the entry `(i, j)` of the matrix will have a `true`
value. The criteria to evaluate distances between data points are defined
based on the keyword arguments.

[^Marwan2007]:
    N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems",
    [Phys. Reports 438*(5-6), 237-329 (2007)](https://doi.org/10.1016/j.physrep.2006.11.001)

[^Marwan2015]:
    N. Marwan & C.L. Webber, *Recurrence Quantification Analysis. Theory and Best Practices*
    [Springer (2015)](https://link.springer.com/book/10.1007/978-3-319-07155-8)
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

"""
    RecurrenceMatrix{FAN}(x, r::Real; metric = Euclidean(), parallel::Bool)

The parametrized constructor `RecurrenceMatrix{FAN}` creates the recurrence matrix
with a fixed number of neighbors for each point in the phase space, i.e. the number
of recurrences is the same for all columns of the recurrence matrix.
In such case, `r` is taken as the target fixed local recurrence rate,
defined as a value between 0 and 1 which means that `k = round(Int, r*N)` recurrences
for each point are identified with `N = length(x)`.
This is often referred to in the literature as the method of
"Fixed Amount of Nearest Neighbors" (or FAN for short).

Note that while recurrence matrices
with neighbors defined within a given range are always symmetric, those defined
by a fixed amount of neighbors can be non-symmetric.
"""
function RecurrenceMatrix{NeighborNumber}(x, ε; metric = DEFAULT_METRIC,
        parallel::Bool = length(x) > 500 && Threads.nthreads() > 1
    )
    m = getmetric(metric)
    s = get_fan_threshold(x, m, ε)
    m = recurrence_matrix(x, x, m, s, Val(parallel)) # to allow for non-symmetry
    return RecurrenceMatrix{NeighborNumber}(m)
end

"""
    CrossRecurrenceMatrix(x, y, ε::Real; kwargs...)
    CrossRecurrenceMatrix{FAN}(x, y, k::Int; kwargs...)

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
function JointRecurrenceMatrix{RT}(x, y, ε; kwargs...) where RT
    n = min(size(x,1), size(y,1))
    if n == size(x,1) && n == size(y,1)
        rm1 = RecurrenceMatrix{RT}(x, ε; kwargs...)
        rm2 = RecurrenceMatrix{RT}(y, ε; kwargs...)
    else
        rm1 = RecurrenceMatrix{RT}(x[1:n,:], ε; kwargs...)
        rm2 = RecurrenceMatrix{RT}(y[1:n,:], ε; kwargs...)
    end
    return JointRecurrenceMatrix{RT}(rm1.data .* rm2.data)
end

JointRecurrenceMatrix(args...; kwargs...) = JointRecurrenceMatrix{WithinRange}(args...; kwargs...)

"""
    JointRecurrenceMatrix(R1, R2; kwargs...)

Create a joint recurrence matrix from given recurrence matrices `R1, R2`.
"""
function JointRecurrenceMatrix(
    R1::AbstractRecurrenceMatrix{RT}, R2::AbstractRecurrenceMatrix{RT}; kwargs...
    ) where RT
    R3 = R1.data .* R2.data
    return JointRecurrenceMatrix{RT}(R3)
end

################################################################################
# Scaling / fixed rate / fixed amount of neighbors
################################################################################
# here args... is (x, y, metric, ε) or just (x, metric, ε)
function resolve_scale(x, y, metric, ε; scale=1, fixedrate=false)
    if fixedrate
        @assert 0 < ε < 1 "ε must be within (0, 1) for fixed rate"
        sfun = (m) -> quantile(vec(m), ε)
        return resolve_scale(x, y, metric, 1.0; scale=sfun, fixedrate=false)
    else
        scale_value = _computescale(scale, x, y, metric)
        return ε*scale_value
    end
end
resolve_scale(x, metric, ε) = resolve_scale(x, x, metric, ε)

# If `scale` is a function, compute the numeric value of the scale based on the
# distance matrix; otherwise return the value of `scale` itself
_computescale(scale::Real, args...) = scale

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
    for i in axes(d, 2)
        fan_threshold[i] = quantile(view(d, : ,i), ε)
    end
    return fan_threshold
end

get_fan_threshold(x, metric, ε) = get_fan_threshold(x, x, metric, ε)


################################################################################
# recurrence_matrix - Low level function
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
    for j in eachindex(y)
        nzcol = 0
        for i in eachindex(x)
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
    for j in eachindex(y)
        nzcol = 0
        for i in eachindex(x)
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
    for j in eachindex(x)
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
    for j in eachindex(x)
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
    Threads.@threads for j in eachindex(y)
        threadn = Threads.threadid()
        nzcol = 0
        for i in eachindex(x)
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
    Threads.@threads for j in eachindex(y)
        threadn = Threads.threadid()
        nzcol = 0
        for i in eachindex(x)
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
