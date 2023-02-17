#=
This file contains:
- Definition of recurrence matrices types and Base functions extensions
- Definition of types of recurrences
- Core computations for creating a recurrence matrix

The low level interface is contained in the function
`recurrence_matrix`, and this is where any specialization should happen.
=#
################################################################################
# AbstractRecurrenceMatrix type hierarchy
################################################################################
abstract type AbstractRecurrenceMatrix{RT} end
const ARM = AbstractRecurrenceMatrix

# The recurrence type is included as a field in all matrix types,
# because we use it in the pretty printing.
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

# Pretty printing:
function Base.summary(R::AbstractRecurrenceMatrix)
    N = nnz(R.data)
    return """
    $(nameof(typeof(R))) of size $(size(R.data)),
    with $(nameof(typeof(R.recurrence_type)))-type recurrences and $N entries.
    """
end
Base.show(io::IO, R::AbstractRecurrenceMatrix) = println(io, summary(R))

function Base.show(io::IO, ::MIME"text/plain", R::ARM)
    # `recurrenceplot` comes from `plot.jl` file.
    a = recurrenceplot(io, R)
    show(io, a)
end

################################################################################
# Definition of ways to find a recurrence
################################################################################
abstract type AbstractRecurrenceType end

"""
    RecurrenceThreshold(ε::Real)
Recurrences are defined as any point with distance `≤ ε` from the referrence point.
See [`RecurrenceMatrix`](@ref) for more.
"""
struct RecurrenceThreshold{T<:Real} <: AbstractRecurrenceType
    ε::T
end

"""
    RecurrenceThresholdScaled(ratio::Real, scale::Function)
Recurrences are defined as any point with distance `≤ d` from the referrence point,
where `d` is a scaled ratio (specified by `ratio, scale`) of the distance matrix.
See [`RecurrenceMatrix`](@ref) for more.
"""
struct RecurrenceThresholdScaled{T<:Real, S<:Function} <: AbstractRecurrenceType
    ratio::T
    scale::S
end

"""
    GlobalRecurrenceRate(rate::Real)
Recurrences are defined as a constant global recurrence rate.
See [`RecurrenceMatrix`](@ref) for more.
"""
struct GlobalRecurrenceRate{T<:Real} <: AbstractRecurrenceType
    rate::T
end


"""
    LocalRecurrenceRate(rate::Real)
Recurrences are defined as a constant local recurrence rate.
See [`RecurrenceMatrix`](@ref) for more.
"""
struct LocalRecurrenceRate{T<:Real} <: AbstractRecurrenceType
    rate::T
end

################################################################################
# Extensions of standard library
################################################################################
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

LinearAlgebra.issymmetric(R::ARM) = issymmetric(R.data)
LinearAlgebra.issymmetric(::RecurrenceMatrix{X}) where
    {X <: Union{RecurrenceThreshold,RecurrenceThresholdScaled,GlobalRecurrenceRate}} = true
LinearAlgebra.issymmetric(::JointRecurrenceMatrix{X}) where
    {X <: Union{RecurrenceThreshold,RecurrenceThresholdScaled,GlobalRecurrenceRate}} = true


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
# Documentation strings and dispatch to `recurrence_matrix`
################################################################################
"""
    RecurrenceMatrix(x, ε::Real; metric = Euclidean(), parallel::Bool)

Create a recurrence matrix from trajectory `x` (either a `Dataset` or a `Vector`)
and with recurrence distance threshold `ε`.
Instead of `ε::Real`, you can specify an `AbstractRecurrenceType`, see the method below.

The keyword `metric`, if given, must be any subtype of `Metric` from
[Distances.jl](https://github.com/JuliaStats/Distances.jl)
and defines the metric used to calculate distances for recurrences.
By default the Euclidean metric is used, typical alternatives are `Chebyshev(), Cityblock()`.

The keyword `parallel` decides if the comptutation should be done in parallel using threads.
Defaults to `length(x) > 500 && Threads.nthreads() > 1`.

## Description
A (cross-)recurrence matrix is a way to quantify *recurrences* that occur in a trajectory.
A recurrence happens when a trajectory visits the same neighborhood on the state space that
it was at some previous time.

The recurrence matrix is a numeric representation of a recurrence plot,
described in detail in [^Marwan2007] and [^Marwan2015]. It represents a
a sparse square matrix of Boolean values that quantifies recurrences in the trajectory,
i.e., points where the trajectory returns close to itself.
Given trajectories ``x, y``, the matrix is defined as:
```math
R[i, j] = \\begin{cases}
1 \\quad \\text{if}\\quad d(x[i], y[j]) \\le \\varepsilon\\
0 \\quad \\text{else}
\\end{cases}
```
with ``d`` is the distance function specified by the given `metric`.
The difference between a `RecurrenceMatrix` and a [`CrossRecurrenceMatrix`](@ref)
is that in the first case `x === y`.

Objects of type `<:AbstractRecurrenceMatrix` are displayed as a [`recurrenceplot`](@ref).

See also: [`CrossRecurrenceMatrix`](@ref), [`JointRecurrenceMatrix`](@ref) and
use [`recurrenceplot`](@ref) to turn the result of these functions into a plottable format.

[^Marwan2007]:
    N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems",
    [Phys. Reports 438*(5-6), 237-329 (2007)](https://doi.org/10.1016/j.physrep.2006.11.001)

[^Marwan2015]:
    N. Marwan & C.L. Webber, *Recurrence Quantification Analysis. Theory and Best Practices*
    [Springer (2015)](https://link.springer.com/book/10.1007/978-3-319-07155-8)
"""
function RecurrenceMatrix(x, ε::Real;
    # DEPRECATED keywords. TODO: Remove them in next stable release.
    scale = nothing, fixedrate = nothing,
    # VALID keywords (propagated)
    kwargs...)
    if !isnothing(scale)
        @warn "Providing keyword `scale` is deprecated. Use `RecurrenceThresholdScaled`."
        rt = RecurrenceThresholdScaled(ε, scale)
    elseif !isnothing(fixedrate)
        @warn "Providing keyword `fixedrate` is deprecated. Use `GlobalRecurrenceRate`."
        rt = GlobalRecurrenceRate(ε)
    elseif ε isa Real
        rt = RecurrenceThreshold(ε)
    else
        rt = ε
    end
    return RecurrenceMatrix(x, rt; kwargs...)
end

"""
    RecurrenceMarix(x, recurrence_type::AbstractRecurrenceType; metric, parallel)
Use this method to specify with more options how recurrences are identified and counted.
`metric, parallel` are as in the above method.

The `recurrence_type` can be:
* `RecurrenceThreshold(ε::Real)`: Recurrences are defined as any point with distance
  `≤ ε` from the referrence point. This is identical to `RecurrenceMatrix(x, ε::Real)`.
* `RecurrenceThresholdScaled(ratio::Real, scale::Function)`: Here `scale` is a function
  of the distance matrix `dm` (see [`distancematrix`](@ref)) that is used
  to scale the value of the recurrence threshold `ε` so that `ε = ratio*scale(dm)`,
  with `ratio` ∈ (0, 1). After the new `ε` is obtained, the method works
  just like the `RecurrenceThreshold`.
* `GlobalRecurrenceRate(ratio::Real)`: Here the number of total recurrence rate over the whole
  matrix (see [`recurrencerate`](@ref) is specified to be a `ratio` ∈ (0,1). This means that
  a distance threshold `ε` will be calculated such that there is a fixed `ratio` of
  recurrences, out of the total possible `N^2` (with `N = length(x)`).
  After the new `ε` is obtained, the method works just like the `RecurrenceThreshold`.
* `LocalRecurrenceRate(r::Real)`: The recurrence threhsold here is point-dependent. It is
  defined so that each point of `x` has a fixed number of `k = r*N` neighbors, with ratio
  `r` out of the total possible `N`. Equivalently, this means that each column of the
  recurrence matrix will have exactly `k` true entries. Notice that `LocalRecurrenceRate`
  does not guarantee that the resulting recurrence matrix will be symmetric.
"""
function RecurrenceMatrix(x, rt::AbstractRecurrenceType; metric = Euclidean(),
        parallel::Bool = length(x) > 500 && Threads.nthreads() > 1
    )
    metric = getmetric(metric) # TODO: Remove this in major update.
    ε = recurrence_threshold(rt, x, metric)
    m = recurrence_matrix(x, metric, ε, Val(parallel))
    return RecurrenceMatrix{typeof(rt)}(m, rt)
end


"""
    CrossRecurrenceMatrix(x, y, ε; kwargs...)

Create a cross recurrence matrix from trajectories `x` and `y`.
See [`RecurrenceMatrix`](@ref) for possible value sfor `ε` and `kwargs`.

The cross recurrence matrix is a bivariate extension of the recurrence matrix.
For the time series `x`, `y`, of length `n` and `m`, respectively, it is a
sparse `n×m` matrix of Boolean values, such that if `d(x[i], y[j]) ≤ ε`,
then the cell `(i, j)` of the matrix will have a `true` value.

Note that cross recurrence matrices are generally not symmetric irrespectively of `ε`.
"""
function CrossRecurrenceMatrix(x, y, ε;
        # DEPRECATED keywords. TODO: Remove them in next stable release.
        scale = nothing, fixedrate = nothing,
        # Normal keywords
        metric = Euclidean(), parallel::Bool = length(x) > 500 && Threads.nthreads() > 1
    )
    if ε isa AbstractRecurrenceType
        rt = ε
    elseif !isnothing(scale)
        @warn "Providing keyword `scale` is deprecated. Use `RecurrenceThresholdScaled`."
        rt = RecurrenceThresholdScaled(ε, scale)
    elseif !isnothing(fixedrate)
        @warn "Providing keyword `fixedrate` is deprecated. Use `GlobalRecurrenceRate`."
        rt = GlobalRecurrenceRate(ε)
    elseif ε isa Real
        rt = RecurrenceThreshold(ε)
    else
        throw(ArgumentError("Unknown specification of recurrence type."))
    end
    metric = getmetric(metric) # TODO: Remove this in major update.
    ε = recurrence_threshold(rt, x, y, metric)
    m = recurrence_matrix(x, y, metric, ε, Val(parallel))
    return CrossRecurrenceMatrix{typeof(rt)}(m, rt)
end

"""
    JointRecurrenceMatrix(x, y, ε; kwargs...)

Create a joint recurrence matrix from trajectories `x` and `y`.
See [`RecurrenceMatrix`](@ref) for possible value sfor `ε` and `kwargs`.

The joint recurrence matrix considers the recurrences of the trajectories
of `x` and `y` separately, and looks for points where both recur
simultaneously. It is calculated by the element-wise multiplication
of the recurrence matrices of `x` and `y`. If `x` and `y` are of different
length, the recurrences are only calculated until the length of the shortest one.

See [`RecurrenceMatrix`](@ref) for details, references and keywords.
"""
function JointRecurrenceMatrix(x, y, ε; kwargs...)
    n = min(length(x), length(y))
    rm1 = RecurrenceMatrix(x[1:n], ε; kwargs...)
    rm2 = RecurrenceMatrix(y[1:n], ε; kwargs...)
    ε = ε isa Real ? RecurrenceThreshold(ε) : ε # to be sure we have recurrence type
    return JointRecurrenceMatrix{typeof(ε)}(rm1.data .* rm2.data, ε)
end

"""
    JointRecurrenceMatrix(R1::AbstractRecurrenceMatrix, R2::AbstractRecurrenceMatrix)

Equivalent with `R1 .* R2`.
"""
function JointRecurrenceMatrix(
        R1::AbstractRecurrenceMatrix{RT}, R2::AbstractRecurrenceMatrix{RT}
    ) where {RT}
    R3 = R1.data .* R2.data
    return JointRecurrenceMatrix{RT}(R3, R1.rt)
end

# Extend methods that just take in directly a boolean matrix
for type in (:RecurrenceMatrix, :CrossRecurrenceMatrix, :JointRecurrenceMatrix)
    @eval begin
        function $(type)(m::AbstractMatrix{Bool})
            x = SparseMatrixCSC(m)
            return $(type){RecurrenceThreshold}(x, RecurrenceThreshold(0))
        end
    end
end

################################################################################
# Resolving the recurrence threshold and/or scaling
################################################################################
"""
    recurrence_threshold(rt::AbstractRecurrenceType, x [, y], metric) → ε
Return the calculated threshold `ε` for `rt`. The output is real, unless
`rt isa LocalRecurrenceRate`, where `ε isa Vector`.
"""
recurrence_threshold(rt, x, metric::Metric) =
    recurrence_threshold(rt, x, x, metric)
recurrence_threshold(rt::Real, x, y, metric) = rt
recurrence_threshold(rt::RecurrenceThreshold, x, y, metric) = rt.ε
function recurrence_threshold(rt::RecurrenceThresholdScaled, x, y, metric)
    scale_value = _computescale(rt.scale, x, y, metric)
    return rt.ratio*scale_value
end
function recurrence_threshold(rt::GlobalRecurrenceRate, x, y, metric)
    # We leverage the code of the scaled version to get the global recurrence rate
    scale = (m) -> quantile(vec(m), rt.rate)
    scale_value = _computescale(scale, x, y, metric)
    return rt.rate*scale_value
end
function recurrence_threshold(rt::LocalRecurrenceRate, x, y, metric)
    rate = rt.rate
    @assert 0 < rate < 1 "Recurrence rate must be ∈ (0, 1)."
    thresholds = zeros(length(y))
    d = distancematrix(x, y, metric)
    if x === y
        rate += 1/length(x) # because of recurrences guaranteed in the diagonal
    end
    for i in axes(d, 2)
        thresholds[i] = quantile(view(d, : ,i), rate)
    end
    return thresholds
end


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


################################################################################
# recurrence_matrix - Low level function
################################################################################
# TODO: increase the efficiency here by not computing everything!

# Core function

# First, we define the non-parallel versions.

# For two datasets

recurrence_matrix(x, metric, ε, parallel) = recurrence_matrix(x, x, metric, ε, parallel)

"""
    recurrence_matrix(x [, y], metric::Metric, ε, parallel::Val)

Return a sparse matrix which encodes recurrence points.
`ε` is the `recurrence_threhsold`, which is a vector for `LocalRecurrenceRate`.

`parallel` may be either `Val(true)` or `Val(false)`.
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
