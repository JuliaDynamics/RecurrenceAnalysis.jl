#=
In this file the core computations for creating a recurrence matrix
are defined (via multiple dispatch).

The low level interface is contained in the function
`recurrence_matrix`, and this is where any specialization should happen.
=#
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
* `fan::Bool=false` : a flag that indicates if `ε` should be
  taken as a target fixed global recurrence rate, but here also ensuring a local
  fixed recurrence rate, i.e. for each column of the recurrence matrix. This
  corresponds to preselcting a Fixed Amount of Neighbors (fan) for each point in
  phase space. If `fan` is set to `true`, `ε` must be a value between 0 and 1,
  and `scale` is ignored. If `fixedrate` is also set to `true` the fan-method is
  used.
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

"""
    JointRecurrenceMatrix(R1, R2; kwargs...)
Create a joint recurrence matrix from given recurrence matrices `R1, R2`.
"""
function JointRecurrenceMatrix(R1::AbstractRecurrenceMatrix, R2::AbstractRecurrenceMatrix; kwargs...)
    R3 = R1.data .* R2.data
    return JointRecurrenceMatrix(R3)
end

################################################################################
# Scaling / fixed rate
################################################################################
# here args... is (x, y, metric, ε) or just (x, metric, ε)
function resolve_scale(args...; scale=1, fixedrate=false, fan=false)
    ε = args[end]
    x = args[1]
    if length(args) > 3
        y = args[2]
        metric = args[3]
    else
        metric = args[2]
        y = x
    end
    # Check fixed recurrence rate - ε must be within (0, 1)
    if fixedrate && (fan == false)
        sfun = (m) -> quantile(vec(m), ε)
        return resolve_scale(Base.front(args)..., 1.0; scale=sfun, fixedrate=false)
    elseif (fixedrate && fan) || (fixedrate == false && fan)
        return get_fan_threshold(x, y, metric, ε)
    else
        scale_value = _computescale(scale, Base.front(args)...)
        return ε*scale_value
    end
end

"""
    get_fan_threshold(x, y, metric, ε) → ε_fan
Compute the adaptive Fixed Amount of Neibours (FAN) `ε_fan`
"""
function get_fan_threshold(x, y, metric, ε::Real)
    @assert length(x) == length(y)
    @assert 0 < ε < 1 "Global recurrence rate must be ∈ (0, 1)"
    fan_threshold = zeros(length(x))
    d = distancematrix(x, y, metric)

    for i in 1:size(d, 1)
        fan_threshold[i] = quantile(view(d, i ,:), ε)
    end
    return fan_threshold
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
function recurrence_matrix(xx::Dataset, yy::Dataset, metric::Metric, ε::Real, parallel::Val{false})
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

# for FAN threshold: not necessarily symmetric RP anymore
function recurrence_matrix(xx::Dataset, yy::Dataset, metric::Metric, ε::Vector, parallel::Val{false})
    x = xx.data
    y = yy.data
    @assert length(ε) == length(x)
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
    for j in 1:length(y)
        nzcol = 0
        for i in 1:length(x)
            @inbounds if evaluate(metric, x[i], y[j]) ≤ ε[i]
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
function recurrence_matrix(x::AbstractVector, y::AbstractVector, metric, ε::Real, parallel::Val{false})
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

# for FAN threshold: not necessarily symmetric RP anymore
function recurrence_matrix(x::AbstractVector, y::AbstractVector, metric, ε::Vector, parallel::Val{false})
    @assert length(ε) == length(x)
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
    for j in 1:length(y)
        nzcol = 0
        for i in 1:length(x)
            if @inbounds abs(x[i] - y[j]) ≤ ε[i]
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
            if @inbounds evaluate(metric, x[i], x[j]) ≤ ε
                push!(rowvals, i)
                nzcol += 1
            end
        end
        append!(colvals, fill(j, (nzcol,)))
    end
    nzvals = fill(true, (length(rowvals),))
    return Symmetric(sparse(rowvals, colvals, nzvals, length(x), length(x)), :U)
end

function recurrence_matrix(x::AbstractVector, metric::Metric, ε::Vector, parallel::Val{false})
    @assert length(ε) == length(x)
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
    for j in 1:length(x)
        nzcol = 0
        for i in 1:length(x)
            if @inbounds evaluate(metric, x[i], x[j]) ≤ ε[i]
                push!(rowvals, i)
                nzcol += 1
            end
        end
        append!(colvals, fill(j, (nzcol,)))
    end
    nzvals = fill(true, (length(rowvals),))
    return sparse(rowvals, colvals, nzvals, length(x), length(x))
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

# for FAN threshold: not necessarily symmetric RP anymore
function recurrence_matrix(xx::Dataset, metric::Metric, ε::Vector, parallel::Val{false})
    x = xx.data
    @assert length(ε) == length(x)
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
    for j in 1:length(x)
        nzcol = 0
        for i in 1:length(x)
            @inbounds if evaluate(metric, x[i], x[j]) ≤ ε[i]
                push!(rowvals, i)
                nzcol += 1
            end
        end
        append!(colvals, fill(j, (nzcol,)))
    end
    nzvals = fill(true, (length(rowvals),))
    return sparse(rowvals, colvals, nzvals, length(x), length(x))
end


# Now, we define the parallel versions of these functions.

# Core function
function recurrence_matrix(xx::Dataset, yy::Dataset, metric::Metric, ε::Real, parallel::Val{true})
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

function recurrence_matrix(xx::Dataset, yy::Dataset, metric::Metric, ε::Vector, parallel::Val{true})
    x = xx.data
    y = yy.data
    @assert length(ε) == length(x)
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
            @inbounds if evaluate(metric, x[i], y[j]) ≤ ε[i]
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
function recurrence_matrix(x::AbstractVector, y::AbstractVector, metric, ε::Real, parallel::Val{true})
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

function recurrence_matrix(x::AbstractVector, y::AbstractVector, metric, ε::Vector, parallel::Val{true})
    @assert length(ε) == length(x)
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
            @inbounds if abs(x[i] - y[j]) ≤ ε[i]
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

function recurrence_matrix(xx::AbstractVector, metric::Metric, ε::Vector, parallel::Val{true})
    x = xx
    @assert length(ε) == length(x)
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
            for i in 1:length(x)
                @inbounds if abs(x[i] - x[j]) ≤ ε[i]
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
    return sparse(finalrows, finalcols, nzvals, length(x), length(x))
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

function recurrence_matrix(xx::Dataset, metric::Metric, ε::Vector, parallel::Val{true})
    x = xx.data
    @assert length(ε) == length(x)
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
            for i in 1:length(x)
                @inbounds if evaluate(metric, x[i], x[j]) ≤ ε[i]
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
    return sparse(finalrows, finalcols, nzvals, length(x), length(x))
end
