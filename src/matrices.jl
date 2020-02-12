using Statistics
using ChaosTools:linear_region
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
distancematrix(x, metric::Union{Metric,String}=DEFAULT_METRIC) = distancematrix(x, x, metric)

# For 1-dimensional arrays (vectors), the distance does not depend on the metric
distancematrix(x::Vector, y::Vector, metric=DEFAULT_METRIC) = abs.(x .- y')

# If the metric is supplied as a string, get the corresponding Metric from Distances
distancematrix(x, y, metric::String) = distancematrix(x, y, getmetric(metric))

const MAXDIM = 9
function distancematrix(x::Tx, y::Ty, metric::Metric=DEFAULT_METRIC) where
         {Tx<:Union{AbstractMatrix, Dataset}} where {Ty<:Union{AbstractMatrix, Dataset}}
    sx, sy = size(x), size(y)
    if sx[2] != sy[2]
        error("the dimensions of `x` and `y` data points must be the equal")
    end
    if sx[2] ≥ MAXDIM && typeof(metric) == Euclidean # Blas optimization
        return _distancematrix(Matrix(x), Matrix(y), metric)
    elseif Tx <: Matrix && Ty <: Matrix && metric == Chebyshev()
        return _distancematrix(x, y, metric)
    else
        return _distancematrix(Dataset(x), Dataset(y), metric)
    end
end

# Core function for Matrices (wrapper of `pairwise` from the Distances package)
_distancematrix(x::AbstractMatrix, y::AbstractMatrix, metric::Metric) =
pairwise(metric, x', y', dims=2)
# Core function for Datasets
function _distancematrix(x::Dataset{S,Tx}, y::Dataset{S,Ty},
    metric::Metric) where {S, Tx, Ty}

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

The recurrence matrix is a numeric representation of a "recurrence plot" [1, 2],
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

See also: [`CrossRecurrenceMatrix`](@ref), [`JointRecurrenceMatrix`](@ref) and
use [`recurrenceplot`](@ref) to turn the result of these functions into a plottable format.

## References
[1] : N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems",
*Phys. Reports 438*(5-6), 237-329 (2007).

[2] : N. Marwan & C.L. Webber, "Mathematical and computational foundations of
recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence
Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).
"""
function RecurrenceMatrix(x, ε; metric = DEFAULT_METRIC, kwargs...)
    m = getmetric(metric)
    if !haskey(kwargs, :scale)
        s = resolve_scale(x, x, m, ε; kwargs...)
    elseif kwargs[:scale] == "powerlaw"
        s = _powerlawepsilon(x, x; kwargs...)
    elseif kwargs[:scale] == "peak"
        s =_significantpeaksepsilon(x, x; kwargs...)
    end

    m = recurrence_matrix(x, m, s)
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
function CrossRecurrenceMatrix(x, y, ε; metric = DEFAULT_METRIC, kwargs...)
    m = getmetric(metric)

    if !haskey(kwargs, :scale)
        s = resolve_scale(x, y, m, ε; kwargs...)
    elseif kwargs[:scale] == "powerlaw"
        s = _powerlawepsilon(x, y; kwargs...)
    elseif kwargs[:scale] == "peak"
        s =_significantpeaksepsilon(x, y; kwargs...)
    end

    m = recurrence_matrix(x, y, m, s)
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

function _computedistances(x, y, metric::Metric)
    if x===y
        distlist = Vector{Real}(undef, Int(length(x)*(length(x)-1)/2))
        index = 1
        @inbounds for i in 1:length(x)-1, j=(i+1):length(x)
            distlist[index] = evaluate(metric, x[i], y[j])
            index += 1
        end
    else
        distlist = Vector{Real}(undef, length(x)*length(y))
        index = 1
        @inbounds for xi in x, yj in y
            distlist[index] += evaluate(metric, xi, yj)
            index += 1
        end
    end
    return distlist
end

function _computedistancematrix(x, y; metric::Metric)
    distmatrix = Matrix{Real}(undef, length(x), length(y))
    @inbounds for i in 1:length(x), j in 1:length(y)
        distmatrix[i, j] = evaluate(metric, x[i], y[j])
    end
    return distmatrix
end

function _computescale(scale::typeof(var), x, y, metric::Metric)
    return Statistics.var(_computedistances(x, y, metric))
end

function _computescale(scale::typeof(median), x, y, metric::Metric)
    return Statistics.mean(_computedistances(x, y, metric))
end


################################################################################
# recurrence_matrix - Low level interface
################################################################################
# TODO: increase the efficiency here by not computing everything:
recurrence_matrix(x, metric::Metric, ε::Real) =
recurrence_matrix(x, x, metric, ε)

# Convert Matrices to Datasets (better performance in all cases)
function recurrence_matrix(x::AbstractMatrix, y::AbstractMatrix,
                           metric::Metric, ε)
    return recurrence_matrix(Dataset(x), Dataset(y), metric, ε)
end

# Core function
function recurrence_matrix(xx::Dataset, yy::Dataset, metric::Metric, ε)
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
function recurrence_matrix(x::AbstractVector, y::AbstractVector, metric, ε)
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

################################################################################
# Power law method for selecting epsilon
################################################################################

function _powerlawepsilon(x, y; quantilelist=[1/2^n for n in 0:0.5:10], metric = DEFAULT_METRIC, kwargs...)
    rr = _computedistancematrix(x, y; metric=metric)
    rrvector = reshape(rr, length(rr))
    ϵlist = quantile(rrvector, quantilelist)

    logquantilelist = log.(quantilelist)
    logϵlist = log.(ϵlist)

    region, _ = linear_region(logϵlist, logquantilelist)

    return middle(region)
end


################################################################################
# Peak-counting method for selecting epsilon
################################################################################

function _diagonalhistogram(arr)
    x, y = size(arr)
    hist = Vector{typeof(arr[1])}(undef, x + y  - 1)
    for i in 1:x
        for j in 1:y
            hist[x+y-1] += arr[i, j]
        end
    end
    return hist
end

# From https://gist.github.com/usmcamp0811/ad3efd069755dc582df53b6cb0de3375
# Julia implimentation of http://stackoverflow.com/a/22640362/6029703
function _SmoothedZscoreAlgo(y, lag, threshold, influence)
    n = length(y)
    signals = zeros(n) # init signal results
    filteredY = copy(y) # init filtered series
    avgFilter = zeros(n) # init average filter
    stdFilter = zeros(n) # init std filter
    avgFilter[lag - 1] = mean(y[1:lag]) # init first value
    stdFilter[lag - 1] = std(y[1:lag]) # init first value

    for i in range(lag, stop=n-1)
        if abs(y[i] - avgFilter[i-1]) > threshold*stdFilter[i-1]
            if y[i] > avgFilter[i-1]
                signals[i] += 1 # postive signal
            else
                signals[i] += -1 # negative signal
            end
            # Make influence lower
            filteredY[i] = influence*y[i] + (1-influence)*filteredY[i-1]
        else
            signals[i] = 0
            filteredY[i] = y[i]
        end
        avgFilter[i] = mean(filteredY[i-lag+1:i])
        stdFilter[i] = std(filteredY[i-lag+1:i])
    end
    return signals
end

function _peakcounter(xs)
    count = 0
    for i in 1:length(xs)
        if xs[i] > 0 && !(xs[i-1] > 0)
            count += 1
        end
    end
    return count
end

function _significantpeaksepsilon(x, y; quantilelist=[1/2^n for n in 0:0.5:10], metric = DEFAULT_METRIC, lag=30, threshold=2, influence=0.7, kwargs...)
    rr = _computedistancematrix(x, y; metric=metric)
    rrvector = reshape(rr, length(rr))
    ϵlist = quantile(rrvector, quantilelist)

    bestscore = 100
    bestϵ = 0
    samplesize = size(rr)[1]
    binaryrr = zeros(size(rr))

    for i in 1:length(quantilelist)
        quantile = quantilelist[i]
        ϵ = ϵlist[i]

        map!(x -> x < ϵ ? 1.0 : 0.0, binaryrr, rr)
        xs = (_diagonalhistogram(binaryrr))
        peaklist = _SmoothedZscoreAlgo(binaryrr, lag, threshold, influence)
        peakcount = _peakcounter(peaklist)

        score = abs(1 - peakcount / (samplesize * quantile))
        if score > bestscore
            bestscore = score
            bestϵ = ϵ
        end
    end

    return bestϵ
end
