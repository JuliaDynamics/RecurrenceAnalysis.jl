# get distance metric of the Distance packages
function getmetric(normtype::AbstractString)
    normtype = lowercase(normtype)
    metrics = Dict(
        "euclidean"=>Euclidean(),
        "max"=>Chebyshev(),
        "inf"=>Chebyshev(),
        "cityblock"=>Cityblock(),
        "manhattan"=>Cityblock(),
        "taxicab"=>Cityblock()
        )
    !haskey(metrics,normtype) && error("incorrect norm type. Accepted values are \""
        *join(keys(metrics),"\", \"", "\" or \"") * "\".")
    metrics[normtype]
end


# Distance matrices
# _maxdimension is used as threshold to choose the fastest method
const _maxdimension = 10
"""
    distancematrix(x[, y, metric])
    
Create a matrix with the distances between each pair of points of the
time series `x`, or between each point of `x` and `y`, with the distance type
defined by `metric`.

The time series `x` and `y` can be `Dataset`s or matrices with data points in rows.
The data point dimensions (or number of columns) must be the same for `x` and `y`.
The returned value is a `n`×`m` matrix, with `n` being the length (or number of rows)
of `x`, and `m` the length of `y`. If `y` is not provided, the result is a square
`n`×`n` matrix.

The metric can be identified by a string, or any of the `Metric`s defined in
the [Distances package](https://github.com/JuliaStats/Distances.jl).
The list of strings available to define the metric are:

* `"max"` or `"inf"` for the maximum or infinity norm 
  (`Chebyshev()` in the Distances package, used by default).
* `"euclidean"` for the L₂ or Euclidean norm
  (`Euclidean()` in Distances).
* `"manhattan"`, `"cityblock"` or `"taxicab"` for the Manhattan norm
  (`Cityblock()in Distances).
"""
function distancematrix(x, metric::Metric=Chebyshev())
    sx = size(x)
    if length(sx) > 1  &&  sx[2] > _maxdimension
        return _distancematrix(Matrix(x), metric)
    else
        return _distancematrix(Dataset(x), metric)
    end
end

function distancematrix(x, y, metric::Metric=Chebyshev())
    sx, sy = size(x), size(y)
    if length(sx) != length(sy)
        error("the dimension of `x` and `y` data points must be the equal")
    end
    if length(sx) > 1  &&  max(sx[2], sy[2]) > _maxdimension
        return _distancematrix(Matrix(x), Matrix(y), metric)
    else
        return _distancematrix(Dataset(x), Dataset(y), metric)
    end
end

# Methods specialized for matrices
_distancematrix(x::AbstractVecOrMat, metric::Metric) = pairwise(metric, x')
_distancematrix(x::AbstractVecOrMat, y::AbstractVecOrMat, metric::Metric) = pairwise(metric, x', y')
# Methods specialized for datasets
_distancematrix(x::Dataset, metric::Metric) = _distancematrix(x, x, metric)
# Core function for two datasets
function _distancematrix(x::D, y::D, metric::Metric) where {D<:Dataset{S,T}} where {S} where {T}
    x = x.data
    y = y.data
    d = zeros(T, length(x), length(y))
    for j in 1:length(y)
        for i in 1:length(x)
            @inbounds d[i,j] = evaluate(metric, x[i], y[j])
        end
    end
    return d
end

# Interface with strings for metrics
distancematrix(x, metric::String) = distancematrix(x, getmetric(metric))
distancematrix(x, y, metric::String) = distancematrix(x, y, getmetric(metric))


"""
    recurrencematrix(x, radius; <keyword arguments>)
    
Create a recurrence matrix from an embedded time series.

# Arguments:
* `x::Any` : embedded time series.
* `radius::Any` : threshold parameter to classify distances as recurrences.

# Keyword arguments:
* `scale::Any` : function of the distance matrix, or fixed number to scale the distances
between points. Typical choices are `maximum` to scale distances into
the unit inteval, or `mean`. Use `1` to keep the distances unscaled (default).
* `fixedrate::Bool` : a flag that indicates if the `radius` should be taken as a
target fixed recurrence rate of the full matrix (`false` by default).
If `fixedrate` is set to `true`, `radius` must be a value between 0 and 1 and `scale` is ignored.
* `metric::String` : metric of the norm, as in `distancematrix`.
"""
function recurrencematrix(x, radius; scale=1, fixedrate=false, kwargs...)
    kwargs = Dict(kwargs)
    # check fixed recurrence rate within (0,1)
    if fixedrate
        sfun = (m) -> quantile(m[:], radius)
        return recurrencematrix(x, 1; scale=sfun, fixedrate=false, kwargs...)
    end
    argsdm = haskey(kwargs,:metric) ? (x, kwargs[:metric]) : (x,)
    dm = distancematrix(argsdm...)
    (typeof(scale) <: Function) && (scale = scale(dm))
    sparse(dm .<= radius*scale)
end

"""
    crossrecurrencematrix(x, y, radius; <keyword arguments>)
    
Create a cross recurrence matrix from two embeded time series.

See `?recurrencematrix` for details.
"""
function crossrecurrencematrix(x, y, radius; scale=1, fixedrate=false, kwargs...)
    kwargs = Dict(kwargs)
    # check fixed recurrence rate within (0,1)
    if fixedrate
        sfun = (m) -> quantile(m[:], p)
        return crossrecurrencematrix(x, y, 1; scale=sfun, fixedrate=false, kwargs...)
    end
    argsdm = haskey(kwargs,:metric) ? (x, y, kwargs[:metric]) : (x, y)
    dm = distancematrix(argsdm...)
    (typeof(scale) <: Function) && (scale = scale(dm))
    sparse(dm .<= radius*scale)
end

"""
    crossrecurrencematrix(x, y, radius; <keyword arguments>)
    
Create a joint recurrence matrix from two embeded time series.

See `?recurrencematrix` for details.
"""
function jointrecurrencematrix(x, y, radius; kwargs...)
    rm1 = recurrencematrix(x, radius, kwargs...)
    rm2 = recurrencematrix(y, radius, kwargs...)
    n = minimum([size(rm1)..., size(rm2)...])
    rm1[1:n, 1:n] .* rm2[1:n, 1:n]
end

