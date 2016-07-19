module RecurrenceAnalysis

using Distances
using StatsBase: fit, Histogram

export embed,
       distancematrix,
       recurrencematrix,
       crossrecurrencematrix,
       jointrecurrencematrix,
       recurrencerate,
       determinism,
       avgdiag,
       maxdiag,
       divergence,
       entropy,
       trend,
       laminarity,
       trappingtime,
       maxvert,
       autocorrelation,
       ami,
       gmi,
       fnn,
       afnn,
       ffnn

# get distance metric of the Distance packages
function getmetric(normtype::AbstractString)
    normtype = lowercase(normtype)
    metrics = Dict(
        "euclidean"=>Euclidean(),
        "max"=>Chebyshev(),
        "inf"=>Chebyshev()
        )
    !haskey(metrics,normtype) && error("incorrect norm type. Accepted values are \""
        *join(keys(metrics),"\", \"", "\" or \"") * "\".")
    metrics[normtype]
end

"""
    embed(x, m, delay)
    
Embed a time series in `m` dimensions with a given delay. 
"""
function embed(x::AbstractVecOrMat, m::Integer, delay::Integer)
    dims = size(x)
    n = dims[1]
    nm = n-delay*(m-1)
    (nm < 2) && warning("the emedded time series has length < 2")
    ix = (1:nm) .+ (0:delay:delay*(m-1))'
    embed_indices(x, ix)
end

embed(x, m, delay) = isinteger(m) && isinteger(delay) ?
    embed(x, m, delay) : error("embedding dimension and delay must be integer values")

embed_indices(x::AbstractVector, indices) = x[indices]

function embed_indices(x::AbstractMatrix, indices)
    dx = size(x)
    dxm = size(indices)
    ix_rep = repeat(indices, inner=[1,dx[2]])
    ix_rep += repeat(dx[1]*(0:dx[2]-1)', outer=[dxm...])
    x[ix_rep]
end

# Recurrence matrix creation

"""
    distancematrix(x, metric="max")
    distancematrix(x, y, metric="max")
    
Create a distance matrix from one or two (possibly embedded) time series.

Available metrics for the distances are `"max"` (default), `"inf"` (same), and `"euclidean"`.
"""
function distancematrix(x::AbstractVecOrMat, metric::AbstractString="max")
    dist = getmetric(metric)
    pairwise(dist, x')
end

function distancematrix(x::AbstractVecOrMat, y::AbstractVecOrMat, metric::AbstractString="max")
    dist = getmetric(metric)
    pairwise(dist, x', y')
end

"""
    recurrencematrix(x, radius; <keyword arguments>)
    
Create a recurrence matrix from an embedded time series.

# Arguments
* `x::Any` : embedded time series.
* `radius::Any` : threshold parameter to classify distances as recurrences.
* `scale::Any` : function of the distance matrix, or fixed number to scale the distances
between points. Typical choices are `maximum` (default) to scale distances into
the unit inteval, or `mean`. Use `1` to keep the distances unscaled.
* `metric::String` : metric of the norm, as in `distancematrix`.
"""
function recurrencematrix(x, radius; scale=maximum, kwargs...)
    kwargs = Dict(kwargs)
    argsdm = haskey(kwargs,:metric) ? (x, kwargs[:metric]) : (x,)
    dm = distancematrix(argsdm...)
    scfac = (typeof(scale) == Function) ? scale(dm[:]) : scale
    dm ./= scfac
    sparse(dm .< radius)
end

"""
    crossrecurrencematrix(x, y, radius; <keyword arguments>)
    
Create a cross recurrence matrix from two embeded time series.

See `?recurrencematrix` for details.
"""
function crossrecurrencematrix(x, y, radius; scale=maximum, kwargs...)
    kwargs = Dict(kwargs)
    argsdm = haskey(kwargs,:metric) ? (x, y, kwargs[:metric]) : (x, y)
    dm = distancematrix(argsdm...)
    scfac = (typeof(scale) == Function) ? scale(dm[:]) : scale
    dm ./= scfac
    sparse(dm .< radius)
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

include("rqa.jl")
include("delay.jl")
include("dimension.jl")
include("radius.jl")

end
