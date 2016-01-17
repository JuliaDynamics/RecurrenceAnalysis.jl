module RecurrenceAnalysis

using Distances

export embed,
       distancematrix,
       recurrencematrix,
       recurrencerate,
       determinism
       avgdiag,
       maxdiag,
       divergence,
       entropy,
       trend,
       laminarity,
       trappingtie,
       maxvert,
       autocorrelation,
       ami,
       gmi

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

# Embed series
function embed(x::AbstractVecOrMat, m::Integer, delay::Integer)
    dims = size(x)
    n = dims[1]
    nm = n-delay*(m-1)
    if nm < 2 &&
        warning("the emedded time series has length < 2")
    end
    ix = (1:nm) .+ (0:delay:delay*(m-1))'
    embed_indices(x, ix)
end

embed(x, m, delay) = isinteger(m) && isinteger(delay) ?
    embed(x,m, delay) : error("embedding dimension and delay must be integer values")

embed_indices(x::AbstractVector, indices) = x[indices]

function embed_indices(x::AbstractMatrix, indices)
    dx = size(x)
    dxm = size(indices)
    ix_rep = repeat(indices, inner=[1,dx[2]])
    ix_rep += repeat(dx[1]*(0:dx[2]-1)', outer=[dxm...])
    x[ix_rep]
end

# Recurrence matrix creation

"""Create a distance matrix from an embedded time series"""
function distancematrix(x, normtype="max")
    dist = getmetric(normtype)
    pairwise(dist, x')
end

"""Create a recurrence matrix from an embedded time series"""
function recurrencematrix(x, radius, args...)
    sparse(distancematrix(x, args...) .< radius)
end

include("rqa.jl")
include("delay.jl")
include("radius.jl")

end
