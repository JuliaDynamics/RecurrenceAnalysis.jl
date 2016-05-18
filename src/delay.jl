"""
    autocorrelation(x)
    
Calculate autocorrelation of a vector for all positive lags.
"""
autocorrelation(x) = xcorr(x,x)[length(x):end]

# Numbers of bin
# Freeman-Diaconis rule
function n_fd(x::AbstractVector)
    first, last = extrema(x)
    iqr = quantile(x,.75) - quantile(x,.25)
    ceil(Integer, (last-first) / (2iqr * length(x)^(-1/3)) )
end

# Sturges
n_sturges(x::AbstractVector) = ceil(Integer, 1+log2(length(x)))

makebins(x, n) = linspace(extrema(x)..., n+1)

"""
    ami(x, delay, nbins)
    
Calculate the average mutual information of a signal for given delays.

The delay can either be an integer (in which case the function returns an integer),
or a series of numbers (returning a vector) specified as an array of ascending integers,
a range, or a tuple with the limits of the range (minimum, maximum).

The number of bins used to estimate the marginal and joint entropies of the
original and delayed signal can be given as a fixed integer, or as a
string that specifies a binning criterion. Currently available criteria are
Sturges (`"Sturges"`) or Freeman-Diaconis (`"FD"`).
"""
function ami(x, delay::Integer, nbins::Integer)
    n = length(x)-delay
    log2n = log2(n)
    edges = makebins(x, nbins)
    b1, b2, pxy = hist2d(x[(1:n) .+ [0 delay]], edges, edges)
    px = sum(pxy, 2)
    py = sum(pxy, 1)
    ret = 0.
    for ix=1:nbins, iy=1:nbins
        if pxy[ix,iy] != 0
            pxpy = px[ix]*py[iy]
            ret += pxy[ix,iy]/n*(log2(pxy[ix,iy]/pxpy)+log2n)
        end
    end
    # Maximum entropy for a given nbins is log2(nbins)
    ret/log2(nbins)
end

function ami(x, delay::Union{Array, Range}, nbins::Integer)
    d2 = maximum(delay)
    miv = zeros(d2+1)
    n = length(x)-d2
    log2n = log2(n)
    # Pre-allocated arrays
    pxy = zeros(nbins, nbins)
    px = zeros(nbins,1)
    py = zeros(1,nbins)
    for d in delay
        edges = makebins(x[1:n+d], nbins)
        b1, b2, pxy = hist2d!(pxy, x[(1:n) .+ [0 d]], edges, edges)
        px = sum!(px, pxy)
        py = sum!(py, pxy)
        for ix=1:nbins, iy=1:nbins
            if pxy[ix,iy] != 0
                pxpy = px[ix]*py[iy]
                miv[d+1] += pxy[ix,iy]/n*(log2(pxy[ix,iy]/pxpy)+log2n)
            end
        end
    end
    miv[1+delay]/log2(nbins)
end

ami(x, delay::Tuple{Integer,Integer}, nbins) = ami(x, colon(delay...), nbins)

function ami(x, delay, nbins="Sturges")
    # Define number of bins
    dict_binfuns = Dict("Sturges" => n_sturges,
        "FD" => n_fd)
    !haskey(dict_binfuns,nbins) && error("incorrect definition of bins.")
    binfun = dict_binfuns[nbins]
    nbins = binfun(x)
    ami(x, delay, nbins)
end

"""
    gmi(x, delay, radius; <keyword arguments>)
    
Calculate the generalized mutual information of a signal for given delays, based
on Rényi entropies.

The delay can either be an integer (in which case the function returns an integer),
or a series of numbers (returning a vector) specified as an array of ascending integers,
a range, or a tuple with the limits of the range (minimum, maximum).

The radius is a fixed value in the absolute scale of the signal,
used for the calculation of recurrence matrices to estimate
Rényi entropies (see `?recurrencematrix` for details).
Maximum norm and no scaling of distances are assumed.
"""
function gmi(x, delay::Integer, radius::Real)
    rmfull = recurrencematrix(x,radius,scale=1)
    n = length(x) - delay
    # Entropy of the non-delayed series
    rm1 = rmfull[1:n, 1:n]
    h2 = -log2(nnz(rm1)/n^2)
    # Joint entropy
    rmj = rm1 .* rmfull[delay+(1:n), delay+(1:n)]
    h2tau = -log2(nnz(rmj)/n^2)
    # Maximum entropy for a given radius corresponds to uniform distribution
    maxh2 = log2((maximum(x)-minimum(x))/radius)
    (2h2 - h2tau)/maxh2
end

function gmi(x, delay::Union{Array, Range}, radius::Real)
    d2 = maximum(delay)
    rmfull = recurrencematrix(x,radius,scale=1)
    h2tau = zeros(d2+1)
    n = length(x) - d2
    # Entropy of the non-delayed series
    rm1 = rmfull[1:n, 1:n]
    h2 = -log2(nnz(rm1)/n^2)
    # Joint entropy
    rmj = spzeros(n, n)
    for d in delay
        rmj = rm1 .* rmfull[d+(1:n), d+(1:n)]
        h2tau[d+1] = -log2(nnz(rmj)/n^2)
    end
    # Maximum entropy for a given radius corresponds to uniform distribution
    maxh2 = log2((maximum(x)-minimum(x))/radius)
    (2h2 - h2tau[1+delay])/maxh2
end

gmi(x, delay::Tuple{Integer,Integer}, radius) = gmi(x, colon(delay...), radius)

# gmi(x::AbstractVector, delay) = gmi(x, delay, radius_mrr(x, .01))
