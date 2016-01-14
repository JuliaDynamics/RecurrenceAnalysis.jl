# Autocorrelation

autocorrelation(x) = xcorr(x,x)[length(x):end]

# Bin according to Freeman-Diaconis
function fdbins(x::AbstractVector)
    first, last = extrema(x)
    iqr = quantile(x,.75) - quantile(x,.25)
    n = ceil( (last-first) / (2iqr * length(x)^(-1/3)) )
    linspace(first, last, n)
end

# Average mutual information

function ami(x, d1, d2)
    edges = fdbins(x)
    nbins = length(edges)-1
    ami = zeros(d2)
    n = length(x)-d2
    log2n = log2(n)
    pxy = zeros(nbins, nbins)
    px = zeros(nbins,1)
    py = zeros(1,nbins)
    for d in d1:d2
        b1, b2, pxy = hist2d!(pxy, x[(1:n) .+ [0 d]], edges, edges)
        px = sum!(px, pxy)
        py = sum!(py, pxy)
        for ix=1:nbins, iy=1:nbins
            if pxy[ix,iy] != 0
                pxpy = px[ix]*py[iy]
                ami[d] += pxy[ix,iy]*(log2(pxy[ix,iy]/pxpy)+log2n)
            end
        end
    end
    ami[d1:d2]
end

# Generalized mutual information
function gmi(x, d1, d2, radius)
    rmfull = recurrencematrix(x,radius)
    h2tau = zeros(d2)
    n = length(x) - d2
    rm1 = rmfull[1:n, 1:n]
    rmj = spzeros(n, n)
    for d in d1:d2 
        rmj = rm1 .* rmfull[d+(1:n), d+(1:n)]
        h2tau[d] = -log(nnz(rmj))
    end
    h2 = -log(recurrencerate(rm1))
    2h2 - h2tau[d1:d2] - 2log(n)
end

        
