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
    nbins = fdbins(x)
    ami = zeros(d2)
    for d in d1:d2
        n = length(x) - d
        logn = log2(n)
        pxy = hist2d(x[(1:n) .+ [0 d]], nbins, nbins)[3]
        pxpy = sum(pxy, 2) * sum(pxy, 1)
        logp = ifelse(pxy.!=0, log2(pxy./pxpy)+logn, 0.)
        ami[d] = sum(pxy.*logp)
    end
    ami[d1:d2]
end

# Generalized mutual information
function gmi(x, d1, d2, radius)
    rmfull = recurrencematrix(x,radius)
    gmi = zeros(d2)
    n = length(x) - d2
    rm1 = rmfull[1:n, 1:n]
    h2 = -log(recurrencerate(rm1))
    for d in d1:d2
        rm2 = rmfull[d+(1:n), d+(1:n)]
        h2tau = -log(recurrencerate(rm1 .* rm2))
        gmi[d] = 2h2 - h2tau
    end
    gmi[d1:d2]
end

        
