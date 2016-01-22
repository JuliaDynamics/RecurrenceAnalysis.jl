# Embed distance matrix
function embedmatrix1(x, delay, metric="max")
    dist = getmetric(metric)
    fun = Dict(Euclidean()=>+, Chebyshev()=>max)[dist]
    nr, nc = size(x)
    x = abs2(x)
    sqrt(fun(x[1:nr-delay, 1:nc-delay], x[1+delay:nr, 1+delay:nc]))
end

function embedmatrix(x, m, delay, metric="max")
    dist = getmetric(metric)
    fun = Dict(Euclidean()=>+, Chebyshev()=>max)[dist]
    nr, nc = size(x)
    nrm = nr - delay*(m-1)
    ncm = nc - delay*(m-1)
    x = abs2(x)
    xm = x[1:nrm, 1:ncm]
    for k = 1:m-1
        xm = fun(xm, x[delay*k+(1:nrm), delay*k+(1:ncm)])
    end
    sqrt(xm)
end

# Kennel's algorithm
function fnn(x, m1, m2, delay, thresholds; metric="max")
    Rtol2 = thresholds[1]^2
    Atol = thresholds[2]
    Ra = std(x)
    nfnn = zeros(m2-m1+1)
    # DM - with maximum value in the diagonal to rule it out in neighbour search
    dm = distancematrix(x)
    n = size(x)[1]
    dm[1:n+1:end] = maximum(dm)
    # Start with matrix embedded at m=m1
    dme = embedmatrix(dm, m1, delay, metric)
    nm = n - delay*m2
    for m = 1:length(nfnn)
        # Look for nearest neighbours within the first nm points
        nnval, nnpos = findmin(dme[1:nm,1:nm], 2)
        dme = embedmatrix1(dme, delay, metric)
        nnv1 = nnval[1:nm]
        nnv2 = dme[1:nm,1:nm][nnpos]
        fnn1 = ( (nnv2./nnv1-1) .> Rtol2 )
        fnn2 = ( (nnv2/Ra) .> Atol )
        nfnn[m] = sum(fnn1 | fnn2)
    end
end

# Cao's algorithm
function afnn(x, m1, m2, delay; metric="max")
end
