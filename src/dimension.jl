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
    nrm = nr - delay*m
    ncm = nc - delay*m
    x = abs2(x)
    xm = x[1:nrm, 1:ncm]
    for k = 1:m
        xm = fun(xm, x[delay*k+(1:nrm), delay*k+(1:ncm)])
    end
    sqrt(xm)
end

# Kennel's algorithm
function fnn(x, mbounds, delay, thresholds; metric="max")
    Rtol2 = thresholds[1]^2
    Atol = thresholds[2]
    Ra = std(x)
    m1, m2 = mbounds
    nfnn = zeros(m2-m1+1)
    # DM - with maximum value in the diagonal to rule it out in neighbour search
    # Start with matrix embedded at m=m1
    dm = distancematrix(embed(x, m1, delay), metric)
    n = size(dm)[1]
    dm[1:n+1:end] = maximum(dm)
    for m = 1:length(nfnn)
        # Look for nearest neighbours within the first n-delay
        n -= delay
        nnval, nnpos = findmin(dm[1:n,1:n], 2)
        dm = embedmatrix1(dm, delay, metric)
        nnv1 = nnval[1:n]
        nnv2 = dm[1:n,1:n][nnpos]
        fnn1 = ( ((nnv2.^2)./(nnv1.^2)-1) .> Rtol2 )
        fnn2 = ( (nnv2/Ra) .> Atol )
        nfnn[m] = sum(fnn1 | fnn2)
    end
    nfnn
end

# Cao's algorithm
function afnn(x, mbounds, delay; metric="max")
    m1, m2 = mbounds
    dm = distancematrix(embed(x, m1, delay), metric)
    n = size(dm)[1]
    dm[1:n+1:end] = maximum(dm)
    mean_increment = zeros(m2-m1+2)
    mean_ratio = zeros(m2-m1+2)
    for m = 1:length(mean_ratio)
        n -= delay
        nnval, nnpos = findmin(dm[1:n,1:n], 2)
        dm = embedmatrix1(dm, delay, metric)
        nnv1 = nnval[1:n]
        nnv2 = dm[1:n,1:n][nnpos]
        mean_ratio[m] = mean(nnv2./nnv1)
        mean_increment[m] = mean(sqrt(nnv2.^2 - nnv1.^2))
    end
    e1 = mean_ratio[2:end]./mean_ratio[1:end-1]
    e2 = mean_increment[2:end]./mean_increment[1:end-1]
    e1, e2
end

# Krakovská's algorithm
function ffnn(x, bounds, delay; metric="max")
    m1, m2 = mbounds
    dm = distancematrix(embed(x, m1, delay), metric)
    n = size(x)[1]
    dm[1:n+1:end] = maximum(dm)
    fnn_ratio = zeros(m2-m1+1)
    for m = 1:length(fnn_ratio)
        n -= delay
        nnval, nnpos1 = findmin(dm[1:n,1:n], 2)
        dm = embedmatrix1(dm, delay, metric)
        nnval, nnpos2 = findmin(dm, 2)
        fnn_ratio[m] = sum(nnpos1 .!= nnpos2)/n
    end
    fnn_ratio
end
