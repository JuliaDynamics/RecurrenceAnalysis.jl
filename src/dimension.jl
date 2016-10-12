# Embed distance matrix
function embedmatrix1{T}(x::Array{T,2}, delay, metric="max")
    dist = getmetric(metric)
    fun = Dict(Euclidean()=>+, Chebyshev()=>max)[dist]
    nr, nc = size(x)
    @compat x .= abs2.(x)
    @compat y = sqrt.(fun.(x[1:nr-delay, 1:nc-delay], x[1+delay:nr, 1+delay:nc]))
    Array{T,2}(y)
end

"""
    fnn(x, mbounds, delay, thresholds; metric="max")
    
Calculate the number of false nearest neighbours (FNN) by Kennel's algorithm.

FNN are calculated for the embedding dimensions between `mbounds[1]` and
`mbounds[2]`, and for the given delay, taking `thresholds[1]` as the parameter
Rtol, and `thresholds[2]` as the parameter Atol. The argument `metric` defines
what kind of distance is calculated between pairs of points, as in
`distancematrix`.
"""
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

"""
    afnn(x, mbounds, delay; metric="max")
    
Calculate the ratios E1 and E2 of Cao's method of "averaged false nearest neighbours".

The ratios are calculated for the embedding dimensions between `mbounds[1]` and
`mbounds[2]`, and for the given delay. The argument `metric` defines
what kind of distance is calculated between pairs of points, as in
`distancematrix`.
"""

function afnn(x, mbounds, delay; metric="max")
    m1, m2 = mbounds
    dm = distancematrix(embed(x, m1, delay), metric)
    n = size(dm)[1]
    dm[1:n+1:end] = maximum(dm)
    mean_ratio = zeros(m2-m1+2)
    mean_increment = zeros(m2-m1+2)
    for m = 1:length(mean_ratio)
        n -= delay
        nnval, nnpos = findmin(dm[1:n,1:n], 2)
        dm = embedmatrix1(dm, delay, metric)
        nnv1 = nnval[1:n]
        nnv2 = dm[1:n,1:n][nnpos]
        mean_ratio[m] = mean(nnv2./nnv1)
        # Project nnpos in original series
        nnx = div(nnpos, n) + 1
        d = (m1+m-1)*delay
        @compat mean_increment[m] = mean(abs.(x[(1:n)+d]-x[nnx+d]))
    end
    e1 = mean_ratio[2:end]./mean_ratio[1:end-1]
    e2 = mean_increment[2:end]./mean_increment[1:end-1]
    e1, e2
end

"""
    ffnn(x, mbounds, delay; metric="max")
    
Calculate the ratio of false first nearest neighbours (FFNN) by Krakovská's algorithm.

The ratios are calculated for the embedding dimensions between `mbounds[1]` and
`mbounds[2]`, and for the given delay. The argument `metric` defines
what kind of distance is calculated between pairs of points, as in
`distancematrix`.
"""
function ffnn(x, mbounds, delay; metric="max")
    m1, m2 = mbounds
    dm = distancematrix(embed(x, m1, delay), metric)
    n = size(dm)[1]
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
