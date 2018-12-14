# Radius for a given median recurrence rate
function radius_mrr(x::AbstractVector, rr::Real)
    lx = length(x)
    nr = ceil(Integer, rr*lx)
    nr == 0 && error("the recurrence rate does not account for more than one sample")
    xs = sort(x)
    d = zeros(lx)
    for i=1:lx
        dxs = xs[max(1,i-nr+1):min(lx,i+nr-1)] - xs[i]
        d[i] = sort(abs.(dxs))[nr]
    end
    median(d)
end

"""
    sorteddistances(x; kwargs...)

Return a tuple with the sorted distances between points of the
embedded time series `x`, and the recurrence rates under those values.

The keyword arguments are the same that should be passed to the functions
`recurrencematrix` and `recurrencerate` to obtain those results. I.e., if
`d,r = sorteddistances(x; kwargs...)`, and
`rmat = recurrencematrix(x, d[i]; kwargs...)`, then
`recurrencerate(rmat; kwargs...) == r[i]`.
"""
function sorteddistances(x; theiler::Integer=0, scale=1, kwargs...)
    # Create distance matrix
    kwargs = Dict(kwargs)
    argsdm = haskey(kwargs,:metric) ? (x, kwargs[:metric]) : (x,)
    dm = distancematrix(argsdm...)
    scale = (typeof(scale) <: Function) ? scale(dm) : scale
    dm /= scale
    # Browse upper triangle after Theiler window
    n = size(x,1)
    nd = (n+1)*n/2
    ntheiler = theiler*n - theiler*(theiler-1)/2
    distarray = zeros(round(Integer,nd-ntheiler))
    pos = 0
    for d = 1:n
        tmp = dm[d+theiler:n,d]
        distarray[pos .+ (1:length(tmp))] .= tmp
        pos += length(tmp)
    end
    sort(distarray), (1.:length(distarray))/length(distarray)
end
