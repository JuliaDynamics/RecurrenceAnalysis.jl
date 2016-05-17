# Radius for a given median recurrence rate
function radius_mrr(x::AbstractVector, rr::Real)
    lx = length(x)
    nr = ceil(Integer, rr*lx)
    nr == 0 && error("the recurrence rate does not account for more than one sample")
    xs = sort(x)
    d = zeros(lx)
    for i=1:lx
        dxs = xs[max(1,i-nr+1):min(lx,i+nr-1)] - xs[i]
        d[i] = sort(abs(dxs))[nr]
    end
    median(d)
end

