# Radius for a given median recurrence rate
function radius_mrr(x::AbstractVector, rr::Real)
    nr = ceil(rr*length(x))
    nr == 0 && error("the recurrence rate does not account for more than one sample")
    xs = sort(x)
    d = xs[nr+1:end] - xs[1:end-nr]
    median(d)
end

