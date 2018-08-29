# Approximation of diagonal line based measures as in Schultz et al. (2015)

# Re-embed a d-dimensional time series one further dimension, using a unit delay
reembed1(x,d) = @views hcat(x[1:end-1,:], x[2:end, (end-d+1):end])

# Number of pairwise proximities in an embedded time series
# Based on histogram of row values
function proximities(x::AbstractVecOrMat)
    bins = [@view x[1,:]]
    freqs = [1]
    for i in 2:size(x,1)
        row = @view x[i,:]
        if row in bins
            pos = findfirst(isequal(row), bins)
            freqs[pos] += 1
        else
            push!(bins, row)
            push!(freqs, 1)
        end
    end
    sum(freqs.^2)
end

# Initial routine: bin x according to the radius
function _bin(x, radius)
    if radius > 0
        delta = 2radius
        minx = minimum(x)
        typeint = typeof(1)
        return typeint.(div.(x.-minx, delta)) .+ 1
    else
        return x
    end
end

"""
    approx_recurrencerate(x, radius; loi=false)

Calculate an approximation to the recurrence rate (RR) of the time series `x`,
using `radius` as the maximum distance between recurrent points.
The LOI is excluded from the analysis by default.
The keyword argument `loi` can be set to `true` to include the LOI.
"""
function approx_recurrencerate(x, radius; loi::Bool=false)
    x = _bin(x, radius)
    pp = proximities(x)
    n = size(x,1)
    loi ? pp/n^2 : (pp-n)/n^2
end

"""
    approx_determinism(x, radius; lmin=2, loi=false)

Calculate an approximation to the determinism (DET) of the time series `x`,
using `radius` as the maximum distance between recurrent points,
and taking into account diagonal structures longer than `lmin`.
The LOI is excluded from the analysis by default.
The keyword argument `loi` can be set to `true` to include the LOI.
"""
function approx_determinism(x, radius; lmin=2, loi::Bool=false)
    x = _bin(x, radius)
    pp1 = proximities(x)
    xe = embed(x,lmin,1)
    ppm = proximities(xe)
    xex = reembed1(xe,size(x,2))
    ppmx = proximities(xex)
    n = size(x,1)
    loi ? (lmin*ppm - (lmin-1)*ppmx)/pp1 : (lmin*ppm - (lmin-1)*ppmx - n)/(pp1-n)
end

"""
    approx_avgdiag(x, radius; lmin=2, loi=false)

Calculate an approximation to the average length of diagonal recurrent structures
(L) of the time series `x`, using `radius` as the maximum distance between recurrent points,
and taking into account diagonal structures longer than `lmin`.
The LOI is excluded from the analysis by default.
The keyword argument `loi` can be set to `true` to include the LOI.
"""
function approx_avgdiag(x, radius; lmin=2, loi::Bool=false)
    x = _bin(x, radius)
    pp1 = proximities(x)
    xe = embed(x,lmin,1)
    ppm = proximities(xe)
    xex = reembed1(xe,size(x,2))
    ppmx = proximities(xex)
    n = size(x,1)
    loi ? (lmin*ppm - (lmin-1)*ppmx)/(ppm-ppmx) : (lmin*ppm - (lmin-1)*ppmx - n)/(ppm-ppmx-1)
end

"""
    approx_rqa(x, radius; lmin=2, loi=false)

Calculate an approximation to the RQA parameters of the time series `x` based
on diagonal structures, using `radius` as the maximum distance between recurrent points,
and taking into account diagonal structures longer than `lmin`.
The LOI is excluded from the analysis by default.
The keyword argument `loi` can be set to `true` to include the LOI.

The returned value is a dictionary with the following keys:

* "RR": recurrence rate (see `approx_recurrencerate`)
* "DET": determinsm (see `approx_determinism`)
* "L": average length of diagonal structures (see `approx_avgdiag`)
"""
function approx_rqa(x, radius; lmin=2, loi::Bool=false)
    x = _bin(x, radius)
    pp1 = proximities(x)
    xe = embed(x,lmin,1)
    ppm = proximities(xe)
    xex = reembed1(xe,size(x,2))
    ppmx = proximities(xex)
    n = size(x,1)
    Dict("RR"  => typeof(0.0)( loi ? pp1/n^2 : (pp1-n)/n^2 ),
        "DET"  => typeof(0.0)( loi ? (lmin*ppm - (lmin-1)*ppmx)/pp1 : (lmin*ppm - (lmin-1)*ppmx - n)/(pp1-n) ),
        "L"    => typeof(0.0)( loi ? (lmin*ppm - (lmin-1)*ppmx)/(ppm-ppmx) : (lmin*ppm - (lmin-1)*ppmx - n)/(ppm-ppmx-1) )
    )
end



###
# The following are not efficient - do not export

"""
    approx_maxdiag
"""
function approx_maxdiag(x, radius; lmin=2, loi::Bool=false)
    n = size(x,1)
    nd = size(x,2)
    # The longest diagonal is the LOI (if it is not excluded)
    loi && (return n)
    x = _bin(x, radius)
    # Search until the number of proximities is equal to the diagonal length
    m = lmin
    xe = embed(x,lmin,1)
    pp = proximities(xe)
    while pp > n-m+1
        m += 1
        xe = reembed1(xe,nd)
        pp = proximities(xe)
    end
    m-1
end

"""
    approx_divergence
"""
approx_divergence(x, radius; kwargs...) = typeof(0.0)( 1/approx_maxdiag(x, radius; kwargs...) )

"""
    approx_entropy
"""
function approx_entropy(x, radius; lmin=2, loi=false)
    n = size(x,1)
    nd = size(x,2)    
    m = lmin
    x = _bin(x, radius)
    # Record number of proximities until they are equal to the diagonal length
    xe = embed(x,lmin,1)
    ppitem = proximities(xe)
    pp = [ppitem]
    while pp[end] > n-m+1
        m += 1
        xe = reembed1(xe,nd)
        ppitem = proximities(xe)
        push!(pp, ppitem)
    end
    # Add proximities one further dimension (0 if m reached n)
    push!(pp, n-m)
    # Make histogram of diagonals, prepended by zeros before lmin
    diag_hist = zeros(eltype(pp), lmin-1)
    # (Formula not in the paper)
    @views append!(diag_hist, pp[1:end-2] .- 2pp[2:end-1] .+ pp[3:end])
    # Add LOI to the histogram if required
    loi && append!(diag_hist, 1)
    entropy(diag_hist; lmin=lmin)
end

"""
    approx_rqa
"""
function approx_rqa_extra(x, radius; lmin=2, loi::Bool=false)
    n = size(x,1)
    nd = size(x,2)    
    m = 1
    x = _bin(x, radius)
    # Record number of proximities until they are equal to the diagonal length 
    ppitem = proximities(x)
    pp = [ppitem]
    while ppitem > n-m+1
        m += 1
        x = reembed1(x,nd)
        ppitem = proximities(x)
        push!(pp, ppitem)
    end
    # Add proximities one further dimension (0 if m reached n)
    push!(pp, n-m)
    @views dhist = pp[1:end-2] .- 2pp[2:end-1] .+ pp[3:end]
    # Add LOI to the histogram if required
    loi && append!(dhist, 1)
    # 
    pp1 = pp[1]
    ppm = pp[lmin]
    ppmx = pp[lmin+1]
    Dict("RR"  => typeof(0.0)( loi ? pp1/n^2 : (pp1-n)/n^2 ),
        "DET"  => typeof(0.0)( loi ? (lmin*ppm - (lmin-1)*ppmx)/pp1 : (lmin*ppm - (lmin-1)*ppmx - n)/(pp1-n) ),
        "L"    => typeof(0.0)( loi ? (lmin*ppm - (lmin-1)*ppmx)/(ppm-ppmx) : (lmin*ppm - (lmin-1)*ppmx - n)/(ppm-ppmx-1) ),
        "Lmax" => typeof(0.0)( loi ? n : length(dhist) ),
        "DIV"  => typeof(0.0)( loi ? 1/n : 1/length(dhist) ),
        "ENT"  => entropy(dhist; lmin=lmin)
    )
end


