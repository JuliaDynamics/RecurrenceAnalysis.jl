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

# Initial routine: bin x according to the radius and set loi according to kwargs
function _approx_common(x, radius, loi; kwargs...)
    kwargs = Dict(kwargs)
    haskey(kwargs, :scale) && (radius *= scale)
    if haskey(kwargs, :theiler)
        if kwargs[:theiler] == 0
            loi = loi && true
        elseif kwargs[:theiler] == 1
            loi = loi && false
        else
            @warn "The Theiler window can only be 0 or 1 for the approximated recurrence rate."
        end
        @info "Use `loi` (`true` or `false`) to declare if the LOI must be included."
        @info "The LOI will $(loi ? "" : "not ")be evaluated."
    end
    if radius > 0
        delta = 2radius
        minx = minimum(x)
        typeint = typeof(1)
        xb = typeint.(div.(x.-minx, delta)) .+ 1
    end
    (xb, loi)
end

"""
    approx_recurrencerate
"""
function approx_recurrencerate(x, radius; loi::Bool=true, kwargs...)
    (x, loi) = _approx_common(x, radius, loi; kwargs...)
    pp = proximities(x)
    n = size(x,1)
    loi ? pp/n^2 : (pp-n)/n^2
end

"""
    approx_determinism
"""
function approx_determinism(x, radius; lmin=2, loi::Bool=true, kwargs...)
    (x, loi) = _approx_common(x, radius, loi; kwargs...)
    pp1 = proximities(x)
    xe = embed(x,lmin,1)
    ppm = proximities(xe)
    xex = reembed1(xe,size(x,2))
    ppmx = proximities(xex)
    n = size(x,1)
    loi ? (lmin*ppm - (lmin-1)*ppmx)/pp1 : (lmin*ppm - (lmin-1)*ppmx - n)/(pp1-n)
end

function approx_avgdiag(x, radius; lmin=2, loi::Bool=true, kwargs...)
    (x, loi) = _approx_common(x, radius, loi; kwargs...)
    pp1 = proximities(x)
    xe = embed(x,lmin,1)
    ppm = proximities(xe)
    xex = reembed1(xe,size(x,2))
    ppmx = proximities(xex)
    n = size(x,1)
    loi ? (lmin*ppm - (lmin-1)*ppmx)/(ppm-ppmx) : (lmin*ppm - (lmin-1)*ppmx - n)/(ppm-ppmx-1)
end

"""
    approx_maxdiag
"""
function approx_maxdiag(x, radius; lmin=2, loi::Bool=true, kwargs...)
    n = size(x,1)
    nd = size(x,2)
    # The longest diagonal is the LOI (if it is not excluded)
    loi && (return n)
    (x, loi) = _approx_common(x, radius, loi; kwargs...)
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
function approx_entropy(x, radius; lmin=2, loi=true, kwargs...)
    n = size(x,1)
    nd = size(x,2)    
    m = lmin
    (x, loi) = _approx_common(x, radius, loi; kwargs...)
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
    entropy(diag_hist; lmin=lmin, kwargs...)
end

"""
    approx_rqa
"""
function approx_rqa(x, radius; lmin=2, loi::Bool=true, kwargs...)
    # Manage custom arguments for diagonals
    kwargs = Dict(kwargs)
    haskey(kwargs, :lmindiag) && (lmin = kwargs[:lmindiag])
    haskey(kwargs, :theilerdiag) && (kwargs[:theiler] = kwargs[:theilerdiag])
    n = size(x,1)
    nd = size(x,2)    
    m = 1
    (x, loi) = _approx_common(x, radius, loi; kwargs...)    
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
        "ENT"  => entropy(dhist; kwargs...)
    )
end


