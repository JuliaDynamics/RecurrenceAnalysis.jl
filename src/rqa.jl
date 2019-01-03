# Recurrence parameters as defined by Marwan et al. (2007)

"""
    recurrencerate(x; theiler=0)

Calculate the recurrence rate (RR) of a recurrence matrix, ruling out
the points within the Theiler window.
"""
function recurrencerate(x::ARM; theiler::Integer=0, kwargs...)::Float64
    theiler < 0 && error("Theiler window length must be greater than or equal to 0")
    if theiler == 0
        return nnz(x)/length(x)
    end
    diags_remove = -(theiler-1):(theiler-1)
    theiler_nz = 0
    for d in diags_remove
        theiler_nz += nnz(diag(x,d))
    end
    return (nnz(x) - theiler_nz)/length(x)
end

function tau_recurrence(x::ARM)
    n = minimum(size(x))
    [count(!iszero, diag(x,d))/(n-d) for d in (0:n-1)]
end

function recurrencestructures(x::ARM; diagonal=true, vertical=true, recurrencetimes=true, kwargs...)
    
    # Parse arguments for diagonal and vertical structures
    histograms = Dict{String,Vector{Int}}()
    if diagonal
        kw_d = Dict(kwargs)
        haskey(kw_d, :theilerdiag) && (kw_d[:theiler] = kw_d[:theilerdiag])
        haskey(kw_d, :lmindiag) && (kw_d[:lmin] = kw_d[:lmindiag])
        histograms["diagonal"] = diagonalhistogram(x; kw_d...)
    end
    if vertical || recurrencetimes
        kw_v = Dict(kwargs)
        haskey(kw_v, :theilervert) && (kw_v[:theiler] = kw_v[:theilervert])
        haskey(kw_v, :lminvert) && (kw_v[:lmin] = kw_v[:lminvert])
        vhist = verticalhistogram(x; kw_v...)
        vertical && (histograms["vertical"] = vhist[1])
        recurrencetimes && (histograms["recurrencetimes"] = vhist[2])
    end
    return histograms
end

# Based on diagonal lines

function diagonalhistogram(x::ARM; theiler::Integer=0, kwargs...)
    theiler < 0 && error("Theiler window length must be greater than or equal to 0")
    m,n=size(x)
    rv = rowvals(x)
    dv = colvals(x) .- rowvals(x)
    loi_hist = Int[]
    if issymmetric(x)
        valid = (dv .>= max(theiler,1))
        f = 2
        # If theiler==0, the LOI is counted separately to avoid duplication
        if theiler == 0
            loi_hist = verticalhistogram(CrossRecurrenceMatrix(hcat(diag(x,0))))[1]
        end
    else
        valid = (abs.(dv) .>= theiler)
        f = 1
    end
    vmat = CrossRecurrenceMatrix(sparse(rv[valid], dv[valid] .+ (m+1), true))
    dh = f .* verticalhistogram(vmat; theiler=0)[1]
    # Add frequencies of LOI if suitable
    if (nbins_loi = length(loi_hist)) > 0
        nbins = length(dh)
        if nbins_loi > nbins
            loi_hist[1:nbins] .+= dh
            dh = loi_hist
        else
            dh[1:nbins_loi] .+= loi_hist
        end
    end
    dh
end

"""
    determinism(x; lmin=2, theiler=0)

Calculate the determinism (DET) of a recurrence matrix, ruling out
the points within the Theiler window and diagonals shorter than a minimum value.
"""
function determinism(diag_hist::Vector; lmin=2, kwargs...)::Float64
    if lmin < 2
        error("lmin must be 2 or greater")
    end
    nbins = length(diag_hist)
    diag_points = collect(1:nbins) .* diag_hist
    return (lmin > nbins) ? 0.0 : sum(diag_points[lmin:nbins])/sum(diag_points)
end

determinism(x::ARM; kwargs...) = determinism(diagonalhistogram(x; kwargs...); kwargs...)

"""
    avgdiag(x; lmin=2, theiler=0)

Calculate the average diagonal length (L) in a recurrence matrix, ruling out
the points within the Theiler window and diagonals shorter than a minimum value.
"""
function avgdiag(diag_hist::Vector; lmin=2, kwargs...)::Float64
    if lmin < 2
        error("lmin must be 2 or greater")
    end
    nbins = length(diag_hist)
    diag_points = collect(1:nbins) .* diag_hist
    return (lmin > nbins) ? 0.0 : sum(diag_points[lmin:nbins])/sum(diag_hist[lmin:nbins])
end

avgdiag(x::ARM; kwargs...) =  avgdiag(diagonalhistogram(x; kwargs...); kwargs...)

"""
    maxdiag(x; theiler=0)

Calculate the longest diagonal (Lmax) in a recurrence matrix, ruling out
the points within the Theiler window.
"""
maxdiag(diag_hist::Vector; kwargs...) = length(diag_hist)
maxdiag(x::ARM; kwargs...) = maxdiag(diagonalhistogram(x; kwargs...))

"""
    divergence(x; lmin=2, theiler=0)

Calculate the divergence of a recurrence matrix
(actually the inverse of `maxdiag`).
"""
divergence(x; kwargs...) = Float64( 1/maxdiag(x; kwargs...) )

"""
    rqaentropy(x; lmin=2, theiler=0)

Calculate the Shannon entropy of diagonal lengths (ENT) of a recurrence matrix, ruling out
the points within the Theiler window and diagonals shorter than a minimum value.
"""
function rqaentropy(diag_hist::Vector; lmin=2, kwargs...)::Float64
    if lmin < 2
        error("lmin must be 2 or greater")
    end
    nbins = length(diag_hist)
    if lmin <= nbins
        prob_bins = diag_hist[lmin:nbins] ./ sum(diag_hist[lmin:nbins])
        prob_bins = prob_bins[findall(!iszero, prob_bins)]
        return -sum(prob_bins .* log.(prob_bins))
    else
        return 0.0
    end
end

rqaentropy(x::ARM; kwargs...) = rqaentropy(diagonalhistogram(x; kwargs...); kwargs...)

"""
    trend(x; theiler=0, border=10)

Calculate the trend of recurrences in recurrence matrix towards its edges, ruling out
the points within the Theiler window and in the outermost diagonals.
"""
function trend(npoints::Vector; theiler=0, border=10, kwargs...)::Float64
    nmax = length(npoints)
    rrk = npoints./collect(nmax:-1:1)
    a = 1+theiler
    b = nmax-border
    w = collect(a:b) .- b/2
    (w'*(rrk[a:b] .- mean(rrk[a:b])) ./ (w'*w))[1]
end

trend(x::ARM; kwargs...) = trend(tau_recurrence(x); kwargs...)

# Number of l-length sequences, based on diagonals
function countsequences(diag_hist::Vector; lmin=2, kwargs...)
    overlap = (1:length(diag_hist))' .- (lmin+1)
    overlap[overlap .< 0] = 0
    overlap * diag_hist
end

countsequences(x::ARM; kwargs...) = countsequences(diagonalhistogram(x; kwargs...); kwargs...)

# 2. Based on vertical lines

# Add one item to position `p` in the histogram `h` that has precalculated length `n`
# - update the histogram and return its new length
@inline function extendhistogram!(h::Vector{Int}, n::Int, p::Int)
    if p > n
        append!(h, zeros(p-n))
        n = p
    end
    h[p] += 1
    return n
end

# macro extend ex # the expression must be of the type h[p] += 1

function verticalhistogram(x::ARM; theiler::Integer=0, whitelines=true, kwargs...)
    m,n=size(x)
    # histogram for "black lines"
    bins = [0]
    nbins = 1
    # histogram for "white lines"
    bins_w = [0]
    nbins_w = 1
    rv = rowvals(x)
    # Iterate over columns
    for c = 1:n
        rvc = rv[nzrange(x,c)]
        # Remove theiler window if needed
        if theiler != 0
            rvc = rvc[(rvc .<= c-theiler) .| (rvc .>= c+theiler)]
        end
        nc = length(rvc)
        if nc>1
            r1 = rvc[1]
            # add white line prior to the first point if it exists
            whitelines && (r1 != 1) && (nbins_w = extendhistogram!(bins_w, nbins_w, r1-1))
            rprev = r1
            for r in rvc[2:end]
                # Look for nonzero that starts a new column fragment
                # (more than one row after the previous one)
                if r-rprev != 1
                    # white line
                    whitelines && (nbins_w = extendhistogram!(bins_w, nbins_w, r-rprev-1))
                    # black line
                    current_vert = rprev-r1+1
                    nbins = extendhistogram!(bins, nbins, current_vert)
                    r1 = r
                end
                rprev = r
            end
            # Last column fragment
            if rprev-rvc[end-1] == 1
                current_vert = rprev-r1+1
                nbins = extendhistogram!(bins, nbins, current_vert)
            else
                bins[1] += 1
            end
            # add white lines after the last point if it exists
            whitelines && (rprev != m) && (nbins_w = extendhistogram!(bins_w, nbins_w, m-rprev))
        elseif nc==1
            # add a single black point, and white lines around it
            bins[1] += 1
            whitelines && begin
               (rvc[1] != 1) && (nbins_w = extendhistogram!(bins_w, nbins_w, rvc[1]-1))
             (rvc[1] != m) && (nbins_w = extendhistogram!(bins_w, nbins_w, m-rvc[1]))
            end
        else # add full column as white line
            whitelines && (nbins_w = extendhistogram!(bins_w, nbins_w, m))
        end
    end
    return (bins, bins_w)
end

"""
    laminarity(x; lmin=2, theiler=0)

Calculate the laminarity (LAM) of a recurrence matrix, ruling out
vertical lines shorter than a minimum value.
"""
laminarity(vert_hist::Vector; kwargs...) = determinism(vert_hist; kwargs...)
laminarity(x::ARM; kwargs...) =
laminarity(verticalhistogram(x; kwargs...)[1]; kwargs...)

"""
    trappingtime(x; lmin=2, theiler=0)

Calculate the trapping time (TT) of a recurrence matrix, ruling out
vertical lines shorter than a minimum value.
"""
trappingtime(vert_hist::Vector; kwargs...) = avgdiag(vert_hist; kwargs...)
trappingtime(x::ARM; kwargs...) =
trappingtime(verticalhistogram(x; kwargs...)[1]; kwargs...)

"""
    maxvert(x; theiler=0)

Calculate the longest vertical line (Vmax) of a recurrence matrix.
"""
maxvert(vert_hist::Vector; kwargs...) = length(vert_hist)
maxvert(x::ARM; kwargs...) = maxvert(verticalhistogram(x; kwargs...)[1])

"""
    rqa(x; kwargs...)

Calculate all RQA parameters of a recurrence matrix. See the functions
`recurrencerate`, `determinism`, `avgdiag`, `maxdiag`, `divergence`, `rqaentropy`,
`trend`, `laminarity`, `trappingtime` and `maxvert` for the definition of
the different parameters and the default values of the arguments.
Using this function is much more efficient than calling all individual functions
one by one.

The keyword arguments `theilerdiag`, `lmindiag` may be used to declare specific values
that override the values of `theiler` and `lmin` in the calculation of
parameters related to diagonal structures. Likewise, `theilervert` and
`lminvert` can be used for the calculation of parameters related to vertical
structures.

The returned value is a dictionary with the following keys:

* "RR": recurrence rate (see `recurrencerate`)
* "DET": determinsm (see `determinism`)
* "L": average length of diagonal structures (see `avgdiag`)
* "Lmax": maximum length of diagonal structures (see `maxdiag`)
* "DIV": divergence (see `divergence`)
* "ENT": entropy of diagonal structures (see `rqaentropy`)
* "TND": trend of recurrences (see `trend`)
* "LAM": laminarity (see `laminarity`)
* "TT": trapping time (see `trappingtime`)
* "Vmax": maximum length of vertical structures (`see `maxvert`)

The keyword argument `onlydiagonal` (`false` by default) can be set to `true`
in order to restrict the analysis to the recurrence rate and the parameters related
to diagonal structures ("RR", "DET", "L", "Lmax", "DIV" and "ENT").
"""
function rqa(x; onlydiagonal=false, kwargs...)
    # Parse arguments for diagonal and vertical structures
    kw_d = Dict(kwargs)
    haskey(kw_d, :theilerdiag) && (kw_d[:theiler] = kw_d[:theilerdiag])
    haskey(kw_d, :lmindiag) && (kw_d[:lmin] = kw_d[:lmindiag])
    dhist = diagonalhistogram(x; kw_d...)
    if onlydiagonal
        return Dict("RR"  => recurrencerate(x; kwargs...),
        "DET"  => determinism(dhist; kw_d...),
        "L"    => avgdiag(dhist; kw_d...),
        "Lmax" => maxdiag(dhist),
        "DIV"  => divergence(dhist),
        "ENT"  => rqaentropy(dhist; kw_d...)
        )
   else
        kw_v = Dict(kwargs)
        haskey(kw_v, :theilervert) && (kw_v[:theiler] = kw_v[:theilervert])
        haskey(kw_v, :lminvert) && (kw_v[:lmin] = kw_v[:lminvert])
        vhist = verticalhistogram(x; kw_v...)[1]
        return Dict("RR"  => recurrencerate(x; kwargs...),
            "DET"  => determinism(dhist; kw_d...),
            "L"    => avgdiag(dhist; kw_d...),
            "Lmax" => maxdiag(dhist),
            "DIV"  => divergence(dhist),
            "ENT"  => rqaentropy(dhist; kw_d...),
            "TND"  => trend(x; kw_d...),
            "LAM"  => laminarity(vhist; kw_v...),
            "TT"   => trappingtime(vhist; kw_v...),
            "Vmax" => maxvert(vhist)
        )
    end
end
