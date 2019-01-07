# Recurrence parameters as defined by Marwan et al. (2007)

# Recurrence rate

"""
    recurrencerate(x; theiler=0)

Calculate the recurrence rate (RR) of the recurrence matrix `x`, ruling out
the points within the Theiler window of size `theiler`.
"""
function recurrencerate(x::ARM; theiler::Integer=0, kwargs...)::Float64
    (theiler < 0) && throw(ErrorException(
        "Theiler window length must be greater than or equal to 0"))
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

# 1. Based on diagonal lines

"""
    determinism(x::AbstractRecurrenceMatrix; lmin=2, theiler=0)

Calculate the determinism (DET) of the recurrence matrix `x`, ruling out
the points within the Theiler window of size `theiler` and diagonals shorter
than `lmin`.


    determinism(h::Vector{<:Integer}, npoints)

Calculate the determinism from the histogram of diagonal lines `h`,
e.g. obtained from [`recurrencestructures`](@ref), and the total number of points
(`npoints`).
"""
function determinism(x::ARM; kwargs...)
    npoints = recurrencerate(x; kwargs...)*length(x)
    return determinism(diagonalhistogram(x; kwargs...), npoints)
end

function determinism(diag_hist::Vector{<:Integer}, npoints)::Float64
    diag_points = (1:length(diag_hist)) .* diag_hist
    return sum(diag_points)/npoints
end

"""
    avgdiag(x; lmin=2, theiler=0)

Calculate the average diagonal length (L) in the recurrence matrix `x`, ruling out
the points within the Theiler window of size `theiler` and diagonals shorter
than `lmin`.

    avgdiag(h::Vector{<:Integer})

Calculate the average diagonal length from the histogram of diagonal lines `h`,
e.g. obtained from [`recurrencestructures`](@ref).
"""
avgdiag(x::ARM; kwargs...) = avgdiag(diagonalhistogram(x; kwargs...))

function avgdiag(diag_hist::Vector{<:Integer})::Float64
    diag_points = collect(1:length(diag_hist)) .* diag_hist
    return sum(diag_points)/sum(diag_hist)
end

"""
    maxdiag(x; theiler=0)

Calculate the longest diagonal (Lmax) in the recurrence matrix `x`, ruling out
the points within the Theiler window of size `theiler`.

    maxdiag(h::Vector{<:Integer})

Calculate the longest diagonal from the histogram of diagonal lines `h`,
e.g. obtained from [`recurrencestructures`](@ref).
"""
maxdiag(x::ARM, kwargs...) = maxdiag(diagonalhistogram(x; kwargs...))
    
maxdiag(diag_hist::Vector{<:Integer}) = length(diag_hist)

"""
    divergence(x; theiler=0)

Calculate the divergence of the recurrence matrix `x`
(actually the inverse of [`maxdiag`](@ref)), ruling out
the points within the Theiler window of size `theiler`.

    divergence(h::Vector{<:Integer})

Calculate the divergence from the histogram of diagonal  lines `h`,
e.g. obtained from [`recurrencestructures`](@ref).
"""
divergence(x; kwargs...) = Float64( 1/maxdiag(x; kwargs...) )

"""
    rqaentropy(x; lmin=2, theiler=0)

Calculate the Shannon entropy of diagonal lengths (ENT) of
the recurrence matrix `x`, ruling out the points within the Theiler window of
size `theiler`.

    rqaentropy(h::Vector{<:Integer})

Calculate the RQA entropy from the histogram of diagonal lines `h`,
e.g. obtained from [`recurrencestructures`](@ref).
"""
rqaentropy(x::ARM, kwargs...) = rqaentropy(diagonalhistogram(x; kwargs...))
    
function rqaentropy(diag_hist::Vector)::Float64
    prob_bins = diag_hist ./ sum(diag_hist)
    prob_bins = prob_bins[findall(!iszero, prob_bins)]
    return -sum(prob_bins .* log.(prob_bins))
end

"""
    trend(x; theiler=0, border=10)

Calculate the trend of recurrences in the recurrence matrix `x`
towards its edges, ruling out the points within the Theiler window of size `theiler`
and in the outermost diagonals (at `border` from the corners of the matrix).
"""
trend(x::ARM; kwargs...) = _trend(tau_recurrence(x); kwargs...)

function tau_recurrence(x::ARM)
    n = minimum(size(x))
    [count(!iszero, diag(x,d))/(n-d) for d in (0:n-1)]
end

function _trend(npoints::Vector; theiler=0, border=10, kwargs...)::Float64
    nmax = length(npoints)
    rrk = npoints./collect(nmax:-1:1)
    a = 1+theiler
    b = nmax-border
    w = collect(a:b) .- b/2
    (w'*(rrk[a:b] .- mean(rrk[a:b])) ./ (w'*w))[1]
end

# Number of l-length sequences, based on diagonals
countsequences(x::ARM; kwargs...) = countsequences(diagonalhistogram(x; kwargs...))

function countsequences(diag_hist::Vector; lmin=2, kwargs...)
    overlap = (1:length(diag_hist))' .- (lmin+1)
    overlap[overlap .< 0] = 0
    overlap * diag_hist
end


# 2. Based on vertical lines

"""
    laminarity(x; lmin=2, theiler=0)

Calculate the laminarity (LAM) of the recurrence matrix `x`, ruling out the
points within the Theiler window of size `theiler` and diagonals shorter
than `lmin`.

    laminarity(h::Vector{<:Integer}, npoints)

Calculate the laminarity from the histogram of diagonal lines `h`,
e.g. obtained from [`recurrencestructures`](@ref), and the total number of points
(`npoints`).
"""
function laminarity(x::ARM; kwargs...)
    npoints = recurrencerate(x)*length(x)
    return laminarity(verticalhistograms(x; kwargs...)[1], npoints)
end

laminarity(vert_hist::Vector{<:Integer}, npoints) = determinism(vert_hist, npoints)

"""
    trappingtime(x; lmin=2, theiler=0)

Calculate the trapping time (TT) of the recurrence matrix `x`, ruling out the
points within the Theiler window of size `theiler` and diagonals shorter
than `lmin`.

    trappingtime(h::Vector{<:Integer})

Calculate the trapping time from the histogram of diagonal lines `h`,
e.g. obtained from [`recurrencestructures`](@ref).
"""
trappingtime(x::ARM; kwargs...) = trappingtime(verticalhistograms(x; kwargs...)[1])
trappingtime(vert_hist::Vector{<:Integer}) = avgdiag(vert_hist)

"""
    maxvert(x; theiler=0)

Calculate the longest vertical line (Vmax) of the recurrence matrix `x`,
ruling out the points within the Theiler window of size `theiler`.

    maxvert(h::Vector{<:Integer})

Calculate the trapping time from the histogram of diagonal lines `h`,
e.g. obtained from [`recurrencestructures`](@ref).

"""
maxvert(x::ARM; kwargs...) = maxvert(verticalhistograms(x; kwargs...)[1])
maxvert(vert_hist::Vector{<:Integer}) = length(vert_hist)


# 3. Based on recurrence times

"""
    meanrecurrencetime(x[; theiler=0])

Calculate the mean recurrence time (MRT) of a recurrence matrix.
"""
meanrecurrencetime(vert_hist::Vector) = avgdiag(vert_hist)

"""
    recurrencetimeentropy(x[; theiler=0])

Calculate the entropy of the distribution of recurrence times (RTE)
of a recurrence matrix.
"""
recurrencetimeentropy(vert_hist::Vector) = rqaentropy(vert_hist)

"""
    maxnumberrecurrencetime(x[; theiler=0])
    
"""
maxnumberrecurrencetime(verth_hist::Vector) = maximum(verthist)

# Functions for AbstractRecurrenceMatrix

for fun in [:meanrecurrencetime, :recurrencetimeentropy, :maxnumberrecurrencetime]
    @eval begin $(fun)(x::ARM; lmin=2, kwargs...) =
        $(fun)(verticalhistograms(x; lmin=lmin, kwargs...)[2])
    end
end


# 4. All in one

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
    rr_d = recurrencerate(x; kw_d...)
    if onlydiagonal
        return Dict("RR"  => recurrencerate(x; kwargs...),
        "DET"  => determinism(dhist, rr_d*length(x)),
        "L"    => avgdiag(dhist),
        "Lmax" => maxdiag(dhist),
        "DIV"  => divergence(dhist),
        "ENT"  => rqaentropy(dhist)
        )
   else
        kw_v = Dict(kwargs)
        haskey(kw_v, :theilervert) && (kw_v[:theiler] = kw_v[:theilervert])
        haskey(kw_v, :lminvert) && (kw_v[:lmin] = kw_v[:lminvert])
        vhist = verticalhistograms(x; kw_v...)[1]
        rr_v = recurrencerate(x; kw_v...)
        return Dict("RR"  => recurrencerate(x; kwargs...),
            "DET"  => determinism(dhist, rr_d*length(x)),
            "L"    => avgdiag(dhist),
            "Lmax" => maxdiag(dhist),
            "DIV"  => divergence(dhist),
            "ENT"  => rqaentropy(dhist),
            "TND"  => trend(x; kw_d...),
            "LAM"  => laminarity(vhist, rr_v*length(x)),
            "TT"   => trappingtime(vhist),
            "Vmax" => maxvert(vhist)
        )
    end
end

### Histograms of recurrence structures ###

## (to delete: old methods)

function verticalhistograms_old(x::ARM; theiler::Integer=0,  distances=true, kwargs...)
    m,n=size(x)
    # histogram for "black lines"
    bins = [0]
    nbins = 1
    # histogram for "white lines"
    bins_d = [0]
    nbins_d = 1
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
            rprev = r1
            for r in rvc[2:end]
                # Look for nonzero that starts a new column fragment
                # (more than one row after the previous one)
                if r-rprev != 1
                    # white line
                     distances && (nbins_d = extendhistogram!(bins_d, nbins_d, r-rprev-1))
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
        elseif nc==1
            bins[1] += 1
        end
    end
    return (bins, bins_d)
end

function diagonalhistogram_old(x::ARM; theiler::Integer=0, kwargs...)
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
            loi_hist = verticalhistograms_old(CrossRecurrenceMatrix(hcat(diag(x,0))))[1]
        end
    else
        valid = (abs.(dv) .>= theiler)
        f = 1
    end
    vmat = CrossRecurrenceMatrix(sparse(rv[valid], dv[valid] .+ (m+1), true))
    dh = f .* verticalhistograms_old(vmat; theiler=0)[1]
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
## (end: to delete) ##


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

# Calculate the histograms of segments and distances between segments
# from the indices of rows and columns/diagonals of the matrix
# `theiler` is used for histograms of vertical structures
# `distances` is used to simplify calculations of the distances are not wanted
function _linehistograms(rows::T, cols::T, lmin::Integer=1, theiler::Integer=0, 
     distances::Bool=true) where {T<:AbstractVector{Int}}
    
    # check bounds
    n = length(rows)
    if length(cols) != n
        throw(ErrorException("mismatch between number of row and column indices"))
    end
    # histogram for segments
    bins = [0]
    nbins = 1
    # histogram for distances between segments
    bins_d = [0]
    nbins_d = 1
    # find the first point outside the Theiler window
    firstindex = 1
    while (firstindex<=n) && (abs(rows[firstindex]-cols[firstindex])<theiler)
        firstindex += 1
    end
    if firstindex > n
        return ([0],[0])
    end
    # Iterate over columns
    cprev = cols[firstindex]
    r1 = rows[firstindex]
    rprev = r1 - 1 # coerce that (a) is not hit in the first iteration
    dist = 0 # placeholder for distances between segments
    @inbounds for i=firstindex:n
        r = rows[i]
        c = cols[i]
        if abs(r-c)>=theiler
            # Search the second and later segments in the column
            if c == cprev
                if r-rprev !=1 # (a): there is a separation between rprev and r
                    # update histogram of segments
                    current_vert = rprev-r1+1
                    if current_vert >= lmin
                        nbins = extendhistogram!(bins, nbins, current_vert)
                    end
                    # update histogram of distances if it there were at least
                    # two previous segments in the column
                    if distances
                        if current_vert >= lmin
                            halfline = div(current_vert, 2)
                                if dist != 0
                                    nbins_d = extendhistogram!(bins_d, nbins_d, dist+halfline)
                                end                        
                            # update the distance
                            dist = r-rprev+halfline-1
                        else
                            dist += r - r1
                        end
                    end
                    r1 = r # update the start of the next segment
                end
                rprev = r  # update the previous position
            else # hit in the first point of a new column
                # process the last fragment of the previous column
                current_vert = rprev-r1+1
                if current_vert >= lmin
                    nbins = extendhistogram!(bins, nbins, current_vert)
                end
                if  distances
                    if dist!= 0 && current_vert >= lmin
                        halfline = div(current_vert, 2)
                        nbins_d = extendhistogram!(bins_d, nbins_d, dist+halfline)
                    end
                    dist = 0
                end
                # initialize values for searching new fragments
                cprev = c
                r1 = r
                rprev = r
            end
        end
    end
    # Process the latest fragment
    current_vert = rprev-r1+1
    if current_vert >= lmin
        nbins = extendhistogram!(bins, nbins, current_vert)
        if  distances && dist != 0
            halfline = div(current_vert, 2)
            nbins_d = extendhistogram!(bins_d, nbins_d, dist+halfline)
        end
    end
    return (bins, bins_d)
end

function diagonalhistogram(x::ARM; lmin::Integer=2, theiler::Integer=0, kwargs...)
    (theiler < 0) && throw(ErrorException(
        "Theiler window length must be greater than or equal to 0"))
    (lmin < 1) && throw(ErrorException("lmin must be 1 or greater"))
    m,n=size(x)
    rv = rowvals(x)[:]
    dv = colvals(x) .- rowvals(x)
    loi_hist = Int[]
    if issymmetric(x)
        # If theiler==0, the LOI is counted separately to avoid duplication
        if theiler == 0
            loi_hist = verticalhistograms(CrossRecurrenceMatrix(hcat(diag(x,0))),lmin=lmin)[1]
        end
        inside = (dv .< max(theiler,1))
        f = 2
    else
        inside = (abs.(dv) .< theiler)
        f = 1
    end
    # remove Theiler window and short by diagonals
    deleteat!(rv, inside)
    deleteat!(dv, inside)
    rv = rv[sortperm(dv)]
    dv = sort!(dv)
    dh = f.*_linehistograms(rv, dv, lmin, 0, false)[1]
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
    return dh
end

function verticalhistograms(x::ARM;
    lmin::Integer=2, theiler::Integer=0, distances=true, kwargs...)
    (theiler < 0) && throw(ErrorException(
        "Theiler window length must be greater than or equal to 0"))
    (lmin < 1) && throw(ErrorException("lmin must be 1 or greater"))
    m,n=size(x)
    rv = rowvals(x)
    cv = colvals(x)
    return _linehistograms(rv,cv,lmin,theiler,distances)
end

"""
    recurrencestructures(x::AbstractRecurrenceMatrix;
                             diagonal=true,
                             vertical=true,
                             recurrencetimes=true,
                             kwargs...)
    
Histograms of the recurrence structures contained in the recurrence matrix `x`.

## Description:

Returns a dictionary with the keys `"diagonal"`, `"vertical"` or
`"recurrencetimes"`, depending on what keyword arguments are given as `true`.
Each item of the dictionary is a vector of integers, such that the `i`-th
element of the vector is the number of lines of length `i` contained in `x`.

* `"diagonal"` counts the diagonal lines, i.e. the recurrent trajectories.
* `"vertical"` counts the vertical lines, i.e. the laminar states.
* `"recurrencetimes"` counts the vertical distances between recurrent states,
    i.e. the recurrence times. These are calculated as the distance between
    the middle points of consecutive vertical lines.

All the points of the matrix are counted by default. Extra keyword arguments can
be passed to rule out the lines shorter than a minimum length or around the main
diagonal. See the arguments of the function [`rqa`](@ref) for further details.

## References:

N. Marwan & C.L. Webber, "Mathematical and computational foundations of
recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence
Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).
"""
function recurrencestructures(x::ARM;
    diagonal=true, vertical=true, recurrencetimes=true, lmin=1, theiler=0, kwargs...)
    
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
        vhist = verticalhistograms(x; kw_v...)
        vertical && (histograms["vertical"] = vhist[1])
        recurrencetimes && (histograms["recurrencetimes"] = vhist[2])
    end
    return histograms
end

