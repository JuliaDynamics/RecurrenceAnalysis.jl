### Histograms of recurrence structures

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
    (dh==[0]) && @warn "no diagonal lines found in the matrix"
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

function vl_histogram(x; kwargs...) 
    h = verticalhistograms(x; kwargs...)[1]
    (h==[0]) && @warn "no vertical lines found in the matrix"
    return h
end
function rt_histogram(x; kwargs...)
    h = verticalhistograms(x; kwargs...)[2]
    (h==[0]) && @warn "no recurrence times found in the matrix"
    return h
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
        histograms["diagonal"] = diagonalhistogram(x; lmin=lmin, theiler=theiler, kw_d...)
    end
    if vertical || recurrencetimes
        kw_v = Dict(kwargs)
        haskey(kw_v, :theilervert) && (kw_v[:theiler] = kw_v[:theilervert])
        haskey(kw_v, :lminvert) && (kw_v[:lmin] = kw_v[:lminvert])
        vhist = verticalhistograms(x; lmin=lmin, theiler=theiler, kw_v...)
        if vertical
            (vhist[1]==[0]) && @warn "no vertical lines found in the matrix"
            (histograms["vertical"] = vhist[1])
        end
        if recurrencetimes
            (vhist[2]==[0]) && @warn "no recurrence times found in the matrix"
            (histograms["recurrencetimes"] = vhist[2])
        end
    end
    return histograms
end

