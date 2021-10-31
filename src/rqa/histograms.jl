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
function _linehistograms(rows::T, cols::T, lmin::Integer, theiler::Integer,
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
    while (firstindex ≤ n) && (abs(rows[firstindex]-cols[firstindex]) < theiler)
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
        if abs(r-c) ≥ theiler
            # Search the second and later segments in the column
            if c == cprev
                if r-rprev != 1 # (a): there is a separation between rprev and r
                    # update histogram of segments
                    current_vert = rprev-r1+1
                    if current_vert ≥ lmin
                        nbins = extendhistogram!(bins, nbins, current_vert)
                    end
                    # update histogram of distances if it there were at least
                    # two previous segments in the column
                    if distances
                        halfline = div(current_vert, 2)
                        if dist != 0
                            nbins_d = extendhistogram!(bins_d, nbins_d, dist+halfline)
                        end
                        # update the distance
                        dist = r-rprev+halfline-1
                    end
                    r1 = r # update the start of the next segment
                end
                rprev = r  # update the previous position
            else # hit in the first point of a new column
                # process the last fragment of the previous column
                current_vert = rprev-r1+1
                if current_vert ≥ lmin
                    nbins = extendhistogram!(bins, nbins, current_vert)
                end
                if  distances
                    if dist!= 0
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
    if current_vert ≥ lmin
        nbins = extendhistogram!(bins, nbins, current_vert)
    end
    if  distances && dist != 0
        halfline = div(current_vert, 2)
        nbins_d = extendhistogram!(bins_d, nbins_d, dist+halfline)
    end
    return (bins, bins_d)
end

deftheiler(x::Union{RecurrenceMatrix,JointRecurrenceMatrix}) = 1
deftheiler(x::Union{CrossRecurrenceMatrix,SparseMatrixCSC}) = 0

function diagonalhistogram(x::Union{ARM,SparseMatrixCSC}; lmin::Integer=2, theiler::Integer=deftheiler(x), kwargs...)
    (theiler < 0) && throw(ErrorException(
        "Theiler window length must be greater than or equal to 0"))
    (lmin < 1) && throw(ErrorException("lmin must be 1 or greater"))
    m,n = size(x)
    rv = rowvals(x)[:]
    dv = colvals(x) .- rowvals(x)
    loi_hist = Int[]
    if issymmetric(x)
        # If theiler==0, the LOI is counted separately to avoid duplication
        if theiler == 0
            if all(isequal(1),x.nzval)==false #CrossRecurrenceMatrix won't work if all values aren't equal to 1 and theiler=0
                xtmp = deepcopy(x)
                xtmp.nzval .= 1
                loi_hist = verticalhistograms(CrossRecurrenceMatrix(hcat(diag(xtmp,0))),lmin=lmin)[1]
            else
                loi_hist = verticalhistograms(CrossRecurrenceMatrix(hcat(diag(x,0))),lmin=lmin)[1]
            end
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


function verticalhistograms(x::Union{ARM,SparseMatrixCSC};
    lmin::Integer=2, theiler::Integer=deftheiler(x), distances=true, kwargs...)
    (theiler < 0) && throw(ErrorException(
        "Theiler window length must be greater than or equal to 0"))
    (lmin < 1) && throw(ErrorException("lmin must be 1 or greater"))
    m,n=size(x)
    rv = rowvals(x)
    cv = colvals(x)
    return _linehistograms(rv,cv,lmin,theiler,distances)
end

vl_histogram(x; kwargs...) = verticalhistograms(x; kwargs...)[1]
rt_histogram(x; kwargs...) = verticalhistograms(x; kwargs...)[2]

"""
    recurrencestructures(x::AbstractRecurrenceMatrix;
                             diagonal=true,
                             vertical=true,
                             recurrencetimes=true,
                             kwargs...)

Return a dictionary with the
histograms of the recurrence structures contained in the recurrence matrix `x`,
with the keys `"diagonal"`, `"vertical"` or
`"recurrencetimes"`, depending on what keyword arguments are given as `true`.

## Description
Each item of the dictionary is a vector of integers, such that the `i`-th
element of the vector is the number of lines of length `i` contained in `x`.

* `"diagonal"` counts the diagonal lines, i.e. the recurrent trajectories.
* `"vertical"` counts the vertical lines, i.e. the laminar states.
* `"recurrencetimes"` counts the vertical distances between recurrent states,
    i.e. the recurrence times.

All the points of the matrix are counted by default. The keyword argument
`theiler` can be passed to rule out the lines around the main
diagonal. See the arguments of the function [`rqa`](@ref) for further details.

"Empty" histograms are represented always as `[0]`.

*Notice*: There is not a unique operational definition of "recurrence times". In the
analysis of recurrence plots, usually the  "second type" of recurrence times as
defined by Gao and Cai [1] are considered, i.e. the distance between
consecutive (but separated) recurrent structures in the vertical direction of
the matrix. But that distance is not uniquely defined when the vertical recurrent
structures are longer than one point. The recurrence times calculated here are
the distance between the midpoints of consecutive lines, which is a balanced
estimator of the Poincaré recurrence times [2].

## References

[1] J. Gao & H. Cai. "On the structures and quantification of recurrence plots".
[*Physics Letters A*, 270(1-2), 75–87 (2000)](https://www.sciencedirect.com/science/article/pii/S0375960100003042?via%3Dihub).

[2] N. Marwan & C.L. Webber, "Mathematical and computational foundations of
recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence
Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).
"""
function recurrencestructures(x::Union{ARM,SparseMatrixCSC};
    diagonal=true, vertical=true, recurrencetimes=true,
    theiler=deftheiler(x), kwargs...)

    # Parse arguments for diagonal and vertical structures
    histograms = Dict{String,Vector{Int}}()
    if diagonal
        kw_d = Dict(kwargs)
        haskey(kw_d, :theilerdiag) && (kw_d[:theiler] = kw_d[:theilerdiag])
        histograms["diagonal"] = diagonalhistogram(x; theiler=theiler, kw_d..., lmin=1)
    end
    if vertical || recurrencetimes
        kw_v = Dict(kwargs)
        haskey(kw_v, :theilervert) && (kw_v[:theiler] = kw_v[:theilervert])
        vhist = verticalhistograms(x; theiler=theiler, kw_v..., lmin=1)
        vertical && (histograms["vertical"] = vhist[1])
        recurrencetimes && (histograms["recurrencetimes"] = vhist[2])
    end
    return histograms
end
