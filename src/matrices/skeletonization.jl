#=
In this file the core computations for skeletonizing a recurrence matrix
are defined.
=#

# Core function, which gets exported

"""
    skeletonize(R) → R_skel

Skeletonizes the `RecurrenceMatrix R` by using the algorithm proposed by
Kraemer & Marwan [^Kraemer2019]. This function returns `R_skel`, a recurrence
matrix, which only consists of diagonal lines of "thickness" one.

[^Kraemer2019]: Kraemer, K.H., Marwan, N. (2019). [Border effect corrections for diagonal line based recurrence quantification analysis measures. Physics Letters A 383(34)](https://doi.org/10.1016/j.physleta.2019.125977).
"""
function skeletonize(X::Union{ARM,SparseMatrixCSC})
    if issymmetric(X)
        symm = true
        # convert lower triangle into a close returns map
        X_cl1 = create_close_returns_map(tril(X))
        # get "vertical" line distribution with position indices
        lines1t, lines1l, lines1c = verticalhisto(X_cl1)
        lines_copy1t = deepcopy(lines1t)
        lines_copy1l = deepcopy(lines1l)
        lines_copy1c = deepcopy(lines1c)

        # create a close returns map with horizontal lines represented by
        # numbers, equal to their lengths
        X_vertical1 = create_close_returns_map(lines_copy1t, lines_copy1l, lines_copy1c, size(X_cl1))

    else
        symm = false
        # convert upper and lower triangles into close returns maps
        X_cl1 = create_close_returns_map(tril(X))
        X_cl2 = create_close_returns_map(sparse(triu(X)'))

        # get "vertical" line distribution with position indices
        lines1t, lines1l, lines1c = verticalhisto(X_cl1)
        lines_copy1t = deepcopy(lines1t)
        lines_copy1l = deepcopy(lines1l)
        lines_copy1c = deepcopy(lines1c)
        lines2t, lines2l, lines2c = verticalhisto(X_cl2)
        lines_copy2t = deepcopy(lines2t)
        lines_copy2l = deepcopy(lines2l)
        lines_copy2c = deepcopy(lines2c)

        # create close returns maps with horizontal lines represented by
        # numbers, equal to their lengths
        X_vertical1 = create_close_returns_map(lines_copy1t, lines_copy1l, lines_copy1c, size(X_cl1))
        X_vertical2 = create_close_returns_map(lines_copy2t, lines_copy2l, lines_copy2c, size(X_cl2))
    end

    # scan the lines, start with the longest one and discard all adjacent lines
    lines_final_t, lines_final_l, lines_final_c = get_final_line_matrix(lines1t, lines1l, lines1c, lines_copy1t, lines_copy1l, lines_copy1c, X_vertical1)

    # if not symmetric input RP, than compute for the upper triangle as well
    if !symm
        lines_final2_t, lines_final2_l, lines_final2_c = get_final_line_matrix(lines2t, lines2l, lines2c, lines_copy2t, lines_copy2l, lines_copy2c, X_vertical2)
        # build RP based on the histogramm of the reduced lines
        X_new = build_skeletonized_RP(lines_final_t, lines_final_l, lines_final_c, lines_final2_t, lines_final2_l, lines_final2_c, size(X_vertical1,1), size(X_vertical1,2))
    else
        # build RP based on the histogramm of the reduced lines
        X_new = build_skeletonized_RP(lines_final_t, lines_final_l, lines_final_c, size(X_vertical1,1), size(X_vertical1,2))
    end
    return X_new
end

# Auxiliary functions

# Transforms the standard RP into a close returns map
function create_close_returns_map(R::SparseMatrixCSC; triangle::Bool = true)
    nr = size(R, 1)
    nrplus = nr+1
    rowvalues = rowvals(R)
    diagvalues = nrplus .+ colvals(R) .- rowvalues # LOI at nrplus (nr+1)
    if triangle
        lower = (diagvalues .<= nrplus)  # indices of the lower triangle
        return sparse(rowvalues[lower], diagvalues[lower], trues(count(lower)), nr, nrplus)
    else
        return sparse(rowvalues, diagvalues, trues(nnz(R)), nr, 2nr+1)
    end
end
# from the 3 vectors, which store all horizontal lines in the reverted RP,
# (close returns map) [obtained by calling first `create_close_returns_map` on a
# recurrence Matrix R, and consecutevely `verticalhisto(R)`] construct and return
# a "dummy" close returns map, where all lines are encoded according to their length.
function create_close_returns_map(lines1, lines2, lines3, dims)
    nzv = sum(lines1)
    rowvalues = zeros(Int, nzv)
    colvalues = zeros(Int, nzv)
    lenvalues = zeros(Int, nzv)
    startspan = 0
    for (i, val) in enumerate(lines1)
        span = startspan .+ (1:val)
        colvalues[span] .= lines3[i]
        rowvalues[span] .= lines2[i]-1 .+ (1:val)
        lenvalues[span] .= val
        startspan += val
    end
    return sparse(rowvalues, colvalues, lenvalues, dims...)
end

# Transforms the reverted RP (close returns map) into a normal RP
function revert_close_returns_map(R::SparseMatrixCSC; triangle::Bool = true)
    nr = size(R, 1)
    rowvalues = rowvals(R)
    columnvalues = rowvalues .+ colvals(R) .- (nr + 1)
    if triangle
        lower = (columnvalues .<= rowvalues)
        return sparse(rowvalues[lower], columnvalues[lower], trues(count(lower)), nr, nr)
    else
        return sparse(rowvalues, columnvalues, trues(nnz(R)), nr, nr)
    end
end


# deletes a line, specified in 'l_vec' (line vector, with first line being
# the total line length, the second line the line-index of the starting point
# and the third line the column-index of the starting point) from the close
# returns "RP" 'R'.
function delete_line_from_cl_ret_RP!(R::Union{ARM,SparseMatrixCSC}, line1::Int, line2::Int, line3::Int)
    R[line2 .+ (1:line1).-1 , line3] .= 0
end

# build the skeletonized RP from the line-Matrix
function build_skeletonized_RP(lines1::Vector{Int}, lines2::Vector{Int}, lines3::Vector{Int}, N::Int, M::Int)

    X_cl_new = spzeros(Bool, N, M)
    # fill up close returns map with lines stored in the new line matrix
    @inbounds for i = 1:length(lines1)
        l_max = lines1[i]
        linei = lines2[i]
        columni = lines3[i]
        for j = 1:l_max
            X_cl_new[linei+j-1, columni] = true
        end
    end
    XX = revert_close_returns_map(X_cl_new) # revert this close returns map into a legal RP
    X_new = XX .+ XX'
    X_new[diagind(X_new)] .= true # LOI

    return X_new
end
function build_skeletonized_RP(lines1::Vector{Int}, lines2::Vector{Int}, lines3::Vector{Int},
    lines11::Vector{Int}, lines22::Vector{Int}, lines33::Vector{Int}, N::Int, M::Int)

    X_cl_new = spzeros(Bool, N, M)
    X_cl2_new = spzeros(Bool, N, M)
    # fill up close returns map with lines stored in the new line matrix
    @inbounds for i = 1:length(lines1)
        l_max = lines1[i]
        linei = lines2[i]
        columni = lines3[i]
        for j = 1:l_max
            X_cl_new[linei+j-1, columni] = true
        end
    end
    @inbounds for i = 1:length(lines11)
        l_max = lines11[i]
        linei = lines22[i]
        columni = lines33[i]
        for j = 1:l_max
            X_cl2_new[linei+j-1,columni] = true
        end
    end

    XX = revert_close_returns_map(X_cl_new) # revert this close returns map into a legal RP
    XXX= revert_close_returns_map(X_cl2_new)
    X_new = XX .+ XXX'
    X_new[diagind(X_new)] .= true # LOI

    return X_new
end

function get_final_line_matrix(lines1t::Vector{Int}, lines1l::Vector{Int},
    lines1c::Vector{Int}, lines_copy1t::Vector{Int}, lines_copy1l::Vector{Int},
    lines_copy1c::Vector{Int}, X_vertical1::SparseMatrixCSC)

    N, M = size(X_vertical1)
    # initialize final lines
    lines_final_t = Int[]
    lines_final_l = Int[]
    lines_final_c = Int[]

    Nlines1 = length(lines1t) # number of found lines

    # go through all lines stored in the sorted line matrix
    @inbounds for l_ind = 1:Nlines1
        # check if line is still in the rendered line matrix
        common_ind = intersect(findall(x-> x==true, lines1t[l_ind] .== lines_copy1t),
                        findall(x-> x==true, lines1l[l_ind] .== lines_copy1l),
                        findall(x-> x==true, lines1c[l_ind] .== lines_copy1c))
        isempty(common_ind) ? continue : nothing
        # get index pair for start of the line
        linei, columni = lines1l[l_ind], lines1c[l_ind]
        # copy this line in the final line matrix
        push!(lines_final_t, lines1t[l_ind])
        push!(lines_final_l, lines1l[l_ind])
        push!(lines_final_c, lines1c[l_ind])
        # delete this line from the RP
        delete_line_from_cl_ret_RP!(X_vertical1, lines1t[l_ind], lines1l[l_ind], lines1c[l_ind])
        # go along each point of the line and check for neighbours
        l_max = lines1t[l_ind]

        @inbounds for l = 1:l_max
            # scan each line twice - above and underneth
            for index = -1:2:1
                # make sure not to exceed RP-boundaries
                (columni+index > M || columni+index == 0) ? break : nothing
                # if there is a neighbouring point, call recursive scan-function
                (X_vertical1[linei+l-1, columni+index] != 0) ? scan_lines!(X_vertical1, lines_copy1t,
                                lines_copy1l, lines_copy1c, linei+l-1, columni+index) : nothing
            end
        end
    end
    return lines_final_t, lines_final_l, lines_final_c
end

# compute a `3-by-N` matrix, which stores all vertical lines from a converted
# recurrence matrix (after applying `create_close_returns_map` to an ARM). The
# line lengths of each found line are stored in the first row, and the starting
# indices of its row and column, are stored in the 2nd and 3rd row, respectively.
# Returns the lines of this sorted matrix as vectors
function verticalhisto(R::SparseMatrixCSC)
    M, N = size(R)
    lengthvector = Int[]
    columnvector = Int[]
    startvector = Int[]
    Rjdiffs = zeros(Int,M+1)
    @inbounds for j = 1:N
        extendeddiff!(Rjdiffs, @view R[:,j])
        starts = findall(isequal(1), Rjdiffs)
        ends = findall(isequal(-1), Rjdiffs)
        linelengths = ends .- starts
        mask = (!iszero).(linelengths)
        append!(lengthvector, view(linelengths, mask))
        append!(startvector, view(starts, mask))
        append!(columnvector, repeat([j], count(mask)))
    end
    ordered = sortperm(lengthvector; rev=true)
    return lengthvector[ordered], startvector[ordered], columnvector[ordered]
end

# in-place extension of diff, adding diffs of initial and final points
function extendeddiff!(d, x::AbstractVector{T}) where T
    n = length(x)
    d[1] = prev = x[1]
    @inbounds for i in 2:n
            next = x[i]
            d[i] = next-prev
            prev = next
    end
    d[n+1] = -prev
    return d
end

function verticalhisto2(R::SparseMatrixCSC)
    rows = rowvals(R)
    # early return
    isempty(rows) && return (Int[], Int[], Int[])
    cols = colvals(R)
    n = length(rows)
    lengthvector = [0]
    columnvector = cols[1:1]
    startvector = rows[1:1]
    # Iterate over columns
    cprev = cols[1]
    r1 = rows[1]
    rprev = r1 - 1 # coerce that (a) is not hit in the first iteration
    @inbounds for i=1:n
        r = rows[i]
        c = cols[i]
        # Search the second and later segments in the column
        if c == cprev
            if r-rprev != 1 # (a): there is a separation between rprev and r
                # update histogram of segments
                extend_skeleton_histogram!(lengthvector, startvector, columnvector, r, c, rprev-r1+1)
                r1 = r # update the start of the next segment
            end
            rprev = r  # update the previous position
        else # hit in the first point of a new column
            # process the last fragment of the previous column
            extend_skeleton_histogram!(lengthvector, startvector, columnvector, r, c, rprev-r1+1)
            # initialize values for searching new fragments
            cprev = c
            r1 = r
            rprev = r
        end
    end
    # Process the latest fragment
    lengthvector[end] = rprev-r1+1

    ordered = sortperm(lengthvector; rev=true)
    return lengthvector[ordered], startvector[ordered], columnvector[ordered]
end

@inline function extend_skeleton_histogram!(lengthvector, startvector, columnvector, r, c, p)
    lengthvector[end] = p
    push!(lengthvector, 0)
    push!(startvector, r)
    push!(columnvector, c)
end

# manipulates XX inplace and returns a view on l_vec
function scan_lines!(XX::SparseMatrixCSC, l_vec1::AbstractVector{<:Integer}, l_vec2::AbstractVector{<:Integer}, l_vec3::AbstractVector{<:Integer}, line::Integer, column::Integer)

    # for the input index tuple look for the start indices
    index = 0
    local del_ind::Int
    while true
        # check whether the input index tuple is a listed index for starting
        # points of line lengths in the line matrix
        loc_line = findall(column .== l_vec3)
        found = findfirst(isequal(line+index), l_vec2[loc_line])
        if found !== nothing
            del_ind = loc_line[found]
            break
        end
        index -= 1
    end
    # delete the line from RP
    delete_line_from_cl_ret_RP!(XX, l_vec1[del_ind], l_vec2[del_ind], l_vec3[del_ind])

    len, li , co =  l_vec1[del_ind], l_vec2[del_ind], l_vec3[del_ind]

    #delete the line from the line matix, i.e. set to zero
    deleteat!(l_vec1, del_ind)
    deleteat!(l_vec2, del_ind)
    deleteat!(l_vec3, del_ind)

    # check for borders of the RP
    for i = 1:len
        newli, newco = neighborindices(XX, i, li, co)
        if newli != 0 && newco != 0
            scan_lines!(XX, l_vec1, l_vec2, l_vec3, newli, newco)
        end
    end
end

function neighborindices(XX, i::T, li::T, co::T) where T<:Integer
    N, M = size(XX)
    for newli in (li+i .+ (-2:0)), newco in (co-1, co+1)
        if (1 ≤ newli ≤ N) && (1 ≤ newco ≤ M)
            (XX[newli, newco] != 0) && return (newli, newco)
        end
    end
    return zero(T), zero(T) # default if the conditions are not met
end
