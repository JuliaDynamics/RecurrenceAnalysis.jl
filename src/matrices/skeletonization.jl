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
        X_cl1 = convert_recurrence_matrix(tril(X))
        # get "horizontal" line distribution with position indices
        lines1t, lines1l, lines1c = horizontalhisto(X_cl1)
        lines_copy1t = deepcopy(lines1t)
        lines_copy1l = deepcopy(lines1l)
        lines_copy1c = deepcopy(lines1c)

        # create a close returns map with horizontal lines represented by
        # numbers, equal to their lengths
        X_hori1 = create_close_returns_map(lines_copy1t, lines_copy1l, lines_copy1c, size(X_cl1))
    else
        symm = false
        # convert upper and lower triangles into close returns maps
        X_cl1 = convert_recurrence_matrix(tril(X))
        X_cl2 = convert_recurrence_matrix(sparse(triu(X)'))

        # get "horizontal" line distribution with position indices
        lines1t, lines1l, lines1c = horizontalhisto(X_cl1)
        lines_copy1t = deepcopy(lines1t)
        lines_copy1l = deepcopy(lines1l)
        lines_copy1c = deepcopy(lines1c)
        lines2t, lines2l, lines2c = horizontalhisto(X_cl2)
        lines_copy2t = deepcopy(lines2t)
        lines_copy2l = deepcopy(lines2l)
        lines_copy2c = deepcopy(lines2c)

        # create close returns maps with horizontal lines represented by
        # numbers, equal to their lengths
        X_hori1 = create_close_returns_map(lines_copy1t, lines_copy1l, lines_copy1c, size(X_cl1))
        X_hori2 = create_close_returns_map(lines_copy2t, lines_copy2l, lines_copy2c, size(X_cl2))
    end

    # scan the lines, start with the longest one and discard all adjacent lines
    lines_final_t, lines_final_l, lines_final_c = get_final_line_matrix(lines1t, lines1l, lines1c, lines_copy1t, lines_copy1l, lines_copy1c, X_hori1)

    # if not symmetric input RP, than compute for the upper triangle as well
    if !symm
        lines_final2_t, lines_final2_l, lines_final2_c = get_final_line_matrix(lines2t, lines2l, lines2c, lines_copy2t, lines_copy2l, lines_copy2c, X_hori2)
        # build RP based on the histogramm of the reduced lines
        X_new = build_skeletonized_RP(lines_final_t, lines_final_l, lines_final_c, lines_final2_t, lines_final2_l, lines_final2_c, size(X_hori1,1), size(X_hori1,2))
    else
        # build RP based on the histogramm of the reduced lines
        X_new = build_skeletonized_RP(lines_final_t, lines_final_l, lines_final_c, size(X_hori1,1), size(X_hori1,2))
    end
    return X_new
end

# Auxiliary functions

# Transforms the standard RP into a close returns map
function convert_recurrence_matrix(R::SparseMatrixCSC; triangle::Bool = true)
    if triangle
        N = size(R)
        # init new matrix
        Y = zeros(Bool, N[1]+1, N[1])
        # fill rows of Y with the diagonals of R
        # lower triangle
        @inbounds for i = 0:N[1]-1
           Y[N[1]-i+1,(1:(N[1]-i)).+i] = view(R, diagind(R, -i))
        end
    else
        N = size(R)
        # init new matrix
        Y = zeros(Bool, 2*N[1]+1,N[1])
        # fill rows of Y with the diagonals of R
        # upper triangle
        @inbounds for i = 0:N[1]-1
           Y[N[1]+i+1,(1:(N[1]-i))] = view(R, diagind(R, i))
        end
        # lower triangle
        @inbounds for i = 0:N[1]-1
           Y[N[1]-i+1,(1:(N[1]-i)).+i] = view(R, diagind(R, -i))
        end
    end
    return sparse(Y)
end

# Transforms the reverted RP (close returns map) into a normal RP
function revert_close_returns_map(R::SparseMatrixCSC; triangle::Bool = true)
    if triangle
        N = size(R)
        # init new matrix
        Y = zeros(Bool, N[2], N[2])
        # make R to a square matrix, fill the new part with zeros
        Z = [R ; zeros(N[1]-1,N[2])]
        Z = Z[end:-1:1,:]

        # fill columns of Y with the diagonals of Z (but only the first N points)
        @inbounds for i = 1:N[2]
            di = diag(Z,-i)
            Y[:,N[2]-i+1] = di[1:N[2]]
        end
    else
        N = size(R)
        # init new matrix
        Y = zeros(Bool, N[2], N[2])
        # make R to a square matrix, fill the new part with zeros
        Z = [R zeros(N[1],N[1]+1)]
        Z = Z[end:-1:1,:]

        # fill columns of Y with the diagonals of Z (but only the first N points)
        @inbounds for i = 1:N[2]
            di = diag(Z,-i)
            Y[:,N[2]-i+1] = di[1:N[2]]
        end
    end
    return sparse(Y)
end

# from the 3-by-N Matrix, which stores all horizontal lines in the reverted RP,
# (close returns map) [obtained by calling first `revert_recurrence_matrix` on a
# recurrence Matrix R, and consecutevely `horizontalhisto(R)`] construct and return
# a "dummy" close returns map, where all lines are encoded according to their length.
function create_close_returns_map(lines1::Vector{Int}, lines2::Vector{Int}, lines3::Vector{Int}, size_of_cl_ret_RP::Tuple{Int, Int})
    X = zeros(Int, size_of_cl_ret_RP)
    for i = 1:length(lines1)
        line_length = lines1[i]
        line_ind = lines2[i]
        column_ind = lines3[i]
        for j = 0:line_length-1
            X[line_ind,column_ind+j] = line_length
        end
    end
    return sparse(X)
end

# deletes a line, specified in 'l_vec' (line vector, with first line being
# the total line length, the second line the line-index of the starting point
# and the third line the column-index of the starting point) from the close
# returns "RP" 'R'.
function delete_line_from_cl_ret_RP!(R::Union{ARM,SparseMatrixCSC}, line1::Int, line2::Int, line3::Int)
    R[line2, line3 .+ (1:line1).-1] .= 0
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
            X_cl_new[linei, columni+j-1] = true
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
            X_cl_new[linei, columni+j-1] = true
        end
    end
    @inbounds for i = 1:length(lines11)
        l_max = lines11[i]
        linei = lines22[i]
        columni = lines33[i]
        for j = 1:l_max
            X_cl2_new[linei,columni+j-1] = true
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
    lines_copy1c::Vector{Int}, X_hori1::SparseMatrixCSC)

    N, M = size(X_hori1)
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
        delete_line_from_cl_ret_RP!(X_hori1, lines1t[l_ind], lines1l[l_ind], lines1c[l_ind])
        # go along each point of the line and check for neighbours
        l_max = lines1t[l_ind]

        @inbounds for l = 1:l_max
            # scan each line twice - above and underneth
            for index = -1:2:1
                # make sure not to exceed RP-boundaries
                (linei+index > N || linei+index == 0) ? break : nothing
                # if there is a neighbouring point, call recursive scan-function
                (X_hori1[linei+index, columni+l-1] != 0) ? scan_lines!(X_hori1, lines_copy1t,
                                lines_copy1l, lines_copy1c, linei+index, columni+l-1) : nothing
            end
        end
    end
    return lines_final_t, lines_final_l, lines_final_c
end

# compute a `3-by-N` matrix, which stores all horizontal lines from a converted
# recurrence matrix (after applying `convert_recurrence_matrix` to a ARM). The
# line lengths of each found line are stored in the first row, and the starting
# indices of its row and column, are stored in the 2nd and 3rd row, respectively.
function horizontalhisto(R::SparseMatrixCSC)

    N = size(R)[1]
    liness = zeros(Int, 3, 1)
    @inbounds for j = 1:N
        starts = findall(diff([0; @view R[j,:]]).==1)
        ends = findall(diff([@view R[j,:]; 0]).==-1)

        if ~isempty(starts)
            lines = zeros(Int, 3, length(starts))
            for n=1:length(starts)
                lines[2,n] = j
                lines[3,n] = starts[n]
                lines[1,n] = ends[n] - starts[n] + 1
            end
            liness = hcat(liness,lines)
        end
    end
    # remove lines of length zero (=no line)
    no_zero_lines = findall(liness[1,:].!=0)
    liness = liness[:, no_zero_lines]
    ls = size(liness)
    liness = liness[:,sortperm(@view liness[1, :]; rev=true)]
    return vec(liness[1,:]), vec(liness[2,:]), vec(liness[3,:])
end

# function horizontalhisto2(R::SparseMatrixCSC)
#     # Transpose R in order to get the "horizontal" lines
#     # rows = colvals(R[2:end,:])
#     # cols = rowvals(R[2:end,:])
#     rows = colvals(R)
#     cols = rowvals(R)
#     p = sortperm(cols)
#     rows = rows[p]
#     cols = cols[p]
#     # check bounds
#     n = length(rows)
#     if length(cols) != n
#         throw(ErrorException("mismatch between number of row and column indices"))
#     end
#     # histogram for lines with start & end indices
#     liness = zeros(Int, 3, 1)
#     # Iterate over columns
#     cprev = cols[1]
#     r1 = rows[1]
#     rprev = r1
#     @inbounds for i=1:n
#         r = rows[i]
#         c = cols[i]
#         # Search the second and later segments in the column
#         if c == cprev
#             if r-rprev != 1 # (a): there is a separation between rprev and r
#                 # update histogram of segments
#                 current_vert = rprev-r1+1
#                 if current_vert ≥ 1
#                     liness = extend_skeleton_histogram!(liness, r, c, current_vert)
#                 end
#                 r1 = r # update the start of the next segment
#             end
#             rprev = r  # update the previous position
#         else # hit in the first point of a new column
#             # process the last fragment of the previous column
#             current_vert = rprev-r1+1
#             if current_vert ≥ 1
#                 liness = extend_skeleton_histogram!(liness, r, c, current_vert)
#             end
#             # initialize values for searching new fragments
#             cprev = c
#             r1 = r
#             rprev = r
#         end
#     end
#     # Process the latest fragment
#     current_vert = rprev-r1+1
#     liness[1,end] = current_vert
#     # process the first dummy-fragment
#     liness[1,1] = 0
#
#     # remove lines of length zero (=no line)
#     no_zero_lines = findall(liness[1,:].!=0)
#     liness = liness[:, no_zero_lines]
#     ls = size(liness)
#     liness = liness[:,sortperm(@view liness[1, :]; rev=true)]
#     return vec(liness[1,:]), vec(liness[2,:]), vec(liness[3,:])
# end

# manipulates XX inplace and returns a view on l_vec
function scan_lines!(XX::SparseMatrixCSC, l_vec1::AbstractVector{<:Integer}, l_vec2::AbstractVector{<:Integer}, l_vec3::AbstractVector{<:Integer}, line::Integer, column::Integer)

    # for the input index tuple look for the start indices
    index = 0
    local del_ind::Int
    while true
        # check whether the input index tuple is a listed index for starting
        # points of line lengths in the line matrix
        loc_line = findall(line .== l_vec2)
        found = findfirst(isequal(column+index), l_vec3[loc_line])
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
    for newli in (li-1, li+1), newco in (co+i .+ (-2:0))
        if (1 ≤ newli ≤ N) && (newco != 0) && (newco ≤ M)
            (XX[newli, newco] != 0) && return (newli, newco)
        end
    end
    return zero(T), zero(T) # default if the conditions are not met
end
