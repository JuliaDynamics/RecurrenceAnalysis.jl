module RecurrenceAnalysis

export embed
export recurrencematrix
export recurrencerate, determinism, avgdiag, maxdiag, divergence, entropy
export trend, laminarity, trappingtime, maxvert

# vnorm function can be redefined
vnorm = x -> norm(x, Inf)

# Embed series
function embed(x::AbstractVecOrMat, delay::Integer, m::Integer)
    dims = size(x)
    n = dims[1]
    nm = n-delay*(m-1)
    if nm < 2 &&
        warning("the emedded time series has length < 2")
    end
    ix = (1:nm) .+ (0:delay:delay*(m-1))'
    embed_indices(x, ix)
end

embed(x, delay, m) = isinteger(delay) && isinteger(m) ?
    embed(x,delay,m) : error("delay and embedding dimension must be integer values")

embed_indices(x::AbstractVector, indices) = x[indices]

function embed_indices(x::AbstractMatrix, indices)
    dx = size(x)
    dxm = size(indices)
    ix_rep = repeat(indices, inner=[1,dx[2]])
    ix_rep += repeat(dx[1]*(0:dx[2]-1)', outer=[dxm...])
    x[ix_rep]
end

# Recurrence matrix creation

"""Create a recurrence matrix from an embedded time series"""
function recurrencematrix(x, radius)
    halfradius = 0.5radius
    n = size(x)[1]
    # Sort by the first dimension
    ix = sortperm(x[:,1])
    x_sort = x[ix,:]
    # Create and fill sparse matrix (upper triangle)
    rmat = map(Bool, speye(n))
    first_neighbour = collect(2:n)
    for p = 1:n-1
        # First iteration at nearest point in the first dimension 
        q0 = q = first_neighbour[p]
        # Exit if all columns were already identified as neighbours
        if q > n
            continue
        end
        dq = vnorm(x_sort[q,:] - x_sort[p,:])
        match_all = false
        if dq < radius
            rmat[p,q] = rmat[q,p] = true
            match_all = (dq < halfradius)
        end
        # Continue with following points
        while ((q+=1) <= n) && (x_sort[q,1] - x_sort[p,1] < radius)
            dq = vnorm(x_sort[q,:] - x_sort[p,:])
            if dq < radius
                rmat[p, q] = rmat[q, p] = true
                # Mark all crossings within half radius
                if match_all && dq < halfradius
                    rmat[q0:q-1, q] = true
                    rmat[q, q0:q-1] = true
                    first_neighbour[q0:q-1] = q+1
                else
                    match_all = false
                end
            else
                match_all = false
            end
        end
    end
    ix_inv = sortperm(ix)
    rmat[ix_inv, ix_inv]
end

"""Create a recurrence matrix from umbedded time series
with specified dimension and delay"""
function recurrencematrix_v(x::AbstractVector, m, delay, radius)
    n = length(x)
    nm = n-delay*(m-1)
    if nm < 2 &&
        warning("the emedded time series has length < 2")
    end
    # Sort by the first dimension
    ix = sortperm(x)
    x_sort = x[ix]
    # Create and fill sparse matrix (upper triangle)
    rmat_full = map(Bool, speye(n))
    last_c = 1
    for p = 1:n-1
        next_c = findnext(x_sort-x_sort[p] .> radius, last_c)
        last_c = next_c == 0 ? n : next_c
        for q = p+1:last_c
            rmat_full[p,q] = true
            rmat_full[q,p] = true
        end
    end
    ix_inv = sortperm(ix)
    rmat_full = rmat_full[ix_inv, ix_inv]
    # Multiply matrices to embed
    rmat = rmat_full[1:nm, 1:nm]
    for k = (1:m-1)
        rmat .*= rmat_full[k*delay+(1:nm), k*delay+(1:nm)]
    end
    rmat
end


"""Brute force algorithm for computation of the recurrence matrix"""
function recurrencematrix_bruteforce(x, radius)
    n = size(x)[1]
    rmat = map(Bool, speye(n))
    for r = 1:n
        for c = (r+1):n
            dist = vnorm(x[r,:] - x[c,:])
            if dist < radius
                rmat[r,c] = rmat[c,r] = true
            end
        end
    end
    rmat
end

# Recurrence parameters as defined by Marwan et al. (2007)

# 1. Based on diagonal lines

function diagonalhistogram(x::SparseMatrixCSC{Bool})
    bins = [0]
    nbins = 1
    current_diag = 0
    n = minimum(size(x))
    # Iterate over diagonals
    for d = 1:n-2
        previous_cell = false
        for c = d+1:n
            # Operate on existing diagonal line
            if previous_cell
                # Extend length with current cell
                (extend = x[c-d,c]) && (current_diag += 1)
                # If arrived to the end of a line
                # add the current length to the corresponding bin
                if (!extend || (c == n)) && (current_diag > 0)
                    # Append new positions to the bins if needed
                    if current_diag > nbins
                        append!(bins, zeros(current_diag - nbins))
                        nbins = current_diag
                    end
                    bins[current_diag] += 1
                    current_diag = 0
                end
            end
            previous_cell = x[c-d,c] 
        end
    end
    # Add isolated points in first bin
    [nnz(triu(x,1)) - collect(2:nbins+1)'*bins; bins]
end

# Alternative definition for full matrices
function diagonalhistogram(x::AbstractMatrix{Bool})
    bins = [0]
    nbins = 1
    n = minimum(size(x))
    # Iterate over diagonals
    for d = 1:n-2
        current_diag = 0
        # Look for continuous sequences (cells where diff .== 1)
        diff_vals = diff(find(diag(x,d)))
        nv = length(diff_vals)
        for v = 1:nv
            # Increment length of current_diag if the diagonal is being extended 
            (extend = (diff_vals[v] == 1)) && (current_diag += 1)
            # If arrived to the end of a line
            # add the current length to the corresponding bin
            if (!extend || (v == nv)) && (current_diag > 0)
                # Make the histogram longer if needed
                if current_diag > nbins
                    append!(bins, zeros(current_diag - nbins))
                    nbins = current_diag
                end
                bins[current_diag] += 1
                current_diag = 0
            end
        end
    end
    # Add isolated points in first bin
    [countnz(triu(x,1)) - collect(2:nbins+1)'*bins; bins]
end

localrecurrence(x::SparseMatrixCSC{Bool}) =
    [nnz(diag(x,d)) for d in (1:minimum(size(x))-1)]

recurrencerate(x::SparseMatrixCSC) = nnz(x)/prod(size(x))

function determinism(diag_hist::AbstractVector, lmin)
    if lmin < 2
        error("lmin must be 2 or higher")
    end
    nbins = length(diag_hist)
    diag_points = collect(1:nbins) .* diag_hist
    (lmin > nbins) ? 0.0 : sum(diag_points[lmin:nbins])/sum(diag_points) 
end

determinism(x::AbstractMatrix, lmin) = determinism(diagonalhistogram(x), lmin)

function avgdiag(diag_hist::AbstractVector, lmin)
    if lmin < 2
        error("lmin must be 2 or higher")
    end
    nbins = length(diag_hist)
    diag_points = collect(1:nbins) .* diag_hist
    (lmin > nbins) ? 0.0 : sum(diag_points[lmin:nbins])/sum(diag_hist[lmin:nbins]) 
end

avgdiag(x::AbstractMatrix, lmin) = avgdiag(diagonalhistogram(x), lmin)

maxdiag(diag_hist::AbstractVector) = length(diag_hist)
maxdiag(x::AbstractMatrix) = maxdiag(diagonalhistogram(x))

divergence(x) = 1/maxdiag(x)

function entropy(diag_hist::AbstractVector, lmin)
    if lmin < 2
        error("lmin must be 2 or higher")
    end
    nbins = length(diag_hist)
    if lmin <= nbins
        prob_bins = diag_hist[lmin:nbins] ./ sum(diag_hist[lmin:nbins])
        prob_bins = prob_bins[find(prob_bins)]
        -sum(prob_bins .* log2(prob_bins))
    else
        0.0
    end
end

entropy(x::AbstractMatrix, lmin) = entropy(diagonalhistogram(x), lmin)

function trend(npoints::AbstractVector, border)
    nmax = length(npoints)
    rrk = npoints./collect(nmax:-1:1)
    m = nmax-border
    w = collect(1:m)-m/2
    w'*(rrk[1:m]-mean(rrk[1:m])) / (w'*w)
end

trend(x::AbstractMatrix, border) = trend(localrecurrence(x), border)

# Number of l-length sequences, based on diagonals
function countsequences(diag_hist::AbstractVector, lmin)
    overlap = (1:length(diag_hist))' - lmin + 1
    overlap[overlap .< 0] = 0
    overlap * diag_hist
end

# 2. Based on vertical lines

function verticalhistogram(x::SparseMatrixCSC{Bool})
    bins = [0]
    nbins = 1
    current_vert = 0
    n = minimum(size(x))
    # Iterate over columns
    for c = 1:n
        previous_cell = false
        for r = 1:n
            if previous_cell
                (extend = x[r,c]) && (current_vert += 1)
                if (!extend || (r == n)) && (current_vert > 0)
                    if current_vert > nbins
                        append!(bins, zeros(current_vert - nbins))
                        nbins = current_vert
                    end
                    bins[current_vert] += 1
                    current_vert = 0
                end
            end
            previous_cell = x[r,c] 
        end
    end
    # Add isolated points in first bin
    [nnz(x) - collect(2:nbins+1)'*bins; bins]
end

# Alternative definition for full matrices

function verticalhistogram(x::AbstractMatrix{Bool})
    bins = [0]
    nbins = 1
    n = minimum(size(x))
    # Iterate over columns
    for c = 1:n
        current_vert = 0
        diff_vals = diff(find(x[r,c]))
        nv = length(diff_vals)
        for v = 1:nv
            (extend = (diff_vals[v] == 1)) && (current_vert += 1)
            if (!extend || (v == nv)) && (current_vert > 0)
                if current_vert > nbins
                    append!(bins, zeros(current_vert - nbins))
                    nbins = current_vert
                end
                bins[current_vert] += 1
                current_vert = 0
            end
        end
    end
    # Add isolated points in first bin
    [countnz(x) - collect(2:nbins+1)'*bins; bins]
end


laminarity(vert_hist::AbstractVector, lmin) = determinism(vert_hist, lmin)
laminarity(x::AbstractMatrix, lmin) = laminarity(verticalhistogram(x), lmin)

trappingtime(vert_hist::AbstractVector, lmin) = avgdiag(vert_hist, lmin)
trappingtime(x::AbstractMatrix, lmin) = trappingtime(verticalhistogram(x), lmin)

maxvert(vert_hist::AbstractVector) = length(vert_hist)
maxvert(x::AbstractMatrix) = maxvert(verticalhistogram(x))

end
