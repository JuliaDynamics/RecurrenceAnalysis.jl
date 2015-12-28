module RecurrenceAnalysis

export RecurrenceMatrix

# vnorm function can be redefined
vnorm = x -> norm(x, Inf)

function RecurrenceMatrix(x::Array{Float64,1}, d::Int, k::Int, radius::Float64)
    n = length(x)
    # Create k-vectors with delay d
    # Each k-vector is in a row
    nk = n-d*(k-1)
    nk < 2 && error("Input is too short: length must be greater than k*d")
    xk = x[collect(1:nk) .+ collect(0:d:d*(k-1))']
    # Sort by the first dimension
    ix = sortperm(xk[:,1])
    xk_sort = xk[ix,:]
    # Create and fill sparse matrix (upper triangle)
    rmat = map(Bool, speye(nk))
    first_neighbour = collect(2:nk)
    for p = 1:nk-1
        # First iteration in nearest point in the first dimension 
        q0 = q = first_neighbour[p]
        dq = vnorm(xk_sort[q,:] - xk_sort[p,:])
        match_all = false
        if dq < radius
            rmat[p,q] = rmat[q,p] = true
            match_all = (dq < 0.5radius)
        end
        # Continue with following points
        while ((q+=1) <= nk) && (xk_sort[q,1] - xk_sort[p,1] < radius)
            dq = vnorm(xk_sort[q,:] - xk_sort[p,:])
            if dq < radius
                rmat[p, q] = rmat[q, p] = true
                # Mark all crossings within half radius
                if match_all && dq < 0.5radius
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

function RecurrenceMatrix(x::AbstractVecOrMat, d::Real, k::Real, radius::Real)
    if length(x) != maximum(size(x))
        error("`x` must be a row or column vector")
    end
    RecurrenceMatrix(float(x[:]), Int(d), Int(k), float(radius))
end

function RecurrenceMatrixBruteForce(x::Array{Float64,1}, d::Int, k::Int, radius::Float64)
    n = length(x)
    # Create k-vectors with delay d
    # Each k-vector is in a row
    nk = n-d*(k-1)
    nk < 2 && error("Input is too short: length must be greater than k*d")
    xk = x[collect(1:nk) .+ collect(0:d:d*(k-1))']
    rmat = map(Bool, speye(nk))
    for r = 1:nk
        for c = (r+1):nk
            dist = vnorm(xk[r,:] - xk[c,:])
            if dist < radius
                rmat[r,c] = rmat[c,r] = true
            end
        end
    end
    rmat
end

# Benchmark: 
# The optimized algorithm is about 3 times faster and takes about 9.5 times less
# memory than the brute force algorithm, for a matrix with sparsity ratio ~ 0.06

function diagonal_bins(x::SparseMatrixCSC{Bool})
    bins = [0]
    nbins = 1
    current_diag = 0
    n = minimum(size(x))
    # Iterate over diagonals
    for d = 1:n-2
        previous_cell = false
        for r = 1:n-d
            # Operate on existing diagonal line
            if previous_cell
                # Extend length with current cell
                (extend = x[r,r+d]) && (current_diag += 1)
                # If arrived to the end of a line
                # add the current length to the corresponding bin
                if (!extend || (r == n-d)) && (current_diag > 0)
                    # Append new positions to the bins if needed
                    if current_diag > nbins
                        append!(bins, zeros(current_diag - nbins))
                        nbins = current_diag
                    end
                    bins[current_diag] += 1
                    current_diag = 0
                end
            end
            previous_cell = x[r,r+d]
        end
    end
    # Add isolated points in first bin
    [nnz(triu(x,1)) - sum(bin); bins]
end

recurrence_rate(x::SparseMatrixCSC{Bool}) = nnz(x)/prod(size(R))

function determinism(diag_bins::Array{Int}, lmin::Int)
    if lmin < 2
        error("lmin must be 2 or higher")
    end
    nbins = length(diag_bins)
    (lmin > nbins) ? 0.0 : collect(lmin:nbins)'*bins[lmin:end]/sum(diag_bins) 
end

determinism(x::SparseMatrixCSC{Bool}, lmin) = determinism(diagonal_bins(x), lmin)

avgdiag(diag_bins::Array{Int}) = mean(diag_bins)
avgdiag(x::SparseMatrixCSC{Bool}) = avgdiag(diagonal_bins(x))

maxdiag(diag_bins::Array{Int}) = length(diag_bins)
maxdiag(x::SparseMatrixCSC{Bool}) = maxdiag(diagonal_bins(x))

divergence(x) = 1/maxdiag(x)

function entropy(diag_bins::Array{Int}, lmin::Int)
    if lmin < 2
        error("lmin must be 2 or higher")
    end
    nbins = length(diag_bins)
    if lmin <= nbins
        prob_bins = diag_bins[lmin:end]/sum(diag_bins[lmin:end])
        -sum(prob_bins .* log2(prob_bins))
    else
        0.0
    end
end
