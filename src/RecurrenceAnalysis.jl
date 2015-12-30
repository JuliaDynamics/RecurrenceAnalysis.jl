module RecurrenceAnalysis

export RecurrenceMatrix

# vnorm function can be redefined
vnorm = x -> norm(x, Inf)

function RecurrenceMatrix(x::Array{Float64,1}, d::Integer, k::Integer, radius::Real)
    halfradius = 0.5radius
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
            match_all = (dq < halfradius)
        end
        # Continue with following points
        while ((q+=1) <= nk) && (xk_sort[q,1] - xk_sort[p,1] < radius)
            dq = vnorm(xk_sort[q,:] - xk_sort[p,:])
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

function RecurrenceMatrix(x::AbstractVecOrMat, d::Real, k::Real, radius::Real)
    if length(x) != maximum(size(x))
        error("`x` must be a row or column vector")
    end
    RecurrenceMatrix(float(x[:]), Int(d), Int(k), radius)
end

function RecurrenceMatrixBruteForce(x::Array{Float64,1}, d::Integer, k::Integer, radius::Real)
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

function diagonalstructure(x::SparseMatrixCSC{Bool})
    bins = [0]
    nbins = 1
    current_diag = 0
    n = minimum(size(x))
    npoints = zeros(n-1)
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
            (previous_cell = x[c-d,c]) && (npoints[d] += 1) 
        end
    end
    # Add isolated points in first bin
    bins = [nnz(triu(x,1)) - collect(2:nbins+1)'*bins; bins]
    npoints[n-1] = (x[1,n] ? 1 : 0)
    (bins, npoints)
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

function entropy(diag_bins::Array{Integer}, lmin::Int)
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

# Marwan
function trend(npoints::Integer, border::Integer)
  nmax = length(npoints)
  rrk = npoints./collect(nmax:-1:1)
  m = n-border
  w = (1:m-m/2)
  w'*(rrk[1:m]-mean(rrk[1:m])) / (w'*w)
end

end
