# Recurrence parameters as defined by Marwan et al. (2007)

recurrencerate(x::SparseMatrixCSC) = nnz(x)/prod(size(x))

function recurrencerate(x::SparseMatrixCSC; theiler::Integer=0)
    theiler_points = theiler > 0 ? length(diag(x)) : 0
    theiler_nz = theiler > 0 ? nnz(diag(x)) : 0
    for d = 1:theiler-1
        theiler_points += 2length(diag(x,d))
        theiler_nz += 2nnz(diag(x,d))
    end
    (nnz(x)-theiler_nz)/(prod(size(x))-theiler_points)
end

localrecurrence(x::SparseMatrixCSC{Bool}) =
    [nnz(diag(x,d)) for d in (1:minimum(size(x))-1)]

# Based on diagonal lines

function diagonalhistogram(x::SparseMatrixCSC{Bool}; theiler::Integer=1, kwargs...)
    bins = [0]
    nbins = 1
    current_diag = 0
    n = minimum(size(x))
    # Iterate over diagonals - excluding LOI and Theiler window
    for d = theiler:n-2
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
function diagonalhistogram(x::AbstractMatrix{Bool}; theiler::Integer=1, kwargs...)
    bins = [0]
    nbins = 1
    n = minimum(size(x))
    # Iterate over diagonals
    for d = theiler:n-2
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

function determinism(diag_hist::AbstractVector; lmin=2, kwargs...)
    if lmin < 2
        error("lmin must be 2 or higher")
    end
    nbins = length(diag_hist)
    diag_points = collect(1:nbins) .* diag_hist
    (lmin > nbins) ? 0.0 : sum(diag_points[lmin:nbins])/sum(diag_points) 
end

function determinism(x::AbstractMatrix; kwargs...)
    determinism(diagonalhistogram(x; kwargs...); kwargs...)
end

function avgdiag(diag_hist::AbstractVector; lmin=2, kwargs...)
    if lmin < 2
        error("lmin must be 2 or higher")
    end
    nbins = length(diag_hist)
    diag_points = collect(1:nbins) .* diag_hist
    (lmin > nbins) ? 0.0 : sum(diag_points[lmin:nbins])/sum(diag_hist[lmin:nbins]) 
end

function avgdiag(x::AbstractMatrix; kwargs...)
    avgdiag(diagonalhistogram(x; kwargs...); kwargs...)
end

maxdiag(diag_hist::AbstractVector) = length(diag_hist)
maxdiag(x::AbstractMatrix; kwargs...) = maxdiag(diagonalhistogram(x; kwargs...))

divergence(x) = 1/maxdiag(x)

function entropy(diag_hist::AbstractVector; lmin=2, kwargs...)
    if lmin < 2
        error("lmin must be 2 or higher")
    end
    nbins = length(diag_hist)
    if lmin <= nbins
        prob_bins = diag_hist[lmin:nbins] ./ sum(diag_hist[lmin:nbins])
        prob_bins = prob_bins[find(prob_bins)]
        -sum(prob_bins .* log(prob_bins))
    else
        0.0
    end
end

function entropy(x::AbstractMatrix; kwargs...)
    entropy(diagonalhistogram(x; kwargs...); kwargs...)
end

function trend(npoints::AbstractVector; theiler=1, border=10, kwargs...)
    nmax = length(npoints)
    rrk = npoints./collect(nmax:-1:1)
    m = nmax-border
    w = collect(theiler:m)-m/2
    w'*(rrk[theiler:m]-mean(rrk[theiler:m])) / (w'*w)
end

function trend(x::AbstractMatrix; kwargs...)
    trend(localrecurrence(x); kwargs...)
end

# Number of l-length sequences, based on diagonals
function countsequences(diag_hist::AbstractVector; lmin=2, kwargs...)
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


laminarity(vert_hist::AbstractVector; kwargs...) = determinism(vert_hist; kwargs...)
laminarity(x::AbstractMatrix; kwargs...) = laminarity(verticalhistogram(x); kwargs...)

trappingtime(vert_hist::AbstractVector; kwargs...) = avgdiag(vert_hist; kwargs...)
trappingtime(x::AbstractMatrix; kwargs...) = trappingtime(verticalhistogram(x); kwargs...)

maxvert(vert_hist::AbstractVector) = length(vert_hist)
maxvert(x::AbstractMatrix) = maxvert(verticalhistogram(x))

