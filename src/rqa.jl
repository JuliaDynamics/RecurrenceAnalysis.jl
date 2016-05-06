# Recurrence parameters as defined by Marwan et al. (2007)

recurrencerate(x::AbstractMatrix) = countnz(x)/prod(size(x))

function recurrencerate(x::AbstractMatrix; theiler::Integer=0)
    theiler < 0 && error("Theiler window length must be greater than 0")
    diags_remove = -(theiler-1):(theiler-1)
    theiler_points = 0
    theiler_nz = 0
    for d in diags_remove
        theiler_points += length(diag(x,d))
        theiler_nz += countnz(diag(x,d))
    end
    (countnz(x)-theiler_nz)/(prod(size(x))-theiler_points)
end

function tau_recurrence(x::AbstractMatrix{Bool})
    n = minimum(size(x))
    [countnz(diag(x,d))/(n-d) for d in (1:n-1)]
end

# Based on diagonal lines

function diagonalhistogram(x::AbstractMatrix{Bool}; theiler::Integer=1, kwargs...)
    theiler < 0 && error("Theiler window length must be greater than 0")
    bins = [0]
    nbins = 1
    current_diag = 0
    m, n = size(x)
    # Iterate over diagonals - excluding LOI and Theiler window
    # If the matrix is symmetric, examine only the upper triangle
    diag_collection = collect(theiler:n-2)
    xsym = issym(x)
    !xsym && prepend!(diag_collection, collect(-(m-2):-max(theiler,1)))
    for d in diag_collection
        increment = (xsym && d > 0) ? 2 : 1
        previous_cell = false
        first_c = max(1, d+1)
        last_c = min(n, m+d)
        for c = first_c:last_c
            # Operate on existing diagonal line
            if previous_cell
                # Extend length with current cell
                (extend = x[c-d,c]) && (current_diag += 1)
                # If arrived to the end of a line
                # add the current length to the corresponding bin
                if (!extend || (c == last_c)) && (current_diag > 0)
                    # Append new positions to the bins if needed
                    if current_diag > nbins
                        append!(bins, zeros(current_diag - nbins))
                        nbins = current_diag
                    end
                    bins[current_diag] += increment
                    current_diag = 0
                end
            end
            previous_cell = x[c-d,c] 
        end
    end
    # Add isolated points in first bin
    allpoints = (theiler == 0) ? countnz(x) : countnz(triu(x, theiler)) + countnz(tril(x,-theiler))
    [allpoints - collect(2:nbins+1)'*bins; bins]
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
    avgdiag(diagonalhistogram(x; kwargs); kwargs...)
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
        -sum(prob_bins .* log2(prob_bins))
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
    trend(tau_recurrence(x); kwargs...)
end

# Number of l-length sequences, based on diagonals
function countsequences(diag_hist::AbstractVector; lmin=2, kwargs...)
    overlap = (1:length(diag_hist))' - lmin + 1
    overlap[overlap .< 0] = 0
    overlap * diag_hist
end

function countsequences(x::AbstractMatrix; kwargs...)
    countsequences(diagonalhistogram(x; kwargs...); kwargs...)
end

# 2. Based on vertical lines

function verticalhistogram(x::AbstractMatrix{Bool})
    bins = [0]
    nbins = 1
    current_vert = 0
    m, n = size(x)
    # Iterate over columns
    for c = 1:n
        previous_cell = false
        for r = 1:m
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
    [countnz(x) - collect(2:nbins+1)'*bins; bins]
end

laminarity(vert_hist::AbstractVector; kwargs...) = determinism(vert_hist; kwargs...)
laminarity(x::AbstractMatrix; kwargs...) = laminarity(verticalhistogram(x); kwargs...)

trappingtime(vert_hist::AbstractVector; kwargs...) = avgdiag(vert_hist; kwargs...)
trappingtime(x::AbstractMatrix; kwargs...) = trappingtime(verticalhistogram(x); kwargs...)

maxvert(vert_hist::AbstractVector) = length(vert_hist)
maxvert(x::AbstractMatrix) = maxvert(verticalhistogram(x))

