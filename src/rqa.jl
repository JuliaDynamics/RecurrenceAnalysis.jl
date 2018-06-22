# Recurrence parameters as defined by Marwan et al. (2007)

"""
    recurrencerate(x, theiler=0)
    
Calculate the recurrence rate (RR) of a recurrence matrix, ruling out
the points within the Theiler window.
"""
function recurrencerate(x::AbstractMatrix; theiler::Integer=0, kwargs...)
    theiler < 0 && error("Theiler window length must be greater than 0")
    if theiler == 0
        return typeof(0.0)( count(!iszero, x)/prod(size(x)) )
    end
    diags_remove = -(theiler-1):(theiler-1)
    theiler_points = 0
    theiler_nz = 0
    for d in diags_remove
        theiler_points += length(diag(x,d))
        theiler_nz += count(!iszero, diag(x,d))
    end
    typeof(0.0)( (count(!iszero, x)-theiler_nz)/(prod(size(x))-theiler_points) )
end

function tau_recurrence(x::AbstractMatrix{Bool})
    n = minimum(size(x))
    [count(!iszero, diag(x,d))/(n-d) for d in (1:n-1)]
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
    xsym = issymmetric(x)
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
    allpoints = (theiler == 0) ? count(!iszero, x) : count(!iszero, triu(x, theiler)) + count(!iszero, tril(x,-theiler))
    [allpoints - collect(2:nbins+1)'*bins; bins]
end

function diagonalhistogram(x::SparseMatrixCSC{Bool}; theiler::Integer=1, kwargs...)
    theiler < 0 && error("Theiler window length must be greater than 0")
    m,n=size(x)
    rv = rowvals(x)
    dv = colvals(x) - rowvals(x)
    if issymmetric(x)
        valid = (dv .>= theiler)
        f = 2
    else
        valid = (abs.(dv) .>= theiler)
        f = 1
    end
    vmat = sparse(rv[valid], dv[valid] .+ (m+1), true)
    f .* verticalhistogram(vmat)
end

"""
    determinism(x; lmin=2, theiler=1)
    
Calculate the determinism (DET) of a recurrence matrix, ruling out
the points within the Theiler window and diagonals shorter than a minimum value.
"""
function determinism(diag_hist::Vector; lmin=2, kwargs...)
    if lmin < 2
        error("lmin must be 2 or higher")
    end
    nbins = length(diag_hist)
    diag_points = collect(1:nbins) .* diag_hist
    (lmin > nbins) ? 0.0 : typeof(0.0)( sum(diag_points[lmin:nbins])/sum(diag_points) )
end

function determinism(x::AbstractMatrix; kwargs...)
    determinism(diagonalhistogram(x; kwargs...); kwargs...)
end

"""
    avgdiag(x; lmin=2, theiler=1)
    
Calculate the average diagonal length (L) in a recurrence matrix, ruling out
the points within the Theiler window and diagonals shorter than a minimum value.
"""
function avgdiag(diag_hist::Vector; lmin=2, kwargs...)
    if lmin < 2
        error("lmin must be 2 or higher")
    end
    nbins = length(diag_hist)
    diag_points = collect(1:nbins) .* diag_hist
    (lmin > nbins) ? 0.0 : typeof(0.0)( sum(diag_points[lmin:nbins])/sum(diag_hist[lmin:nbins]) )
end

function avgdiag(x::AbstractMatrix; kwargs...)
    avgdiag(diagonalhistogram(x; kwargs...); kwargs...)
end

"""
    maxdiag(x; theiler=1)
    
Calculate the longest diagonal (Lmax) in a recurrence matrix, ruling out
the points within the Theiler window.
"""
maxdiag(diag_hist::Vector; kwargs...) = length(diag_hist)
maxdiag(x::AbstractMatrix; kwargs...) = maxdiag(diagonalhistogram(x; kwargs...))

"""
    divergence(x; lmin=2, theiler=1)
    
Calculate the divergence of a recurrence matrix
(actually the inverse of `maxdiag`.
"""
divergence(x; kwargs...) = typeof(0.0)( 1/maxdiag(x; kwargs...) )

"""
    entropy(x; lmin=2, theiler=1)
    
Calculate the entropy of diagonal lengths (ENT) of a recurrence matrix, ruling out
the points within the Theiler window and diagonals shorter than a minimum value.
"""
function entropy(diag_hist::Vector; lmin=2, kwargs...)
    if lmin < 2
        error("lmin must be 2 or higher")
    end
    nbins = length(diag_hist)
    if lmin <= nbins
        prob_bins = diag_hist[lmin:nbins] ./ sum(diag_hist[lmin:nbins])
        prob_bins = prob_bins[findall(!iszero, prob_bins)]
        typeof(0.0)( -sum(prob_bins .* log.(prob_bins)) )
    else
        0.0
    end
end

function entropy(x::AbstractMatrix; kwargs...)
    entropy(diagonalhistogram(x; kwargs...); kwargs...)
end

"""
    trend(x; theiler=1, border=10)
    
Calculate the trend of recurrences in recurrence matrix towards its edges, ruling out
the points within the Theiler window and in the outermost diagonals.
"""
function trend(npoints::Vector; theiler=1, border=10, kwargs...)
    nmax = length(npoints)
    rrk = npoints./collect(nmax:-1:1)
    m = nmax-border
    w = collect(theiler:m) .- m/2
    typeof(0.0)( (w'*(rrk[theiler:m] .- mean(rrk[theiler:m])) ./ (w'*w))[1] )
end

function trend(x::AbstractMatrix; kwargs...)
    trend(tau_recurrence(x); kwargs...)
end

# Number of l-length sequences, based on diagonals
function countsequences(diag_hist::Vector; lmin=2, kwargs...)
    overlap = (1:length(diag_hist))' .- (lmin+1)
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
    [count(!iszero, x) - collect(2:nbins+1)'*bins; bins]
end

function verticalhistogram(x::SparseMatrixCSC{Bool})
    m,n=size(x)
    bins = [0]
    nbins = 1
    rv = rowvals(x)
    # Iterate over columns
    for d = 1:n
        rvd = rv[nzrange(x,d)]
        nd = length(rvd)
        if nd>1
            r1 = rvd[1]
            rprev = r1
            for r in rvd[2:end]
                # Look for nonzero that starts a new column fragment
                # (more than one row after the previous one)
                if r-rprev != 1
                    current_diag = rprev-r1+1
                    if current_diag > nbins
                        append!(bins, zeros(current_diag-nbins))
                        nbins = current_diag
                    end
                    bins[current_diag] += 1
                    r1 = r
                end
                rprev = r
            end
            # Last column fragment
            if rprev-rvd[end-1] == 1
                current_diag = rprev-r1+1
                if current_diag > nbins
                    append!(bins, zeros(current_diag-nbins))
                    nbins = current_diag
                end
                bins[current_diag] += 1
            else
                bins[1] += 1
            end
        elseif nd==1
            bins[1] += 1
        end
    end
    bins
end

"""
    laminarity(x; lmin=2)
    
Calculate the laminarity (LAM) of a recurrence matrix, ruling out
vertical lines shorter than a minimum value.
"""
laminarity(vert_hist::Vector; kwargs...) = determinism(vert_hist; kwargs...)
laminarity(x::AbstractMatrix; kwargs...) = laminarity(verticalhistogram(x); kwargs...)

"""
    trappingtime(x; lmin=2)
    
Calculate the trapping time (TT) of a recurrence matrix, ruling out
vertical lines shorter than a minimum value.
"""
trappingtime(vert_hist::Vector; kwargs...) = avgdiag(vert_hist; kwargs...)
trappingtime(x::AbstractMatrix; kwargs...) = trappingtime(verticalhistogram(x); kwargs...)

"""
    maxvert(x)
    
Calculate the longest vertical line (Vmax) of a recurrence matrix.
"""
maxvert(vert_hist::Vector) = length(vert_hist)
maxvert(x::AbstractMatrix) = maxvert(verticalhistogram(x))

"""
    rqa(x; <keyword arguments>)
"""
function rqa(x; kwargs...)
    dhist = diagonalhistogram(x; kwargs...)
    vhist = verticalhistogram(x)
    Dict("RR"  => recurrencerate(x;kwargs...),
        "DET"  => determinism(dhist;kwargs...),
        "L"    => avgdiag(dhist;kwargs...),
        "Lmax" => maxdiag(dhist),
        "DIV"  => divergence(dhist),
        "ENT"  => entropy(dhist;kwargs...),
        "TND"  => trend(x;kwargs...),
        "LAM"  => laminarity(vhist;kwargs...),
        "TT"   => trappingtime(vhist;kwargs...),
        "Vmax" => maxvert(vhist)
    )
end

