################################################################################
# recurrence_matrix - Low level function
################################################################################
# TODO: increase the efficiency here by not computing everything!

# Core function
# First, we define the non-parallel versions.
# For two datasets

recurrence_matrix(x, metric, ε, parallel) = recurrence_matrix(x, x, metric, ε, parallel)

"""
    recurrence_matrix(x [, y], metric::Metric, ε, parallel::Val)

Return a sparse matrix which encodes recurrence points.
`ε` is the `recurrence_threhsold`, which is a vector for `LocalRecurrenceRate`.

`parallel` may be either `Val(true)` or `Val(false)`.

This is the low-level method that makes the matrices, and it is not part of the public API.
"""
function recurrence_matrix(x::Vector_or_SSSet, y::Vector_or_SSSet, metric::Metric, ε, ::Val{false})
    @assert ε isa Real || length(ε) == length(y)
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
    for j in eachindex(y)
        nzcol = 0
        for i in eachindex(x)
            @inbounds if evaluate(metric, x[i], y[j]) ≤ ( (ε isa Real) ? ε : ε[j] )
                push!(rowvals, i)
                nzcol += 1
            end
        end
        append!(colvals, fill(j, (nzcol,)))
    end
    nzvals = fill(true, (length(rowvals),))
    return sparse(rowvals, colvals, nzvals, length(x), length(y))
end

# For one dataset
function recurrence_matrix(x::Vector_or_SSSet, metric::Metric, ε, ::Val{false})
    @assert ε isa Real || length(ε) == length(x)
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
    for j in eachindex(x)
        nzcol = 0
        for i in 1:j
            @inbounds if evaluate(metric, x[i], x[j]) ≤ ( (ε isa Real) ? ε : ε[j] )
                push!(rowvals, i)
                nzcol += 1
            end
        end
        append!(colvals, fill(j, (nzcol,)))
    end
    nzvals = fill(true, (length(rowvals),))
    return Symmetric(sparse(rowvals, colvals, nzvals, length(x), length(x)), :U)
end


# Now, we define the parallel versions of these functions.

# Core function
function recurrence_matrix(x::Vector_or_SSSet, y::Vector_or_SSSet, metric::Metric, ε, ::Val{true})
    @assert ε isa Real || length(ε) == length(y)
    # We create a Channel for `Array`s of `Array`s, for each thread to have its
    # own array to push to.  This avoids race conditions with
    # multiple threads pushing to the same `Array` (`Array`s are not atomic).
    nbuffers = Threads.nthreads()
    threadchannel = Channel{NTuple{2, Vector{Int}}}(nbuffers) # for rows and columns
    foreach(1:nbuffers) do _
        put!(threadchannel, (Int[], Int[]))
    end

    # This is the same logic as the serial function, but parallelized.
    Threads.@threads for j in eachindex(y)
        rowvals, colvals = take!(threadchannel)
        nzcol = 0
        for i in eachindex(x)
            @inbounds if evaluate(metric, x[i], y[j]) ≤ ( (ε isa Real) ? ε : ε[j] )
                push!(rowvals, i) # push to the thread-specific row array
                nzcol += 1
            end
        end
        append!(colvals, fill(j, (nzcol,)))
        put!(threadchannel, (rowvals, colvals))
    end
    # merge into one array
    finalrows = Int[]
    finalcols = Int[]
    foreach(1:nbuffers) do _
        rowvals, colvals = take!(threadchannel)
        append!(finalrows, rowvals)
        append!(finalcols, colvals)
    end
    nzvals = fill(true, (length(finalrows),))
    return sparse(finalrows, finalcols, nzvals, length(x), length(y))
end

function recurrence_matrix(x::Vector_or_SSSet, metric::Metric, ε, ::Val{true})
    @assert ε isa Real || length(ε) == length(x)
    # We create a Channel for `Array`s of `Array`s, for each thread to have its
    # own array to push to.  This avoids race conditions with
    # multiple threads pushing to the same `Array` (`Array`s are not atomic).
    nbuffers = Threads.nthreads()
    threadchannel = Channel{NTuple{2, Vector{Int}}}(nbuffers) # for rows and columns
    foreach(1:nbuffers) do _
        put!(threadchannel, (Int[], Int[]))
    end

    # This is the same logic as the serial function, but parallelized.
    Threads.@threads for k in partition_indices(length(x))
        rowvals, colvals = take!(threadchannel)
        for j in k
            nzcol = 0
            for i in 1:j
                @inbounds if evaluate(metric, x[i], x[j]) ≤ ( (ε isa Real) ? ε : ε[j] )
                    push!(rowvals, i) # push to the thread-specific row array
                    nzcol += 1
                end
            end
            append!(colvals, fill(j, (nzcol,)))
        end
    end
    # merge into one array
    finalrows = Int[]
    finalcols = Int[]
    foreach(1:nbuffers) do _
        rowvals, colvals = take!(threadchannel)
        append!(finalrows, rowvals)
        append!(finalcols, colvals)
    end
    nzvals = fill(true, (length(finalrows),))
    return Symmetric(sparse(finalrows, finalcols, nzvals, length(x), length(x)), :U)
end
