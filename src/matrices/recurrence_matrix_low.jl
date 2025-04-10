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

function _full_task_worker!(rowvals, colvals, x::Vector_or_SSSet, y::Vector_or_SSSet, k, metric::Metric, ε::E) where E
    for j in k
        nzcol = 0
        for i in eachindex(x)
            @inbounds if evaluate(metric, x[i], y[j]) ≤ ( (ε isa Real) ? ε : ε[j] )
                push!(rowvals, i) # push to the thread-specific row array
                nzcol += 1
            end
        end
        append!(colvals, fill(j, (nzcol,)))
    end
    return nothing
end

function recurrence_matrix(x::Vector_or_SSSet, y::Vector_or_SSSet, metric::Metric, ε, ::Val{true})
    @assert ε isa Real || length(ε) == length(y)
    # We create an `Array` of `Array`s, for each task to have its
    # own array to push to.  This avoids race conditions with
    # multiple threads pushing to the same `Array` (`Array`s are not atomic).

    # For load balancing reasons, we will use 2 tasks per CPU thread that we have
    # available.
    ntasks = Threads.nthreads()
    chunks = ntuple(ntasks) do ichunk
        # Adapted from JuliaFolds2/ChunkSplitters.jl, MIT license, credit to the authors.
        n_per_chunk, n_remaining = divrem(length(y), ntasks)
        first = 1 + (ichunk - 1) * n_per_chunk + ifelse(ichunk <= n_remaining, ichunk - 1, n_remaining)
        last = (first - 1) + n_per_chunk + ifelse(ichunk <= n_remaining, 1, 0)
        first:last
    end

    rowvals = [Vector{Int}() for _ in 1:ntasks]
    colvals = [Vector{Int}() for _ in 1:ntasks]

    # This is the same logic as the serial function, but parallelized.
    # We create all tasks in this array comprehension so we have them stored,
    # but they are launched at the same time approximately.
    tasks = [
        Threads.@spawn _full_task_worker!($(rowvals[i]), $(colvals[i]), $(x), $(y), $(k), $(metric), $(ε))
        for (i, k) in enumerate(chunks)
    ]
    # The array comprehension above scheduled the tasks...now we have to wait on them to make
    # sure they are all complete before we access the results.
    foreach(fetch, tasks)
    # Now that we know all tasks are complete, we can merge the results.
    finalrows = reduce(vcat, rowvals) # merge into one array
    finalcols = reduce(vcat, colvals) # merge into one array
    nzvals = fill(true, (length(finalrows),))
    return sparse(finalrows, finalcols, nzvals, length(x), length(y))
end

function _lower_triangular_task_worker!(rowvals, colvals, x::Vector_or_SSSet, k, metric::Metric, ε::E) where E
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
    return nothing
end

function recurrence_matrix(x::Vector_or_SSSet, metric::Metric, ε, ::Val{true})
    @assert ε isa Real || length(ε) == length(x)
    # We create an `Array` of `Array`s, for each thread to have its
    # own array to push to.  This avoids race conditions with
    # multiple threads pushing to the same `Array` (`Array`s are not atomic).

    ntasks = Threads.nthreads()

    rowvals = [Vector{Int}() for _ in 1:ntasks]
    colvals = [Vector{Int}() for _ in 1:ntasks]

    # This is the same logic as the serial function, but parallelized.
    tasks = [
        Threads.@spawn _lower_triangular_task_worker!($(rowvals[i]), $(colvals[i]), $(x), $(k), $(metric), $(ε))
        for (i, k) in enumerate(partition_indices(length(x), ntasks))
    ]
    foreach(fetch, tasks)
    finalrows = reduce(vcat, rowvals) # merge into one array
    finalcols = reduce(vcat, colvals) # merge into one array
    nzvals = fill(true, (length(finalrows),))
    return Symmetric(sparse(finalrows, finalcols, nzvals, length(x), length(x)), :U)
end
