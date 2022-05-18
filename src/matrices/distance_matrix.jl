################################################################################
# Distance matrix
################################################################################
"""
    distancematrix(x [, y = x], metric = "euclidean")

Create a matrix with the distances between each pair of points of the
time series `x` and `y` using `metric`.

The time series `x` and `y` can be `AbstractDataset`s or vectors or matrices with data points
in rows.
The data point dimensions (or number of columns) must be the same for `x` and `y`.
The returned value is a `n×m` matrix, with `n` being the length (or number of rows)
of `x`, and `m` the length of `y`.

The metric can be identified by a string, or any of the `Metric`s defined in
the [`Distances` package](https://github.com/JuliaStats/Distances.jl).
The list of strings available to define the metric are:

* `"max"` or `"inf"` for the maximum or L∞ norm
  (`Chebyshev()` in the `Distances` package).
* `"euclidean"` for the L2 or Euclidean norm, used by default
  (`Euclidean()` in `Distances`).
* `"manhattan"`, `"cityblock"`, `"taxicab"` or `"min"` for the Manhattan or L1 norm
  (`Cityblock()` in `Distances`).
"""
distancematrix(x, metric::Union{Metric,String}=DEFAULT_METRIC, parallel = size(x)[1] > 500) = _distancematrix(x, getmetric(metric), Val(parallel))

# For 1-dimensional arrays (vectors), the distance does not depend on the metric
distancematrix(x::Vector, y::Vector, metric=DEFAULT_METRIC, parallel = size(x)[1] > 500) = abs.(x .- y')

# If the metric is supplied as a string, get the corresponding Metric from Distances
distancematrix(x, y, metric::String, parallel = size(x)[1] > 500) = distancematrix(x, y, getmetric(metric), parallel)

const MAXDIM = 9
function distancematrix(x::Tx, y::Ty, metric::Metric=DEFAULT_METRIC, parallel = size(x)[1] > 500) where
         {Tx<:Union{AbstractMatrix, AbstractDataset}} where {Ty<:Union{AbstractMatrix, AbstractDataset}}
    sx, sy = size(x), size(y)
    @assert sx[2] == sy[2] """
        The dimensions of the data points in `x` and `y` must be equal!
        Found dim(x)=$(sx[2]), dim(y)=$(sy[2]).
        """

    if sx[2] ≥ MAXDIM && typeof(metric) == Euclidean # Blas optimization
        return _distancematrix(Matrix(x), Matrix(y), metric, Val(parallel))
    elseif Tx <: Matrix && Ty <: Matrix && metric == Chebyshev()
        return _distancematrix(x, y, metric, Val(parallel))
    else
        return _distancematrix(Dataset(x), Dataset(y), metric, Val(parallel))
    end
end

# Core function for Matrices (wrapper of `pairwise` from the Distances package)
_distancematrix(x::AbstractMatrix, y::AbstractMatrix, metric::Metric, ::Val{false}) = pairwise(metric, x', y', dims=2)

# First we define the serial versions of the functions.
# Core function for AbstractDatasets
function _distancematrix(x::AbstractDataset{S,Tx}, y::AbstractDataset{S,Ty},
    metric::Metric, ::Val{false}) where {S, Tx, Ty}

    x = x.data
    y = y.data
    d = zeros(promote_type(Tx,Ty), length(x), length(y))
    for j in eachindex(y), i in eachindex(x)
        @inbounds d[i,j] = evaluate(metric, x[i], y[j])
    end
    return d
end

# Now, we define the parallel versions.

function _distancematrix(x::AbstractDataset{S,Tx}, y::AbstractDataset{S,Ty},
    metric::Metric, ::Val{true}) where {S, Tx, Ty}

    x = x.data
    y = y.data
    d = zeros(promote_type(Tx,Ty), length(x), length(y))
    Threads.@threads for j in eachindex(y)
        for i in eachindex(x)
            @inbounds d[i,j] = evaluate(metric, x[i], y[j])
        end
    end
    return d
end

function _distancematrix(x::Matrix{Tx}, y::Matrix{Ty},
    metric::Metric, ::Val{true}) where {Tx, Ty}

    x = x.data
    y = y.data
    d = zeros(promote_type(Tx,Ty), length(x), length(y))
    Threads.@threads for j in eachindex(y)
        for i in eachindex(x)
            @inbounds d[i,j] = evaluate(metric, x[i, :], y[j, :])
        end
    end
    return d
end


# Now, we define methods for the single trajectory.
# We can speed this up by only calculating the lower triangle,
# which halves the number of computations needed.

# Again, we'll define the serial version first:
function _distancematrix(x::Vector{T}, metric::Metric, ::Val{false}) where T
    d = zeros(T, length(x), length(x))

    for j in eachindex(x)
        for i in 1:j
            @inbounds d[i, j] = abs(x[i] - x[j])
        end
    end

    return Symmetric(d, :U)

end

function _distancematrix(x::AbstractDataset{S, T}, metric::Metric, ::Val{false}) where T where S
    d = zeros(T, length(x), length(x))

    for j in eachindex(x)
        for i in 1:j
            @inbounds d[i, j] = evaluate(metric, x[i], x[j])
        end
    end

    return Symmetric(d, :U)

end

# Now, we define the parallel version.  There's a twist, though.

# The single trajectory case is a little tricky.  If you try it naively,
# using the method we used for the serial computation above, the points are
# partioned out unequally between threads.  This means that performance actually
# **decreases** compared to the full parallel version.  To mitigate this, we need
# to partition the computation equally among all threads.
# The function I've written below does essentially this - given a length,
# it calculates the appropriate partitioning, and returns a list of indices,
# by utilizing the fact that the area is proportional to the square of the height.
# It partitions the "triangle" which needs to be computed into "trapezoids",
# which all have an equal area.
function partition_indices(len)
    indices = Vector{UnitRange{Int}}(undef, Threads.nthreads())
    length = len
    offset = 1

    # Here, we have to iterate in reverse, to get an equal area every time
    for n in Threads.nthreads():-1:1
        partition = round(Int, length / sqrt(n)) # length varies as square root of area
        indices[n] = offset:(partition .+ offset .- 1)
        length -= partition
        offset += partition
    end

    return indices
end

# And now for the actual implementation:
function _distancematrix(x::Vector{T}, metric::Metric, ::Val{true}) where T
    d = zeros(T, length(x), length(x))

    Threads.@threads for k in partition_indices(length(x))
        for j in k
            for i in 1:j
                @inbounds d[i, j] = abs(x[i] - x[j])
            end
        end
    end

    return Symmetric(d, :U)

end

function _distancematrix(x::AbstractDataset{S, T}, metric::Metric, ::Val{true}) where T where S
    d = zeros(T, length(x), length(x))

    Threads.@threads for k in partition_indices(length(x))
        for j in k
            for i in 1:j
                @inbounds d[i, j] = evaluate(metric, x[i], x[j])
            end
        end
    end

    return Symmetric(d, :U)

end
