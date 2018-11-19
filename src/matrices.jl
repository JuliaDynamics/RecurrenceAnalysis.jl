# get distance metric of the Distance package
function getmetric(normtype::AbstractString)
    normtype = lowercase(normtype)
    metrics = Dict(
        "euclidean"=>Euclidean(),
        "max"=>Chebyshev(),
        "inf"=>Chebyshev(),
        "cityblock"=>Cityblock(),
        "manhattan"=>Cityblock(),
        "taxicab"=>Cityblock()
        )
    !haskey(metrics,normtype) && error("incorrect norm type. Accepted values are \""
        *join(keys(metrics),"\", \"", "\" or \"") * "\".")
    metrics[normtype]
end


#### Distance matrix ####

"""
    distancematrix(x[, y, metric])
    
Create a matrix with the distances between each pair of points of the
time series `x`, or between each point of `x` and `y`, with the distance type
defined by `metric`.

The time series `x` and `y` can be `Dataset`s or matrices with data points in rows.
The data point dimensions (or number of columns) must be the same for `x` and `y`.
The returned value is a `n×m` matrix, with `n` being the length (or number of rows)
of `x`, and `m` the length of `y`. If `y` is not provided, the result is a square
`n×n` matrix.

The metric can be identified by a string, or any of the `Metric`s defined in
the [Distances package](https://github.com/JuliaStats/Distances.jl).
The list of strings available to define the metric are [1]:

* `"max"` or `"inf"` for the maximum or L∞ norm 
  (`Chebyshev()` in the `Distances` package, used by default).
* `"euclidean"` for the L2 or Euclidean norm
  (`Euclidean()` in `Distances`).
* `"manhattan"`, `"cityblock"`, `"taxicab"` or `"min"` for the Manhattan or L1 norm
  (`Cityblock()` in `Distances`).
  

# References
[1] : N. Marwan *et al.*. "Recurrence plots for the analysis of complex systems",
*Phys. Reports 438*(5-6), 237-329 (2007).
"""
distancematrix(x, metric::Union{Metric,String}=Chebyshev()) = distancematrix(x, x, metric)

# For 1-dimensional arrays (vectors), the distance does not depend on the metric
distancematrix(x::Vector, y::Vector, metric=Chebyshev()) = abs.(x .- y')

# If the metric is supplied as a string, get the corresponding Metric from Distances
distancematrix(x, y, metric::String) = distancematrix(x, y, getmetric(metric))

# Choose the metod based on Matrices or Datasets, depending on _maxdimension:
# for smaller dimensions, Datasets are more efficient
const _maxdimension = 10
function distancematrix(x::Tx, y::Ty, metric::Metric=Chebyshev()) where
         {Tx<:Union{AbstractMatrix, Dataset}} where {Ty<:Union{AbstractMatrix, Dataset}}
    sx, sy = size(x), size(y)
    if sx[2] != sy[2]
        error("the dimensions of `x` and `y` data points must be the equal")
    end
    if sx[2] > _maxdimension
        return _distancematrix(Matrix(x), Matrix(y), metric)
    else
        return _distancematrix(Dataset(x), Dataset(y), metric)
    end
end

# Core function for Matrices (wrapper of `pairwise` from the Distances package)
_distancematrix(x::AbstractMatrix, y::AbstractMatrix, metric::Metric) = pairwise(metric, x', y')
# Core function for Datasets
function _distancematrix(x::Dx, y::Dy, metric::Metric) where 
         {Dx<:Dataset{S,Tx}} where {Dy<:Dataset{S,Ty}} where {S} where {Tx} where {Ty}
    x = x.data
    y = y.data
    d = zeros(promote_type(Tx,Ty), length(x), length(y))
    for j in 1:length(y)
        for i in 1:length(x)
            @inbounds d[i,j] = evaluate(metric, x[i], y[j])
        end
    end
    return d
end


#### Recurrence matrix ####
# Defined as a wrapper of crossrecurrencematrix

"""
    recurrencematrix(x, ε; <keyword arguments>)
    
Create a recurrence matrix from an embedded time series.

# Description

The "recurrence matrix" is a numeric representation of a "recurrence plot" [1, 2],
in the form of a sparse square matrix of Boolean values.

`x` must be `Dataset` or a Vector or Matrix with data points in rows
(possibly representing and embedded phase, space; see [`embed`](@ref)). 
If ∥`x[i]` – `x[j]`∥ ≤ ε, then the cell `(i, j)` of the matrix will have a `true`
value. The criteria to evaluate distances between data points are defined
by the following keyword arguments:

* `scale` : a function of the distance matrix (see [`distancematrix`](@ref)),
  or a fixed number, used to scale the value of `ε`. Typical choices are
  `maximum` or `mean`, such that the threshold `ε` is defined as a ratio of the
  maximum or the mean distance between data points, respectively.
  Use `1` to keep the distances unscaled (default).
* `fixedrate::Bool=false` : a flag that indicates if `ε` should be
  taken as a target fixed recurrence rate (see [`recurrencerate`](@ref)).
  If `fixedrate` is set to `true`, `ε` must be a value between 0 and 1,
  and `scale` is ignored.
* `metric` : metric of the distances, as in [`distancematrix`](@ref).

See also: [`crossrecurrencematrix`](@ref), [`jointrecurrencematrix`](@ref)

# References
[1] : N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems",
*Phys. Reports 438*(5-6), 237-329 (2007).

[2] : N. Marwan & C.L. Webber, "Mathematical and computational foundations of
recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence
Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).
"""
recurrencematrix(x, ε; kwargs...) = crossrecurrencematrix(x, x, ε; kwargs...) 


#### Cross recurrence matrix ####

"""
    crossrecurrencematrix(x, y, ε; <keyword arguments>)
    
Create a cross recurrence matrix from the time series `x` and `y`.

The cross recurrence matrix is a bivariate extension of the recurrence matrix [1, 2].
For the time series `x`, `y`, of length `n` and `m`, respectively, it is a
sparse `n×m` matrix of Boolean values, such that if ∥`x[i]` – `x[j]`∥ ≤ ε,
then the cell `(i, j)` of the matrix will have a `true` value.

See [`recurrencematrix`](@ref) for details and a description of the arguments.

See also: [`jointrecurrencematrix`](@ref)

# References
[1] : N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems",
*Phys. Reports 438*(5-6), 237-329 (2007).

[2] : N. Marwan & C.L. Webber, "Mathematical and computational foundations of
recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence
Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).
"""
function crossrecurrencematrix(x, y, ε; scale=1, fixedrate=false, metric=Chebyshev())
    # Check fixed recurrence rate - ε must be within (0, 1)
    if fixedrate
        sfun = (m) -> quantile(m[:], ε)
        return crossrecurrencematrix(x, y, 1; scale=sfun, fixedrate=false, metric=metric)
    else
        scale_value = _computescale(scale, x, y, metric)
        return _crossrecurrencematrix(x, y, ε*scale_value, metric)
    end
end

# If `scale` is a function, compute the numeric value of the scale based on the
# distance matrix; otherwise return the value of `scale` itself
_computescale(scale::Function, x, y, metric) = scale(distancematrix(x, y, metric))
_computescale(scale::Real, args...) = scale

# Internal methods to calculate the matrix:
# If the metric is supplied as a string, get the corresponding Metric from Distances 
_crossrecurrencematrix(x, y, ε, metric::String="max") = _crossrecurrencematrix(x, y, ε, getmetric(metric))

# Convert the inputs to Datasets (better performance in all cases)
function _crossrecurrencematrix(x::AbstractVecOrMat, y::AbstractVecOrMat, ε, metric::Metric=Chebyshev())
    return _crossrecurrencematrix(Dataset(x), Dataset(y), ε, metric)
end

# Core function
function _crossrecurrencematrix(x::Dataset, y::Dataset, ε, metric::Metric) 
    x = x.data
    y = y.data
    rowvals = Vector{Int64}()
    colvals = Vector{Int64}()
    for j in 1:length(y)
        nzcol = 0
        for i in 1:length(x)
            @inbounds if evaluate(metric, x[i], y[j]) ≤ ε
                push!(rowvals,i)
                nzcol += 1
            end
        end
        append!(colvals, fill(j, (nzcol,)))
    end
    nzvals = fill(true, (length(rowvals),))
    return sparse(rowvals, colvals, nzvals, length(x), length(y))
end


#### Joint recurrence matrix ####

"""
    jointrecurrencematrix(x, y, ε; <keyword arguments>)
    
Create a joint recurrence matrix from the time series `x` and `y`.

The joint recurrence matrix considers the recurrences of the trajectories
of `x` and `y` separately, and looks for points where both recur
simultaneously [1, 2]. It is calculated by the element-wise multiplication
of the recurrence matrices of `x` and `y`. If `x` and `y` are of different
length, the recurrences are only calculated until the length of the shortest one.

See [`recurrencematrix`](@ref) for details and a description of the arguments.

See also: [`crossrecurrencematrix`](@ref)

# References
[1] : N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems",
*Phys. Reports 438*(5-6), 237-329 (2007).

[2] : N. Marwan & C.L. Webber, "Mathematical and computational foundations of
recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence
Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).
"""
function jointrecurrencematrix(x, y, ε; kwargs...)
    n = min(size(x,1), size(y,1))
    rm1 = recurrencematrix(x[1:n,:], ε, kwargs...)
    rm2 = recurrencematrix(y[1:n,:], ε, kwargs...)
    return rm1 .* rm2
end

