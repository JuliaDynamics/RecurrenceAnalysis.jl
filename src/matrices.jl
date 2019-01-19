const METRICS = Dict(
    "euclidean"=>Euclidean(),
    "max"=>Chebyshev(),
    "inf"=>Chebyshev(),
    "cityblock"=>Cityblock(),
    "manhattan"=>Cityblock(),
    "taxicab"=>Cityblock(),
    "min"=>Cityblock()
)

const DEFAULT_METRIC = Euclidean()

function getmetric(normtype::AbstractString)
    normtype = lowercase(normtype)
    !haskey(METRICS,normtype) && error("incorrect norm type. Accepted values are \""
        *join(keys(METRICS),"\", \"", "\" or \"") * "\".")
    METRICS[normtype]
end


#### Distance matrix ####

"""
    distancematrix(x [, y = x], metric = "euclidean")

Create a matrix with the distances between each pair of points of the
time series `x` and `y` using `metric`.

The time series `x` and `y` can be `Dataset`s or vectors or matrices with data points
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
distancematrix(x, metric::Union{Metric,String}=DEFAULT_METRIC) = distancematrix(x, x, metric)

# For 1-dimensional arrays (vectors), the distance does not depend on the metric
distancematrix(x::Vector, y::Vector, metric=DEFAULT_METRIC) = abs.(x .- y')

# If the metric is supplied as a string, get the corresponding Metric from Distances
distancematrix(x, y, metric::String) = distancematrix(x, y, getmetric(metric))

const MAXDIM = 9
function distancematrix(x::Tx, y::Ty, metric::Metric=DEFAULT_METRIC) where
         {Tx<:Union{AbstractMatrix, Dataset}} where {Ty<:Union{AbstractMatrix, Dataset}}
    sx, sy = size(x), size(y)
    if sx[2] != sy[2]
        error("the dimensions of `x` and `y` data points must be the equal")
    end
    if sx[2] ≥ MAXDIM && typeof(metric) == Euclidean # Blas optimization
        return _distancematrix(Matrix(x), Matrix(y), metric)
    elseif Tx <: Matrix && Ty <: Matrix && metric == Chebyshev()
        return _distancematrix(x, y, metric)
    else
        return _distancematrix(Dataset(x), Dataset(y), metric)
    end
end

# Core function for Matrices (wrapper of `pairwise` from the Distances package)
_distancematrix(x::AbstractMatrix, y::AbstractMatrix, metric::Metric) =
pairwise(metric, x', y')
# Core function for Datasets
function _distancematrix(x::Dataset{S,Tx}, y::Dataset{S,Ty},
    metric::Metric) where {S, Tx, Ty}

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



#######################
# Type
#######################
abstract type AbstractRecurrenceMatrix end
const ARM = AbstractRecurrenceMatrix
struct RecurrenceMatrix <: AbstractRecurrenceMatrix
    data::SparseMatrixCSC{Bool,Int}
end
struct CrossRecurrenceMatrix <: AbstractRecurrenceMatrix
    data::SparseMatrixCSC{Bool,Int}
end
struct JointRecurrenceMatrix <: AbstractRecurrenceMatrix
    data::SparseMatrixCSC{Bool,Int}
end

function Base.summary(R::AbstractRecurrenceMatrix)
    N = nnz(R.data)
    return "$(nameof(typeof(R))) of size $(size(R.data)) with $N entries"
end

# Propagate used functions:
begin
    extentions = [
        (:Base, (:getindex, :size, :length, :view, :iterate)),
        (:LinearAlgebra, (:diag, :triu, :tril, :issymmetric)),
        (:SparseArrays, (:nnz, :rowvals, :nzrange, :nonzeros))
    ]
    for (M, fs) in extentions
        for f in fs
            @eval $M.$(f)(x::ARM, args...) = $(f)(x.data, args...)
        end
    end
end

for operator in [:(==), :(!=)]
    @eval Base.$operator(x::ARM, y::ARM) = $operator(x.data, y.data)
end

LinearAlgebra.issymmetric(::RecurrenceMatrix) = true
# column values in sparse matrix (parallel to rowvals)
function colvals(x::SparseMatrixCSC)
    cv = zeros(Int,nnz(x))
    @inbounds for c=1:size(x,2)
        cv[nzrange(x,c)] .= c
    end
    cv
end
colvals(x::ARM) = colvals(x.data)

@deprecate recurrencematrix RecurrenceMatrix
@deprecate crossrecurrencematrix CrossRecurrenceMatrix
@deprecate jointrecurrencematrix JointRecurrenceMatrix


"""
    RecurrenceMatrix(x, ε; kwargs...)

Create a recurrence matrix from trajectory `x`.

## Description

The recurrence matrix is a numeric representation of a "recurrence plot" [1, 2],
in the form of a sparse square matrix of Boolean values.

`x` must be a `Vector` or a `Dataset` or a `Matrix` with data points in rows
(possibly representing and embedded phase space; see [`embed`](@ref)).
If `d(x[i], x[j]) ≤ ε` (with `d` the distance function),
then the cell `(i, j)` of the matrix will have a `true`
value. The criteria to evaluate distances between data points are defined
by the following keyword arguments:

* `scale=1` : a function of the distance matrix (see [`distancematrix`](@ref)),
  or a fixed number, used to scale the value of `ε`. Typical choices are
  `maximum` or `mean`, such that the threshold `ε` is defined as a ratio of the
  maximum or the mean distance between data points, respectively (using
  `mean` or `maximum` calls specialized versions that are faster than the naive
  approach).  Use `1` to keep the distances unscaled (default).
* `fixedrate::Bool=false` : a flag that indicates if `ε` should be
  taken as a target fixed recurrence rate (see [`recurrencerate`](@ref)).
  If `fixedrate` is set to `true`, `ε` must be a value between 0 and 1,
  and `scale` is ignored.
* `metric="euclidean"` : metric of the distances, either `Metric` or a string,
   as in [`distancematrix`](@ref).

See also: [`CrossRecurrenceMatrix`](@ref), [`JointRecurrenceMatrix`](@ref) and
use [`recurrenceplot`](@ref) to turn the result of these functions into a plottable format.

## References
[1] : N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems",
*Phys. Reports 438*(5-6), 237-329 (2007).

[2] : N. Marwan & C.L. Webber, "Mathematical and computational foundations of
recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence
Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).
"""
function RecurrenceMatrix(x, ε; kwargs...)
    m = recurrence_matrix(x, x, ε; kwargs...)
    return RecurrenceMatrix(m)
end


#### Cross recurrence matrix ####

"""
    CrossRecurrenceMatrix(x, y, ε; kwargs...)

Create a cross recurrence matrix from trajectories `x` and `y`.

The cross recurrence matrix is a bivariate extension of the recurrence matrix.
For the time series `x`, `y`, of length `n` and `m`, respectively, it is a
sparse `n×m` matrix of Boolean values, such that if `d(x[i], y[j]) ≤ ε`,
then the cell `(i, j)` of the matrix will have a `true` value.

See [`RecurrenceMatrix`](@ref) for details, references and keywords.
See also: [`JointRecurrenceMatrix`](@ref).
"""
function CrossRecurrenceMatrix(x, y, ε; kwargs...)
    m = recurrence_matrix(x, y, ε; kwargs...)
    return CrossRecurrenceMatrix(m)
end

function recurrence_matrix(x, y, ε; scale=1, fixedrate=false, metric=DEFAULT_METRIC)
    # Check fixed recurrence rate - ε must be within (0, 1)
    if fixedrate
        sfun = (m) -> quantile(m[:], ε)
        return recurrence_matrix(x, y, 1; scale=sfun, fixedrate=false, metric=metric)
    else
        scale_value = _computescale(scale, x, y, metric)
        spm = _recurrence_matrix(x, y, ε*scale_value, metric)
        return spm
    end
end

# If `scale` is a function, compute the numeric value of the scale based on the
# distance matrix; otherwise return the value of `scale` itself
_computescale(scale::Function, x, y, metric) = scale(distancematrix(x, y, metric))
_computescale(scale::Real, args...) = scale
# specific methods to avoid `distancematrix`
function _computescale(scale::typeof(maximum), x::T, y::T, metric::Metric) where {T}
    maxvalue = zero(eltype(x))
    @inbounds for xi in x, yj in y
        newvalue = evaluate(metric, xi, yj)
        (newvalue > maxvalue) && (maxvalue = newvalue)
    end
    return maxvalue
end
function _computescale(scale::typeof(mean), x, y, metric::Metric)
    meanvalue = 0.0
    @inbounds for xi in x, yj in y
        meanvalue += evaluate(metric, xi, yj)
    end
    return meanvalue/(length(x)*length(y))
end


# Internal methods to calculate the matrix:
# If the metric is supplied as a string, get the corresponding Metric from Distances
_recurrence_matrix(x, y, ε, metric::String="max") =
_recurrence_matrix(x, y, ε, getmetric(metric))

# Convert the inputs to Datasets (better performance in all cases)
function _recurrence_matrix(x::AbstractVecOrMat, y::AbstractVecOrMat,
                                ε, metric::Metric=DEFAULT_METRIC)
    return _recurrence_matrix(Dataset(x), Dataset(y), ε, metric)
end

# Core function
function _recurrence_matrix(x::Dataset, y::Dataset, ε, metric::Metric)
    x = x.data
    y = y.data
    rowvals = Vector{Int}()
    colvals = Vector{Int}()
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
    JointRecurrenceMatrix(x, y, ε; kwargs...)

Create a joint recurrence matrix from `x` and `y`.

The joint recurrence matrix considers the recurrences of the trajectories
of `x` and `y` separately, and looks for points where both recur
simultaneously. It is calculated by the element-wise multiplication
of the recurrence matrices of `x` and `y`. If `x` and `y` are of different
length, the recurrences are only calculated until the length of the shortest one.

See [`RecurrenceMatrix`](@ref) for details, references and keywords.
See also: [`CrossRecurrenceMatrix`](@ref).
"""
function JointRecurrenceMatrix(x, y, ε; kwargs...)
    n = min(size(x,1), size(y,1))
    if n == size(x,1) && n == size(y,1)
        rm1 = RecurrenceMatrix(x, ε, kwargs...)
        rm2 = RecurrenceMatrix(y, ε, kwargs...)
    else
        rm1 = RecurrenceMatrix(x[1:n,:], ε, kwargs...)
        rm2 = RecurrenceMatrix(y[1:n,:], ε, kwargs...)
    end
    return JointRecurrenceMatrix(rm1.data .* rm2.data)
end

#######################
# Pretty printing
#######################
function oldshow(io::IO, R::AbstractRecurrenceMatrix)
    s = sprint(io -> show(IOContext(io, :limit=>true), MIME"text/plain"(), R.data))
    s = split(s, '\n')[2:end]
    s = [replace(line, "=  true"=>"", count=1) for line in s]
    s = join(s, '\n')
    tos = summary(R)*"\n"*s
    println(io, tos)
end

Base.show(io::IO, R::AbstractRecurrenceMatrix) = println(io, summary(R))

"""
    scatterdata(R::ARM) -> xs, ys
Transform the data of a recurrence matrix to scatter data `xs, ys`.
"""
function scatterdata(R::ARM) # TODO: scan every ÷100 elements depending on size?
   rows = rowvals(R)
   is = zeros(Int, nnz(R))
   js = zeros(Int, nnz(R))
   k = 1;
   m, n = size(R)
   for j = 1:n
     for i in nzrange(R, j)
        is[k] = j
        js[k] = rows[i]
        k += 1
     end
   end
   return is, js
end

using UnicodePlots
export unicode_arm

unicode_arm(R::ARM; kwargs...) = unicode_arm(stdout, R::ARM; kwargs...)

"""
    unicode_arm([io,] R; minh = 25, maxh = 0.5, ascii, kwargs...) -> u

Obtain a text-scatterplot representation of a recurrence matrix `R` to be displayed in
`io` (by default `stdout`). The matrix spans at minimum `minh` rows and at maximum
`maxh*displaysize(io)[1]` (i.e. by default half the display).
As we always try to plot in equal aspect ratio, if the width of the plot is even less,
the minimum height is dictated by the width.

The keyword `ascii::Bool` can ensure that all elements of the plot are ASCII characters
(`true`, useful when e.g. wanting to print to a file) or Unicode (`false`).

The rest of the `kwargs` are propagated into `UnicodePlots.scatterplot`.

Internally this function calls `RecurrenceAnalysis.scatterdata` to transform
a recurence matrix to scatter data.
"""
function unicode_arm(io::IO, R::ARM; minh = 25, maxh = 0.5, ascii = nothing, kwargs...)
    @assert maxh ≤ 1
    h, w = displaysize(io)
    h = max(minh, round(Int, maxh * h)) # make matrix as long as half the screen (but not too short)
    s = 2.4 # scale that visually brings width and height to equal aspect ratio
    if w < round(Int, h*s) # ensure equal aspect ratio
        h = round(Int, w/s)
    else
        w = round(Int, h*s)
    end

    is, js = scatterdata(R)
    n, m = size(R)

    if ascii == true
        asciidef = (border = :ascii, canvas = DotCanvas)
    elseif ascii == false
        asciidef = (border = :solid, canvas = BrailleCanvas)
        #default:
    elseif ascii == nothing
        # TODO: always use ascii on Jupyter
        if (isdefined(Main, :Juno) && Main.Juno.isactive()) || !Sys.iswindows()
            asciidef = (border = :solid, canvas = BrailleCanvas)
        else
            asciidef = (border = :ascii, canvas = DotCanvas)
        end
    end

    # TODO: Change limits to tuples instead of arrays
    UnicodePlots.scatterplot(
        is, js; xlim = [1, n], ylim = [1, m], title = summary(R), labels = false,
        color = :cyan, width = w, height = h, asciidef..., kwargs...
    )
end

function Base.show(io::IO, ::MIME"text/plain", R::ARM)
    a = unicode_arm(io, R)
    show(io, a)
end
