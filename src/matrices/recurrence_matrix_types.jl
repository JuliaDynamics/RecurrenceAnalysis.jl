################################################################################
# AbstractRecurrenceMatrix type hierarchy
################################################################################
# TODO: `RT` is `AbstractRecurrenceType`, but due to deprecations, it is allowed
# to be anything at the moment
abstract type AbstractRecurrenceMatrix{RT} <: AbstractMatrix{Bool} end
const ARM = AbstractRecurrenceMatrix

# The recurrence type is included as a field in all matrix types,
# because we use it in the pretty printing.
struct RecurrenceMatrix{RT} <: AbstractRecurrenceMatrix{RT}
    data::SparseMatrixCSC{Bool,Int}
    recurrence_type::RT
end
struct CrossRecurrenceMatrix{RT} <: AbstractRecurrenceMatrix{RT}
    data::SparseMatrixCSC{Bool,Int}
    recurrence_type::RT
end
struct JointRecurrenceMatrix{RT} <: AbstractRecurrenceMatrix{RT}
    data::SparseMatrixCSC{Bool,Int}
    recurrence_type::RT
end

# Pretty printing:
function Base.summary(R::AbstractRecurrenceMatrix)
    N = nnz(R.data)
    return "$(size(R.data)) $(nameof(typeof(R))) "*
    "with $N recurrences of type $(nameof(typeof(R.recurrence_type)))."
end
Base.show(io::IO, R::AbstractRecurrenceMatrix) = println(io, summary(R))

function Base.show(io::IO, ::MIME"text/plain", R::ARM)
    # `recurrenceplot` comes from `plot.jl` file.
    a = recurrenceplot(io, R)
    show(io, a)
end


################################################################################
# Extensions of standard library
################################################################################
begin
    extentions = [
        (:Base, (:getindex, :size, :length, :view, :iterate,
            :eachindex, :axes, :CartesianIndices)),
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

# Ambiguity resolution
LinearAlgebra.diag(X::ARM, k::Integer) = diag(X.data, k)

# Special symmetry cases
LinearAlgebra.issymmetric(::RecurrenceMatrix{X}) where
    {X <: Union{RecurrenceThreshold,RecurrenceThresholdScaled,GlobalRecurrenceRate}} = true
LinearAlgebra.issymmetric(::JointRecurrenceMatrix{X}) where
    {X <: Union{RecurrenceThreshold,RecurrenceThresholdScaled,GlobalRecurrenceRate}} = true


# column values in sparse matrix (parallel to rowvals)
function colvals(x::SparseMatrixCSC)
    cv = zeros(Int, nnz(x))
    @inbounds for c in axes(x,2)
        cv[nzrange(x,c)] .= c
    end
    cv
end
colvals(x::ARM) = colvals(x.data)

# Convert to matrices
Base.Array(R::ARM) = Matrix(R.data)
Base.Matrix(R::ARM) = Matrix(R.data)
SparseArrays.SparseMatrixCSC(R::ARM) = SparseMatrixCSC(R.data)

################################################################################
# Documentation strings and dispatch to `recurrence_matrix`
################################################################################
"""
    RecurrenceMatrix(x, rthres; metric = Euclidean(), parallel::Bool)

Create a recurrence matrix from timeseries or trajectory `x`
and with recurrence threshold `rthres`.
`x` is either a [`StateSpaceSet`](@ref) for multivariate data
or an `AbstractVector{<:Real}` for timeseries.

The variable `rthres` defines how recurrences are estimated.
It can be any subtype of [`AbstractRecurrenceType`](@ref),
and different types can specify recurrences differently.
Alternatively, `rthres` can be a real number, which then becomes
an instance of [`RecurrenceThreshold`](@ref).

The keyword `metric`, if given, must be any subtype of `Metric` from
[Distances.jl](https://github.com/JuliaStats/Distances.jl)
and defines the metric used to calculate distances for recurrences.
By default the Euclidean metric is used, typical alternatives are `Chebyshev(), Cityblock()`.

The keyword `parallel` decides if the comptutation should be done in parallel using threads.
Defaults to `length(x) > 500 && Threads.nthreads() > 1`.

## Description

A (cross-)recurrence matrix is a way to quantify *recurrences* that occur in a trajectory.
A recurrence happens when a trajectory visits the same neighborhood on the state space that
it was at some previous time.

The recurrence matrix is a numeric representation of a recurrence plot,
described in detail in [^Marwan2007] and [^Marwan2015]. It represents a
a sparse square matrix of Boolean values that quantifies recurrences in the trajectory,
i.e., points where the trajectory returns close to itself.
Given trajectories `x, y`, and asumming `ε isa Real`, the matrix is defined as:

```julia
R[i,j] = metric(x[i], y[i]) ≤ ε ? true : false
```
with the `metric` being the distance function.
The difference between a `RecurrenceMatrix` and a [`CrossRecurrenceMatrix`](@ref)
is that in the first case `x === y`.

Objects of type `<:AbstractRecurrenceMatrix` are displayed as a [`recurrenceplot`](@ref).

See also: [`CrossRecurrenceMatrix`](@ref), [`JointRecurrenceMatrix`](@ref) and
use [`recurrenceplot`](@ref) to turn the result of these functions into a plottable format.

[^Marwan2007]:
    N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems",
    [Phys. Reports 438*(5-6), 237-329 (2007)](https://doi.org/10.1016/j.physrep.2006.11.001)

[^Marwan2015]:
    N. Marwan & C.L. Webber, *Recurrence Quantification Analysis. Theory and Best Practices*
    [Springer (2015)](https://link.springer.com/book/10.1007/978-3-319-07155-8)
"""
function RecurrenceMatrix(x, ε::Real;
    # DEPRECATED keywords. TODO: Remove them in next stable release.
    scale = nothing, fixedrate = nothing,
    # VALID keywords (propagated)
    kwargs...)
    if !isnothing(scale)
        @warn "Providing keyword `scale` is deprecated. Use `RecurrenceThresholdScaled`."
        rt = RecurrenceThresholdScaled(ε, scale)
    elseif !isnothing(fixedrate)
        @warn "Providing keyword `fixedrate` is deprecated. Use `GlobalRecurrenceRate`."
        rt = GlobalRecurrenceRate(ε)
    elseif ε isa Real
        rt = RecurrenceThreshold(ε)
    else
        rt = ε
    end
    return RecurrenceMatrix(x, rt; kwargs...)
end

function RecurrenceMatrix(x, rt::AbstractRecurrenceType; metric = Euclidean(),
        parallel::Bool = length(x) > 500 && Threads.nthreads() > 1
    )
    metric = getmetric(metric) # TODO: Remove this in major update.
    ε = recurrence_threshold(rt, x, metric)
    m = recurrence_matrix(x, metric, ε, Val(parallel))
    return RecurrenceMatrix{typeof(rt)}(m, rt)
end


"""
    CrossRecurrenceMatrix(x, y, rthres; kwargs...)

Create a cross recurrence matrix from trajectories `x` and `y`.
See [`RecurrenceMatrix`](@ref) for possible value for `rthres` and `kwargs`.

The cross recurrence matrix is a bivariate extension of the recurrence matrix.
For the time series `x`, `y`, of length `n` and `m`, respectively, it is a
sparse `n×m` matrix of Boolean values.

Note that cross recurrence matrices are generally not symmetric irrespectively of `rthres`.
"""
function CrossRecurrenceMatrix(x, y, ε;
        # DEPRECATED keywords. TODO: Remove them in next stable release.
        scale = nothing, fixedrate = nothing,
        # Normal keywords
        metric = Euclidean(), parallel::Bool = length(x) > 500 && Threads.nthreads() > 1
    )
    if ε isa AbstractRecurrenceType
        rt = ε
    elseif !isnothing(scale)
        @warn "Providing keyword `scale` is deprecated. Use `RecurrenceThresholdScaled`."
        rt = RecurrenceThresholdScaled(ε, scale)
    elseif !isnothing(fixedrate)
        @warn "Providing keyword `fixedrate` is deprecated. Use `GlobalRecurrenceRate`."
        rt = GlobalRecurrenceRate(ε)
    elseif ε isa Real
        rt = RecurrenceThreshold(ε)
    else
        throw(ArgumentError("Unknown specification of recurrence type."))
    end
    metric = getmetric(metric) # TODO: Remove this in major update.
    ε = recurrence_threshold(rt, x, y, metric)
    m = recurrence_matrix(x, y, metric, ε, Val(parallel))
    return CrossRecurrenceMatrix{typeof(rt)}(m, rt)
end

"""
    JointRecurrenceMatrix(x, y, rthres; kwargs...)

Create a joint recurrence matrix from trajectories `x` and `y`.
See [`RecurrenceMatrix`](@ref) for possible values for `rthres` and `kwargs`.

The joint recurrence matrix considers the recurrences of the trajectories
of `x` and `y` separately, and looks for points where both recur
simultaneously. It is calculated by the element-wise multiplication
of the recurrence matrices of `x` and `y`. If `x` and `y` are of different
length, the recurrences are only calculated until the length of the shortest one.

See [`RecurrenceMatrix`](@ref) for details, references and keywords.
"""
function JointRecurrenceMatrix(x, y, ε; kwargs...)
    m, n = length(x), length(y)
    if m > n
        x = x[1:n]
    elseif n > m
        y = y[1:m]
    end
    rm1 = RecurrenceMatrix(x, ε; kwargs...)
    rm2 = RecurrenceMatrix(y, ε; kwargs...)
    r = ε isa Real ? RecurrenceThreshold(ε) : ε # to be sure we have recurrence type
    return JointRecurrenceMatrix{typeof(r)}(rm1.data .* rm2.data, r)
end

"""
    JointRecurrenceMatrix(R1::AbstractRecurrenceMatrix, R2::AbstractRecurrenceMatrix)

Equivalent with `R1 .* R2`.
"""
function JointRecurrenceMatrix(
        R1::AbstractRecurrenceMatrix{RT}, R2::AbstractRecurrenceMatrix{RT}
    ) where {RT}
    R3 = R1.data .* R2.data
    return JointRecurrenceMatrix{RT}(R3, R1.rt)
end

# Extend methods that just take in directly a boolean matrix
for type in (:RecurrenceMatrix, :CrossRecurrenceMatrix, :JointRecurrenceMatrix)
    @eval begin
        function $(type)(m::AbstractMatrix{Bool})
            x = SparseMatrixCSC(m)
            return $(type){RecurrenceThreshold}(x, RecurrenceThreshold(0))
        end
    end
end
