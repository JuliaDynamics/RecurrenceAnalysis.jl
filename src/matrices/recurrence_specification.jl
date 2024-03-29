################################################################################
# Definition of ways to find a recurrence
################################################################################
"""
    AbstractRecurrenceType
Supertype of all recurrence specification types. Instances of subtypes are given to
[`RecurrenceMatrix`](@ref) and similar constructors to specify recurrences.
Use [`recurrence_threshold`](@ref) to extract the numeric distance threshold.

Possible subtypes are:

* `RecurrenceThreshold(ε::Real)`: Recurrences are defined as any point with distance
  `≤ ε` from the referrence point.
* `RecurrenceThresholdScaled(ratio::Real, scale::Function)`: Here `scale` is a function
  of the distance matrix `dm` (see [`distancematrix`](@ref)) that is used
  to scale the value of the recurrence threshold `ε` so that `ε = ratio*scale(dm)`.
  After the new `ε` is obtained, the method works just like the `RecurrenceThreshold`.
  Specialized versions are employed if `scale` is `mean` or `maximum`.
* `GlobalRecurrenceRate(q::Real)`: Here the number of total recurrence rate over the whole
  matrix is specified to be a quantile `q ∈ (0,1)` of the [`distancematrix`](@ref).
  In practice this yields (approximately) a ratio `q` of recurrences out of the total
  `Nx * Ny` for input trajectories `x, y`.
* `LocalRecurrenceRate(r::Real)`: The recurrence threhsold here is point-dependent. It is
  defined so that each point of `x` has a fixed number of `k = r*N` neighbors, with ratio
  `r` out of the total possible `N`. Equivalently, this means that each column of the
  recurrence matrix will have exactly `k` true entries. Notice that `LocalRecurrenceRate`
  does not guarantee that the resulting recurrence matrix will be symmetric.
"""
abstract type AbstractRecurrenceType end

"""
    RecurrenceThreshold(ε::Real)
Recurrences are defined as any point with distance `≤ ε` from the referrence point.
See [`AbstractRecurrenceType`](@ref) for more.
"""
struct RecurrenceThreshold{T<:Real} <: AbstractRecurrenceType
    ε::T
end

"""
    RecurrenceThresholdScaled(ratio::Real, scale::Function)
Recurrences are defined as any point with distance `≤ d` from the referrence point,
where `d` is a scaled ratio (specified by `ratio, scale`) of the distance matrix.
See [`AbstractRecurrenceType`](@ref) for more.
"""
struct RecurrenceThresholdScaled{T<:Real, S<:Function} <: AbstractRecurrenceType
    ratio::T
    scale::S
end

"""
    GlobalRecurrenceRate(rate::Real)
Recurrences are defined as a constant global recurrence rate.
See [`AbstractRecurrenceType`](@ref) for more.
"""
struct GlobalRecurrenceRate{T<:Real} <: AbstractRecurrenceType
    rate::T
end


"""
    LocalRecurrenceRate(rate::Real)
Recurrences are defined as a constant local recurrence rate.
See [`AbstractRecurrenceType`](@ref) for more.
"""
struct LocalRecurrenceRate{T<:Real} <: AbstractRecurrenceType
    rate::T
end


################################################################################
# Resolving the recurrence threshold and/or scaling
################################################################################
# TODO: We can optimize stuff here; for methods that compute the distance matrix,
# we don't have to compute distances again when creating the actual recurrence
# matrices........

"""
    recurrence_threshold(rt::AbstractRecurrenceType, x [, y] [, metric]) → ε
Return the calculated distance threshold `ε` for `rt`. The output is real, unless
`rt isa LocalRecurrenceRate`, where `ε isa Vector`.
"""
recurrence_threshold(rt, x::Array_or_SSSet, metric::Metric) =
    recurrence_threshold(rt, x, x, metric)
recurrence_threshold(rt, x::Array_or_SSSet, y::Array_or_SSSet) =
    recurrence_threshold(rt, x, y, Euclidean())
recurrence_threshold(rt, x::Array_or_SSSet) =
    recurrence_threshold(rt, x, x, Euclidean())

recurrence_threshold(rt::Real, x, y, metric) = rt
recurrence_threshold(rt::RecurrenceThreshold, x, y, metric) = rt.ε

function recurrence_threshold(rt::RecurrenceThresholdScaled, x, y, metric)
    scale_value = scale_of_distmatrix(rt.scale, x, y, metric)
    return rt.ratio*scale_value
end

function recurrence_threshold(rt::GlobalRecurrenceRate, x, y, metric)
    # We leverage the code of the scaled version to get the global recurrence rate
    scale = (m) -> quantile(vec(m), rt.rate)
    scale_value = scale_of_distmatrix(scale, x, y, metric)
    return scale_value
end

function recurrence_threshold(rt::LocalRecurrenceRate, x, y, metric)
    rate = rt.rate
    @assert 0 < rate < 1 "Recurrence rate must be ∈ (0, 1)."
    thresholds = zeros(length(y))
    d = distancematrix(x, y, metric)
    if x === y
        rate += 1/length(x) # because of recurrences guaranteed in the diagonal
    end
    for i in axes(d, 2)
        thresholds[i] = quantile(view(d, :, i), rate)
    end
    return thresholds
end

function scale_of_distmatrix(scale::Function, x, y, metric)
    return scale(distancematrix(x, y, metric))
end

# specific methods to avoid making a `distancematrix`
function scale_of_distmatrix(::typeof(maximum), x, y, metric::Metric)
    maxvalue = zero(eltype(x))
    if x===y
        @inbounds for i in 1:length(x)-1, j=(i+1):length(x)
            newvalue = evaluate(metric, x[i], y[j])
            (newvalue > maxvalue) && (maxvalue = newvalue)
        end
    else
        @inbounds for xi in x, yj in y
            newvalue = evaluate(metric, xi, yj)
            (newvalue > maxvalue) && (maxvalue = newvalue)
        end
    end
    return maxvalue
end
function scale_of_distmatrix(::typeof(mean), x, y, metric::Metric)
    meanvalue = 0.0
    if x===y
        @inbounds for i in 1:length(x)-1, j=(i+1):length(x)
            meanvalue += evaluate(metric, x[i], y[j])
        end
        denominator = length(x)*(length(x)-1)/2
    else
        @inbounds for xi in x, yj in y
            meanvalue += evaluate(metric, xi, yj)
        end
        denominator = Float64(length(x)*length(y))
    end
    return meanvalue/denominator
end
