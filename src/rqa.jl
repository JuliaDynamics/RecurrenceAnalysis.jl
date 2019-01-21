# Recurrence parameters as defined by Marwan et al. (2007)

# Recurrence rate

"""
    recurrencerate(x[; theiler])

Calculate the recurrence rate of the recurrence matrix `x`.

## Description

The recurrence rate is calculated as:

```math
RR = \\frac{1}{S} \\sum_{i,j} x[i,j]
```

where *S* is the size of `x` or the region of `x` with potential recurrent points.
There is not a unique definition of that denominator, which is defined as the
full size of the matrix in many sources (e.g. [1]), whereas
in others it is adjusted to remove the points of the LOI when they are
excluded from the count [2,3].

For matrices of type `RecurrenceMatrix` or `JointRecurrenceMatrix`, where the
points around the central diagonal are usually excluded, the denominator is
adjusted to the size of the matrix outside the Theiler window
(by default equal to the LOI, and adjustable with the keyword argument `theiler`;
see [`rqa`](@ref) for details). For matrices of type `CrossRecurrenceMatrix`,
where normally all points are analyzed, the denominator is always the full
size of the matrix, regardless of the Theiler window that might be defined
(none by default). 

**Hint**: to reproduce the calculations done following the formulas that use
the full size of the matrix in the denominator, use
`CrossRecurrenceMatrix(s,s,ε)` to define the recurrence matrix, instead of
`RecurrenceMatrix(s,ε), setting `theiler=1` (or `theiler=n` in general) to
explicitly exclude the LOI or other diagonals around it.

## References

[1] : N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems",
*Phys. Reports 438*(5-6), 237-329 (2007).

[2] : C.L. Webber & J.P. Zbilut, "Recurrence Quantification Analysis of Nonlinear
Dynamical Systems", in: Riley MA & Van Orden GC, Tutorials in Contemporary
Nonlinear Methods for the Behavioral Sciences, 26-94 (2005).
URL: https://www.nsf.gov/pubs/2005/nsf05057/nmbs/nmbs.pdf

[3] : N. Marwan & C.L. Webber, "Mathematical and computational foundations of
recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence
Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).
"""
function recurrencerate(x::ARM; theiler::Integer=deftheiler(x), kwargs...)::Float64
    (theiler < 0) && throw(ErrorException(
        "Theiler window length must be greater than or equal to 0"))
    if theiler == 0
        return nnz(x)/length(x)
    end
    diags_remove = -(theiler-1):(theiler-1)
    theiler_nz = 0
    for d in diags_remove
        theiler_nz += nnz(diag(x,d))
    end
    return (nnz(x) - theiler_nz)/_rrdenominator(x; theiler=theiler, kwargs...)
end

# Calculate the denominator for the recurrence rate
_rrdenominator(x::ARM; theiler=0, kwargs...) = length(x)

function _rrdenominator(x::M; theiler=0, kwargs...) where
    M<:Union{RecurrenceMatrix,JointRecurrenceMatrix}
    
    (theiler == 0) && (return length(x))
    k = size(x,1) - theiler
    return k*(k+1)
end

# Generic parameters of histograms: mean, average and entropy

macro histogram_params(keyword, description, hist_fun)
    combined_descriptions = Dict(:average => "average of the $(description)s",
                                 :max     => "longest $(description)",
                                 :entropy => "entropy of the $(description)s")
    function_bodies = Dict(
        :average => quote
                (hist==[0]) && return 0.0
                points = (1:length(hist)) .* hist
                return sum(points)/sum(hist)
            end,
        :max => quote
                (hist==[0]) && return 0
                return length(hist)
            end,
        :entropy => quote
                (hist==[0]) && return NaN
                prob_bins = hist ./ sum(hist)
                prob_bins = prob_bins[findall(!iszero, prob_bins)]
                return -sum(prob_bins .* log.(prob_bins))
            end
        )

    ret = quote end
    for (name, param) in combined_descriptions
        fname = Symbol("$(keyword)_$(name)")
        _fname = Symbol("_$(keyword)_$(name)")
        fbody = function_bodies[name]
        doc = """
                $fname(x[; lmin=2, theiler])

            Calculate the $(param) contained in the recurrence matrix `x`,
            ruling out the lines shorter than `lmin` (2 by default) and all the
            points inside the Theiler window (see [`rqa`](@ref) for the
            default values and usage of the keyword argument `theiler`).
            """
        push!(ret.args, quote
            @doc $doc ->
            $fname(x::ARM; kwargs...) = $_fname($hist_fun(x; kwargs...))

            function $_fname(hist::Vector{<:Integer})
                $fbody
            end
        end)
    end
    return esc(ret)
end


# 1. Based on diagonal lines

@histogram_params dl "diagonal line" diagonalhistogram

"""
    determinism(x[; lmin=2, theiler])

Calculate the determinism of the recurrence matrix `x`, ruling out
the diagonal lines shorter than `lmin` (2 by default) and all the
points inside the Theiler window (see [`rqa`](@ref) for the
default values and usage of the keyword argument `theiler`).
"""
function determinism(x::ARM; kwargs...)
    npoints = recurrencerate(x; kwargs...)*_rrdenominator(x; kwargs...)
    return _determinism(diagonalhistogram(x; kwargs...), npoints)
end

function _determinism(diag_hist::Vector{<:Integer}, npoints)::Float64
    (diag_hist==[0]) && return 0.0
    diag_points = (1:length(diag_hist)) .* diag_hist
    return sum(diag_points)/npoints
end


"""
    divergence(x[; theiler])

Calculate the divergence of the recurrence matrix `x`
(actually the inverse of [`dl_max`](@ref)).
"""
divergence(x::ARM; kwargs...) = ( 1.0/dl_max(x; kwargs...) )


"""
    trend(x[; border=10, theiler])

Calculate the trend of recurrences in the recurrence matrix `x`.

The 10 outermost diagonals (counting from the corners of the matrix)
are excluded by default to avoid "border effects". Use the keyword argument
`border` to define a different number of excluded lines, and `theiler`
to exclude the central points inside the Theiler window
(see [`rqa`](@ref) for details).
"""
trend(x::ARM; theiler=deftheiler(x), kwargs...) =
    _trend(tau_recurrence(x); theiler=theiler, kwargs...)

function tau_recurrence(x::ARM)
    n = minimum(size(x))
    rr_τ1 = [count(!iszero, diag(x,d))/length(diag(x,d)) for d in (0:n-1)]
    rr_τ2 = [count(!iszero, diag(x,-d))/length(diag(x,-d)) for d in (0:n-1)]
    return (rr_τ1 + rr_τ2)
end

function _trend(rr_τ::Vector; theiler=1, border=10, kwargs...)::Float64
    nmax = length(rr_τ)
    a = 1+theiler
    b = nmax-border
    w = collect(a:b) .- b/2
    (w'*(rr_τ[a:b] .- mean(rr_τ[a:b])) ./ (w'*w))[1]
end

# Number of l-length sequences, based on diagonals
# countsequences(x::ARM; kwargs...) = countsequences(diagonalhistogram(x; kwargs...))
#
# function countsequences(diag_hist::Vector; lmin=2, kwargs...)
#     overlap = (1:length(diag_hist))' .- (lmin+1)
#     overlap[overlap .< 0] = 0
#     overlap * diag_hist
# end


# 2. Based on vertical lines

@histogram_params vl "vertical line" vl_histogram


"""
    laminarity(x[; lmin=2, theiler])

Calculate the laminarity of the recurrence matrix `x`, ruling out the
lines shorter than `lmin` (2 by default) and all the
points inside the Theiler window (see [`rqa`](@ref) for the
default values and usage of the keyword argument `theiler`).
"""
function laminarity(x::ARM; kwargs...)
    npoints = recurrencerate(x)*_rrdenominator(x; kwargs...)
    return _laminarity(verticalhistograms(x; kwargs...)[1], npoints)
end

_laminarity(vert_hist::Vector{<:Integer}, npoints) = _determinism(vert_hist, npoints)

"""
    trappingtime(x[; lmin=2, theiler])

Calculate the trapping time of the recurrence matrix `x`, ruling out the
lines shorter than `lmin` (2 by default) and all the
points inside the Theiler window (see [`rqa`](@ref) for the
default values and usage of the keyword argument `theiler`).

The trapping time is the average of the vertical line structures and thus equal
to [`vl_average`](@ref).
"""
trappingtime(x::ARM; kwargs...) = vl_average(x; kwargs...)


# 3. Based on recurrence times

@histogram_params rt "recurrence time" rt_histogram

"""
    meanrecurrencetime(x[; lmin=2, theiler])

Calculate the mean recurrence time of the recurrence matrix `x`, ruling out the
lines shorter than `lmin` (2 by default) and all the
points inside the Theiler window (see [`rqa`](@ref) for the
default values and usage of the keyword argument `theiler`).

Equivalent to [`rt_average`](@ref).
"""
meanrecurrencetime(x::ARM; kwargs...) = rt_average(x; kwargs...)


"""
    nmprt(x[; lmin=2, theiler])

Calculate the number of the most probable recurrence time (NMPRT), ruling out the
lines shorter than `lmin` (2 by default) and all the
points inside the Theiler window (see [`rqa`](@ref) for the
default values and usage of the keyword argument `theiler`).
"""
nmprt(x::ARM; kwargs) = maximum(verticalhistograms(x; kwargs...)[2])


# 4. All in one

"""
    rqa(x; kwargs...)

Calculate all RQA parameters of a recurrence matrix. See the functions
referred to below for the definition of
the different parameters and the default values of the arguments.
Using this function is much more efficient than calling all individual functions
one by one.

## Return
The returned value is a dictionary with the following keys:

* "RR": recurrence rate (see [`recurrencerate`](@ref))
* "DET": determinsm (see [`determinism`](@ref))
* "L": average length of diagonal structures (see [`dl_average`](@ref))
* "Lmax": maximum length of diagonal structures (see [`dl_max`](@ref))
* "DIV": divergence (see [`divergence`](@ref))
* "ENTR": entropy of diagonal structures (see [`dl_entropy`](@ref))
* "TREND": trend of recurrences (see [`trend`](@ref))
* "LAM": laminarity (see [`laminarity`](@ref))
* "TT": trapping time (see [`trappingtime`](@ref))
* "Vmax": maximum length of vertical structures (see [`vl_max`](@ref))
* "VENTR": entropy of vertical structures (see [`vl_entropy`](@ref))
* "MRT": mean recurrence time (see [`meanrecurrencetime`](@ref))
* "RTE" recurrence time entropy (see [`rt_entropy`](@ref))
* "NMPRT": number of the most probable recurrence time (see [`nmprt`](@ref))

In the case of empty histograms (e.g. no existing vertical lines
less than the keyword `lminvert`) the average and maximum values
("L", "Lmax", "TT", "Vmax", "MRT")
are returned as `0.0` but their respective entropies ("ENTR", "VENTR", "RTE")
are returned as `NaN`.

## Keyword Arguments

Standard keyword arguments are the ones accepted by the functions listed below,
i.e. `theiler`, `lmin`, and `border`:

* `theiler` is used to define a "Theiler window" around the central diagonal or
  "line of identity" (LOI): a region of points that are excluded in the calculation
  of RQA parameters, in order to rule out self-recurrences and apparent recurrences
  for smooth or high resolution data. The LOI is excluded by default for matrices
  of the types `RecurrenceMatrix` or `JointRecurrenceMatrix`, but it is included
  for matrices of the type `CrossRecurrenceMatrix`. `theiler=0` means that the
  whole matrix is scanned for lines. `theiler=1` means that the LOI is excluded.
  In general, `theiler=n` means that the `n` central diagonals are excluded
  (at both sides of the LOI, i.e. actually `2n-1` diagonals are excluded).
  
* `lmin` is used to define the minimum line length in the parameters that
  describe the distributions of diagonal or vertical lines (it is set as 2 by
  default).
  
* `border` is used to avoid border effects in the calculation of "TREND"
  (cf. [`trend`](@ref)).

In addition `theilerdiag`, `lmindiag` may be used to
declare specific values that override the values of `theiler` and `lmin` in the
calculation of parameters related to diagonal structures. Likewise, `theilervert` and
`lminvert` can be used for the calculation of parameters related to vertical
structures.

The keyword argument `onlydiagonal` (`false` by default) can be set to `true`
in order to restrict the analysis to the recurrence rate and the parameters related
to diagonal structures ("RR", "DET", "L", "Lmax", "DIV" and "ENTR"), which makes
this function slightly faster.
"""
function rqa(x; onlydiagonal=false, kwargs...)
    # Parse arguments for diagonal and vertical structures
    kw_d = Dict(kwargs)
    haskey(kw_d, :theilerdiag) && (kw_d[:theiler] = kw_d[:theilerdiag])
    haskey(kw_d, :lmindiag) && (kw_d[:lmin] = kw_d[:lmindiag])
    dhist = diagonalhistogram(x; kw_d...)
    rr_d = recurrencerate(x; kw_d...)
    if onlydiagonal
        return Dict("RR"  => recurrencerate(x; kwargs...),
        "DET"  => _determinism(dhist, rr_d*_rrdenominator(x; kw_d...)),
        "L"    => _dl_average(dhist),
        "Lmax" => _dl_max(dhist),
        "DIV"  => 1.0/_dl_max(dhist),
        "ENTR"  => _dl_entropy(dhist)
        )
   else
        kw_v = Dict(kwargs)
        haskey(kw_v, :theilervert) && (kw_v[:theiler] = kw_v[:theilervert])
        haskey(kw_v, :lminvert) && (kw_v[:lmin] = kw_v[:lminvert])
        vhist, rthist = verticalhistograms(x; kw_v...)
        rr_v = recurrencerate(x; kw_v...)
        return Dict("RR"  => recurrencerate(x; kwargs...),
            "DET"  => _determinism(dhist, rr_d*_rrdenominator(x; kw_v...)),
            "L"    => _dl_average(dhist),
            "Lmax" => _dl_max(dhist),
            "DIV"  => 1.0/_dl_max(dhist),
            "ENTR"  => _dl_entropy(dhist),
            "TREND" => trend(x; kw_d...),
            "LAM"  => _laminarity(vhist, rr_v*_rrdenominator(x; kw_v...)),
            "TT"   => _vl_average(vhist),
            "Vmax" => _vl_max(vhist),
            "VENTR" => _vl_entropy(vhist),
            "MRT"  => _rt_average(rthist),
            "RTE" => _rt_entropy(rthist),
            "NMPRT" => maximum(rthist)
        )
    end
end
