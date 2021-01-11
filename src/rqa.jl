# Recurrence quantifaction analysis measures

###########################################################################################
# 0. Recurrence rate based
###########################################################################################
"""
    recurrencerate(R[; theiler])

Calculate the recurrence rate of the recurrence matrix `R`.

## Description

The recurrence rate is calculated as:

```math
RR = \\frac{1}{S} \\sum R
```

where ``S`` is the size of `R` or the region of `R` with potential recurrent points.
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

*Hint*: to reproduce the calculations done following the formulas that use
the full size of the matrix in the denominator, use
`CrossRecurrenceMatrix(s,s,ε)` to define the recurrence matrix, instead of
`RecurrenceMatrix(s,ε)`, setting `theiler=1` (or `theiler=n` in general) to
explicitly exclude the LOI or other diagonals around it.

## References

[1] : N. Marwan *et al.*, "Recurrence plots for the analysis of complex systems",
*Phys. Reports 438*(5-6), 237-329 (2007).
[DOI:10.1016/j.physrep.2006.11.001](https://doi.org/10.1016/j.physrep.2006.11.001)

[2] : C.L. Webber & J.P. Zbilut, "Recurrence Quantification Analysis of Nonlinear
Dynamical Systems", in: Riley MA & Van Orden GC, Tutorials in Contemporary
Nonlinear Methods for the Behavioral Sciences, 26-94 (2005).
URL: https://www.nsf.gov/pubs/2005/nsf05057/nmbs/nmbs.pdf

[3] : N. Marwan & C.L. Webber, "Mathematical and computational foundations of
recurrence quantifications", in: Webber, C.L. & N. Marwan (eds.), *Recurrence
Quantification Analysis. Theory and Best Practices*, Springer, pp. 3-43 (2015).
"""
function recurrencerate(R::ARM; theiler::Integer=deftheiler(R), kwargs...)::Float64
    (theiler < 0) && throw(ErrorException(
        "Theiler window length must be greater than or equal to 0"))
    if theiler == 0
        return nnz(R)/length(R)
    end
    diags_remove = -(theiler-1):(theiler-1)
    theiler_nz = 0
    for d in diags_remove
        theiler_nz += nnz(diag(R,d))
    end
    return (nnz(R) - theiler_nz)/_rrdenominator(R; theiler=theiler, kwargs...)
end

# Calculate the denominator for the recurrence rate
_rrdenominator(R::ARM; theiler=0, kwargs...) = length(R)

function _rrdenominator(R::M; theiler=0, kwargs...) where
    M<:Union{RecurrenceMatrix,JointRecurrenceMatrix}

    (theiler == 0) && (return length(R))
    k = size(R,1) - theiler
    return k*(k+1)
end

###########################################################################################
# 0. Histograms
###########################################################################################
macro histogram_params(keyword, description, hist_fun)
    combined_descriptions = Dict(:average => "average of the $(description)s",
                                 :max     => "longest $(description)",
                                 :entropy => "Shannon entropy of the $(description)s")
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
                $fname(R[; lmin=2, theiler])

            Calculate the $(param) contained in the recurrence matrix `R`,
            ruling out the lines shorter than `lmin` (2 by default) and all the
            points inside the Theiler window (see [`rqa`](@ref) for the
            default values and usage of the keyword argument `theiler`).
            """
        push!(ret.args, quote
            @doc $doc ->
            $fname(R::ARM; kwargs...) = $_fname($hist_fun(R; kwargs...))

            function $_fname(hist::Vector{<:Integer})
                $fbody
            end
        end)
    end
    return esc(ret)
end

###########################################################################################
# 1. Based on diagonal lines
###########################################################################################
@histogram_params dl "diagonal line" diagonalhistogram

"""
    determinism(R[; lmin=2, theiler])

Calculate the determinism of the recurrence matrix `R`:

## Description

The determinism is calculated as:

```math
DET = \\frac{\\sum_{l=lmin}{l P(l)}}{\\sum_{l=1}{l P(l)}} =
\\frac{\\sum_{l=lmin}{l P(l)}}{\\sum R}
```

where ``l`` stands for the lengths of diagonal lines in the matrix, and ``P(l)``
is the number of lines of length equal to ``l``.

`lmin` is set to 2 by default, and this calculation rules out all the
points inside the Theiler window (see [`rqa`](@ref) for the
default values and usage of the keyword argument `theiler`).
"""
function determinism(R::ARM; kwargs...)
    npoints = recurrencerate(R; kwargs...)*_rrdenominator(R; kwargs...)
    return _determinism(diagonalhistogram(R; kwargs...), npoints)
end

function _determinism(diag_hist::Vector{<:Integer}, npoints)::Float64
    (diag_hist==[0]) && return 0.0
    diag_points = (1:length(diag_hist)) .* diag_hist
    return sum(diag_points)/npoints
end


"""
    divergence(R[; theiler])

Calculate the divergence of the recurrence matrix `R`
(actually the inverse of [`dl_max`](@ref)).
"""
divergence(R::ARM; kwargs...) = ( 1.0/dl_max(R; kwargs...) )


"""
    trend(R[; border=10, theiler])

Calculate the trend of recurrences in the recurrence matrix `R`.

## Description

The trend is the slope of the linear regression that relates the density of
recurrent points in the diagonals parallel to the LOI and the distance between
those diagonals and the LOI. It quantifies the degree of system stationarity,
such that in recurrence plots where points "fade away" from the central diagonal,
the trend will have a negative value.

It is calculated as:

```math
TREND = 10^3\\frac{\\sum_{d=\\tau}^{\\tilde{N}}\\delta[d]\\left(RR[d]-\\langle RR[d]\\rangle\\right)}{\\sum_{d=\\tau}^{\\tilde{N}}\\delta[d]^2}
```

where ``RR[d]`` is the local recurrence rate of the diagonal ``d``,
``\\delta[d]`` is a balanced measure of the distance between that diagonal and the LOI,
``\\tau`` is the Theiler window (number of central diagonals that are excluded), and
``\\tilde{N}`` is the number of the outmost diagonal that is included.

This parameter is expressed in units of variation recurrence rate every
1000 data points, hence the factor ``10^3`` in the formula [1].

The 10 outermost diagonals (counting from the corners of the matrix)
are excluded by default to avoid "border effects". Use the keyword argument
`border` to define a different number of excluded lines, and `theiler`
to define the size of the Theiler window (see [`rqa`](@ref) for details).

*Note*: In rectangular cross-recurrence plots (i.e. when the time series that
originate them are not of the same length), the limits of the formula for TREND
are not clearly defined. For the sake of consistency, this function limits the
calculations to the biggest square matrix that contains the LOI.

## References
[1] C.L. Webber & J.P. Zbilut, "Recurrence Quantification Analysis of Nonlinear
Dynamical Systems", in: Riley MA & Van Orden GC, *Tutorials in Contemporary
Nonlinear Methods for the Behavioral Sciences*, 2005, 26-94.
https://www.nsf.gov/pubs/2005/nsf05057/nmbs/nmbs.pdf
"""
trend(R::ARM; theiler=deftheiler(R), kwargs...) =
    _trend(tau_recurrence(R); theiler=theiler, kwargs...)

function tau_recurrence(R::ARM)
    n = minimum(size(R))
    rv = rowvals(R)
    rr_τ = zeros(n)
    @inbounds for col=1:n
        for i in nzrange(R, col)
            if (r=rv[i]) ≤ n
                d = abs(r-col)
                if d==0
                    rr_τ[1] += 1.0/(n-d)
                else
                    rr_τ[d+1] += 0.5/(n-d)
                end
            end
        end
    end
    return rr_τ
end

function _trend(rr_τ::Vector; theiler=1, border=10, kwargs...)::Float64
    nmax = length(rr_τ)
    a = 1+theiler
    b = nmax-border
    numerator = denominator = 0.0
    mean_rr = mean(@view rr_τ[a:b])
    for d = a:b
        δ = d - b/2
        numerator += δ*(rr_τ[d] - mean_rr)
        denominator += δ*δ
    end
    return 1000.0*numerator/denominator
end

# Number of l-length sequences, based on diagonals
# countsequences(R::ARM; kwargs...) = countsequences(diagonalhistogram(R; kwargs...))
#
# function countsequences(diag_hist::Vector; lmin=2, kwargs...)
#     overlap = (1:length(diag_hist))' .- (lmin+1)
#     overlap[overlap .< 0] = 0
#     overlap * diag_hist
# end

###########################################################################################
# 2. Based on vertical lines
###########################################################################################

@histogram_params vl "vertical line" vl_histogram


"""
    laminarity(R[; lmin=2, theiler])

Calculate the laminarity of the recurrence matrix `R`.

## Description

The laminarity is calculated as:

```math
LAM = \\frac{\\sum_{v=lmin}{v P(l)}}{\\sum_{v=1}{v P(v)}} =
\\frac{\\sum_{v=lmin}{v P(l)}}{\\sum R}
```

where ``v`` stands for the lengths of vertical lines in the matrix, and ``P(v)``
is the number of lines of length equal to ``v``.

`lmin` is set to 2 by default, and this calculation rules out all the
points inside the Theiler window (see [`rqa`](@ref) for the
default values and usage of the keyword argument `theiler`).
"""
function laminarity(R::ARM; kwargs...)
    npoints = recurrencerate(R)*_rrdenominator(R; kwargs...)
    return _laminarity(verticalhistograms(R; kwargs...)[1], npoints)
end

_laminarity(vert_hist::Vector{<:Integer}, npoints) = _determinism(vert_hist, npoints)

"""
    trappingtime(R[; lmin=2, theiler])

Calculate the trapping time of the recurrence matrix `R`, ruling out the
lines shorter than `lmin` (2 by default) and all the
points inside the Theiler window (see [`rqa`](@ref) for the
default values and usage of the keyword argument `theiler`).

The trapping time is the average of the vertical line structures and thus equal
to [`vl_average`](@ref).
"""
trappingtime(R::ARM; kwargs...) = vl_average(R; kwargs...)

###########################################################################################
# 3. Based on recurrence times
###########################################################################################


@histogram_params rt "recurrence time" rt_histogram

"""
    meanrecurrencetime(R[; lmin=2, theiler])

Calculate the mean recurrence time of the recurrence matrix `R`, ruling out the
lines shorter than `lmin` (2 by default) and all the
points inside the Theiler window (see [`rqa`](@ref) for the
default values and usage of the keyword argument `theiler`).

Equivalent to [`rt_average`](@ref).
"""
meanrecurrencetime(R::ARM; kwargs...) = rt_average(R; kwargs...)


"""
    nmprt(R[; lmin=2, theiler])

Calculate the number of the most probable recurrence time (NMPRT), ruling out the
lines shorter than `lmin` (2 by default) and all the
points inside the Theiler window (see [`rqa`](@ref) for the
default values and usage of the keyword argument `theiler`).

This number indicates how many times the system has recurred using the recurrence
time that appears most frequently, i.e it is the maximum value of the histogram
of recurrence times [1].

## References

[1] : E.J. Ngamga *et al.* "Recurrence analysis of strange nonchaotic dynamics",
*Physical Review E*, 75(3), 036222(1-8) (2007)
[DOI:10.1103/physreve.75.036222](https://doi.org/10.1103/physreve.75.036222)

"""
nmprt(R::ARM; kwargs) = maximum(verticalhistograms(R; kwargs...)[2])

###########################################################################################
# 4. All in one
###########################################################################################

# Transient type
struct RQA
    data::Dict
end

function Base.getproperty(result::RQA, name::Symbol)
    if name === :data
        return getfield(result, :data)
    end
    @warn "x.$name is deprecated for results of `rqa`; use x[:$name] instead"
    getindex(result.data, name)
end

for f in (:getindex, :get, :length, :keys, :values, :collect, :iterate)
    @eval Base.$f(result::RQA, args...; kwargs...) = $f(result.data, args...; kwargs...)
end

Dict(rqa::RQA) = rqa.data

function Base.show(io::IO, mime::MIME"text/plain", result::RQA)
    print(io, "RQA parameters in ")
    show(io, mime, result.data)
end

"""
    rqa(R; kwargs...)

Calculate all RQA parameters of a recurrence matrix `R`. See the functions
referred to below for the definition of
the different parameters and the default values of the arguments.
Using this function is much more efficient than calling all individual functions
one by one.

## Return
The returned value contains the following entries,
which can be retrieved as from a dictionary (e.g. `results[:RR]`, etc.):

* `:RR`: recurrence rate (see [`recurrencerate`](@ref))
* `:DET`: determinsm (see [`determinism`](@ref))
* `:L`: average length of diagonal structures (see [`dl_average`](@ref))
* `:Lmax`: maximum length of diagonal structures (see [`dl_max`](@ref))
* `:DIV`: divergence (see [`divergence`](@ref))
* `:ENTR`: entropy of diagonal structures (see [`dl_entropy`](@ref))
* `:TREND`: trend of recurrences (see [`trend`](@ref))
* `:LAM`: laminarity (see [`laminarity`](@ref))
* `:TT`: trapping time (see [`trappingtime`](@ref))
* `:Vmax`: maximum length of vertical structures (see [`vl_max`](@ref))
* `:VENTR`: entropy of vertical structures (see [`vl_entropy`](@ref))
* `:MRT`: mean recurrence time (see [`meanrecurrencetime`](@ref))
* `:RTE` recurrence time entropy (see [`rt_entropy`](@ref))
* `:NMPRT`: number of the most probable recurrence time (see [`nmprt`](@ref))

All the parameters returned by `rqa` are `Float64` numbers,
even for parameters like `:Lmax`, `:Vmax` or `:NMPRT` which are integer values.
In the case of empty histograms (e.g. no existing vertical lines
less than the keyword `lminvert`) the average and maximum values
(`:L`, `:Lmax`, `:TT`, `:Vmax`, `:MRT`)
are returned as `0.0` but their respective entropies (`:ENTR`, `:VENTR`, `:RTE`)
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

* `border` is used to avoid border effects in the calculation of `:TREND`
  (cf. [`trend`](@ref)).

In addition `theilerdiag`, `lmindiag` may be used to
declare specific values that override the values of `theiler` and `lmin` in the
calculation of parameters related to diagonal structures. Likewise, `theilervert` and
`lminvert` can be used for the calculation of parameters related to vertical
structures.

The keyword argument `onlydiagonal` (`false` by default) can be set to `true`
in order to restrict the analysis to the recurrence rate and the parameters related
to diagonal structures (`:RR`, `:DET`, `:L`, `:Lmax`, `:DIV` and `:ENTR`), which makes
this function slightly faster.

## Transitional note on the returned type

In older versions, the `rqa` function returned a `NamedTuple`,
and in future versions it is planned to return a `Dict` instead.
In both cases, the results can be indexed with square brackets and `Symbol` keys,
as `result[:RR]`, `result[:DET]`, etc.
However, named tuples can also be indexed with "dot syntax", e.g. `result.RR`,
whereas this will not be possible with dictionaries, and there are other
differences in the indexing and iteration of those two types.

In order to facilitate the transition between versions, this function currently
returns a `RQA` object that essentially works as a dictionary, but can
also be indexed with the dot syntax (logging a deprecation warning).
The returned type can also be specified as a first argument of `rqa`
in order to replicate the output of different versions:

* `rqa(NamedTuple, R...)` to obtain the output of the older version (as in 1.3).
* `rqa(Dict, R...)` to obtain the output of the planned future version.
* `rqa(RQA, R...)` to obtain the default current output (same as `rqa(R...)`)

"""
rqa(R; kwargs...) = rqa(RQA, R; kwargs...)

function rqa(::Type{RQA}, R; kwargs...)
    rqa_dict = rqa(Dict, R; kwargs...)
    RQA(rqa_dict)
end

function rqa(::Type{Dict}, R; onlydiagonal=false, kwargs...)
    # Parse arguments for diagonal and vertical structures
    kw_d = Dict(kwargs)
    haskey(kw_d, :theilerdiag) && (kw_d[:theiler] = kw_d[:theilerdiag])
    haskey(kw_d, :lmindiag) && (kw_d[:lmin] = kw_d[:lmindiag])
    dhist = diagonalhistogram(R; kw_d...)
    rr_d = recurrencerate(R; kw_d...)
    if onlydiagonal
        return Dict{Symbol, Float64}(
            :RR    => recurrencerate(R; kwargs...),
            :DET   => _determinism(dhist, rr_d*_rrdenominator(R; kw_d...)),
            :L     => _dl_average(dhist),
            :Lmax  => _dl_max(dhist),
            :DIV   => 1.0/_dl_max(dhist),
            :ENTR  => _dl_entropy(dhist)
        )
   else
        kw_v = Dict(kwargs)
        haskey(kw_v, :theilervert) && (kw_v[:theiler] = kw_v[:theilervert])
        haskey(kw_v, :lminvert) && (kw_v[:lmin] = kw_v[:lminvert])
        vhist, rthist = verticalhistograms(R; kw_v...)
        rr_v = recurrencerate(R; kw_v...)
        rqa_dict = Dict{Symbol, Float64}(
            :RR    => rr_d,
            :DET   => _determinism(dhist, rr_d*_rrdenominator(R; kw_v...)),
            :L     => _dl_average(dhist),
            :Lmax  => _dl_max(dhist),
            :DIV   => 1.0/_dl_max(dhist),
            :ENTR  => _dl_entropy(dhist),
            :TREND => trend(R; kw_d...),
            :LAM   => _laminarity(vhist, rr_v*_rrdenominator(R; kw_v...)),
            :TT    => _vl_average(vhist),
            :Vmax  => _vl_max(vhist),
            :VENTR => _vl_entropy(vhist),
            :MRT   => _rt_average(rthist),
            :RTE   => _rt_entropy(rthist),
            :NMPRT => maximum(rthist)
        )
        return rqa_dict
    end
end
