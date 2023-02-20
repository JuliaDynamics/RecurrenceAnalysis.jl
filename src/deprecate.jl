export @windowed
macro windowed(ex, options...)
    error("""
    The `@windowed` macro is removed. It was some incredibly complicated 250 lines of code
    that offer little to no benefit. Instead of using this macro, write a trivial loop over
    a view of a recurrence matrix. E.g., replace
    ```julia
    rmat = RecurrenceMatrix(...)
    @windowed determinism(rmat, theiler=2, lmin=3) width=1000 step=100
    ```
    with
    ```julia
    width = 1000
    step = 100
    windows = 1:step:(size(rmat, 1)-width)
    map(1:length(windows)) do i
        rmat_view = view(rmat, windows[i]:(windows[i]+width), windows[i]:(windows[i]+width))
        determinism(rmat_view; theiler=2, lmin=3)
    end
    ```
    """)
end

for T in (:WithinRange, :NeighborNumber)
    @eval function RecurrenceMatrix{$T}(x::AbstractMatrix, ε; kwargs...)
        @warn string("`RecurrenceMatrix{", $T, "}(x::AbstractMatrix, ε; kwargs...)` is deprecated, use `RecurrenceMatrix{", $T, "}(Dataset(x), ε; kwargs...)`")
        RecurrenceMatrix{$T}(Dataset(x), ε; kwargs...)
    end

    for call in (:CrossRecurrenceMatrix, :JointRecurrenceMatrix)
        @eval function $call{$T}(x::AbstractMatrix, y, ε; kwargs...)
            @warn string("`", $call, "{", $T, "}(x::AbstractMatrix, y, ε; kwargs...)` is deprecated, use `", $call, "{", $T, "}(Dataset(x), y, ε; kwargs...)`")
            $call{$T}(Dataset(x), y, ε; kwargs...)
        end

        @eval function $call{$T}(x, y::AbstractMatrix, ε; kwargs...)
            @warn string("`", $call, "{", $T, "}(x, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `", $call, "{", $T, "}(x, Dataset(y), ε; kwargs...)`")
            $call{$T}(x, Dataset(y), ε; kwargs...)
        end

        @eval function $call{$T}(x::AbstractMatrix, y::AbstractMatrix, ε; kwargs...)
            @warn string("`", $call, "{", $T, "}(x::AbstractMatrix, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `", $call, "{", $T, "}(Dataset(x), Dataset(y), ε; kwargs...)`")
            $call{$T}(Dataset(x), Dataset(y), ε; kwargs...)
        end
    end
end

function rqa(::Type{NamedTuple}, R; onlydiagonal=false, kwargs...)
    @warn "RQA with `NamedTuple` output is deprecated, and will be removed in later versions"
    # Parse arguments for diagonal and vertical structures
    kw_d = Dict(kwargs)
    haskey(kw_d, :theilerdiag) && (kw_d[:theiler] = kw_d[:theilerdiag])
    haskey(kw_d, :lmindiag) && (kw_d[:lmin] = kw_d[:lmindiag])
    dhist = diagonalhistogram(R; kw_d...)
    rr_d = recurrencerate(R; kw_d...)
    if onlydiagonal
        return (
            RR  = recurrencerate(R; kwargs...),
            DET   = _determinism(dhist, rr_d*_rrdenominator(R; kw_d...)),
            L     = _dl_average(dhist),
            Lmax  = _dl_max(dhist),
            DIV   = 1.0/_dl_max(dhist),
            ENTR  = _dl_entropy(dhist)
        )
   else
        kw_v = Dict(kwargs)
        haskey(kw_v, :theilervert) && (kw_v[:theiler] = kw_v[:theilervert])
        haskey(kw_v, :lminvert) && (kw_v[:lmin] = kw_v[:lminvert])
        vhist, rthist = verticalhistograms(R; kw_v...)
        rr_v = recurrencerate(R; kw_v...)
        return (
            RR  = rr_d,
            TRANS = transitivity(R),
            DET  = _determinism(dhist, rr_d*_rrdenominator(R; kw_v...)),
            L    = _dl_average(dhist),
            Lmax = _dl_max(dhist),
            DIV  = 1.0/_dl_max(dhist),
            ENTR  = _dl_entropy(dhist),
            TREND = trend(R; kw_d...),
            LAM  = _laminarity(vhist, rr_v*_rrdenominator(R; kw_v...)),
            TT   = _vl_average(vhist),
            Vmax = _vl_max(vhist),
            VENTR = _vl_entropy(vhist),
            MRT  = _rt_average(rthist),
            RTE = _rt_entropy(rthist),
            NMPRT = maximum(rthist)
        )
    end
end

function transitivity(R::AbstractRecurrenceMatrix)
    @warn "`transitivity(x::AbstractRecurrenceMatrix)` is deprecated`, use `rna` to analyse network parameters"
    if size(R, 1) ≠ size(R, 2)
        @warn "Computing network transitivity of a non-square adjacency matrix is impossible"
        return NaN
    end
    R² = R.data * R.data
    numerator = zero(eltype(R²))
    for col in axes(R,2)
        rows = view(rowvals(R), nzrange(R,col))
        for r in rows
            numerator += R²[r, col]
        end
    end
    trans = numerator / (sum(R²) - LinearAlgebra.tr(R²))
end

export transitivity


import Base.Meta.parse
const METRICS = Dict(
    "euclidean"=>Euclidean(),
    "max"=>Chebyshev(),
    "inf"=>Chebyshev(),
    "cityblock"=>Cityblock(),
    "manhattan"=>Cityblock(),
    "taxicab"=>Cityblock(),
    "min"=>Cityblock(),
)
getmetric(m::Metric) = m
function getmetric(normtype::AbstractString)
    @warn "Specifying metric with strings is deprecated! "*
    "Use a formal instance of a `Metric` from Distances.jl, e.g., `Euclidean()`."
    normtype = lowercase(normtype)
    !haskey(METRICS, normtype) && error("incorrect norm type. Accepted values are \""
        *join(keys(METRICS),"\", \"", "\" or \"") * "\".")
    METRICS[normtype]
end

################################################################################
# TODO: Deprecations of old recurrence matrix interface
################################################################################
const FAN = NeighborNumber
export WithinRange, NeighborNumber, FAN

# OLD KEYWORD ARGUMENTS in `RecurrenceMatrix`: scale, fixedrate.
function RecurrenceMatrix{FAN}(x, ε; fixedrate = true, metric = Euclidean(), parallel = false)
    @warn "Specifying `RecurrenceMatrix{FAN}` is deprecated! Use `LocalRecurrenceRate`."
    rt = LocalRecurrenceRate(ε)
    return RecurrenceMatrix(x, rt; metric, parallel)
end

function CrossRecurrenceMatrix{FAN}(x, y, ε; fixedrate = true, metric = Euclidean(), parallel = false)
    @warn "Specifying `CrossRecurrenceMatrix{FA}` is deprecated! Use `LocalRecurrenceRate`."
    rt = LocalRecurrenceRate(ε)
    return CrossRecurrenceMatrix(x, y, rt; metric, parallel)
end

function _computescale(scale::Real, args...)
    @warn "specifying `scale` as a number is deprecated because its pointless. "*
    "Use a modified `ε = ε*scale` instead..."
    return scale
end
