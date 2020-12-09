function RecurrenceMatrix(x::AbstractMatrix, ε; kwargs...)
    @warn "`RecurrenceMatrix(x::AbstractMatrix, ε; kwargs...)` is deprecated, use `RecurrenceMatrix(Dataset(x), ε; kwargs...)`"
    RecurrenceMatrix(Dataset(x), ε; kwargs...)
end

for call in (:CrossRecurrenceMatrix, :JointRecurrenceMatrix)
    @eval function ($call)(x::AbstractMatrix, y, ε; kwargs...)
        @warn string("`", $call, "(x::AbstractMatrix, y, ε; kwargs...)` is deprecated, use `", $call, "(Dataset(x), y, ε; kwargs...)`")
        ($call)(Dataset(x), y, ε; kwargs...)
    end
    
    @eval function ($call)(x, y::AbstractMatrix, ε; kwargs...)
        @warn string("`", $call, "(x, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `", $call, "(x, Dataset(y), ε; kwargs...)`")
        ($call)(x, Dataset(y), ε; kwargs...)
    end

    @eval function ($call)(x::AbstractMatrix, y::AbstractMatrix, ε; kwargs...)
        @warn string("`", $call, "(x::AbstractMatrix, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `", $call, "(Dataset(x), Dataset(y), ε; kwargs...)`")
        ($call)(Dataset(x), Dataset(y), ε; kwargs...)
    end
end

function rqa(::Type{NamedTuple}, R; onlydiagonal=false, kwargs...)
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