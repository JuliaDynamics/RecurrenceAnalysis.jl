for T in (:FixedRange, :FixedAmount)
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
