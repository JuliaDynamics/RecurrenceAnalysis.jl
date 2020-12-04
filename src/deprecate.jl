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

