@deprecate(RecurrenceMatrix(x::AbstractMatrix, ε; kwargs...), 
           RecurrenceMatrix(Dataset(x), ε; kwargs...))

for call in (:CrossRecurrenceMatrix, :JointRecurrenceMatrix)
    eval(quote
        @deprecate(($call)(x::AbstractMatrix, y, ε; kwargs...), 
                   ($call)(Dataset(x), y, ε; kwargs...))
        
        @deprecate(($call)(x, y::AbstractMatrix, ε; kwargs...), 
                   ($call)(x, Dataset(y), ε; kwargs...))
        
        @deprecate(($call)(x::AbstractMatrix, y::AbstractMatrix, ε; kwargs...),
                   ($call)(Dataset(x), Dataset(y), ε; kwargs...))
    end)
end