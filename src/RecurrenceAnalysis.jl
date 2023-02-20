module RecurrenceAnalysis

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end RecurrenceAnalysis

using Reexport
@reexport using StateSpaceSets

Array_or_SSSet = Union{AbstractArray{<:Real}, AbstractStateSpaceSet}
Vector_or_SSSet = Union{AbstractVector{<:Real}, AbstractStateSpaceSet}

using Distances, Statistics, LinearAlgebra, SparseArrays

const DEFAULT_METRIC = Euclidean()

# recurrence_specification.jl
export AbstractRecurrenceType, RecurrenceThreshold, RecurrenceThresholdScaled,
    GlobalRecurrenceRate, LocalRecurrenceRate, recurrence_threshold


export RecurrenceMatrix, CrossRecurrenceMatrix, JointRecurrenceMatrix

export embed,
       reconstruct,
       recurrence_threshold,
       Dataset,
       distancematrix,
       textrecurrenceplot,
       recurrenceplot,
       recurrencerate,
       recurrencestructures,
       dl_average, dl_max, dl_entropy,
       vl_average, vl_max, vl_entropy,
       rt_average, rt_max, rt_entropy,
       determinism,
       divergence,
       trend,
       laminarity,
       trappingtime,
       meanrecurrencetime,
       nmprt,
       rqa,
       sorteddistances,
       skeletonize,
       rna


include("matrices/distance_matrix.jl")
include("matrices/recurrence_specification.jl")
include("matrices/recurrence_matrix_types.jl")
include("matrices/recurrence_matrix_low.jl")
include("matrices/plot.jl")
include("matrices/skeletonization.jl")
include("rqa/histograms.jl")
include("rqa/rqa.jl")
include("rqa/radius.jl")
include("rna/graphs.jl")
include("rna/rna.jl")
# include("rqa/windowed.jl") # this is removed!!!
include("deprecate.jl")

end
