module RecurrenceAnalysis

using Reexport
@reexport using StateSpaceSets

using Distances, Statistics, LinearAlgebra, SparseArrays

const DEFAULT_METRIC = Euclidean()

export RecurrenceMatrix, CrossRecurrenceMatrix, JointRecurrenceMatrix
export RecurrenceThreshold, RecurrenceThresholdScaled,
    GlobalRecurrenceRate, LocalRecurrenceRate

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
       @windowed,
       rna


include("matrices/distance_matrix.jl")
include("matrices/matrices.jl")
include("matrices/plot.jl")
include("matrices/skeletonization.jl")
include("rqa/histograms.jl")
include("rqa/rqa.jl")
include("rqa/radius.jl")
include("rna/graphs.jl")
include("rna/rna.jl")
include("rqa/windowed.jl")
include("deprecate.jl")

end
