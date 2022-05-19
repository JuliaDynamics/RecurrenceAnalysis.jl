module RecurrenceAnalysis

using Distances, Statistics, LinearAlgebra, SparseArrays, DelayEmbeddings, StaticArrays
import Base.Meta.parse

const DEFAULT_METRIC = Euclidean()


export RecurrenceMatrix, CrossRecurrenceMatrix, JointRecurrenceMatrix,
       AbstractRecurrenceMatrix, WithinRange, NeighborNumber, FAN

export embed,
       reconstruct,
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
       rna,


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
