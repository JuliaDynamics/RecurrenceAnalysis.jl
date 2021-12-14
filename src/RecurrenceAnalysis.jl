module RecurrenceAnalysis

using Distances, Statistics, LinearAlgebra, SparseArrays, DelayEmbeddings, StaticArrays
import Base.Meta.parse

const METRICS = Dict(
    "euclidean"=>Euclidean(),
    "max"=>Chebyshev(),
    "inf"=>Chebyshev(),
    "cityblock"=>Cityblock(),
    "manhattan"=>Cityblock(),
    "taxicab"=>Cityblock(),
    "min"=>Cityblock()
)
const DEFAULT_METRIC = Euclidean()
getmetric(m::Metric) = m
function getmetric(normtype::AbstractString)
    normtype = lowercase(normtype)
    !haskey(METRICS,normtype) && error("incorrect norm type. Accepted values are \""
        *join(keys(METRICS),"\", \"", "\" or \"") * "\".")
    METRICS[normtype]
end


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
       # deprecated:
       transitivity

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
