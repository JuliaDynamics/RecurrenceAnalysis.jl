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
       AbstractRecurrenceMatrix, WithinRange, NeighborNumber, FAN,
       SimpleGraph

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
       transitivity,
       @windowed,
       # reexported from LightGraphs:
       degree, density, local_clustering,
       closeness_centrality, betweenness_centrality,
       laplacian_matrix, laplacian_spectrum,
       # based on (but not exported from) LightGraphs
       averageclustering, transitivity, averagepath
       

include("distance_matrix.jl")
include("matrices.jl")
include("plot.jl")
include("histograms.jl")
include("rqa.jl")
include("radius.jl")
include("graphs.jl")
include("rna.jl")
include("windowed.jl")
include("deprecate.jl")

end
