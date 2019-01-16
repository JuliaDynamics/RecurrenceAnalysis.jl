module RecurrenceAnalysis

using Distances, Statistics, LinearAlgebra, SparseArrays, DelayEmbeddings, StaticArrays
import Base.Meta.parse

export RecurrenceMatrix, CrossRecurrenceMatrix, JointRecurrenceMatrix
export embed,
       reconstruct,
       Dataset,
       distancematrix,
       recurrencematrix,
       crossrecurrencematrix,
       jointrecurrencematrix,
       recurrenceplot,
       recurrencerate,
       recurrencestructures,
       dl_average, dl_max, dl_entropy,
       vl_average, vl_max, vl_entropy,
       rt_average, rt_max, rt_entropy,
       determinism,
       avgdiag,
       maxdiag,
       divergence,
       rqaentropy,
       trend,
       laminarity,
       trappingtime,
       maxvert,
       meanrecurrencetime,
       nmprt,
       rqa,
       sorteddistances,
       @windowed

include("matrices.jl")
include("plot.jl")
include("histograms.jl")
include("rqa.jl")
include("radius.jl")
include("windowed.jl")

end
