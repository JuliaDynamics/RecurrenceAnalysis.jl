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
       determinism,
       avgdiag,
       maxdiag,
       divergence,
       rqaentropy,
       trend,
       laminarity,
       trappingtime,
       maxvert,
       rqa,
       sorteddistances,
       @windowed

include("matrices.jl")
include("plot.jl")
include("rqa.jl")
include("radius.jl")
include("windowed.jl")

end
