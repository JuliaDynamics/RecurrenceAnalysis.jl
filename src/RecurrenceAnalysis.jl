module RecurrenceAnalysis

using Distances, Statistics, LinearAlgebra, SparseArrays, DelayEmbeddings, StaticArrays
import Base.Meta.parse

export embed,
       distancematrix,
       recurrencematrix,
       crossrecurrencematrix,
       jointrecurrencematrix,
       recurrenceplot,
       recurrencerate,
       determinism,
       avgdiag,
       maxdiag,
       divergence,
       entropy,
       trend,
       laminarity,
       trappingtime,
       maxvert,
       rqa,
       autocorrelation,
       ami,
       gmi,
       fnn,
       afnn,
       ffnn,
       sorteddistances,
       @windowed

# column values in sparse matrix (parallel to rowvals)
function colvals(x::SparseMatrixCSC)
    cv = zeros(Int,nnz(x))
    @inbounds for c=1:size(x,2)
        cv[nzrange(x,c)] .= c
    end
    cv
end

include("matrices.jl")
include("plot.jl")
include("rqa.jl")
include("delay.jl")
include("dimension.jl")
include("radius.jl")
include("windowed.jl")

end
