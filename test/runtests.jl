using RecurrenceAnalysis
using Test


ti = time()

@testset "RecurrenceAnalysis tests" begin
include("dynamicalsystems.jl")
include("smallmatrix.jl")
include("deprecations.jl")
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti, digits=3), " seconds or ", round(ti/60, digits=3), " minutes")
