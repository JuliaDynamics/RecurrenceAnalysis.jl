using RecurrenceAnalysis
using Test
using Downloads

# Download some test timeseries
tsfolder = joinpath(@__DIR__, "timeseries")
todownload1 = ["$n.csv" for n in 1:4]
todownload = ["test_time_series_lorenz_standard_N_10000_multivariate.csv", "test_time_series_roessler_N_10000_multivariate.csv"]
append!(todownload, todownload1)
repo = "https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/timeseries"
for a in todownload
    Downloads.download(repo*"/"*a, joinpath(tsfolder, a))
end

ti = time()

@testset "RecurrenceAnalysis tests" begin
    include("dynamicalsystems.jl")
    include("smallmatrix.jl")
    include("skeletontest.jl")
    include("deprecations.jl")
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti, digits=3), " seconds or ", round(ti/60, digits=3), " minutes")
