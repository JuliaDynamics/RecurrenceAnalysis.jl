using RecurrenceAnalysis
using Test

# Download some test timeseries
# Old:
repo = "https://raw.githubusercontent.com/JuliaDynamics/NonlinearDynamicsTextbook/master/exercise_data"
tsfolder = joinpath(@__DIR__, "timeseries")
todownload = ["$n.csv" for n in 1:4]

mkpath(tsfolder)
for a in todownload
    download(repo*"/"*a, joinpath(tsfolder, a))
end

#New:
todownload = ["test_time_series_lorenz_standard_N_10000_multivariate.csv", "test_time_series_roessler_N_10000_multivariate.csv"]
repo = "https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/timeseries"
for a in todownload
    download(repo*"/"*a, joinpath(tsfolder, a))
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
