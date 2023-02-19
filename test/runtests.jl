using RecurrenceAnalysis
using Test

function testfile(file, testname=defaultname(file))
    println("running test file $(file)")
    @testset "$testname" begin; include(file); end
    return
end
defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))


# Download some test timeseries
using Downloads
tsfolder = joinpath(@__DIR__, "timeseries")
todownload1 = ["$n.csv" for n in 1:4]
todownload = ["test_time_series_lorenz_standard_N_10000_multivariate.csv", "test_time_series_roessler_N_10000_multivariate.csv"]
append!(todownload, todownload1)
repo = "https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/timeseries"
mkpath(tsfolder)
for a in todownload
    Downloads.download(repo*"/"*a, joinpath(tsfolder, a))
end

ti = time()

@testset "RecurrenceAnalysis tests" begin
    include("dynamicalsystems.jl")
    include("rqa_rna_analytic.jl")
    include("skeletontest.jl")
    include("deprecations.jl")
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti, digits=3), " seconds or ", round(ti/60, digits=3), " minutes")
