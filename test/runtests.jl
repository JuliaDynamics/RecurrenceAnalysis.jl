using RecurrenceAnalysis
using Test

function testfile(file, testname=defaultname(file))
    println("running test file $(file)")
    @testset "$testname" begin; include(file); end
    return
end
defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))

@testset "RecurrenceAnalysis tests" begin
    testfile("rmatrix_analytic.jl")
    testfile("rqa_rna_analytic.jl")
    testfile("skeletontest.jl")
end
