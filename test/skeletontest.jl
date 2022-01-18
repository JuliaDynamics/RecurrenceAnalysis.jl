using RecurrenceAnalysis, SparseArrays
using Test, Statistics, LinearAlgebra, DelimitedFiles

### Skeletonization of RPs
@testset "Skeletonized sinusoidal recurrence structures" begin
    data = sin.(2*π .* (0:1000)./ 60)
    Y = embed(data, 3, 15)

    RP = RecurrenceMatrix{FAN}(Y, 0.2)
    RP2 = RecurrenceMatrix(Y, 0.25; fixedrate=true)
    RP3 = RecurrenceMatrix(Y, 0.9)

    RP_skel = skeletonize(RP)
    RP_skel2= skeletonize(RP2)
    RP_skel3= skeletonize(RP3)

    @test RP_skel2 == RP_skel3

    tauRR1 = RecurrenceAnalysis.tau_recurrence(RP_skel)
    tauRR2 = RecurrenceAnalysis.tau_recurrence(RP_skel2)
    tauRR3 = RecurrenceAnalysis.tau_recurrence(RP_skel3)

    @test tauRR2 == tauRR3

    c1 = findall(x->x!=0, diff(tauRR1))
    c2 = findall(x->x!=0, diff(tauRR2))

    @test mean(diff(c1[3:2:end])) == 60
    @test mean(diff(c2[3:2:end])) == 60
end

@testset "Skeletonized chaotic recurrence structures" begin
    data = readdlm(joinpath(tsfolder, "test_time_series_lorenz_standard_N_10000_multivariate.csv"))

    RP = RecurrenceMatrix(Dataset(data[1:250,:]), 0.05; fixedrate=true)
    RP_skel = skeletonize(RP)

    @test RP_skel[135,61] == 1
    @test RP_skel[135,60] == 0
    @test RP_skel[135,60] == 0
    @test RP_skel[136,61] == 0
    @test RP_skel[134,61] == 0

    @test isempty(findall(!iszero,diag(RP_skel,1))) == true
    @test isempty(findall(!iszero,diag(RP_skel,-1))) == true

    d = Bool.(RP_skel[61:79,134:153])

    @test length(findall(diag(d,1))) == 18
    @test isempty(findall(diag(d))) == true
end
