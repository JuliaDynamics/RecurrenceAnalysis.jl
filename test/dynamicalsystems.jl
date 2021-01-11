using RecurrenceAnalysis
using DynamicalSystemsBase, Random, Statistics, SparseArrays
using LightGraphs, LinearAlgebra
using Test

RA = RecurrenceAnalysis

rng = Random.seed!(194)
# Trajectories of 200 points
# Examples of the Hénon map based on:
#     Webber CL & Zbilut JP, "Recurrence Quantification Analysis of Nonlinear
#     Dynamical Systems", in: Riley MA & Van Orden GC, Tutorials in Contemporary
#     Nonlinear Methods for the Behavioral Sciences, 2005, 26-94.
#     https://www.nsf.gov/pubs/2005/nsf05057/nmbs/nmbs.pdf
#
trajectories = Dict(
    "Sine wave" => RA.Dataset(map(x->[sin.(x) cos.(x)], StepRangeLen(0.0,0.2,200))),
    "White noise" => RA.Dataset(randn!(zeros(200,2))),
    "Hénon (chaotic)" => trajectory(Systems.henon(a=1.4, b=0.3), 199, Ttr=1000),
    "Hénon (periodic)" => trajectory(Systems.henon(a=1.054, b=0.3), 199, Ttr=1000)
)
embed_params = Dict(       #(m, τ)
    "Sine wave"   => (9, 7),
    "White noise" => (1, 1),
    "Hénon (chaotic)" => (3, 1),
    "Hénon (periodic)" => (3, 1)
)
rqa_threshold = Dict(
    "Sine wave"   => 0.15,
    "White noise" => 0.15,
    "Hénon (chaotic)" => 0.15,
    "Hénon (periodic)" => 0.15
)

dict_keys = ["Sine wave","White noise","Hénon (chaotic)","Hénon (periodic)"]
@testset "$k" for k in dict_keys
    data = trajectories[k]
    x = data[:,1].*0.00000001.*randn(length(data[:,1]))
    y = data[1:100,2]
    if k ≠ "White noise"
        xe = embed(x, embed_params[k]...)
        ye = embed(y, embed_params[k]...)
    else
        xe = x
        ye = copy(x)
    end

    # Distance and recurrence matrices
    ε = rqa_threshold[k]
    dmat = distancematrix(xe, ye)
    crmat = CrossRecurrenceMatrix(xe, ye, ε)
    @test Matrix(crmat) == (dmat .≤ ε)
    rmat = RecurrenceMatrix(xe, ε; parallel = false)
    rmat_p = RecurrenceMatrix(xe, ε; parallel = true)
    jrmat = JointRecurrenceMatrix(xe, ye, ε; parallel = false)
    sz = size(jrmat)
    @test (jrmat.data .& rmat.data[1:sz[1],1:sz[2]]) == jrmat.data
    @test (jrmat.data .& rmat_p.data[1:sz[1],1:sz[2]]) == jrmat.data
    # Compare metrics
    m_euc = rmat
    m_euc_p = rmat_p
    m_max = RecurrenceMatrix(xe, ε, metric="max")
    m_min = RecurrenceMatrix(xe, ε, metric="manhattan")
    @test (m_max.data .& m_euc.data) == m_euc.data
    @test (m_euc.data .& m_min.data) == m_min.data
    @test (m_max.data .& m_euc_p.data) == m_euc.data
    @test (m_euc_p.data .& m_min.data) == m_min.data
    # Compare scales
    m_scalemax = CrossRecurrenceMatrix(xe, ye, 0.1, scale=maximum)
    m_scalemean = CrossRecurrenceMatrix(xe, ye, 0.1, scale=mean)
    m_scalemax_p = CrossRecurrenceMatrix(xe, ye, 0.1, scale=maximum, parallel = true)
    m_scalemean_p = CrossRecurrenceMatrix(xe, ye, 0.1, scale=mean, parallel = true)
    @test (m_scalemax.data .& m_scalemean.data) == m_scalemean.data
    @test (m_scalemax_p.data .& m_scalemean_p.data) == m_scalemean.data
    # Fixed rate for recurrence matrix
    crmat_fixed = CrossRecurrenceMatrix(xe, ye, 0.05; fixedrate=true)
    crmat_fixed_p = CrossRecurrenceMatrix(xe, ye, 0.05; fixedrate=true, parallel = true)
    @test .04 < recurrencerate(crmat_fixed) < .06
    @test .04 < recurrencerate(crmat_fixed_p) < .06
    @test recurrencerate(crmat_fixed) ≈ recurrencerate(crmat_fixed_p)
    # fan method for recurrence threshold
    cr_fan = CrossRecurrenceMatrix{FAN}(xe, ye, 0.05; fixedrate=true, parallel = false)
    cr_fan_p = CrossRecurrenceMatrix{FAN}(xe, ye, 0.05; parallel = true)
    rp_fan = RecurrenceMatrix{FAN}(xe, 0.05; fixedrate=true, parallel = false)
    rp_fan_p = RecurrenceMatrix{FAN}(Dataset(xe), 0.05; parallel = true)
    n = length(xe)
    @test all(.04 < nnz(cr_fan[:,i])/n < .06 for i=1:size(cr_fan, 2))
    @test all(.04 < (nnz(rp_fan[:,i])-1)/n < .06 for i=1:size(rp_fan, 2))
    @test .04 < recurrencerate(cr_fan) < .06
    @test .04 < recurrencerate(cr_fan_p) < .06
    @test .04 < recurrencerate(rp_fan) < .06
    @test .04 < recurrencerate(rp_fan_p) < .06
    @test cr_fan == cr_fan_p
    @test rp_fan == rp_fan_p
    @test recurrencerate(cr_fan) ≈ recurrencerate(cr_fan_p)
    @test recurrencerate(rp_fan_p; theiler=0) > recurrencerate(rp_fan_p)

    # Recurrence plot
    crp = grayscale(crmat, width=125)
    szplot = size(crp)
    szmat  = size(crmat)
    @test szplot[1] ≈ szplot[2]*szmat[2]/szmat[1] atol = 1

    # RQA
    rqapar = rqa(rmat, theiler=1, lmin=3, border=20)
    rqadiag = rqa(rmat, theiler=1, lmin=3, border=20, onlydiagonal=true)
    for p in keys(rqadiag)
        @test rqapar[p]==rqadiag[p]
    end

    # Windowed RQA
    rmatw = @windowed RecurrenceMatrix(xe, ε, metric=RecurrenceAnalysis.Chebyshev()) 50
    @windowed RecurrenceMatrix{FAN}(xe, ε) 50 # not meaningful, only to check that it does not error
    crmatw = @windowed(CrossRecurrenceMatrix(xe, ye, ε),30)
    @windowed jrmatw = JointRecurrenceMatrix(xe, ye, ε) 30
    @test jrmatw[3 .+ (1:30), 3 .+ (1:30)] == jrmat[3 .+ (1:30), 3 .+ (1:30)]
    @windowed(rrw = recurrencerate(rmatw), width=50, step=40)
    @windowed rqaw = rqa(rmatw) width=50 step=40
    @test rqaw[:RR] == rrw

    # Recurrence networks
    graph = SimpleGraph(rmat)
    amat = adjacency_matrix(graph)
    @test amat == Matrix(rmat) - I
    rna_dict = rna(rmat)
    rna_dict[:density] == recurrencerate(rmat)
end
