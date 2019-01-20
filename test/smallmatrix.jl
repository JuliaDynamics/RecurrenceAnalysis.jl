using RecurrenceAnalysis, SparseArrays
using Test

cells = [(1,1),(2,2),(3,2),(6,2),(7,2),(8,2),(6,3),(7,3),(10,3),(2,4),(4,4),
    (1,5),(2,5),(3,5),(5,5),(9,5),(10,5),(2,6),(3,6),(6,6),(7,6),(10,6),
    (3,7),(8,7),(4,8),(1,9),(5,9),(7,9),(9,9),(9,10),(10,10),(3,11),(7,11)]
i = [x[1] for x in cells]
j = [x[2] for x in cells]
rmat = CrossRecurrenceMatrix(sparse(i,j,trues(length(i))))
### Graph of the recurrence matrix:
## · : non-recurrent (white) points
## x : recurrent (black) points
## o : recurrent points within Theiler window = 2
#
#    1 2 3 4 5 6 7 8 9 A B
#  1 o · · · x · · · x · ·
#  2 · o · x x x · · · · ·
#  3 · o · · x x x · · · x
#  4 · · · o · · · x · · ·
#  5 · · · · o · · · x · ·
#  6 · x x · · o · · · · ·
#  7 · x x · · o · · x · x
#  8 · x · · · · o · · · ·
#  9 · · · · x · · · o o ·
#  A · · x · x x · · · o · 
#
###

@testset "Recurrence structures" begin
    @testset "Default parameters" begin
        histograms = recurrencestructures(rmat)
        true_histograms = Dict("diagonal" => [11,7,1,0,1],
                               "vertical" => [15,6,2],
                               "recurrencetimes" => [3,1,4,4])
        for k in keys(histograms)
            @test histograms[k] == true_histograms[k]
        end
        rqa_params = rqa(rmat, theiler=0, lmin=1, border=1)
        @test rqa_params["RR"] == 33/110
        @test rqa_params["DET"] == 1.0
        @test rqa_params["L"] == 33/20
        @test rqa_params["Lmax"] == 5
        @test rqa_params["DIV"] == 0.2
        @test rqa_params["ENTR"] ≈ 0.996 atol=0.001
        @test rqa_params["TREND"] ≈ -5.450e-2 atol=0.001
        @test rqa_params["LAM"] == 1.0
        @test rqa_params["TT"] == 33/23
        @test rqa_params["Vmax"] == 3
        @test rqa_params["MRT"] == 33/12
        @test rqa_params["RTE"] ≈ 1.286 atol=0.001
        @test rqa_params["NMPRT"] == 4
    end
    @testset "With Theiler window" begin
        histograms = recurrencestructures(rmat, theiler=2)
        true_histograms = Dict("diagonal" => [9,4,0,0,1],
                               "vertical" => [10,3,2],
                               "recurrencetimes" => [1,0,3,0,0,0,2])
        for k in keys(histograms)
            @test histograms[k] == true_histograms[k]
        end
        rqa_params = rqa(rmat, theiler=2, lmin=1, border=1)
        @test rqa_params["RR"] == 22/110
        @test rqa_params["DET"] == 1.0
        @test rqa_params["L"] == 22/14
        @test rqa_params["Lmax"] == 5
        @test rqa_params["DIV"] == 0.2
        @test rqa_params["ENTR"] ≈ 0.830 atol=0.001
        @test rqa_params["TREND"] ≈ -8.98e-3 atol=1e-6
        @test rqa_params["LAM"] == 1.0
        @test rqa_params["TT"] == 22/15
        @test rqa_params["Vmax"] == 3
        @test rqa_params["MRT"] == 24/6
        @test rqa_params["RTE"] ≈ 1.011 atol=0.001
        @test rqa_params["NMPRT"] == 3
    end
    @testset "With minimum line" begin
        histograms = recurrencestructures(rmat, lmin=2)
        true_histograms = Dict("diagonal" => [11,7,1,0,1],
                               "vertical" => [15,6,2],
                               "recurrencetimes" => [3,1,4,4])
        for k in keys(histograms)
            @test histograms[k] == true_histograms[k]
        end
        rqa_params = rqa(rmat, theiler=0, lmin=2, border=1)
        @test rqa_params["RR"] == 33/110
        @test rqa_params["DET"] == 22/33
        @test rqa_params["L"] == 22/9
        @test rqa_params["Lmax"] == 5
        @test rqa_params["DIV"] == 0.2
        @test rqa_params["ENTR"] ≈ 0.684 atol=0.001
        @test rqa_params["TREND"] ≈ -5.450e-2 atol=0.001
        @test rqa_params["LAM"] == 18/33
        @test rqa_params["TT"] == 18/8
        @test rqa_params["Vmax"] == 3
        @test rqa_params["MRT"] == 33/12
        @test rqa_params["RTE"] ≈ 1.286 atol=0.001
        @test rqa_params["NMPRT"] == 4
    end
    @testset "Theiler and minimum line" begin
        histograms = recurrencestructures(rmat, theiler=2, lmin=2)
        true_histograms = Dict("diagonal" => [9,4,0,0,1],
                               "vertical" => [10,3,2],
                               "recurrencetimes" => [1,0,3,0,0,0,2])
        for k in keys(histograms)
            @test histograms[k] == true_histograms[k]
        end
        rqa_params = rqa(rmat, theiler=2, lmin=2, border=1)
        @test rqa_params["RR"] == 22/110
        @test rqa_params["DET"] == 13/22
        @test rqa_params["L"] == 13/5
        @test rqa_params["Lmax"] == 5
        @test rqa_params["DIV"] == 0.2
        @test rqa_params["ENTR"] ≈ 0.500 atol=0.001
        @test rqa_params["TREND"] ≈ -8.98e-3 atol=1e-6
        @test rqa_params["LAM"] == 12/22
        @test rqa_params["TT"] == 12/5
        @test rqa_params["Vmax"] == 3
        @test rqa_params["MRT"] == 24/6
        @test rqa_params["RTE"] ≈ 1.011 atol=0.001
        @test rqa_params["NMPRT"] == 3
    end
end


