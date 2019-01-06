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
    @test recurrencerate(rmat) == 33/110
    histograms = recurrencestructures(rmat)
    @test histograms["diagonal"] == [11,7,1,0,1]
    @test histograms["vertical"] == [15,6,2]
    @test histograms["recurrencetimes"] == [3,1,4,4]
    # with theiler window
    @test recurrencerate(rmat, theiler=2) == 22/110
    histograms = recurrencestructures(rmat, theiler=2)
    @test histograms["diagonal"] == [9,4,0,0,1]
    @test histograms["vertical"] == [10,3,2]
    @test histograms["recurrencetimes"] == [1,0,3,0,0,0,2]
    # lmin = 2
    histograms = recurrencestructures(rmat, lmin=2)
    @test histograms["diagonal"] == [0,7,1,0,1]
    @test histograms["vertical"] == [0,6,2]
    @test histograms["recurrencetimes"] == [0,0,0,2,0,0,1]
    # with theiler window
    histograms = recurrencestructures(rmat, theiler=2, lmin=2)
    @test histograms["diagonal"] == [0,4,0,0,1]
    @test histograms["vertical"] == [0,3,2]
    @test histograms["recurrencetimes"] == [0,0,0,0,0,0,1]   
end
