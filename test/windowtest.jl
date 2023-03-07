using RecurrenceAnalysis, Test

t = range(0, 2π; length = 300) # length is crucial and decides distance and thresholds
c = cos.(t)
s = sin.(t)
X = StateSpaceSet(c, s)

dmat = distancematrix(X)
threshold = dmat[1, 2] + 0.01 # distance between two points

@testset "window function" begin
    
    rec = RecurrenceMatrix(X, threshold)

    #classical

    @test length(windowed(rec, recurrencerate, 30, 30)) == 10
    @test length(windowed(rec, determinism, 30, 30)) == 10
    @test length(windowed(rec, dl_average, 30, 30)) == 10
    @test length(windowed(rec, dl_max, 30, 30)) == 10
    @test length(windowed(rec, dl_entropy, 30, 30)) == 10
    @test length(windowed(rec, divergence, 30, 30)) == 10
    @test length(windowed(rec, trend, 30, 30)) == 10

    #extended

    @test length(windowed(rec, laminarity, 30, 30)) == 10
    @test length(windowed(rec, trappingtime, 30, 30)) == 10
    @test length(windowed(rec, vl_average, 30, 30)) == 10
    @test length(windowed(rec, vl_max, 30, 30)) == 10
    @test length(windowed(rec, vl_entropy, 30, 30)) == 10

    #recurrence time

    @test length(windowed(rec, meanrecurrencetime, 30, 30)) == 10
    @test length(windowed(rec, nmprt, 30, 30)) == 10
    @test length(windowed(rec, rt_entropy, 30, 30)) == 10
    @test length(windowed(rec, rt_average, 30, 30)) == 10

    #check if reacts to changes in width and step
    @test length(windowed(rec, recurrencerate, 30, 1)) == 271
    @test length(windowed(rec, recurrencerate, 10, 1)) == 291
end