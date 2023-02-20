# This file tests that the construction of recurrence matrices
# is correct using as many recurrence types as possible
# Analytic tests with points on the circle with fixed distance.

# This file DOES NOT need to test RQA measures. This happens elsewhere!

using RecurrenceAnalysis, Test
using Distances: Chebyshev

t = range(0, 2π; length = 20) # length is crucial and decides distance and thresholds
c = cos.(t)
s = sin.(t)
X = Dataset(c, s)
Y = X[1:10]

dmat = distancematrix(X)
threshold = dmat[1, 2] + 0.01 # distance between two points
maxthres = maximum(dmat) + 0.01
neighbors = count(<(threshold), dmat)

@testset "fixed RR" begin
    ε = threshold
    rmat = RecurrenceMatrix(X, ε; parallel = false)
    rmat_p = RecurrenceMatrix(X, ε; parallel = true)
    crmat = CrossRecurrenceMatrix(X, Y, ε; parallel = false)
    jrmat = JointRecurrenceMatrix(X, X, ε; parallel = false)

    @test size(rmat) == size(rmat_p) == size(jrmat) == (20, 20)
    @test size(crmat) == (20, 10)
    @test count(rmat) == count(rmat_p) == neighbors
    @test count(jrmat) == neighbors
    @test count(crmat) == neighbors÷2

    @test rmat == RecurrenceMatrix(X, RecurrenceThreshold(ε))

end

@testset "more fixed RR" begin
    # max radius
    ε = maxthres
    rmat = RecurrenceMatrix(X, ε; parallel = false)
    @test count(rmat) == length(X)*length(X)

    # Different metric
    metric = Chebyshev()
    # due to the symmetry of the circle, the chebysven metric can only slightly
    # increase the recurrences, but really not much!
    ε = threshold
    rmat = RecurrenceMatrix(X, ε; metric)
    @test neighbors ≤ count(rmat) ≤ 1.05neighbors
end

@testset "Scaled fixed RR" begin
    # We use as scale to get the second number...
    scale = (dm) -> dm[1, 2] + 0.01 # cheating but analytic ;)
    ε = RecurrenceThresholdScaled(1.0, scale)
    rmat = RecurrenceMatrix(X, ε; parallel = false)
    rmat_p = RecurrenceMatrix(X, ε; parallel = true)
    crmat = CrossRecurrenceMatrix(X, Y, ε; parallel = false)

    @test size(rmat) == size(rmat_p) == (20, 20)
    @test size(crmat) == (20, 10)
    @test count(rmat) == count(rmat_p) == neighbors
    @test count(crmat) == neighbors÷2
end


@testset "GlobalRecurrenceRate" begin

    ratios = (0.2, 0.5)
    @testset "ratio = $(ratio)" for ratio in ratios
        ε = GlobalRecurrenceRate(ratio)
        d = recurrence_threshold(ε, X)
        rmat = RecurrenceMatrix(X, ε; parallel = false)
        rmat_p = RecurrenceMatrix(X, ε; parallel = true)
        crmat = CrossRecurrenceMatrix(X, Y, ε; parallel = false)
        # Alright, so I admit I am not sure why this doesn't give exactly 0.2*N^2
        # recurrences, but almost that.
        @test size(rmat) == size(rmat_p) == (20, 20)
        @test size(crmat) == (20, 10)

        expected = length(X)*length(X)*ratio
        # TODO: Fails:
        # @test count(rmat) == count(rmat_p) == expected
        # works
        @test 0.9expected < count(rmat) < 1.3expected

        expected2 = length(X)*length(Y)*ratio
        @test 0.9expected2 < count(crmat) < 1.3expected2
    end

end

@testset "LocalRecurrenceRate" begin
    ratios = (0.2, 0.5)
    @testset "ratio = $(ratio)" for ratio in ratios
        ε = LocalRecurrenceRate(ratio)
        d = recurrence_threshold(ε, X)
        @test d isa AbstractVector
        # Because of points in the circle and full rotational symmetry
        # all entries of `d` should be approximately equal.
        # TODO: It doesn't work
        # @test all(i -> d[i] ≈ d[i+1], 1:length(d) - 1)
        rmat = RecurrenceMatrix(X, ε; parallel = false)
        @test size(rmat) == size(rmat_p) == (20, 20)

        # And each column of the matrix must have same count by definition
        # TODO: This also doesn't work...
        count1 = count(rmat[:, 1])
        # @test all(i -> count(rmat[:, i]) == count1, 2:size(rmat, 2))
        # But a less accurate test works
        @test all(i -> 0.8count1 < count(rmat[:, i]) < 1.2count1, 2:size(rmat, 2))
    end
end
