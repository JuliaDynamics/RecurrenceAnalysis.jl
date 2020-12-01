using Random

xm = rand(10, 2)
ym = rand(10, 2)
v = Dataset(rand(10, 2))
@testset "Deprecations" begin
    @test_deprecated RecurrenceMatrix(xm, 0.5)
    @test_deprecated CrossRecurrenceMatrix(xm, v, 0.5)
    @test_deprecated CrossRecurrenceMatrix(v, ym, 0.5)
    @test_deprecated CrossRecurrenceMatrix(xm, ym, 0.5)
    @test_deprecated JointRecurrenceMatrix(xm, v, 0.5)
    @test_deprecated JointRecurrenceMatrix(v, ym, 0.5)
    @test_deprecated JointRecurrenceMatrix(xm, ym, 0.5)
end