using Random

xm = rand(10, 2)
ym = rand(10, 2)
v = Dataset(rand(10, 2))
@testset "Deprecations" begin
    @test_logs (:warn, "`RecurrenceMatrix{FixedRange}(x::AbstractMatrix, ε; kwargs...)` is deprecated, use `RecurrenceMatrix{FixedRange}(Dataset(x), ε; kwargs...)`") RecurrenceMatrix(xm, 0.5)
    @test_logs (:warn, "`CrossRecurrenceMatrix{FixedRange}(x::AbstractMatrix, y, ε; kwargs...)` is deprecated, use `CrossRecurrenceMatrix{FixedRange}(Dataset(x), y, ε; kwargs...)`") CrossRecurrenceMatrix(xm, v, 0.5)
    @test_logs (:warn, "`CrossRecurrenceMatrix{FixedRange}(x, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `CrossRecurrenceMatrix{FixedRange}(x, Dataset(y), ε; kwargs...)`") CrossRecurrenceMatrix(v, ym, 0.5)
    @test_logs (:warn, "`CrossRecurrenceMatrix{FixedRange}(x::AbstractMatrix, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `CrossRecurrenceMatrix{FixedRange}(Dataset(x), Dataset(y), ε; kwargs...)`") CrossRecurrenceMatrix(xm, ym, 0.5)
    @test_logs (:warn, "`JointRecurrenceMatrix{FixedRange}(x::AbstractMatrix, y, ε; kwargs...)` is deprecated, use `JointRecurrenceMatrix{FixedRange}(Dataset(x), y, ε; kwargs...)`") JointRecurrenceMatrix(xm, v, 0.5)
    @test_logs (:warn, "`JointRecurrenceMatrix{FixedRange}(x, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `JointRecurrenceMatrix{FixedRange}(x, Dataset(y), ε; kwargs...)`") JointRecurrenceMatrix(v, ym, 0.5)
    @test_logs (:warn, "`JointRecurrenceMatrix{FixedRange}(x::AbstractMatrix, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `JointRecurrenceMatrix{FixedRange}(Dataset(x), Dataset(y), ε; kwargs...)`") JointRecurrenceMatrix(xm, ym, 0.5)
end
