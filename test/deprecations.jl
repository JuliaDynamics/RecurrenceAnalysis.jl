using Random

xm = rand(10, 2)
ym = rand(10, 2)
v = Dataset(rand(10, 2))
@testset "Deprecations" begin
    @test_logs (:warn, "`RecurrenceMatrix(x::AbstractMatrix, ε; kwargs...)` is deprecated, use `RecurrenceMatrix(Dataset(x), ε; kwargs...)`") RecurrenceMatrix(xm, 0.5)
    @test_logs (:warn, "`CrossRecurrenceMatrix(x::AbstractMatrix, y, ε; kwargs...)` is deprecated, use `CrossRecurrenceMatrix(Dataset(x), y, ε; kwargs...)`") CrossRecurrenceMatrix(xm, v, 0.5)
    @test_logs (:warn, "`CrossRecurrenceMatrix(x, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `CrossRecurrenceMatrix(x, Dataset(y), ε; kwargs...)`") CrossRecurrenceMatrix(v, ym, 0.5)
    @test_logs (:warn, "`CrossRecurrenceMatrix(x::AbstractMatrix, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `CrossRecurrenceMatrix(Dataset(x), Dataset(y), ε; kwargs...)`") CrossRecurrenceMatrix(xm, ym, 0.5)
    @test_logs (:warn, "`JointRecurrenceMatrix(x::AbstractMatrix, y, ε; kwargs...)` is deprecated, use `JointRecurrenceMatrix(Dataset(x), y, ε; kwargs...)`") JointRecurrenceMatrix(xm, v, 0.5)
    @test_logs (:warn, "`JointRecurrenceMatrix(x, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `JointRecurrenceMatrix(x, Dataset(y), ε; kwargs...)`") JointRecurrenceMatrix(v, ym, 0.5)
    @test_logs (:warn, "`JointRecurrenceMatrix(x::AbstractMatrix, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `JointRecurrenceMatrix(Dataset(x), Dataset(y), ε; kwargs...)`") JointRecurrenceMatrix(xm, ym, 0.5)
end
