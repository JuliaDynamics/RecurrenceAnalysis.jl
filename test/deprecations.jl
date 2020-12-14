using Random

xm = rand(10, 2)
ym = rand(10, 2)
v = Dataset(rand(10, 2))
@testset "Deprecations" begin
    @test_logs (:warn, "`RecurrenceMatrix{WithinRange}(x::AbstractMatrix, ε; kwargs...)` is deprecated, use `RecurrenceMatrix{WithinRange}(Dataset(x), ε; kwargs...)`") RecurrenceMatrix(xm, 0.5)
    @test_logs (:warn, "`CrossRecurrenceMatrix{WithinRange}(x::AbstractMatrix, y, ε; kwargs...)` is deprecated, use `CrossRecurrenceMatrix{WithinRange}(Dataset(x), y, ε; kwargs...)`") CrossRecurrenceMatrix(xm, v, 0.5)
    @test_logs (:warn, "`CrossRecurrenceMatrix{WithinRange}(x, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `CrossRecurrenceMatrix{WithinRange}(x, Dataset(y), ε; kwargs...)`") CrossRecurrenceMatrix(v, ym, 0.5)
    @test_logs (:warn, "`CrossRecurrenceMatrix{WithinRange}(x::AbstractMatrix, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `CrossRecurrenceMatrix{WithinRange}(Dataset(x), Dataset(y), ε; kwargs...)`") CrossRecurrenceMatrix(xm, ym, 0.5)
    @test_logs (:warn, "`JointRecurrenceMatrix{WithinRange}(x::AbstractMatrix, y, ε; kwargs...)` is deprecated, use `JointRecurrenceMatrix{WithinRange}(Dataset(x), y, ε; kwargs...)`") JointRecurrenceMatrix(xm, v, 0.5)
    @test_logs (:warn, "`JointRecurrenceMatrix{WithinRange}(x, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `JointRecurrenceMatrix{WithinRange}(x, Dataset(y), ε; kwargs...)`") JointRecurrenceMatrix(v, ym, 0.5)
    @test_logs (:warn, "`JointRecurrenceMatrix{WithinRange}(x::AbstractMatrix, y::AbstractMatrix, ε; kwargs...)` is deprecated, use `JointRecurrenceMatrix{WithinRange}(Dataset(x), Dataset(y), ε; kwargs...)`") JointRecurrenceMatrix(xm, ym, 0.5)
end
