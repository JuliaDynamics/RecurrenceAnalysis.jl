#=
This file benchmarks the parallel and serial computation of a recurrence matrix
=#
using DynamicalSystemsBase, RecurrenceAnalysis
using Statistics, BenchmarkTools, Base.Threads

Ns = round.(Int, 10.0 .^ (2:0.5:4.5))
ro = Systems.roessler()

println("I am using $(nthreads()) threads. I am reporting results as:")
println("time of (parallel/serial) for 3D and 1D trajectories.")
for N in Ns
    tr = trajectory(ro, N*0.1; dt = 0.1, Tr = 10.0)
    println("For N = $(length(tr))...")
    x = tr[:, 1]
    b_serial   = @benchmark RecurrenceMatrix($(tr), 5.0; metric = Euclidean(), parallel = false)
    b_parallel = @benchmark RecurrenceMatrix($(tr), 5.0; metric = Euclidean(), parallel = true)
    b_s_1 = @benchmark RecurrenceMatrix($(x), 2.0; parallel = false)
    b_p_1 = @benchmark RecurrenceMatrix($(x), 2.0; parallel = true)
    threeD = median(b_serial.times)/median(b_parallel.times)
    oneD = median(b_s_1.times)/median(b_p_1.times)
    threeD, oneD = round.((threeD, oneD); sigdigits = 3)
    println("3D=$(threeD), 1D=$(oneD)")
end
