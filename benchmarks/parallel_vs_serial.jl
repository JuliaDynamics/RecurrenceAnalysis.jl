#=
This file benchmarks the parallel and serial computation of a recurrence matrix
=#
using DynamicalSystemsBase, RecurrenceAnalysis
using Statistics, BenchmarkTools, Base.Threads

function pretty_time(b::BenchmarkTools.Trial)
    med = median(b)

    return string(
        BenchmarkTools.prettytime(BenchmarkTools.time(med)),
        " (",
        BenchmarkTools.prettypercent(BenchmarkTools.gcratio(med)),
        " GC)"
    )

end

Ns = round.(Int, 10.0 .^ (2:0.5:4.5))
ro = Systems.roessler()

println("I am using $(nthreads()) threads.")
for N in Ns
    # set up datasets
    tr = trajectory(ro, N*0.1; dt = 0.1, Tr = 10.0)
    printstyled("For N = $(length(tr))..."; color = :blue, bold = true)
    println()
    x = tr[:, 1]

    # calculate only the lower triangle
    b_serial   = @benchmark RecurrenceMatrix($(tr), 5.0; metric = Euclidean(), parallel = false)
    b_parallel = @benchmark RecurrenceMatrix($(tr), 5.0; metric = Euclidean(), parallel = true)
    b_s_1 = @benchmark RecurrenceMatrix($(x), 2.0; parallel = false)
    b_p_1 = @benchmark RecurrenceMatrix($(x), 2.0; parallel = true)

    # calculate triangular metrics
    threeD = median(b_serial.times)/median(b_parallel.times)
    oneD = median(b_s_1.times)/median(b_p_1.times)
    threeD, oneD = round.((threeD, oneD); sigdigits = 3)
    printstyled("Lower triangle only"; color = :green)
    println("\nSpeedups: \n    3D=$threeD\n    1D=$oneD")
    println(
        """
        Raw timings:
            3D serial:   $(pretty_time(b_serial))
            3D parallel: $(pretty_time(b_parallel))
            1D serial:   $(pretty_time(b_s_1))
            1D parallel: $(pretty_time(b_p_1))
        """
    )

    # calculate the full matrix
    b_serial_full   = @benchmark CrossRecurrenceMatrix($(tr), $(tr), 5.0; metric = Euclidean(), parallel = false)
    b_parallel_full = @benchmark CrossRecurrenceMatrix($(tr), $(tr), 5.0; metric = Euclidean(), parallel = true)
    b_s_1_full = @benchmark CrossRecurrenceMatrix($(x), $(x), 2.0; parallel = false)
    b_p_1_full = @benchmark CrossRecurrenceMatrix($(x), $(x), 2.0; parallel = true)

    # calculate full metrics
    threeD = median(b_serial_full.times)/median(b_parallel_full.times)
    oneD = median(b_s_1_full.times)/median(b_p_1_full.times)
    threeD, oneD = round.((threeD, oneD); sigdigits = 3)
    printstyled("Full matrix"; color = :green)
    println("\nSpeedups: \n    3D=$threeD\n    1D=$oneD")
    println(
        """
        Raw timings:
            3D serial:   $(pretty_time(b_serial_full))
            3D parallel: $(pretty_time(b_parallel_full))
            1D serial:   $(pretty_time(b_s_1_full))
            1D parallel: $(pretty_time(b_p_1_full))
        """
    )
    println()
end
