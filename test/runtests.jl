using RecurrenceAnalysis
using Test, DynamicalSystemsBase

systems = (Systems.lorenz(),)

testvalues = (
[426/234^2, 364/426, 364/50, 23, 2.02058292, 3.21794972e-7, 92/426, 92/30, 4],
)

u0s = (ones(3),)

dt = 0.01
for (ds, vals, u0) in zip(systems, testvalues, u0s)

    data = trajectory(ds, dt*1000, u0; dt = dt)
    x = data[501:2:end,1]

    # Look for optimal threshold
    dd, rr = sorteddistances(x, theiler=1)
    # Distance and recurrence matrices
    xe = embed(x, 3, 8)
    rmat = RecurrenceMatrix(xe, 1.5, metric="max")
    y = data[701:2:end,3]
    crmat = CrossRecurrenceMatrix(x, y, 1.5)
    jrmat = JointRecurrenceMatrix(x, y, 1.5)
    # Recurrence plot
    crp = recurrenceplot(crmat, width=125)
    @test size(crp)[1] == 75
    # RQA
    rqapar = rqa(rmat, theiler=2, lmin=3, border=20)
    tol = 1e-5
    @test rqapar["RR"] ≈ vals[1]    atol = tol
    @test rqapar["DET"] ≈ vals[2]   atol = tol
    @test rqapar["L"] ≈ vals[3]     atol = tol
    @test rqapar["Lmax"] ≈ vals[4]  atol = tol
    @test rqapar["ENT"] ≈ vals[5]   atol = tol
    @test rqapar["TND"] ≈ vals[6]   atol = tol
    @test rqapar["LAM"] ≈ vals[7]   atol = tol
    @test rqapar["TT"] ≈ vals[8]    atol = tol
    @test rqapar["Vmax"] ≈ vals[9]  atol = tol
    rqadiag = rqa(rmat, theiler=2, lmin=3, border=20, onlydiagonal=true)
    @test all([rqapar[p]==rqadiag[p] for p in keys(rqadiag)])
    # Fixed rate for recurrence matrix
    rmat2 = RecurrenceMatrix(xe[1:3:end,:], 0.05; fixedrate=true)
    @test .049 < recurrencerate(rmat2) < .051
    # Windowed RQA
    rmatw = @windowed RecurrenceMatrix(xe, 1.5, metric=RecurrenceAnalysis.Chebyshev()) 50
    crmatw = @windowed(CrossRecurrenceMatrix(x, y, 1.5),30)
    @windowed jrmatw = JointRecurrenceMatrix(x, y, 1.5) 30
    @test jrmatw[33 .+ (1:30), 33 .+ (1:30)] == jrmat[33 .+ (1:30), 33 .+ (1:30)]
    @windowed(rrw = recurrencerate(rmatw), width=50, step=40)
    @windowed rqaw = rqa(rmatw) width=50 step=40
    @test rqaw["RR"] == rrw
end
