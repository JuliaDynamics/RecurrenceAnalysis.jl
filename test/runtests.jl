using RecurrenceAnalysis
using Base.Test

# Functions for dynamic systems
# Increment of x[ix] through RK4-integration
# dxdh is the array of derivatives, dt is the time increment
function dx_rk4(dxdt, x, ix, dt)
    x = collect(x)
    incr = zeros(length(x))
    k1 = dxdt[ix](x...)
    incr[ix] = dt*k1/2
    k2 = dxdt[ix]((x+incr)...)
    incr[ix] = dt*k2/2
    k3 = dxdt[ix]((x+incr)...)
    incr[ix] = dt*k3
    k4 = dxdt[ix]((x+incr)...)
    dt/6*(k1 + 2k2 + 2k3 + k4)
end
# Integrate dynamical system through n intervals of length dt
function dynamical_system(x0, dxdt, dt, n)
    m = length(x0)
    x = zeros(n,m)
    x[1,:] = collect(x0)
    for t=1:n-1
        dx = zeros(1,m)
        for ix=1:m
            dx[ix] = dx_rk4(dxdt, x[t,:], ix, dt)
        end
        x[t+1,:] = x[t,:] + dx
    end
    x
end

# Test with Lorenz system
lorenz_eq(sigma, rho, beta) = (
  (x,y,z) -> sigma*(y - x),
  (x,y,z) -> x*(rho - z) - y,
  (x,y,z) -> x*y - beta*z)

sigma=10.; rho=28.; b=8/3;
x0 = ones(3)
lorenz_data = dynamical_system(x0, lorenz_eq(sigma,rho,b), .01, 1000)
x = lorenz_data[:,1]
# Look for optimal delay
ami_def = ami(x, (1,30))
@test findmin(ami_def)[2] == 19
ami_st = ami(x, (1,30), "Sturges")
@test ami_st == ami_def
ami_fd = ami(x, (1,30), "FD")
@test findmin(ami_fd)[2] == 16
ami_15 = ami(x, (1,30), 15)
gmi_10 = gmi(x, (1,30), 0.1)
# Distance and recurrence matrices
xe = embed(x, 3, 17)
dmat = distancematrix(xe)
dmat_euc = distancematrix(xe,"euclidean")
@test all(dmat .<= dmat_euc)
rmat = recurrencematrix(xe, 0.05)
crmat = crossrecurrencematrix(x, lorenz_data[:,3], 0.05)
jrmat = jointrecurrencematrix(x, lorenz_data[:,3], 0.05)
# RQA
rr = recurrencerate(rmat)
@test recurrencerate(rmat,theiler=1) < rr
det = determinism(rmat)
@test determinism(rmat, theiler=5, lmin=5) < det
lmean = avgdiag(rmat)
@test avgdiag(rmat, lmin=5) > lmean
@test maxdiag(rmat) > lmean
@test divergence(rmat, theiler=5) == 1/maxdiag(rmat, theiler=5)
ent = entropy(rmat)
tnd = trend(rmat, border=20, theiler=5)
lam = laminarity(rmat)
tt = trappingtime(rmat, lmin=5)
@test maxvert(rmat) > tt
