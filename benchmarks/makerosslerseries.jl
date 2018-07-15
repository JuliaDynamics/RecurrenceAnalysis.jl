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
        dx = zeros(m)
        for ix=1:m
            dx[ix] = dx_rk4(dxdt, x[t,:], ix, dt)
        end
        x[t+1,:] = x[t,:][:] + dx
    end
    x
end
# Rössler system
rossler_eq(a, b, c) = (
    (x,y,z) -> -(y+z),
    (x,y,z) -> x + a*y,
    (x,y,z) -> b + z*(x-c))

a=.25; b=.25; c=4
x0 = zeros(3)
m = zeros(3000,40)
for i=0:11
    x = randn(3)
    rossler_data = dynamical_system(x0, rossler_eq(a,b,c), .05, 4000)
    m[:,2i.+(1:2)] .= rossler_data[1001:end,1:2]
end
writedlm("rossler.txt",m,'\t')

