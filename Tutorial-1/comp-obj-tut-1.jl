using DifferentialEquations
using Printf
using Plots
using LaTeXStrings
gr()
# define the set of odes
function stellar_structure!(du, u, params, m)
    a, b, n, s, nu, lambda, A = params
    p, x, t, l = u
    du[1] = - (m / x^4)
    du[2] = (t^b / (x^2 * p^a))
    du[3] = - (p^(a * n) * l) / (x^4 * t^(3 + s + b * n))
    du[4] = A * p^(a * lambda) * t^(nu - (b * lambda))
    # @printf("m = %.f \t dp/dm = %.f dx/dm = %.f \t dt/dm = %.f \t dl/dm = %.f \t p = %.f\n", m, du[1], du[2], du[3], du[4], p)
end

# We define the condition to find the event where p becomes less than some tolerace set by user.
condition(u, m, integrator) =  u[1] - 1e-4
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition, affect!)

# Initial conditions, at m = 0, as per eq 1.2
p0 = 1.0
x0 = 1e-3
t0 = 1.0
l0 = 0.0

u0 = [p0, x0, t0, l0]

# parameters
a = 1.0
b = 1.0
n = 0.0
s = 0.0
nu = 4.0
lambda = 1.0
A = 0.5

params = [a, b, n, s, nu, lambda, A]
mspan = (1e-9, 20.0)
prob1 = ODEProblem(stellar_structure!, u0, mspan, params)
sol1 = solve(prob1, alg_hints=[:stiff], callback=cb)
plot(sol1, vars=(0, 1), lw = 2, xaxis = (L"m"), yaxis = (L"p(m)"))
savefig("p_vs_m.pdf")
