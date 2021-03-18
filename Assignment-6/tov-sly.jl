using DifferentialEquations, StaticArrays
using HTTP, CSV
using Interpolations
using Unitful, UnitfulAstro

# Downlaod the  data for SLy equation of state
df = CSV.File(HTTP.get("http://www.ioffe.ru/astro/NSG/NSEOS/sly4.dat").body;
    header=false,skipto=7,delim=" ",ignorerepeated=true)

ρ = df.Column3
p = df.Column4

# Interpolate and then extrapolate the data
p_interp = interpolate((ρ, ), p, Gridded(Linear()))
ρ_interp = interpolate((p,), ρ, Gridded(Linear()))
ρ_exterp = extrapolate(ρ_interp, Line())

# Define the tov odes
function tov(u, params, r)
    m, p = u
    G, C = params
    ρ = ρ_exterp(p)
    dmdr = 4 * pi * r^2 * ρ
    dpdr = - (G * m * ρ / r^2) * (1 + (p/(C^2 * ρ))) * (1 + (4* pi * r^3 * p/(C^2 * m))) * (1 - (2 * G * m/(C^2 * r)))^(-1)
    @SVector [dmdr, dpdr]
end

condition(u, r, integrator) = u[2]
affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition, affect!)

ρc = 1e15
u0 = @SVector [1e-10, p_interp(ρc)]
rspan = [1e-10, 1e10]
G = 6.6743 * 1e-8
c = 3 * 1e10

problem = ODEProblem(tov, u0, rspan, [G, c])
sol = solve(problem, callback=cb)

Msun = uconvert(NoUnits, 1 * u"Msun"/u"g")
sol.t[end]/1e5
sol.u[end][2]
sol.u[end][1]/Msun

# We solve for different values of the central density and plot the mass and radius of the star
ρcs = 10 .^ LinRange(14.5, log10(ρ[end]), 50)

Masses = zeros(length(ρcs))
Radii = zeros(length(ρcs))

for (idx, ρc) in enumerate(ρcs)
    u0 = @SVector [1e-10, p_interp(ρc)]
    problem = ODEProblem(tov, u0, rspan, [G, c])
    sol = solve(problem, callback=cb)
    Radii[idx] = sol.t[end]/1e5
    Masses[idx] = sol.u[end][1]/Msun
end

using Plots, LaTeXStrings
scatter(Radii, Masses, xlabel=L"$r \textrm{ [km]}$", ylabel=L"$\textrm{Mass }[M_\odot]$", label=L"\textrm{SLy}")
savefig("mass_vs_radius_SLy.png")
