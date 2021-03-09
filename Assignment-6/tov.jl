# Solve TOV equation for neutron star
# Md Arif Shaikh, Postdoc, ICTS-TIFR, March 1 2021"

using CSV, DataFrames, HTTP
using Plots, LaTeXStrings
using Interpolations 
using DifferentialEquations, StaticArrays

# First we solve Newtonian stellar structure equation for polytropic equations
function StellarStrunctureNewt(u, params, r)
	m, ρ = u
	G, K, ρ0 = params
	dpdρ = K * abs(ρ)^(2/3)*(5 - 4 * (ρ0/abs(ρ))^(1/3))/3
	dpdr = - G * m * ρ / r^2
	dρdr = dpdr/dpdρ
	dmdr = 4 * pi * r^2 * ρ
	@SVector [dmdr, dρdr]
end

# Define the constants
G = 6.6743 * 1e-8
μe = 2
Z = 6
ρ0 = 0.4 * Z^2 * μe
K = 1e13 * μe^(-5/3)

# Callback function
condition(u, r, integrator) = u[2] - ρ0
# Stop the integration when we hit the surface of the star, i.e
# when ρ becomes equal to ρ0

affect!(integrator) = terminate!(integrator)
cb = ContinuousCallback(condition, affect!)

# Define the ODE problem
u0 = @SVector [0, 1e7]  # initial values
rspan = [1e-10, 1e10]  # range of integration
params = [G, K, ρ0]  # parameters in the ODE
probNewt = ODEProblem(StellarStrunctureNewt, u0, rspan, params)

# Solve the ODEProblem
solNewt = solve(probNewt, callback=cb)

# Plot the solution
Msun = 1.989 * 1e33
plot(solNewt.t ./ 1e5, solNewt[1, :] ./Msun, xlabel=L"\textrm{radius [km]}", ylabel=L"$\textrm{Mass }[M_\odot]$", label=L"$\rho_c = 10^{7}$", legend=:bottomright)

# Now let's vary the central density and plot the final mass and radius
ρcs = 10 .^ LinRange(6, 9, 20)
Masses = zeros(size(ρcs))
radii = zeros(size(ρcs))

# We now solve the Newtonian structure equation for different central densities
# We want to save only the events and therefore we modify the callback function slightly
cb2 = ContinuousCallback(condition, affect!, save_positions=(false,true))
for (idx, ρc) in enumerate(ρcs)
	u0 = @SVector [0.0, ρc]
	probNewt = ODEProblem(StellarStrunctureNewt, u0, rspan, params)
	# since we want to save only the events we add few more otions in our solver
	solNewt = solve(probNewt, callback=cb, save_everystep=false, save_start=false, save_end=false)
	println(idx, "\t", ρc, "\t", solNewt.t[end], "\t" , solNewt.u[end])
	# Save the events
	Masses[idx] = solNewt.u[end][1]
	radii[idx] = solNewt.t[end]
end

# Now we plot the masses vs radii for different central densities
using ColorSchemes
scatter(radii ./ 1e5, Masses ./ Msun, zcolor=log10.(ρcs),
	xlabel=L"$\textrm{radius [km]}$",
	ylabel=L"$\textrm{Mass }[M_\odot]$",
	label=L"$\textrm{central density}$")