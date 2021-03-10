### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 3f352894-81a6-11eb-3fd7-0333460c112f
using DifferentialEquations, StaticArrays

# ╔═╡ e88be666-81a7-11eb-12ed-67300d719779
using Plots, LaTeXStrings

# ╔═╡ 708f50ea-81a6-11eb-1cd8-477238b61afe
md"# Solve TOV equation for neutron star
Md Arif Shaikh, Postdoc, ICTS-TIFR, March 10 2021
"

# ╔═╡ 929b5bca-81a6-11eb-3793-dffff4210e82
md"## Newtonian Case:
First we solve for Newtonian stellar structure equations"

# ╔═╡ a9558cf2-81a6-11eb-0afc-9b1133f43885
function StellarStrunctureNewt(u, params, r)
	m, ρ = u
	G, K, ρ0 = params
	dpdρ = K * abs(ρ)^(2/3)*(5 - 4 * (ρ0/abs(ρ))^(1/3))/3
	dpdr = - G * m * ρ / r^2
	dρdr = dpdr/dpdρ
	dmdr = 4 * pi * r^2 * ρ
	@SVector [dmdr, dρdr]
end

# ╔═╡ b3e63b26-81a6-11eb-3536-cfc5c3cd70e2
md"Define the constants"

# ╔═╡ ba55efcc-81a6-11eb-3d16-43035b98aa9f
begin
	G = 6.6743 * 1e-8
	μe = 2
	Z = 6
	ρ0 = 0.4 * Z^2 * μe
	K = 1e13 * μe^(-5/3)
end

# ╔═╡ cf0ca384-81a6-11eb-3be4-dbd44e9b276c
md"Define the callback condition"

# ╔═╡ dbd5633a-81a6-11eb-2c33-e1555b27a642
condition(u, r, integrator) = u[2] - ρ0

# ╔═╡ f2042eca-81a6-11eb-12c3-abe61c337809
md" Stop the integration when we hit the surface of the star, i.e  when $\rho$ becomes equal to $\rho_0$"

# ╔═╡ e6737ce6-81a6-11eb-0e25-d3304d934ea7
affect!(integrator) = terminate!(integrator)

# ╔═╡ 06468d24-81a7-11eb-150f-1b0190fb30cc
cb = ContinuousCallback(condition, affect!)

# ╔═╡ 0f690b34-81a7-11eb-2a6a-b5e92506ba43
md"Now we define Newtonian ODE problem."

# ╔═╡ 73a5c95c-81a7-11eb-3011-2f51c93c10b2
begin
	u0 = @SVector [0, 1e8]  # initial values
	rspan = [1e-10, 1e10]  # range of integration
	params = [G, K, ρ0]  # parameters in the ODE
	probNewt = ODEProblem(StellarStrunctureNewt, u0, rspan, params)
end

# ╔═╡ b13fef66-81a7-11eb-38ff-491f793cde05
md"Solution"

# ╔═╡ b7f2a512-81a7-11eb-1e24-ddd9b8bd0f90
solNewt = solve(probNewt, callback=cb)

# ╔═╡ d50944bc-81a7-11eb-2f41-c1de1d753953
md"Plot the mass vs radius"

# ╔═╡ fced62ba-81a7-11eb-3e7f-4b95cf280ef5
begin
	Msun = 1.989 * 1e33
	plot(solNewt.t ./ 1e5, solNewt[1, :] ./Msun, xlabel=L"\textrm{radius [km]}", 			ylabel=L"$\textrm{Mass }[M_\odot]$", label=L"$\rho_c = 10^{8}$", 					legend=:bottomright)
end

# ╔═╡ 3001faa6-81a8-11eb-2c91-43731d360e7c
md"Now we solve for Mass of the star vs radius of the star for different values of the central density."

# ╔═╡ 534fd89a-81a8-11eb-156a-e3cf27733d96
begin
	ρcs = 10 .^ LinRange(6, 9, 25)
	MassesNewt = zeros(size(ρcs))
	radiiNewt = zeros(size(ρcs))
end

# ╔═╡ 64eea234-81a8-11eb-032a-15ffc27dcac3
md"We now solve the Newtonian structure equation for different central densities.
We want to save only the events and therefore we modify the callback function slightly."

# ╔═╡ 7932ec20-81a8-11eb-2a03-c70394f794d5
cbMassRadiiNewt = ContinuousCallback(condition, affect!, save_positions=(false,true))

# ╔═╡ 8e840ac6-81a8-11eb-32aa-2b7dfb553f1d
for (idx, ρc) in enumerate(ρcs)
	u0 = @SVector [0.0, ρc]
	probMassRadiiNewt = ODEProblem(StellarStrunctureNewt, u0, rspan, params)
	# since we want to save only the events we add few more otions in our solver
	solMassRadiiNewt = solve(probMassRadiiNewt, callback=cbMassRadiiNewt, 	save_everystep=false, save_start=false, save_end=false)
	# Save the events
	MassesNewt[idx] = solMassRadiiNewt.u[end][1]
	radiiNewt[idx] = solMassRadiiNewt.t[end]
end

# ╔═╡ f0c763b6-81a8-11eb-32ba-0304d36c7fe6
scatter(radiiNewt ./ 1e5, MassesNewt ./ Msun, zcolor=log10.(ρcs),
	xlabel=L"$\textrm{radius [km]}$",
	ylabel=L"$\textrm{Mass }[M_\odot]$",
	label=L"$\textrm{central density}$")

# ╔═╡ 36e57fc4-81a9-11eb-21fe-654f42968a86
md"## Solve TOV"

# ╔═╡ fd2bffbe-81a9-11eb-294b-3f09a5cb31aa
md"Write down the TOV equations with $G$ and $c$ restored.

$$m\to Gm/c^2$$

$$\rho\to G\rho/c^2$$

$$p\to Gp/c^4$$

So the TOV equations become

$$\frac{dp}{dr} = -\frac{Gm\rho}{r^2}\left[1 + \frac{1}{c^2}\frac{p}{\rho}\right]\left[1 + \frac{1}{c^2}\frac{4\pi r^3 p}{m}\right]\left[1 - \frac{2Gm}{c^2r}\right]^{-1}$$"

# ╔═╡ 1b29c084-81aa-11eb-1788-4946c721fc33
function tov(u, params, r)
    m, ρ = u
	G, C, K, ρ0 = params
    p = K * abs(ρ)^(5/3) * (1 - (ρ0/abs(ρ))^(1/3))
    dmdr = 4 * pi * r^2 * ρ
    dpdr = - (G * m * ρ / r^2) * (1 + (p/(C^2 * ρ))) * (1 + (4* pi * r^3 * p/(C^2 * m))) * (1 - (2 * G * m/(C^2 * r)))^(-1)
    dpdρ = K * abs(ρ)^(2/3) * (5 - 4 * (ρ0/abs(ρ))^(1/3)) / 3
    dρdr = dpdr/dpdρ
    
	@SVector [dmdr, dρdr]
end

# ╔═╡ eccb8058-81aa-11eb-3bd5-ffb497b470f4
md"First solve for a single central density and compare with the Newtonian solution"

# ╔═╡ 3728e3e8-81ab-11eb-0ffb-f1a9360d6ee9
begin
	c = 1e10
	u0TOV = @SVector [1e-10, u0[2]]
	paramsTOV = [G, c, K, ρ0]
	ProbTOV = ODEProblem(tov, u0TOV, rspan, paramsTOV)
end

# ╔═╡ aed24844-81ab-11eb-0a02-4f3c1af140ff
solveTOV = solve(ProbTOV, callback=cb)

# ╔═╡ 712b33ba-81ac-11eb-2b07-eb59d01b946d
md"Plot both Newtonian and TOV solution"

# ╔═╡ 80b123b4-81ac-11eb-2c41-a3e4826ef900
begin
	p = plot(xlabel=L"$r \textrm{ [km]}$", ylabel=L"$m [M_\odot]$")
	plot!(p, solNewt.t ./ 1e5, solNewt[1, :] ./Msun, label=L"$\textrm{Newtonian}$")
	plot!(p, solveTOV.t ./ 1e5, solveTOV[1, :] ./Msun, label=L"$\textrm{TOV}$", legend=:bottomright)
end

# ╔═╡ c0a059ba-81ad-11eb-0d10-c510ec43175e
md"As before we can vary the central density and plot mass vs radius"

# ╔═╡ d04ac76a-81ad-11eb-36b2-df50ae9fd5e1
begin
	MassesTOV = zeros(size(ρcs))
	radiiTOV = zeros(size(ρcs))
end

# ╔═╡ 049d4e8e-81ae-11eb-3563-e986cf2198cd
for (idx, ρc) in enumerate(ρcs)
	u0 = @SVector [1e-10, ρc]
	probMassRadiiTOV = ODEProblem(tov, u0, rspan, paramsTOV)
	# since we want to save only the events we add few more otions in our solver
	solMassRadiiTOV = solve(probMassRadiiTOV, callback=cbMassRadiiNewt, 	save_everystep=false, save_start=false, save_end=false)
	# Save the events
	MassesTOV[idx] = solMassRadiiTOV.u[end][1]
	radiiTOV[idx] = solMassRadiiTOV.t[end]
end

# ╔═╡ 4cfe94a8-81ae-11eb-2099-39022c741d8e
md"Plot mass vs radii for both Newtonian and TOV"

# ╔═╡ 5dce1b8c-81ae-11eb-17fc-dbd92026d51b
begin
	pTOVvsNewt = plot(xlabel=L"$r \textrm{ [km]}$", ylabel=L"$m [M_\odot]$")
	scatter!(pTOVvsNewt, radiiNewt ./ 1e5, MassesNewt ./ Msun, zcolor=log10.(ρcs), label=L"$\textrm{Newtonian}$")
	scatter!(pTOVvsNewt, radiiTOV ./ 1e5, MassesTOV ./ Msun, zcolor=log10.(ρcs), markershape =:square, label=L"$\textrm{TOV}$")
end

# ╔═╡ 18e64b4c-81af-11eb-0757-47bd91603d06
md"TODO: Solve TOV for Neutron star equation of state: SLy"

# ╔═╡ Cell order:
# ╟─708f50ea-81a6-11eb-1cd8-477238b61afe
# ╠═3f352894-81a6-11eb-3fd7-0333460c112f
# ╟─929b5bca-81a6-11eb-3793-dffff4210e82
# ╠═a9558cf2-81a6-11eb-0afc-9b1133f43885
# ╟─b3e63b26-81a6-11eb-3536-cfc5c3cd70e2
# ╠═ba55efcc-81a6-11eb-3d16-43035b98aa9f
# ╟─cf0ca384-81a6-11eb-3be4-dbd44e9b276c
# ╠═dbd5633a-81a6-11eb-2c33-e1555b27a642
# ╟─f2042eca-81a6-11eb-12c3-abe61c337809
# ╠═e6737ce6-81a6-11eb-0e25-d3304d934ea7
# ╠═06468d24-81a7-11eb-150f-1b0190fb30cc
# ╟─0f690b34-81a7-11eb-2a6a-b5e92506ba43
# ╠═73a5c95c-81a7-11eb-3011-2f51c93c10b2
# ╟─b13fef66-81a7-11eb-38ff-491f793cde05
# ╠═b7f2a512-81a7-11eb-1e24-ddd9b8bd0f90
# ╟─d50944bc-81a7-11eb-2f41-c1de1d753953
# ╠═e88be666-81a7-11eb-12ed-67300d719779
# ╠═fced62ba-81a7-11eb-3e7f-4b95cf280ef5
# ╟─3001faa6-81a8-11eb-2c91-43731d360e7c
# ╠═534fd89a-81a8-11eb-156a-e3cf27733d96
# ╟─64eea234-81a8-11eb-032a-15ffc27dcac3
# ╠═7932ec20-81a8-11eb-2a03-c70394f794d5
# ╠═8e840ac6-81a8-11eb-32aa-2b7dfb553f1d
# ╠═f0c763b6-81a8-11eb-32ba-0304d36c7fe6
# ╟─36e57fc4-81a9-11eb-21fe-654f42968a86
# ╟─fd2bffbe-81a9-11eb-294b-3f09a5cb31aa
# ╠═1b29c084-81aa-11eb-1788-4946c721fc33
# ╟─eccb8058-81aa-11eb-3bd5-ffb497b470f4
# ╠═3728e3e8-81ab-11eb-0ffb-f1a9360d6ee9
# ╠═aed24844-81ab-11eb-0a02-4f3c1af140ff
# ╟─712b33ba-81ac-11eb-2b07-eb59d01b946d
# ╠═80b123b4-81ac-11eb-2c41-a3e4826ef900
# ╟─c0a059ba-81ad-11eb-0d10-c510ec43175e
# ╟─d04ac76a-81ad-11eb-36b2-df50ae9fd5e1
# ╠═049d4e8e-81ae-11eb-3563-e986cf2198cd
# ╟─4cfe94a8-81ae-11eb-2099-39022c741d8e
# ╠═5dce1b8c-81ae-11eb-17fc-dbd92026d51b
# ╟─18e64b4c-81af-11eb-0757-47bd91603d06
