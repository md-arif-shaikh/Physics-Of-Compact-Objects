### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ b213914c-6fa0-11eb-0ba6-950195900c41
using DifferentialEquations, CairoMakie, StaticArrays

# ╔═╡ 73341cde-6fab-11eb-1c14-6b12fc1e85b8
using Unitful, UnitfulAstro

# ╔═╡ 853114f6-6fa0-11eb-39a8-edca713a7c1a
md"### Assignment 5

**Course on Physics of Compact Objects**

Instructor: Prof. P. Ajith

Author: Md Arif Shaikh, ICTS-TIFR, Feb 15, 2021"

# ╔═╡ a5da3e0a-6fa0-11eb-2518-ef262ae935d9
md"### Problem 2
Numerically solve the stellar structure equations

$\frac{dm}{dr} = 4\pi r^2 \rho, \frac{dP}{dr} = - \frac {G m \rho}{r^2}$

assuming the equation of state of non-relativistic degenerate electron gas *including Coulomb correction*s:

$ P(\rho) = K \rho^{5/3}\left[ 1 − (\rho_0 /\rho)^{1/3}\right],$

where $K \simeq 10^{13} \mu_e^{-5/3}$ (cgs units) and $\rho_0 \simeq 0.4 Z^2 \mu_e$ g/cm$^3$. Assume the atomic numeber to be $Z = 6$ (Carbon) an the mean molecular weight per free electron to be $\mu_e = 2$. The obvious boundary conditions
are $m(r = 0) = 0$ and $P(r = 0) = P(\rho_c)$, where $\rho_c$ is the central density. Compute the mass-radius of white dwarfs for central densities ranging $\rho_c = 10^6 − 10^9$ g/cm$^3$ . Compare the mass-radius relation with the same computed earlier
neglecting Coulomb corrections.
"

# ╔═╡ fe34ad36-6fa0-11eb-1187-11d611ecdbb7
begin
	const G = 6.6743 * 1e-8
	const μe = 2.0
	const K = 1e13 * μe^(-5.0 / 3.0) 
	const Z = 2.0
	const ρ0 = 0.4 * Z^2 * μe
end

# ╔═╡ 0f31b1e2-6fa1-11eb-3a90-19e6e49800c4
function stellarStructure(u, params, r)
    m, ρ = u
    ρ0, G, K = params
    dmdr = 4.0 * pi * r^2 * abs(ρ)
    dρdr = - (3.0 * G * m * abs(ρ)^(1.0 /3.0 )) / (K * r^2 * (5.0 - 4.0 * abs((ρ0 / ρ))^(1.0 / 3.0)))
    @SVector [dmdr, dρdr]
end

# ╔═╡ 1c2ef3ca-6fa1-11eb-3946-453e872da2ca
begin
	condition(u, r, integrator) =  u[2] - ρ0;
	affect!(integrator) = terminate!(integrator);
	cb = ContinuousCallback(condition, affect!);
end

# ╔═╡ 371db5f2-6fa1-11eb-3546-cf1fa66b48cf
begin
	rspan = [1e-7, 1e10]
	params = [ρ0, G, K]
	ρc = 1e8

	u0 = @SVector [0.0, ρc]
	problem = ODEProblem(stellarStructure, u0, rspan, params)
	solution = solve(problem, callback = cb)
end

# ╔═╡ b115c610-6fa6-11eb-3bf9-eb593e3f7d78
begin
	fig = Figure(resolution = (800, 450), fonts = "Times Roman", fontsize = 12)
	ax = fig[1, 1] = Axis(fig, xlabel = "radius [cgs]", ylabel = "density [cgs]")
	lines!(solution.t, solution[2, :], color = :red)
	fig
end

# ╔═╡ 83c00ae4-6fa7-11eb-121a-a76e67b869e3
begin
	masses = Float64[]
	radii = Float64[]
	ρcs = LinRange(1.e6, 1.e8, 100)
	for ρc in ρcs
		u0 = @SVector [0.0, ρc]
		problem = ODEProblem(stellarStructure, u0, rspan, params)
		sol = solve(problem, callback = cb)
		push!(radii, sol.t[end])
		push!(masses, sol[1, :][end])
	end
end

# ╔═╡ 7bdf068c-6fab-11eb-31fd-b5b610054ccb
Msun = 1.989 * 1e33 

# ╔═╡ eb191c4a-6fab-11eb-18a9-fdad808bddc5
Rsun = 696340 * 1e5 

# ╔═╡ 0eb3412c-6fab-11eb-12d0-29d1a35e11ba
begin
	figMR = Figure(resolution = (800, 450), fonts = "Times Roman", fontsize = 12)
	axMR = figMR[1, 1] = Axis(figMR, 
		ylabel = "radius [km]", 
		xlabel = "mass [Msun]")
	lines!(masses / Msun, radii / 1e5, color = :red)
	xlims!(axMR, [masses[1], masses[end]]/Msun)
	ylims!(axMR, [radii[end], radii[1]]/1e5)
	figMR
end

# ╔═╡ Cell order:
# ╟─853114f6-6fa0-11eb-39a8-edca713a7c1a
# ╟─a5da3e0a-6fa0-11eb-2518-ef262ae935d9
# ╠═b213914c-6fa0-11eb-0ba6-950195900c41
# ╠═fe34ad36-6fa0-11eb-1187-11d611ecdbb7
# ╠═0f31b1e2-6fa1-11eb-3a90-19e6e49800c4
# ╠═1c2ef3ca-6fa1-11eb-3946-453e872da2ca
# ╠═371db5f2-6fa1-11eb-3546-cf1fa66b48cf
# ╠═b115c610-6fa6-11eb-3bf9-eb593e3f7d78
# ╠═83c00ae4-6fa7-11eb-121a-a76e67b869e3
# ╠═73341cde-6fab-11eb-1c14-6b12fc1e85b8
# ╠═7bdf068c-6fab-11eb-31fd-b5b610054ccb
# ╠═eb191c4a-6fab-11eb-18a9-fdad808bddc5
# ╠═0eb3412c-6fab-11eb-12d0-29d1a35e11ba
