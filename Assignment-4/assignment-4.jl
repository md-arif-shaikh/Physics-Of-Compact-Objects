### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ fa72a230-6612-11eb-2c8c-41f072db62c8
using DifferentialEquations

# ╔═╡ 19d62798-6613-11eb-3432-fd10668c8447
using StaticArrays

# ╔═╡ dc2f69a6-6613-11eb-39a0-456ca3cfcafe
using Plots; gr()

# ╔═╡ 5f050fe8-6614-11eb-38e1-3b878533a01c
using LaTeXStrings

# ╔═╡ ab527e9c-66a2-11eb-3531-8b666bb82451
using DataFrames

# ╔═╡ c7ffa252-660e-11eb-381b-2f2ceec373c1
md"## white dwarfs"

# ╔═╡ b98c5b72-6611-11eb-3001-71f3dddb7f23
md"### Lane-Emden equation"

# ╔═╡ e1480d1c-660e-11eb-1c29-5decbd0775cd
md"$\frac{dm(r)}{dr} = 4\pi r^2\rho(r), \frac{dp(r)}{dr}=-\frac{Gm(r)\rho(r)}{r^2}$. These two equations can be combined to give

$\frac{1}{r^2}\frac{d}{dr}\left(\frac{r^2}{\rho(r)}\frac{dp(r)}{dr}\right) = -4\pi G \rho(r)$

This can be further written in terms of a dimensionless form using

$\rho = \rho_c \theta^n, r = a\xi, \Gamma = 1 + \frac{1}{n}$

Then we get the following equation

$\frac{1}{\xi^2}\frac{d}{d\xi}\left(\xi^2 \frac{d\theta}{d\xi}\right) = -\theta^n$

This the Lane-Emden equation.
"

# ╔═╡ ce805c54-6611-11eb-1c6b-ff532c7973b6
md"### Solve Lane-Emden Equation"

# ╔═╡ dad519b6-6611-11eb-20f6-df1d57ee4fcb
md"To solve the lane-emden equation we write it as two first order odes by defining $\psi = \frac{d\theta}{d\xi}$

$\frac{d\theta}{d\xi} = \psi, \frac{d\psi}{d\xi} = - \frac{2\psi}{\xi} - \theta^n$

With the following initial values

$\theta(0) = 1, \psi(0) = 0$"

# ╔═╡ 9ca0c614-6612-11eb-2dce-1b8e502da696
function laneEmden(u, n, ξ)
	θ, ψ = u
	dθ = ψ
	dψ = - (2 * ψ / ξ) - θ^n
	@SVector [dθ, dψ]
end

# ╔═╡ 38699210-6613-11eb-2fae-1502e73646d3
u0 = @SVector [1., 0.]

# ╔═╡ 59092ad0-6613-11eb-3cca-0be9b3b034f6
ξspan1 = (1e-10, 10)

# ╔═╡ 71d03ac2-6613-11eb-31a9-376c57d17fbe
problem1 = ODEProblem(laneEmden, u0, ξspan1, 1.)

# ╔═╡ b568cdc6-6613-11eb-3bf6-0f0fce82c1d2
sol1 = solve(problem1)

# ╔═╡ 033df3a0-6614-11eb-17e8-d175880db41d
plot(sol1, vars = (0, 1), xlabel = L"\xi", ylabel = L"\theta", label = L"n = 1")

# ╔═╡ b386900e-6629-11eb-1a5e-7b5c82483b88
md"Solutions for different integer $n$"

# ╔═╡ f0587620-662a-11eb-1f28-5f6a0938d300


# ╔═╡ 3eb68340-662b-11eb-1465-1b4e43a9750d
condition(u, ξ, integrator) =  u[1];

# ╔═╡ c22a1caa-662b-11eb-2fe0-99791c25951d
affect!(integrator) = terminate!(integrator);

# ╔═╡ cdfd4ed8-662b-11eb-28e3-21f8f05045e8
cb = ContinuousCallback(condition, affect!);

# ╔═╡ eb305878-6629-11eb-224c-0775dbfbbe01
ξspan = (1e-10, 20)

# ╔═╡ c2483324-6629-11eb-096f-892e2660ddbe
begin
	p = plot();
	for n in 1:4
		problem = ODEProblem(laneEmden, u0, ξspan, n)
		sol = solve(problem, callback = cb)
		plot!(p, sol, vars = (0, 1), lw = 2, label = L"n = %$n")
	end
	plot!(xlabel = L"\xi", ylabel = L"\theta", ylim = (0, 1))
	p
end

# ╔═╡ 75295168-6640-11eb-2f4d-79dbafd7e74c
md"for fractional $n$ the solver would throw error as negative value of $\theta$ would give complex value for $\theta^n$. Callback function does not work for some reason. Therefore we manually check at every step whether $\theta$ is negative. We stop whenever $\theta < \epsilon$. where $\epsilon$ is a very small number"

# ╔═╡ c57492b8-6640-11eb-1f31-29cc2440c6e0
begin
	pfrac = plot();
	ns = 0.5:1.0:4.5
	ξ1s = []
	for n in ns
		θfracsol = []
		ψfracsol = []
		ξs = []
		θfrac = 1.
		ψfrac = 0.
		ξ = 1e-8
		dξ = 1e-4
		while θfrac > dξ
			push!(θfracsol, θfrac)
			push!(ψfracsol, ψfrac)
			push!(ξs, ξ)
			θfracini = θfrac
			ψfracini = ψfrac
			u0frac = @SVector [θfracini, ψfracini]
			ξspanfrac = (ξ, ξ+dξ)
			problem = ODEProblem(laneEmden, u0frac, ξspanfrac, n)
			sol = solve(problem)
			θfrac = sol.u[end][1]
			ψfrac = sol.u[end][2]
			println(sol)
			ξ += dξ
		end
		plot!(pfrac, ξs, θfracsol, lw = 2, xlabel = L"\xi", ylabel = L"\theta", label = L"n = %$n")
		push!(ξ1s, (ξ + dξ/2.))
	end
end

# ╔═╡ 1d4b843c-66a1-11eb-188c-871b357bfdd7
ξ1s

# ╔═╡ 72c6e67c-66a1-11eb-2171-852cad67d86f
pfrac

# ╔═╡ 913a247c-66b8-11eb-296a-b92827dff6c8
DataFrame(n = ns, ξ1 = ξ1s)

# ╔═╡ 774b990a-66b9-11eb-2f71-ef6e4b08ec6f
md"### Mass of white dwarf"

# ╔═╡ a1b993b2-66bf-11eb-3df9-ddd3a3363f1b


# ╔═╡ Cell order:
# ╟─c7ffa252-660e-11eb-381b-2f2ceec373c1
# ╟─b98c5b72-6611-11eb-3001-71f3dddb7f23
# ╟─e1480d1c-660e-11eb-1c29-5decbd0775cd
# ╟─ce805c54-6611-11eb-1c6b-ff532c7973b6
# ╟─dad519b6-6611-11eb-20f6-df1d57ee4fcb
# ╠═fa72a230-6612-11eb-2c8c-41f072db62c8
# ╠═19d62798-6613-11eb-3432-fd10668c8447
# ╠═9ca0c614-6612-11eb-2dce-1b8e502da696
# ╠═38699210-6613-11eb-2fae-1502e73646d3
# ╠═59092ad0-6613-11eb-3cca-0be9b3b034f6
# ╠═71d03ac2-6613-11eb-31a9-376c57d17fbe
# ╠═b568cdc6-6613-11eb-3bf6-0f0fce82c1d2
# ╠═dc2f69a6-6613-11eb-39a0-456ca3cfcafe
# ╠═5f050fe8-6614-11eb-38e1-3b878533a01c
# ╠═033df3a0-6614-11eb-17e8-d175880db41d
# ╟─b386900e-6629-11eb-1a5e-7b5c82483b88
# ╠═f0587620-662a-11eb-1f28-5f6a0938d300
# ╠═3eb68340-662b-11eb-1465-1b4e43a9750d
# ╠═c22a1caa-662b-11eb-2fe0-99791c25951d
# ╠═cdfd4ed8-662b-11eb-28e3-21f8f05045e8
# ╠═eb305878-6629-11eb-224c-0775dbfbbe01
# ╠═c2483324-6629-11eb-096f-892e2660ddbe
# ╟─75295168-6640-11eb-2f4d-79dbafd7e74c
# ╠═c57492b8-6640-11eb-1f31-29cc2440c6e0
# ╠═1d4b843c-66a1-11eb-188c-871b357bfdd7
# ╠═72c6e67c-66a1-11eb-2171-852cad67d86f
# ╠═ab527e9c-66a2-11eb-3531-8b666bb82451
# ╠═913a247c-66b8-11eb-296a-b92827dff6c8
# ╟─774b990a-66b9-11eb-2f71-ef6e4b08ec6f
# ╠═a1b993b2-66bf-11eb-3df9-ddd3a3363f1b
