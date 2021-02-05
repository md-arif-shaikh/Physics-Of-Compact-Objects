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

# ╔═╡ 3c0653ec-67ba-11eb-1bce-e5ec4c92989e
using Unitful, UnitfulAstro

# ╔═╡ c7ffa252-660e-11eb-381b-2f2ceec373c1
md"## white dwarfs"

# ╔═╡ de5358ce-67c0-11eb-1c87-49efb42f32d7
md"### Course on Compact Objects: Assignment 4

Instructor: Prof. P. Ajith
Author: Md Arif Shaikh, Postdoc, ICTS-TIFR
Date: Feb 5, 2021"

# ╔═╡ 73e9d7ec-6793-11eb-3bca-bd8d40ef03d6
md"### Problem 1

Show that the pressure exerted by a gas of particle with isotropic momentum distribution $n(p)$ is given by

$P = \frac{1}{3}\int_0^\infty p v_p n(p) dp$

where $v_p$ is the velocity associated with momentum $p$.
"

# ╔═╡ b71be1a4-6793-11eb-0012-dba495d16df2
md"Let us consider the pressure in the $z$ direction. The component of velocity of along the $z$ direction would be given by $v \cos\theta$ where $\theta$ is the angle of the velocity vector with the $z$ axis.


$v_z = v \cos\theta$

and therefore,

$p_z = p \cos\theta$

Now the number of particles in a volume element $d^3p$ in momentum space with momentum between $p$ and $p +dp$ is

$d^3p = p^2 dp\sin\theta d\theta d\phi$

Now the rate of change of momentum is given by

$\frac{p_z}{\Delta z/v_z} = \frac{p_z v_z}{\Delta z} = \frac{p v \cos^2\theta}{\Delta z}$

So the pressure per particle is 

$\frac{p v \cos^2\theta}{\Delta z \Delta x \Delta y} = \frac{p v \cos^2\theta}{V}$

So the total pressure would be given by 

$P = \int_0^\infty \frac{p v \cos^2\theta}{V} N(p) dp$

where $N(p)$ is the total number of particles with momentum between $p$ and $p+dp$ and with angle $\theta$ to $\theta + d\theta$ which is $2\pi p^2 dp\sin\theta d\theta$ (which is obtained by integrating $d^3p$ over $\phi$)

Thus,

$P = \int_0^\infty \frac{p v \cos^2\theta}{V} 2\pi \sin\theta d\theta p^2 dp$

or

$P = \int_0^\pi \cos^2\theta d\theta\sin\theta \int_0^\infty \frac{p v }{V} 2\pi p^2 dp = \frac{1}{3} \int_0^\infty \frac{p v }{V} 4\pi p^2 dp = \frac{1}{3}\int_0^\infty p v n(p) dp$


where $n(p) = 4\pi p^2/V$
"

# ╔═╡ 6dab1184-6798-11eb-01e2-af1bc7dbdd1d
md"### Problem 2

Argue why we are justified in using a 'cold' degenerate equation of state to describe a white dwarf with a temperature
$T \sim 10^4$K (Hint: Show that the degeneracy parameter $\mu/kT \gg 0$, where $\mu$ is the chemical potential and $k$ the Boltzmann
constant. The density of the white dwarf is $\sim 10^6 g/cm^3$ and the chemical potential $\sim$ the Fermi energy. Assume that
$\mu_e = 2$)."

# ╔═╡ 93d8e972-6799-11eb-0149-ed6bd5c87857
md"### Problem 3

Derive Lane-Emden equation"

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

# ╔═╡ add3d718-6799-11eb-0a12-398c88b47440
md"### Problem 4"

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
	ψs = []
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
		push!(ψs, ψfrac)
	end
end

# ╔═╡ 1d4b843c-66a1-11eb-188c-871b357bfdd7
ξ1s

# ╔═╡ 72c6e67c-66a1-11eb-2171-852cad67d86f
pfrac

# ╔═╡ 913a247c-66b8-11eb-296a-b92827dff6c8
data = DataFrame(n = ns, ξ1 = ξ1s, ψ1 = ψs)

# ╔═╡ 774b990a-66b9-11eb-2f71-ef6e4b08ec6f
md"### Problem 5

Radius and Mass of white dwarf"

# ╔═╡ 27a64b22-679d-11eb-0763-6d5931f2d0bb
md"Mass of the star is given by

$M = \int_0^\infty 4\pi\rho r^2 dr$

now,

$r = \xi  \left(\frac{4 \pi  G {\rho_c}^{1-\frac{1}{n}}}{K (n+1)}\right)^{-\frac{1}{2}}$

$\rho=\rho_c \theta^n$

$\theta^n = - \frac{1}{\xi^2}\frac{d}{d\xi}\left(\xi^2 \frac{d\theta}{d\xi}\right)$

Thus,

$M = 4\pi \rho_c\left(\frac{4 \pi  G {\rho_c}^{1-\frac{1}{n}}}{K (n+1)}\right)^{-3/2}\int_0^{\xi_1} d\left(\xi^2 \frac{d\theta}{d\xi}\right) = 4\pi \rho_c\left(\frac{4 \pi  G {\rho_c}^{1-\frac{1}{n}}}{K (n+1)}\right)^{-3/2}\left(\xi^2 \frac{d\theta}{d\xi}\right)_{\xi_1}$

or 

$M = 4\pi \rho_c^{(3-n)/2n}\left(\frac{4 \pi G}{K (n+1)}\right)^{-3/2}\left(\xi^2 \frac{d\theta}{d\xi}\right)_{\xi_1}$

Now, the radius would be given by value of $r$ at $\xi = \xi_1$ which is 

$R_\star = \xi_1  \left(\frac{4 \pi  G {\rho_c}^{1-\frac{1}{n}}}{K (n+1)}\right)^{-\frac{1}{2}}$

From this we can express $\rho_c$ in terms of $R_\star$ and $\xi_1$ as

$\rho_c = \left(\frac{K (n+1) \xi_1 ^2}{4\pi G R_\star^2}\right)^{\frac{n}{n-1}}$

So,

$\boxed{M = 4\pi R_\star^{(3-n)/(1-n)} \left(\frac{K (n+1)}{4\pi G}\right)^{n/(n-1)} \xi_1^{-(3-n)/(1-n)} \xi_1^2\left.\frac{d\theta}{d\xi}\right|_{\xi_1}}$
"

# ╔═╡ d299d78a-67b7-11eb-2bf5-192e423764d7
md"### Problem 6

Compute mass and radius of white dwarfs
"

# ╔═╡ ea8073d8-67b7-11eb-3a59-99769f0574f9
md"From problem 4, we get that for $n=3/2$, $\xi_1 \approx 3.655$. Now given that

$\rho_c = 10^6 - 10^9 g/cm^3, K \approx 10^{13} \mu_e^{-5/3}, \mu_e = 2.$"

# ╔═╡ a1b993b2-66bf-11eb-3df9-ddd3a3363f1b
function MassWhiteDwarf(R, K, G, ξ1, ψ1, n)
	4 * pi * R^((3. -n)/(1. - n))*((K*(n+1))/(4*pi*G))^(n/(n-1))*ξ1^(-(3. - n)/(1. -n))*ξ1^2 * abs(ψ1)
end

# ╔═╡ bb649190-67b9-11eb-012d-2bc5d384ddba
function Rstar(G, ξ1, ρc, K, n)
	return ξ1 * (4 * pi * G * ρc^(1.0 - (1.0 / n))/(K * (n + 1)))^(-1.0 / 2.0)
end

# ╔═╡ 056ca994-67ba-11eb-2459-5bf59f7db5c6
ξ1 = data.ξ1[2]

# ╔═╡ 2f4ab814-67ba-11eb-34e6-db97da22b6ad
ψ1 = data.ψ1[2]

# ╔═╡ 47d0d1c0-67ba-11eb-2d35-f32177efc1ab
G = u"GMsun"/u"Msun"

# ╔═╡ 6b1ef04e-67ba-11eb-3c96-d15015059f76
μ = 2.0

# ╔═╡ dbb73808-67be-11eb-12d0-916378cf0c7f
n = 3/2

# ╔═╡ 58b6b5ea-67ba-11eb-3bca-1d89747d69bd
K = 1e13 * μ^(-5.0/3.0) * u"dyn/cm^2"/(u"g/cm^3")^(1+(1/n))

# ╔═╡ 7205031c-67ba-11eb-1e45-75ac2b20c76f
ρcs = [1e6, 1e7, 1e8, 1e9] * u"g/cm^3"

# ╔═╡ 7233ae8c-67bb-11eb-13dc-716bb11add02
ρcs[2]

# ╔═╡ be271b86-67ba-11eb-2828-517b586bc4ad
r1 = uconvert(u"km", Rstar(G, ξ1, ρcs[1], K, n))

# ╔═╡ 3e871518-67bd-11eb-100b-9dce3872fa40
mass1 = MassWhiteDwarf(r1, K, G, ξ1, ψ1, n)

# ╔═╡ 93838510-67bd-11eb-1ea0-7fcb8db20c6f
uconvert(u"Msun", mass1)

# ╔═╡ 171b5cca-67c0-11eb-3c51-c525076f9966
begin
	masses = []
	radii = []
	
	for ρc in ρcs
		radius = uconvert(u"km", Rstar(G, ξ1, ρc, K, n))
		mass = uconvert(u"Msun", MassWhiteDwarf(radius, K, G, ξ1, ψ1, n))
		push!(masses, mass)
		push!(radii, radius)
	end
end

# ╔═╡ 7299b18c-67c0-11eb-1feb-2533c617d02d
massRadius = DataFrame(central_density = ρcs, radius = radii, mass = masses)

# ╔═╡ Cell order:
# ╟─c7ffa252-660e-11eb-381b-2f2ceec373c1
# ╟─de5358ce-67c0-11eb-1c87-49efb42f32d7
# ╟─73e9d7ec-6793-11eb-3bca-bd8d40ef03d6
# ╟─b71be1a4-6793-11eb-0012-dba495d16df2
# ╟─6dab1184-6798-11eb-01e2-af1bc7dbdd1d
# ╟─93d8e972-6799-11eb-0149-ed6bd5c87857
# ╟─b98c5b72-6611-11eb-3001-71f3dddb7f23
# ╟─e1480d1c-660e-11eb-1c29-5decbd0775cd
# ╟─add3d718-6799-11eb-0a12-398c88b47440
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
# ╟─27a64b22-679d-11eb-0763-6d5931f2d0bb
# ╟─d299d78a-67b7-11eb-2bf5-192e423764d7
# ╟─ea8073d8-67b7-11eb-3a59-99769f0574f9
# ╠═a1b993b2-66bf-11eb-3df9-ddd3a3363f1b
# ╠═bb649190-67b9-11eb-012d-2bc5d384ddba
# ╠═056ca994-67ba-11eb-2459-5bf59f7db5c6
# ╠═2f4ab814-67ba-11eb-34e6-db97da22b6ad
# ╠═3c0653ec-67ba-11eb-1bce-e5ec4c92989e
# ╠═47d0d1c0-67ba-11eb-2d35-f32177efc1ab
# ╠═6b1ef04e-67ba-11eb-3c96-d15015059f76
# ╠═dbb73808-67be-11eb-12d0-916378cf0c7f
# ╠═58b6b5ea-67ba-11eb-3bca-1d89747d69bd
# ╠═7205031c-67ba-11eb-1e45-75ac2b20c76f
# ╠═7233ae8c-67bb-11eb-13dc-716bb11add02
# ╠═be271b86-67ba-11eb-2828-517b586bc4ad
# ╠═3e871518-67bd-11eb-100b-9dce3872fa40
# ╠═93838510-67bd-11eb-1ea0-7fcb8db20c6f
# ╠═171b5cca-67c0-11eb-3c51-c525076f9966
# ╠═7299b18c-67c0-11eb-1feb-2533c617d02d
