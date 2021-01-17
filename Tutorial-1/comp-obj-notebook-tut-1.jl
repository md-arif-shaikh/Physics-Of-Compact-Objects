### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ 17ce391a-58a1-11eb-0cc3-27ab49468914
begin
	using DifferentialEquations
	using ODEInterfaceDiffEq
	using ODEInterface
	using Printf
	using Plots
	using LaTeXStrings
	gr()
end

# ╔═╡ 50faac0c-58a2-11eb-0b45-913e36bec28c
md"**ODEs**

We solve the following set of first order odes which describes **structure of homologus star**

$$\frac{dp}{dm} = - \frac{m}{x^4}, \frac{dx}{dm} = \frac{t^b}{x^2p^a}, \frac{dt}{dm} = - \frac{p^{an}l}{x^4 t^{3+s+bn}}, \frac{dl}{dm} = A p^{a\lambda}t^{\nu-b\lambda}$$
"

# ╔═╡ abe0e024-58a1-11eb-3d89-45c25419d716
function stellar_structure!(du, u, params, m)
    a, b, n, s, ν, λ, A = params
    p, x, t, l = u
    du[1] = - (m / x^4)
    du[2] = (t^b / (x^2 * p^a))
    du[3] = - (p^(a * n) * l) / (x^4 * t^(3 + s + b * n))
    du[4] = A * p^(a * λ) * t^(ν - (b * λ))
end

# ╔═╡ 1863311a-58a3-11eb-18ca-bfa5e090d38e
md"**Callback function**

We want to know the value of $m$ where $p$ goes to zero. For this we need to use the callback function. First we set the condtion i.e $p$ goes to zero at the event. In practical we set a value below which say that p is sufficiently small."

# ╔═╡ fd2312b8-58a1-11eb-2e96-a94a16e0e8ca
condition(u, m, integrator) =  u[1] - 1e-4;

# ╔═╡ 6d9108e4-58a3-11eb-2bf7-e78d60fe87f4
md"We then need to specify what to do when the event is detected, i.e., the condition is satisfied. We here want to just stop the integration. This is done using the following"

# ╔═╡ 0bd12e44-58a2-11eb-0268-8b77422861ff
affect!(integrator) = terminate!(integrator);

# ╔═╡ 2a177554-58a5-11eb-21fb-f3f3112aa6d5
cb = ContinuousCallback(condition, affect!);

# ╔═╡ aa0096a8-58a3-11eb-0151-a333ecc27932
md"**Parameters**

Now to solve the set of odes we specify the values of the parameters as well as the initial values. First we set the values of the parameters"

# ╔═╡ c6f7d528-58a3-11eb-0286-67213e4687f7
begin
	a = 1.0
	b = 1.0
	n = 0.0
	s = 0.0
	ν = 4.0
	λ = 1.0
	A = 0.5
	params = [a, b, n, s, ν, λ, A]
end

# ╔═╡ fc9efeae-58a3-11eb-28ec-857b69eb26bb
md"**Initial values**

Now we set the initial values. If we take a carefull look at the odes, then we see that $\frac{dp}{dm}$ and $\frac{dx}{dm}$ have $x^4$ and $x^2$, respectively, in the denominator. We start our integration from $m=0$, where $x$ and $l$ are zero. However this makes the gardient diverge. Therefore we take small non-zero values for $x$ and m instead of zero"

# ╔═╡ a0fa6b00-58a4-11eb-15c8-1540d56d1781
begin
	p0 = 1.0
	x0 = 1e-3
	t0 = 1.0
	l0 = 0.0
	u0 = [p0, x0, t0, l0]
end

# ╔═╡ be3c61d2-58a4-11eb-0b28-e7092c004530
md"We also need specify the domain of integration"

# ╔═╡ cb8fa8bc-58a4-11eb-066f-915f0c57d6fa
mspan = (1e-9, 100.0)

# ╔═╡ dea60fd6-58a4-11eb-1656-6fe62b037ed4
md"**Solve ODE**

Now we are ready to set the ODE problem"

# ╔═╡ ec10984e-58a4-11eb-14c6-c5b8b047b955
problem1 = ODEProblem(stellar_structure!, u0, mspan, params);

# ╔═╡ 1141ef96-58a5-11eb-347e-b1a3fc705b0c
solution1 = solve(problem1, alg=RadauIIA5(), callback=cb)

# ╔═╡ 93817b8e-58aa-11eb-3116-4f4be1645bc3
solution1.t[end]

# ╔═╡ 9d037e8c-58aa-11eb-0cd2-818a684914fe
solution1.u[end]

# ╔═╡ 61558144-58a7-11eb-3b11-9dd48ad009a5
md"**Plot Solution**"

# ╔═╡ 3b294794-58a7-11eb-1380-df7639e7dd2f
begin
	plot(solution1, vars=(0, 1), lw = 2, label = L"p")
	plot!(solution1, vars=(0, 3), lw = 2, label = L"t")
	plot!(xaxis = (L"m"))
end

# ╔═╡ Cell order:
# ╠═17ce391a-58a1-11eb-0cc3-27ab49468914
# ╠═50faac0c-58a2-11eb-0b45-913e36bec28c
# ╠═abe0e024-58a1-11eb-3d89-45c25419d716
# ╟─1863311a-58a3-11eb-18ca-bfa5e090d38e
# ╠═fd2312b8-58a1-11eb-2e96-a94a16e0e8ca
# ╟─6d9108e4-58a3-11eb-2bf7-e78d60fe87f4
# ╠═0bd12e44-58a2-11eb-0268-8b77422861ff
# ╠═2a177554-58a5-11eb-21fb-f3f3112aa6d5
# ╟─aa0096a8-58a3-11eb-0151-a333ecc27932
# ╠═c6f7d528-58a3-11eb-0286-67213e4687f7
# ╟─fc9efeae-58a3-11eb-28ec-857b69eb26bb
# ╠═a0fa6b00-58a4-11eb-15c8-1540d56d1781
# ╟─be3c61d2-58a4-11eb-0b28-e7092c004530
# ╠═cb8fa8bc-58a4-11eb-066f-915f0c57d6fa
# ╟─dea60fd6-58a4-11eb-1656-6fe62b037ed4
# ╠═ec10984e-58a4-11eb-14c6-c5b8b047b955
# ╠═1141ef96-58a5-11eb-347e-b1a3fc705b0c
# ╠═93817b8e-58aa-11eb-3116-4f4be1645bc3
# ╠═9d037e8c-58aa-11eb-0cd2-818a684914fe
# ╟─61558144-58a7-11eb-3b11-9dd48ad009a5
# ╠═3b294794-58a7-11eb-1380-df7639e7dd2f
