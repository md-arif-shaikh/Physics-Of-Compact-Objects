### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 17ce391a-58a1-11eb-0cc3-27ab49468914
using DifferentialEquations

# ╔═╡ 45f4efec-65f0-11eb-1a0a-3788d5efdc78
using StaticArrays

# ╔═╡ 402dc84e-65f1-11eb-2539-37f0fe257197
using Plots

# ╔═╡ 50faac0c-58a2-11eb-0b45-913e36bec28c
md"**ODEs**

We solve the following set of first order odes which describes **structure of homologus star**

$$\frac{dp}{dm} = - \frac{m}{x^4}, \frac{dx}{dm} = \frac{t^b}{x^2p^a}, \frac{dt}{dm} = - \frac{p^{an}l}{x^4 t^{3+s+bn}}, \frac{dl}{dm} = A p^{a\lambda}t^{\nu-b\lambda}$$
"

# ╔═╡ abe0e024-58a1-11eb-3d89-45c25419d716
function stellar_structure(u, params, m)
    a, b, n, s, ν, λ, A = params
    p, x, t, l = u
    dp = - (m / x^4)
    dx = (t^b / (x^2 * p^a))
    dt = - (p^(a * n) * l) / (x^4 * t^(3 + s + b * n))
    dl = A * p^(a * λ) * t^(ν - (b * λ))
	@SVector [dp, dx, dt, dl]
end

# ╔═╡ 1863311a-58a3-11eb-18ca-bfa5e090d38e
md"**Callback function**

We want to know the value of $m$ where $p$ goes to zero. For this we need to use the callback function. First we set the condtion i.e $p$ goes to zero at the event. In practical we set a value below which say that p is sufficiently small."

# ╔═╡ fd2312b8-58a1-11eb-2e96-a94a16e0e8ca
condition(u, m, integrator) =  u[1];

# ╔═╡ 6d9108e4-58a3-11eb-2bf7-e78d60fe87f4
md"We then need to specify what to do when the event is detected, i.e., the condition is satisfied. We here want to just stop the integration. This is done using the following"

# ╔═╡ 0bd12e44-58a2-11eb-0268-8b77422861ff
affect!(integrator) = terminate!(integrator);

# ╔═╡ 2a177554-58a5-11eb-21fb-f3f3112aa6d5
cb = ContinuousCallback(condition, affect!);

# ╔═╡ aa0096a8-58a3-11eb-0151-a333ecc27932
md"**Parameters**

Now to solve the set of odes we specify the values of the parameters as well as the initial values. First we set the values of the parameters. We will set A later."

# ╔═╡ c6f7d528-58a3-11eb-0286-67213e4687f7
begin
	a = 1.0
	b = 1.0
	n = 0.0
	s = 0.0
	ν = 4.0
	λ = 1.0
end

# ╔═╡ fc9efeae-58a3-11eb-28ec-857b69eb26bb
md"**Initial values**

Now we set the initial values. If we take a carefull look at the odes, then we see that $\frac{dp}{dm}$ and $\frac{dx}{dm}$ have $x^4$ and $x^2$, respectively, in the denominator. We start our integration from $m=0$, where $x$ and $l$ are zero. However this makes the gardient diverge. Therefore we take small non-zero values for $x$ instead of zero"

# ╔═╡ a0fa6b00-58a4-11eb-15c8-1540d56d1781
begin
	p0 = 1.0
	x0 = 1e-4
	t0 = 1.0
	l0 = 0.0
	u0 = @SVector [p0, x0, t0, l0]
end

# ╔═╡ be3c61d2-58a4-11eb-0b28-e7092c004530
md"We also need specify the domain of integration"

# ╔═╡ cb8fa8bc-58a4-11eb-066f-915f0c57d6fa
mspan = (0, 100.0)

# ╔═╡ dea60fd6-58a4-11eb-1656-6fe62b037ed4
md"**Solve ODE**

Now we are ready to set the ODE problem"

# ╔═╡ ec10984e-58a4-11eb-14c6-c5b8b047b955
begin
	params1 = [a, b, n, s, ν, λ, 0.5]
	problem1 = ODEProblem(stellar_structure, u0, mspan, params1);
end

# ╔═╡ 1141ef96-58a5-11eb-347e-b1a3fc705b0c
solution1 = solve(problem1, callback=cb)

# ╔═╡ 93817b8e-58aa-11eb-3116-4f4be1645bc3
solution1.t[end]

# ╔═╡ 9d037e8c-58aa-11eb-0cd2-818a684914fe
solution1.u[end]

# ╔═╡ 61558144-58a7-11eb-3b11-9dd48ad009a5
md"**Plot Solution**"

# ╔═╡ 9cdff58c-65f4-11eb-0704-47631de5b2a4
gr()

# ╔═╡ cf7cd588-65f3-11eb-0f95-97a0bc0c586d
begin
	colors = [:red, :green, :blue, :orange]
	labels = ["p", "x", "t", "m"]
	styles = [:solid, :dash, :dot, :dashdot]
end

# ╔═╡ 7c9d4d50-65f2-11eb-282b-f55027ec1222
begin
	global p = plot()
	for idx in 1:4
		global p = plot!(solution1, vars = (0, idx), color = colors[idx], ls = styles[idx], label = labels[idx])
	end
	plot!(xlim = (0, 100), ylim = (0, 1))
end

# ╔═╡ 4aec840e-58d7-11eb-2f48-93558ced3c13
md"Let us now find A for which p and t goes zero at the same value of m"

# ╔═╡ b835d0fa-65f5-11eb-34b7-8f18d3a7ef68
begin
	params2 = [a, b, n, s, ν, λ, 0.534]
	problem2 = ODEProblem(stellar_structure, u0, mspan, params2);
end

# ╔═╡ d1872548-65f5-11eb-02b5-c184e4141842
solution2 = solve(problem2, callback = cb)

# ╔═╡ e8eda41e-65f5-11eb-133b-e7140ded270b
begin
	global p2 = plot()
	for idx in 1:4
		global p2 = plot!(solution2, vars = (0, idx), color = colors[idx], ls = styles[idx], label = labels[idx])
	end
	plot!(xlim = (0, 12), ylim = (0, 1))
end

# ╔═╡ Cell order:
# ╠═17ce391a-58a1-11eb-0cc3-27ab49468914
# ╟─50faac0c-58a2-11eb-0b45-913e36bec28c
# ╠═45f4efec-65f0-11eb-1a0a-3788d5efdc78
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
# ╠═402dc84e-65f1-11eb-2539-37f0fe257197
# ╠═9cdff58c-65f4-11eb-0704-47631de5b2a4
# ╠═cf7cd588-65f3-11eb-0f95-97a0bc0c586d
# ╠═7c9d4d50-65f2-11eb-282b-f55027ec1222
# ╟─4aec840e-58d7-11eb-2f48-93558ced3c13
# ╠═b835d0fa-65f5-11eb-34b7-8f18d3a7ef68
# ╠═d1872548-65f5-11eb-02b5-c184e4141842
# ╠═e8eda41e-65f5-11eb-133b-e7140ded270b
