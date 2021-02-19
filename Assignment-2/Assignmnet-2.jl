# Problem 3

using CairoMakie


const G = 6.67 * 1e-8
const M = 1.989 * 1e33
const R = 8.3 * 1e7
const μe = 2.0

function Tc(ρc, K, γ, C)
	(μe/R) * (C * G * M^(2/3) * ρc^(1/3) - (K * ρc^(γ - 1)/μe^γ))
end


ρc = 10 .^LinRange(6, 9, 50)
C = 10
ks = [1.24 * 1e15, 1e13]
γs = [4/3, 5/3]
colors = [:orange, :cyan]
fig = Figure(resolution=(600, 400), font = "Times New Roman", fontsize = 12)
ax = fig[1,1] = Axis(fig, xlabel = "rho_c", ylabel = "T_c")
lins = [lines!(ρc, Tc.(ρc, kval, γval, C), color = color, linewidth = 2) for (kval, γval, color) in zip(ks, γs, colors)]
legs = ["gamma = $γval, K = $kval" for (kval, γval) in zip(ks, γs)]
leg = fig[1,1] = Legend(fig, lins, legs)
fig
