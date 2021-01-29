### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ ce069b56-6212-11eb-359b-d567223d2ca7
using UnitfulAstro, Unitful

# ╔═╡ dabacae4-6211-11eb-24e6-0f01be7534c2
md"Assignment3: Stellar Collapse"

# ╔═╡ 13fff7ca-6212-11eb-077c-7947cfdc39c8
md"## Problem 1."

# ╔═╡ 2d4e997a-6212-11eb-3671-c55e57a35aa6
md"The energy released in a core collapse supernovae could be estimated by the difference in the gravitational energy pre and post the collapse."

# ╔═╡ 79329336-6213-11eb-08b6-21052f99e3d8
rini = 1e4 * u"km"

# ╔═╡ ec734190-6213-11eb-245d-49b11cbb0cdc
rfinal = 10 * u"km"

# ╔═╡ f4789400-6213-11eb-23f9-690d151ca69b
mass = 1.4 * u"Msun"

# ╔═╡ fe84e7fa-6213-11eb-3024-b50806285bd5
G = 6.67430 * 1e-11 * u"m^3/kg/s^2"

# ╔═╡ 4cac0352-6212-11eb-2dd2-3370a67a7b01
function EngCoreCollapse(rini, rfinal, mass)
	G * mass^2 * ((1/rfinal) - (1/rini))
end

# ╔═╡ 2b5b3ef0-6214-11eb-1f2a-397bc6558586
Eng = uconvert(u"erg", EngCoreCollapse(rini, rfinal, mass))

# ╔═╡ 2ff74464-6216-11eb-3109-8b1189247127
md"## Problem 2"

# ╔═╡ 39137b50-6216-11eb-038d-932576307c40
md"We calculate the energy lost in the ejecta and photons and neutrons"

# ╔═╡ 50eb550e-6216-11eb-3b07-9f78ee7ea244
Mprogenitor = 10 * u"Msun"

# ╔═╡ cfb61ed2-6216-11eb-1422-b7ad3f259e40
Mejecta = Mprogenitor - mass

# ╔═╡ e21a9b20-6216-11eb-1115-c9e34718d03a
Vejecta = 10^4 * u"km/s"

# ╔═╡ f3ee18c2-6216-11eb-1c6b-6fdd191655e3
KE = uconvert(u"erg", (1/2)* Mejecta * Vejecta^2)

# ╔═╡ 1f000c50-6217-11eb-3b29-35e08d006c79
fractionKE = KE/Eng

# ╔═╡ 6261a06c-6217-11eb-0ed4-591915acb820
Luminosity = 2 * 10^8 * u"Lsun" 

# ╔═╡ ab933340-6217-11eb-2388-9db82ff944a7
Tvisible = 60 * 24 * 3600 * u"s" 

# ╔═╡ 90837fa6-6217-11eb-1801-839ebe2789a6
EMEnergy = uconvert(u"erg", Luminosity * u"s")

# ╔═╡ 3b24e44a-6218-11eb-396c-2932821c89e7
md"## Problem 3"

# ╔═╡ 42f1582a-6218-11eb-2999-bd0511099c52
Etotnutrino = Eng - KE - EMEnergy

# ╔═╡ 6644c924-6218-11eb-00ea-bda8f6c78709
Enutrino = uconvert(u"erg", 5 * 1e6 * u"eV")

# ╔═╡ 8899849c-6218-11eb-2773-d90d69839c9b
Nnutrino = Etotnutrino / Enutrino

# ╔═╡ 2985585e-6219-11eb-0a1d-13516cc0943a
dL = uconvert(u"m", 1e4 * u"pc")

# ╔═╡ 184cf092-6219-11eb-2d3a-75ef0c0f5580
Flux = Nnutrino / (4*pi*dL^2)

# ╔═╡ Cell order:
# ╠═dabacae4-6211-11eb-24e6-0f01be7534c2
# ╟─13fff7ca-6212-11eb-077c-7947cfdc39c8
# ╟─2d4e997a-6212-11eb-3671-c55e57a35aa6
# ╠═4cac0352-6212-11eb-2dd2-3370a67a7b01
# ╠═ce069b56-6212-11eb-359b-d567223d2ca7
# ╠═79329336-6213-11eb-08b6-21052f99e3d8
# ╠═ec734190-6213-11eb-245d-49b11cbb0cdc
# ╠═f4789400-6213-11eb-23f9-690d151ca69b
# ╠═fe84e7fa-6213-11eb-3024-b50806285bd5
# ╠═2b5b3ef0-6214-11eb-1f2a-397bc6558586
# ╟─2ff74464-6216-11eb-3109-8b1189247127
# ╟─39137b50-6216-11eb-038d-932576307c40
# ╠═50eb550e-6216-11eb-3b07-9f78ee7ea244
# ╠═cfb61ed2-6216-11eb-1422-b7ad3f259e40
# ╠═e21a9b20-6216-11eb-1115-c9e34718d03a
# ╠═f3ee18c2-6216-11eb-1c6b-6fdd191655e3
# ╠═1f000c50-6217-11eb-3b29-35e08d006c79
# ╠═6261a06c-6217-11eb-0ed4-591915acb820
# ╠═ab933340-6217-11eb-2388-9db82ff944a7
# ╠═90837fa6-6217-11eb-1801-839ebe2789a6
# ╟─3b24e44a-6218-11eb-396c-2932821c89e7
# ╠═42f1582a-6218-11eb-2999-bd0511099c52
# ╠═6644c924-6218-11eb-00ea-bda8f6c78709
# ╠═8899849c-6218-11eb-2773-d90d69839c9b
# ╠═2985585e-6219-11eb-0a1d-13516cc0943a
# ╠═184cf092-6219-11eb-2d3a-75ef0c0f5580
