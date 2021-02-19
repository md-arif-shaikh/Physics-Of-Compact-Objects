import numpy as np
import matplotlib.pyplot as plt
import figStyle
from figStyle import plot_for

G = 6.67 * 1e-8
M = 1.989 * 1e33
R = 8.3 * 1e7
mue = 2.0


def Tc(rhoc, gamma, K, C):
    return ((mue/R) * (C * G * M**(2/3) * rhoc**(1/3)
                       - (K * rhoc**(gamma - 1)/mue**gamma)))


rhoc = 10.0 ** np.linspace(6, 9, 50)
C = 10.0
gammas = [4/3, 5/3]
gamma_labels = [r"$\frac{4}{3}$", r"$\frac{5}{3}$"]
ks = [1.24 * 1e15, 1e13]
k_labels = [r"$1.24 \times 10^{15}$", r"$10^{13}$"]

plot_for("PRD")
fig, ax = plt.subplots(nrows=1, ncols=1)
for gamma, k, gamma_label, k_label in zip(gammas, ks, gamma_labels, k_labels):
    ax.plot(rhoc, Tc(rhoc, gamma, k, C),
            label=fr"$\gamma$ =  {gamma_label}, $k$ = {k_label}, $C = 10$")
ax.legend(loc="lower right")
ax.grid(ls='dashed')
ax.set_xlabel(r"$\rho_c$")
ax.set_ylabel(r"$T_c$")
fig.savefig("Tc_vs_rhoc.pdf")
