import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from figStyle import plot_for


def stellar_structure(m, y, params):
    a, b, n, s, nu, lam, A = params
    p, x, t, ell = y

    dpdm = - (m / x**4)
    dxdm = (t**b / (x**2 * p**a))
    dtdm = - (p**(a*n) * ell) / (x**4 * t**(3 + s + b * n))
    dldm = A * p**(a * lam) * t**(nu - (b * lam))

    return [dpdm, dxdm, dtdm, dldm]


# In the first part of the problem we set A to 1.
# The initial values at the center of the start could be read off from eq 1.2
a = 1.0
b = 1.0
n = 0.0
s = 0.0
lam = 1.0
nu = 4.0
A = 0.534
params = [a, b, n, s, nu, lam, A]
mspan = [0, 100.0]
m_eval = np.arange(mspan[0], mspan[-1], 0.001)
y0 = [1.0, 1e-7, 1.0, 0.0]


# We want to find where p(m) goes to zero. So we define an event function
def hit_surface(m, y, params):
    return y[0] - 1e-3


hit_surface.terminal = True
hit_surface.direction = -1

sol1 = solve_ivp(stellar_structure, mspan, y0, t_eval=m_eval, args=(params,),
                 events=hit_surface)
mstar = sol1.t_events
# print(mstar)
print(sol1.t[-1], sol1.y[0][-1], sol1.y[3][-1])
# Let's also plot p(m) vs m to see where it becomes zero.
fig, ax = plt.subplots(nrows=1, ncols=1)
variables = ["$p$", "$x$", "$t$", "$l$"]
styles = ["solid", "dashed", "dotted", "dashdot"]
for idx, (style, label) in enumerate(zip(styles, variables)):
    ax.plot(sol1.t, sol1.y[idx], ls=style, label=rf"{label}")

ax.set_xlim(0, round(sol1.t[-1]) + 1)
ax.set_ylim(0, 2)
ax.set_xlabel(r"$m$")
ax.grid()
ax.legend()
fig.savefig(f"./p_vs_m_A_{A}.pdf")
