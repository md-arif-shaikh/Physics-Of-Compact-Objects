using DifferentialEquations

function stellar_structure(m, u, du)
    du[1] = - m / x^4
    du[2] = - t^b / (x^2 * p^a)
    du[3] = - p^(a * n) * l / (x^4 * t^(3 + s + b*n))
    du[4] = A * p^(a * lambda) * t^(nu - b * lambda)
end
