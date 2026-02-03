# problems.jl

"""
    exact_solution(x, y, z, alpha)

Analytical solution for the verification case (f = 0, steady state).
"""
function exact_solution(x::Real, y::Real, z::Real, alpha::Real)
    sxy = sin(pi * x) * sin(pi * y)
    denom = sinh(sqrt(2) * pi)
    return sxy / denom * (sinh(sqrt(2) * pi * z) + alpha * sinh(sqrt(2) * pi * (1 - z)))
end

"""
    source_term(x, y, z)

Source term for the Poisson equation. Default is zero.
"""
source_term(x::Real, y::Real, z::Real) = zero(promote_type(typeof(x), typeof(y), typeof(z)))

"""
    dirichlet_bc(x, y, z, alpha)

Dirichlet boundary condition for the default verification case.
"""
function dirichlet_bc(x::Real, y::Real, z::Real, alpha::Real)
    if z == 0
        return alpha * sin(pi * x) * sin(pi * y)
    elseif z == 1
        return sin(pi * x) * sin(pi * y)
    else
        return zero(promote_type(typeof(x), typeof(y), typeof(z), typeof(alpha)))
    end
end

"""
    boundary_from_prob(prob::ProblemSpec)

Build face-wise Dirichlet boundary functions from the problem definition.

Face function signatures:
- x-faces: g_xlo(y, z), g_xhi(y, z)
- y-faces: g_ylo(x, z), g_yhi(x, z)
- z-faces: g_zlo(x, y), g_zhi(x, y)
"""
function boundary_from_prob(prob::ProblemSpec)
    f = prob.dirichlet
    alpha = prob.alpha
    Lx = prob.Lx
    Ly = prob.Ly
    Lz = prob.Lz

    g_xlo = (y, z) -> f(0, y, z, alpha)
    g_xhi = (y, z) -> f(Lx, y, z, alpha)
    g_ylo = (x, z) -> f(x, 0, z, alpha)
    g_yhi = (x, z) -> f(x, Ly, z, alpha)
    g_zlo = (x, y) -> f(x, y, 0, alpha)
    g_zhi = (x, y) -> f(x, y, Lz, alpha)

    return BoundaryConditions(g_xlo, g_xhi, g_ylo, g_yhi, g_zlo, g_zhi)
end
