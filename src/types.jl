# types.jl

"""
    ProblemSpec

Defines physical parameters and boundary/source definitions for the Poisson problem.
"""
struct ProblemSpec{T<:Real,Fsrc,Fbc}
    Lx::T
    Ly::T
    Lz::T
    alpha::T
    source::Fsrc
    dirichlet::Fbc
end

"""
    BoundaryConditions

Dirichlet boundary condition functions on each face.
"""
struct BoundaryConditions{Fxlo,Fxhi,Fylo,Fyhi,Fzlo,Fzhi}
    g_xlo::Fxlo
    g_xhi::Fxhi
    g_ylo::Fylo
    g_yhi::Fyhi
    g_zlo::Fzlo
    g_zhi::Fzhi
end

"""
    TaylorBuffers3D

Buffers for ping-pong Taylor coefficient computation.
"""
struct TaylorBuffers3D{T<:Real}
    bufA::Array{T,3}
    bufB::Array{T,3}
    acc::Array{T,3}
end

"""
    TaylorArrays3D

Optional full storage of Taylor coefficients for verification.
"""
struct TaylorArrays3D{T<:Real}
    U::Array{T,4}
end

"""
    SolverConfig

Numerical parameters for the solver.
"""
struct SolverConfig{T<:Real}
    nx::Int
    ny::Int
    nz::Int
    M::Int
    dt::T
    max_steps::Int
    epsilon::T
end

const DEFAULT_CLI_OPTS = Dict{String,Any}(
    "nx" => 16,
    "ny" => 16,
    "nz" => 16,
    "M" => 10,
    "dt" => 2e-4,
    "max_steps" => 10000,
    "epsilon" => 1e-10,
    "alpha" => 1.0,
    "bc_order" => "spec",
    "solver" => "taylor",
    "cg_precond" => "none",
)

"""
    default_cli_options()

Return a fresh dictionary of CLI default options.
"""
function default_cli_options()
    return Dict(k => v for (k, v) in DEFAULT_CLI_OPTS)
end

"""
    Solution

Results container for the computed solution.
"""
struct Solution{T<:Real}
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    u::Array{T,3}
    t::T
    iter::Int
end
