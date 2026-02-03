module ADPoisson

using LinearAlgebra
using Printf
using Statistics

include("types.jl")
include("core.jl")
include("boundary.jl")
include("problems.jl")
include("visualization.jl")
include("factory.jl")

export ProblemSpec, BoundaryConditions, TaylorBuffers3D, TaylorArrays3D
export SolverConfig, Solution
export default_cli_options
export solve, make_problem, boundary_from_prob
export exact_solution, source_term, dirichlet_bc
export plot_slice

end # module
