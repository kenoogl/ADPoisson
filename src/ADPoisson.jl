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
include("sor.jl")
include("cg.jl")
include("mg.jl")

export ProblemSpec, BoundaryConditions, TaylorBuffers3D, TaylorArrays3D
export SolverConfig, Solution
export default_cli_options
export solve, solve_with_runtime, make_problem, boundary_from_prob
export exact_solution, source_term, dirichlet_bc
export plot_slice
export sor_solve, sor_solve!, sor_solve_with_runtime
export ssor_solve, ssor_solve!, ssor_solve_with_runtime
export cg_solve, cg_solve!, cg_solve_with_runtime
export pseudo_mg_correction!, restrict_full_weighting!, prolong_trilinear!, vcycle!, two_level_mg_correction!

end # module
