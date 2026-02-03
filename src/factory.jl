# factory.jl

"""
    make_problem(config; alpha=1, Lx=1, Ly=1, Lz=1)

Construct the default Poisson verification problem.
"""
function make_problem(config::SolverConfig; alpha=1.0, Lx=1.0, Ly=1.0, Lz=1.0)
    prob = ProblemSpec(Lx, Ly, Lz, alpha, source_term, dirichlet_bc)
    bc = boundary_from_prob(prob)
    return prob, bc
end
