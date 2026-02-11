using Test
using ADPoisson

@testset "solvers" begin
    config = SolverConfig(8, 8, 8, 4, 1.0e-4, 1000, 1.0e-6)
    prob, _ = ADPoisson.make_problem(config; alpha=1.0)
    bc = ADPoisson.boundary_from_prob(prob)

    u0 = ADPoisson.initialize_solution(config, prob).u
    f = zeros(eltype(u0), size(u0))
    ADPoisson.compute_source!(f, prob, config)

    # SOR
    sol_sor = ADPoisson.initialize_solution(config, prob)
    results_dir = joinpath(@__DIR__, "results")
    converged_sor, result_sor = ADPoisson.sor_solve!(sol_sor, f, bc, prob, config; output_dir=results_dir)
    @test converged_sor
    @test result_sor.iter <= config.max_steps

    # SSOR
    sol_ssor = ADPoisson.initialize_solution(config, prob)
    converged_ssor, result_ssor = ADPoisson.ssor_solve!(sol_ssor, f, bc, prob, config; output_dir=results_dir)
    @test converged_ssor
    @test result_ssor.iter <= config.max_steps

    # CG (no precond)
    sol_cg = ADPoisson.initialize_solution(config, prob)
    converged_cg, result_cg = ADPoisson.cg_solve!(sol_cg, f, bc, prob, config; precond=:none, output_dir=results_dir)
    @test converged_cg
    @test result_cg.iter <= config.max_steps

    # CG (SSOR precond)
    sol_cg_ssor = ADPoisson.initialize_solution(config, prob)
    converged_cg_ssor, result_cg_ssor = ADPoisson.cg_solve!(sol_cg_ssor, f, bc, prob, config; precond=:ssor, output_dir=results_dir)
    @test converged_cg_ssor
    @test result_cg_ssor.iter <= config.max_steps
end
