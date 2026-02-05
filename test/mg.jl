using Test
using ADPoisson

@testset "mg" begin
    config = SolverConfig(8, 8, 8, 4, 1.0e-4, 1000, 1.0e-6)
    prob, _ = ADPoisson.make_problem(config; alpha=1.0)
    bc = ADPoisson.boundary_from_prob(prob)

    sol = ADPoisson.initialize_solution(config, prob)
    f = zeros(eltype(sol.u), config.nx + 2, config.ny + 2, config.nz + 2)
    ADPoisson.compute_source!(f, prob, config)

    r = similar(sol.u)
    ADPoisson.apply_bc!(sol.u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=:spec)
    r0 = ADPoisson.compute_residual_norm!(r, sol.u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)

    ADPoisson.pseudo_mg_correction!(sol.u, f, bc, config, prob; dt_scale=2.0, M=2, bc_order=:spec)
    r1 = ADPoisson.compute_residual_norm!(r, sol.u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)

    @test r1 <= r0

    sol2 = ADPoisson.initialize_solution(config, prob)
    ADPoisson.apply_bc!(sol2.u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=:spec)
    r2_0 = ADPoisson.compute_residual_norm!(r, sol2.u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    ADPoisson.two_level_mg_correction!(sol2.u, f, bc, config, prob;
                                       nu1=1, nu2=1, dt_scale=2.0, M=2, bc_order=:spec)
    r2_1 = ADPoisson.compute_residual_norm!(r, sol2.u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    @test r2_1 <= r2_0
end
