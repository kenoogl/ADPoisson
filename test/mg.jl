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

    sol3 = ADPoisson.initialize_solution(config, prob)
    ADPoisson.apply_bc!(sol3.u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=:spec)
    r3_0 = ADPoisson.compute_residual_norm!(r, sol3.u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    ADPoisson.vcycle!(sol3.u, f, bc, config, prob;
                      nu1=1, nu2=1, dt_scale=2.0, M=2, bc_order=:spec)
    r3_1 = ADPoisson.compute_residual_norm!(r, sol3.u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    @test r3_1 <= r3_0

    sol4 = ADPoisson.initialize_solution(config, prob)
    ADPoisson.apply_bc!(sol4.u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=:spec)
    r4_0 = ADPoisson.compute_residual_norm!(r, sol4.u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    ADPoisson.vcycle!(sol4.u, f, bc, config, prob;
                      nu1=1, nu2=1, dt_scale=2.0, M=2,
                      correction_mode=:correction_taylor,
                      corr_M=2, corr_dt_scale=1.0, corr_steps=1,
                      bc_order=:spec)
    r4_1 = ADPoisson.compute_residual_norm!(r, sol4.u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    @test r4_1 <= r4_0
    @test all(isfinite, sol4.u)
end
