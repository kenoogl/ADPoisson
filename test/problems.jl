using Test
using ADPoisson

const RUN_IMPL_TESTS = get(ENV, "ADPOISSON_IMPL_TEST", "0") == "1"

if RUN_IMPL_TESTS
    @testset "problems" begin
        alpha = 1.0
        @test exact_solution(0.25, 0.5, 0.0, alpha) ≈ alpha * sin(pi * 0.25) * sin(pi * 0.5)
        @test exact_solution(0.25, 0.5, 1.0, alpha) ≈ sin(pi * 0.25) * sin(pi * 0.5)
        @test source_term(0.1, 0.2, 0.3) == 0.0

        prob = ProblemSpec(1.0, 1.0, 1.0, alpha, source_term, dirichlet_bc)
        bc = boundary_from_prob(prob)
        @test bc.g_zlo(0.25, 0.5) ≈ alpha * sin(pi * 0.25) * sin(pi * 0.5)
        @test bc.g_zhi(0.25, 0.5) ≈ sin(pi * 0.25) * sin(pi * 0.5)
        @test bc.g_xlo(0.2, 0.3) == 0.0
        @test bc.g_yhi(0.2, 0.3) == 0.0
    end
end
