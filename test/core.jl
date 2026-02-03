using Test
using ADPoisson
import ADPoisson: laplacian!, taylor_step!, apply_bc!, boundary_from_prob, grid_spacing
import ADPoisson: exact_solution, source_term, dirichlet_bc

function l2_norm_interior(a, config)
    s = 0.0
    @inbounds for k in 2:config.nz+1
        for j in 2:config.ny+1
            for i in 2:config.nx+1
                s += a[i, j, k]^2
            end
        end
    end
    return sqrt(s)
end

function l2_norm_interior_weighted(a, config, prob)
    dx, dy, dz = grid_spacing(config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    return l2_norm_interior(a, config) * sqrt(dx * dy * dz)
end

function l2_error_solution(sol, prob, config)
    dx = prob.Lx / config.nx
    dy = prob.Ly / config.ny
    dz = prob.Lz / config.nz
    s = 0.0
    sref = 0.0
    @inbounds for k in 1:config.nz
        for j in 1:config.ny
            for i in 1:config.nx
                x = sol.x[i]
                y = sol.y[j]
                z = sol.z[k]
                uex = exact_solution(x, y, z, prob.alpha)
                u = sol.u[i + 1, j + 1, k + 1]
                s += (u - uex)^2
                sref += uex^2
            end
        end
    end
    return sqrt(s * dx * dy * dz) / sqrt(sref * dx * dy * dz)
end

@testset "laplacian!" begin
    config = SolverConfig(8, 8, 8, 2, 1e-3, 1.0, 1e-10)
    u = fill(2.5, config.nx + 2, config.ny + 2, config.nz + 2)
    Lu = similar(u)
    laplacian!(Lu, u, config)
    @test maximum(abs.(Lu[2:end-1, 2:end-1, 2:end-1])) == 0.0
end

@testset "taylor_step!" begin
    config = SolverConfig(6, 6, 6, 2, 1e-3, 1.0, 1e-10)
    u = fill(1.0, config.nx + 2, config.ny + 2, config.nz + 2)
    f = fill(0.25, config.nx + 2, config.ny + 2, config.nz + 2)
    next = similar(u)
    taylor_step!(next, u, f, 0, config)
    @test maximum(abs.(next[2:end-1, 2:end-1, 2:end-1] .+ 0.25)) < 1e-12
end

@testset "convergence_order" begin
    alpha = 1.0
    ns = [16, 32, 64]
    errs = Float64[]

    for n in ns
        config = SolverConfig(n, n, n, 2, 1e-3, 1.0, 1e-10)
        prob = ProblemSpec(1.0, 1.0, 1.0, alpha, source_term, dirichlet_bc)
        bc = boundary_from_prob(prob)
        u = zeros(n + 2, n + 2, n + 2)

        dx, dy, dz = grid_spacing(config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
        @inbounds for k in 2:n+1
            z = (k - 1.5) * dz
            for j in 2:n+1
                y = (j - 1.5) * dy
                for i in 2:n+1
                    x = (i - 1.5) * dx
                    u[i, j, k] = exact_solution(x, y, z, alpha)
                end
            end
        end

        apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=:high)
        Lu = similar(u)
        laplacian!(Lu, u, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)

        err = l2_norm_interior_weighted(Lu, config, prob) /
              l2_norm_interior_weighted(u, config, prob)
        push!(errs, err)
    end

    order1 = log2(errs[1] / errs[2])
    order2 = log2(errs[2] / errs[3])
    @test order1 > 1.5
    @test order2 > 1.5
end

@testset "solver_error" begin
    alpha = 1.0
    ns = get(ENV, "ADPOISSON_FULL_TEST", "0") == "1" ? [16, 32, 64] : [16, 32]
    errs = Float64[]

    for n in ns
        prob = ProblemSpec(1.0, 1.0, 1.0, alpha, source_term, dirichlet_bc)
        dt = 0.1 / (3 * n^2)
        config = SolverConfig(n, n, n, 4, dt, 0.5, 1e-6)
        sol = solve(config, prob; bc_order=:high)
        err = l2_error_solution(sol, prob, config)
        push!(errs, err)
    end

    @test errs[end] <= 1e-3
    if length(errs) >= 2
        order = log2(errs[1] / errs[2])
        @test order > 1.5
    end
end
