using Test
using ADPoisson
using Printf
import ADPoisson: laplacian!, laplacian4!, taylor_step!, apply_bc!, boundary_from_prob, grid_spacing
import ADPoisson: exact_solution, source_term, dirichlet_bc

const RUN_IMPL_TESTS = get(ENV, "ADPOISSON_IMPL_TEST", "0") == "1"

function l2_norm_interior(a, config)
    s = 0.0
    gx = (size(a, 1) - config.nx) ÷ 2
    gy = (size(a, 2) - config.ny) ÷ 2
    gz = (size(a, 3) - config.nz) ÷ 2
    @inbounds for k in gz+1:gz+config.nz
        for j in gy+1:gy+config.ny
            for i in gx+1:gx+config.nx
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
    gx = (size(sol.u, 1) - config.nx) ÷ 2
    gy = (size(sol.u, 2) - config.ny) ÷ 2
    gz = (size(sol.u, 3) - config.nz) ÷ 2
    @inbounds for k in 1:config.nz
        for j in 1:config.ny
            for i in 1:config.nx
                x = sol.x[i]
                y = sol.y[j]
                z = sol.z[k]
                uex = exact_solution(x, y, z, prob.alpha)
                u = sol.u[i + gx, j + gy, k + gz]
                s += (u - uex)^2
                sref += uex^2
            end
        end
    end
    return sqrt(s * dx * dy * dz) / sqrt(sref * dx * dy * dz)
end

if RUN_IMPL_TESTS
    @testset "laplacian!" begin
        config = SolverConfig(8, 8, 8, 2, 1e-3, 1, 1e-10)
        ghost = 2
        u = fill(2.5, config.nx + 2 * ghost, config.ny + 2 * ghost, config.nz + 2 * ghost)
        Lu = similar(u)
        laplacian!(Lu, u, config)
        gx = ghost
        @test maximum(abs.(Lu[gx+1:end-gx, gx+1:end-gx, gx+1:end-gx])) == 0.0
    end

    @testset "taylor_step!" begin
        config = SolverConfig(6, 6, 6, 2, 1e-3, 1, 1e-10)
        ghost = 2
        u = fill(1.0, config.nx + 2 * ghost, config.ny + 2 * ghost, config.nz + 2 * ghost)
        f = fill(0.25, config.nx + 2 * ghost, config.ny + 2 * ghost, config.nz + 2 * ghost)
        next = similar(u)
        taylor_step!(next, u, f, 0, config)
        gx = ghost
        @test maximum(abs.(next[gx+1:end-gx, gx+1:end-gx, gx+1:end-gx] .+ 0.25)) < 1e-12
    end

    @testset "convergence_order" begin
        alpha = 1.0
        ns = [16, 32, 64]
        errs = Float64[]

        for n in ns
            config = SolverConfig(n, n, n, 2, 1e-3, 1, 1e-10)
            prob = ProblemSpec(1.0, 1.0, 1.0, alpha, source_term, dirichlet_bc)
            bc = boundary_from_prob(prob)
            ghost = 2
            u = zeros(n + 2 * ghost, n + 2 * ghost, n + 2 * ghost)

            dx, dy, dz = grid_spacing(config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
            @inbounds for k in ghost+1:ghost+n
                z = (k - (ghost + 0.5)) * dz
                for j in ghost+1:ghost+n
                    y = (j - (ghost + 0.5)) * dy
                    for i in ghost+1:ghost+n
                        x = (i - (ghost + 0.5)) * dx
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
end

@testset "laplacian4!" begin
    ns = [16, 32, 64]
    errs = Float64[]

    for n in ns
        config = SolverConfig(n, n, n, 2, 1e-3, 1, 1e-10)
        dx = 1.0 / n
        dy = 1.0 / n
        dz = 1.0 / n
        u = zeros(n + 4, n + 4, n + 4)
        Lu = similar(u)

        @inbounds for k in 1:n+4
            z = (k - 2.5) * dz
            for j in 1:n+4
                y = (j - 2.5) * dy
                for i in 1:n+4
                    x = (i - 2.5) * dx
                    u[i, j, k] = sin(2 * pi * x) * sin(2 * pi * y) * sin(2 * pi * z)
                end
            end
        end

        laplacian4!(Lu, u, config)
        s = 0.0
        @inbounds for k in 3:n+2
            z = (k - 2.5) * dz
            for j in 3:n+2
                y = (j - 2.5) * dy
                for i in 3:n+2
                    x = (i - 2.5) * dx
                    exact = -12 * pi^2 * sin(2 * pi * x) * sin(2 * pi * y) * sin(2 * pi * z)
                    d = Lu[i, j, k] - exact
                    s += d * d
                end
            end
        end
        push!(errs, sqrt(s / (n^3)))
    end

    order1 = log2(errs[1] / errs[2])
    order2 = log2(errs[2] / errs[3])
    @test order1 > 3.5
    @test order2 > 3.5
end

@testset "apply_bc ghost2" begin
    config = SolverConfig(6, 6, 6, 2, 1e-3, 1, 1e-10)
    g_lo = 1.2
    g_hi = -0.7
    g0 = (a, b) -> 0.0
    bc = BoundaryConditions((y, z) -> g_lo, (y, z) -> g_hi, g0, g0, g0, g0)
    ghost = 2
    i_lo = ghost + 1
    i_hi = ghost + config.nx
    j = ghost + 1
    k = ghost + 1

    u = zeros(config.nx + 2 * ghost, config.ny + 2 * ghost, config.nz + 2 * ghost)
    u[i_lo, j, k] = 0.4
    u[i_lo + 1, j, k] = -0.1
    u[i_lo + 2, j, k] = 0.7
    u[i_lo + 3, j, k] = 1.5
    u[i_hi, j, k] = -0.3
    u[i_hi - 1, j, k] = 0.6
    u[i_hi - 2, j, k] = -0.8
    u[i_hi - 3, j, k] = 2.2

    apply_bc!(u, bc, 0, config; Lx=1.0, Ly=1.0, Lz=1.0, order=:spec)
    @test u[i_lo - 1, j, k] ≈ 2 * g_lo - u[i_lo, j, k]
    @test u[i_lo - 2, j, k] ≈ 2 * g_lo - u[i_lo + 1, j, k]
    @test u[i_hi + 1, j, k] ≈ 2 * g_hi - u[i_hi, j, k]
    @test u[i_hi + 2, j, k] ≈ 2 * g_hi - u[i_hi - 1, j, k]

    apply_bc!(u, bc, 1, config; Lx=1.0, Ly=1.0, Lz=1.0, order=:spec)
    @test u[i_lo - 1, j, k] ≈ -u[i_lo, j, k]
    @test u[i_lo - 2, j, k] ≈ -u[i_lo + 1, j, k]
    @test u[i_hi + 1, j, k] ≈ -u[i_hi, j, k]
    @test u[i_hi + 2, j, k] ≈ -u[i_hi - 1, j, k]

    u .= 0
    u[i_lo, j, k] = 0.4
    u[i_lo + 1, j, k] = -0.1
    u[i_lo + 2, j, k] = 0.7
    u[i_lo + 3, j, k] = 1.5
    u[i_hi, j, k] = -0.3
    u[i_hi - 1, j, k] = 0.6
    u[i_hi - 2, j, k] = -0.8
    u[i_hi - 3, j, k] = 2.2

    apply_bc!(u, bc, 0, config; Lx=1.0, Ly=1.0, Lz=1.0, order=:high)
    @test u[i_lo - 1, j, k] ≈ (128 / 35) * g_lo - 4 * u[i_lo, j, k] + 2 * u[i_lo + 1, j, k] -
                              (4 / 5) * u[i_lo + 2, j, k] + (1 / 7) * u[i_lo + 3, j, k]
    @test u[i_lo - 2, j, k] ≈ (128 / 7) * g_lo - 30 * u[i_lo, j, k] + 20 * u[i_lo + 1, j, k] -
                              9 * u[i_lo + 2, j, k] + (12 / 7) * u[i_lo + 3, j, k]
    @test u[i_hi + 1, j, k] ≈ (128 / 35) * g_hi - 4 * u[i_hi, j, k] + 2 * u[i_hi - 1, j, k] -
                              (4 / 5) * u[i_hi - 2, j, k] + (1 / 7) * u[i_hi - 3, j, k]
    @test u[i_hi + 2, j, k] ≈ (128 / 7) * g_hi - 30 * u[i_hi, j, k] + 20 * u[i_hi - 1, j, k] -
                              9 * u[i_hi - 2, j, k] + (12 / 7) * u[i_hi - 3, j, k]

    apply_bc!(u, bc, 1, config; Lx=1.0, Ly=1.0, Lz=1.0, order=:high)
    @test u[i_lo - 1, j, k] ≈ -4 * u[i_lo, j, k] + 2 * u[i_lo + 1, j, k] -
                              (4 / 5) * u[i_lo + 2, j, k] + (1 / 7) * u[i_lo + 3, j, k]
    @test u[i_lo - 2, j, k] ≈ -30 * u[i_lo, j, k] + 20 * u[i_lo + 1, j, k] -
                              9 * u[i_lo + 2, j, k] + (12 / 7) * u[i_lo + 3, j, k]
end

@testset "solver_error" begin
    alpha = 1.0
    ns = get(ENV, "ADPOISSON_FULL_TEST", "0") == "1" ? [16, 32, 64] : [16, 32]
    errs = Float64[]
    do_plot = get(ENV, "ADPOISSON_TEST_PLOT", "0") == "1"
    defaults = default_cli_options()
    use_safe_dt = get(ENV, "ADPOISSON_TEST_USE_SAFE_DT", "1") == "1"

    bc_order = :high
    for n in ns
        prob = ProblemSpec(1.0, 1.0, 1.0, alpha, source_term, dirichlet_bc)
        dt = defaults["dt"]
        if use_safe_dt
            safety = bc_order === :high ? 0.25 : 0.5
            dt_safe = safety / (3 * n^2)
            dt = min(dt, dt_safe)
        end
        max_steps = defaults["max_steps"]
        println("test config:")
        @printf("  n=%d dt=%.3e max_steps=%d epsilon=%.3e bc_order=%s use_safe_dt=%s\n",
                n, dt, max_steps, 1e-6, string(bc_order), string(use_safe_dt))
        config = SolverConfig(n, n, n, 4, dt, max_steps, 1e-6)
        sol = solve(config, prob; bc_order=bc_order)
        if do_plot
            plot_slice(sol, prob, config)
        end
        err = l2_error_solution(sol, prob, config)
        push!(errs, err)
    end

    @test errs[end] <= 1e-3
    if length(errs) >= 2
        order = log2(errs[1] / errs[2])
        @test order > 1.5
    end
end
