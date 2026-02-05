# sor.jl

"""
    sor_solve_with_runtime(prob, config; omega=1.0, output_dir="results", bc_order=:spec)

Solve Poisson equation using RB-SOR and return (Solution, runtime_s).
"""
function sor_solve_with_runtime(prob::ProblemSpec, config::SolverConfig;
                                omega::Real=1.0, output_dir::AbstractString="results",
                                bc_order=:spec)
    sol = initialize_solution(config, prob)
    bc = boundary_from_prob(prob)
    f = zeros(eltype(sol.u), config.nx + 2, config.ny + 2, config.nz + 2)
    compute_source!(f, prob, config)
    omega_t = convert(eltype(sol.u), omega)
    _, sol_out, runtime = sor_solve_with_runtime!(sol, f, bc, prob, config;
                                                  omega=omega_t, output_dir=output_dir,
                                                  bc_order=bc_order)
    return sol_out, runtime
end

"""
    ssor_solve_with_runtime(prob, config; omega=1.0, output_dir="results", bc_order=:spec)

Solve Poisson equation using RBSSOR (symmetric RB-SOR) and return (Solution, runtime_s).
"""
function ssor_solve_with_runtime(prob::ProblemSpec, config::SolverConfig;
                                 omega::Real=1.0, output_dir::AbstractString="results",
                                 bc_order=:spec)
    sol = initialize_solution(config, prob)
    bc = boundary_from_prob(prob)
    f = zeros(eltype(sol.u), config.nx + 2, config.ny + 2, config.nz + 2)
    compute_source!(f, prob, config)
    omega_t = convert(eltype(sol.u), omega)
    _, sol_out, runtime = ssor_solve_with_runtime!(sol, f, bc, prob, config;
                                                   omega=omega_t, output_dir=output_dir,
                                                   bc_order=bc_order)
    return sol_out, runtime
end

"""
    ssor_solve(prob, config; omega=1.0, output_dir="results", bc_order=:spec)

Solve Poisson equation using RBSSOR and return Solution.
"""
function ssor_solve(prob::ProblemSpec, config::SolverConfig;
                    omega::Real=1.0, output_dir::AbstractString="results",
                    bc_order=:spec)
    sol_out, _ = ssor_solve_with_runtime(prob, config;
                                         omega=omega, output_dir=output_dir,
                                         bc_order=bc_order)
    return sol_out
end

"""
    sor_solve(prob, config; omega=1.0, output_dir="results", bc_order=:spec)

Solve Poisson equation using RB-SOR and return Solution.
"""
function sor_solve(prob::ProblemSpec, config::SolverConfig;
                   omega::Real=1.0, output_dir::AbstractString="results",
                   bc_order=:spec)
    sol_out, _ = sor_solve_with_runtime(prob, config;
                                        omega=omega, output_dir=output_dir,
                                        bc_order=bc_order)
    return sol_out
end

"""
    sor_solve!(sol, f, bc, prob, config; omega=1.0, output_dir="results", bc_order=:spec)

RB-SOR solver with residual history output.
"""
function sor_solve!(sol::Solution{T}, f::Array{T,3}, bc::BoundaryConditions,
                    prob::ProblemSpec, config::SolverConfig;
                    omega::T=one(T), output_dir::AbstractString="results",
                    bc_order=:spec) where {T<:Real}
    converged, result, _ = sor_solve_with_runtime!(sol, f, bc, prob, config;
                                                   omega=omega, output_dir=output_dir,
                                                   bc_order=bc_order)
    return converged, result
end

function sor_solve_with_runtime!(sol::Solution{T}, f::Array{T,3}, bc::BoundaryConditions,
                                 prob::ProblemSpec, config::SolverConfig;
                                 omega::T=one(T), output_dir::AbstractString="results",
                                 bc_order=:spec) where {T<:Real}
    u = sol.u
    r = similar(u)
    u_exact = exact_solution_array(sol, prob, config)

    dx, dy, dz = grid_spacing(config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    inv_dx2 = one(T) / (dx * dx)
    inv_dy2 = one(T) / (dy * dy)
    inv_dz2 = one(T) / (dz * dz)
    diag = 2 * (inv_dx2 + inv_dy2 + inv_dz2)

    apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
    r0 = compute_residual_norm!(r, u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    denom = max(r0, one(r0))

    history = IOBuffer()
    println(history, "# step err_l2 res_l2")
    err_l2 = l2_error_exact_precomputed(u, u_exact, prob, config)
    @printf(history, "%d %.6e %.6e\n", 0, err_l2, r0 / denom)

    converged = false
    iter = 0
    t_start = time()
    for step in 1:config.max_steps
        apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        sor_sweep_color!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 0)
        apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        sor_sweep_color!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 1)

        res = compute_residual_norm!(r, u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz) / denom
        err_l2 = l2_error_exact_precomputed(u, u_exact, prob, config)
        @printf(history, "%d %.6e %.6e\n", step, err_l2, res)
        iter = step
        if res <= config.epsilon
            converged = true
            break
        end
    end
    runtime = time() - t_start

    output_dir = string(output_dir)
    isdir(output_dir) || mkpath(output_dir)
    tag = "nx$(config.nx)_ny$(config.ny)_nz$(config.nz)_steps$(iter)"
    history_path = joinpath(output_dir, "history_sor_$(tag).txt")
    open(history_path, "w") do io
        write(io, String(take!(history)))
    end

    result = Solution(sol.x, sol.y, sol.z, sol.u, zero(T), iter)
    return converged, result, runtime
end

function ssor_solve_with_runtime!(sol::Solution{T}, f::Array{T,3}, bc::BoundaryConditions,
                                  prob::ProblemSpec, config::SolverConfig;
                                  omega::T=one(T), output_dir::AbstractString="results",
                                  bc_order=:spec) where {T<:Real}
    u = sol.u
    r = similar(u)
    u_exact = exact_solution_array(sol, prob, config)

    dx, dy, dz = grid_spacing(config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    inv_dx2 = one(T) / (dx * dx)
    inv_dy2 = one(T) / (dy * dy)
    inv_dz2 = one(T) / (dz * dz)
    diag = 2 * (inv_dx2 + inv_dy2 + inv_dz2)

    apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
    r0 = compute_residual_norm!(r, u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    denom = max(r0, one(r0))

    history = IOBuffer()
    println(history, "# step err_l2 res_l2")
    err_l2 = l2_error_exact_precomputed(u, u_exact, prob, config)
    @printf(history, "%d %.6e %.6e\n", 0, err_l2, r0 / denom)

    converged = false
    iter = 0
    t_start = time()
    for step in 1:config.max_steps
        # RBSSOR: forward R→B, backward B→R, forward B→R, backward R→B
        apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        sor_sweep_color!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 0)
        sor_sweep_color!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 1)

        apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        sor_sweep_color_backward!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 1)
        sor_sweep_color_backward!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 0)

        apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        sor_sweep_color!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 1)
        sor_sweep_color!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 0)

        apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        sor_sweep_color_backward!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 0)
        sor_sweep_color_backward!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 1)

        res = compute_residual_norm!(r, u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz) / denom
        err_l2 = l2_error_exact_precomputed(u, u_exact, prob, config)
        @printf(history, "%d %.6e %.6e\n", step, err_l2, res)
        iter = step
        if res <= config.epsilon
            converged = true
            break
        end
    end
    runtime = time() - t_start

    output_dir = string(output_dir)
    isdir(output_dir) || mkpath(output_dir)
    tag = "nx$(config.nx)_ny$(config.ny)_nz$(config.nz)_steps$(iter)"
    history_path = joinpath(output_dir, "history_ssor_$(tag).txt")
    open(history_path, "w") do io
        write(io, String(take!(history)))
    end

    result = Solution(sol.x, sol.y, sol.z, sol.u, zero(T), iter)
    return converged, result, runtime
end

function sor_sweep_color!(u::Array{T,3}, f::Array{T,3}, config::SolverConfig,
                          inv_dx2::T, inv_dy2::T, inv_dz2::T, diag::T,
                          omega::T, color::Int) where {T<:Real}
    @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in (2 + ((color - ((j + k) & 1)) & 1)):2:config.nx+1
        rhs = f[i, j, k]
        sum_nb = (u[i+1, j, k] + u[i-1, j, k]) * inv_dx2 +
                 (u[i, j+1, k] + u[i, j-1, k]) * inv_dy2 +
                 (u[i, j, k+1] + u[i, j, k-1]) * inv_dz2
        u_star = (sum_nb - rhs) / diag
        u[i, j, k] = (one(T) - omega) * u[i, j, k] + omega * u_star
    end
end

function sor_sweep_color_backward!(u::Array{T,3}, f::Array{T,3}, config::SolverConfig,
                                   inv_dx2::T, inv_dy2::T, inv_dz2::T, diag::T,
                                   omega::T, color::Int) where {T<:Real}
    @inbounds for k in (config.nz+1):-1:2, j in (config.ny+1):-1:2,
        i in (2 + ((color - ((j + k) & 1)) & 1) + 2 * ((config.nx + 1 - (2 + ((color - ((j + k) & 1)) & 1))) ÷ 2)):-2:2
        rhs = f[i, j, k]
        sum_nb = (u[i+1, j, k] + u[i-1, j, k]) * inv_dx2 +
                 (u[i, j+1, k] + u[i, j-1, k]) * inv_dy2 +
                 (u[i, j, k+1] + u[i, j, k-1]) * inv_dz2
        u_star = (sum_nb - rhs) / diag
        u[i, j, k] = (one(T) - omega) * u[i, j, k] + omega * u_star
    end
end
