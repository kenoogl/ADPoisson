# sor.jl

function _sor_diag_terms(config::SolverConfig{T}, prob::ProblemSpec{T}) where {T<:Real}
    dx, dy, dz = grid_spacing(config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    inv_dx2 = one(T) / (dx * dx)
    inv_dy2 = one(T) / (dy * dy)
    inv_dz2 = one(T) / (dz * dz)
    diag = 2 * (inv_dx2 + inv_dy2 + inv_dz2)
    return inv_dx2, inv_dy2, inv_dz2, diag
end

function sor_sweep_point_forward!(u::Array{T,3}, f::Array{T,3}, config::SolverConfig,
                                  inv_dx2::T, inv_dy2::T, inv_dz2::T, diag::T,
                                  omega::T) where {T<:Real}
    i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(u, config)
    @inbounds for k in k_lo:k_hi, j in j_lo:j_hi, i in i_lo:i_hi
        rhs = f[i, j, k]
        sum_nb = (u[i+1, j, k] + u[i-1, j, k]) * inv_dx2 +
                 (u[i, j+1, k] + u[i, j-1, k]) * inv_dy2 +
                 (u[i, j, k+1] + u[i, j, k-1]) * inv_dz2
        u_star = (sum_nb - rhs) / diag
        u[i, j, k] = (one(T) - omega) * u[i, j, k] + omega * u_star
    end
    return u
end

function sor_sweep_point_backward!(u::Array{T,3}, f::Array{T,3}, config::SolverConfig,
                                   inv_dx2::T, inv_dy2::T, inv_dz2::T, diag::T,
                                   omega::T) where {T<:Real}
    i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(u, config)
    @inbounds for k in k_hi:-1:k_lo, j in j_hi:-1:j_lo, i in i_hi:-1:i_lo
        rhs = f[i, j, k]
        sum_nb = (u[i+1, j, k] + u[i-1, j, k]) * inv_dx2 +
                 (u[i, j+1, k] + u[i, j-1, k]) * inv_dy2 +
                 (u[i, j, k+1] + u[i, j, k-1]) * inv_dz2
        u_star = (sum_nb - rhs) / diag
        u[i, j, k] = (one(T) - omega) * u[i, j, k] + omega * u_star
    end
    return u
end

function sor_sweep_color!(u::Array{T,3}, f::Array{T,3}, config::SolverConfig,
                          inv_dx2::T, inv_dy2::T, inv_dz2::T, diag::T,
                          omega::T, color::Int) where {T<:Real}
    i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(u, config)
    @inbounds for k in k_lo:k_hi, j in j_lo:j_hi
        i_start = i_lo + ((color - ((i_lo + j + k) & 1)) & 1)
        for i in i_start:2:i_hi
            rhs = f[i, j, k]
            sum_nb = (u[i+1, j, k] + u[i-1, j, k]) * inv_dx2 +
                     (u[i, j+1, k] + u[i, j-1, k]) * inv_dy2 +
                     (u[i, j, k+1] + u[i, j, k-1]) * inv_dz2
            u_star = (sum_nb - rhs) / diag
            u[i, j, k] = (one(T) - omega) * u[i, j, k] + omega * u_star
        end
    end
    return u
end

function sor_sweep_color_backward!(u::Array{T,3}, f::Array{T,3}, config::SolverConfig,
                                   inv_dx2::T, inv_dy2::T, inv_dz2::T, diag::T,
                                   omega::T, color::Int) where {T<:Real}
    i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(u, config)
    @inbounds for k in k_hi:-1:k_lo, j in j_hi:-1:j_lo
        i_start = i_hi - ((i_hi + j + k - color) & 1)
        for i in i_start:-2:i_lo
            rhs = f[i, j, k]
            sum_nb = (u[i+1, j, k] + u[i-1, j, k]) * inv_dx2 +
                     (u[i, j+1, k] + u[i, j-1, k]) * inv_dy2 +
                     (u[i, j, k+1] + u[i, j, k-1]) * inv_dz2
            u_star = (sum_nb - rhs) / diag
            u[i, j, k] = (one(T) - omega) * u[i, j, k] + omega * u_star
        end
    end
    return u
end

function _run_sor_iteration!(u::Array{T,3}, f::Array{T,3}, bc::BoundaryConditions,
                             config::SolverConfig, prob::ProblemSpec;
                             omega::T, bc_order::Symbol, mode::Symbol,
                             inv_dx2::T, inv_dy2::T, inv_dz2::T, diag::T) where {T<:Real}
    if mode === :point_sor
        apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        sor_sweep_point_forward!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega)
    elseif mode === :rb_sor
        apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        sor_sweep_color!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 0)
        apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        sor_sweep_color!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 1)
    elseif mode === :point_ssor
        apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        sor_sweep_point_forward!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega)
        apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        sor_sweep_point_backward!(u, f, config, inv_dx2, inv_dy2, inv_dz2, diag, omega)
    elseif mode === :rb_ssor
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
    else
        error("unknown SOR mode: $(mode)")
    end
    return u
end

function _sor_family_with_runtime!(sol::Solution{T}, f::Array{T,3}, bc::BoundaryConditions,
                                   prob::ProblemSpec, config::SolverConfig;
                                   omega::T=one(T), output_dir::AbstractString="results",
                                   bc_order::Symbol=:spec, mode::Symbol) where {T<:Real}
    u = sol.u
    r = similar(u)
    u_exact = exact_solution_array(sol, prob, config)
    inv_dx2, inv_dy2, inv_dz2, diag = _sor_diag_terms(config, prob)

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
        _run_sor_iteration!(u, f, bc, config, prob;
                            omega=omega, bc_order=bc_order, mode=mode,
                            inv_dx2=inv_dx2, inv_dy2=inv_dy2, inv_dz2=inv_dz2, diag=diag)
        res = compute_residual_norm!(r, u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz) / denom
        err_l2 = l2_error_exact_precomputed(u, u_exact, prob, config)
        @printf(history, "%d %.6e %.6e\n", step, err_l2, res)
        iter = step
        if !isfinite(res) || !isfinite(err_l2)
            break
        end
        if res <= config.epsilon
            converged = true
            break
        end
    end
    runtime = time() - t_start

    output_dir = string(output_dir)
    isdir(output_dir) || mkpath(output_dir)
    tag = "nx$(config.nx)_ny$(config.ny)_nz$(config.nz)_steps$(iter)"
    prefix = (mode === :point_sor || mode === :rb_sor) ? "history_sor_" : "history_ssor_"
    history_path = joinpath(output_dir, "$(prefix)$(tag).txt")
    open(history_path, "w") do io
        write(io, String(take!(history)))
    end

    result = Solution(sol.x, sol.y, sol.z, sol.u, zero(T), iter)
    return converged, result, runtime
end

function _solve_family_with_runtime(prob::ProblemSpec, config::SolverConfig;
                                    omega::Real=1.0, output_dir::AbstractString="results",
                                    bc_order=:spec, mode::Symbol)
    sol = initialize_solution(config, prob)
    bc = boundary_from_prob(prob)
    f = zeros(eltype(sol.u), size(sol.u))
    compute_source!(f, prob, config)
    omega_t = convert(eltype(sol.u), omega)
    _, sol_out, runtime = _sor_family_with_runtime!(sol, f, bc, prob, config;
                                                    omega=omega_t, output_dir=output_dir,
                                                    bc_order=bc_order, mode=mode)
    return sol_out, runtime
end

function _solve_family(prob::ProblemSpec, config::SolverConfig;
                       omega::Real=1.0, output_dir::AbstractString="results",
                       bc_order=:spec, mode::Symbol)
    sol_out, _ = _solve_family_with_runtime(prob, config;
                                            omega=omega, output_dir=output_dir,
                                            bc_order=bc_order, mode=mode)
    return sol_out
end

function _solve_family!(sol::Solution{T}, f::Array{T,3}, bc::BoundaryConditions,
                        prob::ProblemSpec, config::SolverConfig;
                        omega::T=one(T), output_dir::AbstractString="results",
                        bc_order=:spec, mode::Symbol) where {T<:Real}
    converged, result, _ = _sor_family_with_runtime!(sol, f, bc, prob, config;
                                                     omega=omega, output_dir=output_dir,
                                                     bc_order=bc_order, mode=mode)
    return converged, result
end

function sor_solve_with_runtime(prob::ProblemSpec, config::SolverConfig;
                                omega::Real=1.0, output_dir::AbstractString="results",
                                bc_order=:spec)
    return _solve_family_with_runtime(prob, config;
                                      omega=omega, output_dir=output_dir,
                                      bc_order=bc_order, mode=:point_sor)
end

function rbsor_solve_with_runtime(prob::ProblemSpec, config::SolverConfig;
                                  omega::Real=1.0, output_dir::AbstractString="results",
                                  bc_order=:spec)
    return _solve_family_with_runtime(prob, config;
                                      omega=omega, output_dir=output_dir,
                                      bc_order=bc_order, mode=:rb_sor)
end

function ssor_solve_with_runtime(prob::ProblemSpec, config::SolverConfig;
                                 omega::Real=1.0, output_dir::AbstractString="results",
                                 bc_order=:spec)
    return _solve_family_with_runtime(prob, config;
                                      omega=omega, output_dir=output_dir,
                                      bc_order=bc_order, mode=:point_ssor)
end

function rbssor_solve_with_runtime(prob::ProblemSpec, config::SolverConfig;
                                   omega::Real=1.0, output_dir::AbstractString="results",
                                   bc_order=:spec)
    return _solve_family_with_runtime(prob, config;
                                      omega=omega, output_dir=output_dir,
                                      bc_order=bc_order, mode=:rb_ssor)
end

function sor_solve(prob::ProblemSpec, config::SolverConfig;
                   omega::Real=1.0, output_dir::AbstractString="results",
                   bc_order=:spec)
    return _solve_family(prob, config;
                         omega=omega, output_dir=output_dir,
                         bc_order=bc_order, mode=:point_sor)
end

function rbsor_solve(prob::ProblemSpec, config::SolverConfig;
                     omega::Real=1.0, output_dir::AbstractString="results",
                     bc_order=:spec)
    return _solve_family(prob, config;
                         omega=omega, output_dir=output_dir,
                         bc_order=bc_order, mode=:rb_sor)
end

function ssor_solve(prob::ProblemSpec, config::SolverConfig;
                    omega::Real=1.0, output_dir::AbstractString="results",
                    bc_order=:spec)
    return _solve_family(prob, config;
                         omega=omega, output_dir=output_dir,
                         bc_order=bc_order, mode=:point_ssor)
end

function rbssor_solve(prob::ProblemSpec, config::SolverConfig;
                      omega::Real=1.0, output_dir::AbstractString="results",
                      bc_order=:spec)
    return _solve_family(prob, config;
                         omega=omega, output_dir=output_dir,
                         bc_order=bc_order, mode=:rb_ssor)
end

function sor_solve!(sol::Solution{T}, f::Array{T,3}, bc::BoundaryConditions,
                    prob::ProblemSpec, config::SolverConfig;
                    omega::T=one(T), output_dir::AbstractString="results",
                    bc_order=:spec) where {T<:Real}
    return _solve_family!(sol, f, bc, prob, config;
                          omega=omega, output_dir=output_dir,
                          bc_order=bc_order, mode=:point_sor)
end

function rbsor_solve!(sol::Solution{T}, f::Array{T,3}, bc::BoundaryConditions,
                      prob::ProblemSpec, config::SolverConfig;
                      omega::T=one(T), output_dir::AbstractString="results",
                      bc_order=:spec) where {T<:Real}
    return _solve_family!(sol, f, bc, prob, config;
                          omega=omega, output_dir=output_dir,
                          bc_order=bc_order, mode=:rb_sor)
end

function ssor_solve!(sol::Solution{T}, f::Array{T,3}, bc::BoundaryConditions,
                     prob::ProblemSpec, config::SolverConfig;
                     omega::T=one(T), output_dir::AbstractString="results",
                     bc_order=:spec) where {T<:Real}
    return _solve_family!(sol, f, bc, prob, config;
                          omega=omega, output_dir=output_dir,
                          bc_order=bc_order, mode=:point_ssor)
end

function rbssor_solve!(sol::Solution{T}, f::Array{T,3}, bc::BoundaryConditions,
                       prob::ProblemSpec, config::SolverConfig;
                       omega::T=one(T), output_dir::AbstractString="results",
                       bc_order=:spec) where {T<:Real}
    return _solve_family!(sol, f, bc, prob, config;
                          omega=omega, output_dir=output_dir,
                          bc_order=bc_order, mode=:rb_ssor)
end
