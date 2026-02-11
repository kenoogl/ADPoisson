# cg.jl

"""
    cg_solve_with_runtime(prob, config; precond=:none, omega_ssor=1.0, output_dir="results", bc_order=:spec)

Solve Poisson equation using PCG and return (Solution, runtime_s).
"""
function cg_solve_with_runtime(prob::ProblemSpec, config::SolverConfig;
                               precond::Symbol=:none, omega_ssor::Real=1.0,
                               output_dir::AbstractString="results", bc_order=:spec)
    sol = initialize_solution(config, prob)
    bc = boundary_from_prob(prob)
    f = zeros(eltype(sol.u), size(sol.u))
    compute_source!(f, prob, config)
    omega_t = convert(eltype(sol.u), omega_ssor)
    _, sol_out, runtime = cg_solve_with_runtime!(sol, f, bc, prob, config;
                                                 precond=precond, omega_ssor=omega_t,
                                                 output_dir=output_dir, bc_order=bc_order)
    return sol_out, runtime
end

"""
    cg_solve(prob, config; precond=:none, omega_ssor=1.0, output_dir="results", bc_order=:spec)

Solve Poisson equation using PCG and return Solution.
"""
function cg_solve(prob::ProblemSpec, config::SolverConfig;
                  precond::Symbol=:none, omega_ssor::Real=1.0,
                  output_dir::AbstractString="results", bc_order=:spec)
    sol_out, _ = cg_solve_with_runtime(prob, config;
                                       precond=precond, omega_ssor=omega_ssor,
                                       output_dir=output_dir, bc_order=bc_order)
    return sol_out
end

"""
    cg_solve!(sol, f, bc, prob, config; precond=:none, omega_ssor=1.0, output_dir="results", bc_order=:spec)

PCG solver with optional preconditioner.
"""
function cg_solve!(sol::Solution{T}, f::Array{T,3}, bc::BoundaryConditions,
                   prob::ProblemSpec, config::SolverConfig;
                   precond::Symbol=:none, omega_ssor::T=one(T),
                   output_dir::AbstractString="results",
                   bc_order=:spec) where {T<:Real}
    converged, result, _ = cg_solve_with_runtime!(sol, f, bc, prob, config;
                                                  precond=precond, omega_ssor=omega_ssor,
                                                  output_dir=output_dir, bc_order=bc_order)
    return converged, result
end

function cg_solve_with_runtime!(sol::Solution{T}, f::Array{T,3}, bc::BoundaryConditions,
                                prob::ProblemSpec, config::SolverConfig;
                                precond::Symbol=:none, omega_ssor::T=one(T),
                                output_dir::AbstractString="results",
                                bc_order=:spec) where {T<:Real}
    precond === :ssor || precond === :none || error("precond must be :ssor or :none")
    u = sol.u
    r = similar(u)
    p = similar(u)
    z = similar(u)
    q = similar(u)
    u_exact = exact_solution_array(sol, prob, config)

    bc0 = zero_boundary_conditions(T)

    apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
    r0 = compute_residual_norm!(r, u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    denom = max(r0, one(r0))

    history = IOBuffer()
    println(history, "# step err_l2 res_l2")
    err_l2 = l2_error_exact_precomputed(u, u_exact, prob, config)
    @printf(history, "%d %.6e %.6e\n", 0, err_l2, r0 / denom)

    if precond === :ssor
        ssor_precond!(z, r, bc0, config, prob; omega=omega_ssor)
        copy_interior!(p, z, config)
        rho = dot_interior(r, z, config)
    else
        copy_interior!(z, r, config)
        copy_interior!(p, r, config)
        rho = dot_interior(r, r, config)
    end

    converged = false
    iter = 0
    t_start = time()
    for step in 1:config.max_steps
        apply_A!(q, p, bc0, config, prob)
        denom_cg = dot_interior(p, q, config)
        if abs(denom_cg) < eps(denom_cg)
            break
        end
        alpha = rho / denom_cg
        axpy_interior!(u, p, alpha, config)
        axpy_interior!(r, q, -alpha, config)

        res = l2_norm_interior(r, config) / denom
        err_l2 = l2_error_exact_precomputed(u, u_exact, prob, config)
        @printf(history, "%d %.6e %.6e\n", step, err_l2, res)
        iter = step
        if res <= config.epsilon
            converged = true
            break
        end

        if precond === :ssor
            ssor_precond!(z, r, bc0, config, prob; omega=omega_ssor)
            rho_new = dot_interior(r, z, config)
        else
            copy_interior!(z, r, config)
            rho_new = dot_interior(r, r, config)
        end
        if abs(rho_new) < eps(rho_new)
            break
        end
        beta = rho_new / rho
        update_p!(p, z, beta, config)
        rho = rho_new
    end
    runtime = time() - t_start

    output_dir = string(output_dir)
    isdir(output_dir) || mkpath(output_dir)
    tag = "nx$(config.nx)_ny$(config.ny)_nz$(config.nz)_steps$(iter)"
    history_path = joinpath(output_dir, "history_cg_$(tag).txt")
    open(history_path, "w") do io
        write(io, String(take!(history)))
    end

    result = Solution(sol.x, sol.y, sol.z, sol.u, zero(T), iter)
    return converged, result, runtime
end

function zero_boundary_conditions(::Type{T}) where {T<:Real}
    g0 = (a, b) -> zero(T)
    return BoundaryConditions(g0, g0, g0, g0, g0, g0)
end

function apply_A!(q::Array{T,3}, p::Array{T,3}, bc0::BoundaryConditions,
                  config::SolverConfig, prob::ProblemSpec) where {T<:Real}
    apply_bc!(p, bc0, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=:spec)
    laplacian!(q, p, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(q, config)
    @inbounds for k in k_lo:k_hi, j in j_lo:j_hi, i in i_lo:i_hi
        q[i, j, k] = -q[i, j, k]
    end
end

function ssor_precond!(z::Array{T,3}, r::Array{T,3}, bc0::BoundaryConditions,
                       config::SolverConfig, prob::ProblemSpec;
                       omega::T=one(T)) where {T<:Real}
    z .= zero(T)
    dx, dy, dz = grid_spacing(config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    inv_dx2 = one(T) / (dx * dx)
    inv_dy2 = one(T) / (dy * dy)
    inv_dz2 = one(T) / (dz * dz)
    diag = 2 * (inv_dx2 + inv_dy2 + inv_dz2)

    # RBSSOR: forward R→B, backward B→R, forward B→R, backward R→B
    # NOTE: order=:spec は ghost が隣接内点のみ依存のため、色ごとの再適用は不要。
    #       order=:high を使う場合は ghost が複数内点に依存するため、
    #       色ごとに境界更新するか、事前に :spec に固定すること。
    apply_bc!(z, bc0, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=:spec)
    sor_sweep_color_rhs_forward!(z, r, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 0)
    sor_sweep_color_rhs_forward!(z, r, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 1)

    apply_bc!(z, bc0, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=:spec)
    sor_sweep_color_rhs_backward!(z, r, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 1)
    sor_sweep_color_rhs_backward!(z, r, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 0)

    apply_bc!(z, bc0, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=:spec)
    sor_sweep_color_rhs_forward!(z, r, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 1)
    sor_sweep_color_rhs_forward!(z, r, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 0)

    apply_bc!(z, bc0, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=:spec)
    sor_sweep_color_rhs_backward!(z, r, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 0)
    sor_sweep_color_rhs_backward!(z, r, config, inv_dx2, inv_dy2, inv_dz2, diag, omega, 1)
    return z
end

function sor_sweep_color_rhs_forward!(u::Array{T,3}, r::Array{T,3}, config::SolverConfig,
                                      inv_dx2::T, inv_dy2::T, inv_dz2::T, diag::T,
                                      omega::T, color::Int) where {T<:Real}
    i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(u, config)
    @inbounds for k in k_lo:k_hi, j in j_lo:j_hi
        i_start = i_lo + ((color - ((i_lo + j + k) & 1)) & 1)
        for i in i_start:2:i_hi
        rhs = -r[i, j, k]
        sum_nb = (u[i+1, j, k] + u[i-1, j, k]) * inv_dx2 +
                 (u[i, j+1, k] + u[i, j-1, k]) * inv_dy2 +
                 (u[i, j, k+1] + u[i, j, k-1]) * inv_dz2
        u_star = (sum_nb - rhs) / diag
        u[i, j, k] = (one(T) - omega) * u[i, j, k] + omega * u_star
        end
    end
end

function sor_sweep_color_rhs_backward!(u::Array{T,3}, r::Array{T,3}, config::SolverConfig,
                                       inv_dx2::T, inv_dy2::T, inv_dz2::T, diag::T,
                                       omega::T, color::Int) where {T<:Real}
    i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(u, config)
    @inbounds for k in k_hi:-1:k_lo, j in j_hi:-1:j_lo
        i_start = i_hi - ((i_hi + j + k - color) & 1)
        for i in i_start:-2:i_lo
        rhs = -r[i, j, k]
        sum_nb = (u[i+1, j, k] + u[i-1, j, k]) * inv_dx2 +
                 (u[i, j+1, k] + u[i, j-1, k]) * inv_dy2 +
                 (u[i, j, k+1] + u[i, j, k-1]) * inv_dz2
        u_star = (sum_nb - rhs) / diag
        u[i, j, k] = (one(T) - omega) * u[i, j, k] + omega * u_star
        end
    end
end

function dot_interior(a::Array{T,3}, b::Array{T,3}, config::SolverConfig) where {T<:Real}
    s = zero(T)
    i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(a, config)
    @inbounds for k in k_lo:k_hi, j in j_lo:j_hi, i in i_lo:i_hi
        s += a[i, j, k] * b[i, j, k]
    end
    return s
end

function copy_interior!(dest::Array{T,3}, src::Array{T,3}, config::SolverConfig) where {T<:Real}
    i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(dest, config)
    @inbounds for k in k_lo:k_hi, j in j_lo:j_hi, i in i_lo:i_hi
        dest[i, j, k] = src[i, j, k]
    end
    return dest
end

function update_p!(p::Array{T,3}, z::Array{T,3}, beta::T, config::SolverConfig) where {T<:Real}
    i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(p, config)
    @inbounds for k in k_lo:k_hi, j in j_lo:j_hi, i in i_lo:i_hi
        p[i, j, k] = z[i, j, k] + beta * p[i, j, k]
    end
    return p
end

function axpy_interior!(x::Array{T,3}, y::Array{T,3}, a::T, config::SolverConfig) where {T<:Real}
    i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(x, config)
    @inbounds for k in k_lo:k_hi, j in j_lo:j_hi, i in i_lo:i_hi
        x[i, j, k] += a * y[i, j, k]
    end
    return x
end
