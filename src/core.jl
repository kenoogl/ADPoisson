# core.jl

"""
    grid_spacing(config; Lx=1, Ly=1, Lz=1)

Compute uniform grid spacing for the domain.
"""
function grid_spacing(config::SolverConfig; Lx=1, Ly=1, Lz=1)
    dx = Lx / config.nx
    dy = Ly / config.ny
    dz = Lz / config.nz
    return dx, dy, dz
end

"""
    make_grid(config; Lx=1, Ly=1, Lz=1)

Generate cell-centered coordinates for interior points only.
"""
function make_grid(config::SolverConfig; Lx=1, Ly=1, Lz=1)
    dx, dy, dz = grid_spacing(config; Lx=Lx, Ly=Ly, Lz=Lz)
    T = promote_type(typeof(dx), typeof(dy), typeof(dz))
    x = Vector{T}(undef, config.nx)
    y = Vector{T}(undef, config.ny)
    z = Vector{T}(undef, config.nz)
    for i in 1:config.nx
        x[i] = (i - 0.5) * dx
    end
    for j in 1:config.ny
        y[j] = (j - 0.5) * dy
    end
    for k in 1:config.nz
        z[k] = (k - 0.5) * dz
    end
    return x, y, z
end

"""
    initialize_solution(config, prob)

Create grid and zero initial condition with ghost cells.
"""
function initialize_solution(config::SolverConfig, prob::ProblemSpec)
    x, y, z = make_grid(config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    T = promote_type(typeof(prob.Lx), typeof(prob.alpha), typeof(config.dt))
    u = zeros(T, config.nx + 2, config.ny + 2, config.nz + 2)
    return Solution(x, y, z, u, zero(T), 0)
end

"""
    laplacian!(Lu, u, config; Lx=1, Ly=1, Lz=1)

Compute 7-point Laplacian for interior points. Ghost cells are assumed updated.
"""
function laplacian!(Lu::Array{T,3}, u::Array{T,3}, config::SolverConfig; Lx=1, Ly=1, Lz=1) where {T}
    dx, dy, dz = grid_spacing(config; Lx=Lx, Ly=Ly, Lz=Lz)
    inv_dx2 = one(T) / (dx * dx)
    inv_dy2 = one(T) / (dy * dy)
    inv_dz2 = one(T) / (dz * dz)
    @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
        Lu[i, j, k] = (u[i+1, j, k] - 2 * u[i, j, k] + u[i-1, j, k]) * inv_dx2 +
                     (u[i, j+1, k] - 2 * u[i, j, k] + u[i, j-1, k]) * inv_dy2 +
                     (u[i, j, k+1] - 2 * u[i, j, k] + u[i, j, k-1]) * inv_dz2
    end
    return Lu
end

"""
    taylor_step!(next, curr, f, m, config; Lx=1, Ly=1, Lz=1)

Compute Taylor coefficient u_{m+1} from u_m.
If f is time-invariant: for m >= 1, the source term is zero.
"""
function taylor_step!(next::Array{T,3}, curr::Array{T,3}, f::Array{T,3}, m::Int,
                      config::SolverConfig; Lx=1, Ly=1, Lz=1) where {T}
    laplacian!(next, curr, config; Lx=Lx, Ly=Ly, Lz=Lz)
    if m == 0
        @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
            next[i, j, k] = (next[i, j, k] - f[i, j, k]) / (m + 1)
        end
    else
        @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
            next[i, j, k] = next[i, j, k] / (m + 1)
        end
    end
    return next
end

"""
    accumulate_taylor!(acc, coeff, dt_pow)

Accumulate Taylor series sum: acc += coeff * dt_pow (interior only).
"""
function accumulate_taylor!(acc::Array{T,3}, coeff::Array{T,3}, dt_pow::T,
                            config::SolverConfig) where {T}
    @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
        acc[i, j, k] += coeff[i, j, k] * dt_pow
    end
    return acc
end

"""
    horner_update!(u_new, coeffs, dt, M)

Horner evaluation for full coefficient storage (verification use).
"""
function horner_update!(u_new::Array{T,3}, coeffs::TaylorArrays3D{T}, dt::T, M::Int,
                        config::SolverConfig) where {T}
    @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
        u_new[i, j, k] = coeffs.U[i, j, k, M + 1]
    end
    for m in M-1:-1:0
        @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
            u_new[i, j, k] = u_new[i, j, k] * dt + coeffs.U[i, j, k, m + 1]
        end
    end
    return u_new
end

"""
    taylor_series_update!(u_next, buffers, f, bc, config, prob)

Compute Taylor series update using ping-pong buffers.
"""
function taylor_series_update!(u_next::Array{T,3}, buffers::TaylorBuffers3D{T},
                               f::Array{T,3}, bc::BoundaryConditions,
                               config::SolverConfig, prob::ProblemSpec;
                               bc_order=:spec) where {T}
    bufA = buffers.bufA
    bufB = buffers.bufB
    acc = buffers.acc

    bufA .= u_next
    acc .= u_next
    dt_pow = one(T)

    for m in 0:config.M-1
        taylor_step!(bufB, bufA, f, m, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
        apply_bc!(bufB, bc, m + 1, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        dt_pow *= config.dt
        accumulate_taylor!(acc, bufB, dt_pow, config)
        bufA, bufB = bufB, bufA
    end

    u_next .= acc
    return u_next
end

"""
    compute_source!(f, prob, config)

Fill source term on interior points.
"""
function compute_source!(f::Array{T,3}, prob::ProblemSpec, config::SolverConfig) where {T}
    dx, dy, dz = grid_spacing(config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
        z = (k - 1.5) * dz
        y = (j - 1.5) * dy
        x = (i - 1.5) * dx
        f[i, j, k] = prob.source(x, y, z)
    end
    return f
end

"""
    compute_residual!(r, u, f, config; Lx=1, Ly=1, Lz=1)

Compute residual r = Lu - f on interior points.
"""
function compute_residual!(r::Array{T,3}, u::Array{T,3}, f::Array{T,3},
                           config::SolverConfig; Lx=1, Ly=1, Lz=1) where {T}
    laplacian!(r, u, config; Lx=Lx, Ly=Ly, Lz=Lz)
    @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
        r[i, j, k] -= f[i, j, k]
    end
    return r
end

"""
    compute_residual_norm!(r, u, f, config; Lx=1, Ly=1, Lz=1)

Compute residual r = Lu - f and return its L2 norm (interior only) in one pass.
"""
function compute_residual_norm!(r::Array{T,3}, u::Array{T,3}, f::Array{T,3},
                                config::SolverConfig; Lx=1, Ly=1, Lz=1) where {T}
    dx, dy, dz = grid_spacing(config; Lx=Lx, Ly=Ly, Lz=Lz)
    inv_dx2 = one(T) / (dx * dx)
    inv_dy2 = one(T) / (dy * dy)
    inv_dz2 = one(T) / (dz * dz)
    s = zero(T)
    @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
        Lu = (u[i+1, j, k] - 2 * u[i, j, k] + u[i-1, j, k]) * inv_dx2 +
             (u[i, j+1, k] - 2 * u[i, j, k] + u[i, j-1, k]) * inv_dy2 +
             (u[i, j, k+1] - 2 * u[i, j, k] + u[i, j, k-1]) * inv_dz2
        val = Lu - f[i, j, k]
        r[i, j, k] = val
        s += val * val
    end
    return sqrt(s)
end

"""
    l2_norm_interior(a, config)

Compute L2 norm over interior points.
"""
function l2_norm_interior(a::Array{T,3}, config::SolverConfig) where {T}
    s = zero(T)
    @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
        s += a[i, j, k]^2
    end
    return sqrt(s)
end

"""
    l2_error_exact(sol, prob, config)

Compute absolute L2 norm of error against the analytical solution.
"""
function l2_error_exact(sol::Solution, prob::ProblemSpec, config::SolverConfig)
    dx = prob.Lx / config.nx
    dy = prob.Ly / config.ny
    dz = prob.Lz / config.nz
    s = 0.0
    @inbounds for k in 1:config.nz, j in 1:config.ny, i in 1:config.nx
        uex = exact_solution(sol.x[i], sol.y[j], sol.z[k], prob.alpha)
        u = sol.u[i + 1, j + 1, k + 1]
        s += (u - uex)^2
    end
    return sqrt(s * dx * dy * dz)
end

"""
    exact_solution_array(sol, prob, config)

Precompute exact solution on interior points.
"""
function exact_solution_array(sol::Solution, prob::ProblemSpec, config::SolverConfig)
    uex = Array{eltype(sol.u)}(undef, config.nx, config.ny, config.nz)
    @inbounds for k in 1:config.nz, j in 1:config.ny, i in 1:config.nx
        uex[i, j, k] = exact_solution(sol.x[i], sol.y[j], sol.z[k], prob.alpha)
    end
    return uex
end

"""
    l2_error_exact_precomputed(u, uex, prob, config)

Compute absolute L2 norm of error using precomputed exact solution.
"""
function l2_error_exact_precomputed(u::Array{T,3}, uex::Array{T,3},
                                    prob::ProblemSpec, config::SolverConfig) where {T}
    dx = prob.Lx / config.nx
    dy = prob.Ly / config.ny
    dz = prob.Lz / config.nz
    s = zero(T)
    @inbounds for k in 1:config.nz, j in 1:config.ny, i in 1:config.nx
        diff = u[i + 1, j + 1, k + 1] - uex[i, j, k]
        s += diff^2
    end
    return sqrt(s * dx * dy * dz)
end

"""
    error_stats_precomputed(u, uex, prob, config)

Compute absolute L2 norm and max error using precomputed exact solution.
"""
function error_stats_precomputed(u::Array{T,3}, uex::Array{T,3},
                                 prob::ProblemSpec, config::SolverConfig) where {T}
    dx = prob.Lx / config.nx
    dy = prob.Ly / config.ny
    dz = prob.Lz / config.nz
    s = zero(T)
    err_max = zero(T)
    @inbounds for k in 1:config.nz, j in 1:config.ny, i in 1:config.nx
        diff = u[i + 1, j + 1, k + 1] - uex[i, j, k]
        s += diff^2
        ad = abs(diff)
        if ad > err_max
            err_max = ad
        end
    end
    return sqrt(s * dx * dy * dz), err_max
end

function solve_core(config::SolverConfig, prob::ProblemSpec; bc_order=:spec, output_dir="results")
    sol = initialize_solution(config, prob)
    bc = boundary_from_prob(prob)

    f = zeros(eltype(sol.u), config.nx + 2, config.ny + 2, config.nz + 2)
    compute_source!(f, prob, config)

    buffers = TaylorBuffers3D(similar(sol.u), similar(sol.u), similar(sol.u))
    r = similar(sol.u)
    u_exact = exact_solution_array(sol, prob, config)

    dx, dy, dz = grid_spacing(config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    Fo = config.dt * (1 / dx^2 + 1 / dy^2 + 1 / dz^2)
    if Fo > 0.5
        @warn "diffusion number exceeds 0.5" Fo=Fo dt=config.dt
    end

    iter = 0
    apply_bc!(sol.u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
    r0 = compute_residual_norm!(r, sol.u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    denom = max(r0, one(r0))
    history = IOBuffer()
    println(history, "# step err_l2 res_l2")

    rnorm = r0 / denom
    t_start = time()
    while true
        err_l2 = l2_error_exact_precomputed(sol.u, u_exact, prob, config)
        @printf(history, "%d %.6e %.6e\n", iter, err_l2, rnorm)
        if rnorm <= config.epsilon || iter >= config.max_steps
            break
        end
        taylor_series_update!(sol.u, buffers, f, bc, config, prob; bc_order=bc_order)
        iter += 1
        apply_bc!(sol.u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        res_l2 = compute_residual_norm!(r, sol.u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
        rnorm = res_l2 / denom
    end
    runtime = time() - t_start

    t = config.dt * iter
    result = Solution(sol.x, sol.y, sol.z, sol.u, t, iter)
    err_l2, err_max = error_stats_precomputed(result.u, u_exact, prob, config)
    output_dir = string(output_dir)
    isdir(output_dir) || mkpath(output_dir)
    tag = "nx$(config.nx)_ny$(config.ny)_nz$(config.nz)_M$(config.M)_steps$(iter)"
    history_path = joinpath(output_dir, "history_$(tag).txt")
    open(history_path, "w") do io
        write(io, String(take!(history)))
    end
    @info "summary" Fo=Fo dt=config.dt err_l2=@sprintf("%.3e", err_l2) err_max=@sprintf("%.3e", err_max) steps=iter runtime_s=@sprintf("%.3f", runtime)
    return result, runtime
end

"""
    solve_with_runtime(config, prob; bc_order=:spec, output_dir="results")

Main solver loop using Taylor series pseudo-time stepping.
Returns (Solution, runtime_s) where runtime is the solve loop only.
"""
function solve_with_runtime(config::SolverConfig, prob::ProblemSpec; bc_order=:spec, output_dir="results")
    return solve_core(config, prob; bc_order=bc_order, output_dir=output_dir)
end

"""
    solve(config, prob; bc_order=:spec, output_dir="results")

Main solver loop using Taylor series pseudo-time stepping.
"""
function solve(config::SolverConfig, prob::ProblemSpec; bc_order=:spec, output_dir="results")
    sol, _ = solve_core(config, prob; bc_order=bc_order, output_dir=output_dir)
    return sol
end
