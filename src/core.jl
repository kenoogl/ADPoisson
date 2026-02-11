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

@inline function ghost_layers(a::Array{T,3}, config::SolverConfig) where {T}
    gx = (size(a, 1) - config.nx) ÷ 2
    gy = (size(a, 2) - config.ny) ÷ 2
    gz = (size(a, 3) - config.nz) ÷ 2
    if gx < 1 || gy < 1 || gz < 1 ||
       size(a, 1) != config.nx + 2 * gx ||
       size(a, 2) != config.ny + 2 * gy ||
       size(a, 3) != config.nz + 2 * gz
        error("array size does not match config with ghost layers")
    end
    return gx, gy, gz
end

@inline function interior_bounds(a::Array{T,3}, config::SolverConfig) where {T}
    gx, gy, gz = ghost_layers(a, config)
    return gx + 1, gx + config.nx, gy + 1, gy + config.ny, gz + 1, gz + config.nz
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
    ghost = 2
    u = zeros(T, config.nx + 2 * ghost, config.ny + 2 * ghost, config.nz + 2 * ghost)
    return Solution(x, y, z, u, zero(T), 0)
end

"""
    laplacian!(Lu, u, config; Lx=1, Ly=1, Lz=1, order=:second)

Compute Laplacian for interior points. Ghost cells are assumed updated.
"""
function laplacian!(Lu::Array{T,3}, u::Array{T,3}, config::SolverConfig;
                    Lx=1, Ly=1, Lz=1, order::Symbol=:second) where {T}
    if order === :second
        return laplacian2!(Lu, u, config; Lx=Lx, Ly=Ly, Lz=Lz)
    elseif order === :fourth
        return laplacian4!(Lu, u, config; Lx=Lx, Ly=Ly, Lz=Lz)
    end
    throw(ArgumentError("unknown laplacian order: $(order). use :second or :fourth"))
end

"""
    laplacian2!(Lu, u, config; Lx=1, Ly=1, Lz=1)

Compute 2nd-order 7-point Laplacian for interior points.
"""
function laplacian2!(Lu::Array{T,3}, u::Array{T,3}, config::SolverConfig; Lx=1, Ly=1, Lz=1) where {T}
    dx, dy, dz = grid_spacing(config; Lx=Lx, Ly=Ly, Lz=Lz)
    inv_dx2 = one(T) / (dx * dx)
    inv_dy2 = one(T) / (dy * dy)
    inv_dz2 = one(T) / (dz * dz)
    ngx = (size(u, 1) - config.nx) ÷ 2
    ngy = (size(u, 2) - config.ny) ÷ 2
    ngz = (size(u, 3) - config.nz) ÷ 2
    @inbounds for k in ngz+1:ngz+config.nz, j in ngy+1:ngy+config.ny, i in ngx+1:ngx+config.nx
        Lu[i, j, k] = (u[i+1, j, k] - 2 * u[i, j, k] + u[i-1, j, k]) * inv_dx2 +
                     (u[i, j+1, k] - 2 * u[i, j, k] + u[i, j-1, k]) * inv_dy2 +
                     (u[i, j, k+1] - 2 * u[i, j, k] + u[i, j, k-1]) * inv_dz2
    end
    return Lu
end

"""
    laplacian4!(Lu, u, config; Lx=1, Ly=1, Lz=1)

Compute 4th-order Laplacian (radius=2, cross terms omitted) for interior points.
Requires at least two ghost layers in every direction.
"""
function laplacian4!(Lu::Array{T,3}, u::Array{T,3}, config::SolverConfig; Lx=1, Ly=1, Lz=1) where {T}
    dx, dy, dz = grid_spacing(config; Lx=Lx, Ly=Ly, Lz=Lz)
    inv_12dx2 = one(T) / (12 * dx * dx)
    inv_12dy2 = one(T) / (12 * dy * dy)
    inv_12dz2 = one(T) / (12 * dz * dz)
    ngx = (size(u, 1) - config.nx) ÷ 2
    ngy = (size(u, 2) - config.ny) ÷ 2
    ngz = (size(u, 3) - config.nz) ÷ 2
    if ngx < 2 || ngy < 2 || ngz < 2
        throw(ArgumentError("laplacian4! requires two ghost layers; got (ngx,ngy,ngz)=($(ngx),$(ngy),$(ngz))"))
    end

    @inbounds for k in ngz+1:ngz+config.nz, j in ngy+1:ngy+config.ny, i in ngx+1:ngx+config.nx
        d2x = (-u[i+2, j, k] + 16 * u[i+1, j, k] - 30 * u[i, j, k] +
               16 * u[i-1, j, k] - u[i-2, j, k]) * inv_12dx2
        d2y = (-u[i, j+2, k] + 16 * u[i, j+1, k] - 30 * u[i, j, k] +
               16 * u[i, j-1, k] - u[i, j-2, k]) * inv_12dy2
        d2z = (-u[i, j, k+2] + 16 * u[i, j, k+1] - 30 * u[i, j, k] +
               16 * u[i, j, k-1] - u[i, j, k-2]) * inv_12dz2
        Lu[i, j, k] = d2x + d2y + d2z
    end
    return Lu
end

"""
    taylor_step!(next, curr, f, m, config; Lx=1, Ly=1, Lz=1)

Compute Taylor coefficient u_{m+1} from u_m.
If f is time-invariant: for m >= 1, the source term is zero.
"""
function taylor_step!(next::Array{T,3}, curr::Array{T,3}, f::Array{T,3}, m::Int,
                      config::SolverConfig; Lx=1, Ly=1, Lz=1,
                      lap_order::Symbol=:second) where {T}
    laplacian!(next, curr, config; Lx=Lx, Ly=Ly, Lz=Lz, order=lap_order)
    if m == 0
        i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(next, config)
        @inbounds for k in k_lo:k_hi, j in j_lo:j_hi, i in i_lo:i_hi
            next[i, j, k] = (next[i, j, k] - f[i, j, k]) / (m + 1)
        end
    else
        i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(next, config)
        @inbounds for k in k_lo:k_hi, j in j_lo:j_hi, i in i_lo:i_hi
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
    i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(acc, config)
    @inbounds for k in k_lo:k_hi, j in j_lo:j_hi, i in i_lo:i_hi
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
    i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(u_new, config)
    @inbounds for k in k_lo:k_hi, j in j_lo:j_hi, i in i_lo:i_hi
        u_new[i, j, k] = coeffs.U[i, j, k, M + 1]
    end
    for m in M-1:-1:0
        @inbounds for k in k_lo:k_hi, j in j_lo:j_hi, i in i_lo:i_hi
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
                               bc_order=:spec, lap_order::Symbol=:second) where {T}
    bufA = buffers.bufA
    bufB = buffers.bufB
    acc = buffers.acc

    bufA .= u_next
    acc .= u_next
    dt_pow = one(T)

    for m in 0:config.M-1
        taylor_step!(bufB, bufA, f, m, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz,
                     lap_order=lap_order)
        apply_bc!(bufB, bc, m + 1, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        dt_pow *= config.dt
        accumulate_taylor!(acc, bufB, dt_pow, config)
        bufA, bufB = bufB, bufA
    end

    u_next .= acc
    apply_bc!(u_next, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
    return u_next
end

"""
    taylor_series_update_reuse!(u_next, buffers, r0, f, bc, config, prob; bc_order=:spec)

Compute Taylor series update using precomputed residual r0 = Lu - f as the m=0 coefficient.
"""
function taylor_series_update_reuse!(u_next::Array{T,3}, buffers::TaylorBuffers3D{T},
                                     r0::Array{T,3}, f::Array{T,3}, bc::BoundaryConditions,
                                     config::SolverConfig, prob::ProblemSpec;
                                     bc_order=:spec, lap_order::Symbol=:second) where {T}
    if config.M == 0
        apply_bc!(u_next, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        return u_next
    end
    bufA = buffers.bufA
    bufB = buffers.bufB
    acc = buffers.acc

    bufA .= u_next
    acc .= u_next
    dt_pow = one(T)

    bufB .= r0
    apply_bc!(bufB, bc, 1, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
    dt_pow *= config.dt
    accumulate_taylor!(acc, bufB, dt_pow, config)
    bufA, bufB = bufB, bufA

    for m in 1:config.M-1
        taylor_step!(bufB, bufA, f, m, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz,
                     lap_order=lap_order)
        apply_bc!(bufB, bc, m + 1, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        dt_pow *= config.dt
        accumulate_taylor!(acc, bufB, dt_pow, config)
        bufA, bufB = bufB, bufA # not copy, pointer exchange
    end

    u_next .= acc
    apply_bc!(u_next, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
    return u_next
end

"""
    compute_source!(f, prob, config)

Fill source term on interior points.
"""
function compute_source!(f::Array{T,3}, prob::ProblemSpec, config::SolverConfig) where {T}
    dx, dy, dz = grid_spacing(config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    gx, gy, gz = ghost_layers(f, config)
    coord_offset_x = gx + 0.5
    coord_offset_y = gy + 0.5
    coord_offset_z = gz + 0.5
    @inbounds for k in gz+1:gz+config.nz, j in gy+1:gy+config.ny, i in gx+1:gx+config.nx
        z = (k - coord_offset_z) * dz
        y = (j - coord_offset_y) * dy
        x = (i - coord_offset_x) * dx
        f[i, j, k] = prob.source(x, y, z)
    end
    return f
end

"""
    compute_residual!(r, u, f, config; Lx=1, Ly=1, Lz=1)

Compute residual r = Lu - f on interior points.
"""
function compute_residual!(r::Array{T,3}, u::Array{T,3}, f::Array{T,3},
                           config::SolverConfig; Lx=1, Ly=1, Lz=1,
                           lap_order::Symbol=:second) where {T}
    laplacian!(r, u, config; Lx=Lx, Ly=Ly, Lz=Lz, order=lap_order)
    ngx = (size(u, 1) - config.nx) ÷ 2
    ngy = (size(u, 2) - config.ny) ÷ 2
    ngz = (size(u, 3) - config.nz) ÷ 2
    @inbounds for k in ngz+1:ngz+config.nz, j in ngy+1:ngy+config.ny, i in ngx+1:ngx+config.nx
        r[i, j, k] -= f[i, j, k]
    end
    return r
end

"""
    compute_residual_norm!(r, u, f, config; Lx=1, Ly=1, Lz=1)

Compute residual r = Lu - f and return its L2 norm (interior only) in one pass.
"""
function compute_residual_norm!(r::Array{T,3}, u::Array{T,3}, f::Array{T,3},
                                config::SolverConfig; Lx=1, Ly=1, Lz=1,
                                lap_order::Symbol=:second) where {T}
    laplacian!(r, u, config; Lx=Lx, Ly=Ly, Lz=Lz, order=lap_order)
    ngx = (size(u, 1) - config.nx) ÷ 2
    ngy = (size(u, 2) - config.ny) ÷ 2
    ngz = (size(u, 3) - config.nz) ÷ 2
    s = zero(T)
    @inbounds for k in ngz+1:ngz+config.nz, j in ngy+1:ngy+config.ny, i in ngx+1:ngx+config.nx
        val = r[i, j, k] - f[i, j, k]
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
    i_lo, i_hi, j_lo, j_hi, k_lo, k_hi = interior_bounds(a, config)
    @inbounds for k in k_lo:k_hi, j in j_lo:j_hi, i in i_lo:i_hi
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
    gx = (size(sol.u, 1) - config.nx) ÷ 2
    gy = (size(sol.u, 2) - config.ny) ÷ 2
    gz = (size(sol.u, 3) - config.nz) ÷ 2
    @inbounds for k in 1:config.nz, j in 1:config.ny, i in 1:config.nx
        uex = exact_solution(sol.x[i], sol.y[j], sol.z[k], prob.alpha)
        u = sol.u[i + gx, j + gy, k + gz]
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
    gx = (size(u, 1) - config.nx) ÷ 2
    gy = (size(u, 2) - config.ny) ÷ 2
    gz = (size(u, 3) - config.nz) ÷ 2
    @inbounds for k in 1:config.nz, j in 1:config.ny, i in 1:config.nx
        diff = u[i + gx, j + gy, k + gz] - uex[i, j, k]
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
    gx = (size(u, 1) - config.nx) ÷ 2
    gy = (size(u, 2) - config.ny) ÷ 2
    gz = (size(u, 3) - config.nz) ÷ 2
    @inbounds for k in 1:config.nz, j in 1:config.ny, i in 1:config.nx
        diff = u[i + gx, j + gy, k + gz] - uex[i, j, k]
        s += diff^2
        ad = abs(diff)
        if ad > err_max
            err_max = ad
        end
    end
    return sqrt(s * dx * dy * dz), err_max
end

function solve_core(config::SolverConfig, prob::ProblemSpec;
                    bc_order=:spec, output_dir="results",
                    lap_order::Symbol=:second,
                    mg_interval::Int=0, mg_dt_scale::Real=2.0, mg_M::Int=4,
                    mg_vcycle::Bool=false, mg_nu1::Int=1, mg_nu2::Int=1,
                    mg_level_Ms::Union{Nothing,AbstractVector{Int}}=nothing,
                    mg_level_dt_scales::Union{Nothing,AbstractVector{<:Real}}=nothing,
                    mg_correction::Symbol=:classic,
                    mg_corr_M::Int=2, mg_corr_dt_scale::Real=1.0, mg_corr_steps::Int=1,
                    mg_corr_nu1::Union{Nothing,Int}=nothing,
                    mg_corr_nu2::Union{Nothing,Int}=nothing,
                    debug_residual::Bool=false, debug_vcycle::Bool=false,
                    mg_workspace=nothing)
    sol = initialize_solution(config, prob)
    bc = boundary_from_prob(prob)

    f = zeros(eltype(sol.u), size(sol.u))
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
    r0 = compute_residual_norm!(r, sol.u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz,
                                lap_order=lap_order)
    denom = max(r0, one(r0))
    history = IOBuffer()
    println(history, "# step err_l2 res_l2")
    debug_io = debug_residual ? IOBuffer() : nothing
    r_check = debug_residual ? similar(sol.u) : nothing
    if debug_residual
        println(debug_io, "# step res_l2 res_l2_check diff")
    end
    vcycle_io = debug_vcycle ? IOBuffer() : nothing
    if debug_vcycle
        println(vcycle_io, "# level nx ny res_l2_norm")
    end
    if mg_vcycle && mg_interval > 0
        if mg_workspace !== nothing &&
           size(mg_workspace.levels[1].r, 1) == size(sol.u, 1) &&
           size(mg_workspace.levels[1].r, 2) == size(sol.u, 2) &&
           size(mg_workspace.levels[1].r, 3) == size(sol.u, 3)
            mg_ws = mg_workspace
        else
            mg_ws = build_mg_workspace(sol.u, config)
        end
    else
        mg_ws = nothing
    end
    rnorm = r0 / denom
    t_start = time()
    while true
        err_l2 = l2_error_exact_precomputed(sol.u, u_exact, prob, config)
        @printf(history, "%d %.6e %.6e\n", iter, err_l2, rnorm)
        if debug_residual
            res_check = compute_residual_norm!(r_check, sol.u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz,
                                               lap_order=lap_order) / denom
            @printf(debug_io, "%d %.6e %.6e %.6e\n", iter, rnorm, res_check, abs(rnorm - res_check))
        end
        if rnorm <= config.epsilon || iter >= config.max_steps
            break
        end
        taylor_series_update_reuse!(sol.u, buffers, r, f, bc, config, prob;
                                    bc_order=bc_order, lap_order=lap_order)
        iter += 1
        if mg_vcycle && mg_interval > 0 && (iter % mg_interval == 0)
            vcycle!(sol.u, f, bc, config, prob, mg_ws;
                    nu1=mg_nu1, nu2=mg_nu2,
                    dt_scale=mg_dt_scale, M=mg_M,
                    level_dt_scales=mg_level_dt_scales, level_Ms=mg_level_Ms,
                    correction_mode=mg_correction,
                    corr_M=mg_corr_M, corr_dt_scale=mg_corr_dt_scale, corr_steps=mg_corr_steps,
                    corr_nu1=mg_corr_nu1, corr_nu2=mg_corr_nu2,
                    bc_order=bc_order, debug_io=vcycle_io, debug_denom=denom)
        end
        res_l2 = compute_residual_norm!(r, sol.u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz,
                                        lap_order=lap_order)
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
    if debug_residual
        debug_path = joinpath(output_dir, "residual_check_$(tag).txt")
        open(debug_path, "w") do io
            write(io, String(take!(debug_io)))
        end
    end
    if debug_vcycle
        debug_path = joinpath(output_dir, "vcycle_residual_$(tag).txt")
        open(debug_path, "w") do io
            write(io, String(take!(vcycle_io)))
        end
    end
    @info "summary" Fo=Fo dt=config.dt err_l2=@sprintf("%.3e", err_l2) err_max=@sprintf("%.3e", err_max) steps=iter runtime_s=@sprintf("%.3f", runtime)
    return result, runtime
end

"""
    solve_with_runtime(config, prob; bc_order=:spec, output_dir="results")

Main solver loop using Taylor series pseudo-time stepping.
Returns (Solution, runtime_s) where runtime is the solve loop only.
"""
function solve_with_runtime(config::SolverConfig, prob::ProblemSpec;
                            bc_order=:spec, output_dir="results",
                            lap_order::Symbol=:second,
                            mg_interval::Int=0, mg_dt_scale::Real=2.0, mg_M::Int=4,
                            mg_vcycle::Bool=false, mg_nu1::Int=1, mg_nu2::Int=1,
                            mg_level_Ms::Union{Nothing,AbstractVector{Int}}=nothing,
                            mg_level_dt_scales::Union{Nothing,AbstractVector{<:Real}}=nothing,
                            mg_correction::Symbol=:classic,
                            mg_corr_M::Int=2, mg_corr_dt_scale::Real=1.0, mg_corr_steps::Int=1,
                            mg_corr_nu1::Union{Nothing,Int}=nothing,
                            mg_corr_nu2::Union{Nothing,Int}=nothing,
                            debug_residual::Bool=false, debug_vcycle::Bool=false,
                            mg_workspace=nothing)
    return solve_core(config, prob; bc_order=bc_order, output_dir=output_dir,
                      lap_order=lap_order,
                      mg_interval=mg_interval, mg_dt_scale=mg_dt_scale, mg_M=mg_M,
                      mg_vcycle=mg_vcycle, mg_nu1=mg_nu1, mg_nu2=mg_nu2,
                      mg_level_Ms=mg_level_Ms, mg_level_dt_scales=mg_level_dt_scales,
                      mg_correction=mg_correction,
                      mg_corr_M=mg_corr_M, mg_corr_dt_scale=mg_corr_dt_scale, mg_corr_steps=mg_corr_steps,
                      mg_corr_nu1=mg_corr_nu1, mg_corr_nu2=mg_corr_nu2,
                      debug_residual=debug_residual, debug_vcycle=debug_vcycle,
                      mg_workspace=mg_workspace)
end

"""
    solve(config, prob; bc_order=:spec, output_dir="results")

Main solver loop using Taylor series pseudo-time stepping.
"""
function solve(config::SolverConfig, prob::ProblemSpec; bc_order=:spec, output_dir="results",
               lap_order::Symbol=:second,
               mg_interval::Int=0, mg_dt_scale::Real=2.0, mg_M::Int=4,
               mg_vcycle::Bool=false, mg_nu1::Int=1, mg_nu2::Int=1,
               mg_level_Ms::Union{Nothing,AbstractVector{Int}}=nothing,
               mg_level_dt_scales::Union{Nothing,AbstractVector{<:Real}}=nothing,
               mg_correction::Symbol=:classic,
               mg_corr_M::Int=2, mg_corr_dt_scale::Real=1.0, mg_corr_steps::Int=1,
               mg_corr_nu1::Union{Nothing,Int}=nothing,
               mg_corr_nu2::Union{Nothing,Int}=nothing,
               debug_residual::Bool=false, debug_vcycle::Bool=false,
               mg_workspace=nothing)
    sol, _ = solve_core(config, prob; bc_order=bc_order, output_dir=output_dir,
                        lap_order=lap_order,
                        mg_interval=mg_interval, mg_dt_scale=mg_dt_scale, mg_M=mg_M,
                        mg_vcycle=mg_vcycle, mg_nu1=mg_nu1, mg_nu2=mg_nu2,
                        mg_level_Ms=mg_level_Ms, mg_level_dt_scales=mg_level_dt_scales,
                        mg_correction=mg_correction,
                        mg_corr_M=mg_corr_M, mg_corr_dt_scale=mg_corr_dt_scale, mg_corr_steps=mg_corr_steps,
                        mg_corr_nu1=mg_corr_nu1, mg_corr_nu2=mg_corr_nu2,
                        debug_residual=debug_residual, debug_vcycle=debug_vcycle,
                        mg_workspace=mg_workspace)
    return sol
end
