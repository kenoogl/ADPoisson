# mg.jl

"""
    pseudo_mg_correction!(u, f, bc, config, prob; dt_scale=2.0, M=2, bc_order=:spec)

Apply a pseudo multigrid-style correction by taking a stronger Taylor step
with a larger dt (clipped by Fo) and low order M.
"""
function pseudo_mg_correction!(u::Array{T,3}, f::Array{T,3}, bc::BoundaryConditions,
                               config::SolverConfig, prob::ProblemSpec;
                               dt_scale::T=convert(T, 2.0), M::Int=2,
                               bc_order::Symbol=:spec) where {T<:Real}
    r = similar(u)
    buffers = TaylorBuffers3D(similar(u), similar(u), similar(u))

    apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
    compute_residual_norm!(r, u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)

    dx, dy, dz = grid_spacing(config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    denom = (one(T) / (dx * dx) + one(T) / (dy * dy) + one(T) / (dz * dz))
    dt_corr = config.dt * dt_scale
    if dt_corr * denom > 0.5
        dt_corr = convert(T, 0.5) / denom
    end

    corr_config = SolverConfig(config.nx, config.ny, config.nz, M, dt_corr,
                               config.max_steps, config.epsilon)
    taylor_series_update_reuse!(u, buffers, r, f, bc, corr_config, prob; bc_order=bc_order)
    return u
end

function taylor_smoother!(u::Array{T,3}, buffers::TaylorBuffers3D{T}, r::Array{T,3},
                          f::Array{T,3}, bc::BoundaryConditions,
                          config::SolverConfig, prob::ProblemSpec;
                          steps::Int=1, bc_order::Symbol=:spec) where {T<:Real}
    for _ in 1:steps
        apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)
        compute_residual_norm!(r, u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
        taylor_series_update_reuse!(u, buffers, r, f, bc, config, prob; bc_order=bc_order)
    end
    return u
end

function zero_ghost!(a::Array{T,3}, config::SolverConfig) where {T<:Real}
    nx = config.nx
    ny = config.ny
    nz = config.nz
    @inbounds for k in 1:nz+2, j in 1:ny+2
        a[1, j, k] = zero(T)
        a[nx+2, j, k] = zero(T)
    end
    @inbounds for k in 1:nz+2, i in 1:nx+2
        a[i, 1, k] = zero(T)
        a[i, ny+2, k] = zero(T)
    end
    @inbounds for j in 1:ny+2, i in 1:nx+2
        a[i, j, 1] = zero(T)
        a[i, j, nz+2] = zero(T)
    end
    return a
end

"""
    restrict_full_weighting!(rc, rf, cfg_f, cfg_c)

3D full-weighting restriction from fine residual rf to coarse rc (cell-centered).
Assumes rc and rf include ghost cells.
"""
function restrict_full_weighting!(rc::Array{T,3}, rf::Array{T,3},
                                  cfg_f::SolverConfig, cfg_c::SolverConfig) where {T<:Real}
    w = (convert(T, 1 / 8), convert(T, 3 / 8), convert(T, 3 / 8), convert(T, 1 / 8))
    @inbounds for K in 1:cfg_c.nz, J in 1:cfg_c.ny, I in 1:cfg_c.nx
        ic = I + 1
        jc = J + 1
        kc = K + 1
        ifi = 2 * I
        jfi = 2 * J
        kfi = 2 * K
        acc = zero(T)
        for dk in -1:2, dj in -1:2, di in -1:2
            acc += w[di + 2] * w[dj + 2] * w[dk + 2] * rf[ifi + di, jfi + dj, kfi + dk]
        end
        rc[ic, jc, kc] = acc
    end
    return rc
end

"""
    prolong_trilinear!(ef, ec, cfg_f, cfg_c)

3D trilinear prolongation from coarse ec to fine ef (cell-centered).
Assumes arrays include ghost cells.
"""
function prolong_trilinear!(ef::Array{T,3}, ec::Array{T,3},
                            cfg_f::SolverConfig, cfg_c::SolverConfig) where {T<:Real}
    @inbounds for kf in 1:cfg_f.nz, jf in 1:cfg_f.ny, ifine in 1:cfg_f.nx
        x = (ifine + 0.5) / 2
        y = (jf + 0.5) / 2
        z = (kf + 0.5) / 2
        I0 = clamp(Int(floor(x)), 0, cfg_c.nx)
        J0 = clamp(Int(floor(y)), 0, cfg_c.ny)
        K0 = clamp(Int(floor(z)), 0, cfg_c.nz)
        I1 = I0 + 1
        J1 = J0 + 1
        K1 = K0 + 1
        tx = x - I0
        ty = y - J0
        tz = z - K0

        v000 = ec[I0 + 1, J0 + 1, K0 + 1]
        v100 = ec[I1 + 1, J0 + 1, K0 + 1]
        v010 = ec[I0 + 1, J1 + 1, K0 + 1]
        v110 = ec[I1 + 1, J1 + 1, K0 + 1]
        v001 = ec[I0 + 1, J0 + 1, K1 + 1]
        v101 = ec[I1 + 1, J0 + 1, K1 + 1]
        v011 = ec[I0 + 1, J1 + 1, K1 + 1]
        v111 = ec[I1 + 1, J1 + 1, K1 + 1]

        v00 = v000 * (1 - tx) + v100 * tx
        v10 = v010 * (1 - tx) + v110 * tx
        v01 = v001 * (1 - tx) + v101 * tx
        v11 = v011 * (1 - tx) + v111 * tx
        v0 = v00 * (1 - ty) + v10 * ty
        v1 = v01 * (1 - ty) + v11 * ty
        ef[ifine + 1, jf + 1, kf + 1] = v0 * (1 - tz) + v1 * tz
    end
    return ef
end

function mg_max_levels(config::SolverConfig; min_n::Int=4)
    levels = 1
    nx = config.nx
    ny = config.ny
    nz = config.nz
    while (nx % 2 == 0) && (ny % 2 == 0) && (nz % 2 == 0)
        nx2 = nx ÷ 2
        ny2 = ny ÷ 2
        nz2 = nz ÷ 2
        if nx2 < min_n || ny2 < min_n || nz2 < min_n
            break
        end
        levels += 1
        nx = nx2
        ny = ny2
        nz = nz2
    end
    return levels
end

struct MGLevelWorkspace{T}
    r::Array{T,3}
    rhs::Array{T,3}
    e::Array{T,3}
    taylor::TaylorBuffers3D{T}
    tmp::Array{T,3}
end

mutable struct MGWorkspace{T}
    levels::Vector{MGLevelWorkspace{T}}
    coarse_A::Union{Nothing,Array{T,2}}
    coarse_fact::Union{Nothing,Factorization{T}}
    coarse_b::Union{Nothing,Vector{T}}
    coarse_x::Union{Nothing,Vector{T}}
    coarse_dims::NTuple{3,Int}
end

function mg_level_workspace(::Type{T}, nx::Int, ny::Int, nz::Int) where {T<:Real}
    r = Array{T}(undef, nx + 2, ny + 2, nz + 2)
    rhs = Array{T}(undef, nx + 2, ny + 2, nz + 2)
    e = Array{T}(undef, nx + 2, ny + 2, nz + 2)
    tmp = Array{T}(undef, nx + 2, ny + 2, nz + 2)
    taylor = TaylorBuffers3D(Array{T}(undef, nx + 2, ny + 2, nz + 2),
                             Array{T}(undef, nx + 2, ny + 2, nz + 2),
                             Array{T}(undef, nx + 2, ny + 2, nz + 2))
    return MGLevelWorkspace{T}(r, rhs, e, taylor, tmp)
end

function build_mg_workspace(u::Array{T,3}, config::SolverConfig; min_n::Int=4) where {T<:Real}
    levels = mg_max_levels(config; min_n=min_n)
    work_levels = Vector{MGLevelWorkspace{T}}(undef, levels)
    nx = config.nx
    ny = config.ny
    nz = config.nz
    for level in 1:levels
        work_levels[level] = mg_level_workspace(T, nx, ny, nz)
        if level < levels
            nx = nx ÷ 2
            ny = ny ÷ 2
            nz = nz ÷ 2
        end
    end
    return MGWorkspace{T}(work_levels, nothing, nothing, nothing, nothing, (0, 0, 0))
end

@inline function mg_level_value(values::Nothing, level::Int, default)
    return default
end

@inline function mg_level_value(values::AbstractVector, level::Int, default)
    if level <= length(values)
        return values[level]
    end
    return values[end]
end

function mg_dt_clipped(config::SolverConfig, prob::ProblemSpec, dt_scale::T) where {T<:Real}
    dx, dy, dz = grid_spacing(config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    denom = (one(T) / (dx * dx) + one(T) / (dy * dy) + one(T) / (dz * dz))
    dt_corr = config.dt * dt_scale
    if dt_corr * denom > 0.5
        dt_corr = convert(T, 0.5) / denom
    end
    return dt_corr
end

"""
    correction_taylor_solve!(e, r, bc, config, prob; M=2, dt_scale=1.0, steps=1, bc_order=:spec, buffers=nothing, work=nothing)

Solve the correction equation L e = r using Taylor pseudo-time integration
for a fixed number of steps, starting from e = 0.
"""
function correction_taylor_solve!(e::Array{T,3}, r::Array{T,3}, bc::BoundaryConditions,
                                  config::SolverConfig, prob::ProblemSpec;
                                  M::Int=2, dt_scale::Real=1.0, steps::Int=1,
                                  bc_order::Symbol=:spec,
                                  buffers::Union{Nothing,TaylorBuffers3D{T}}=nothing,
                                  work::Union{Nothing,Array{T,3}}=nothing) where {T<:Real}
    buffers = buffers === nothing ? TaylorBuffers3D(similar(e), similar(e), similar(e)) : buffers
    work = work === nothing ? similar(e) : work
    fill!(e, zero(T))
    dt_corr = mg_dt_clipped(config, prob, convert(T, dt_scale))
    M_corr = max(M, 1)
    steps_corr = max(steps, 1)
    cfg_corr = SolverConfig(config.nx, config.ny, config.nz, M_corr, dt_corr,
                            config.max_steps, config.epsilon)
    taylor_smoother!(e, buffers, work, r, bc, cfg_corr, prob;
                     steps=steps_corr, bc_order=bc_order)
    return e
end

function direct_solve_poisson!(u::Array{T,3}, f::Array{T,3},
                               config::SolverConfig, prob::ProblemSpec,
                               ws::MGWorkspace{T}) where {T<:Real}
    nx = config.nx
    ny = config.ny
    nz = config.nz
    n = nx * ny * nz
    dims = (nx, ny, nz)
    rebuild = ws.coarse_A === nothing || ws.coarse_fact === nothing || ws.coarse_dims != dims
    if rebuild
        A = zeros(T, n, n)
        dx, dy, dz = grid_spacing(config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
        inv_dx2 = one(T) / (dx * dx)
        inv_dy2 = one(T) / (dy * dy)
        inv_dz2 = one(T) / (dz * dz)
        diag = -2 * (inv_dx2 + inv_dy2 + inv_dz2)
        stride_xy = nx * ny
        @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
            idx = i + (j - 1) * nx + (k - 1) * stride_xy
            A[idx, idx] = diag
            if i > 1
                A[idx, idx - 1] = inv_dx2
            end
            if i < nx
                A[idx, idx + 1] = inv_dx2
            end
            if j > 1
                A[idx, idx - nx] = inv_dy2
            end
            if j < ny
                A[idx, idx + nx] = inv_dy2
            end
            if k > 1
                A[idx, idx - stride_xy] = inv_dz2
            end
            if k < nz
                A[idx, idx + stride_xy] = inv_dz2
            end
        end
        ws.coarse_A = A
        ws.coarse_fact = lu(A)
        ws.coarse_b = Vector{T}(undef, n)
        ws.coarse_x = Vector{T}(undef, n)
        ws.coarse_dims = dims
    elseif ws.coarse_b === nothing || ws.coarse_x === nothing || length(ws.coarse_b) != n
        ws.coarse_b = Vector{T}(undef, n)
        ws.coarse_x = Vector{T}(undef, n)
        ws.coarse_dims = dims
    end
    b = ws.coarse_b::Vector{T}
    x = ws.coarse_x::Vector{T}
    stride_xy = nx * ny
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        idx = i + (j - 1) * nx + (k - 1) * stride_xy
        b[idx] = f[i + 1, j + 1, k + 1]
    end
    ldiv!(x, ws.coarse_fact, b)
    @inbounds for k in 1:nz, j in 1:ny, i in 1:nx
        idx = i + (j - 1) * nx + (k - 1) * stride_xy
        u[i + 1, j + 1, k + 1] = x[idx]
    end
    return u
end

"""
    vcycle!(u, f, bc, config, prob, ws;
            level=1, max_level=0, min_n=4,
            nu1=1, nu2=1, dt_scale=2.0, M=2,
            level_dt_scales=nothing, level_Ms=nothing,
            correction_mode=:classic, corr_M=2, corr_dt_scale=1.0,
            corr_steps=1, corr_nu1=nothing, corr_nu2=nothing,
            is_correction=false, bc_order=:spec)

Apply a recursive V-cycle to solve Lu = f using Taylor smoothing.
"""
function vcycle!(u::Array{T,3}, f::Array{T,3}, bc::BoundaryConditions,
                 config::SolverConfig, prob::ProblemSpec, ws::MGWorkspace{T};
                 level::Int=1, max_level::Int=0, min_n::Int=4,
                 nu1::Int=1, nu2::Int=1, dt_scale::Real=2.0, M::Int=2,
                 level_dt_scales=nothing, level_Ms=nothing,
                 correction_mode::Symbol=:classic,
                 corr_M::Int=2, corr_dt_scale::Real=1.0, corr_steps::Int=1,
                 corr_nu1::Union{Nothing,Int}=nothing,
                 corr_nu2::Union{Nothing,Int}=nothing,
                 is_correction::Bool=false,
                 bc_order::Symbol=:spec,
                 debug_io::Union{Nothing,IO}=nothing,
                 debug_denom::Union{Nothing,Real}=nothing) where {T<:Real}
    if max_level == 0
        max_level = length(ws.levels)
    end
    if max_level > length(ws.levels)
        max_level = length(ws.levels)
    end
    ws_level = ws.levels[level]
    bc_order_level = level == 1 ? bc_order : :spec
    corr_nu1_eff = (corr_nu1 === nothing) ? corr_steps : corr_nu1
    corr_nu2_eff = (corr_nu2 === nothing) ? corr_steps : corr_nu2
    use_corr_params = is_correction && (correction_mode === :correction_taylor)
    dt_scale_level = use_corr_params ? corr_dt_scale : mg_level_value(level_dt_scales, level, dt_scale)
    M_level = use_corr_params ? corr_M : mg_level_value(level_Ms, level, M)
    M_level = max(Int(M_level), 1)
    nu1_level = use_corr_params ? corr_nu1_eff : nu1
    nu2_level = use_corr_params ? corr_nu2_eff : nu2
    dt_corr = mg_dt_clipped(config, prob, convert(T, dt_scale_level))
    cfg_level = SolverConfig(config.nx, config.ny, config.nz, M_level,
                             dt_corr, config.max_steps, config.epsilon)

    nx_c = config.nx ÷ 2
    ny_c = config.ny ÷ 2
    nz_c = config.nz ÷ 2
    can_coarsen = (level < max_level) &&
                  (config.nx % 2 == 0) && (config.ny % 2 == 0) && (config.nz % 2 == 0) &&
                  (nx_c >= min_n) && (ny_c >= min_n) && (nz_c >= min_n)

    if !can_coarsen
        if is_correction && (correction_mode === :correction_taylor)
            correction_taylor_solve!(u, f, bc, config, prob;
                                     M=corr_M, dt_scale=corr_dt_scale,
                                     steps=corr_nu1_eff + corr_nu2_eff,
                                     bc_order=:spec,
                                     buffers=ws_level.taylor,
                                     work=ws_level.tmp)
        else
            direct_solve_poisson!(u, f, config, prob, ws)
        end
        apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order_level)
        return u
    end

    r = ws_level.r
    buffers = ws_level.taylor
    taylor_smoother!(u, buffers, r, f, bc, cfg_level, prob; steps=nu1_level, bc_order=bc_order_level)

    compute_residual_norm!(r, u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    if debug_io !== nothing
        res = l2_norm_interior(r, config)
        if debug_denom !== nothing
            res /= debug_denom
        end
        @printf(debug_io, "%d %d %d %.6e\n", level, config.nx, config.ny, res)
    end
    @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
        r[i, j, k] = -r[i, j, k]
    end
    zero_ghost!(r, config)

    ws_coarse = ws.levels[level + 1]
    rc = ws_coarse.rhs
    fill!(rc, zero(T))
    cfg_c_base = SolverConfig(nx_c, ny_c, nz_c, config.M, config.dt,
                              config.max_steps, config.epsilon)
    restrict_full_weighting!(rc, r, config, cfg_c_base)

    ec = ws_coarse.e
    fill!(ec, zero(T))
    bc0 = zero_boundary_conditions(T)
    vcycle!(ec, rc, bc0, cfg_c_base, prob, ws;
            level=level + 1, max_level=max_level, min_n=min_n,
            nu1=nu1, nu2=nu2, dt_scale=dt_scale, M=M,
            level_dt_scales=level_dt_scales, level_Ms=level_Ms,
            correction_mode=correction_mode,
            corr_M=corr_M, corr_dt_scale=corr_dt_scale,
            corr_steps=corr_steps, corr_nu1=corr_nu1, corr_nu2=corr_nu2,
            is_correction=true,
            bc_order=bc_order, debug_io=debug_io, debug_denom=debug_denom)

    ef = ws_level.e
    prolong_trilinear!(ef, ec, config, cfg_c_base)
    @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
        u[i, j, k] += ef[i, j, k]
    end
    apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order_level)

    taylor_smoother!(u, buffers, r, f, bc, cfg_level, prob; steps=nu2_level, bc_order=bc_order_level)
    return u
end

function vcycle!(u::Array{T,3}, f::Array{T,3}, bc::BoundaryConditions,
                 config::SolverConfig, prob::ProblemSpec;
                 level::Int=1, max_level::Int=0, min_n::Int=4,
                 nu1::Int=1, nu2::Int=1, dt_scale::Real=2.0, M::Int=2,
                 level_dt_scales=nothing, level_Ms=nothing,
                 correction_mode::Symbol=:classic,
                 corr_M::Int=2, corr_dt_scale::Real=1.0, corr_steps::Int=1,
                 corr_nu1::Union{Nothing,Int}=nothing,
                 corr_nu2::Union{Nothing,Int}=nothing,
                 is_correction::Bool=false,
                 bc_order::Symbol=:spec,
                 debug_io::Union{Nothing,IO}=nothing,
                 debug_denom::Union{Nothing,Real}=nothing) where {T<:Real}
    ws = build_mg_workspace(u, config; min_n=min_n)
    return vcycle!(u, f, bc, config, prob, ws;
                   level=level, max_level=max_level, min_n=min_n,
                   nu1=nu1, nu2=nu2, dt_scale=dt_scale, M=M,
                   level_dt_scales=level_dt_scales, level_Ms=level_Ms,
                   correction_mode=correction_mode,
                   corr_M=corr_M, corr_dt_scale=corr_dt_scale, corr_steps=corr_steps,
                   corr_nu1=corr_nu1, corr_nu2=corr_nu2,
                   is_correction=is_correction,
                   bc_order=bc_order, debug_io=debug_io, debug_denom=debug_denom)
end

"""
    two_level_mg_correction!(u, f, bc, config, prob;
                             nu1=1, nu2=1, dt_scale=2.0, M=2, bc_order=:spec)

Apply a 2-level MG correction using full-weighting restriction and trilinear prolongation.
"""
function two_level_mg_correction!(u::Array{T,3}, f::Array{T,3}, bc::BoundaryConditions,
                                  config::SolverConfig, prob::ProblemSpec;
                                  nu1::Int=1, nu2::Int=1, dt_scale::Real=2.0, M::Int=2,
                                  bc_order::Symbol=:spec) where {T<:Real}
    nx_c = Int(floor(config.nx / 2))
    ny_c = Int(floor(config.ny / 2))
    nz_c = Int(floor(config.nz / 2))
    if nx_c < 4 || ny_c < 4 || nz_c < 4
        return u
    end

    r = similar(u)
    buffers_f = TaylorBuffers3D(similar(u), similar(u), similar(u))
    taylor_smoother!(u, buffers_f, r, f, bc, config, prob; steps=nu1, bc_order=bc_order)

    compute_residual_norm!(r, u, f, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
        r[i, j, k] = -r[i, j, k]
    end
    zero_ghost!(r, config)

    rc = zeros(T, nx_c + 2, ny_c + 2, nz_c + 2)
    cfg_c_base = SolverConfig(nx_c, ny_c, nz_c, M, config.dt,
                              config.max_steps, config.epsilon)
    restrict_full_weighting!(rc, r, config, cfg_c_base)

    ec = zeros(T, nx_c + 2, ny_c + 2, nz_c + 2)
    buffers_c = TaylorBuffers3D(similar(ec), similar(ec), similar(ec))
    bc0 = zero_boundary_conditions(T)

    dx, dy, dz = grid_spacing(cfg_c_base; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    denom = (one(T) / (dx * dx) + one(T) / (dy * dy) + one(T) / (dz * dz))
    dt_corr = convert(T, dt_scale) * config.dt
    if dt_corr * denom > 0.5
        dt_corr = convert(T, 0.5) / denom
    end
    cfg_c = SolverConfig(nx_c, ny_c, nz_c, M, dt_corr, config.max_steps, config.epsilon)

    rc_buf = similar(ec)
    taylor_smoother!(ec, buffers_c, rc_buf, rc, bc0, cfg_c, prob; steps=1, bc_order=bc_order)

    ef = zeros(T, config.nx + 2, config.ny + 2, config.nz + 2)
    prolong_trilinear!(ef, ec, config, cfg_c)

    @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
        u[i, j, k] += ef[i, j, k]
    end
    apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)

    taylor_smoother!(u, buffers_f, r, f, bc, config, prob; steps=nu2, bc_order=bc_order)
    return u
end
