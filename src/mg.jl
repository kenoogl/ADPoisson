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

3D full-weighting restriction from fine residual rf to coarse rc.
Assumes rc and rf include ghost cells.
"""
function restrict_full_weighting!(rc::Array{T,3}, rf::Array{T,3},
                                  cfg_f::SolverConfig, cfg_c::SolverConfig) where {T<:Real}
    w = (convert(T, 0.25), convert(T, 0.5), convert(T, 0.25))
    @inbounds for K in 1:cfg_c.nz, J in 1:cfg_c.ny, I in 1:cfg_c.nx
        ic = I + 1
        jc = J + 1
        kc = K + 1
        ifi = 2 * I
        jfi = 2 * J
        kfi = 2 * K
        acc = zero(T)
        for dk in -1:1, dj in -1:1, di in -1:1
            acc += w[di + 2] * w[dj + 2] * w[dk + 2] * rf[ifi + di, jfi + dj, kfi + dk]
        end
        rc[ic, jc, kc] = acc
    end
    return rc
end

"""
    prolong_trilinear!(ef, ec, cfg_f, cfg_c)

3D trilinear prolongation from coarse ec to fine ef.
Assumes arrays include ghost cells.
"""
function prolong_trilinear!(ef::Array{T,3}, ec::Array{T,3},
                            cfg_f::SolverConfig, cfg_c::SolverConfig) where {T<:Real}
    @inbounds for kf in 1:cfg_f.nz, jf in 1:cfg_f.ny, ifine in 1:cfg_f.nx
        x = ifine / 2
        y = jf / 2
        z = kf / 2
        I0 = clamp(Int(floor(x)), 1, cfg_c.nx)
        J0 = clamp(Int(floor(y)), 1, cfg_c.ny)
        K0 = clamp(Int(floor(z)), 1, cfg_c.nz)
        I1 = clamp(I0 + 1, 1, cfg_c.nx)
        J1 = clamp(J0 + 1, 1, cfg_c.ny)
        K1 = clamp(K0 + 1, 1, cfg_c.nz)
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
    zero_ghost!(r, config)

    rc = zeros(T, nx_c + 2, ny_c + 2, nz_c + 2)
    cfg_c_base = SolverConfig(nx_c, ny_c, nz_c, M, config.dt,
                              config.max_steps, config.epsilon)
    restrict_full_weighting!(rc, r, config, cfg_c_base)

    ec = zeros(T, nx_c + 2, ny_c + 2, nz_c + 2)
    buffers_c = TaylorBuffers3D(similar(ec), similar(ec), similar(ec))

    dx, dy, dz = grid_spacing(cfg_c_base; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz)
    denom = (one(T) / (dx * dx) + one(T) / (dy * dy) + one(T) / (dz * dz))
    dt_corr = convert(T, dt_scale) * config.dt
    if dt_corr * denom > 0.5
        dt_corr = convert(T, 0.5) / denom
    end
    cfg_c = SolverConfig(nx_c, ny_c, nz_c, M, dt_corr, config.max_steps, config.epsilon)

    rc_buf = similar(ec)
    taylor_smoother!(ec, buffers_c, rc_buf, rc, bc, cfg_c, prob; steps=1, bc_order=bc_order)

    ef = zeros(T, config.nx + 2, config.ny + 2, config.nz + 2)
    prolong_trilinear!(ef, ec, config, cfg_c)

    @inbounds for k in 2:config.nz+1, j in 2:config.ny+1, i in 2:config.nx+1
        u[i, j, k] += ef[i, j, k]
    end
    apply_bc!(u, bc, 0, config; Lx=prob.Lx, Ly=prob.Ly, Lz=prob.Lz, order=bc_order)

    taylor_smoother!(u, buffers_f, r, f, bc, config, prob; steps=nu2, bc_order=bc_order)
    return u
end
