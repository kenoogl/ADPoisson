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
