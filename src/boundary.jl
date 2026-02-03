# boundary.jl

"""
    apply_bc!(u, bc, m, config; Lx=1, Ly=1, Lz=1, order=:spec)

Apply Dirichlet boundary conditions to ghost cells.
- order=:spec (default)
  - m == 0: u_ghost = 2g - u_adj
  - m >= 1: u_ghost = -u_adj
- order=:high (m==0 only): cubic extrapolation to improve boundary accuracy.

Ghost updates are applied on faces only.
"""
function apply_bc!(u::Array{T,3}, bc::BoundaryConditions, m::Int, config::SolverConfig;
                   Lx=1, Ly=1, Lz=1, order=:spec) where {T}
    dx = Lx / config.nx
    dy = Ly / config.ny
    dz = Lz / config.nz

    use_high = (order === :high) && (m == 0) &&
               (config.nx >= 3) && (config.ny >= 3) && (config.nz >= 3)

    if use_high
        # x-min / x-max faces (varying y,z)
        @inbounds for k in 2:config.nz+1
            z = (k - 1.5) * dz
            for j in 2:config.ny+1
                y = (j - 1.5) * dy
                g_lo = bc.g_xlo(y, z)
                g_hi = bc.g_xhi(y, z)
                u[1, j, k] = (16 / 5) * g_lo - 3 * u[2, j, k] + u[3, j, k] - (1 / 5) * u[4, j, k]
                u[config.nx+2, j, k] = (16 / 5) * g_hi - 3 * u[config.nx+1, j, k] +
                                       u[config.nx, j, k] - (1 / 5) * u[config.nx-1, j, k]
            end
        end

        # y-min / y-max faces (varying x,z)
        @inbounds for k in 2:config.nz+1
            z = (k - 1.5) * dz
            for i in 2:config.nx+1
                x = (i - 1.5) * dx
                g_lo = bc.g_ylo(x, z)
                g_hi = bc.g_yhi(x, z)
                u[i, 1, k] = (16 / 5) * g_lo - 3 * u[i, 2, k] + u[i, 3, k] - (1 / 5) * u[i, 4, k]
                u[i, config.ny+2, k] = (16 / 5) * g_hi - 3 * u[i, config.ny+1, k] +
                                       u[i, config.ny, k] - (1 / 5) * u[i, config.ny-1, k]
            end
        end

        # z-min / z-max faces (varying x,y)
        @inbounds for j in 2:config.ny+1
            y = (j - 1.5) * dy
            for i in 2:config.nx+1
                x = (i - 1.5) * dx
                g_lo = bc.g_zlo(x, y)
                g_hi = bc.g_zhi(x, y)
                u[i, j, 1] = (16 / 5) * g_lo - 3 * u[i, j, 2] + u[i, j, 3] - (1 / 5) * u[i, j, 4]
                u[i, j, config.nz+2] = (16 / 5) * g_hi - 3 * u[i, j, config.nz+1] +
                                       u[i, j, config.nz] - (1 / 5) * u[i, j, config.nz-1]
            end
        end
    else
        scale = m == 0 ? one(T) : zero(T)

        # x-min / x-max faces (varying y,z)
        @inbounds for k in 2:config.nz+1
            z = (k - 1.5) * dz
            for j in 2:config.ny+1
                y = (j - 1.5) * dy
                g_lo = bc.g_xlo(y, z)
                g_hi = bc.g_xhi(y, z)
                u[1, j, k] = 2 * scale * g_lo - u[2, j, k]
                u[config.nx+2, j, k] = 2 * scale * g_hi - u[config.nx+1, j, k]
            end
        end

        # y-min / y-max faces (varying x,z)
        @inbounds for k in 2:config.nz+1
            z = (k - 1.5) * dz
            for i in 2:config.nx+1
                x = (i - 1.5) * dx
                g_lo = bc.g_ylo(x, z)
                g_hi = bc.g_yhi(x, z)
                u[i, 1, k] = 2 * scale * g_lo - u[i, 2, k]
                u[i, config.ny+2, k] = 2 * scale * g_hi - u[i, config.ny+1, k]
            end
        end

        # z-min / z-max faces (varying x,y)
        @inbounds for j in 2:config.ny+1
            y = (j - 1.5) * dy
            for i in 2:config.nx+1
                x = (i - 1.5) * dx
                g_lo = bc.g_zlo(x, y)
                g_hi = bc.g_zhi(x, y)
                u[i, j, 1] = 2 * scale * g_lo - u[i, j, 2]
                u[i, j, config.nz+2] = 2 * scale * g_hi - u[i, j, config.nz+1]
            end
        end
    end

    return u
end
