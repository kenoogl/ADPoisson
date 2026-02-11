# boundary.jl

"""
    apply_bc!(u, bc, m, config; Lx=1, Ly=1, Lz=1, order=:spec)

Apply Dirichlet boundary conditions to ghost cells.
- order=:spec (default)
  - m == 0: u_ghost = 2g - u_adj
  - m >= 1: u_ghost = -u_adj
- order=:high (ghost 2-layer only): 4th-order boundary extrapolation.

Ghost updates are applied on faces only.
"""
function apply_bc!(u::Array{T,3}, bc::BoundaryConditions, m::Int, config::SolverConfig;
                   Lx=1, Ly=1, Lz=1, order=:spec) where {T}
    dx = Lx / config.nx
    dy = Ly / config.ny
    dz = Lz / config.nz

    ghost = (size(u, 1) - config.nx) ÷ 2
    if ghost < 1 || ghost > 2 ||
       size(u, 1) != config.nx + 2 * ghost ||
       size(u, 2) != config.ny + 2 * ghost ||
       size(u, 3) != config.nz + 2 * ghost
        error("apply_bc! expects ghost layers = 1 or 2 with consistent array sizes.")
    end
    i_lo = ghost + 1
    i_hi = config.nx + ghost
    j_lo = ghost + 1
    j_hi = config.ny + ghost
    k_lo = ghost + 1
    k_hi = config.nz + ghost
    coord_offset = ghost + 0.5

    use_high4 = (order === :high) && (ghost == 2) &&
                (config.nx >= 4) && (config.ny >= 4) && (config.nz >= 4)
    use_high_legacy = (order === :high) && (ghost == 1) && (m == 0) &&
                      (config.nx >= 4) && (config.ny >= 4) && (config.nz >= 4)

    if use_high4
        g_scale = m == 0 ? one(T) : zero(T)

        # x-min / x-max faces (varying y,z)
        @inbounds for k in k_lo:k_hi
            z = (k - coord_offset) * dz
            for j in j_lo:j_hi
                y = (j - coord_offset) * dy
                g_lo = bc.g_xlo(y, z) * g_scale
                g_hi = bc.g_xhi(y, z) * g_scale
                u[i_lo-1, j, k] = (128 / 35) * g_lo - 4 * u[i_lo, j, k] + 2 * u[i_lo+1, j, k] -
                                  (4 / 5) * u[i_lo+2, j, k] + (1 / 7) * u[i_lo+3, j, k]
                u[i_lo-2, j, k] = (128 / 7) * g_lo - 30 * u[i_lo, j, k] + 20 * u[i_lo+1, j, k] -
                                  9 * u[i_lo+2, j, k] + (12 / 7) * u[i_lo+3, j, k]
                u[i_hi+1, j, k] = (128 / 35) * g_hi - 4 * u[i_hi, j, k] + 2 * u[i_hi-1, j, k] -
                                  (4 / 5) * u[i_hi-2, j, k] + (1 / 7) * u[i_hi-3, j, k]
                u[i_hi+2, j, k] = (128 / 7) * g_hi - 30 * u[i_hi, j, k] + 20 * u[i_hi-1, j, k] -
                                  9 * u[i_hi-2, j, k] + (12 / 7) * u[i_hi-3, j, k]
            end
        end

        # y-min / y-max faces (varying x,z)
        @inbounds for k in k_lo:k_hi
            z = (k - coord_offset) * dz
            for i in i_lo:i_hi
                x = (i - coord_offset) * dx
                g_lo = bc.g_ylo(x, z) * g_scale
                g_hi = bc.g_yhi(x, z) * g_scale
                u[i, j_lo-1, k] = (128 / 35) * g_lo - 4 * u[i, j_lo, k] + 2 * u[i, j_lo+1, k] -
                                  (4 / 5) * u[i, j_lo+2, k] + (1 / 7) * u[i, j_lo+3, k]
                u[i, j_lo-2, k] = (128 / 7) * g_lo - 30 * u[i, j_lo, k] + 20 * u[i, j_lo+1, k] -
                                  9 * u[i, j_lo+2, k] + (12 / 7) * u[i, j_lo+3, k]
                u[i, j_hi+1, k] = (128 / 35) * g_hi - 4 * u[i, j_hi, k] + 2 * u[i, j_hi-1, k] -
                                  (4 / 5) * u[i, j_hi-2, k] + (1 / 7) * u[i, j_hi-3, k]
                u[i, j_hi+2, k] = (128 / 7) * g_hi - 30 * u[i, j_hi, k] + 20 * u[i, j_hi-1, k] -
                                  9 * u[i, j_hi-2, k] + (12 / 7) * u[i, j_hi-3, k]
            end
        end

        # z-min / z-max faces (varying x,y)
        @inbounds for j in j_lo:j_hi
            y = (j - coord_offset) * dy
            for i in i_lo:i_hi
                x = (i - coord_offset) * dx
                g_lo = bc.g_zlo(x, y) * g_scale
                g_hi = bc.g_zhi(x, y) * g_scale
                u[i, j, k_lo-1] = (128 / 35) * g_lo - 4 * u[i, j, k_lo] + 2 * u[i, j, k_lo+1] -
                                  (4 / 5) * u[i, j, k_lo+2] + (1 / 7) * u[i, j, k_lo+3]
                u[i, j, k_lo-2] = (128 / 7) * g_lo - 30 * u[i, j, k_lo] + 20 * u[i, j, k_lo+1] -
                                  9 * u[i, j, k_lo+2] + (12 / 7) * u[i, j, k_lo+3]
                u[i, j, k_hi+1] = (128 / 35) * g_hi - 4 * u[i, j, k_hi] + 2 * u[i, j, k_hi-1] -
                                  (4 / 5) * u[i, j, k_hi-2] + (1 / 7) * u[i, j, k_hi-3]
                u[i, j, k_hi+2] = (128 / 7) * g_hi - 30 * u[i, j, k_hi] + 20 * u[i, j, k_hi-1] -
                                  9 * u[i, j, k_hi-2] + (12 / 7) * u[i, j, k_hi-3]
            end
        end
    elseif use_high_legacy
        # Legacy (ghost 1-layer): cubic extrapolation at m == 0
        @inbounds for k in k_lo:k_hi
            z = (k - coord_offset) * dz
            for j in j_lo:j_hi
                y = (j - coord_offset) * dy
                g_lo = bc.g_xlo(y, z)
                g_hi = bc.g_xhi(y, z)
                u[i_lo-1, j, k] = (16 / 5) * g_lo - 3 * u[i_lo, j, k] + u[i_lo+1, j, k] -
                                  (1 / 5) * u[i_lo+2, j, k]
                u[i_hi+1, j, k] = (16 / 5) * g_hi - 3 * u[i_hi, j, k] + u[i_hi-1, j, k] -
                                  (1 / 5) * u[i_hi-2, j, k]
            end
        end

        @inbounds for k in k_lo:k_hi
            z = (k - coord_offset) * dz
            for i in i_lo:i_hi
                x = (i - coord_offset) * dx
                g_lo = bc.g_ylo(x, z)
                g_hi = bc.g_yhi(x, z)
                u[i, j_lo-1, k] = (16 / 5) * g_lo - 3 * u[i, j_lo, k] + u[i, j_lo+1, k] -
                                  (1 / 5) * u[i, j_lo+2, k]
                u[i, j_hi+1, k] = (16 / 5) * g_hi - 3 * u[i, j_hi, k] + u[i, j_hi-1, k] -
                                  (1 / 5) * u[i, j_hi-2, k]
            end
        end

        @inbounds for j in j_lo:j_hi
            y = (j - coord_offset) * dy
            for i in i_lo:i_hi
                x = (i - coord_offset) * dx
                g_lo = bc.g_zlo(x, y)
                g_hi = bc.g_zhi(x, y)
                u[i, j, k_lo-1] = (16 / 5) * g_lo - 3 * u[i, j, k_lo] + u[i, j, k_lo+1] -
                                  (1 / 5) * u[i, j, k_lo+2]
                u[i, j, k_hi+1] = (16 / 5) * g_hi - 3 * u[i, j, k_hi] + u[i, j, k_hi-1] -
                                  (1 / 5) * u[i, j, k_hi-2]
            end
        end
    else
        scale = m == 0 ? one(T) : zero(T)

        # x-min / x-max faces (varying y,z)
        @inbounds for k in k_lo:k_hi
            z = (k - coord_offset) * dz
            for j in j_lo:j_hi
                y = (j - coord_offset) * dy
                g_lo = bc.g_xlo(y, z)
                g_hi = bc.g_xhi(y, z)
                u[i_lo-1, j, k] = 2 * scale * g_lo - u[i_lo, j, k]
                u[i_hi+1, j, k] = 2 * scale * g_hi - u[i_hi, j, k]
                if ghost == 2
                    u[i_lo-2, j, k] = 2 * scale * g_lo - u[i_lo+1, j, k]
                    u[i_hi+2, j, k] = 2 * scale * g_hi - u[i_hi-1, j, k]
                end
            end
        end

        # y-min / y-max faces (varying x,z)
        @inbounds for k in k_lo:k_hi
            z = (k - coord_offset) * dz
            for i in i_lo:i_hi
                x = (i - coord_offset) * dx
                g_lo = bc.g_ylo(x, z)
                g_hi = bc.g_yhi(x, z)
                u[i, j_lo-1, k] = 2 * scale * g_lo - u[i, j_lo, k]
                u[i, j_hi+1, k] = 2 * scale * g_hi - u[i, j_hi, k]
                if ghost == 2
                    u[i, j_lo-2, k] = 2 * scale * g_lo - u[i, j_lo+1, k]
                    u[i, j_hi+2, k] = 2 * scale * g_hi - u[i, j_hi-1, k]
                end
            end
        end

        # z-min / z-max faces (varying x,y)
        @inbounds for j in j_lo:j_hi
            y = (j - coord_offset) * dy
            for i in i_lo:i_hi
                x = (i - coord_offset) * dx
                g_lo = bc.g_zlo(x, y)
                g_hi = bc.g_zhi(x, y)
                u[i, j, k_lo-1] = 2 * scale * g_lo - u[i, j, k_lo]
                u[i, j, k_hi+1] = 2 * scale * g_hi - u[i, j, k_hi]
                if ghost == 2
                    u[i, j, k_lo-2] = 2 * scale * g_lo - u[i, j, k_lo+1]
                    u[i, j, k_hi+2] = 2 * scale * g_hi - u[i, j, k_hi-1]
                end
            end
        end
    end

    return u
end
