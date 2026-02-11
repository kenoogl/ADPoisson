# visualization.jl

using Plots

"""
    plot_slice(sol, prob, config; y_index=nothing, output_dir="results")

Plot XZ slice at y=0.5 (nearest index) and save figures.
Saves the analytical solution slice and the error slice.
"""
function plot_slice(sol::Solution, prob::ProblemSpec, config::SolverConfig;
                    y_index=nothing, output_dir="results")
    nx = length(sol.x)
    ny = length(sol.y)
    nz = length(sol.z)

    if y_index === nothing
        y_target = 0.5 * prob.Ly
        y_index = argmin(abs.(sol.y .- y_target))
    end
    gx = (size(sol.u, 1) - nx) ÷ 2
    gy = (size(sol.u, 2) - ny) ÷ 2
    gz = (size(sol.u, 3) - nz) ÷ 2

    Uex = Array{eltype(sol.u)}(undef, nx, nz)
    Uerr = Array{eltype(sol.u)}(undef, nx, nz)
    for k in 1:nz
        z = sol.z[k]
        for i in 1:nx
            x = sol.x[i]
            u_num = sol.u[i + gx, y_index + gy, k + gz]
            u_ex = exact_solution(x, sol.y[y_index], z, prob.alpha)
            Uex[i, k] = u_ex
            Uerr[i, k] = abs(u_num - u_ex)
        end
    end

    isdir(output_dir) || mkpath(output_dir)
    p_exact = contourf(sol.x, sol.z, Uex';
                       aspect_ratio=:equal,
                       xlims=(-0.1, 1.1),
                       right_margin=10Plots.mm,
                       colorbar_formatter=:scientific,
                       title="u_exact(x,z) slice")
    p_err = contourf(sol.x, sol.z, Uerr';
                     aspect_ratio=:equal,
                     xlims=(-0.1, 1.1),
                     right_margin=10Plots.mm,
                     colorbar_formatter=:scientific,
                     title="|u - u_exact|")

    tag = "nx$(nx)_ny$(ny)_nz$(nz)_M$(config.M)_steps$(sol.iter)"
    exact_tag = "nx$(nx)_ny$(ny)_nz$(nz)"
    png(p_exact, joinpath(output_dir, "exact_$(exact_tag).png"))
    png(p_err, joinpath(output_dir, "error_$(tag).png"))
    return p_exact, p_err
end
