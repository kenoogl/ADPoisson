# visualization.jl

using Plots

"""
    plot_slice(sol, prob, config; y_index=nothing, output_dir="results")

Plot XZ slice at y=0.5 (nearest index) and save figures.
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

    U = Array{eltype(sol.u)}(undef, nx, nz)
    Uerr = Array{eltype(sol.u)}(undef, nx, nz)
    for k in 1:nz
        z = sol.z[k]
        for i in 1:nx
            x = sol.x[i]
            u_num = sol.u[i + 1, y_index + 1, k + 1]
            u_ex = exact_solution(x, sol.y[y_index], z, prob.alpha)
            U[i, k] = u_num
            Uerr[i, k] = abs(u_num - u_ex)
        end
    end

    isdir(output_dir) || mkpath(output_dir)
    p1 = contourf(sol.x, sol.z, U';
                  aspect_ratio=:equal,
                  xlims=(-0.1, 1.1),
                  right_margin=10Plots.mm,
                  colorbar_formatter=:scientific,
                  title="u(x,z) slice")
    p2 = contourf(sol.x, sol.z, Uerr';
                  aspect_ratio=:equal,
                  xlims=(-0.1, 1.1),
                  right_margin=10Plots.mm,
                  colorbar_formatter=:scientific,
                  title="|u - u_exact|")

    tag = "nx$(nx)_ny$(ny)_nz$(nz)_M$(config.M)_t$(config.tend)"
    png(p1, joinpath(output_dir, "solution_$(tag).png"))
    png(p2, joinpath(output_dir, "error_$(tag).png"))
    return p1, p2
end
