#!/usr/bin/env julia

using ADPoisson
using DelimitedFiles
using Printf
using Plots
using Dates
using TOML

function make_run_dir(output_dir; prefix="run")
    mkpath(output_dir)
    ts = Dates.format(now(), dateformat"yyyymmdd_HHMMSS")
    run_dir = joinpath(output_dir, "$(prefix)_$(ts)")
    i = 1
    while isdir(run_dir)
        run_dir = joinpath(output_dir, "$(prefix)_$(ts)_$(i)")
        i += 1
    end
    mkpath(run_dir)
    return run_dir
end

function parse_alpha_list(s::String)
    parts = split(s, ",")
    alphas = Float64[]
    for p in parts
        t = strip(p)
        isempty(t) && continue
        push!(alphas, parse(Float64, t))
    end
    isempty(alphas) && error("alpha list is empty")
    return alphas
end

function format_alpha_tag(alpha::Real)
    s = @sprintf("%.3g", alpha)
    s = replace(s, "." => "p", "-" => "m", "+" => "")
    return s
end

function parse_args(args)
    opts = default_cli_options()
    opts["alphas"] = "0.5,1.0,1.5"
    opts["output_dir"] = "results"
    opts["Fo"] = nothing

    i = 1
    while i <= length(args)
        if startswith(args[i], "--")
            key = args[i][3:end]
            if i == length(args)
                error("Missing value for --$key")
            end
            if key == "bc-order" || key == "bc_order"
                opts["bc_order"] = args[i + 1]
            elseif key == "max-steps" || key == "max_steps"
                opts["max_steps"] = parse(Int, args[i + 1])
            elseif key == "alphas" || key == "alpha-list" || key == "alpha_list"
                opts["alphas"] = args[i + 1]
            elseif key == "output-dir" || key == "output_dir"
                opts["output_dir"] = args[i + 1]
            elseif key == "Fo" || key == "fo"
                opts["Fo"] = parse(Float64, args[i + 1])
            else
                if key == "nx" || key == "ny" || key == "nz" || key == "M"
                    opts[key] = parse(Int, args[i + 1])
                else
                    opts[key] = parse(Float64, args[i + 1])
                end
            end
            i += 2
        else
            i += 1
        end
    end
    return opts
end

function load_history(path)
    data = readdlm(path, comments=true, comment_char='#')
    if ndims(data) == 1
        data = reshape(data, 1, :)
    end
    steps = Float64.(data[:, 1])
    res = Float64.(data[:, 3])
    return steps, res
end

function main()
    opts = parse_args(ARGS)
    alphas = parse_alpha_list(string(opts["alphas"]))

    nx = Int(opts["nx"])
    ny = Int(opts["ny"])
    nz = Int(opts["nz"])
    M = Int(opts["M"])
    dt = opts["dt"]
    max_steps = Int(opts["max_steps"])
    epsilon = opts["epsilon"]
    bc_order = Symbol(opts["bc_order"])
    output_dir = string(opts["output_dir"])
    run_dir = make_run_dir(output_dir)

    dx = 1.0 / nx
    dy = 1.0 / ny
    dz = 1.0 / nz
    denom = 1 / dx^2 + 1 / dy^2 + 1 / dz^2
    dt_source = "dt"
    fo_requested = nothing
    if opts["Fo"] !== nothing
        fo_requested = opts["Fo"]
        dt = fo_requested / denom
        dt_source = "Fo"
    end
    fo = dt * denom
    dt_clipped = false
    if fo > 0.5
        dt_requested = dt
        dt = 0.5 / denom
        fo = 0.5
        dt_clipped = true
        @printf("Fo > 0.5; dt clipped: %.3e -> %.3e (Fo=0.5)\n", dt_requested, dt)
    end

    println("run config:")
    @printf("  nx=%d ny=%d nz=%d M=%d\n", nx, ny, nz, M)
    @printf("  dt=%.3e max_steps=%d epsilon=%.3e\n", dt, max_steps, epsilon)
    @printf("  bc_order=%s\n", string(bc_order))
    @printf("  alphas=%s\n", join(alphas, ","))
    @printf("  output_dir=%s\n", output_dir)
    @printf("  run_dir=%s\n", run_dir)

    run_config = Dict(
        "timestamp" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SSzzzz"),
        "script" => "scripts/compare_alpha.jl",
        "output_dir" => output_dir,
        "run_dir" => run_dir,
        "nx" => nx,
        "ny" => ny,
        "nz" => nz,
        "M" => M,
        "dt" => dt,
        "Fo" => fo,
        "Fo_requested" => fo_requested,
        "dt_source" => dt_source,
        "dt_clipped" => dt_clipped,
        "max_steps" => max_steps,
        "epsilon" => epsilon,
        "alphas" => alphas,
        "bc_order" => string(bc_order),
    )
    open(joinpath(run_dir, "run_config.toml"), "w") do io
        TOML.print(io, run_config)
    end

    results = Vector{Tuple{Float64, Float64, Int, Float64, Float64, Float64, String}}()
    for alpha in alphas
        config = SolverConfig(nx, ny, nz, M, dt, max_steps, epsilon)
        prob, _ = make_problem(config; alpha=alpha)
        t_start = time()
        sol = solve(config, prob; bc_order=bc_order, output_dir=run_dir)
        runtime = time() - t_start
        u_exact = ADPoisson.exact_solution_array(sol, prob, config)
        err_l2, err_max = ADPoisson.error_stats_precomputed(sol.u, u_exact, prob, config)
        tag = "nx$(nx)_ny$(ny)_nz$(nz)_M$(M)_steps$(sol.iter)"
        history_path = joinpath(run_dir, "history_$(tag).txt")
        alpha_tag = format_alpha_tag(alpha)
        history_path_alpha = joinpath(run_dir, "history_alpha$(alpha_tag)_$(tag).txt")
        if history_path != history_path_alpha
            mv(history_path, history_path_alpha; force=true)
        end
        push!(results, (alpha, dt, sol.iter, err_l2, err_max, runtime, history_path_alpha))
    end

    alpha_tag_list = join(map(format_alpha_tag, alphas), "-")
    summary_tag = "nx$(nx)_ny$(ny)_nz$(nz)_M$(M)_alphas$(alpha_tag_list)"
    summary_path = joinpath(run_dir, "compare_alpha_$(summary_tag).txt")
    open(summary_path, "w") do io
        println(io, "# alpha dt steps err_l2 err_max runtime_s")
        for (alpha, dt, steps, err_l2, err_max, runtime, _) in results
            @printf(io, "%.6g %.6e %d %.6e %.6e %.6f\n", alpha, dt, steps, err_l2, err_max, runtime)
        end
    end

    p = plot(xlabel="step", ylabel="res_l2 (relative)", yscale=:log10,
             title="Convergence history (res_l2)")
    for (alpha, _, _, _, _, _, history_path) in results
        steps, res = load_history(history_path)
        plot!(p, steps, res; label=@sprintf("alpha=%.3g", alpha))
    end
    plot_path = joinpath(run_dir, "history_compare_alpha_$(summary_tag).png")
    png(p, plot_path)

    run_summary = Dict(
        "summary_file" => basename(summary_path),
        "plot_file" => basename(plot_path),
        "history_files" => [basename(r[7]) for r in results],
    )
    open(joinpath(run_dir, "run_summary.toml"), "w") do io
        TOML.print(io, run_summary)
    end

    println("summary:")
    for (alpha, dt, steps, err_l2, err_max, runtime, _) in results
        @printf("  alpha=%.3g dt=%.3e steps=%d err_l2=%.3e err_max=%.3e runtime_s=%.3f\n",
                alpha, dt, steps, err_l2, err_max, runtime)
    end
    @info "outputs" summary=summary_path plot=plot_path
end

main()
