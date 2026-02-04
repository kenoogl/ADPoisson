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

function warmup_solve(config::SolverConfig, prob::ProblemSpec, bc_order::Symbol, output_dir::String)
    warm_dir = joinpath(output_dir, "_warmup")
    isdir(warm_dir) || mkpath(warm_dir)
    warm_config = SolverConfig(config.nx, config.ny, config.nz, config.M, config.dt, 1, 0.0)
    solve(warm_config, prob; bc_order=bc_order, output_dir=warm_dir)
    rm(warm_dir; recursive=true, force=true)
end

function parse_fo_list(s::String)
    parts = split(s, ",")
    Fos = Float64[]
    for p in parts
        t = strip(p)
        isempty(t) && continue
        push!(Fos, parse(Float64, t))
    end
    isempty(Fos) && error("Fo list is empty")
    return sort(unique(Fos))
end

function format_fo_tag(fo::Real)
    s = @sprintf("%.3g", fo)
    s = replace(s, "." => "p", "-" => "m", "+" => "")
    return s
end

function parse_args(args)
    opts = default_cli_options()
    opts["Fos"] = "0.1,0.2,0.3,0.4,0.5"
    opts["output_dir"] = "results"

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
            elseif key == "Fos" || key == "fos"
                opts["Fos"] = args[i + 1]
            elseif key == "output-dir" || key == "output_dir"
                opts["output_dir"] = args[i + 1]
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
    Fos = parse_fo_list(string(opts["Fos"]))

    nx = Int(opts["nx"])
    ny = Int(opts["ny"])
    nz = Int(opts["nz"])
    M = Int(opts["M"])
    max_steps = Int(opts["max_steps"])
    epsilon = opts["epsilon"]
    alpha = opts["alpha"]
    bc_order = Symbol(opts["bc_order"])
    output_dir = string(opts["output_dir"])
    run_dir = make_run_dir(output_dir)

    dx = 1.0 / nx
    dy = 1.0 / ny
    dz = 1.0 / nz
    denom = 1 / dx^2 + 1 / dy^2 + 1 / dz^2

    println("run config:")
    @printf("  nx=%d ny=%d nz=%d M=%d\n", nx, ny, nz, M)
    @printf("  max_steps=%d epsilon=%.3e\n", max_steps, epsilon)
    @printf("  alpha=%.6f bc_order=%s\n", alpha, string(bc_order))
    @printf("  Fos=%s\n", join(Fos, ","))
    @printf("  output_dir=%s\n", output_dir)
    @printf("  run_dir=%s\n", run_dir)

    run_config = Dict(
        "timestamp" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SSzzzz"),
        "script" => "scripts/compare_fo.jl",
        "output_dir" => output_dir,
        "run_dir" => run_dir,
        "nx" => nx,
        "ny" => ny,
        "nz" => nz,
        "M" => M,
        "Fos" => Fos,
        "max_steps" => max_steps,
        "epsilon" => epsilon,
        "alpha" => alpha,
        "bc_order" => string(bc_order),
        "warmup" => true,
    )
    open(joinpath(run_dir, "run_config.toml"), "w") do io
        TOML.print(io, run_config)
    end

    warm_fo = Fos[1]
    if warm_fo > 0.5
        warm_fo = 0.5
    end
    warm_dt = warm_fo / denom
    warm_config = SolverConfig(nx, ny, nz, M, warm_dt, 1, 0.0)
    warm_prob, _ = make_problem(warm_config; alpha=alpha)
    warmup_solve(warm_config, warm_prob, bc_order, run_dir)

    results = Vector{Tuple{Float64, Float64, Int, Float64, Float64, Float64, String}}()
    for fo in Fos
        dt = fo / denom
        if fo > 0.5
            @printf("Fo > 0.5; running with Fo=%.3g (dt=%.3e)\n", fo, dt)
        end
        config = SolverConfig(nx, ny, nz, M, dt, max_steps, epsilon)
        prob, _ = make_problem(config; alpha=alpha)
        sol, runtime = solve_with_runtime(config, prob; bc_order=bc_order, output_dir=run_dir)
        u_exact = ADPoisson.exact_solution_array(sol, prob, config)
        err_l2, err_max = ADPoisson.error_stats_precomputed(sol.u, u_exact, prob, config)
        tag = "nx$(nx)_ny$(ny)_nz$(nz)_M$(M)_steps$(sol.iter)"
        history_path = joinpath(run_dir, "history_$(tag).txt")
        fo_tag = format_fo_tag(fo)
        history_path_fo = joinpath(run_dir, "history_Fo$(fo_tag)_$(tag).txt")
        if history_path != history_path_fo
            mv(history_path, history_path_fo; force=true)
        end
        push!(results, (fo, dt, sol.iter, err_l2, err_max, runtime, history_path_fo))
    end

    fo_tag_list = join(map(format_fo_tag, Fos), "-")
    summary_tag = "nx$(nx)_ny$(ny)_nz$(nz)_M$(M)_Fos$(fo_tag_list)"
    summary_path = joinpath(run_dir, "compare_Fo_$(summary_tag).txt")
    open(summary_path, "w") do io
        println(io, "# Fo dt steps err_l2 err_max runtime_s")
        for (fo, dt, steps, err_l2, err_max, runtime, _) in results
            @printf(io, "%.6g %.6e %d %.6e %.6e %.6f\n", fo, dt, steps, err_l2, err_max, runtime)
        end
    end

    p = plot(xlabel="step", ylabel="res_l2 (relative)", yscale=:log10,
             title="Convergence history (res_l2)")
    for (fo, _, _, _, _, _, history_path) in results
        steps, res = load_history(history_path)
        plot!(p, steps, res; label=@sprintf("Fo=%.3g", fo))
    end
    plot_path = joinpath(run_dir, "history_compare_Fo_$(summary_tag).png")
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
    for (fo, dt, steps, err_l2, err_max, runtime, _) in results
        @printf("  Fo=%.3g dt=%.3e steps=%d err_l2=%.3e err_max=%.3e runtime_s=%.3f\n",
                fo, dt, steps, err_l2, err_max, runtime)
    end
    @info "outputs" summary=summary_path plot=plot_path
end

main()
