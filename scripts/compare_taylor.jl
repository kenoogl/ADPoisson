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

function warmup_solve(config::SolverConfig, prob::ProblemSpec, bc_order::Symbol,
                      output_dir::String, solver::Symbol, cg_precond::Symbol)
    warm_dir = joinpath(output_dir, "_warmup")
    isdir(warm_dir) || mkpath(warm_dir)
    warm_config = SolverConfig(config.nx, config.ny, config.nz, config.M, config.dt, 1, 0.0)
    if solver === :taylor
        solve(warm_config, prob; bc_order=bc_order, output_dir=warm_dir)
    elseif solver === :sor
        sor_solve(prob, warm_config; bc_order=bc_order, output_dir=warm_dir)
    elseif solver === :ssor
        ssor_solve(prob, warm_config; bc_order=bc_order, output_dir=warm_dir)
    elseif solver === :cg
        cg_solve(prob, warm_config; precond=cg_precond, bc_order=bc_order, output_dir=warm_dir)
    else
        error("unknown solver: $(solver)")
    end
    rm(warm_dir; recursive=true, force=true)
end

function parse_ms_list(s::String)
    parts = split(s, ",")
    Ms = Int[]
    for p in parts
        t = strip(p)
        isempty(t) && continue
        push!(Ms, parse(Int, t))
    end
    isempty(Ms) && error("Ms list is empty")
    return sort(unique(Ms))
end

function parse_args(args)
    opts = default_cli_options()
    opts["Ms"] = "2,4,6,8,10"
    opts["output_dir"] = "results"
    opts["Fo"] = nothing
    opts["solver"] = "taylor"
    opts["cg_precond"] = "none"

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
            elseif key == "Ms" || key == "ms"
                opts["Ms"] = args[i + 1]
            elseif key == "output-dir" || key == "output_dir"
                opts["output_dir"] = args[i + 1]
            elseif key == "Fo" || key == "fo"
                opts["Fo"] = parse(Float64, args[i + 1])
            elseif key == "n"
                opts["n"] = parse(Int, args[i + 1])
            elseif key == "solver"
                opts["solver"] = lowercase(args[i + 1])
            elseif key == "cg-precond" || key == "cg_precond"
                opts["cg_precond"] = lowercase(args[i + 1])
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
    Ms = parse_ms_list(string(opts["Ms"]))
    solver = Symbol(lowercase(opts["solver"]))
    cg_precond = Symbol(lowercase(opts["cg_precond"]))
    (solver === :taylor || solver === :sor || solver === :ssor || solver === :cg) ||
        error("solver must be taylor/sor/ssor/cg")
    (cg_precond === :ssor || cg_precond === :none) ||
        error("cg-precond must be ssor/none")

    if opts["n"] !== nothing
        nx = Int(opts["n"])
        ny = Int(opts["n"])
        nz = Int(opts["n"])
    else
        nx = Int(opts["nx"])
        ny = Int(opts["ny"])
        nz = Int(opts["nz"])
    end
    dt = opts["dt"]
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

    if solver !== :taylor
        @printf("solver=%s: M list is ignored; --M is optional and ignored for iterative solvers\n", string(solver))
        Ms = [Int(opts["M"])]
    end

    println("run config:")
    @printf("  nx=%d ny=%d nz=%d\n", nx, ny, nz)
    @printf("  dt=%.3e max_steps=%d epsilon=%.3e\n", dt, max_steps, epsilon)
    @printf("  alpha=%.6f bc_order=%s\n", alpha, string(bc_order))
    @printf("  solver=%s\n", string(solver))
    if solver === :cg
        @printf("  cg_precond=%s\n", string(cg_precond))
    end
    @printf("  Ms=%s\n", join(Ms, ","))
    @printf("  output_dir=%s\n", output_dir)
    @printf("  run_dir=%s\n", run_dir)

    run_config = Dict(
        "timestamp" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SSzzzz"),
        "script" => "scripts/compare_taylor.jl",
        "output_dir" => output_dir,
        "run_dir" => run_dir,
        "nx" => nx,
        "ny" => ny,
        "nz" => nz,
        "max_steps" => max_steps,
        "epsilon" => epsilon,
        "alpha" => alpha,
        "bc_order" => string(bc_order),
        "solver" => string(solver),
        "cg_precond" => string(cg_precond),
        "warmup" => true,
    )
    if solver === :taylor
        run_config["Ms"] = Ms
        run_config["dt"] = dt
        run_config["Fo"] = fo
        run_config["dt_source"] = dt_source
        run_config["dt_clipped"] = dt_clipped
    end
    if fo_requested !== nothing
        run_config["Fo_requested"] = fo_requested
    end
    open(joinpath(run_dir, "run_config.toml"), "w") do io
        TOML.print(io, run_config)
    end

    warm_config = SolverConfig(nx, ny, nz, Ms[1], dt, max_steps, epsilon)
    warm_prob, _ = make_problem(warm_config; alpha=alpha)
    warmup_solve(warm_config, warm_prob, bc_order, run_dir, solver, cg_precond)

    results = Vector{Tuple{Int, Int, Float64, Float64, Float64, String}}()
    for M in Ms
        config = SolverConfig(nx, ny, nz, M, dt, max_steps, epsilon)
        prob, _ = make_problem(config; alpha=alpha)
        sol, runtime = if solver === :taylor
            solve_with_runtime(config, prob; bc_order=bc_order, output_dir=run_dir)
        elseif solver === :sor
            sor_solve_with_runtime(prob, config; bc_order=bc_order, output_dir=run_dir)
        elseif solver === :ssor
            ssor_solve_with_runtime(prob, config; bc_order=bc_order, output_dir=run_dir)
        else
            cg_solve_with_runtime(prob, config; precond=cg_precond, bc_order=bc_order, output_dir=run_dir)
        end
        u_exact = ADPoisson.exact_solution_array(sol, prob, config)
        err_l2, err_max = ADPoisson.error_stats_precomputed(sol.u, u_exact, prob, config)
        tag = "nx$(nx)_ny$(ny)_nz$(nz)_M$(M)_steps$(sol.iter)"
        history_path = if solver === :taylor
            joinpath(run_dir, "history_$(tag).txt")
        elseif solver === :sor
            joinpath(run_dir, "history_sor_nx$(nx)_ny$(ny)_nz$(nz)_steps$(sol.iter).txt")
        elseif solver === :ssor
            joinpath(run_dir, "history_ssor_nx$(nx)_ny$(ny)_nz$(nz)_steps$(sol.iter).txt")
        else
            joinpath(run_dir, "history_cg_nx$(nx)_ny$(ny)_nz$(nz)_steps$(sol.iter).txt")
        end
        push!(results, (M, sol.iter, err_l2, err_max, runtime, history_path))
    end

    summary_tag = "nx$(nx)_ny$(ny)_nz$(nz)_Ms$(join(Ms, "-"))_solver$(solver)"
    summary_path = joinpath(run_dir, "compare_M_$(summary_tag).txt")
    open(summary_path, "w") do io
        println(io, "# M steps err_l2 err_max runtime_s")
        for (M, steps, err_l2, err_max, runtime, _) in results
            @printf(io, "%d %d %.6e %.6e %.6f\n", M, steps, err_l2, err_max, runtime)
        end
    end

    p = plot(xlabel="step", ylabel="res_l2 (relative)", yscale=:log10,
             title="Convergence history (res_l2)")
    for (M, _, _, _, _, history_path) in results
        steps, res = load_history(history_path)
        plot!(p, steps, res; label="M=$(M)")
    end
    plot_path = joinpath(run_dir, "history_compare_$(summary_tag).png")
    png(p, plot_path)

    run_summary = Dict(
        "summary_file" => basename(summary_path),
        "plot_file" => basename(plot_path),
        "history_files" => [basename(r[6]) for r in results],
    )
    open(joinpath(run_dir, "run_summary.toml"), "w") do io
        TOML.print(io, run_summary)
    end

    println("summary:")
    for (M, steps, err_l2, err_max, runtime, _) in results
        @printf("  M=%d steps=%d err_l2=%.3e err_max=%.3e runtime_s=%.3f\n",
                M, steps, err_l2, err_max, runtime)
    end
    @info "outputs" summary=summary_path plot=plot_path
end

main()
