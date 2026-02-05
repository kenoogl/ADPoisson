#!/usr/bin/env julia

using ADPoisson
using DelimitedFiles
using Printf
using Plots
using Dates
using TOML

function last_res_l2(path::AbstractString)
    last = nothing
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            if isempty(s) || startswith(s, "#")
                continue
            end
            last = s
        end
    end
    if last === nothing
        return nothing
    end
    parts = split(last)
    if length(parts) < 3
        return nothing
    end
    try
        return parse(Float64, parts[3])
    catch
        return nothing
    end
end

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

function parse_solvers_list(s::String)
    parts = split(s, ",")
    solvers = Symbol[]
    for p in parts
        t = lowercase(strip(p))
        isempty(t) && continue
        push!(solvers, Symbol(t))
    end
    isempty(solvers) && error("solvers list is empty")
    return solvers
end

function parse_args(args)
    opts = default_cli_options()
    opts["output_dir"] = "results"
    opts["Fo"] = nothing
    opts["solver_list"] = "taylor,sor,ssor,cg"
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
            elseif key == "output-dir" || key == "output_dir"
                opts["output_dir"] = args[i + 1]
            elseif key == "Fo" || key == "fo"
                opts["Fo"] = parse(Float64, args[i + 1])
            elseif key == "n"
                opts["n"] = parse(Int, args[i + 1])
            elseif key == "solvers" || key == "solver-list" || key == "solver_list"
                opts["solver_list"] = args[i + 1]
            elseif key == "cg-precond" || key == "cg_precond"
                opts["cg_precond"] = lowercase(args[i + 1])
            elseif key == "max-steps" || key == "max_steps"
                opts["max_steps"] = parse(Int, args[i + 1])
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
    solvers = parse_solvers_list(string(opts["solver_list"]))
    cg_precond = Symbol(lowercase(opts["cg_precond"]))
    (cg_precond === :ssor || cg_precond === :none) || error("cg-precond must be ssor/none")

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

    println("run config:")
    @printf("  nx=%d ny=%d nz=%d\n", nx, ny, nz)
    @printf("  dt=%.3e max_steps=%d epsilon=%.3e\n", dt, max_steps, epsilon)
    @printf("  alpha=%.6f bc_order=%s\n", alpha, string(bc_order))
    @printf("  solvers=%s\n", join(string.(solvers), ","))
    if :cg in solvers
        @printf("  cg_precond=%s\n", string(cg_precond))
    end
    @printf("  output_dir=%s\n", output_dir)
    @printf("  run_dir=%s\n", run_dir)

    run_config = Dict(
        "timestamp" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SSzzzz"),
        "script" => "scripts/compare_solvers.jl",
        "output_dir" => output_dir,
        "run_dir" => run_dir,
        "nx" => nx,
        "ny" => ny,
        "nz" => nz,
        "M" => Int(opts["M"]),
        "max_steps" => max_steps,
        "epsilon" => epsilon,
        "alpha" => alpha,
        "bc_order" => string(bc_order),
        "solvers" => string.(solvers),
        "cg_precond" => string(cg_precond),
        "warmup" => true,
    )
    run_config["dt"] = dt
    run_config["Fo"] = fo
    run_config["dt_source"] = dt_source
    run_config["dt_clipped"] = dt_clipped
    if fo_requested !== nothing
        run_config["Fo_requested"] = fo_requested
    end
    open(joinpath(run_dir, "run_config.toml"), "w") do io
        TOML.print(io, run_config)
    end

    warm_config = SolverConfig(nx, ny, nz, Int(opts["M"]), dt, max_steps, epsilon)
    warm_prob, _ = make_problem(warm_config; alpha=alpha)
    for solver in solvers
        bc_order_eff = bc_order
        if solver !== :taylor && bc_order_eff !== :spec
            @warn "bc-order is only valid for taylor; forcing to spec for iterative solvers" solver=solver bc_order=bc_order_eff
            bc_order_eff = :spec
        end
        warmup_solve(warm_config, warm_prob, bc_order_eff, run_dir, solver, cg_precond)
    end

    results = Vector{Tuple{Symbol, Int, Float64, Float64, Float64, Float64, String}}()
    for solver in solvers
        bc_order_eff = bc_order
        if solver !== :taylor && bc_order_eff !== :spec
            @warn "bc-order is only valid for taylor; forcing to spec for iterative solvers" solver=solver bc_order=bc_order_eff
            bc_order_eff = :spec
        end
        config = SolverConfig(nx, ny, nz, Int(opts["M"]), dt, max_steps, epsilon)
        prob, _ = make_problem(config; alpha=alpha)
        sol, runtime = if solver === :taylor
            solve_with_runtime(config, prob; bc_order=bc_order_eff, output_dir=run_dir)
        elseif solver === :sor
            sor_solve_with_runtime(prob, config; bc_order=bc_order_eff, output_dir=run_dir)
        elseif solver === :ssor
            ssor_solve_with_runtime(prob, config; bc_order=bc_order_eff, output_dir=run_dir)
        elseif solver === :cg
            cg_solve_with_runtime(prob, config; precond=cg_precond, bc_order=bc_order_eff, output_dir=run_dir)
        else
            error("unknown solver: $(solver)")
        end
        u_exact = ADPoisson.exact_solution_array(sol, prob, config)
        err_l2, err_max = ADPoisson.error_stats_precomputed(sol.u, u_exact, prob, config)
        history_path = if solver === :taylor
            tag = "nx$(nx)_ny$(ny)_nz$(nz)_M$(config.M)_steps$(sol.iter)"
            joinpath(run_dir, "history_$(tag).txt")
        elseif solver === :sor
            joinpath(run_dir, "history_sor_nx$(nx)_ny$(ny)_nz$(nz)_steps$(sol.iter).txt")
        elseif solver === :ssor
            joinpath(run_dir, "history_ssor_nx$(nx)_ny$(ny)_nz$(nz)_steps$(sol.iter).txt")
        elseif solver === :cg
            joinpath(run_dir, "history_cg_nx$(nx)_ny$(ny)_nz$(nz)_steps$(sol.iter).txt")
        else
            error("unknown solver: $(solver)")
        end
        res_l2 = last_res_l2(history_path)
        push!(results, (solver, sol.iter, err_l2, err_max, res_l2, runtime, history_path))
    end

    summary_tag = "nx$(nx)_ny$(ny)_nz$(nz)_solver$(join(string.(solvers), "-"))"
    summary_path = joinpath(run_dir, "compare_solvers_$(summary_tag).txt")
    open(summary_path, "w") do io
        println(io, "# solver steps err_l2 err_max res_l2 runtime_s")
        for (solver, steps, err_l2, err_max, res_l2, runtime, _) in results
            @printf(io, "%s %d %.6e %.6e %.6e %.6f\n",
                    string(solver), steps, err_l2, err_max, res_l2, runtime)
        end
    end

    p = plot(xlabel="step", ylabel="res_l2 (relative)", yscale=:log10,
             yticks=10.0 .^ (-16:1:0),
             title="Convergence history (res_l2)")
    for (solver, _, _, _, _, _, history_path) in results
        steps, res = load_history(history_path)
        plot!(p, steps, res; label=string(solver))
    end
    plot_path = joinpath(run_dir, "history_compare_solvers_$(summary_tag).png")
    png(p, plot_path)

    run_summary = Dict(
        "summary_file" => basename(summary_path),
        "plot_file" => basename(plot_path),
        "history_files" => [basename(r[7]) for r in results],
        "res_l2_list" => [r[5] for r in results],
    )
    open(joinpath(run_dir, "run_summary.toml"), "w") do io
        TOML.print(io, run_summary)
    end

    println("summary:")
    for (solver, steps, err_l2, err_max, res_l2, runtime, _) in results
        @printf("  solver=%s steps=%d err_l2=%.3e err_max=%.3e res_l2=%.3e runtime_s=%.3f\n",
                string(solver), steps, err_l2, err_max, res_l2, runtime)
    end
    @info "outputs" summary=summary_path plot=plot_path
end

main()
