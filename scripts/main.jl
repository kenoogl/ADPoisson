#!/usr/bin/env julia

using ADPoisson
using Printf
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

function parse_args(args)
    opts = default_cli_options()
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

function main()
    opts = parse_args(ARGS)
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
    if opts["Fo"] !== nothing
        dx = 1.0 / nx
        dy = 1.0 / ny
        dz = 1.0 / nz
        denom = 1 / dx^2 + 1 / dy^2 + 1 / dz^2
        dt = opts["Fo"] / denom
    end
    dx = 1.0 / nx
    dy = 1.0 / ny
    dz = 1.0 / nz
    denom = 1 / dx^2 + 1 / dy^2 + 1 / dz^2
    fo = dt * denom
    dt_source = "dt"
    if opts["Fo"] !== nothing
        fo = opts["Fo"]
        dt = fo / denom
        dt_source = "Fo"
    end
    config = SolverConfig(nx, ny, nz,
                          Int(opts["M"]), dt, Int(opts["max_steps"]), opts["epsilon"])
    prob, _ = make_problem(config; alpha=opts["alpha"])
    bc_order = Symbol(opts["bc_order"])
    solver = Symbol(lowercase(opts["solver"]))
    cg_precond = Symbol(lowercase(opts["cg_precond"]))
    (solver === :taylor || solver === :sor || solver === :ssor || solver === :cg) ||
        error("solver must be taylor/sor/cg")
    (cg_precond === :ssor || cg_precond === :none) ||
        error("cg-precond must be ssor/none")
    output_dir = string(opts["output_dir"])
    run_dir = make_run_dir(output_dir)
    println("run config:")
    @printf("  nx=%d ny=%d nz=%d M=%d\n", config.nx, config.ny, config.nz, config.M)
    @printf("  dt=%.3e max_steps=%d epsilon=%.3e\n", config.dt, config.max_steps, config.epsilon)
    @printf("  alpha=%.6f bc_order=%s\n", prob.alpha, string(bc_order))
    @printf("  solver=%s\n", string(solver))
    if solver === :cg
        @printf("  cg_precond=%s\n", string(cg_precond))
    end
    @printf("  output_dir=%s\n", output_dir)
    @printf("  run_dir=%s\n", run_dir)
    run_config = Dict(
        "timestamp" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SSzzzz"),
        "script" => "scripts/main.jl",
        "output_dir" => output_dir,
        "run_dir" => run_dir,
        "nx" => config.nx,
        "ny" => config.ny,
        "nz" => config.nz,
        "max_steps" => config.max_steps,
        "epsilon" => config.epsilon,
        "alpha" => prob.alpha,
        "bc_order" => string(bc_order),
        "solver" => string(solver),
        "warmup" => true,
    )
    if solver === :taylor
        run_config["M"] = config.M
        run_config["dt"] = config.dt
        run_config["Fo"] = fo
        run_config["dt_source"] = dt_source
    elseif solver === :cg
        run_config["cg_precond"] = string(cg_precond)
    end
    open(joinpath(run_dir, "run_config.toml"), "w") do io
        TOML.print(io, run_config)
    end
    warmup_solve(config, prob, bc_order, run_dir, solver, cg_precond)
    sol, runtime = if solver === :taylor
        solve_with_runtime(config, prob; bc_order=bc_order, output_dir=run_dir)
    elseif solver === :sor
        sor_solve_with_runtime(prob, config; bc_order=bc_order, output_dir=run_dir)
    elseif solver === :ssor
        ssor_solve_with_runtime(prob, config; bc_order=bc_order, output_dir=run_dir)
    else
        cg_solve_with_runtime(prob, config; precond=cg_precond, bc_order=bc_order, output_dir=run_dir)
    end
    plot_slice(sol, prob, config; output_dir=run_dir)
    u_exact = ADPoisson.exact_solution_array(sol, prob, config)
    err_l2, err_max = ADPoisson.error_stats_precomputed(sol.u, u_exact, prob, config)
    tag = "nx$(config.nx)_ny$(config.ny)_nz$(config.nz)_M$(config.M)_steps$(sol.iter)"
    history_file = if solver === :taylor
        "history_$(tag).txt"
    elseif solver === :sor
        "history_sor_nx$(config.nx)_ny$(config.ny)_nz$(config.nz)_steps$(sol.iter).txt"
    elseif solver === :ssor
        "history_ssor_nx$(config.nx)_ny$(config.ny)_nz$(config.nz)_steps$(sol.iter).txt"
    else
        "history_cg_nx$(config.nx)_ny$(config.ny)_nz$(config.nz)_steps$(sol.iter).txt"
    end
    run_summary = Dict(
        "steps" => sol.iter,
        "err_l2" => err_l2,
        "err_max" => err_max,
        "runtime_s" => runtime,
        "history_file" => history_file,
        "error_plot" => "error_$(tag).png",
        "exact_plot" => "exact_nx$(config.nx)_ny$(config.ny)_nz$(config.nz).png",
    )
    open(joinpath(run_dir, "run_summary.toml"), "w") do io
        TOML.print(io, run_summary)
    end
    @info "done" t=sol.t iter=sol.iter
end

main()
