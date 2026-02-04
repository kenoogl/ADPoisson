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
    nx = Int(opts["nx"])
    ny = Int(opts["ny"])
    nz = Int(opts["nz"])
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
    output_dir = string(opts["output_dir"])
    run_dir = make_run_dir(output_dir)
    println("run config:")
    @printf("  nx=%d ny=%d nz=%d M=%d\n", config.nx, config.ny, config.nz, config.M)
    @printf("  dt=%.3e max_steps=%d epsilon=%.3e\n", config.dt, config.max_steps, config.epsilon)
    @printf("  alpha=%.6f bc_order=%s\n", prob.alpha, string(bc_order))
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
        "M" => config.M,
        "dt" => config.dt,
        "Fo" => fo,
        "dt_source" => dt_source,
        "max_steps" => config.max_steps,
        "epsilon" => config.epsilon,
        "alpha" => prob.alpha,
        "bc_order" => string(bc_order),
    )
    open(joinpath(run_dir, "run_config.toml"), "w") do io
        TOML.print(io, run_config)
    end
    t_start = time()
    sol = solve(config, prob; bc_order=bc_order, output_dir=run_dir)
    runtime = time() - t_start
    plot_slice(sol, prob, config; output_dir=run_dir)
    u_exact = ADPoisson.exact_solution_array(sol, prob, config)
    err_l2, err_max = ADPoisson.error_stats_precomputed(sol.u, u_exact, prob, config)
    tag = "nx$(config.nx)_ny$(config.ny)_nz$(config.nz)_M$(config.M)_steps$(sol.iter)"
    run_summary = Dict(
        "steps" => sol.iter,
        "err_l2" => err_l2,
        "err_max" => err_max,
        "runtime_s" => runtime,
        "history_file" => "history_$(tag).txt",
        "error_plot" => "error_$(tag).png",
        "exact_plot" => "exact_nx$(config.nx)_ny$(config.ny)_nz$(config.nz).png",
    )
    open(joinpath(run_dir, "run_summary.toml"), "w") do io
        TOML.print(io, run_summary)
    end
    @info "done" t=sol.t iter=sol.iter
end

main()
