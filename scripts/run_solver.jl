#!/usr/bin/env julia

using ADPoisson
using Printf
using Dates
using TOML
using YAML

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
    return output_dir
end

function warmup_solve(config::SolverConfig, prob::ProblemSpec, bc_order::Symbol,
                      lap_order::Symbol, output_dir::String, solver::Symbol, cg_precond::Symbol)
    warm_dir = joinpath(output_dir, "_warmup")
    isdir(warm_dir) || mkpath(warm_dir)
    warm_config = SolverConfig(config.nx, config.ny, config.nz, config.M, config.dt, 1, 0.0)
    if solver === :taylor
        solve(warm_config, prob; bc_order=bc_order, lap_order=lap_order, output_dir=warm_dir)
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

function parse_int_list(value::AbstractString)
    parts = split(value, ',')
    vals = Int[]
    for part in parts
        t = strip(part)
        isempty(t) && continue
        push!(vals, parse(Int, t))
    end
    return isempty(vals) ? nothing : vals
end

function parse_float_list(value::AbstractString)
    parts = split(value, ',')
    vals = Float64[]
    for part in parts
        t = strip(part)
        isempty(t) && continue
        push!(vals, parse(Float64, t))
    end
    return isempty(vals) ? nothing : vals
end

function toml_scalar(v)
    if v === nothing
        return "\"\""
    elseif v isa AbstractString
        s = replace(v, "\\" => "\\\\", "\"" => "\\\"")
        return "\"$(s)\""
    elseif v isa Bool
        return v ? "true" : "false"
    elseif v isa Integer
        return string(v)
    elseif v isa AbstractFloat
        if isnan(v)
            return "nan"
        elseif isinf(v)
            return signbit(v) ? "-inf" : "inf"
        end
        return repr(v)
    elseif v isa Tuple || v isa AbstractVector
        return "[" * join((toml_scalar(x) for x in v), ", ") * "]"
    end
    return toml_scalar(string(v))
end

function write_run_summary(path::AbstractString; steps, runtime_s, res_l2, err_l2, err_max,
                           history_file, error_plot, exact_plot,
                           mg_levels_used=nothing, mg_coarsest_grid=nothing)
    open(path, "w") do io
        println(io, "# Run")
        println(io, "steps = ", toml_scalar(steps))
        println(io, "runtime_s = ", toml_scalar(runtime_s))
        println(io)
        println(io, "# Error & residual")
        println(io, "err_l2 = ", toml_scalar(err_l2))
        println(io, "err_max = ", toml_scalar(err_max))
        println(io, "res_l2 = ", toml_scalar(res_l2))
        println(io)
        println(io, "# Artifacts")
        println(io, "history_file = ", toml_scalar(history_file))
        println(io, "error_plot = ", toml_scalar(error_plot))
        println(io, "exact_plot = ", toml_scalar(exact_plot))
        if mg_levels_used !== nothing && mg_coarsest_grid !== nothing
            println(io)
            println(io, "# Multigrid")
            println(io, "mg_levels_used = ", toml_scalar(mg_levels_used))
            println(io, "mg_coarsest_grid = ", toml_scalar(mg_coarsest_grid))
        end
    end
end

function parse_bool(value::AbstractString)
    v = lowercase(strip(value))
    if v == "1" || v == "true" || v == "yes" || v == "on"
        return true
    elseif v == "0" || v == "false" || v == "no" || v == "off"
        return false
    end
    error("invalid boolean: $(value)")
end

function dict_get(dict::AbstractDict, key)
    if haskey(dict, key)
        return dict[key]
    end
    if key isa String && haskey(dict, Symbol(key))
        return dict[Symbol(key)]
    end
    if key isa Symbol && haskey(dict, String(key))
        return dict[String(key)]
    end
    return nothing
end

function cfg_get(cfg, keys...)
    v = cfg
    for key in keys
        if !(v isa AbstractDict)
            return nothing
        end
        v = dict_get(v, key)
        if v === nothing
            return nothing
        end
    end
    return v
end

function to_int(x)
    return x isa Integer ? Int(x) : parse(Int, string(x))
end

function to_float(x)
    return x isa AbstractFloat ? Float64(x) : parse(Float64, string(x))
end

function to_string(x)
    return x isa AbstractString ? String(x) : string(x)
end

function to_int_list(x)
    if x isa AbstractVector
        return [to_int(v) for v in x]
    end
    return parse_int_list(to_string(x))
end

function to_float_list(x)
    if x isa AbstractVector
        return [to_float(v) for v in x]
    end
    return parse_float_list(to_string(x))
end

function apply_run_config!(opts, run, seen)
    for (raw_key, raw_val) in run
        key = to_string(raw_key)
        val = raw_val
        if key == "n" || key == "nx" || key == "ny" || key == "nz" || key == "M" || key == "max_steps"
            opts[key] = to_int(val)
        elseif key == "Fo" || key == "dt" || key == "alpha" || key == "epsilon" ||
               key == "mg_dt_scale" || key == "mg_corr_dt_scale"
            opts[key] = to_float(val)
        elseif key == "solver" || key == "cg_precond" || key == "bc_order" || key == "lap_order"
            opts[key] = lowercase(to_string(val))
        elseif key == "mg_level_Ms"
            opts[key] = to_int_list(val)
        elseif key == "mg_level_dt_scales"
            opts[key] = to_float_list(val)
        elseif key == "mg_interval" || key == "mg_M" || key == "mg_nu1" || key == "mg_nu2" ||
               key == "mg_corr_M" || key == "mg_corr_steps" || key == "mg_corr_nu1" || key == "mg_corr_nu2"
            opts[key] = to_int(val)
        elseif key == "output_dir"
            opts[key] = to_string(val)
        elseif key == "debug_residual" || key == "debug_vcycle"
            opts[key] = Bool(val)
        end
        push!(seen, key)
    end
end

function apply_config!(opts, cfg)
    seen = Set{String}()
    run = cfg_get(cfg, "run")
    if run !== nothing
        apply_run_config!(opts, run, seen)
    end

    results_dir = cfg_get(cfg, "output", "results_dir")
    if results_dir !== nothing && !("output_dir" in seen)
        opts["output_dir"] = to_string(results_dir)
    end

    grid = cfg_get(cfg, "discretization", "grid")
    if grid !== nothing
        for (k, optkey) in (("nx", "nx"), ("ny", "ny"), ("nz", "nz"))
            v = dict_get(grid, k)
            if v !== nothing && !(optkey in seen)
                opts[optkey] = to_int(v)
            end
        end
    end

    tol = cfg_get(cfg, "convergence", "tolerance")
    if tol !== nothing && !("epsilon" in seen)
        opts["epsilon"] = to_float(tol)
    end
    maxit = cfg_get(cfg, "convergence", "max_iterations")
    if maxit !== nothing && !("max_steps" in seen)
        opts["max_steps"] = to_int(maxit)
    end

    solver_cfg = cfg_get(cfg, "solver")
    if solver_cfg !== nothing
        method = dict_get(solver_cfg, "method")
        stype = dict_get(solver_cfg, "type")
        if method !== nothing && !("solver" in seen)
            m = lowercase(to_string(method))
            t = stype === nothing ? "" : lowercase(to_string(stype))
            if m == "taylor"
                if t == "uniform"
                    opts["solver"] = "mg-uniform-taylor"
                elseif t == "hierarchical"
                    opts["solver"] = "mg-hierarchical-taylor"
                elseif t == "correction-taylor" || t == "correction"
                    opts["solver"] = "mg-correction-taylor"
                else
                    opts["solver"] = "taylor"
                end
            else
                opts["solver"] = m
            end
        end
        if dict_get(solver_cfg, "Fo") !== nothing && !("Fo" in seen)
            opts["Fo"] = to_float(dict_get(solver_cfg, "Fo"))
        end
        if dict_get(solver_cfg, "dt") !== nothing && !("dt" in seen)
            opts["dt"] = to_float(dict_get(solver_cfg, "dt"))
        end
        if dict_get(solver_cfg, "M") !== nothing && !("M" in seen)
            opts["M"] = to_int(dict_get(solver_cfg, "M"))
        end
        smoother = cfg_get(solver_cfg, "smoother", "type")
        if smoother !== nothing && !("cg_precond" in seen)
            opts["cg_precond"] = lowercase(to_string(smoother))
        end
    end
end

function mg_levels_and_coarsest(config::SolverConfig)
    levels = ADPoisson.mg_max_levels(config; min_n=4)
    cx = config.nx
    cy = config.ny
    cz = config.nz
    for _ in 2:levels
        cx ÷= 2
        cy ÷= 2
        cz ÷= 2
    end
    return levels, Dict("nx" => cx, "ny" => cy, "nz" => cz)
end

function parse_args(args)
    opts = default_cli_options()
    opts["config_path"] = nothing
    opts["output_dir"] = "results"
    opts["Fo"] = nothing
    opts["mg_interval"] = 0
    opts["mg_dt_scale"] = 2.0
    opts["mg_M"] = 4
    opts["mg_nu1"] = 1
    opts["mg_nu2"] = 1
    opts["mg_level_Ms"] = nothing
    opts["mg_level_dt_scales"] = nothing
    opts["mg_corr_M"] = 2
    opts["mg_corr_dt_scale"] = 1.0
    opts["mg_corr_steps"] = 1
    opts["mg_corr_nu1"] = nothing
    opts["mg_corr_nu2"] = nothing
    opts["debug_residual"] = false
    opts["debug_vcycle"] = false

    config_path = nothing
    i = 1
    while i <= length(args)
        if args[i] == "--config"
            if i == length(args)
                error("Missing value for --config")
            end
            config_path = args[i + 1]
            i += 2
        elseif startswith(args[i], "--config=")
            config_path = split(args[i], "=", limit=2)[2]
            i += 1
        else
            i += 1
        end
    end
    if config_path !== nothing
        cfg = YAML.load_file(config_path)
        apply_config!(opts, cfg)
        opts["config_path"] = config_path
    end

    i = 1
    while i <= length(args)
        if startswith(args[i], "--")
            key = args[i][3:end]
            if key == "config"
                i += 2
                continue
            elseif startswith(key, "config=")
                i += 1
                continue
            end
            if i == length(args)
                error("Missing value for --$key")
            end
            if key == "bc-order" || key == "bc_order"
                opts["bc_order"] = args[i + 1]
            elseif key == "lap-order" || key == "lap_order"
                opts["lap_order"] = lowercase(args[i + 1])
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
            elseif key == "mg-interval" || key == "mg_interval"
                opts["mg_interval"] = parse(Int, args[i + 1])
            elseif key == "mg-dt-scale" || key == "mg_dt_scale"
                opts["mg_dt_scale"] = parse(Float64, args[i + 1])
            elseif key == "mg-M" || key == "mg_M"
                opts["mg_M"] = parse(Int, args[i + 1])
            elseif key == "mg-nu1" || key == "mg_nu1"
                opts["mg_nu1"] = parse(Int, args[i + 1])
            elseif key == "mg-nu2" || key == "mg_nu2"
                opts["mg_nu2"] = parse(Int, args[i + 1])
            elseif key == "mg-vcycle" || key == "mg_vcycle"
                error("--mg-vcycle is removed. Use --solver mg-uniform-taylor or mg-hierarchical-taylor.")
            elseif key == "mg-vcycle-mode" || key == "mg_vcycle_mode"
                error("--mg-vcycle-mode is removed. Use --solver mg-hierarchical-taylor.")
            elseif key == "mg-correction" || key == "mg_correction"
                error("--mg-correction is removed. Use --solver mg-correction-taylor.")
            elseif key == "mg-corr-M" || key == "mg_corr_M"
                opts["mg_corr_M"] = parse(Int, args[i + 1])
            elseif key == "mg-corr-dt-scale" || key == "mg_corr_dt_scale"
                opts["mg_corr_dt_scale"] = parse(Float64, args[i + 1])
            elseif key == "mg-corr-steps" || key == "mg_corr_steps"
                opts["mg_corr_steps"] = parse(Int, args[i + 1])
            elseif key == "mg-corr-nu1" || key == "mg_corr_nu1"
                opts["mg_corr_nu1"] = parse(Int, args[i + 1])
            elseif key == "mg-corr-nu2" || key == "mg_corr_nu2"
                opts["mg_corr_nu2"] = parse(Int, args[i + 1])
            elseif key == "mg-level-Ms" || key == "mg_level_Ms" || key == "mg-level-ms" || key == "mg_level_ms"
                opts["mg_level_Ms"] = parse_int_list(args[i + 1])
            elseif key == "mg-level-dt-scales" || key == "mg_level_dt_scales"
                opts["mg_level_dt_scales"] = parse_float_list(args[i + 1])
            elseif key == "debug-residual" || key == "debug_residual"
                v = lowercase(args[i + 1])
                opts["debug_residual"] = (v == "1" || v == "true" || v == "yes" || v == "on")
            elseif key == "debug-vcycle" || key == "debug_vcycle"
                v = lowercase(args[i + 1])
                opts["debug_vcycle"] = (v == "1" || v == "true" || v == "yes" || v == "on")
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
    lap_order = Symbol(lowercase(string(opts["lap_order"])))
    (lap_order === :second || lap_order === :fourth) ||
        error("lap-order must be second/fourth")
    solver_mode = lowercase(opts["solver"])
    solver = :taylor
    cg_precond = Symbol(lowercase(opts["cg_precond"]))
    mg_interval = Int(opts["mg_interval"])
    mg_M = Int(opts["mg_M"])
    mg_dt_scale = Float64(opts["mg_dt_scale"])
    mg_nu1 = Int(opts["mg_nu1"])
    mg_nu2 = Int(opts["mg_nu2"])
    mg_level_Ms = opts["mg_level_Ms"]
    mg_level_dt_scales = opts["mg_level_dt_scales"]
    mg_corr_M = Int(opts["mg_corr_M"])
    mg_corr_dt_scale = Float64(opts["mg_corr_dt_scale"])
    mg_corr_steps = Int(opts["mg_corr_steps"])
    corr_nu1_set = opts["mg_corr_nu1"] !== nothing
    corr_nu2_set = opts["mg_corr_nu2"] !== nothing
    mg_corr_nu1 = corr_nu1_set ? Int(opts["mg_corr_nu1"]) : mg_corr_steps
    mg_corr_nu2 = corr_nu2_set ? Int(opts["mg_corr_nu2"]) : mg_corr_steps
    mg_vcycle = false
    mg_vcycle_mode = "uniform"
    mg_correction = "classic"
    if solver_mode == "taylor" || solver_mode == "sor" || solver_mode == "ssor" || solver_mode == "cg"
        solver = Symbol(solver_mode)
    elseif solver_mode == "mg-uniform-taylor"
        solver = :taylor
        mg_vcycle = true
        mg_vcycle_mode = "uniform"
        mg_correction = "classic"
    elseif solver_mode == "mg-hierarchical-taylor"
        solver = :taylor
        mg_vcycle = true
        mg_vcycle_mode = "hierarchical"
        mg_correction = "classic"
    elseif solver_mode == "mg-correction-taylor"
        solver = :taylor
        mg_vcycle = true
        mg_vcycle_mode = "uniform"
        mg_correction = "correction-taylor"
    else
        error("solver must be taylor/sor/ssor/cg/mg-uniform-taylor/mg-hierarchical-taylor/mg-correction-taylor")
    end
    use_hierarchical_taylor = (mg_vcycle_mode == "hierarchical")
    if mg_vcycle && mg_interval == 0
        @warn "mg solver selected; defaulting --mg-interval to 5"
        mg_interval = 5
    end
    debug_residual = Bool(opts["debug_residual"])
    debug_vcycle = Bool(opts["debug_vcycle"])
    (cg_precond === :ssor || cg_precond === :none) ||
        error("cg-precond must be ssor/none")
    if !mg_vcycle &&
       ((mg_level_Ms !== nothing) || (mg_level_dt_scales !== nothing))
        @warn "mg-level-Ms/mg-level-dt-scales are ignored unless --solver mg-hierarchical-taylor"
        mg_level_Ms = nothing
        mg_level_dt_scales = nothing
    end
    if mg_vcycle && !use_hierarchical_taylor &&
       ((mg_level_Ms !== nothing) || (mg_level_dt_scales !== nothing))
        @warn "mg-level-Ms/mg-level-dt-scales are ignored unless --solver mg-hierarchical-taylor"
        mg_level_Ms = nothing
        mg_level_dt_scales = nothing
    end
    if mg_correction == "classic" &&
       (mg_corr_M != 2 || mg_corr_dt_scale != 1.0 || mg_corr_steps != 1 ||
        corr_nu1_set || corr_nu2_set)
        @warn "mg-corr-* options are ignored unless --solver mg-correction-taylor"
    end
    output_dir = string(opts["output_dir"])
    run_dir = make_run_dir(output_dir)
    if solver !== :taylor && bc_order !== :spec
        @warn "bc-order is only valid for taylor; forcing to spec for iterative solvers" solver=solver bc_order=bc_order
        bc_order = :spec
    end
    if solver !== :taylor && lap_order !== :second
        @warn "lap-order is only valid for taylor; forcing to second for iterative solvers" solver=solver lap_order=lap_order
        lap_order = :second
    end
    println("run config:")
    @printf("  nx=%d ny=%d nz=%d M=%d\n", config.nx, config.ny, config.nz, config.M)
    @printf("  dt=%.3e max_steps=%d epsilon=%.3e\n", config.dt, config.max_steps, config.epsilon)
    @printf("  alpha=%.6f bc_order=%s lap_order=%s\n", prob.alpha, string(bc_order), string(lap_order))
    @printf("  solver=%s\n", solver_mode)
    if solver === :cg
        @printf("  cg_precond=%s\n", string(cg_precond))
    end
    if mg_vcycle && mg_interval > 0
        @printf("  mg_vcycle=%s mg_vcycle_mode=%s mg_interval=%d mg_dt_scale=%.3f mg_M=%d mg_nu1=%d mg_nu2=%d\n",
                string(mg_vcycle), mg_vcycle_mode, mg_interval, mg_dt_scale, mg_M, mg_nu1, mg_nu2)
        @printf("  mg_correction=%s mg_corr_M=%d mg_corr_dt_scale=%.3f mg_corr_steps=%d mg_corr_nu1=%d mg_corr_nu2=%d\n",
                mg_correction, mg_corr_M, mg_corr_dt_scale, mg_corr_steps, mg_corr_nu1, mg_corr_nu2)
    end
    @printf("  output_dir=%s\n", output_dir)
    @printf("  run_dir=%s\n", run_dir)
    run_config = Dict(
        "timestamp" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SSzzzz"),
        "script" => "scripts/run_solver.jl",
        "output_dir" => output_dir,
        "run_dir" => run_dir,
        "nx" => config.nx,
        "ny" => config.ny,
        "nz" => config.nz,
        "max_steps" => config.max_steps,
        "epsilon" => config.epsilon,
        "alpha" => prob.alpha,
        "bc_order" => string(bc_order),
        "lap_order" => string(lap_order),
        "solver" => solver_mode,
        "warmup" => true,
    )
    if opts["config_path"] !== nothing
        run_config["config_path"] = opts["config_path"]
    end
    if solver === :taylor
        run_config["M"] = config.M
        run_config["dt"] = config.dt
        run_config["Fo"] = fo
        run_config["dt_source"] = dt_source
        run_config["mg_interval"] = mg_interval
        run_config["mg_dt_scale"] = mg_dt_scale
        run_config["mg_M"] = mg_M
        run_config["mg_nu1"] = mg_nu1
        run_config["mg_nu2"] = mg_nu2
        run_config["mg_correction"] = mg_correction
        run_config["mg_corr_M"] = mg_corr_M
        run_config["mg_corr_dt_scale"] = mg_corr_dt_scale
        run_config["mg_corr_steps"] = mg_corr_steps
        run_config["mg_corr_nu1"] = mg_corr_nu1
        run_config["mg_corr_nu2"] = mg_corr_nu2
        if mg_level_Ms !== nothing
            run_config["mg_level_Ms"] = mg_level_Ms
        end
        if mg_level_dt_scales !== nothing
            run_config["mg_level_dt_scales"] = mg_level_dt_scales
        end
        run_config["debug_residual"] = debug_residual
        run_config["debug_vcycle"] = debug_vcycle
        if mg_vcycle && mg_interval > 0
            levels, coarsest = mg_levels_and_coarsest(config)
            run_config["mg_levels_used"] = levels
            run_config["mg_coarsest_grid"] = coarsest
        end
    elseif solver === :cg
        run_config["cg_precond"] = string(cg_precond)
    end
    open(joinpath(run_dir, "run_config.toml"), "w") do io
        TOML.print(io, run_config)
    end
    warmup_solve(config, prob, bc_order, lap_order, run_dir, solver, cg_precond)
    sol, runtime = if solver === :taylor
        correction_mode = (mg_correction == "correction-taylor") ? :correction_taylor : :classic
        solve_with_runtime(config, prob; bc_order=bc_order, lap_order=lap_order, output_dir=run_dir,
                           mg_interval=mg_interval, mg_dt_scale=mg_dt_scale, mg_M=mg_M,
                           mg_vcycle=mg_vcycle, mg_nu1=mg_nu1, mg_nu2=mg_nu2,
                           mg_level_Ms=(use_hierarchical_taylor ? mg_level_Ms : nothing),
                           mg_level_dt_scales=(use_hierarchical_taylor ? mg_level_dt_scales : nothing),
                           mg_correction=correction_mode,
                           mg_corr_M=mg_corr_M, mg_corr_dt_scale=mg_corr_dt_scale, mg_corr_steps=mg_corr_steps,
                           mg_corr_nu1=mg_corr_nu1, mg_corr_nu2=mg_corr_nu2,
                           debug_residual=debug_residual, debug_vcycle=debug_vcycle)
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
    res_l2 = last_res_l2(joinpath(run_dir, history_file))
    mg_levels_used = nothing
    mg_coarsest_grid = nothing
    if mg_vcycle && mg_interval > 0
        levels, coarsest = mg_levels_and_coarsest(config)
        mg_levels_used = levels
        mg_coarsest_grid = coarsest
    end
    write_run_summary(joinpath(run_dir, "run_summary.toml");
                      steps=sol.iter, runtime_s=runtime, res_l2=res_l2,
                      err_l2=err_l2, err_max=err_max,
                      history_file=history_file,
                      error_plot="error_$(tag).png",
                      exact_plot="exact_nx$(config.nx)_ny$(config.ny)_nz$(config.nz).png",
                      mg_levels_used=mg_levels_used, mg_coarsest_grid=mg_coarsest_grid)
    @info "done" t=sol.t iter=sol.iter
end

main()
