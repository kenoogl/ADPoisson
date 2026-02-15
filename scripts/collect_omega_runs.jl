#!/usr/bin/env julia

using YAML
using Printf

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

function parse_args(args)
    input_dir = "results"
    output_path = ""
    output_file = "omega_runs_summary.md"
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
        input_cfg = cfg_get(cfg, "collect", "input_dir")
        output_cfg = cfg_get(cfg, "collect", "output_path")
        output_file_cfg = cfg_get(cfg, "collect", "output_file")
        results_dir = cfg_get(cfg, "output", "results_dir")
        if input_cfg !== nothing
            input_dir = string(input_cfg)
        end
        if output_cfg !== nothing
            output_path = string(output_cfg)
        elseif results_dir !== nothing
            output_path = joinpath(string(results_dir), output_file)
        end
        if output_file_cfg !== nothing
            output_file = string(output_file_cfg)
            if output_path == "" && results_dir !== nothing
                output_path = joinpath(string(results_dir), output_file)
            end
        end
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
            if key == "input-dir" || key == "input_dir"
                input_dir = args[i + 1]
            elseif key == "output" || key == "output_path"
                output_path = args[i + 1]
            elseif key == "output-file" || key == "output_file"
                output_file = args[i + 1]
            else
                error("Unknown option: --$key")
            end
            i += 2
        else
            i += 1
        end
    end
    if output_path == ""
        output_path = joinpath(input_dir, output_file)
    end
    return input_dir, output_path
end

function safe_get(d::AbstractDict, key::String, default=nothing)
    if haskey(d, key)
        return d[key]
    elseif haskey(d, Symbol(key))
        return d[Symbol(key)]
    end
    return default
end

function try_parse_omega_from_dir(dir::String)
    m = match(r"omega([0-9]+(?:\\.[0-9]+)?)$", dir)
    return m === nothing ? nothing : parse(Float64, m.captures[1])
end

function format_residual(x)
    if x === nothing
        return ""
    end
    if x isa Number
        if isnan(x)
            return "nan"
        elseif isinf(x)
            return "inf"
        end
        return @sprintf("%.6e", x)
    end
    s = lowercase(string(x))
    if s == "nan" || s == "inf" || s == "+inf" || s == "-inf"
        return s
    end
    return string(x)
end

function is_diverged_residual(x)
    if x === nothing
        return false
    elseif x isa Number
        return !isfinite(x)
    else
        s = lowercase(string(x))
        return s == "nan" || s == "inf" || s == "+inf" || s == "-inf"
    end
end

function format_runtime(x)
    if x === nothing
        return ""
    end
    if x isa Number
        return @sprintf("%.6e", x)
    end
    return string(x)
end

function collect_rows(input_dir::String)
    if !isdir(input_dir)
        error("input directory not found: $(input_dir)")
    end

    rows = Dict("sor" => Vector{Dict{String,Any}}(), "ssor" => Vector{Dict{String,Any}}())
    for entry in sort(readdir(input_dir))
        run_dir = joinpath(input_dir, entry)
        if !isdir(run_dir)
            continue
        end
        sum_path = joinpath(run_dir, "run_summary.json")
        if !isfile(sum_path)
            continue
        end
        sum = YAML.load_file(sum_path)
        cfg = Dict{String,Any}()
        cfg_path = safe_get(sum, "config_path", nothing)
        if cfg_path !== nothing && cfg_path != ""
            full_cfg_path = isabspath(cfg_path) ? cfg_path : joinpath(pwd(), cfg_path)
            if isfile(full_cfg_path)
                loaded = YAML.load_file(full_cfg_path)
                run_section = cfg_get(loaded, "run")
                if run_section isa AbstractDict
                    cfg = Dict{String,Any}(string(k) => v for (k, v) in run_section)
                end
            end
        end
        solver = lowercase(string(safe_get(cfg, "solver", "")))
        if !(solver in ("sor", "ssor"))
            continue
        end
        nx = safe_get(cfg, "nx", safe_get(cfg, "n", nothing))
        omega = safe_get(cfg, "omega", nothing)
        if omega === nothing
            omega = try_parse_omega_from_dir(entry)
        end
        if nx === nothing || omega === nothing
            continue
        end

        res_l2 = safe_get(sum, "residual_l2", nothing)
        steps = safe_get(sum, "iterations", nothing)
        steps_display = is_diverged_residual(res_l2) ? "∞" : string(steps)

        row = Dict{String,Any}()
        row["dir"] = entry
        row["nx"] = Int(nx)
        row["omega"] = Float64(omega)
        row["steps_display"] = steps_display
        row["res_l2"] = res_l2
        row["runtime_s"] = safe_get(sum, "runtime_sec", nothing)
        push!(rows[solver], row)
    end

    for solver in ("sor", "ssor")
        sort!(rows[solver], by=r -> (r["nx"], r["omega"]))
    end
    return rows
end

function write_markdown(path::String, rows)
    open(path, "w") do io
        println(io, "# Omega Sweep Summary")
        println(io)
        for solver in ("sor", "ssor")
            println(io, "## ", uppercase(solver))
            println(io)
            println(io, "| nx | omega | steps | res_l2 | runtime_s | dir |")
            println(io, "| --- | --- | --- | --- | --- | --- |")
            for r in rows[solver]
                omega_str = @sprintf("%.1f", r["omega"])
                println(io,
                    "| ", r["nx"],
                    " | ", omega_str,
                    " | ", r["steps_display"],
                    " | ", format_residual(r["res_l2"]),
                    " | ", format_runtime(r["runtime_s"]),
                    " | ", r["dir"], " |")
            end
            println(io)
        end
    end
end

function main()
    input_dir, output_path = parse_args(ARGS)
    rows = collect_rows(input_dir)
    mkpath(dirname(output_path))
    write_markdown(output_path, rows)
    println("wrote $(output_path)")
end

main()
