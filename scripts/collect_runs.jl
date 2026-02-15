#!/usr/bin/env julia

using YAML
using Dates
using Printf

function parse_args(args)
    input_dir = "results"
    output_path = ""
    i = 1
    while i <= length(args)
        if startswith(args[i], "--")
            key = args[i][3:end]
            if i == length(args)
                error("Missing value for --$key")
            end
            if key == "input-dir" || key == "input_dir"
                input_dir = args[i + 1]
            elseif key == "output" || key == "output_path"
                output_path = args[i + 1]
            else
                error("Unknown option: --$key")
            end
            i += 2
        else
            i += 1
        end
    end
    return input_dir, output_path
end

function safe_get(d::AbstractDict, key::String, default="")
    if haskey(d, key)
        return d[key]
    elseif haskey(d, Symbol(key))
        return d[Symbol(key)]
    end
    return default
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

function format_float(x)
    if x === "" || x === nothing
        return ""
    end
    if x isa Number
        return @sprintf("%.6e", x)
    end
    return string(x)
end

function format_int(x)
    if x === "" || x === nothing
        return ""
    end
    return string(x)
end

function collect_rows(input_dir::String)
    rows = Vector{Dict{String,Any}}()
    if !isdir(input_dir)
        error("input directory not found: $(input_dir)")
    end
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
        run_cfg = Dict{String,Any}()
        cfg_path = safe_get(sum, "config_path", nothing)
        if cfg_path !== nothing && cfg_path != ""
            full_cfg_path = isabspath(cfg_path) ? cfg_path : joinpath(pwd(), cfg_path)
            if isfile(full_cfg_path)
                loaded = YAML.load_file(full_cfg_path)
                run_section = cfg_get(loaded, "run")
                if run_section isa AbstractDict
                    run_cfg = Dict{String,Any}(string(k) => v for (k, v) in run_section)
                end
            end
        end
        row = Dict{String,Any}()
        row["dir"] = entry
        row["nx"] = safe_get(run_cfg, "nx", safe_get(run_cfg, "n", ""))
        row["ny"] = safe_get(run_cfg, "ny", safe_get(run_cfg, "n", ""))
        row["nz"] = safe_get(run_cfg, "nz", safe_get(run_cfg, "n", ""))
        row["solver"] = safe_get(run_cfg, "solver", "")
        row["cg_precond"] = safe_get(run_cfg, "cg_precond", "")
        row["alpha"] = safe_get(run_cfg, "alpha", "")
        row["epsilon"] = safe_get(run_cfg, "epsilon", "")
        row["bc_order"] = safe_get(run_cfg, "bc_order", "")
        row["steps"] = safe_get(sum, "iterations", "")
        row["err_l2"] = safe_get(sum, "error_l2", "")
        row["err_max"] = safe_get(sum, "error_max", "")
        row["runtime_s"] = safe_get(sum, "runtime_sec", "")
        row["res_l2"] = safe_get(sum, "residual_l2", "")
        push!(rows, row)
    end
    return rows
end

function write_markdown(path::String, rows)
    open(path, "w") do io
        println(io, "# Run Summary")
        println(io)
        println(io, "| dir | nx | ny | nz | solver | precond | alpha | epsilon | bc_order | steps | err_l2 | err_max | res_l2 | runtime_s |")
        println(io, "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for r in rows
            println(io, "| $(r["dir"]) | $(format_int(r["nx"])) | $(format_int(r["ny"])) | $(format_int(r["nz"])) | $(r["solver"]) | $(r["cg_precond"]) | $(format_float(r["alpha"])) | $(format_float(r["epsilon"])) | $(r["bc_order"]) | $(format_int(r["steps"])) | $(format_float(r["err_l2"])) | $(format_float(r["err_max"])) | $(format_float(r["res_l2"])) | $(format_float(r["runtime_s"])) |")
        end
    end
end

function main()
    input_dir, output_path = parse_args(ARGS)
    rows = collect_rows(input_dir)
    if output_path == ""
        output_path = joinpath(input_dir, "run_summary.md")
    end
    write_markdown(output_path, rows)
    println("wrote $(output_path)")
end

main()
