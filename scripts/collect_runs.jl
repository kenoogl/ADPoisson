#!/usr/bin/env julia

using TOML
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

function safe_get(d::Dict{String,Any}, key::String, default="")
    return haskey(d, key) ? d[key] : default
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
        cfg_path = joinpath(run_dir, "run_config.toml")
        sum_path = joinpath(run_dir, "run_summary.toml")
        if !isfile(cfg_path) || !isfile(sum_path)
            continue
        end
        cfg = TOML.parsefile(cfg_path)
        sum = TOML.parsefile(sum_path)
        row = Dict{String,Any}()
        row["dir"] = entry
        row["nx"] = safe_get(cfg, "nx", "")
        row["ny"] = safe_get(cfg, "ny", "")
        row["nz"] = safe_get(cfg, "nz", "")
        row["solver"] = safe_get(cfg, "solver", "")
        row["cg_precond"] = safe_get(cfg, "cg_precond", "")
        row["alpha"] = safe_get(cfg, "alpha", "")
        row["epsilon"] = safe_get(cfg, "epsilon", "")
        row["bc_order"] = safe_get(cfg, "bc_order", "")
        row["steps"] = safe_get(sum, "steps", "")
        row["err_l2"] = safe_get(sum, "err_l2", "")
        row["err_max"] = safe_get(sum, "err_max", "")
        row["runtime_s"] = safe_get(sum, "runtime_s", "")
        row["res_l2"] = safe_get(sum, "res_l2", "")
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
