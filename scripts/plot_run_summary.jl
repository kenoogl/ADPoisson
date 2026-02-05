#!/usr/bin/env julia

using TOML
using Printf
using Plots

function parse_args(args)
    input_dir = "results"
    output_dir = ""
    i = 1
    while i <= length(args)
        if startswith(args[i], "--")
            key = args[i][3:end]
            if i == length(args)
                error("Missing value for --$key")
            end
            if key == "input-dir" || key == "input_dir"
                input_dir = args[i + 1]
            elseif key == "output-dir" || key == "output_dir"
                output_dir = args[i + 1]
            else
                error("Unknown option: --$key")
            end
            i += 2
        else
            i += 1
        end
    end
    return input_dir, output_dir
end

function safe_get(d::Dict{String,Any}, key::String, default=nothing)
    return haskey(d, key) ? d[key] : default
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
        steps = safe_get(sum, "steps", nothing)
        if steps == 20000
            continue
        end
        row = Dict{String,Any}()
        row["dir"] = entry
        row["solver"] = safe_get(cfg, "solver", "unknown")
        row["precond"] = safe_get(cfg, "cg_precond", "none")
        row["nx"] = safe_get(cfg, "nx", 0)
        row["err_l2"] = safe_get(sum, "err_l2", NaN)
        row["res_l2"] = safe_get(sum, "res_l2", NaN)
        row["runtime_s"] = safe_get(sum, "runtime_s", NaN)
        row["steps"] = steps
        push!(rows, row)
    end
    return rows
end

function make_labels(rows)
    # group by solver+precond, then sort within by nx asc
    groups = Dict{Tuple{String,String}, Vector{Dict{String,Any}}}()
    for r in rows
        key = (string(r["solver"]), string(r["precond"]))
        if !haskey(groups, key)
            groups[key] = Vector{Dict{String,Any}}()
        end
        push!(groups[key], r)
    end
    keys_sorted = sort(collect(keys(groups)); by = k -> (k[1], k[2]))
    ordered = Vector{Dict{String,Any}}()
    labels = String[]
    for key in keys_sorted
        rs = groups[key]
        sort!(rs, by=r -> r["nx"])
        for r in rs
            push!(ordered, r)
            solver_sym = r["solver"] == "taylor" ? "T" :
                         r["solver"] == "sor" ? "SOR" :
                         r["solver"] == "ssor" ? "SGS" :
                         r["solver"] == "cg" ? "CG" : "?"
            precond_sym = r["precond"] == "ssor" ? "SGS" :
                          r["precond"] == "none" ? "" : "?"
            if precond_sym == ""
                push!(labels, @sprintf("%s-%d", solver_sym, r["nx"]))
            else
                push!(labels, @sprintf("%s-%s-%d", solver_sym, precond_sym, r["nx"]))
            end
        end
    end
    return ordered, labels
end

function plot_err(rows, labels, output_dir)
    x = 1:length(rows)
    err_l2 = [r["err_l2"] for r in rows]
    res_l2 = [r["res_l2"] for r in rows]
    p = plot(x, err_l2; label="err_l2", ylabel="error", yscale=:log10,
             yticks=10.0 .^ (-16:1:0),
             xticks=(x, labels), xrotation=60, xtickfontsize=8, bottom_margin=8Plots.mm,
             legend=:right, marker=:circle)
    plot!(p, x, res_l2; label="res_l2", marker=:diamond)
    out = joinpath(output_dir, "compare_errors.png")
    png(p, out)
end

function plot_runtime_steps(rows, labels, output_dir)
    x = 1:length(rows)
    runtime = [r["runtime_s"] for r in rows]
    steps = [r["steps"] for r in rows]
    rt_min = minimum(runtime)
    rt_max = maximum(runtime)
    rt_pad = rt_max == rt_min ? rt_max * 0.1 + 1e-12 : (rt_max - rt_min) * 0.1
    p = plot(x, runtime; label="Runtime", ylabel="Runtime",
             xticks=(x, labels), xrotation=60, xtickfontsize=8, bottom_margin=8Plots.mm,
             legend=(0.12, 0.75), legendfontsize=8, background_color_legend=:transparent,
             marker=:circle, color=:blue, ylims=(rt_min - rt_pad, rt_max + rt_pad))
    p2 = twinx(p)
    plot!(p2, x, steps; label="Iteration", ylabel="Iteration", marker=:diamond, color=:red)
    out = joinpath(output_dir, "compare_runtime_steps.png")
    png(p, out)
end

function main()
    input_dir, output_dir = parse_args(ARGS)
    if output_dir == ""
        output_dir = input_dir
    end
    rows = collect_rows(input_dir)
    if isempty(rows)
        error("no valid run_config.toml/run_summary.toml found in $(input_dir)")
    end
    ordered, labels = make_labels(rows)
    plot_err(ordered, labels, output_dir)
    plot_runtime_steps(ordered, labels, output_dir)
    println("wrote $(joinpath(output_dir, "compare_errors.png"))")
    println("wrote $(joinpath(output_dir, "compare_runtime_steps.png"))")
end

main()
