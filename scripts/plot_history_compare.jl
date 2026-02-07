#!/usr/bin/env julia

using TOML
using Printf
using Plots

function parse_args(args)
    input_dirs = String[]
    input_dir = "results"
    output_dir = ""
    recursive = false
    input_globs = String[]
    i = 1
    while i <= length(args)
        if startswith(args[i], "--")
            key = args[i][3:end]
            if i == length(args)
                error("Missing value for --$key")
            end
            if key == "input-dir" || key == "input_dir"
                input_dir = args[i + 1]
            elseif key == "input-dirs" || key == "input_dirs"
                input_dirs = String.(split(args[i + 1], ","))
            elseif key == "input-glob" || key == "input_glob"
                input_globs = String.(split(args[i + 1], ","))
            elseif key == "output-dir" || key == "output_dir"
                output_dir = args[i + 1]
            elseif key == "recursive"
                recursive = lowercase(args[i + 1]) in ("1", "true", "yes")
            else
                error("Unknown option: --$key")
            end
            i += 2
        else
            i += 1
        end
    end
    if isempty(input_dirs)
        input_dirs = [input_dir]
    end
    return input_dirs, output_dir, recursive, input_globs
end

function label_from_config(run_dir::String)
    cfg_path = joinpath(run_dir, "run_config.toml")
    if !isfile(cfg_path)
        return basename(run_dir)
    end
    cfg = TOML.parsefile(cfg_path)
    solver = get(cfg, "solver", "unknown")
    precond = get(cfg, "cg_precond", "none")
    nx = get(cfg, "nx", "")
    solver_sym = solver == "taylor" ? "T" :
                 solver == "sor" ? "SOR" :
                 solver == "ssor" ? "SGS" :
                 solver == "cg" ? "CG" :
                 solver == "mg-uniform-taylor" ? "MGU" :
                 solver == "mg-hierarchical-taylor" ? "MGH" :
                 solver == "mg-correction-taylor" ? "MGC" : "?"
    precond_sym = precond == "ssor" ? "SGS" :
                  precond == "none" ? "" : "?"
    if precond_sym == ""
        return @sprintf("%s-%s", solver_sym, nx)
    end
    return @sprintf("%s-%s-%s", solver_sym, precond_sym, nx)
end

function load_history(path::String)
    steps = Float64[]
    err = Float64[]
    res = Float64[]
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            if isempty(s) || startswith(s, "#")
                continue
            end
            parts = split(s)
            if length(parts) < 3
                continue
            end
            push!(steps, parse(Float64, parts[1]))
            push!(err, parse(Float64, parts[2]))
            push!(res, parse(Float64, parts[3]))
        end
    end
    return steps, err, res
end

function glob_to_regex(pattern::String)
    buf = IOBuffer()
    for c in pattern
        if c == '*'
            print(buf, ".*")
        elseif c == '?'
            print(buf, ".")
        elseif c in ('.', '+', '(', ')', '[', ']', '{', '}', '^', '$', '|', '\\')
            print(buf, "\\", c)
        else
            print(buf, c)
        end
    end
    return Regex("^" * String(take!(buf)) * "\$")
end

function expand_globs(globs::Vector{String})
    dirs = String[]
    for g in globs
        parent = dirname(g)
        pat = basename(g)
        if parent == "" || parent == "."
            parent = "."
        end
        if !isdir(parent)
            continue
        end
        re = glob_to_regex(pat)
        for entry in readdir(parent)
            if !occursin(re, entry)
                continue
            end
            full = joinpath(parent, entry)
            if isdir(full)
                push!(dirs, full)
            end
        end
    end
    return dirs
end

function main()
    input_dirs, output_dir, recursive, input_globs = parse_args(ARGS)
    if output_dir == ""
        output_dir = input_dirs[1]
    end
    if !isempty(input_globs)
        input_dirs = expand_globs(input_globs)
        if isempty(input_dirs)
            error("no directories matched input-glob")
        end
    end
    for d in input_dirs
        if !isdir(d)
            error("input directory not found: $(d)")
        end
    end

    histories = Vector{Tuple{String, Vector{Float64}, Vector{Float64}, Vector{Float64}}}()
    if recursive
        for root in input_dirs
            for (dirpath, dirnames, filenames) in walkdir(root)
                history_files = filter(f -> startswith(f, "history_") && endswith(f, ".txt"), filenames)
                if isempty(history_files)
                    continue
                end
                label = label_from_config(dirpath)
                for hf in sort(history_files)
                    steps, err, res = load_history(joinpath(dirpath, hf))
                    if isempty(steps)
                        continue
                    end
                    push!(histories, (label, steps, err, res))
                end
            end
        end
    else
        for run_dir in input_dirs
            history_files = filter(f -> startswith(f, "history_") && endswith(f, ".txt"), readdir(run_dir))
            if isempty(history_files)
                continue
            end
            label = label_from_config(run_dir)
            for hf in sort(history_files)
                steps, err, res = load_history(joinpath(run_dir, hf))
                if isempty(steps)
                    continue
                end
                push!(histories, (label, steps, err, res))
            end
        end
    end

    if isempty(histories)
        error("no history files found under $(input_dir)")
    end

    p_err = plot(xlabel="step", ylabel="err_l2", yscale=:log10,
                 yticks=10.0 .^ (-16:1:0),
                 title="History: err_l2", legend=:right)
    for (label, steps, err, _) in histories
        plot!(p_err, steps, err; label=label)
    end
    err_path = joinpath(output_dir, "history_err_l2.png")
    png(p_err, err_path)

    p_res = plot(xlabel="step", ylabel="res_l2", yscale=:log10,
                 yticks=10.0 .^ (-16:1:0),
                 title="History: res_l2", legend=:right)
    for (label, steps, _, res) in histories
        plot!(p_res, steps, res; label=label)
    end
    res_path = joinpath(output_dir, "history_res_l2.png")
    png(p_res, res_path)

    p_cmp = plot(xlabel="step", ylabel="err_l2 / res_l2", yscale=:log10,
                 yticks=10.0 .^ (-16:1:0),
                 title="History: err_l2 vs res_l2", legend=:right)
    for (label, steps, err, res) in histories
        plot!(p_cmp, steps, err; label="$(label) err", linestyle=:solid)
        plot!(p_cmp, steps, res; label="$(label) res", linestyle=:dash)
    end
    cmp_path = joinpath(output_dir, "history_err_vs_res.png")
    png(p_cmp, cmp_path)

    println("wrote $(err_path)")
    println("wrote $(res_path)")
    println("wrote $(cmp_path)")
end

main()
