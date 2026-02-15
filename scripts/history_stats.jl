#!/usr/bin/env julia

using YAML
using Printf

function parse_args(args)
    input_dir = "results"
    output_path = ""
    i = 1
    while i <= length(args)
        if startswith(args[i], "--")
            key = args[i][3:end]
            i == length(args) && error("Missing value for --$key")
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
    if output_path == ""
        output_path = joinpath(input_dir, "history_stats.json")
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

function read_history_residuals(path::AbstractString)
    vals = Float64[]
    open(path, "r") do io
        for line in eachline(io)
            s = strip(line)
            if isempty(s) || startswith(s, "#")
                continue
            end
            parts = split(s)
            length(parts) < 3 && continue
            r = tryparse(Float64, parts[3])
            r === nothing && continue
            push!(vals, r)
        end
    end
    isempty(vals) && error("no residual rows in history: $(path)")
    return vals
end

function is_monotonic_nonincreasing(vals::Vector{Float64})
    for i in 2:length(vals)
        if vals[i] > vals[i - 1]
            return false
        end
    end
    return true
end

function oscillation_detected(vals::Vector{Float64})
    if length(vals) < 3
        return false
    end
    prev_sign = 0
    for i in 2:length(vals)
        d = vals[i] - vals[i - 1]
        sign = d > 0 ? 1 : (d < 0 ? -1 : 0)
        if sign == 0
            continue
        end
        if prev_sign != 0 && sign != prev_sign
            return true
        end
        prev_sign = sign
    end
    return false
end

function convergence_rate_estimate(vals::Vector{Float64})
    pts = Tuple{Float64,Float64}[]
    for i in eachindex(vals)
        v = vals[i]
        if v > 0 && isfinite(v)
            push!(pts, (i - 1.0, log(v)))
        end
    end
    if length(pts) < 2
        return 0.0
    end
    n = length(pts)
    sx = 0.0
    sy = 0.0
    sxx = 0.0
    sxy = 0.0
    for (x, y) in pts
        sx += x
        sy += y
        sxx += x * x
        sxy += x * y
    end
    denom = n * sxx - sx * sx
    abs(denom) < eps(Float64) && return 0.0
    slope = (n * sxy - sx * sy) / denom
    return -slope
end

function json_bool(x::Bool)
    return x ? "true" : "false"
end

function json_num(x)
    if x isa AbstractFloat
        if isnan(x)
            return "\"nan\""
        elseif isinf(x)
            return signbit(x) ? "\"-inf\"" : "\"inf\""
        end
    end
    return repr(x)
end

function write_stats_json(path::AbstractString, stats::Dict{String,Any})
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "{")
        println(io, "  \"history_stats\": {")
        println(io, "    \"monotonic\": ", json_bool(stats["monotonic"]), ",")
        println(io, "    \"oscillation_detected\": ", json_bool(stats["oscillation_detected"]), ",")
        println(io, "    \"converged\": ", json_bool(stats["converged"]), ",")
        println(io, "    \"initial_residual\": ", json_num(stats["initial_residual"]), ",")
        println(io, "    \"min_residual\": ", json_num(stats["min_residual"]), ",")
        println(io, "    \"final_residual\": ", json_num(stats["final_residual"]), ",")
        println(io, "    \"convergence_rate_estimate\": ", json_num(stats["convergence_rate_estimate"]))
        println(io, "  }")
        println(io, "}")
    end
end

function main()
    input_dir, output_path = parse_args(ARGS)
    summary_path = joinpath(input_dir, "run_summary.json")
    isfile(summary_path) || error("run_summary.json not found: $(summary_path)")

    run_summary = YAML.load_file(summary_path)
    artifacts = safe_get(run_summary, "artifacts", Dict{String,Any}())
    history_rel = safe_get(artifacts, "history", nothing)
    history_rel === nothing && error("artifacts.history not found in run_summary.json")
    history_path = joinpath(input_dir, string(history_rel))
    isfile(history_path) || error("history file not found: $(history_path)")

    residuals = read_history_residuals(history_path)
    converged = Bool(safe_get(run_summary, "converged", false))
    stats = Dict{String,Any}(
        "monotonic" => is_monotonic_nonincreasing(residuals),
        "oscillation_detected" => oscillation_detected(residuals),
        "converged" => converged,
        "initial_residual" => residuals[1],
        "min_residual" => minimum(residuals),
        "final_residual" => residuals[end],
        "convergence_rate_estimate" => convergence_rate_estimate(residuals),
    )
    write_stats_json(output_path, stats)
    @printf("wrote %s\n", output_path)
end

main()
