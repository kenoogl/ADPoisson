#!/usr/bin/env julia

using ADPoisson

function parse_args(args)
    opts = Dict{String,Float64}()
    opts["nx"] = 16
    opts["ny"] = 16
    opts["nz"] = 16
    opts["M"] = 10
    opts["dt"] = 1e-4
    opts["tend"] = 1.0
    opts["epsilon"] = 1e-10
    opts["alpha"] = 1.0

    i = 1
    while i <= length(args)
        if startswith(args[i], "--")
            key = args[i][3:end]
            if i == length(args)
                error("Missing value for --$key")
            end
            opts[key] = parse(Float64, args[i + 1])
            i += 2
        else
            i += 1
        end
    end
    return opts
end

function main()
    opts = parse_args(ARGS)
    config = SolverConfig(Int(opts["nx"]), Int(opts["ny"]), Int(opts["nz"]),
                          Int(opts["M"]), opts["dt"], opts["tend"], opts["epsilon"])
    prob, _ = make_problem(config; alpha=opts["alpha"])
    sol = solve(config, prob)
    plot_slice(sol, prob, config)
    @info "done" t=sol.t iter=sol.iter
end

main()
