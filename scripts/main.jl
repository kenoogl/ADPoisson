#!/usr/bin/env julia

using ADPoisson

function parse_args(args)
    opts = Dict{String,Any}()
    opts["nx"] = 16
    opts["ny"] = 16
    opts["nz"] = 16
    opts["M"] = 10
    opts["dt"] = 1e-4
    opts["tend"] = 1.0
    opts["epsilon"] = 1e-10
    opts["alpha"] = 1.0
    opts["bc_order"] = "spec"

    i = 1
    while i <= length(args)
        if startswith(args[i], "--")
            key = args[i][3:end]
            if i == length(args)
                error("Missing value for --$key")
            end
            if key == "bc-order" || key == "bc_order"
                opts["bc_order"] = args[i + 1]
            else
                opts[key] = parse(Float64, args[i + 1])
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
    config = SolverConfig(Int(opts["nx"]), Int(opts["ny"]), Int(opts["nz"]),
                          Int(opts["M"]), opts["dt"], opts["tend"], opts["epsilon"])
    prob, _ = make_problem(config; alpha=opts["alpha"])
    bc_order = Symbol(opts["bc_order"])
    sol = solve(config, prob; bc_order=bc_order)
    plot_slice(sol, prob, config)
    @info "done" t=sol.t iter=sol.iter
end

main()
