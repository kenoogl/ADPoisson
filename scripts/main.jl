#!/usr/bin/env julia

using ADPoisson
using Printf

function parse_args(args)
    opts = default_cli_options()

    i = 1
    while i <= length(args)
        if startswith(args[i], "--")
            key = args[i][3:end]
            if i == length(args)
                error("Missing value for --$key")
            end
            if key == "bc-order" || key == "bc_order"
                opts["bc_order"] = args[i + 1]
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
    config = SolverConfig(Int(opts["nx"]), Int(opts["ny"]), Int(opts["nz"]),
                          Int(opts["M"]), opts["dt"], Int(opts["max_steps"]), opts["epsilon"])
    prob, _ = make_problem(config; alpha=opts["alpha"])
    bc_order = Symbol(opts["bc_order"])
    println("run config:")
    @printf("  nx=%d ny=%d nz=%d M=%d\n", config.nx, config.ny, config.nz, config.M)
    @printf("  dt=%.3e max_steps=%d epsilon=%.3e\n", config.dt, config.max_steps, config.epsilon)
    @printf("  alpha=%.6f bc_order=%s\n", prob.alpha, string(bc_order))
    sol = solve(config, prob; bc_order=bc_order)
    plot_slice(sol, prob, config)
    @info "done" t=sol.t iter=sol.iter
end

main()
