using Test
using YAML

@testset "cli mg level arrays" begin
    mktempdir() do dir
        script = normpath(joinpath(@__DIR__, "..", "scripts", "run_solver.jl"))
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) $(script) --solver mg-hierarchical-taylor --n 8 --Fo 0.1 --M 2 --max-steps 1 --epsilon 1e-6 --alpha 1.0 --bc-order spec --mg-interval 1 --mg-dt-scale 2.0 --mg-M 2 --mg-nu1 1 --mg-nu2 1 --mg-level-Ms 4,2,2 --mg-level-dt-scales 2.0,4.0 --output-dir $(dir)`
        output = read(addenv(cmd, "GKSwstype" => "100"), String)
        m = match(r"run_dir=([^\n]+)", output)
        @test m !== nothing
        run_dir = strip(m.captures[1])
        summary_path = joinpath(run_dir, "run_summary.json")
        @test isfile(summary_path)
        sum = YAML.load_file(summary_path)
        @test haskey(sum, "config_path")
    end
end

@testset "cli mg correction config" begin
    mktempdir() do dir
        script = normpath(joinpath(@__DIR__, "..", "scripts", "run_solver.jl"))
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) $(script) --solver mg-correction-taylor --n 8 --Fo 0.1 --M 2 --max-steps 1 --epsilon 1e-6 --alpha 1.0 --bc-order spec --mg-interval 1 --mg-dt-scale 2.0 --mg-M 2 --mg-nu1 1 --mg-nu2 1 --mg-corr-M 3 --mg-corr-dt-scale 1.5 --mg-corr-steps 2 --mg-corr-nu1 1 --mg-corr-nu2 3 --output-dir $(dir)`
        output = read(addenv(cmd, "GKSwstype" => "100"), String)
        m = match(r"run_dir=([^\n]+)", output)
        @test m !== nothing
        run_dir = strip(m.captures[1])
        sum = YAML.load_file(joinpath(run_dir, "run_summary.json"))
        @test haskey(sum, "converged")
        @test haskey(sum, "residual_l2")
    end
end

@testset "cli mg-uniform ignores level arrays" begin
    mktempdir() do dir
        script = normpath(joinpath(@__DIR__, "..", "scripts", "run_solver.jl"))
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) $(script) --solver mg-uniform-taylor --n 8 --Fo 0.1 --M 2 --max-steps 1 --epsilon 1e-6 --alpha 1.0 --bc-order spec --mg-interval 1 --mg-dt-scale 2.0 --mg-M 2 --mg-nu1 1 --mg-nu2 1 --mg-level-Ms 4,2,2 --mg-level-dt-scales 2.0,4.0 --output-dir $(dir)`
        output = read(addenv(cmd, "GKSwstype" => "100"), String)
        m = match(r"run_dir=([^\n]+)", output)
        @test m !== nothing
        run_dir = strip(m.captures[1])
        sum = YAML.load_file(joinpath(run_dir, "run_summary.json"))
        @test haskey(sum, "artifacts")
    end
end

@testset "cli lap-order fourth" begin
    mktempdir() do dir
        script = normpath(joinpath(@__DIR__, "..", "scripts", "run_solver.jl"))
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) $(script) --solver taylor --n 8 --Fo 0.1 --M 2 --max-steps 1 --epsilon 1e-6 --alpha 1.0 --bc-order spec --lap-order fourth --output-dir $(dir)`
        output = read(addenv(cmd, "GKSwstype" => "100"), String)
        m = match(r"run_dir=([^\n]+)", output)
        @test m !== nothing
        run_dir = strip(m.captures[1])
        sum = YAML.load_file(joinpath(run_dir, "run_summary.json"))
        @test haskey(sum, "iterations")
    end
end

@testset "run_exp config-driven" begin
    root = normpath(joinpath(@__DIR__, ".."))
    run_exp = normpath(joinpath(root, "bin", "run_exp"))

    function write_config(path, name; command, extra_run="")
        open(path, "w") do io
            write(io, """
experiment:
  name: $(name)
  description: test

output:
  results_dir: results/$(name)

execution:
  command: >
    $(command)

run:
  solver: taylor
  n: 8
  Fo: 0.1
  M: 2
  max_steps: 1
  epsilon: 1.0e-6
  alpha: 1.0
  bc_order: spec
  lap_order: second
""")
            if !isempty(extra_run)
                write(io, extra_run)
            end
        end
    end

    # 正常系
    name_ok = "_tmp_config_ok"
    exp_dir_ok = joinpath(root, "experiments", name_ok)
    mkpath(exp_dir_ok)
    config_ok = joinpath(exp_dir_ok, "config.yaml")
    cmd_ok = "julia --project=. scripts/run_solver.jl --config experiments/$(name_ok)/config.yaml"
    write_config(config_ok, name_ok; command=cmd_ok)
    output = read(addenv(Cmd(`$(run_exp) $(name_ok)`; dir=root), "GKSwstype" => "100"), String)
    @test occursin("Experiment completed", output)
    run_sum = YAML.load_file(joinpath(root, "results", name_ok, "run_summary.json"))
    @test run_sum["config_path"] == "experiments/$(name_ok)/config.yaml"

    # 異常系: 旧形式 command（実行パラメータ直書き）
    name_bad = "_tmp_config_bad"
    exp_dir_bad = joinpath(root, "experiments", name_bad)
    mkpath(exp_dir_bad)
    config_bad = joinpath(exp_dir_bad, "config.yaml")
    cmd_bad = "julia --project=. scripts/run_solver.jl --n 8 --Fo 0.1"
    write_config(config_bad, name_bad; command=cmd_bad)
    @test_throws Base.ProcessFailedException run(addenv(Cmd(`$(run_exp) $(name_bad)`; dir=root), "GKSwstype" => "100"))

    # 異常系: --config の実験名不一致
    name_mismatch = "_tmp_config_mismatch"
    exp_dir_mm = joinpath(root, "experiments", name_mismatch)
    mkpath(exp_dir_mm)
    config_mm = joinpath(exp_dir_mm, "config.yaml")
    cmd_mm = "julia --project=. scripts/run_solver.jl --config experiments/other/config.yaml"
    write_config(config_mm, name_mismatch; command=cmd_mm)
    @test_throws Base.ProcessFailedException run(addenv(Cmd(`$(run_exp) $(name_mismatch)`; dir=root), "GKSwstype" => "100"))

    rm(joinpath(root, "results", name_ok); recursive=true, force=true)
    rm(joinpath(root, "logs", "$(name_ok).json"); force=true)
    rm(exp_dir_ok; recursive=true, force=true)
    rm(exp_dir_bad; recursive=true, force=true)
    rm(exp_dir_mm; recursive=true, force=true)
end
