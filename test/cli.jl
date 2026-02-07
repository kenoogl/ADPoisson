using Test
using TOML

@testset "cli mg level arrays" begin
    mktempdir() do dir
        script = normpath(joinpath(@__DIR__, "..", "scripts", "main.jl"))
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) $(script) --solver mg-hierarchical-taylor --n 8 --Fo 0.1 --M 2 --max-steps 1 --epsilon 1e-6 --alpha 1.0 --bc-order spec --mg-interval 1 --mg-dt-scale 2.0 --mg-M 2 --mg-nu1 1 --mg-nu2 1 --mg-level-Ms 4,2,2 --mg-level-dt-scales 2.0,4.0 --output-dir $(dir)`
        output = read(setenv(cmd, "GKSwstype" => "100"), String)
        m = match(r"run_dir=([^\n]+)", output)
        @test m !== nothing
        run_dir = strip(m.captures[1])
        config_path = joinpath(run_dir, "run_config.toml")
        @test isfile(config_path)
        cfg = TOML.parsefile(config_path)
        @test cfg["mg_level_Ms"] == [4, 2, 2]
        @test cfg["mg_level_dt_scales"] == [2.0, 4.0]
    end
end

@testset "cli mg correction config" begin
    mktempdir() do dir
        script = normpath(joinpath(@__DIR__, "..", "scripts", "main.jl"))
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) $(script) --solver mg-correction-taylor --n 8 --Fo 0.1 --M 2 --max-steps 1 --epsilon 1e-6 --alpha 1.0 --bc-order spec --mg-interval 1 --mg-dt-scale 2.0 --mg-M 2 --mg-nu1 1 --mg-nu2 1 --mg-corr-M 3 --mg-corr-dt-scale 1.5 --mg-corr-steps 2 --mg-corr-nu1 1 --mg-corr-nu2 3 --output-dir $(dir)`
        output = read(setenv(cmd, "GKSwstype" => "100"), String)
        m = match(r"run_dir=([^\n]+)", output)
        @test m !== nothing
        run_dir = strip(m.captures[1])
        cfg = TOML.parsefile(joinpath(run_dir, "run_config.toml"))
        @test cfg["mg_correction"] == "correction-taylor"
        @test cfg["mg_corr_M"] == 3
        @test cfg["mg_corr_dt_scale"] == 1.5
        @test cfg["mg_corr_steps"] == 2
        @test cfg["mg_corr_nu1"] == 1
        @test cfg["mg_corr_nu2"] == 3
    end
end

@testset "cli mg-uniform ignores level arrays" begin
    mktempdir() do dir
        script = normpath(joinpath(@__DIR__, "..", "scripts", "main.jl"))
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) $(script) --solver mg-uniform-taylor --n 8 --Fo 0.1 --M 2 --max-steps 1 --epsilon 1e-6 --alpha 1.0 --bc-order spec --mg-interval 1 --mg-dt-scale 2.0 --mg-M 2 --mg-nu1 1 --mg-nu2 1 --mg-level-Ms 4,2,2 --mg-level-dt-scales 2.0,4.0 --output-dir $(dir)`
        output = read(setenv(cmd, "GKSwstype" => "100"), String)
        m = match(r"run_dir=([^\n]+)", output)
        @test m !== nothing
        run_dir = strip(m.captures[1])
        cfg = TOML.parsefile(joinpath(run_dir, "run_config.toml"))
        @test cfg["solver"] == "mg-uniform-taylor"
        @test !haskey(cfg, "mg_level_Ms")
        @test !haskey(cfg, "mg_level_dt_scales")
    end
end
