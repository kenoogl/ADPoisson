using Test
using TOML

@testset "cli mg level arrays" begin
    mktempdir() do dir
        script = normpath(joinpath(@__DIR__, "..", "scripts", "main.jl"))
        cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) $(script) --solver taylor --n 8 --Fo 0.1 --M 2 --max-steps 1 --epsilon 1e-6 --alpha 1.0 --bc-order spec --mg-level 3 --mg-interval 1 --mg-dt-scale 2.0 --mg-M 2 --mg-nu1 1 --mg-nu2 1 --mg-level-Ms 4,2,2 --mg-level-dt-scales 2.0,4.0 --output-dir $(dir)`
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
