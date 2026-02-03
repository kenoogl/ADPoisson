# ADPoisson

Taylor級数による擬似時間発展法で 3D Poisson 方程式を解くソルバーです。

**要件**
1. Julia 1.10 以上

**セットアップ**
```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
```

**実行**
```bash
julia --project scripts/main.jl --nx 32 --ny 32 --nz 32 --M 10 --dt 1e-4 --tend 1.0 --epsilon 1e-10 --alpha 1.0
```

終了時に `Fo`、解析解との **L2誤差（絶対値）**、ステップ数、実行時間をまとめて出力します。`Fo > 0.5` の場合は警告を表示します。

**コマンドライン引数**
- `--nx`, `--ny`, `--nz`: 各方向の分割数（デフォルト: 16）
- `--M`: Taylor 展開次数（デフォルト: 10）
- `--dt`: 擬似時間刻み幅（デフォルト: 1e-4）
- `--tend`: 反復の上限時間（デフォルト: 1.0）
- `--epsilon`: 収束判定の相対残差閾値（デフォルト: 1e-10）
- `--alpha`: 境界条件パラメータ（デフォルト: 1.0）
- `--bc-order`: 境界条件の次数（`spec` または `high`、デフォルト: `spec`）

**推奨設定**
- 拡散数 `Fo = Δt(1/Δx^2 + 1/Δy^2 + 1/Δz^2)` を `0.5` 以下にする
- 立方体格子（`nx=ny=nz=n`）では `Δt <= 0.5 / (3 n^2)` が目安

**高次境界（任意）**
- 境界精度を上げたい場合は `--bc-order high` を指定

例:
```bash
julia --project -e 'using ADPoisson; n=32; dt=0.1/(3n^2); config=SolverConfig(n,n,n,4,dt,0.5,1e-6); prob,_=make_problem(config; alpha=1.0); sol=solve(config, prob; bc_order=:high)'
```

出力は `results/` に保存されます:
- `solution_nx{nx}_ny{ny}_nz{nz}_M{M}_t{tend}.png`
- `error_nx{nx}_ny{ny}_nz{nz}_M{M}_t{tend}.png`

**テスト**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

`N=64` まで含めたフルテストを行う場合:
```bash
ADPOISSON_FULL_TEST=1 julia --project -e 'using Pkg; Pkg.test()'
```

テスト中に可視化を出力する場合:
```bash
ADPOISSON_TEST_PLOT=1 julia --project -e 'using Pkg; Pkg.test()'
```
