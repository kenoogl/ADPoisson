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
julia --project scripts/main.jl --nx 32 --ny 32 --nz 32 --M 10 --dt 1e-4 --max-steps 10000 --epsilon 1e-10 --alpha 1.0
```

終了時に `Fo`、解析解との **L2誤差（絶対値）**、**最大誤差**、ステップ数、実行時間をまとめて出力します。`Fo > 0.5` の場合は警告を表示します。

**コマンドライン引数**
- `--nx`, `--ny`, `--nz`: 各方向の分割数（デフォルト: 16）
- `--M`: Taylor 展開次数（デフォルト: 10）
- `--dt`: 擬似時間刻み幅（デフォルト: 1e-4）
- `--max-steps`: 擬似時間積分の最大ステップ数（デフォルト: 10000）
- `--epsilon`: 収束判定の相対残差閾値（デフォルト: 1e-10）
  - 相対残差は $\|r\|_2 / \max(\|r_0\|_2, 1)$（$r=Lu-f$、$r_0$ は初期残差、内点のみ評価）
- `--alpha`: 境界条件パラメータ（デフォルト: 1.0）
- `--bc-order`: 境界条件の次数（`spec` または `high`、デフォルト: `spec`）

**推奨設定**
- 拡散数 `Fo = Δt(1/Δx^2 + 1/Δy^2 + 1/Δz^2)` を `0.5` 以下にする
- 立方体格子（`nx=ny=nz=n`）では `Δt <= 0.5 / (3 n^2)` が目安

**高次境界（任意）**
- 境界精度を上げたい場合は `--bc-order high` を指定

例:
```bash
julia --project -e 'using ADPoisson; n=32; dt=0.1/(3n^2); max_steps=Int(ceil(0.5/dt)); config=SolverConfig(n,n,n,4,dt,max_steps,1e-6); prob,_=make_problem(config; alpha=1.0); sol=solve(config, prob; bc_order=:high)'
```

出力は `results/` に保存されます:
- `exact_nx{nx}_ny{ny}_nz{nz}.png`（解析解のため格子情報のみ）
- `error_nx{nx}_ny{ny}_nz{nz}_M{M}_steps{steps}.png`
- `history_nx{nx}_ny{ny}_nz{nz}_M{M}_steps{steps}.txt`
  - 擬似時間ステップの履歴（`step`, `err_l2`, `res_l2`。`res_l2` は初期残差で相対化した残差L2、`step` は更新回数で初期状態は 0）

**テスト**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

デフォルトのテストは **問題を解いて精度確認** のみを実行します。
（`solver_error` で擬似時間ソルバを実行し、解析解との相対 L2 誤差を評価）
テストで使用する `dt` と `max_steps` は `scripts/main.jl` のデフォルト値に一致します。
ただし **フルテスト（N=64）を安定に通すため、テストでは `dt` を上限 `0.5/(3 n^2)` で安全側にクリップ**します（デフォルト有効）。
実行時と同じ `dt` を強制する場合は以下を使います。
```bash
ADPOISSON_TEST_USE_SAFE_DT=0 julia --project -e 'using Pkg; Pkg.test()'
```

実装確認（`laplacian!`, `taylor_step!`, `convergence_order`, `problems`）を行う場合:
```bash
ADPOISSON_IMPL_TEST=1 julia --project -e 'using Pkg; Pkg.test()'
```

`N=64` まで含めたフルテストを行う場合:
```bash
ADPOISSON_FULL_TEST=1 julia --project -e 'using Pkg; Pkg.test()'
```

テスト中に可視化を出力する場合:
```bash
ADPOISSON_TEST_PLOT=1 julia --project -e 'using Pkg; Pkg.test()'
```
