# ADPoisson

このコードで扱える解法（実験可能）:

共通オプション（全ソルバー共通、以下の表から除外）:
`--nx/--ny/--nz` または `--n`, `--max-steps`, `--epsilon`, `--alpha`, `--output-dir`
補足: `--M`, `--dt`, `--Fo`, `--bc-order` の有効/無効は次のとおり。

| solver | `--M` | `--dt/--Fo` | `--bc-order` |
| --- | --- | --- | --- |
| `taylor` | 有効 | 有効 | 有効 |
| `sor` | 無視 | 無視 | 無視（`spec` 固定） |
| `ssor` | 無視 | 無視 | 無視（`spec` 固定） |
| `cg` | 無視 | 無視 | 無視（`spec` 固定） |
| `mg-uniform-taylor` | 有効 | 有効 | 有効 |
| `mg-hierarchical-taylor` | 有効 | 有効 | 有効 |
| `mg-correction-taylor` | 有効 | 有効 | 有効 |

| 分類 | solver | 概要 | solver固有オプション（共通を除外） |
| --- | --- | --- | --- |
| Taylor | `taylor` | 擬似時間の Taylor 展開で Poisson を反復的に解く基本ソルバー | `--M`, `--dt`/`--Fo`, `--bc-order` |
| SOR | `sor` | 赤黒 SOR による反復解法 | （なし） |
| SSOR | `ssor` | 赤黒（RB）SSOR による反復解法 | （なし） |
| CG | `cg` | 共役勾配法（必要に応じて SSOR 前処理） | `--cg-precond` |
| MG | `mg-uniform-taylor` | Taylor スムーザを用いた V-cycle | `--M`, `--dt`/`--Fo`, `--bc-order`, `--mg-interval`, `--mg-M`, `--mg-dt-scale`, `--mg-nu1`, `--mg-nu2` |
| MG | `mg-hierarchical-taylor` | レベルごとに `M` と `dt` を変えるスムーザ設定 | `--M`, `--dt`/`--Fo`, `--bc-order`, `--mg-interval`, `--mg-M`, `--mg-dt-scale`, `--mg-nu1`, `--mg-nu2`, `--mg-level-Ms`, `--mg-level-dt-scales` |
| MG | `mg-correction-taylor` | coarse 補正方程式 $L e = -r$ を Taylor 擬似時間で解く | `--M`, `--dt`/`--Fo`, `--bc-order`, `--mg-interval`, `--mg-M`, `--mg-dt-scale`, `--mg-nu1`, `--mg-nu2`, `--mg-corr-M`, `--mg-corr-dt-scale`, `--mg-corr-steps`, `--mg-corr-nu1`, `--mg-corr-nu2` |

Taylor級数による擬似時間発展法で 3D Poisson 方程式を解くソルバーです。

**要件**
1. Julia 1.10 以上

**セットアップ**
```bash
julia --project -e 'using Pkg; Pkg.instantiate()'
```

### **実行**
```bash
julia --project scripts/run_solver.jl --nx 32 --ny 32 --nz 32 --M 10 --dt 1e-4 --max-steps 10000 --epsilon 1e-10 --alpha 1.0 --output-dir results

julia --project scripts/run_solver.jl --n 32 --M 10 --dt 1e-4 --max-steps 10000 --epsilon 1e-10 --alpha 1.0 --output-dir results --solver taylor
```

終了時に `Fo`、解析解との **L2誤差（絶対値）**、**最大誤差**、ステップ数、実行時間をまとめて出力します。`Fo > 0.5` の場合は警告を表示します。
`run_summary.toml` の `runtime_s` は **JIT コンパイルを除外するためのウォームアップ実行後**の計測値です。

**コマンドライン引数**
- `--nx`, `--ny`, `--nz`: 各方向の分割数（デフォルト: 16）
- `--n`: 等方格子の分割数（`nx=ny=nz=n`、`--nx/--ny/--nz` より優先）
- `--M`: Taylor 展開次数（デフォルト: 10）
- `--dt`: 擬似時間刻み幅（デフォルト: 1e-4）
- `--Fo`: 拡散数による指定（`Fo = Δt(1/Δx^2 + 1/Δy^2 + 1/Δz^2)`、`--dt` より優先）
- `--max-steps`: 擬似時間積分の最大ステップ数（デフォルト: 10000）
- `--epsilon`: 収束判定の相対残差閾値（デフォルト: 1e-10）
  - 相対残差は $\|r\|_2 / \max(\|r_0\|_2, 1)$（$r=Lu-f$、$r_0$ は初期残差、内点のみ評価）
- `--alpha`: 境界条件パラメータ（デフォルト: 1.0）
- `--bc-order`: 境界条件の次数（`spec` または `high`、デフォルト: `spec`）
  - Taylor 系（`taylor` / `mg-*`）のみ有効（反復解法では `spec` に固定）
- `--lap-order`: ラプラシアン次数（`second` または `fourth`、デフォルト: `second`）
  - Taylor 系（`taylor` / `mg-*`）のみ有効（反復解法では `second` に固定）
  - `fourth` は 4次差分（半径2、各軸5点）を使用
- `--output-dir`: 出力ディレクトリ（デフォルト: `results`。存在しない場合は作成）
- `--solver`: 実行するソルバー（`taylor` / `sor` / `ssor` / `cg` / `mg-uniform-taylor` / `mg-hierarchical-taylor` / `mg-correction-taylor`。デフォルト: `taylor`）
- `--cg-precond`: CG の前処理（`ssor` / `none`、デフォルト: `none`）
- `--mg-interval`: MG補正の適用間隔（0 で無効、デフォルト: 0）
  - `--solver mg-*` の場合、未指定なら 5 に自動設定される（非 MG では無視される）
- レベル1/2（疑似MG/2-level MG）は実験済みで不採用のため、MG加速は V-cycle のみ対応
- `--mg-dt-scale`: MG補正ステップ用 `dt` 係数（デフォルト: 2.0）
- `--mg-M`: MG補正ステップに使う Taylor 次数（デフォルト: 4）
- `--mg-nu1`: **主方程式（u）側**の前スムージング回数（デフォルト: 1）
- `--mg-nu2`: **主方程式（u）側**の後スムージング回数（デフォルト: 1）
- `--solver mg-uniform-taylor`: 全レベルで `mg_M`, `mg_dt_scale` を共通使用
- `--solver mg-hierarchical-taylor`: レベル別に `mg_level_Ms`, `mg_level_dt_scales` を使用
- `--solver mg-correction-taylor`: coarse 補正方程式 $L e = -r$（$r=Lu-f$）を Taylor 擬似時間で解く
  - coarse へは $-r$（$f-Lu$）を制限して右辺に使う
- `--mg-corr-M`: `correction-taylor` 用の Taylor 次数（デフォルト: 2）
- `--mg-corr-dt-scale`: `correction-taylor` 用の `dt` スケール（デフォルト: 1.0）
- `--mg-corr-steps`: `correction-taylor` の反復回数（デフォルト: 1）
  - `--mg-corr-nu1/--mg-corr-nu2` 未指定時に `nu1=nu2=mg-corr-steps` として利用
- `--mg-corr-nu1`: **補正方程式（e）側**の前スムージング回数（デフォルト: `mg-corr-steps`）
- `--mg-corr-nu2`: **補正方程式（e）側**の後スムージング回数（デフォルト: `mg-corr-steps`）
- `--mg-level-Ms`: 階層 Taylor のレベル別 `M`（例: `4,4,2,2`。未指定時は全レベルで `mg_M`）
- `--mg-level-dt-scales`: 階層 Taylor のレベル別 `dt` スケール（例: `2.0,2.0,4.0,4.0`。未指定時は全レベルで `mg_dt_scale`）
  - `mg-level-Ms` と `mg-level-dt-scales` は同じレベル番号（1=最細）で対応づけて使用する
  - 配列長がレベル数より短い場合は最後の値を繰り返して適用する
  - 2つは独立指定だが実効的には連成する（同一レベルで `M` を下げるほど `dt` を大きくしすぎない方が安定しやすい）
  - `dt` は各レベルで Fo 条件によりクリップされるため、`dt-scale` を上げても効果が頭打ちになる場合がある

**MG オプションの関係**
- MG を使う場合は `--solver mg-*` を指定する。
- **主方程式（u）側**: `--mg-M` / `--mg-dt-scale` / `--mg-level-*` / `--mg-nu1/--mg-nu2` が、
  V-cycle の **u 更新スムーザ**（pre/post）を決める。
- **補正方程式（e）側**: `--mg-corr-*` が、
  coarse 補正 $L e = -r$ の **解き方と反復回数**を決める。

**推奨設定**
- 拡散数 `Fo = Δt(1/Δx^2 + 1/Δy^2 + 1/Δz^2)` を `0.5` 以下にする
- 立方体格子（`nx=ny=nz=n`）では `Δt <= 0.5 / (3 n^2)` が目安

**高次境界（任意）**
- 境界精度を上げたい場合は `--bc-order high` を指定

例:
```bash
julia --project -e 'using ADPoisson; n=32; dt=0.1/(3n^2); max_steps=Int(ceil(0.5/dt)); config=SolverConfig(n,n,n,4,dt,max_steps,1e-6); prob,_=make_problem(config; alpha=1.0); sol=solve(config, prob; bc_order=:high)'
```

出力は `--output-dir` で指定したディレクトリ配下の `run_YYYYMMDD_HHMMSS/` に保存されます（デフォルト: `results/`）。
実行条件と結果の確認用に `run_config.toml` と `run_summary.toml` を出力します。
`--solver mg-*` の場合、`mg_levels_used`（使用レベル数）と `mg_coarsest_grid`（最粗格子の `nx,ny,nz`）を記録します。
- `exact_nx{nx}_ny{ny}_nz{nz}.png`（解析解のため格子情報のみ）
- `error_nx{nx}_ny{ny}_nz{nz}_M{M}_steps{steps}.png`
- `history_nx{nx}_ny{ny}_nz{nz}_M{M}_steps{steps}.txt`
  - 擬似時間ステップの履歴（`step`, `err_l2`, `res_l2`。`res_l2` は初期残差で相対化した残差L2、`step` は更新回数で初期状態は 0）

**CG 実行例**
```bash
julia --project scripts/run_solver.jl --solver cg --cg-precond ssor --nx 32 --ny 32 --nz 32 --max-steps 2000 --epsilon 1e-8 --alpha 1.0 --output-dir results
```

**反復解法（SOR）実行例**
```bash
julia --project scripts/run_solver.jl --solver sor --nx 32 --ny 32 --nz 32 --max-steps 2000 --epsilon 1e-8 --alpha 1.0 --output-dir results
```

**Uniform Taylor (V-cycle MG) 実行例**
```bash
julia --project scripts/run_solver.jl --solver mg-uniform-taylor --n 64 --Fo 0.5 --M 4 --max-steps 20000 --epsilon 1e-8 --alpha 1.0 --bc-order high --mg-interval 5 --mg-dt-scale 2.0 --mg-M 4 --mg-nu1 1 --mg-nu2 1 --output-dir results
```
**Hierarchical Taylor (V-cycle MG) 実行例**
```bash
julia --project scripts/run_solver.jl --solver mg-hierarchical-taylor --n 64 --Fo 0.5 --M 4 --max-steps 20000 --epsilon 1e-8 --alpha 1.0 --bc-order high --mg-interval 5 --mg-dt-scale 2.0 --mg-M 4 --mg-nu1 1 --mg-nu2 1 --mg-level-Ms 4,4,2,2 --mg-level-dt-scales 2.0,2.0,4.0,4.0 --output-dir results
```
**Correction-Taylor (V-cycle MG) 実行例**
```bash
julia --project scripts/run_solver.jl --solver mg-correction-taylor --n 64 --Fo 0.5 --M 4 --max-steps 20000 --epsilon 1e-8 --alpha 1.0 --bc-order high --mg-corr-M 2 --mg-corr-dt-scale 1.0 --mg-corr-nu1 1 --mg-corr-nu2 1 --mg-interval 5 --mg-dt-scale 2.0 --mg-M 4 --mg-nu1 1 --mg-nu2 1 --output-dir results
```
補足:
- V-cycle は `--mg-interval` ごとに適用されます（`0` で無効）。
- `nx,ny,nz` が偶数で段階的に 1/2 にできる場合にのみ多段化します。
- `mg-uniform-taylor` / `mg-hierarchical-taylor` の場合、最粗格子は密行列を組み、`A \ b`（Julia の直接解法/LU ベース）で解きます。
- `mg-correction-taylor` の場合、最粗格子も correction‑taylor（Taylor 擬似時間）で解きます。

### **テスト**
```bash
julia --project -e 'using Pkg; Pkg.test()'
```

デフォルトのテストは **問題を解いて精度確認** のみを実行します。
（`solver_error` で擬似時間ソルバを実行し、解析解との相対 L2 誤差を評価）
テストで使用する `dt` と `max_steps` は `scripts/run_solver.jl` のデフォルト値に一致します。
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

### **Taylor次数比較**
指定パラメータのまま Taylor 展開次数 `M` のみを変化させ、収束解と履歴を比較します。
```bash
julia --project scripts/compare_taylor.jl --solver taylor --nx 32 --ny 32 --nz 32 --dt 1e-4 --max-steps 10000 --epsilon 1e-6 --alpha 1.0 --bc-order high --Ms 2,4,6,8,10 --output-dir results
```
反復解法を使う場合は `--solver` を指定します。`--solver sor/cg` の場合、`Ms` は無視され、`--M` の値のみが使われます。
```bash
julia --project scripts/compare_taylor.jl --solver sor --nx 32 --ny 32 --nz 32 --max-steps 10000 --epsilon 1e-6 --alpha 1.0 --bc-order high --output-dir results
```
`--Fo` を指定した場合は `dt` より優先されます。`Fo > 0.5` の場合は、比較スクリプト内で `dt` を `Fo=0.5` になるようにクリップします。
`runtime_s` はウォームアップ実行後に計測します。
出力は `--output-dir` で指定したディレクトリ配下の `run_YYYYMMDD_HHMMSS/` に保存されます（存在しない場合は作成）。
`run_config.toml` と `run_summary.toml` に実行条件と結果を記録します。
- `compare_M_nx{nx}_ny{ny}_nz{nz}_Ms{Mlist}.txt`（列: `M`, `steps`, `err_l2`, `err_max`, `runtime_s`）
- `history_compare_nx{nx}_ny{ny}_nz{nz}_Ms{Mlist}.png`
- `history_nx{nx}_ny{ny}_nz{nz}_M{M}_steps{steps}.txt`（各 M の履歴）

### **ソルバー比較（Taylor/SOR/SSOR/CG）**
同一条件でソルバーを切り替えて、収束履歴（`res_l2`）を比較します。
```bash
julia --project scripts/compare_solvers.jl --solvers taylor,sor,ssor,cg --nx 32 --ny 32 --nz 32 --M 10 --dt 1e-4 --max-steps 10000 --epsilon 1e-6 --alpha 1.0 --bc-order high --output-dir results
```
`--cg-precond` は CG を含む場合のみ使用します（`none`/`ssor`）。
```bash
julia --project scripts/compare_solvers.jl --solvers taylor,cg --cg-precond ssor --nx 32 --ny 32 --nz 32 --M 10 --dt 1e-4 --max-steps 10000 --epsilon 1e-6 --alpha 1.0 --bc-order high --output-dir results
```
出力は `--output-dir` 配下の `run_YYYYMMDD_HHMMSS/` に保存されます。
- `compare_solvers_nx{nx}_ny{ny}_nz{nz}_solver{solverlist}.txt`（列: `solver`, `steps`, `err_l2`, `err_max`, `res_l2`, `runtime_s`）
- `history_compare_solvers_nx{nx}_ny{ny}_nz{nz}_solver{solverlist}.png`
- 各ソルバーの履歴ファイル（`history_*`）

### **Fo比較**
指定パラメータのまま拡散数 `Fo` のみを変化させ、収束解と履歴を比較します。
```bash
julia --project scripts/compare_fo.jl --solver taylor --nx 32 --ny 32 --nz 32 --M 10 --max-steps 10000 --epsilon 1e-6 --alpha 1.0 --bc-order high --Fos 0.1,0.2,0.3,0.4,0.5 --output-dir results
```
反復解法を使う場合は `--solver` を指定します。`--solver sor/cg` の場合、`Fo` は `dt` にしか影響しないため結果が変わらない可能性があります。
```bash
julia --project scripts/compare_fo.jl --solver cg --cg-precond none --nx 32 --ny 32 --nz 32 --max-steps 10000 --epsilon 1e-6 --alpha 1.0 --bc-order high --Fos 0.1,0.2,0.3 --output-dir results
```
出力は `--output-dir` で指定したディレクトリ配下の `run_YYYYMMDD_HHMMSS/` に保存されます（存在しない場合は作成）。
`run_config.toml` と `run_summary.toml` に実行条件と結果を記録します。
`runtime_s` はウォームアップ実行後に計測します。
- `compare_Fo_nx{nx}_ny{ny}_nz{nz}_M{M}_Fos{Folist}.txt`（列: `Fo`, `dt`, `steps`, `err_l2`, `err_max`, `runtime_s`）
- `history_compare_Fo_nx{nx}_ny{ny}_nz{nz}_M{M}_Fos{Folist}.png`
- `history_Fo{Fo}_nx{nx}_ny{ny}_nz{nz}_M{M}_steps{steps}.txt`（各 Fo の履歴）

### **alpha比較**
指定パラメータのまま境界条件パラメータ `alpha` のみを変化させ、収束解と履歴を比較します。
```bash
julia --project scripts/compare_alpha.jl --solver taylor --nx 32 --ny 32 --nz 32 --M 10 --dt 1e-4 --max-steps 10000 --epsilon 1e-6 --bc-order high --alphas 0.5,1.0,1.5 --output-dir results
```
反復解法を使う場合は `--solver` を指定します。
```bash
julia --project scripts/compare_alpha.jl --solver sor --nx 32 --ny 32 --nz 32 --dt 1e-4 --max-steps 10000 --epsilon 1e-6 --bc-order high --alphas 0.5,1.0,1.5 --output-dir results
```
出力は `--output-dir` で指定したディレクトリ配下の `run_YYYYMMDD_HHMMSS/` に保存されます（存在しない場合は作成）。
`run_config.toml` と `run_summary.toml` に実行条件と結果を記録します。
`runtime_s` はウォームアップ実行後に計測します。
- `compare_alpha_nx{nx}_ny{ny}_nz{nz}_M{M}_alphas{Alist}.txt`（列: `alpha`, `dt`, `steps`, `err_l2`, `err_max`, `runtime_s`）
- `history_compare_alpha_nx{nx}_ny{ny}_nz{nz}_M{M}_alphas{Alist}.png`
- `history_alpha{alpha}_nx{nx}_ny{ny}_nz{nz}_M{M}_steps{steps}.txt`（各 alpha の履歴）

### **実行結果の一覧化**
指定ディレクトリ配下の `run_*/run_config.toml` と `run_summary.toml` を読み、Markdown の一覧表を生成します。
```bash
julia --project scripts/collect_runs.jl --input-dir results
```
デフォルト出力: `results/run_summary.md`

### **実行結果の比較プロット**
指定ディレクトリ配下の `run_*/` を走査し、`solver+precond+nx` を項目名として比較プロットを出力します。
`steps=20000` の項目は除外されます。
```bash
julia --project scripts/plot_run_summary.jl --input-dir results
```
出力:
- `results/compare_errors.png`（`err_l2` と `err_max` の double-Y）
- `results/compare_runtime_steps.png`（`runtime` と `steps` の double-Y）

### **履歴のまとめ描画**
`history_*.txt` から `err_l2` と `res_l2` の履歴をまとめて描画します。

必要な分だけ指定（非再帰）:
```bash
julia --project scripts/plot_history_compare.jl --input-dirs results/run_20260205_164455,results/run_20260205_170955 --output-dir results
```

指定ディレクトリ配下をすべて対象（再帰）:
```bash
julia --project scripts/plot_history_compare.jl --input-dir results --recursive true
```

ワイルドカード指定:
```bash
julia --project scripts/plot_history_compare.jl --input-glob "results/run_*"
```

出力:
- `results/history_err_l2.png`
- `results/history_res_l2.png`
