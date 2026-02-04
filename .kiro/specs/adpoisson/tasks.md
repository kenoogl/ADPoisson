# 実装計画

## 凡例
- (P): Phase内で並列実行可能

## Phase 1: 基盤構築とデータ構造
- [x] 1. プロジェクト構造のセットアップ
  - `src/ADPoisson.jl`, `src/types.jl`, `src/core.jl` 等の空ファイル作成
  - `ADburgers` を参考にエントリポイントを記述
- [x] 1a. Project.toml の作成
  - name = "ADPoisson"
  - deps: Plots, Printf, Statistics, Test, LinearAlgebra
- [x] 2. データ構造の実装 (`src/types.jl`)
  - `ProblemSpec` の定義 (Lx, Ly, Lz, alpha, source, dirichlet)
  - `BoundaryConditions` の定義 (g_xlo, g_xhi, g_ylo, g_yhi, g_zlo, g_zhi)
  - `TaylorBuffers3D` の定義 (bufA, bufB, acc)
  - `TaylorArrays3D` の定義（検証用途のみ、任意）
  - AD対応の型パラメトリック化（`T<:Real`）を適用
  - `SolverConfig` の定義 (グリッド数, M, dt, max_steps, epsilon)
  - `Solution` の定義 (物理座標, u配列)
  - _Requirements: コア機能-ソルバー, 計算格子_
  - _Design: BoundaryConditions, TaylorBuffers3D, TaylorArrays3D_

## Phase 2: コアアルゴリズム実装
> Phase 1 完了後に開始
- [x] 3. グリッド生成と初期化 (`src/core.jl`, `src/problems.jl`)
  - `make_grid(config)` 実装（内点のみ）
  - `initialize_solution(config, prob)` 実装
  - 物理座標配列 $x, y, z$ の計算
  - 初期条件 $u=0$ の設定
  - _Requirements: 計算格子, 初期条件_
- [x] 4. 解析解とソース項の実装 (`src/problems.jl`)
  - `exact_solution(x, y, z, alpha)` 実装
  - `source_term(x, y, z)` 実装 ($f=0$ だが)
  - `dirichlet_bc(x, y, z, alpha)` 実装
  - `boundary_from_prob(prob)` 実装（Dirichlet関数を6面に分解）
  - _Requirements: 解析解_
- [x] 5a. Laplacian 演算子の実装 (`src/core.jl`) (P)
  - `laplacian!(Lu, u, config)`
  - 7点差分の実装
  - _Requirements: Taylor級数漸化式_
- [x] 5b. Taylor係数計算ルーチンの実装 (`src/core.jl`) (P)
  - `taylor_step!(next, curr, f, m, config)`
  - `accumulate_taylor!(acc, coeff, dt_pow)`
  - `horner_update!` 実装（検証用途のみ）
  - 漸化式の実装
  - depends: [5a]
  - _Requirements: 時間積分手法, Taylor級数漸化式_
- [x] 5c. Ping-pong バッファリングの実装 (`src/core.jl`) (P)
  - `TaylorBuffers3D` を用いた逐次係数生成
  - depends: [5a, 5b]
  - _Requirements: メモリ効率_
- [x] 6. 境界条件の実装 (`src/boundary.jl`) (P)
  - `apply_bc!(u, bc, m, config)`
  - Dirichlet条件 (ghost cell更新) の実装
  - ghost更新は面のみ（エッジ・コーナーは更新しない）
  - 高次境界（任意, m=0 のみ、`nx,ny,nz>=4` が条件）
  - _Requirements: 領域・境界条件, Taylor級数漸化式_
- [x] 7. factory.jl の実装 (`src/factory.jl`)
  - `make_problem(config, alpha)` 等のヘルパーを実装
  - `boundary_from_prob(prob)` を利用して `BoundaryConditions` を生成
  - _Design: factory.jl_

## Phase 3: ソルバー統合と検証
- [x] 8. メインループの実装 (`src/core.jl`)
  - `solve(config, prob)` 実装
  - 時間ステップループ、収束判定 ($r=Lu-f$ 相対残差)
  - `compute_residual!` 実装（内点のみ）
  - 相対残差 $\|r\|_2 / \max(\|r_0\|_2, 1)$ を採用
  - 擬似時間ステップ履歴を `run_YYYYMMDD_HHMMSS/` 配下に保存（`history_nx{nx}_ny{ny}_nz{nz}_M{M}_steps{steps}.txt`）
  - 拡散数 $Fo$ の推奨条件チェック（警告）
  - depends: [5a, 5b, 5c, 6]
  - _Requirements: ソルバー, 時間積分手法_
- [x] 9. 可視化機能の実装 (`src/visualization.jl`)
  - `plot_slice(sol, ...)` 実装 (XZ面)
  - 解析解との差分（|u-u_exact|）の可視化
  - `results/` ディレクトリを作成して保存
  - 保存の命名規則: `exact_nx{nx}_ny{ny}_nz{nz}.png`, `error_nx{nx}_ny{ny}_nz{nz}_M{M}_steps{steps}.png`
  - depends: [8]
  - _Requirements: Julia実装-可視化_
- [x] 10. CLIと実行スクリプト (`scripts/main.jl`)
  - コマンドライン引数処理
  - `--nx --ny --nz --M --dt --Fo --max-steps --epsilon --alpha --bc-order --output-dir` に対応
  - `factory.jl` を利用して問題生成
  - depends: [7, 8]
  - _Requirements: Julia実装-パラメータ_
- [x] 11. 検証とテスト (`test/runtests.jl`)
  - $N=16, 32, 64$ でのL2誤差収束確認
  - 相対 $L2$ 誤差 $\|u-u_{\text{exact}}\|_2/\|u_{\text{exact}}\|_2 \le 10^{-3}$ を満たすこと
  - `test/core.jl`, `test/problems.jl` の構成で分割
  - `@testset "laplacian!"`, `@testset "taylor_step!"`, `@testset "convergence_order"` を追加
  - `laplacian!` テスト: 一様場 $u=c$ で $L=0$ を確認
  - `taylor_step!` テスト: $m=0$ で $(u)_1 = Lu - f$ を確認
  - 最大誤差も計算して出力（合否判定には使わない）
  - depends: [8]
  - 単体テスト（laplacian!, taylor_step!）は Phase 2 完了時点で実行可能
  - _Requirements: 検証機能_

## Phase 4: 線形ソルバー（SOR/CG）
> Phase 3 完了後に開始
- [ ] 12. SOR ソルバーの実装 (`src/sor.jl`)（内点のみ、Dirichlet境界の寄与は RHS に取り込み）
  - RB-SOR 反復
  - 相対残差 $\|r\|_2/\max(\|r_0\|_2,1)$ による収束判定
  - 収束履歴の出力（`step`, `err_l2`, `res_l2`）
  - 反復ループ内は `if` 分岐なし
  - depends: [5a, 6]
  - _Requirements: 線形ソルバー_
  - _Design: 線形ソルバー（SOR）_
- [ ] 13. CG ソルバーの実装 (`src/cg.jl`)（前処理: `:none` / `:ssor`）
  - 明示行列を組まず `laplacian!` による行列作用
  - SSOR 前処理（RBSSOR 4 スイープ、対称）
  - 相対残差の収束判定と履歴出力
  - 反復ループ内は `if` 分岐なし
  - depends: [5a, 6, 12]
  - _Requirements: 線形ソルバー_
  - _Design: 線形ソルバー（CG/SSOR）_
- [ ] 14. 収束性の比較
  - 擬似時間法（Taylor）と SOR/CG の収束履歴を同一指標で比較
  - 例題は共通の検証問題を使用
  - 比較結果を `run_YYYYMMDD_HHMMSS/` に保存
  - depends: [8, 12, 13]
  - _Requirements: 検証機能_
  - _Design: 比較スクリプト_
- [ ] 15. Phase 4 テスト (`test/solvers.jl`)
  - SOR: 小規模問題で既知解との一致を確認
  - CG: SORと同一問題で収束を確認
  - depends: [12, 13]
