# テスト報告

実行日: 2026-02-03  
コマンド: `julia --project -e 'using Pkg; Pkg.test()'`

## 目的
- ソルバー実行時の解析解誤差評価（相対 L2 誤差）

## 結果
- `solver_error`: **PASS**
  - 解析解との相対 L2 誤差が閾値以下
  - デフォルトは `N=16,32` で実行（`ADPOISSON_FULL_TEST=1` で `N=64` まで実行）
  - パラメータ: `M=4`, `epsilon=1e-6`, `bc_order=:high`
  - `dt` と `max_steps` は `scripts/main.jl` のデフォルト値を使用
  - `ADPOISSON_TEST_PLOT=1` の場合は可視化画像を `results/` に保存

## 補足
- `ADPOISSON_IMPL_TEST=1` の場合は実装確認テストも実行。
