# テスト報告

実行日: 2026-02-04  
コマンド: `ADPOISSON_FULL_TEST=1 ADPOISSON_TEST_PLOT=1 julia --project -e 'using Pkg; Pkg.test()'`

## 目的
- ソルバー実行時の解析解誤差評価（相対 L2 誤差）

## 結果
- `solver_error`: **PASS**
  - 解析解との相対 L2 誤差が閾値以下
  - デフォルトは `N=16,32` で実行（`ADPOISSON_FULL_TEST=1` で `N=64` まで実行）
  - パラメータ: `M=4`, `epsilon=1e-6`, `bc_order=:high`
  - `dt` と `max_steps` は `scripts/main.jl` のデフォルト値を基準とするが、
    **フルテスト安定化のため `dt <= 0.5/(3 n^2)` にクリップ**（`ADPOISSON_TEST_USE_SAFE_DT=1` がデフォルト）
  - `ADPOISSON_TEST_PLOT=1` の場合は可視化画像を `results/` に保存
  - 出力例には `exact_*.png`, `error_*.png`, `history_*.txt` が含まれる

## 補足
- `ADPOISSON_IMPL_TEST=1` の場合は実装確認テストも実行。
