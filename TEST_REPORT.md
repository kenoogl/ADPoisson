# テスト報告

実行日: 2026-02-03  
コマンド: `julia --project -e 'using Pkg; Pkg.test()'`

## 目的
- `laplacian!` の正当性（定数場で 0）
- `taylor_step!` の正当性（`m=0` で `(u)_1 = Lu - f`）
- 空間差分の収束次数（2次精度の確認）
- ソルバー実行時の解析解誤差評価（相対 L2 誤差）

## 結果
- `laplacian!`: **PASS**
- `taylor_step!`: **PASS**
- `convergence_order`: **PASS**（2次精度を確認）
- `solver_error`: **PASS**
  - 解析解との相対 L2 誤差が閾値以下
  - デフォルトは `N=16,32` で実行（`ADPOISSON_FULL_TEST=1` で `N=64` まで実行）
  - パラメータ: `M=4`, `dt=0.1/(3 n^2)`, `tend=0.5`, `epsilon=1e-6`, `bc_order=:high`

## 補足
- `convergence_order` は空間差分の収束次数のみを確認し、絶対誤差しきい値は `solver_error` に移動。
