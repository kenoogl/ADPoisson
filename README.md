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

**コマンドライン引数**
- `--nx`, `--ny`, `--nz`: 各方向の分割数（デフォルト: 16）
- `--M`: Taylor 展開次数（デフォルト: 10）
- `--dt`: 擬似時間刻み幅（デフォルト: 1e-4）
- `--tend`: 反復の上限時間（デフォルト: 1.0）
- `--epsilon`: 収束判定の相対残差閾値（デフォルト: 1e-10）
- `--alpha`: 境界条件パラメータ（デフォルト: 1.0）

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
