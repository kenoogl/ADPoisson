# 設計要件書

---
**目的**: Implement a 3D Poisson solver using Taylor series pseudo-time stepping, following the ADburgers architecture.

**アプローチ**:
- `ADburgers.jl` のモジュール構成 (`types`, `core`, `problems` 等) を踏襲
- `requirements.md` で定義された数学的仕様を実装に落とし込む
- 将来的なAD対応を見据え、構造体を適切に定義する

---

## 概要
**目的**: 三次元ポアソン方程式 $\nabla^2 u = f$ を、擬似時間発展方程式 $\frac{\partial u}{\partial t} = \nabla^2 u - f$ とみなし、Taylor級数法を用いて定常解を求めるソルバーを提供する。
**ユーザー**: 数値流体力学の研究者、特に高次精度時間積分や自動微分に関心のある開発者。
**インパクト**: 3Dポアソン問題に対する高精度かつAD適用可能なソルバーの基盤となる。

### Goals
- Taylor級数展開（$M=10$次）による擬似時間積分ソルバーの実装
- 解析解（Method of Manufactured Solutionsまたは具体的境界値問題）との比較による精度検証
- コマンドラインからのパラメータ制御 ($N_x,N_y,N_z,M,\Delta t,\text{max\_steps},\epsilon,\alpha,\text{bc-order}$)

### Non-Goals
- 高速化（並列化・GPU化）は現時点での主目的ではない（将来拡張）
- 複雑な境界形状（直方体領域のみ対象）
- Neumann境界の実装（仕様上は定義されているが、今回の必須実装範囲外としてDirichletのみ扱う）

## アーキテクチャ

### アーキテクチャパターン
`ADburgers` と同様に、Juliaのモジュールシステムを利用した機能分割を行う。
`src/ADPoisson.jl` をエントリポイントとし、以下のサブコンポーネントに分割する。

### ディレクトリ構成
```
ADPoisson/
├── src/
│   ├── ADPoisson.jl       # メインモジュール
│   ├── types.jl           # データ構造定義 (ProblemSpec, SolverConfig, Solution 等)
│   ├── core.jl            # コアアルゴリズム (Taylor係数計算, 時間積分)
│   ├── boundary.jl        # 境界条件適用 (ghost cell 更新)
│   ├── problems.jl        # 問題設定 (解析解, ソース項定義)
│   ├── visualization.jl   # 可視化機能 (Plots.jl / Heatmap)
│   └── factory.jl         # インスタンス生成ヘルパー
├── scripts/
│   └── main.jl            # 実行スクリプト (CLI引数処理)
└── test/
    ├── runtests.jl        # テストエントリ
    ├── core.jl            # Taylor係数/更新テスト
    └── problems.jl        # 解析解/ソース項テスト
```

### モジュール依存関係
- `ADPoisson.jl` が全モジュールを `include` し、公開APIを定義
- `core.jl` → `types.jl`, `boundary.jl`, `problems.jl` に依存
- `boundary.jl` → `types.jl` に依存
- `problems.jl` → `types.jl` に依存
- `visualization.jl` → `types.jl`, `problems.jl` に依存
- `factory.jl` → `types.jl`, `problems.jl`, `boundary.jl` に依存

### 技術スタック
| Layer | Choice / Version | Role in Feature | Notes |
|-------|------------------|-----------------|-------|
| Language | Julia 1.10+ | Core Logic | |
| Visualization | Plots.jl | Result Visualization | 2D断面コンター図 |
| CLI | ARGS (Base) | Argument Parsing | 簡易的な引数処理 |

## コンポーネント詳細

### 1. データ構造 (`src/types.jl`)

#### `ProblemSpec`
問題設定（領域、境界条件パラメータ、ソース項）を保持。
```julia
struct ProblemSpec{T<:Real,Fsrc,Fbc}
    Lx::T; Ly::T; Lz::T  # デフォルトは 1.0（[0,1]^3）
    alpha::T # 境界条件用パラメータ
    source::Fsrc      # f(x,y,z)（時間一定を仮定）
    dirichlet::Fbc    # u|∂Ω = g(x,y,z,alpha)
end
```

#### `BoundaryConditions`
Dirichlet境界の6面関数を保持（Neumannは将来拡張）。
```julia
struct BoundaryConditions{Fxlo,Fxhi,Fylo,Fyhi,Fzlo,Fzhi}
    g_xlo::Fxlo; g_xhi::Fxhi
    g_ylo::Fylo; g_yhi::Fyhi
    g_zlo::Fzlo; g_zhi::Fzhi
end
```
`ProblemSpec.dirichlet` は境界条件の元関数（`(x,y,z,alpha)->u`）として受け取り、
`boundary_from_prob(prob)` で6面関数に分解して `BoundaryConditions` を生成する。
境界条件の正規形は `BoundaryConditions` とする。

**マッピング例（要件定義書の境界条件）**:
- $z=0$: $u=\alpha \sin(\pi x)\sin(\pi y)$ → `g_zlo(x,y)=alpha*sin(pi*x)*sin(pi*y)`
- $z=1$: $u=\sin(\pi x)\sin(\pi y)$ → `g_zhi(x,y)=sin(pi*x)*sin(pi*y)`
- 他4面: $u=0$ → `g_xlo=g_xhi=g_ylo=g_yhi=(x,y)->0`

#### `TaylorBuffers3D`
Taylor係数の逐次生成用バッファ（メモリ削減のため全次数保持しない）。
```julia
struct TaylorBuffers3D{T<:Real}
    bufA::Array{T,3}  # u_m
    bufB::Array{T,3}  # u_{m+1}
    acc::Array{T,3}   # Taylor和
end
```

#### `TaylorArrays3D`（任意、検証用）
係数全保持が必要な検証時にのみ使用（通常は非推奨）。
```julia
struct TaylorArrays3D{T<:Real}
    U::Array{T,4}  # (nx+2, ny+2, nz+2, M+1)
end
```

#### `SolverConfig`
数値計算パラメータを保持。
```julia
struct SolverConfig{T<:Real}
    nx::Int; ny::Int; nz::Int
    M::Int         # Taylor展開次数
    dt::T
    max_steps::Int
    epsilon::T # 収束判定閾値（相対残差、基準は初期残差）
end
```

#### `Solution`
計算結果（グリッド、解配列）を保持。
```julia
struct Solution{T<:Real}
    x::Vector{T}; y::Vector{T}; z::Vector{T}
    u::Array{T, 3} # ghost cell込み: (nx+2, ny+2, nz+2)
    t::T
    iter::Int
end
```

### 2. コアアルゴリズム (`src/core.jl`, `src/boundary.jl`)

#### `solve(config::SolverConfig, prob::ProblemSpec; bc_order=:spec, output_dir="results")`
1. グリッド初期化
2. 初期条件設定 ($u=0$)
3. 時間ループ
    - 残差計算: $r=Lu-f$（内点のみ、初期残差で相対化）
    - Taylor係数計算 (再帰的定義)
    - 解更新
4. 結果返却
   - まとめて `Fo`, `dt`, `err_l2`, `err_max`, `steps`, `runtime` を表示
5. 擬似時間ステップ履歴を `results/` に保存（`history_nx{nx}_ny{ny}_nz{nz}_M{M}_steps{steps}.txt`）
   - 出力列: `step`, `err_l2`, `res_l2`（`res_l2` は初期残差で相対化）

#### 関数シグネチャ（主要）
```julia
apply_bc!(u::Array{T,3}, bc::BoundaryConditions, m::Int, config::SolverConfig; order=:spec) where {T}
laplacian!(Lu::Array{T,3}, u::Array{T,3}, config::SolverConfig) where {T}
compute_residual!(r::Array{T,3}, u::Array{T,3}, f::Array{T,3}, config::SolverConfig) where {T}

make_grid(config::SolverConfig) -> (x::Vector, y::Vector, z::Vector)

taylor_step!(next::Array{T,3}, curr::Array{T,3}, f::Array{T,3}, m::Int, config::SolverConfig) where {T}
accumulate_taylor!(acc::Array{T,3}, coeff::Array{T,3}, dt_pow::T) where {T}

# 係数を保持する場合の評価（検証用途のみ）
horner_update!(u_new::Array{T,3}, coeffs::TaylorArrays3D{T}, dt::T, M::Int) where {T}
```

#### Taylor係数計算と更新
`requirements.md` の数式に基づき、メモリ効率を考慮した実装を行う。
- **メモリ管理**: 3D配列のメモリ消費を抑えるため、全次数分の `u_m` を保持しない。
  - `u_sum`: 解の更新用累積配列
  - `u_prev`: ひとつ前の次数の係数 ($u_m$)
  - `u_curr`: 計算中の次数の係数 ($u_{m+1}$)
  - これらを使い回し（ping-pong）、逐次的に `u_sum` に加算していくことで、必要な3D配列を最小限（2~3枚）に抑える。
- **Ghost Cell更新 (`boundary.jl`)**: 各次数 $m$ の計算直後に境界条件を適用する。
  - ghost更新は面のみ（エッジ・コーナーは更新しない）
  - **高次境界（任意）**: `order=:high` の場合、$m=0$ のみ 3次外挿を用いて精度改善。
    - 例（x-min）: $u_{1}=\frac{16}{5}g-3u_{2}+u_{3}-\frac{1}{5}u_{4}$
    - x-max, y-min/y-max, z-min/z-max も同様に内側3点を用いる
    - $m\ge1$ は仕様通り $u_{\text{ghost}}=-u_{\text{adj}}$
    - `nx,ny,nz>=3` を満たす場合にのみ使用可能
  - `solve(...; bc_order=:high)` で高次境界を有効化する
- **終了条件**: 相対残差 $\|r\|_2 / \max(\|r_0\|_2, 1) \le \epsilon$（$r_0$ は初期残差）または反復回数が最大ステップ数に到達した時点の早い方（最大ステップ数のデフォルトは 10000）。

#### メモリ効率（Taylor係数の保持）
ADburgersでは `TaylorArrays` で全次数の係数を保持するが、3Dではメモリが支配的になるため、**係数を逐次生成し、保存せずに和に加算する**方式を採用する。

- **方針**: 係数 $u_m$ と $u_{m+1}$ を **2枚のバッファで交互に保持**（ping-pong）し、Taylor和は別配列に逐次加算する。
- **必要メモリ**: 係数バッファ2枚 + 和の配列1枚（合計3枚、ghost込み）。`f` が必要なら別途1枚。
- **利点**: `TaylorArrays` のように $M+1$ 枚を保持せずに済むため、メモリ使用量が大幅に減る。

**実装イメージ**:
1. `u_next .= u` で和の初期化（$m=0$）
2. `bufA` を $u_m$、`bufB` を $u_{m+1}$ として `m=0..M-1` を反復
3. `bufB` に `laplacian!(bufB, bufA)` を計算し、漸化式で $u_{m+1}$ を作成
4. `apply_bc!(bufB, m+1, ...)` 後に `u_next .+= bufB * dt^m`
5. `bufA` と `bufB` を入れ替えて次の次数へ

この方式で **全次数の係数保持を避けつつ** Taylor和を構築できる。

#### Horner法による時間更新（検証用途のみ）
係数を全保持する場合（`TaylorArrays3D`、検証用途のみ）は、Horner法で
$u^{n+1} = (((u_M)\Delta t + u_{M-1})\Delta t + \cdots + u_0)$
として評価し、乗算回数を削減する。

さらにメモリ削減が必要な場合は、**解配列 `u` を係数バッファと共用**し、和の配列を `u_next` のみとする。
- 例: `u` を $u_m$ の一時格納に使い、`u_next` に和を蓄積
- 必要メモリは係数バッファ1枚 + 和の配列1枚（合計2枚）まで削減可能
- 収束判定や可視化で `u` の旧値が必要な場合は、最小限のコピーで対応する

### 3. 問題設定 (`src/problems.jl`)
解析解およびソース項 $f$ を提供する関数群。
- `exact_solution(x, y, z, alpha)`
- `source_term(x, y, z)` (今回は $f=0$ だが一般化のため用意)
- `dirichlet_bc(x, y, z, alpha)`（境界値の定義）

### 4. 可視化 (`src/visualization.jl`)
結果 `Solution` を受け取り、指定された断面 ($y=0.5$等) の分布図と誤差図を描画し、PNG保存する。
解析解の断面図も保存する（同一解像度なら同一になるため格子情報のみ付与）。
出力先は `output_dir` で指定し、デフォルトは `results/` とする。命名は以下を基本とする。
- `exact_nx{nx}_ny{ny}_nz{nz}.png`
- `error_nx{nx}_ny{ny}_nz{nz}_M{M}_steps{steps}.png`
- `history_nx{nx}_ny{ny}_nz{nz}_M{M}_steps{steps}.txt`

## データモデル
- **グリッド**: Cell-centered。インデックス $i=2 \dots N_x+1$ が内点。
- **座標**: $x_i=(i-1.5)\Delta x,\ y_j=(j-1.5)\Delta y,\ z_k=(k-1.5)\Delta z$
- **配列**: 3次元配列 `Array{T, 3}`（AD対応のため `T` をパラメトリックにする）。
  - `x,y,z` は内点のみの長さ `nx,ny,nz` を持つ（ghostは保持しない）
  - `make_grid` の戻り値も内点のみとする

## CLI引数（scripts/main.jl）
- 形式: `julia scripts/main.jl --nx=32 --ny=32 --nz=32 --M=10 --dt=1e-3 --Fo=0.3 --max-steps=10000 --epsilon=1e-10 --alpha=1.0 --output-dir results`
- 必須: `--nx,--ny,--nz`
- 任意: `--M,--dt,--Fo,--max-steps,--epsilon,--alpha,--output-dir`（`--Fo` があれば `--dt` より優先、デフォルトは requirements.md に準拠）

## エラーハンドリング
- パラメータチェック: $N_x, N_y, N_z > 0$, $M \ge 1$ 等。
- 拡散数 $Fo$ のチェック: 推奨条件 $Fo \le 0.5$ を超える場合に警告を出力。
- 拡散数の定義: $Fo=\Delta t\left(\frac{1}{\Delta x^2}+\frac{1}{\Delta y^2}+\frac{1}{\Delta z^2}\right)$
- 残差評価は内点のみ。ghost（面以外）やエッジ/コーナーは評価対象外とする。

## テスト戦略
- **単体テスト**: Taylor係数計算ルーチンの正当性（低次で手計算と比較）。
- **統合テスト**: $N=16, 32, 64$ で解析解と比較し、L2誤差の収束率が2次精度であることを確認する (`test/runtests.jl` または検証スクリプト)。
- **収束テスト**: 相対残差が $\epsilon$ に達すること、または最大ステップ数打ち切りで収束しない場合は警告する。
  - 例: `@testset "laplacian!"`, `@testset "taylor_step!"`, `@testset "convergence_order"`
  - 検証基準: 相対 $L2$ 誤差 $\|u-u_{\text{exact}}\|_2/\|u_{\text{exact}}\|_2 \le 10^{-3}$
  - 最大誤差も計算して出力（合否判定には使わない）
