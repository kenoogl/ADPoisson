# 要件定義書

## プロジェクト概要
本プロジェクト（ADPoisson）は、**三次元ポアソン方程式**を対象とし、**Taylor級数を用いた時間積分手法**の実装を行うものである。
実装にあたっては、先行する `/Users/Daily/Development/ADTM/ADburgers` プロジェクト（一次元Burgers方程式向け実装）を参考にする。
検証として、解析解を持つ問題を設定し、数値解との誤差評価を行う。

## 機能要件
### コア機能
- [ ] **三次元ポアソン方程式ソルバー**
  - 未知変数: スカラー場 $u(x, y, z, t)$ （※tは擬似時間で時間発展形式で扱う）
  - ポアソン方程式は通常時間独立であるが、本プロジェクトでは時間発展形式で扱う
  - 方程式: $\frac{\partial u}{\partial t} = \nabla^2 u - f$ を想定（擬似時間）
  - 空間離散化: 差分法二次精度中心差分
  - 検証ケースでは $f=0$ を用いるが、一般性のために $f$ は残す
- [ ] **時間積分手法**
  - Taylor級数展開を用いた高次精度時間積分
  - Taylor展開次数は $M=10$（$m=0..M$）まで考える
  - 高次導関数の計算はADburgersの実装を参考にする
  - 3D Poisson式に対するTaylor級数展開の導出
  - ソース項fの扱いの検討。コード実装ではfを考慮する。実際の問題ではf=0として実行
  - 擬似時間の停止条件は残差 $r=Lu-f$ を用いる（内点のみ評価）。相対残差 $\|r\|_2 / \max(\|r_0\|_2, 1) \le \epsilon$（$r_0$ は初期残差）。$\epsilon$ は指定可能とする。デフォルト値: 1e-10
  - 終了条件は「残差が閾値を満たす」または「反復回数が最大ステップ数に到達」のいずれか早い方
  - 最大ステップ数は指定可能とする。デフォルト値: 10000
  - 擬似時間ステップの履歴をファイルに出力する
    - ファイル名: `history_nx{nx}_ny{ny}_nz{nz}_M{M}_steps{steps}.txt`
    - 出力列: `step`, `err_l2`, `res_l2`（`res_l2` は初期残差で相対化した残差L2）
    - 出力先は指定可能（デフォルト: `results/`）。実行ごとに `run_YYYYMMDD_HHMMSS/` を作成し、その配下に保存する
    - 実行条件と結果を `run_config.toml` / `run_summary.toml` に記録する
      - `run_summary.toml` には最終反復時の `res_l2` を記録する
  - 擬似時間刻みの条件は拡散数で管理する（推奨条件として採用）。拡散数（$\nu=1$）は
    $$Fo=\Delta t\left(\frac{1}{\Delta x^2}+\frac{1}{\Delta y^2}+\frac{1}{\Delta z^2}\right)$$
    と定義し、明示的安定性の推奨条件として $Fo \le 1/2$ を満たす
  - コマンドラインでは $\Delta t$ での指定に加えて $Fo$ での指定を許可する（$Fo$ 指定がある場合は優先）
- [ ] **線形ソルバー（SOR/CG）**
  - $\nabla^2 u = f$ を直接解く SOR 法と CG 法を実装する
  - 例題は擬似時間法と同一の検証問題を用いる
  - 収束判定は相対残差 $\|r\|_2 / \max(\|r_0\|_2, 1) \le \epsilon$（$r=Lu-f$、内点のみ）で統一する
  - 線形系は内点のみの方程式として構成し、Dirichlet境界の寄与は RHS に取り込む
  - 収束履歴（`step`, `err_l2`, `res_l2`）を出力し、収束性を比較できるようにする
    - SOR: `history_sor_nx{nx}_ny{ny}_nz{nz}_steps{steps}.txt`
    - SSOR: `history_ssor_nx{nx}_ny{ny}_nz{nz}_steps{steps}.txt`
    - CG: `history_cg_nx{nx}_ny{ny}_nz{nz}_steps{steps}.txt`
  - 反復ループ内（内点更新）では `if` 分岐を使わず、事前に条件分岐を外へ出す
  - SOR の緩和係数 $\omega$ は 1.0 を既定値とする（後で変更可能な形で実装する）
  - CG の前処理はオプション指定とし、既定は **none**
    - 指定: `:none`（既定） / `:ssor`
    - SSOR 使用時の緩和係数 $\omega$ は 1.0 を使用（後で変更可能な形で実装する）
    - **RBSSOR の対称スイープ順**: 前進 R→B、後退 B→R、前進 B→R、後退 R→B（R=red, B=black）
    - SSOR 前処理では `order=:spec` を固定する（`order=:high` は使わない）
  - SSOR ソルバー（RBSSOR）を選択可能とする
  - CLI は `--solver` で `taylor/sor/ssor/cg/mg-uniform-taylor/mg-hierarchical-taylor/mg-correction-taylor` を指定し、`--cg-precond` は `--solver=cg` のときのみ有効
- [ ] **加速（マルチグリッド的アプローチ）**
  - Taylor 擬似時間法の残差履歴が「高周波が早く減衰し低周波が残る」挙動であるため、マルチグリッド的加速を検討する
  - **レベル1/2は実験済みで不採用**のため、実装対象外とする
  - **レベル1（疑似MG, 不採用）**:
    - 数ステップごとに残差 $r=Lu-f$ を計算し、強めの Taylor ステップ（大きめ $\Delta t$ / 低次）で補正を入れる簡易方式
    - 粗格子を明示的に持たず、低周波成分の緩和を狙う
    - 補正ステップでも Fo クリップを適用し、$\Delta t$ は既存の Fo 上限を超えない
    - 適用条件は Taylor 法と同一（高次境界を使う場合は `nx,ny,nz>=4`）
  - **レベル2（2-level MG, 不採用）**:
    - coarse grid を 1段（$N/2$）構築し、制限 $R$ / 補間 $P$ による補正を行う
    - スムーザは Taylor 擬似時間ステップ（$\nu_1,\nu_2$ 回）を用いる
    - 粗格子サイズは内点数で $(n_x,n_y,n_z)\mapsto(\lfloor n_x/2\rfloor,\lfloor n_y/2\rfloor,\lfloor n_z/2\rfloor)$ とし、
      いずれかが 4 未満になった時点で粗格子生成を終了する（2-level MG の適用条件は $n_x,n_y,n_z\ge 8$）
    - 制限 $R$ は 3D full-weighting、補間 $P$ は trilinear を基本とする
    - 粗格子境界値は fine からの転送ではなく、境界条件関数から直接評価する
  - **レベル3（V-cycle MG）**:
    - 多段の V-cycle（必要なら W-cycle）を構成し、最粗格子では直接解法または厳密解を用いる
    - 直接解法は密行列を組み、`A \\ b`（LU ベース）で解く
    - レベル依存で $\Delta t$ や Taylor 次数 $M$ を変えることを許容する
  - **階層 Taylor（Level-dependent Taylor）**:
    - レベル番号は **1=最細（fine）**、以降レベルが増えるほど粗格子とする
    - `mg_level_Ms` / `mg_level_dt_scales` は CLI から指定可能とする
      - `--mg-level-Ms 4,4,2,2`、`--mg-level-dt-scales 2.0,2.0,4.0,4.0`
      - `run_config.toml` に設定値を記録する
    - 配列長がレベル数未満の場合は **最後の値を繰り返して使用**する
    - `mg_level_dt_scales` 未指定時は **`mg_dt_scale` を全レベルに適用**する
    - レベル $\ell$ の Taylor 更新は、格子幅に対応した $\Delta t_\ell,\ M_\ell$ を用いる
    - $\Delta t_\ell = \Delta t \cdot s_\ell$（$s_\ell$ は `mg_level_dt_scales` で指定）
    - $M_\ell$ は `mg_level_Ms` で指定（未指定は全レベルで `mg_M`）
    - 漸化式は各レベルで同一:
      $$(u^{(\ell)})_{m+1}=\frac{1}{m+1}\Bigl((L_\ell u^{(\ell)}_m)-(f_\ell)_m\Bigr)$$
      $f$ が擬似時間一定なら $(f_\ell)_0=f_\ell,\ (f_\ell)_m=0\ (m\ge1)$
    - **誤差方程式** $L e = -r$（$r=Lu-f$）を解く場合は、粗格子で **境界条件はゼロ（Dirichlet）** とする
    - `mg_level_Ms` / `mg_level_dt_scales` は `run_config.toml` に保存する
    - `mg_M` のデフォルトは **4**、`mg_dt_scale` のデフォルトは **2.0**
    - 最粗格子は内点数がいずれか 4 未満になった時点で停止する
  - **補正方程式の Taylor 化（Correction-Taylor）**:
    - 目的: V-cycle の coarse 補正で、誤差方程式 $L e = -r$ を Taylor 擬似時間積分で解く
      - 適用範囲: V-cycle の **全 coarse レベル**に適用
    - 定義:
      - 残差は **全体で** $r=Lu-f$ とする
      - 補正方程式は $L e = -r$ とする
      - 補正擬似時間方程式を $e_t = L e + r$ とする（定常で $L e = -r$）
      - coarse へは $-r$（$f-Lu$）を制限した右辺を用いる
    - レベル $\ell$ の補正漸化式:
      - 初期値は $e^{(\ell)}_0=0$
      - $$(e^{(\ell)})_{m+1}=\frac{1}{m+1}\Bigl((L_\ell e^{(\ell)}_m)+\delta_{m0}r_\ell\Bigr)$$
      - $e^{(\ell)}_{\text{new}} = e^{(\ell)} + \sum_{m=1}^{M_\ell} (\Delta t_\ell)^m (e^{(\ell)})_m$
    - 境界条件:
      - 補正方程式では coarse grid の境界をゼロ Dirichlet とする
      - fine 側の解更新後は元問題の境界条件を再適用する
    - 適用オプション:
      - `--solver mg-correction-taylor` で Correction-Taylor を有効化する
      - `--mg-corr-M` / `--mg-corr-dt-scale` は `mg-correction-taylor` の場合のみ有効
      - 既定値: `mg_corr_M=2`, `mg_corr_dt_scale=1.0`, `mg_corr_steps=1`
      - `--mg-corr-nu1` / `--mg-corr-nu2` で **補正方程式（e）側**の pre/post を個別指定できる
        - 未指定時は `mg_corr_steps` を **pre/post 両方**に適用する
        - 最粗格子は `mg_corr_nu1 + mg_corr_nu2` 回の Taylor 反復で近似する
    - 収束判定:
      - 全体収束判定は既存どおり $\|r\|_2 / \max(\|r_0\|_2,1)$
      - 補正ステップ内は固定回数反復を既定とし、将来拡張で閾値停止を許容
    - 出力:
      - `run_config.toml` に `mg_correction`, `mg_corr_M`, `mg_corr_dt_scale`, `mg_corr_steps`,
        `mg_corr_nu1`, `mg_corr_nu2` を保存する
  - CG が適用可能であること（係数行列の対称性・正定性）を確認する
  - 出力は実行ごとの `run_YYYYMMDD_HHMMSS/` 配下に保存し、`run_config.toml` / `run_summary.toml` を記録する
  - MG の履歴ファイル命名（レベル3のみ）:
    - レベル3: `history_vcycle_nx{nx}_ny{ny}_nz{nz}_steps{steps}.txt`
- [ ] **Taylor級数漸化式（3D Poisson, 擬似時間）**
  - 詳細は `/Users/Daily/Development/ADTM/ADPoisson漸化式.md` に準拠
  - Poisson: $\nabla^2 u = f$ を擬似時間で $u_t = \nabla^2 u - f$ として解く
  - 3D cell-centered、Julia 1-origin、ghost 1層を前提
  - Taylor係数: $(u_{i,j,k})_m=\frac{1}{m!}\left.\frac{\partial^m u_{i,j,k}}{\partial t^m}\right|_{t=t_0}$
  - 更新: $u_{i,j,k}(t_0+\Delta t)\approx \sum_{m=0}^{M}(u_{i,j,k})_m(\Delta t)^m$
  - 7点差分ラプラシアン:
    $$(L u)_{i,j,k}= \frac{u_{i+1,j,k}-2u_{i,j,k}+u_{i-1,j,k}}{\Delta x^2}+
    \frac{u_{i,j+1,k}-2u_{i,j,k}+u_{i,j-1,k}}{\Delta y^2}+
    \frac{u_{i,j,k+1}-2u_{i,j,k}+u_{i,j,k-1}}{\Delta z^2}$$
  - 漸化式（内点）:
    $$(u_{i,j,k})_{m+1}=\frac{1}{m+1}\Bigl((L(u_m))_{i,j,k}-(f_{i,j,k})_m\Bigr)$$
    $f$ が擬似時間一定なら $(f)_0=f,\ (f)_m=0\ (m\ge1)$ とし、$m\ge1$ では
    $$(u_{i,j,k})_{m+1}= \frac{1}{m+1}(L(u_m))_{i,j,k}$$
    なお $f$ の時間依存は将来拡張とする
  - 境界条件（ghost 係数式）:
    - Dirichlet $u=g$: $m=0$ で $u_{\text{ghost}}=2g-u_{\text{adj}}$、$m\ge1$ で $u_{\text{ghost}}=-u_{\text{adj}}$
    - **高次境界（任意, m=0 のみ）**: 境界近傍の精度改善のため、$m=0$ のみ 3次外挿を許可する。
      - 一般式: $u_{\text{ghost}}=\frac{16}{5}g-3u_1+u_2-\frac{1}{5}u_3$
      - ここで $u_1,u_2,u_3$ は境界から内側方向に並ぶ 1,2,3 番目の内点値
      - 適用条件: 各方向で内点が 4 点以上（`nx,ny,nz >= 4`）
    - Neumann $\partial u/\partial n=h$: $m=0$ で $u_{\text{ghost}}=u_{\text{adj}}\pm\Delta n\,h$、$m\ge1$ で $u_{\text{ghost}}=u_{\text{adj}}$
  - 6面の ghost 係数式（Julia 1-origin）:
    - Dirichlet（面データ: `g_xlo/g_xhi/g_ylo/g_yhi/g_zlo/g_zhi`）
      - x-min: $(u_{1,j,k})_0=2g_{xlo}[j-1,k-1]-(u_{2,j,k})_0$, $(u_{1,j,k})_m=-(u_{2,j,k})_m$
      - x-max: $(u_{N_x+2,j,k})_0=2g_{xhi}[j-1,k-1]-(u_{N_x+1,j,k})_0$, $(u_{N_x+2,j,k})_m=-(u_{N_x+1,j,k})_m$
      - y-min: $(u_{i,1,k})_0=2g_{ylo}[i-1,k-1]-(u_{i,2,k})_0$, $(u_{i,1,k})_m=-(u_{i,2,k})_m$
      - y-max: $(u_{i,N_y+2,k})_0=2g_{yhi}[i-1,k-1]-(u_{i,N_y+1,k})_0$, $(u_{i,N_y+2,k})_m=-(u_{i,N_y+1,k})_m$
      - z-min: $(u_{i,j,1})_0=2g_{zlo}[i-1,j-1]-(u_{i,j,2})_0$, $(u_{i,j,1})_m=-(u_{i,j,2})_m$
      - z-max: $(u_{i,j,N_z+2})_0=2g_{zhi}[i-1,j-1]-(u_{i,j,N_z+1})_0$, $(u_{i,j,N_z+2})_m=-(u_{i,j,N_z+1})_m$
      - **高次境界（任意, m=0 のみ）**:
        - x-min: $(u_{1,j,k})_0=\frac{16}{5}g_{xlo}[j-1,k-1]-3(u_{2,j,k})_0+(u_{3,j,k})_0-\frac{1}{5}(u_{4,j,k})_0$
        - x-max: $(u_{N_x+2,j,k})_0=\frac{16}{5}g_{xhi}[j-1,k-1]-3(u_{N_x+1,j,k})_0+(u_{N_x,j,k})_0-\frac{1}{5}(u_{N_x-1,j,k})_0$
        - y-min: $(u_{i,1,k})_0=\frac{16}{5}g_{ylo}[i-1,k-1]-3(u_{i,2,k})_0+(u_{i,3,k})_0-\frac{1}{5}(u_{i,4,k})_0$
        - y-max: $(u_{i,N_y+2,k})_0=\frac{16}{5}g_{yhi}[i-1,k-1]-3(u_{i,N_y+1,k})_0+(u_{i,N_y,k})_0-\frac{1}{5}(u_{i,N_y-1,k})_0$
        - z-min: $(u_{i,j,1})_0=\frac{16}{5}g_{zlo}[i-1,j-1]-3(u_{i,j,2})_0+(u_{i,j,3})_0-\frac{1}{5}(u_{i,j,4})_0$
        - z-max: $(u_{i,j,N_z+2})_0=\frac{16}{5}g_{zhi}[i-1,j-1]-3(u_{i,j,N_z+1})_0+(u_{i,j,N_z})_0-\frac{1}{5}(u_{i,j,N_z-1})_0$
    - Neumann（面データ: `h_xlo/h_xhi/h_ylo/h_yhi/h_zlo/h_zhi`）
      - x-min: $(u_{1,j,k})_0=(u_{2,j,k})_0+\Delta x\,h_{xlo}[j-1,k-1]$, $(u_{1,j,k})_m=(u_{2,j,k})_m$
      - x-max: $(u_{N_x+2,j,k})_0=(u_{N_x+1,j,k})_0-\Delta x\,h_{xhi}[j-1,k-1]$, $(u_{N_x+2,j,k})_m=(u_{N_x+1,j,k})_m$
      - y-min: $(u_{i,1,k})_0=(u_{i,2,k})_0+\Delta y\,h_{ylo}[i-1,k-1]$, $(u_{i,1,k})_m=(u_{i,2,k})_m$
      - y-max: $(u_{i,N_y+2,k})_0=(u_{i,N_y+1,k})_0-\Delta y\,h_{yhi}[i-1,k-1]$, $(u_{i,N_y+2,k})_m=(u_{i,N_y+1,k})_m$
      - z-min: $(u_{i,j,1})_0=(u_{i,j,2})_0+\Delta z\,h_{zlo}[i-1,j-1]$, $(u_{i,j,1})_m=(u_{i,j,2})_m$
      - z-max: $(u_{i,j,N_z+2})_0=(u_{i,j,N_z+1})_0-\Delta z\,h_{zhi}[i-1,j-1]$, $(u_{i,j,N_z+2})_m=(u_{i,j,N_z+1})_m$
  - Neumann境界は将来拡張として扱い、本仕様での必須実装はDirichletのみ（Neumannの6面式は参考記載）
- [ ] **検証機能**
  - 解析解の定義（Method of Manufactured Solutions等）
  - 数値解と解析解の比較は、相対 $L2$ 誤差 $\|u-u_{\text{exact}}\|_2/\|u_{\text{exact}}\|_2 \le 10^{-3}$ を満たすこと
- 計算格子
  - 直交格子で、分割数 $N_x, N_y, N_z$ を入力し、各辺長=1を割り、格子幅 $\Delta x,\Delta y,\Delta z$ を計算
  - 3D cell-centered 配置で ghost 1層を含む配列として保持し、境界条件はghost層に反映して適用する
  - 配列サイズは $(N_x+2, N_y+2, N_z+2)$
  - 内部点の添字は $i=2..N_x+1,\ j=2..N_y+1,\ k=2..N_z+1$
  - cell-centeredの物理座標は $x_i=(i-1.5)\Delta x,\ y_j=(j-1.5)\Delta y,\ z_k=(k-1.5)\Delta z$
- **領域・境界条件・初期条件**
  - 計算領域: $[0,1]^3$
  - 境界条件（Dirichlet）:
    - $z=0$: $u(x,y,0)= \alpha \sin(\pi x)\sin(\pi y)$
    - $z=1$: $u(x,y,1)= \sin(\pi x)\sin(\pi y)$
    - その他境界: $u(x,y,z)=0$
  - 将来拡張としてNeumann境界も実装する
  - $\alpha$ は難易度制御パラメータ（$z$ 方向の勾配の大きさに対応）。通常は $\alpha=1$
  - 初期条件: $u(x,y,z,0)=0$
  - 評価格子サイズ: $N_x=N_y=N_z=16, 32, 64$ をコマンドラインパラメータで指定
  - 合否基準: 格子精細化で 2次精度が確認できること（$N_x=16/32/64$ で収束率2を確認）
  - 解析解（$f=0$ かつ時間一定の検証ケース）:
    $$u(x,y,z)=\frac{\sin(\pi x)\sin(\pi y)}{\sinh(\sqrt{2}\pi)}\Bigl(\sinh(\sqrt{2}\pi z)+\alpha\,\sinh(\sqrt{2}\pi(1-z))\Bigr)$$
  - 最大誤差も評価し、終了時に出力する（合否判定には使わず参考値として扱う）
  - ghost更新は面のみで、エッジ、コーナーは更新しない（評価・出力も内点のみ）
  - $\alpha$ の値はコマンドラインで指定。デフォルト1.0

## Julia実装
- パラメータはコマンドラインで指定
  - $N_x, N_y, N_z$（または `n` で等方格子指定）, Taylor展開次数 $M$, $\Delta t$ または $Fo$, 最大ステップ数, 境界条件次数（`spec`/`high`）, 出力ディレクトリ
  - ソルバー指定: `--solver taylor|sor|ssor|cg|mg-uniform-taylor|mg-hierarchical-taylor|mg-correction-taylor`（既定 `taylor`）
  - 前処理指定: `--cg-precond none|ssor`（`--solver=cg` の場合のみ有効）
  - 可視化は $y=0.5$ の断面（XZ面）をヒートマップ/コンターで表示し、解析解との誤差も可視化（出力先は指定可能）
  - 擬似時間ステップ履歴のファイル出力を行う

## 非機能要件
- [ ] **実装言語**: Julia
- [ ] **参照実装**: データ構造や関数は`ADburgers` の設計思想・コード構造を踏襲する
- [ ] **拡張性**: 将来的な自動微分（AD）適用を考慮した設計とする
- [ ] **メモリ効率**: 3Dのメモリ増加を考慮し、Taylor係数は全次数保持せず、係数バッファの使い回し（ping-pong）と逐次加算で実装する。必要配列は係数2枚＋和1枚を基本とし、必要なら解配列との共用で2枚まで削減する。
- [ ] **MGバッファ再利用**: V-cycle 内で `zeros/similar` による配列再確保を行わない。レベルごとの残差・補正・右辺・作業バッファは初回にワークスペースとして確保し、以降は **配列全域（境界/ゴースト領域を含む）を `fill!`** で初期化して再利用する。

## 制約事項
- 開発環境: macOS
- 既存の `.kiro` エコシステムに従う
