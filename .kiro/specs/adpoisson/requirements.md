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
  - 擬似時間刻みの条件は拡散数で管理する（推奨条件として採用）。拡散数（$\nu=1$）は
    $$Fo=\Delta t\left(\frac{1}{\Delta x^2}+\frac{1}{\Delta y^2}+\frac{1}{\Delta z^2}\right)$$
    と定義し、明示的安定性の推奨条件として $Fo \le 1/2$ を満たす
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

- Julia実装
  - パラメータはコマンドラインで指定
  - $N_x, N_y, N_z$, Taylor展開次数 $M$, $\Delta t$, 最大ステップ数, 境界条件次数（`spec`/`high`）
  - 可視化は $y=0.5$ の断面（XZ面）をヒートマップ/コンターで表示し、解析解との誤差も可視化
  - 擬似時間ステップ履歴のファイル出力を行う

## 非機能要件
- [ ] **実装言語**: Julia
- [ ] **参照実装**: データ構造や関数は`ADburgers` の設計思想・コード構造を踏襲する
- [ ] **拡張性**: 将来的な自動微分（AD）適用を考慮した設計とする
- [ ] **メモリ効率**: 3Dのメモリ増加を考慮し、Taylor係数は全次数保持せず、係数バッファの使い回し（ping-pong）と逐次加算で実装する。必要配列は係数2枚＋和1枚を基本とし、必要なら解配列との共用で2枚まで削減する。

## 制約事項
- 開発環境: macOS
- 既存の `.kiro` エコシステムに従う
