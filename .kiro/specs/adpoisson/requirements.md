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
  - 方程式: $\frac{\partial u}{\partial t} = f-\nabla^2 u$ を想定（擬似時間）
  - 空間離散化: 差分法二次精度中心差分
  - 検証ケースでは $f=0$ を用いるが、一般性のために $f$ は残す
- [ ] **時間積分手法**
  - Taylor級数展開を用いた高次精度時間積分
  - Taylor展開次数は $M=10$（$m=0..M$）まで考える
  - 高次導関数の計算はADburgersの実装を参考にする
  - 3D Poisson式に対するTaylor級数展開の導出
  - ソース項fの扱いの検討
- [ ] **Taylor級数漸化式（3D Poisson, 擬似時間）**
  - 詳細は `/Users/Daily/Development/ADTM/ADPoisson漸化式.md` に準拠
  - Poisson: $\nabla^2 u = f$ を擬似時間で $u_t = f-\nabla^2 u$ として解く
  - 3D cell-centered、Julia 1-origin、ghost 1層を前提
  - Taylor係数: $(u_{i,j,k})_m=\frac{1}{m!}\left.\frac{\partial^m u_{i,j,k}}{\partial t^m}\right|_{t=t_0}$
  - 更新: $u_{i,j,k}(t_0+\Delta t)\approx \sum_{m=0}^{M}(u_{i,j,k})_m(\Delta t)^m$
  - 7点差分ラプラシアン:
    $$(L u)_{i,j,k}= \frac{u_{i+1,j,k}-2u_{i,j,k}+u_{i-1,j,k}}{\Delta x^2}+
    \frac{u_{i,j+1,k}-2u_{i,j,k}+u_{i,j-1,k}}{\Delta y^2}+
    \frac{u_{i,j,k+1}-2u_{i,j,k}+u_{i,j,k-1}}{\Delta z^2}$$
  - 漸化式（内点）:
    $$(u_{i,j,k})_{m+1}=\frac{1}{m+1}\Bigl((f_{i,j,k})_m-(L(u_m))_{i,j,k}\Bigr)$$
    $f$ が擬似時間一定なら $(f)_0=f,\ (f)_m=0\ (m\ge1)$ とし、$m\ge1$ では
    $$(u_{i,j,k})_{m+1}= -\frac{1}{m+1}(L(u_m))_{i,j,k}$$
  - 境界条件（ghost 係数式）:
    - Dirichlet $u=g$: $m=0$ で $u_{\text{ghost}}=2g-u_{\text{adj}}$、$m\ge1$ で $u_{\text{ghost}}=-u_{\text{adj}}$
    - Neumann $\partial u/\partial n=h$: $m=0$ で $u_{\text{ghost}}=u_{\text{adj}}\pm\Delta n\,h$、$m\ge1$ で $u_{\text{ghost}}=u_{\text{adj}}$
  - 6面の ghost 係数式（Julia 1-origin）:
    - Dirichlet（面データ: `g_xlo/g_xhi/g_ylo/g_yhi/g_zlo/g_zhi`）
      - x-min: $(u_{1,j,k})_0=2g_{xlo}[j-1,k-1]-(u_{2,j,k})_0$, $(u_{1,j,k})_m=-(u_{2,j,k})_m$
      - x-max: $(u_{N_x+2,j,k})_0=2g_{xhi}[j-1,k-1]-(u_{N_x+1,j,k})_0$, $(u_{N_x+2,j,k})_m=-(u_{N_x+1,j,k})_m$
      - y-min: $(u_{i,1,k})_0=2g_{ylo}[i-1,k-1]-(u_{i,2,k})_0$, $(u_{i,1,k})_m=-(u_{i,2,k})_m$
      - y-max: $(u_{i,N_y+2,k})_0=2g_{yhi}[i-1,k-1]-(u_{i,N_y+1,k})_0$, $(u_{i,N_y+2,k})_m=-(u_{i,N_y+1,k})_m$
      - z-min: $(u_{i,j,1})_0=2g_{zlo}[i-1,j-1]-(u_{i,j,2})_0$, $(u_{i,j,1})_m=-(u_{i,j,2})_m$
      - z-max: $(u_{i,j,N_z+2})_0=2g_{zhi}[i-1,j-1]-(u_{i,j,N_z+1})_0$, $(u_{i,j,N_z+2})_m=-(u_{i,j,N_z+1})_m$
    - Neumann（面データ: `h_xlo/h_xhi/h_ylo/h_yhi/h_zlo/h_zhi`）
      - x-min: $(u_{1,j,k})_0=(u_{2,j,k})_0+\Delta x\,h_{xlo}[j-1,k-1]$, $(u_{1,j,k})_m=(u_{2,j,k})_m$
      - x-max: $(u_{N_x+2,j,k})_0=(u_{N_x+1,j,k})_0-\Delta x\,h_{xhi}[j-1,k-1]$, $(u_{N_x+2,j,k})_m=(u_{N_x+1,j,k})_m$
      - y-min: $(u_{i,1,k})_0=(u_{i,2,k})_0+\Delta y\,h_{ylo}[i-1,k-1]$, $(u_{i,1,k})_m=(u_{i,2,k})_m$
      - y-max: $(u_{i,N_y+2,k})_0=(u_{i,N_y+1,k})_0-\Delta y\,h_{yhi}[i-1,k-1]$, $(u_{i,N_y+2,k})_m=(u_{i,N_y+1,k})_m$
      - z-min: $(u_{i,j,1})_0=(u_{i,j,2})_0+\Delta z\,h_{zlo}[i-1,j-1]$, $(u_{i,j,1})_m=(u_{i,j,2})_m$
      - z-max: $(u_{i,j,N_z+2})_0=(u_{i,j,N_z+1})_0-\Delta z\,h_{zhi}[i-1,j-1]$, $(u_{i,j,N_z+2})_m=(u_{i,j,N_z+1})_m$
- [ ] **検証機能**
  - 解析解の定義（Method of Manufactured Solutions等）
  - 数値解と解析解の比較は、相対 $L2$ 誤差 $\|u-u_{\text{exact}}\|_2/\|u_{\text{exact}}\|_2 \le 10^{-3}$ を満たすこと
- 計算格子
  - 直交格子で、分割数 $N_x, N_y, N_z$ を入力し、各辺長=1を割り、格子幅 $\Delta x,\Delta y,\Delta z$ を計算
  - 3D cell-centered 配置で ghost 1層を含む配列として保持し、境界条件はghost層に反映して適用する
  - 配列サイズは $(N_x+2, N_y+2, N_z+2)$
  - 内部点の添字は $i=2..N_x+1,\ j=2..N_y+1,\ k=2..N_z+1$
- **領域・境界条件・初期条件**
  - 計算領域: $[0,1]^3$
  - 境界条件（Dirichlet）:
    - $z=0$: $\phi(x,y,0)= \alpha \sin(\pi x)\sin(\pi y)$
    - $z=1$: $\phi(x,y,1)= \sin(\pi x)\sin(\pi y)$
    - その他境界: $\phi(x,y,z)=0$
  - $\alpha$ は難易度制御パラメータ（$z$ 方向の勾配の大きさに対応）。通常は $\alpha=1$
  - 初期条件: $u(x,y,z,0)=0$
  - 合否基準: 格子精細化で 2次精度が確認できること
  - 解析解（$f=0$ の検証ケース）:
    $$u(x,y,z)=\frac{\sin(\pi x)\sin(\pi y)}{\sinh(\sqrt{2}\pi)}\Bigl(\sinh(\sqrt{2}\pi z)+\alpha\,\sinh(\sqrt{2}\pi(1-z))\Bigr)$$

- Julia実装
  - パラメータはコマンドラインで指定
  - $N_x, N_y, N_z$, Taylor展開次数 $M$, $\Delta t$, $t_{\text{end}}$
  - 可視化は $y=0.5$ の断面（XZ面）をヒートマップ/コンターで表示し、解析解との誤差も可視化

## 非機能要件
- [ ] **実装言語**: Julia
- [ ] **参照実装**: データ構造や関数は`ADburgers` の設計思想・コード構造を踏襲する
- [ ] **拡張性**: 将来的な自動微分（AD）適用を考慮した設計とする

## 制約事項
- 開発環境: macOS
- 既存の `.kiro` エコシステムに従う
