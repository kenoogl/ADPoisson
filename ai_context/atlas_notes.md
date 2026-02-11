# ADPoisson (Auto Differential Poisson)

## 本プロジェクトの論点

### 全体としての強み（査読者に刺さる点）

- Taylor 法を **時間積分ではなく誤差伝播制御として再解釈**
- MG スムーザーとしての **新しい自由度（次数・Δt）**
- Correction-Taylor による **V-cycle の再定義**
- Matrix-free / AD 由来の Taylor という **系譜的一貫性**

------

### 論文感が出ること

1. **残差のフーリエスペクトル**（Taylor が高周波をどう潰すか）**安定領域（Δt vs M）を Fourier モードで解析**
2. SSOR vs Taylor の **smoothing factor 比較**
3. 「1 Taylor step ≒ ν SOR steps」の経験則

これが入ると
 👉 **“Taylor is just another smoother” ではなく
 “Taylor is a controllable spectral filter”**
 という主張が成立します。





## 解法のコア

解く問題を$u_t=∇^2u-f$とおき、定常状態で$∇^2u-f=0$、つまり$∇^2u=f$の解が得られるとする。セルセンターの離散化で、ゴーストセル2層、juliaの1オリジン記法でコードを書く。3Dpoisson方程式のTaylor 級数を用いた擬似時間解法において、Taylor 係数と基本漸化式、およびディリクレ（二次精度と四次精度）、ノイマン境界条件を詳しく導出する。

ここでは **擬似時間方程式**
$$
u_t=\nabla^2u-f
$$
を、**3D・セルセンター配置・ゴースト2層・Julia 1-origin**で解くために、

1. **Taylor 係数の定義**
2. **基本漸化式（係数再帰）**
3. **境界条件のゴースト生成**

- Dirichlet：**2次精度**／**4次精度**
- Neumann：**（外向き法線微分）2次精度**（＋参考として4次も付記）

を、導出→そのまま実装できる形でまとめます。

------

### 0. 格子・インデクス（ghost 2層, Julia 1-origin）

内点セル数を ($N_x,N_y,N_z$)、格子幅を ($\Delta x,\Delta y,\Delta z$) とする。

ゴースト2層の配列サイズ：$
(N_x+4,;N_y+4,;N_z+4)
$

- **内点（物理セル中心）**：$
  i=3..N_x+2,\quad j=3..N_y+2,\quad k=3..N_z+2
  $
- **ゴースト層**：
  - x-min 側：$i=1,2$（外側2層）
  - x-max 側：$i=N_x+3,\,N_x+4$
  - y,z も同様に $j=1,2,N_y+3,N_y+4,\,k=1,2,N_z+3,N_z+4$

------

### 1. 空間離散（ラプラシアン）

以降、$(Lu)_{i,j,k}$ が離散ラプラシアンを表すとして議論します。

#### (A) 2次精度（7点, radius=1）

$$
(L^{(2)}u)_{i,j,k}=
\frac{u_{i+1,j,k}-2u_{i,j,k}+u_{i-1,j,k}}{\Delta x^2}+
\frac{u_{i,j+1,k}-2u_{i,j,k}+u_{i,j-1,k}}{\Delta y^2}+
\frac{u_{i,j,k+1}-2u_{i,j,k}+u_{i,j,k-1}}{\Delta z^2}
$$

#### (B) 4次精度（各方向5点, radius=2：クロス項なし 19点）

必要なら（ghost2 があるので可能）：
$$
u_{xx}\approx\frac{-u_{i+2}+16u_{i+1}-30u_i+16u_{i-1}-u_{i-2}}{12\Delta x^2}
$$
を各方向に足して $L^{(4)}$ を作れます（以前提示した形）。

> 以降の Taylor 漸化式は **$L$ を $L^{(2)}$でも $L^{(4)}$ でも同じ形**で成り立ちます。違いは $L(\cdot)$ の計算（ステンシル）だけです。

------

### 2. Taylor 係数と基本漸化式（擬似時間）

#### 2.1 Taylor 係数

展開点 $t=t_0$で
$$
(u_{i,j,k})_m := \frac{1}{m!}\left.\frac{\partial^m u_{i,j,k}}{\partial t^m}\right|_{t=t_0}
\qquad (m=0,1,2,\dots)
$$
と定義すると
$$
u_{i,j,k}(t_0+\Delta t) \approx \sum_{m=0}^{M}(u_{i,j,k})_m(\Delta t)^m.
$$

#### 2.2 基本関係

定義より
$$
(u_{i,j,k})_{m+1}=\frac{1}{m+1}\left((u_t)_{i,j,k}\right)_m.
$$
ここで PDE は
$$
u_t = Lu - f.
$$
$L$ は空間線形演算子なので
$$
(Lu)_m = L(u_m)
$$
が成立。

また通常 $f$ は擬似時間に依存しないので
$$
(f_{i,j,k})_0=f_{i,j,k},\qquad (f_{i,j,k})_m=0\ \ (m\ge1).
$$
よって **内点** $i=3..N_x+2$ 等で：

#### 一般形（全 $m\ge0$）

$$
\boxed{
(u_{i,j,k})_{m+1}
=\frac{1}{m+1}\Bigl( (L(u_m))_{i,j,k}-(f_{i,j,k})_m\Bigr)
}
$$

#### $f$ が擬似時間一定の簡約

- $m=0$：
  $$
  \boxed{
  (u_{i,j,k})_{1}=(L(u_0))_{i,j,k}-f_{i,j,k}
  }
  $$

- $m\ge1$：
  $$
  \boxed{
  (u_{i,j,k})_{m+1}=\frac{1}{m+1}(L(u_m))_{i,j,k}
  }
  $$

更新：
$$
\boxed{
u^{new}_{i,j,k}=\sum_{m=0}^{M}(u_{i,j,k})_m(\Delta t)^m
}
$$

------

### 3. 境界条件 → ghost 値（Taylor 係数ごとに埋める）

以下では **面上データ**を持つとします（内点に対応する 1-origin 配列）：

- Dirichlet：面値 (g)
  - `g_xlo[1:Ny,1:Nz]`, `g_xhi[1:Ny,1:Nz]`
  - `g_ylo[1:Nx,1:Nz]`, `g_yhi[1:Nx,1:Nz]`
  - `g_zlo[1:Nx,1:Ny]`, `g_zhi[1:Nx,1:Ny]`
- Neumann：外向き法線微分 (h=\partial u/\partial n)
  - 同形サイズで `h_*`

各次数 (m) の係数配列 ((u_m))（実装では `Um`）に対して **ghost を埋める**のが境界処理です。

------

#### 3.1 Dirichlet（(u=g)）— 2次精度（ghost2層の埋め方）

cell-centered で境界面はセル中心から半格子なので、面値を線形補間：
$$
u_{\text{face}}\approx\frac{u_{\text{adj}}+u_{\text{ghost1}}}{2}=g
\Rightarrow
u_{\text{ghost1}}=2g-u_{\text{adj}}.
$$
x-min でいうと

- adj：$i=3$（$x=+\Delta x/2$）
- ghost1：$i=2$（$x=-\Delta x/2$）
- ghost2：$i=1$（$x=-3\Delta x/2$）

2層目の ghost は、同じ考えで「一つ奥の内点 $i=4$（$x=3\Delta x/2$）」と対にして
$$
u(x=-3h/2)\approx 2g-u(x=+3h/2)
\Rightarrow u_{\text{ghost2}}=2g-u_{\text{adj2}}.
$$
（これは“反射拡張”で、2次の境界整合を壊しません）

#### Taylor 係数で（(g) が擬似時間一定）

- $m=0$：$(g)_0=g$
- $m\ge1$：$(g)_m=0$

x-min 面の例（$j=3..N_y+2,\, k=3..N_z+2$）：
$$
\boxed{(u_{2,j,k})_m=
\begin{cases}
2g_{x-}[j-2,k-2]-(u_{3,j,k})_0,& m=0\\
-(u_{3,j,k})_m,& m\ge1
\end{cases}}
$$

$$
\boxed{
(u_{1,j,k})_m=
\begin{cases}
2g_{x-}[j-2,k-2]-(u_{4,j,k})_0,& m=0\\
-(u_{4,j,k})_m,& m\ge1
\end{cases}}
$$
他の 5 面も同様に、法線方向に沿って「ghost1 は最外内点」「ghost2 はその次の内点」を使って `2g - u` を作ります。

------

#### 3.2 Dirichlet（(u=g)）— 4次精度（ghost2層）

4次精度で境界を整合させたい場合、**境界面 $x=0$ の値 $u(0)=g$** を含む補間多項式で、ghost 点の値を評価します。

x-min で、次の 5 点（quartic）を使います：

- 境界面：$x=0$ で $u=g$
- 内点セル中心：$x=\frac12h,\frac32h,\frac52h,\frac72h$（対応 $i=3,4,5,6$）

この 5 点で作る4次補間多項式 $p(x)$ に対して、ghost 位置

- ghost1：$x=-\frac12 h$（$i=2$）
- ghost2：$x=-\frac32 h$（$i=1$）
  を評価すると、**係数が一意に決まります**：

#### 4次 Dirichlet：x-min の ghost 値（$m=0$）

$$
\boxed{(u_{2,j,k})_0=
\frac{128}{35}\,g_{x-}[j-2,k-2]
-4(u_{3,j,k})_0+2(u_{4,j,k})_0-\frac{4}{5}(u_{5,j,k})_0+\frac{1}{7}(u_{6,j,k})_0
}
$$

$$
\boxed{
(u_{1,j,k})_0=
\frac{128}{7}\, g_{x-}[j-2,k-2]
-30(u_{3,j,k})_0 + 20(u_{4,j,k})_0 - 9(u_{5,j,k})_0+\frac{12}{7}(u_{6,j,k})_0
}
$$

#### 4次 Dirichlet：Taylor 高次係数（$m\ge1$）

$(g)_m=0$ なので、単に $g$ 項が落ちて
$$
\boxed{
(u_{2,j,k})_m=
-4(u_{3,j,k})_m+2(u_{4,j,k})_m-\frac{4}{5}(u_{5,j,k})_m+\frac{1}{7}(u_{6,j,k})_m
}
$$

$$
\boxed{
(u_{1,j,k})_m=
-30(u_{3,j,k})_m+20(u_{4,j,k})_m-9(u_{5,j,k})_m+\frac{12}{7}(u_{6,j,k})_m
}
$$

x-max 面（右側）も同じ係数で、添字を内側に向けて反転させればOKです（例：ghost1 は $i=N_x+3$、内点は $i=N_x+2,N_x+1,N_x,N_x-1$ を使う）。

> 重要：この4次 Dirichlet ghost は「ghost2層が必要な 4次ラプラシアン（radius=2）」と整合します。
> 内点が少ない場合（たとえば $N_x<4$）は 2次へフォールバックするのが実務的です。

------

### 3.3 Neumann（外向き法線微分 $\partial u/\partial n=h$）

ここでは **$h$ を“外向き法線微分”として与える**（符号はデータ側に吸収）とします。

#### 2次精度（ghost2層の埋め方）

cell-centered で境界面における法線微分を 2次精度で
$$
\left.\frac{\partial u}{\partial n}\right|_{\text{face}}
\approx \frac{u_{\text{ghost1}}-u_{\text{adj}}}{\Delta n}
$$

（ $\Delta n=\Delta x,\Delta y,\Delta z$）と置くと、
$$
\boxed{u_{\text{ghost1}}=u_{\text{adj}}+\Delta n\,h.}
$$
ghost2 は同じ傾きで線形外挿すると（距離が 3 倍）：
$$
u_{\text{ghost2}} \approx u(x=+3h/2)+3h\,u_x(0)
$$
に対応する形になり、
$$
\boxed{u_{\text{ghost2}}=u_{\text{adj2}}+3\,\Delta n\,h.}
$$
（x-min なら adj2 は $i=4$、ghost2 は $i=1$）

Taylor 係数で（$h$ 擬似時間一定）：

- $m=0$：$(h)_0=h$
- $m\ge1$：$(h)_m=0$

x-min 面の例：
$$
\boxed{
(u_{2,j,k})_m=
\begin{cases}
(u_{3,j,k})_0+\Delta x\,h_{x-}[j-2,k-2],& m=0\\
(u_{3,j,k})_m,& m\ge1
\end{cases}}
$$

$$
\boxed{
(u_{1,j,k})_m=
\begin{cases}
(u_{4,j,k})_0+3\Delta x\,h_{x-}[j-2,k-2],& m=0\\
(u_{4,j,k})_m,& m\ge1
\end{cases}}
$$

> 全境界 Neumann の場合は Poisson 解が定数分不定なので、平均値固定などゲージ条件が必要です。

#### （参考）4次 Neumann ghost（ghost2層）

4次ラプラシアンと整合させたい場合、境界での **法線微分条件**と内点4点（(i=3..6)）から quartic を作ると、x-min 側で
$$
(u_{2,j,k})_0=\frac{17u_3+9u_4-5u_5+u_6}{22}+\frac{12}{11}\Delta x\,h
$$
$$
(u_{1,j,k})_0=\frac{-135u_3+265u_4-135u_5+27u_6}{22}+\frac{60}{11}\Delta x\,h
$$
（ここで $u_q$ は $u_{q,j,k}$の略、$h=h_{x-}[j-2,k-2]$）

$m\ge1$ では $h$ 項が落ちて線形結合だけになります。

------

# 

- **空間ラプラシアンが2次（7点）**のままだと、全体の解の収束次数は基本2次に支配されがちです。
- ただし **境界の高次化（Dirichlet 4次）**は、境界誤差が支配的なケースや MG の粗格子誤差の抑制で効くことがあります。
- **空間も4次（radius=2）**にするなら、今の ghost2 はちょうど必要条件です（その場合 `laplacian4!` を追加して差し替えればOK）。



- dt の安定上限（拡散数）と $M$ の選び方



## Taylor 擬似時間法は「SOR的」な収束挙動を示す

今の問題は
$$
u_t = \nabla^2 u - f
$$
で、定常解が Poisson 方程式です。
 これは本質的に **拡散方程式の時間積分による緩和法**です。

スペクトルで見ると

離散ラプラシアン $L$ の固有値を $\lambda_k<0$ とすると、

- 誤差モード $e_k$ は

$$
\frac{d e_k}{dt} = \lambda_k e_k
$$

- 高周波モード（$|\lambda_k|$ が大きい）は **急速に減衰**
- 低周波モード（$|\lambda_k|$ が小さい）は **非常に遅い**

Taylor 法は高次精度の時間積分ですが、**本質は陽的時間積分**なので、

> 高周波は一瞬で消える
>  低周波だけがいつまでも残る

という挙動になります。

👉 これは

- Jacobi
- SOR
- 擬似時間緩和
   すべてに共通する **smoother 的挙動**です。

Taylor 擬似時間法は：

- **極めて優秀なスムーザ**
- 高周波減衰率が高い（特に高次 Taylor）

つまり：

> **Taylor 法 = マルチグリッドの理想的スムーザ**

という位置づけになります。



## V-cycle MG

### (1) 高次 Taylor = 周波数選択性が良い

Taylor 法は次数 $M$ を上げると
$$
e^{\Delta t L} \approx \sum_{m=0}^M \frac{(\Delta t L)^m}{m!}
$$
に近づくので、

- 高周波（大きな $|\lambda|$）ほど **強く抑制**
- 低周波はあまり動かない

👉 **理想的な smoother 特性**

### (2) dt と M を「スムーザ設計パラメータ」にできる

MG 的には：

- fine grid：
  - 小さめ dt
  - 高次 Taylor（M=4〜6）
- coarse grid：
  - 大きめ dt
  - 低次 Taylor（M=2〜3）

という **レベル依存設計**が可能です。

### (3) 行列を一切組まない

- 完全に matrix-free
- 幾何 MG と相性抜群
- GPU / 並列にも強い

これは Krylov + MG より実装が軽い場合もあります。

### (4) SOR との違い（重要）

SOR:

- パラメータ ω 調整がシビア
- 高周波と低周波の分離が弱い

Taylor:

- **次数 M が「自動スペクトルフィルタ」**
- ω チューニング不要
- MG にした瞬間に「最適級」

👉 **SOR + MG より、Taylor + MG の方が設計自由度が高い**



## Taylor 展開そのものに MG 的な階層効果を組み込むアイデア

------

前述は「MGのスムーサーとして Taylor を使う」ですが、**Taylor 展開そのものに MG 的な階層効果を組み込む**設計にできます。方向性を3つに整理します。

### 1. 階層 Taylor（Level-dependent Taylor）

各レベルで **格子幅に応じた Taylor 展開次数と Δt を変える**。
目的は「細かい格子は高次＋小Δt、粗い格子は低次＋大Δt」。

- fine: M_f, dt_f
- coarse: M_c, dt_c
- dt_c = dt_f * 4（格子幅2倍なら安定条件は4倍まで許容）

**MG の効果が Taylor の系列近似自体に埋め込まれる**ため、
V-cycle の補正が「Taylor の階層化」と一致します。

------

### 2. 補正方程式の Taylor 化（Correction-Taylor）

通常 MG は L e = r を解くが、
これを Taylor で **疑似時間の補正方程式として直接解く**。

```
e_t = L e - r 
```

を各レベルで Taylor 展開して **補正量 e を直接更新**。
このとき coarse grid は **e だけを進める**ので、
「補正の Taylor」が MG の本質（誤差伝播）に一致。

------

### 3. スペクトル整形 Taylor（Preconditioned Taylor）

Taylor 展開に **coarse grid の逆作用を事前に混ぜる**：

```
u_{new} = u + P (L_c^{-1}) R (f - L u) + Taylor_smoother(...) 
```

つまり Taylor 更新式の中に **粗格子補正項を組み込む**。
これにより V-cycle を明示的に呼ばずに
「一回の Taylor 更新が MG の V-cycle に相当」になる。

------



# 階層 Taylor（Level-dependent Taylor）

以下に、あなたの意図（**fine: 高次+小Δt / coarse: 低次+大Δt、dt は h²スケールで増やす**）をそのまま使って、**階層 Taylor（Level-dependent Taylor）**を **MG の V-cycle と整合**する形でアルゴリズムとして示します。

ここでは問題は一貫して
$$
u_t = L_\ell u - f_\ell \qquad (L_\ell \approx \nabla^2)
$$
で、各レベルで Taylor 擬似時間ステップを **スムーザー**として使います（＝「Taylor の系列近似に MG の効果を埋め込む」）。

------

### 0) レベル依存パラメータ設計（あなたの前提を明示）

レベル $\ell=1$ を最細（fine）、$\ell=L$ を最粗（coarse）とする。

- 格子幅：$h_\ell = h_1\cdot 2^{\ell-1}$

- 安定条件（陽的拡散）に合わせて
  $$
  \Delta t_\ell = \Delta t_1 \cdot 4^{\ell-1}
  $$
  （= あなたの「dt_c = dt_f * 4」を一般化）

- Taylor 次数：例として
  $$
  M_\ell = \text{Ms}[\ell]
  $$
  （細かいほど大きく、粗いほど小さく）

------

### 1) 核心：レベル依存 Taylor スムーザー

**Taylor 擬似時間 1 step（次数 (M_\ell)）**を
$$
u \leftarrow \mathcal{T}_\ell(u; f_\ell,\Delta t_\ell,M_\ell)
$$
と書く。中身はあなたが実装済みの係数再帰：

- $u_0=u$
- $u_1=L_\ell u_0 - f_\ell$
- $u_{m+1}=\frac{1}{m+1}L_\ell u_m\quad (m\ge1)$
- $u^{new}=\sum_{m=0}^{M_\ell}(\Delta t_\ell)^m u_m$

これを **ν 回**まわすのが smoother（pre/post）です。

------

### 2) 階層 Taylor を含む V-cycle アルゴリズム（Level-dependent Taylor MG）

MG の構造は普通の V-cycle ですが、違いは

- 各レベルの smoother が **Taylor $M_\ell, dt_\ell$**
- coarse solve も同じ Taylor（または direct）

という点です。

------

#### Algorithm: V-cycle with Level-dependent Taylor smoother

**入力**

- レベル集合 $\ell=1..L$
- 作用素 $L_\ell$
- 右辺 $f_\ell$
- 制限 $R_\ell^{\ell+1}$、補間 $P_{\ell+1}^{\ell}$
- Taylor パラメータ $(M_\ell,\Delta t_\ell)$
- pre/post 回数 $\nu_1,\nu_2$、最粗反復回数 $\nu_c$

**出力**

- 更新された $u_1$

------

#### `Vcycle` ($ℓ, u_ℓ, f_ℓ$) （再帰）

1. **Pre-smoothing (Level-dependent Taylor)**
   $$
   u_\ell \leftarrow \mathcal{T}_\ell^{(\nu_1)}(u_\ell; f_\ell, \Delta t_\ell, M_\ell)
   $$

2. **Residual**
   $$
   r_\ell \leftarrow f_\ell - L_\ell u_\ell
   $$

3. **Coarse-grid residual**
   $$
   r_{\ell+1} \leftarrow R_\ell^{\ell+1}, r_\ell
   $$

4. **Coarse-grid correction problem**

   - 初期値 $e_{\ell+1}=0$

   - もし $\ell+1=L$（最粗）なら：
     $$
     e_L \leftarrow \mathcal{T}_L^{(\nu_c)}(e_L;\ r_L,\ \Delta t_L,\ M_L)
     $$

   - ここで **補正方程式**は
     $$
     e_t = L_L e - r_L
     $$
     を Taylor で回していると解釈できます（Correction-Taylor）。

   - そうでなければ：
     $$
     e_{\ell+1}\leftarrow \text{Vcycle}(\ell+1,\ e_{\ell+1},\ r_{\ell+1})
     $$

5. **Prolongation and correction**

   
   $$
   e_\ell \leftarrow P_{\ell+1}^{\ell}, e_{\ell+1}
   $$
   $$
   u_\ell \leftarrow u_\ell + e_\ell
   $$

6. **Post-smoothing (Level-dependent Taylor)**
   $$
   u_\ell \leftarrow \mathcal{T}_\ell^{(\nu_2)}(u_\ell; f_\ell, \Delta t_\ell, M_\ell)
   $$

7. return (u_\ell)

------

## 3) “階層 Taylor が MG を埋め込む”の解釈

- fine grid の Taylor（高次・小dt）は **高周波誤差を激減**（スムーザー）
- coarse grid では同じ誤差が相対的に高周波化するので、**低次・大dt**でも効率よく減衰
- したがって
  - V-cycle の「周波数分離」と
  - Taylor の「スペクトルフィルタ」
    が同じ役割分担になり、**MG の効果が Taylor パラメータ列 ${M_\ell,\Delta t_\ell}$に現れる**

------

## 4) 実装時の最小仕様（重要）

- **dt のクリップ**：$\Delta t_\ell$ は理論上 $h_\ell^2$ スケールで増えるが、Taylor 次数 $M_\ell$ と合わせて安定域が変わるので、実装では

  
  $$
  \Delta t_\ell \leftarrow \min(\Delta t_\ell, \Delta t_{\ell,\max})
  $$
  を入れるのが安全。

- **coarse の BC は誤差の同次 BC**（Dirichlet→0、Neumann→0）

- **operator は各レベルで一致した離散化**（Galerkin か 幾何）

- **係数バッファの使い回し**  
  `bufA, bufB, acc` を levelごとに再利用

- **m=0 の係数は residual reuse**  
  既存の `taylor_series_update_reuse!` をそのまま用いる 



### 期待される効果

#### 1. **高周波誤差の除去が加速（fine レベル）**

fine では `M` を大きくし、`dt` を小さくすることで  
Taylor の高周波減衰が強くなります。  
→ **V-cycle 前の残差がよりクリーンになる**

#### 2. **低周波補正が強化（coarse レベル）**

coarse では `dt` を大きく、`M` を小さくしても  
安定性が保たれます（格子が粗いほど安定条件が緩い）。  
→ **粗格子補正が効きやすくなる**

#### 3. **収束の停滞を減らせる**

現在見えている  
「残差は下がるが誤差が停滞」  
という挙動は、fine と coarse の smoothing のバランスが崩れている兆候です。  
→ レベル依存でバランスを調整できる

#### 4. **計算コストの最適化**

- fine だけ高次（高コスト）  
- coarse は低次（低コスト）  
  → **全体の計算量を抑えつつ収束率を上げる**



##  Taylor 係数・漸化式・更新式

 階層 Taylor（Level-dependent Taylor)を
「各レベルで独立した擬似時間 ODE を Taylor 展開で解く」という立場から、
レベル番号を明示した Taylor 係数・漸化式・更新式を完全に書き下します。

問題設定・符号は一貫して **あなたの指定どおり**
$$
u_t = \nabla^2 u - f
\quad\Longleftrightarrow\quad
u_t = L_\ell u - f_\ell
$$


------

### 1. レベル付き問題設定

レベル $\ell = 1,\dots,L$

- 格子幅
  $$
  h_\ell = h_1,2^{\ell-1}
  $$

- 離散ラプラシアン
  $$
  L_\ell \approx \nabla^2 \quad (\text{2次 or 4次、ghost 2層})
  $$

- 擬似時間方程式
  $$
  \boxed{
  \partial_t u^{(\ell)} = L_\ell u^{(\ell)} - f_\ell
  }
  $$

------

### 2. Taylor 係数の定義（各レベル共通）

各レベル (\ell) で、展開点 $t=t_0$ において
$$
\boxed{
(u^{(\ell)}_{i,j,k})_m
:= \frac{1}{m!}
\left.\frac{\partial^m u^{(\ell)}_{i,j,k}}{\partial t^m}\right|_{t=t_0}
}
\qquad m=0,1,2,\dots
$$
これにより
$$
u^{(\ell)}_{i,j,k}(t_0+\Delta t_\ell)
\;\approx\;
\sum_{m=0}^{M_\ell}
(\Delta t_\ell)^m (u^{(\ell)}_{i,j,k})_m
$$

------

### 3. Taylor 係数の基本関係（一般公式）

定義から必ず成り立つ恒等式：
$$
\boxed{ (u^{(\ell)})_{m+1}

\frac{1}{m+1}
\left( \partial_t u^{(\ell)} \right)_m
}
$$
ここに PDE を代入します。

------

### 4. レベル付き Taylor 漸化式の導出

擬似時間 PDE：
$$
\partial_t u^{(\ell)} = L_\ell u^{(\ell)} - f_\ell
$$

### (1) 線形演算子の性質

$L_\ell$ は空間線形演算子なので
$$
(L_\ell u^{(\ell)})_m = L_\ell (u^{(\ell)}_m)
$$

### (2) RHS (f_\ell) の Taylor 係数

通常 $f_\ell$ は擬似時間に依存しない：
$$
(f_\ell)_0 = f_\ell,
\qquad
(f_\ell)_m = 0 \quad (m\ge1)
$$

------

### 4.1 一般形（全 (m\ge0)）

$$
\boxed{ (u^{(\ell)})_{m+1}

\frac{1}{m+1} \Bigl( L_\ell (u^{(\ell)}_m)

(f_\ell)_m
\Bigr)
}
$$

------

### 4.2 実際に使う形（(f) 擬似時間一定）

#### 初期係数

$$
\boxed{
(u^{(\ell)})_0 = u^{(\ell)}
}
$$

#### 第1係数

$$
\boxed{
(u^{(\ell)})_1 = L_\ell (u^{(\ell)}_0) - f_\ell
}
$$

#### 高次係数（(m\ge1)）

$$
\boxed{ (u^{(\ell)})_{m+1}

\frac{1}{m+1}\,
L_\ell (u^{(\ell)}_m)
}
$$

👉 **ここが重要**

- 各レベルで **同一構造**

- 違うのは
  $$
  L_\ell,\quad \Delta t_\ell,\quad M_\ell
  $$
  だけ

------

### 5. レベル付き Taylor 更新式

各レベル $\ell$ における **1 Taylor step**：
$$
\boxed{ u^{(\ell),\,new}
\sum_{m=0}^{M_\ell}
(\Delta t_\ell)^m
(u^{(\ell)})_m
}
$$
これを **(\nu) 回繰り返す**と smoother：
$$
u^{(\ell)}
\;\leftarrow\;
\bigl(\mathcal T_\ell\bigr)^{\nu}
\bigl(u^{(\ell)}\bigr)
$$

------

### 6. 階層 Taylor の具体像（fine / coarse の違い）

### fine level（(\ell=1)）

$$
\begin{aligned}
(u^{(1)})_1 &= L_1 u^{(1)} - f_1 \\
(u^{(1)})_2 &= \tfrac12 L_1 (u^{(1)}_1) \\
(u^{(1)})_3 &= \tfrac13 L_1 (u^{(1)}_2) \\
&\vdots \\
u^{(1),new} &= u^{(1)} + \Delta t_1 u^{(1)}_1 + \Delta t_1^2 u^{(1)}_2 + \cdots
\end{aligned}
$$

- $M_1$：大
- $\Delta t_1$：小
- **高周波を強烈に減衰（スムーザー）**

------

### coarse level（(\ell=2)）

$$
\begin{aligned}
(u^{(2)})_1 &= L_2 u^{(2)} - f_2 \\
(u^{(2)})_2 &= \tfrac12 L_2 (u^{(2)}_1) \\
&\vdots \\
u^{(2),new} &= u^{(2)} + \Delta t_2 u^{(2)}_1 + \cdots
\end{aligned}
$$

- $h_2 = 2h_1$
- $\Delta t_2 \approx 4\Delta t_1$
- $M_2 < M_1$

👉 fine で残った **低周波誤差が coarse では高周波化**し、
少ない Taylor 次数でもよく減衰。

------

### 7. Correction-Taylor（補正方程式）との完全一致

coarse で解く補正方程式
$$
L_\ell e_\ell = r_\ell
$$
を擬似時間化：
$$
\boxed{
\partial_t e_\ell = L_\ell e_\ell - r_\ell
}
$$
に対しても **全く同じ Taylor 漸化式**：
$$
\boxed{ (e^{(\ell)})_{m+1}

\frac{1}{m+1} \Bigl( L_\ell (e^{(\ell)}_m)

(r_\ell)_m
\Bigr)
}
$$

- 初期値：$(e^{(\ell)})_0=0$
- $(r_\ell)_0=r_\ell,\ (r_\ell)_m=0\ (m\ge1)$

👉 **解と補正が同一アルゴリズム**
👉 MG が「Taylor の階層化」と一致する理由

------

### 8. アルゴリズム上の最重要ポイント（整理）

1. **Taylor 漸化式の形は全レベルで同一**
   1. レベル差は
      $(L_\ell,\ \Delta t_\ell,\ M_\ell)$
      のみ
2. coarse ほど
   - $\Delta t_\ell \uparrow$
   - $M_\ell \downarrow$
3. MG の
   - smoother
   - correction
     を **Taylor 1本で統一**

------

### まとめ（核心）

- 階層 Taylor は
  **「レベルごとに異なる Taylor 多項式で同一の拡散 ODE を解く」**
- MG の V-cycle は
  **Taylor 展開パラメータ列 ${M_\ell,\Delta t_\ell}$ の切り替え**として再解釈できる
- 数式的にも、実装的にも **完全に首尾一貫**

------



# Correction-Taylor

Correction-Taylor では、各レベル $\ell$ で **補正方程式**
$$
L_\ell e_\ell = r_\ell
$$
を、擬似時間 ODE
$$
\boxed{ \partial_t e_\ell = L_\ell e_\ell - r_\ell }
$$
として解き、定常状態 $\partial_t e_\ell=0$ で $L_\ell e_\ell=r_\ell$ を満たすようにします。

以下、**レベル番号付き**で、**Taylor 係数の漸化式**を明示します。

------

## 1) Taylor 係数の定義（補正 (e_\ell)）

展開点 $t=t_0$ で
$$
\boxed{
(e^{(\ell)}_{i,j,k})_m
:=\frac{1}{m!}\left.\frac{\partial^m e^{(\ell)}_{i,j,k}}{\partial t^m}\right|_{t=t_0}
}
\qquad m=0,1,2,\dots
$$
同様に（必要なら）残差 $r_\ell$ の係数を
$$
(r^{(\ell)}_{i,j,k})_m
:=\frac{1}{m!}\left.\frac{\partial^m r^{(\ell)}_{i,j,k}}{\partial t^m}\right|_{t=t_0}
$$
と置けますが、通常 $r_\ell$ は “その時点の残差” を固定して coarse で解くので **擬似時間一定**として扱います。

------

## 2) 基本恒等式

定義より
$$
\boxed{
(e^{(\ell)})_{m+1}
=\frac{1}{m+1}\left( \partial_t e^{(\ell)}\right)_m
}
$$

------

## 3) `Correction-Taylor` の漸化式（一般形）

ODE：
$$
\partial_t e^{(\ell)} = L_\ell e^{(\ell)} - r_\ell
$$
線形性から
$$
(L_\ell e^{(\ell)})_m = L_\ell (e^{(\ell)}_m)
$$
よって
$$
\boxed{
(e^{(\ell)})_{m+1}
=\frac{1}{m+1}\Bigl( L_\ell (e^{(\ell)}_m) - (r_\ell)_m \Bigr)
}
\qquad (m=0,1,2,\dots)
$$

------

## 4) 残差 (r_\ell) を擬似時間一定とする場合（実用形）

coarse-grid で補正を作るときは $r_\ell$ を固定するので
$$
(r_\ell)_0=r_\ell,\qquad (r_\ell)_m=0\ (m\ge1)
$$
したがって：

- 初期値（典型）
  $$
  \boxed{ (e^{(\ell)})_0 = 0 }
  $$

- 1次係数
  $$
  \boxed{
  (e^{(\ell)})_1 = L_\ell (e^{(\ell)}_0) - r_\ell = -r_\ell
  }
  $$
  （初期値を 0 にした場合）

- 高次係数（$m\ge 1$）
  $$
  \boxed{
  (e^{(\ell)})_{m+1}
  =\frac{1}{m+1}\,L_\ell\,(e^{(\ell)}_m)
  }
  $$

------

## 5) Taylor 更新式（補正 (e_\ell) の 1 step）

レベル $\ell$ の Taylor 次数 $M_\ell$、擬似時間ステップ $\Delta t_\ell$ で
$$
\boxed{ e_\ell^{new}

\sum_{m=0}^{M_\ell}
(\Delta t_\ell)^m (e^{(\ell)})_m
}
$$
（通常 $(e^{(\ell)})_0=0$ なので $m=1..M_\ell$ の和だけでよい）

これを $\nu$ 回回して coarse 補正を作り、
$$
u_\ell \leftarrow u_\ell + P_{\ell+1}^{\ell} e_{\ell+1}
$$
に使います。

------

## 6) “MG の誤差伝播”との一致（短い補足）

$u^{*}$を真の離散解とし、誤差 $E_\ell=u^*_\ell-u_\ell$ とすると、
$$
r_\ell=f_\ell-L_\ell u_\ell = L_\ell(u^*_\ell-u_\ell)=L_\ell E_\ell
$$
なので補正方程式 $L_\ell e_\ell=r_\ell$ は **誤差そのもの** $e_\ell=E_\ell$ を求める式です。
従って `Correction-Taylor` は「粗格子で誤差を直接緩和する」MG の本質と一致します。

------

必要なら、この漸化式を **ghost2層・Dirichlet/Neumann（誤差は同次BC）**で実装する際の注意点（特に Neumann-only のゲージ固定）も、コード断片込みで整理します。



## Correction-Taylorの境界条件

`Correction-Taylor` は **「誤差方程式を擬似時間で解く」**ので、
境界条件は **元の問題の BC をそのまま使うのではなく、「誤差の BC」**になります。
ここを曖昧にすると MG が壊れます。

以下、**原問題 → 補正方程式 → Correction-Taylor** の流れで整理します。

------

### 0. 前提の整理

原問題（各レベル $\ell$）：
$$
\begin{cases}
L_\ell u_\ell = f_\ell & \text{in }\Omega_\ell\\
\mathcal B u_\ell = g_\ell & \text{on }\partial\Omega_\ell
\end{cases}
$$
補正を
$$
u_\ell^{new} = u_\ell + e_\ell
$$
と置き、補正方程式

$$
L_\ell e_\ell = r_\ell,\qquad r_\ell=f_\ell-L_\ell u_\ell
$$

を解く。

------

### 1. 誤差の境界条件はどう決まるか（原理）

境界演算子 $\mathcal B$ が線形なら
$$
\mathcal B(u_\ell + e_\ell)=g_\ell
$$
より
$$
\boxed{
\mathcal B e_\ell = g_\ell - \mathcal B u_\ell
}
$$
つまり

- **Dirichlet**：境界値の誤差
- **Neumann**：法線微分の誤差

になります。

------

### 2. MG での基本方針（重要）

MG の補正段階では

> coarse grid は **「誤差だけを伝播させる」**

ので、**原則として誤差の境界条件は同次（ゼロ）**にします。

これは

- fine grid では BC はすでに満たされている（ghost 埋めで強制）
- coarse では **境界誤差を増やしたくない**

からです。

------

### 3. `Correction-Taylor` の境界条件（結論）

#### 3.1 Dirichlet 境界（最重要）

原問題：
$$
u=g \quad\text{on }\partial\Omega
$$
誤差：
$$
e = g - u \quad\text{on }\partial\Omega
$$
しかし **MG では coarse では**
$$
\boxed{
e=0 \quad\text{(homogeneous Dirichlet)}
}
$$
を使います。

#### 理由

- fine grid の (u) は毎回 Dirichlet BC を満たす
- 境界の誤差は **理論上ゼロ**
- coarse で非ゼロにすると補正が boundary-driven になり不安定

👉 **Correction-Taylor では**

```text
Dirichlet BC → e = 0
```

#### ghost 実装（例：x-min）

- 2次：
  $$
  e_{ghost} = - e_{adj}
  $$

- 4次：
  $$
  e_{ghost} = \text{(内点線形結合)} \quad\text{※ g 項なし}
  $$

（あなたが既に導出した **Dirichlet 4次 ghost 式から g を落とした形**）

------

#### 3.2 Neumann 境界

原問題：
$$
\frac{\partial u}{\partial n}=h
$$
誤差：
$$
\frac{\partial e}{\partial n}=h-\frac{\partial u}{\partial n}
$$
しかし MG 補正では
$$
\boxed{
\frac{\partial e}{\partial n}=0
}
$$
を使います（**同次 Neumann**）。

#### ghost（2次例）

x-min：
$$
e_{ghost}=e_{adj}
$$
（4次でも同様に **h 項を落とすだけ**）

------

#### 3.3 全 Neumann 問題の注意（超重要）

全境界が Neumann の場合：

- 原問題：解は定数分不定
- 誤差方程式：**nullspace を持つ**

👉 Correction-Taylor では必ず **ゲージ条件**を入れる：

例：
$$
\boxed{
\sum_{i,j,k} e_{i,j,k} = 0
}
$$
または
$$
e_{i_0,j_0,k_0}=0
$$
（各 Taylor step 後に平均値を引くのが簡単）

------

### 4. Correction-Taylor 用 境界条件まとめ（実装規約）

| 原問題の BC              | Correction-Taylor の BC      |
| ------------------------ | ---------------------------- |
| Dirichlet (u=g)          | **Dirichlet (e=0)**          |
| Neumann (\partial_n u=h) | **Neumann (\partial_n e=0)** |
| Mixed                    | 各面で上記を適用             |
| 全 Neumann               | 同次 Neumann + ゲージ固定    |

------

### 5. 擬似時間 Taylor 漸化式との整合

Correction-Taylor ODE：
$$
e_t = L e - r
$$

- 境界条件は **各 Taylor 係数 (e_m) に同じ同次 BC**
- RHS $r$ は境界では使わない（内点のみ）

よって、あなたの Taylor 漸化式
$$
(e)_{m+1}=\frac{1}{m+1}(L e_m - \delta_{m0}r)
$$
は **境界で完全に閉じる**。

------

### 6. なぜこれが「MG の本質」と一致するか

- MG の coarse correction は **誤差伝播**
- 誤差は境界でゼロ（Dirichlet）または流束ゼロ（Neumann）
- `Correction-Taylor` は **誤差 ODE をそのまま時間積分**

👉 境界条件を同次にすることで

> **Taylor 展開 = 誤差伝播演算子**

になります。

------

### 最終まとめ（短く）

- `Correction-Taylor` の BC は **必ず誤差の同次 BC**
- Dirichlet → $e=0$
- Neumann → $\partial_n e=0$
- 全 Neumann → ゲージ固定必須
- 実装上は **既存 ghost 式から g / h 項を消すだけ**





## Correction-Taylor を smoother として使う場合の安定条件

Correction-Taylor を **smoother（= 高周波誤差を減衰させる反復）**として使うときの安定条件は、結局
$$
e_t = L e - r \quad(\text{r は固定})
$$
の **誤差方程式**に帰着します。定常解 $e^*$ まわりで $\tilde e=e-e^*$ と置くと
$$
\tilde e_t = L\,\tilde e
$$
なので、**安定性は “(L)” のスペクトルと、あなたの Taylor 多項式（次数 $M$、$\Delta t$）が作る増幅因子**で決まります。RHS の (r) は安定性に影響しません。

------

### 1) 1モードに対する増幅因子（核心）

離散ラプラシアン $L$ の固有値を $\lambda \le 0$ とすると、モードは
$$
\tilde e_t = \lambda \tilde e
$$
で、Taylor $M$ 次の 1 step は（あなたの係数再帰の結果と同値）
$$
\tilde e^{new} = P_M(\lambda \Delta t)\,\tilde e
$$


$$
\boxed{
P_M(z)=\sum_{m=0}^{M}\frac{z^m}{m!}
}
\qquad (z=\lambda\Delta t\le 0)
$$


したがって安定条件は
$$
\boxed{
\max_{\lambda\in\sigma(L)}\;|P_M(\lambda\Delta t)|\le 1
}
$$
smootherとしてはさらに「高周波領域で十分小さい」こと（$|P_M|\ll 1$）が望ましいですが、最低限は上式です。

------

### 2) 実務で使う安定条件（上限 Δt の形）

欲しいのは「最も負の固有値 $\lambda_{\min}$」に対して安定な $\Delta t$。
つまり
$$
z_{\min}=\lambda_{\min}\Delta t
$$
が、$P_M$ の安定域（負の実軸上で $|P_M(z)|\le1$ が成り立つ区間）に入ること。

#### 2.1 $\lambda_{\min}$ の上界（よく使う評価）

- **2次精度 7点（dx,dy,dz）**なら
  $$
  |\lambda_{\min}|\;\le\; 4\left(\frac{1}{\Delta x^2}+\frac{1}{\Delta y^2}+\frac{1}{\Delta z^2}\right)
  $$
  （等方格子 $h$ なら $|\lambda_{\min}|\approx 12/h^2$ が目安）

- **4次（±2）**は係数が変わるので $|\lambda_{\min}|$ の定数が少し増えます（安全側に見るなら 2次より厳しく置く）。

#### 2.2 したがって Δt の安全上限は

$$
\boxed{
\Delta t \;\le\; \frac{\alpha_M}{|\lambda_{\min}|}
}
$$
ここで $\alpha_M$ は
$$
|P_M(-\alpha_M)|=1
$$
を満たす「負の実軸上の安定限界」です（次数 $M$ で決まる）。

------

### 3) 具体的な運用ルール（MGで使うならこれが強い）

Taylor truncated は “厳密な A-stable” ではないので、実務では「厳密な $\alpha_M$ を毎回求める」より **安全なスケーリング則**で十分強いです。

#### 推奨ルール（拡散CFL型）

$$
\boxed{
\Delta t_\ell = \theta\;\frac{h_\ell^2}{2d}
}
\quad(d=3,\;\theta\in(0,1])
$$

- 2次 7点・等方格子なら $2d=6$ を分母に置くのが古典的
- 実際は $L$ の定数係数を吸収して $\theta$を調整（例：$\theta=0.2\sim0.8$）

そしてレベル間はあなたの言う通り
$$
h_{\ell+1}=2h_\ell \Rightarrow \Delta t_{\ell+1}\approx 4\Delta t_\ell
$$
が自然に出ます。

#### smoother用途の追加条件（減衰を確保）

高周波の代表として $\lambda\approx \lambda_{\min}$ を狙い、
$$
|P_M(\lambda_{\min}\Delta t)| \le \rho_s
$$
（$\rho_s=0.1\sim0.3$ くらい）を目標に $\Delta t$ を選ぶと、スムーザーとして強いです。

------

### 4) 重要な注意点（MG＋Correction-Taylor特有）

#### (A) 安定性は “境界条件” でほぼ変わらないが、全Neumannはゲージ必須

固有値0（定数モード）があると
$$
P_M(0)=1
$$
で減衰しないので、smootherとしては「平均値固定」などの投影が必要です。

#### (B) coarse ほど Δt を大きくできるが、M を下げすぎると安定域が狭くなる

「粗格子は低次＋大Δt」は理にかないますが、
**M を下げると $\alpha_M$ が小さくなりがち**で、同じ 4倍則が破綻する場合があります。

実務では：

- $\Delta t_{\ell+1}=4\Delta t_\ell)$をまず置く
- 発散／振動したら **$\Delta t_\ell$ をクリップ**する
- もしくは coarse で **M を最低 2〜3 以上に固定**する

が安定です。

------

### 5) すぐ使える “設計レシピ”

（等方格子 $h_\ell$、2次7点を想定）

- fine で $\Delta t_1 = \theta \, h_1^2/6$（$\theta=0.3$ くらいから開始）

- coarse は $\Delta t_\ell = 4^{\ell-1}\Delta t_1$ だが
  $$
  \Delta t_\ell \leftarrow \min(\Delta t_\ell,\;\Delta t_{\max}(M_\ell))
  $$

- $M_\ell$ は例として

  - fine: 4〜6
  - mid: 3〜4
  - coarse: 2〜3（最低2）

------

### まとめ

- `Correction-Taylor` の安定性は **(\tilde e_t=L\tilde e)** のみで決まり、残差 (r) には依存しない

- 1モード増幅は
  $$
  P_M(\lambda\Delta t)=\sum_{m=0}^M \frac{(\lambda\Delta t)^m}{m!}
  $$

- 安定には
  $$
  |P_M(\lambda_{\min}\Delta t)|\le1
  $$
  を満たすよう (\Delta t) を選ぶ

- MG では基本的に $\Delta t \propto h^2$（よって 2倍粗で 4倍 dt）が自然だが、**低次数化しすぎるとクリップが必要**