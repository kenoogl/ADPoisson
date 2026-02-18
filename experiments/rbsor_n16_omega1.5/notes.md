# Experiment Notes: sor_n16_omega1.5

---

## 1. AI Summary (Facts)
### Run Overview
- 実験名: `sor_n16_omega1.5`。
- ソルバー設定: `solver=sor`, `n=16`, `omega=1.5`, `Fo=0.5`, `epsilon=1.0e-8`, `max_steps=20000`。
- 実行結果: `iterations=131`, `runtime_sec=0.0038568973541259766`, `converged=true`。
- 残差・誤差: `residual_l2=9.638116e-9`, `error_l2=0.0017228532668104649`, `error_max=0.007576460042261868`。
- 出力履歴: `history_sor_nx16_ny16_nz16_steps131.txt`。

### Observed Behavior
- 最大反復到達前に収束して終了。
- `run_summary.json` 上で NaN/Inf などの異常値は確認されない。

##### ※ ここでは解釈しない。run_summary.json の内容を整理するのみ。


## 2. AI Analysis (Evaluation)
### Alignment with Intent
- Baseline（SOR）として、比較用の基準データを提供しており Intent と整合。
- 収束判定は `converged=true` で、設定閾値 `epsilon=1.0e-8` を満たしている。
- ただし Intent の「他手法比で高速化」の判定には、SSOR/CG/Taylor 系との横比較が必要。

### Hypotheses 仮説候補
- `n=16` では `omega=1.5` が `omega=1.0` より反復回数を減らす可能性が高い。
- 同じ格子で `omega` をさらに上げると、反復回数がさらに減るか、逆に不安定化する境界が現れる可能性がある。
- 低格子では runtime 差が小さいため、最適 `omega` の評価は iterations と収束安定性の両方で見る必要がある。

### Counterarguments 他実験との比較
- 既知の `sor_n16_omega1.0`（406反復）と比較すると、本ケース（131反復）は大幅に少ない。
- ただし error 指標の差は小さく、収束後精度は同程度とみなせる。
- `n=32/64/128` でも同傾向かは未確認で、格子依存性は追加検証が必要。

### Minimal Next Experiment 次実験提案
- `n=16` 固定で `omega=1.6` を1本追加し、`iterations` と `converged` の変化を確認する。

---

## 3. Human Thoughts (Decision)
(手書き)
- 採用 / 保留 / 却下
- なぜそう判断したか
- 次の実験名
- 気になる点
- 投資するかどうか

---
