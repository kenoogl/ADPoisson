# Experiment Notes: sor_n16_omega1.6

---

## 1. AI Summary (Facts)
### Run Overview
- 実験名: `sor_n16_omega1.6`。
- ソルバー設定: `solver=sor`, `n=16`, `omega=1.6`, `Fo=0.5`, `epsilon=1.0e-8`, `max_steps=20000`。
- 実行結果: `iterations=91`, `runtime_sec=0.002624988555908203`, `converged=true`。
- 残差・誤差: `residual_l2=8.900214e-9`, `error_l2=0.0017228523926655897`, `error_max=0.007576459174722272`。
- 出力履歴: `history_sor_nx16_ny16_nz16_steps91.txt`。

### Observed Behavior
- 最大反復到達前に収束して終了。
- `run_summary.json` 上で NaN/Inf などの異常値は確認されない。

##### ※ ここでは解釈しない。run_summary.json の内容を整理するのみ。


## 2. AI Analysis (Evaluation)
### Alignment with Intent
- Baseline（SOR）として、比較用の基準データを提供しており Intent と整合。
- 収束判定は `converged=true` で、設定閾値 `epsilon=1.0e-8` を満たしている。
- `n=16` 条件で `omega` 探索の進展を示すデータとして有効。

### Hypotheses 仮説候補
- `n=16` では `omega=1.6` が `omega=1.5` よりさらに反復回数を削減できる可能性がある。
- これ以上 `omega` を増やすと、収束は維持されても誤差や安定性に悪影響が出る境界がある可能性がある。
- 低格子では runtime が非常に短いため、差分評価は iterations と収束余裕度で見るのが妥当。

### Counterarguments 他実験との比較
- 既知の `sor_n16_omega1.5`（131反復）と比較すると、本ケース（91反復）はさらに少ない。
- `sor_n16_omega1.0`（406反復）から見ると大幅改善傾向が継続している。
- 一方で `n=16` は小規模格子なので、同傾向が `n=32/64/128` でも成立するかは未確認。

### Minimal Next Experiment 次実験提案
- `n=16` 固定で `omega=1.7` を1本追加し、収束性（converged）と反復回数の改善継続可否を確認する。

---

## 3. Human Thoughts (Decision)
(手書き)
- 採用 / 保留 / 却下
- なぜそう判断したか
- 次の実験名
- 気になる点
- 投資するかどうか

---
