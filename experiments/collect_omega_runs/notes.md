# Experiment Notes: collect_omega_runs

---

## 1. AI Summary (Facts)
### Run Overview
- 対象は `results/collect_omega_runs/omega_runs_summary.md` の集計結果（SOR/SSOR の omega 探索一覧）。
- 格子サイズは `nx=16,32,64,128`、omega は `1.0` から `1.8`。
- SOR は全格子で広い範囲が収束し、`nx=16` では `omega=1.8` が発散（`res_l2=inf`）。
- SSOR は `omega=1.4` 以上で多くのケースが発散（全格子で `res_l2=inf` が多数）。
- 集計表は `results/collect_omega_runs/omega_runs_summary.md` に出力済み。

### Observed Behavior
- SOR は omega 増加に伴い反復回数が減少する傾向が確認される。
- SSOR は低omega側（概ね `1.0-1.3`）で収束し、高omega側で発散が顕著。
- `nx=128` の SOR は `omega=1.0` で 20000 step 到達、`omega=1.1` 以上で収束。

##### ※ ここでは解釈しない。集計表（omega_runs_summary.md）の内容を整理するのみ。


## 2. AI Analysis (Evaluation)
### Alignment with Intent
- Baseline（SOR/SSOR）比較データの整理という目的に整合している。
- Intent の比較基盤として、収束可否・反復回数・runtime の横比較が可能な形になっている。

### Hypotheses 仮説候補
- SOR は格子ごとに最適omega帯があり、`n=16` では `1.6-1.7` 付近が有望。
- SSOR は安定領域が狭く、`omega>=1.4` で不安定化する可能性が高い。
- 格子が大きくなるほど、SOR の最適omegaは高め側にシフトする可能性がある。

### Counterarguments 他実験との比較
- SOR は SSOR より発散領域が狭く、探索レンジ全体で実用的な候補点が多い。
- SSOR は低omegaで収束するが、高omega側で一貫して発散しており探索効率が悪い。
- ただし比較は同一停止条件での結果に依存するため、別停止条件では順位が変わる余地がある。

### 最小反復回数での評価
- `nx=16`: SOR 最小 `91`（`omega=1.6`） vs SSOR 最小 `102`（`omega=1.3`）。
- `nx=32`: SOR 最小 `273`（`omega=1.7`） vs SSOR 最小 `384`（`omega=1.3`）。
- `nx=64`: SOR 最小 `684`（`omega=1.8`） vs SSOR 最小 `1432`（`omega=1.3`）。
- `nx=128`: SOR 最小 `2645`（`omega=1.8`） vs SSOR 最小 `5302`（`omega=1.3`）。
- 最小反復回数のみで比較すると、全格子で SOR が優位。

### Minimal Next Experiment 次実験提案
- `nx=128` の SOR を `omega=1.65` で1本追加し、`omega=1.6` と `1.7` の間で反復回数最小点を確認する。

---

## 3. Human Thoughts (Decision)
- SORはn=16でomega=1.6、n=32でomega=1.7、n=64,128ではomega=1.8が反復回数が最小
- SSORではn=16,32,64,128の全ケースでomega=1.3が反復回数最小、omega=1.4以上で発散している
- SORとSSORでomega最適値の傾向が異なる。
- SORでは格子数により最適値は変化する。
- SSORではomega=1.3が最適であるが、少し大きくなると発散する
- SSORの反復回数がSORよりも多い。これはSSORはRBSORの実装で、ポイントSORベースではないため、スイープ順が適切でない可能性がある。

---
