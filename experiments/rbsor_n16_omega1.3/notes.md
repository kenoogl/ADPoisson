# Experiment Notes: <exp>

---

## 1. AI Summary (Facts)
・対象実験: `sor_n16_omega1.3`（solver=`sor`, `n=16`, `omega=1.3`, `epsilon=1.0e-8`）  
・問題設定: 3次元Poisson方程式、Dirichlet境界条件、格子 `16×16×16`  
・実行結果: `iterations=219`、`converged=true`、`residual_l2=9.723001e-9`  
・誤差指標: `error_l2=0.0017228547248401929`、`error_max=0.007576461485101493`  
・runtime: `0.006287097930908203 sec`  
・異常有無: `run_summary.json`上で異常を示すフラグや失敗記録はなく、収束完了

---

## 2. AI Analysis (Evaluation)

### Alignment with Intent
・本実験はベースライン（SOR）の単独計測として、今後のTaylor系手法比較に必要な基準データという位置づけでIntentと整合する  
・成功基準（反復数・時間・残差履歴比較）のうち、本実験単体では「比較評価」自体は未完了で、基準点の取得段階と評価できる  
・収束性の観点では、所定許容誤差まで到達しており、安定性確認の最小条件は満たしている

### Hypotheses
・`n=16`では`omega=1.3`のSORは少ない反復回数で安定収束する可能性が高い  
・同一格子で`omega`を調整すると、219反復よりさらに短縮できる最適域が存在する可能性がある  
・Taylor系手法が有効であるためには、少なくとも本条件（219反復・約0.0063秒）を上回る収束効率を示す必要がある

### Comparison
・現時点では他実験（SSOR, CG, Taylor系, multigrid）数値が未提示のため、定量比較は未実施  
・位置づけとしては「Poisson 3D, n=16におけるSOR基準ケース」であり、以後の比較軸（反復回数、runtime、誤差）の参照点になる

### Minimal Next Experiment
・`sor_n16_omega1.5` を1本だけ追加実行し、`omega=1.3`との差（iterations/runtime/residual到達性）を比較する

---

## 3. Human Thoughts (Decision)
(手書き)
- 採用 / 保留 / 却下
- なぜそう判断したか
- 次の実験名
- 気になる点
- 投資するかどうか

---