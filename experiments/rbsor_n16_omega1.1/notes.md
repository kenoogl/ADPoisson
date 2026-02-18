# Experiment Notes: <exp>

---

## 1. AI Summary (Facts)
・対象実験: sor_n16_omega1.1（solver: sor, n=16, omega=1.1, max_steps=20000, epsilon=1.0e-8）。
・run_summary.jsonでは、反復回数は333回、runtimeは0.00980997085571289秒。
・convergedはtrue、最終residual_l2は9.759974e-9。
・誤差指標はerror_l2=0.0017228561237746108、error_max=0.007576462867579181。
・履歴ファイルはhistory_sor_nx16_ny16_nz16_steps333.txt。
・異常終了や未収束は記録されていない。

## 2. AI Analysis (Evaluation)

### Alignment with Intent
・成功基準のうち、所定許容誤差までの収束達成と短い実行時間の確認には整合している。
・ただし、Intentにある「ベースライン比較」や「Taylor法の有効性評価」は単独実験だけでは判定できない。

### Hypotheses
・同一格子（n=16）では、omega=1.1はSORの安定収束域にあり、反復333回で収束する設定である可能性が高い。
・収束性は確保できているため、性能改善余地はomega最適化または別solver比較に残っている可能性がある。

### Comparison
・本結果単体では、SSOR/CG/Taylor系列との優劣は未確定。
・比較位置づけを明確化するには、同一条件（n, Fo, epsilon）でのrun_summary横断比較が必要。

### Minimal Next Experiment
・最小追加実験として、同条件で`omega=1.3`のSORを1本実行し、iterationsとruntime_secを`sor_n16_omega1.1`と比較する。

---

## 3. Human Thoughts (Decision)
(手書き)
- 採用 / 保留 / 却下
- なぜそう判断したか
- 次の実験名
- 気になる点
- 投資するかどうか

---