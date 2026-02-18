# Experiment Notes: sor_n16_omega1.0

---

## 1. AI Summary (Facts)
・対象実験は `sor_n16_omega1.0`、ソルバは SOR（`omega=1.0`）、格子は `16×16×16`、収束判定は残差 L2 ノルム `epsilon=1.0e-8`。  
・実行結果は `iterations=406`、`runtime_sec=0.012026071548461914`。  
・`converged=true`、最終 `residual_l2=9.80302e-9`、`error_l2=0.0017228567799728253`、`error_max=0.007576463515522103`。  
・履歴ファイルは `history_sor_nx16_ny16_nz16_steps406.txt`。  
・異常有無: 収束フラグは true で、max_steps（20000）到達前に終了。実行失敗を示す事実は記録されていない。  

---

## 2. AI Analysis (Evaluation)

### Alignment with Intent
・ベースライン（SOR）として、所定許容誤差まで収束しており、Intent の比較対象データとして有効。  
・ただし本実験単体では「Taylor 系手法が SOR/SSOR より高速」という成功基準の達成可否は判断できない（比較実験が未提示）。  

### Hypotheses
・`n=16` では SOR（`omega=1.0`）が安定に収束するため、同一条件下で Taylor 系を適用すると、反復回数または時間で優位性の有無が明確に現れる可能性がある。  
・本ケースの反復数（406）は、平滑化特性比較（高周波誤差減衰）に使う基準点として機能する。  

### Comparison
・現時点では他実験値（SSOR, CG, Taylor, multigrid）が与えられておらず、相対順位は未確定。  
・位置づけとしては「Poisson 3D・`16^3`・SOR(ω=1.0) の基準ラン」。今後の比較軸（反復回数、runtime、残差履歴、誤差分布）の参照点になる。  

### Minimal Next Experiment
・**同一設定でソルバのみ SSOR に変更した 1 実験**（`n=16`, `epsilon=1e-8` を固定）を実施し、SOR 基準との反復回数・runtime・残差履歴を直接比較する。

---

## 3. Human Thoughts (Decision)

(ここは人間が書く)

---