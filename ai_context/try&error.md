### 実験に応じたconfig.yamlの作り方

~~~
./bin/new_config base

実行すると、次ができる
experiments/base/config.yaml
~~~

~~~
experiment:
  name: base
  description: user's comment

output:
  results_dir: results/base

execution:
  command: >
    julia --project=. src/main.jl
    --output-dir results/base
~~~

これに具体的な内容を書き足す。以下は、Taylor uniformの例

~~~
experiment:
  name: taylor_64_uniform
  description: 以前の結果と修正コードの結果が一致するかを確認する

output:
  results_dir: results/taylor_64_uniform
  save_residual_history: true

execution:
  command: >
    julia --project scripts/main.jl 
    --solver mg-uniform-taylor 
    --n 64 --Fo 0.5 --M 4 --max-steps 20000 
    --epsilon 1e-8 --alpha 1.0 
    --bc-order spec 
    --mg-interval 5 
    --mg-dt-scale 2.0 
    --mg-M 2 
    --mg-nu1 1 --mg-nu2 1 
    --lap-order second 
    --output-dir results/taylor_64_uniform

problem:
  equation: poisson
  dimension: 3
  domain:
    type: cartesian
    size: [1.0, 1.0, 1.0]
  boundary_condition:
    type: dirichlet
    solution: analytic

discretization:
  grid:
    nx: 64
    ny: 64
    nz: 64
  layout: cell-centered
  ghost_cells: 2

solver:
  method: taylor
  type: uniform

convergence:
  residual_norm: L2
  tolerance: 1.0e-8
  max_iterations: 5000
~~~

`bin/run_exp taylor_64_uniform`で実行し、`logs/taylor_64_uniform.json`と`rsults/taylor_64_uniform`ディレクトリ下に結果が生成されることを確認。



- yamlで記述されたパラメータとコードのコマンドライン引数で二重にパラメータが記述され、両者が一致しない場合の不具合が懸念される。
- シングルソースにすべき。
- yamlパラメータを正とするようにJuliaコードを修正







格子数違いのconfig.yamlの作成。次のひな形の`__NN__`, `__EXP_NAME__`をパターンで置換し、ディレクトリにconfig.yamlを配置させる

~~~
experiment:
  name: __EXP_NAME__
  description: cg without preconditioner

output:
  results_dir: results/__EXP_NAME__
  save_residual_history: true

run:
  solver: cg
  cg_precond: none
  n: __NN__
  Fo: 0.5
  M: 4
  max_steps: 20000
  epsilon: 1.0e-8
  alpha: 1.0
  bc_order: spec
  output_dir: results/__EXP_NAME__

problem:
  equation: poisson
  dimension: 3
  domain:
    type: cartesian
    size: [1.0, 1.0, 1.0]
  boundary_condition:
    type: dirichlet
    solution: analytic

discretization:
  grid:
    nx: __NN__
    ny: __NN__
    nz: __NN__
  layout: cell-centered
  ghost_cells: 2

convergence:
  residual_norm: L2
~~~



##### CG

プロンプトを次のように与え、16,32,64,128,256の各ケースを作成。

~~~
以下のプレビューだけ出してください。

目的:
- CG法（前処理なし）の実験設定ファイルを、NN=16,32,64,128,256 の5ケース分作成する。

入力テンプレート:
- experiments/base/config.yaml

置換ルール:
- __EXP_NAME__ -> cg_none_<NN>
- __NN__ -> <NN>
- <NN> は 16, 32, 64, 128, 256

出力先（作成時）:
- experiments/cg_none_<NN>/config.yaml

重要:
- まだファイルは作成しないこと（dry-run）。
- まず以下のみ提示すること:
  1) 生成予定ファイル一覧
  2) 書き込み先ディレクトリ
  3) 各 config.yaml の最終内容（5件すべて）

上記を確認して私が「実行して」と返答したら、実際に作成してください。
~~~

##### SOR/SSOR

sor, ssorの加速係数omegaの最適値を探すため、config.yamlを作成

~~~
以下を実施してください。

目的:
- `experiments/base/config_sor.yaml` をテンプレートとして、
  緩和係数 `omega` の最適値を探索する実験設定を作成する。
- 探索条件:
  - 格子数 `NN = 16, 32, 64, 128`
  - 緩和係数 `OMEGA = 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6`
  - ソルバー `SOLVER = sor`

置換ルール:
- `__EXP_NAME__` -> `<SOLVER>_n<NN>_omega<OMEGA>`
- `__NN__` -> `<NN>`
- `__OMEGA__` -> `<OMEGA>`
- `__SOLVER__` -> `<SOLVER>`

作成ルール:
- 各組み合わせ `(NN, OMEGA)` ごとに 1 つの実験ディレクトリを作成し、
  `config.yaml` を生成する。
- 生成先:
  - `experiments/<SOLVER>_n<NN>_omega<OMEGA>/config.yaml`
- `OMEGA` は小数1桁表記を維持する（例: `1.0`, `1.1`, ..., `1.6`）。

必須反映項目（テンプレート内）:
- `experiment.name = <SOLVER>_n<NN>_omega<OMEGA>`
- `output.results_dir = results/<SOLVER>_n<NN>_omega<OMEGA>`
- `execution.command` は config参照型で、
  `experiments/<SOLVER>_n<NN>_omega<OMEGA>/config.yaml` を指すこと
- `run.n = <NN>`
- `run.omega = <OMEGA>`
- `discretization.grid.nx = <NN>`
- `discretization.grid.ny = <NN>`
- `discretization.grid.nz = <NN>`

それ以外の項目は `experiments/base/config_sor.yaml` をそのまま保持する。

出力:
1) 作成した実験名一覧
2) 作成した `config.yaml` のパス一覧
3) 代表として `NN=16` の全 `OMEGA` ケースの `config.yaml` 内容

~~~



##### 実験スクリプト

文字列に合致するディレクトリ名を探し、run_expで実行

~~~
bin/run_exp_patterns

#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./bin/run_exp_patterns '<pattern>'

Examples:
  ./bin/run_exp_patterns 'sor_n*_omega*'
EOF
}

if [ $# -ne 1 ]; then
  usage
  exit 1
fi

PATTERN="$1"

EXPS_RAW="$(find experiments -mindepth 1 -maxdepth 1 -type d -name "$PATTERN" -exec basename {} \; | sort)"
if [ -z "$EXPS_RAW" ]; then
  echo "No experiments matched: experiments/$PATTERN"
  exit 0
fi

echo "Pattern: $PATTERN"
echo "Matched experiments:"
echo "$EXPS_RAW" | sed 's/^/  - /'

FAILED=""
SKIPPED=""
TOTAL=0
FAILED_N=0
SKIPPED_N=0

while IFS= read -r exp; do
  [ -z "$exp" ] && continue
  TOTAL=$((TOTAL + 1))
  cfg="experiments/${exp}/config.yaml"
  if [ ! -f "$cfg" ]; then
    echo "[skip] missing config: $cfg"
    SKIPPED="${SKIPPED}${exp}"$'\n'
    SKIPPED_N=$((SKIPPED_N + 1))
    continue
  fi

  echo "[run] $exp"
  if ./bin/run_exp "$exp"; then
    echo "[ok]  $exp"
  else
    echo "[ng]  $exp"
    FAILED="${FAILED}${exp}"$'\n'
    FAILED_N=$((FAILED_N + 1))
  fi
done <<EOF
$EXPS_RAW
EOF

echo
echo "===== Summary ====="
echo "pattern : $PATTERN"
echo "total   : $TOTAL"
echo "failed  : $FAILED_N"
echo "skipped : $SKIPPED_N"

if [ "$FAILED_N" -gt 0 ]; then
  echo
  echo "Failed experiments:"
  printf "%s" "$FAILED" | sed '/^$/d; s/^/  - /'
fi

if [ "$SKIPPED_N" -gt 0 ]; then
  echo
  echo "Skipped experiments (missing config):"
  printf "%s" "$SKIPPED" | sed '/^$/d; s/^/  - /'
fi

~~~

この後実験結果をみて、omega=1.6では反復回数が下がりきっていないので1.7, 1.8を追加するため、config.yamlの作成を依頼

~~~
sor_n*_omega1.7のconfig.yamlを作成してください
sor_n*_omega1.8のconfig.yamlを作成してください
~~~



##### sor/ssorで加速係数の最適化

sorとssorのomegaパラメータを変えた結果,

~~~
`sor_n*_omega*`, `ssor_n*_omega*`にあるrun_summary.tomlから格子数（nxを代表値として扱う）、反復回数steps、残差res_l2、実行時間runtime_sの一覧mdファイルを作成するスクリプトをcollect_runs.jlを参考に作成。res_l2がinf or nanの場合には発散しているとして反復回数は∞とする。カテゴリ毎にまとめて表示。README.mdに追記。
~~~

`collect_omega_runs.jl`が作成されるが、出力ディレクトリの指定ができないので、下記プロンプトを投げる

~~~
collect_omega_runs.jlもrun_expで実行する。実行のconfig.yamlを作成し、結果をフレームワークのルールに従って出力する
~~~

`experiments/collect_omega_runs/config.yaml`が作成され、一覧の結果がでた。







------

了解。あなたの期待は「AIを“自動実行者”にする」ではなく、

> **実験という事実（data）から、思考（analysis → hypothesis → debate → next exp）を加速する**

に最大化したい、ということですよね。
そのために一番効くやり方は、AIが“考える材料”を毎回同じ形で受け取れるようにして、**思考のループを儀式化**することです。

以下、「最小の複雑さ」で最大の思考支援が出る設計を提案します。

------

## 1) AI思考支援を最大化する鍵は「入力の固定化」

AIが強いのは、**比較・要約・論点抽出・反証・代替案生成**です。
弱いのは、入力が毎回バラバラで「何を見ればいいか」から推測させられる状況。

なので、各実験についてAIに渡す“パケット”を固定します。

### ✅ 実験パケット（AIが必ず読む順）

1. `ai_context/intent.md`（意図・成功条件）
2. `experiments/<exp>/config.yaml`（条件）
3. `logs/<exp>.json`（実行の事実：commit/dirty/command）
4. `results/<exp>/`（図表・csv・数値）
5. `ai_context/codex_results.md`（事実の文章化：AIが書記として整える）

この5点セットが揃うと、AIは“思考”に全振りできます。

------

## 2) 5〜7を「AIの価値が出る順」に再定義する

あなたが言う価値（図表化→分析→仮説→ディベート→新アイデア）は、5〜7をこう分けると最大化します。

### 5️⃣ **Facts（事実の固定）** = AIの書記能力

- **トリガー**：`run_exp` が終わって `logs/<exp>.json` と `results/<exp>/` が揃った
- **AIが見る**：logs + results + config
- **AIが出す**：`ai_context/codex_results.md`
- **ルール**：ここでは“評価”しない（事実だけ）

> ここでAIの強み：自動要約、表/グラフ案、異常値検知、比較軸抽出

------

### 6️⃣ **Judge（批判的評価・反証・比較）** = AIの査読能力

- **トリガー**：`codex_results.md` ができた
- **AIが見る**：plan + config + logs + results + codex_results
- **AIが出す**：`ai_context/judge_reviews.md`
- **出力の必須**：
  - PASS / REVISE / STOP（いずれか）
  - Reviewer #2 が突く3点
  - 次実験の「最小一手」（1〜3個）

> ここでAIの強み：反証、限界、一般化、比較不足の指摘、追加実験の提案

------

### 7️⃣ **PI Decision（意思決定）** = 人間の責任 + AIのディベート補助

- **トリガー**：`judge_reviews.md` が出た
- **人が決める**：採用/却下/保留、次の実験名、投資する時間
- **書く**：`ai_context/atlas_notes.md`（短くてOK）

> ここでAIの強み：ディベート相手、意思決定の理由の論理化、代替戦略の提示

------



------

## 6) 最小の運用ルール（これだけ守れば回る）

1. 実験前：`codex_plan.md` を更新（雑でOK）
2. 実験作成：`new_config exp_name` → config追記
3. 実行：`run_exp exp_name`（logs/results生成）
4. **AIに Facts を作らせる** → `codex_results.md`
5. **AIに Judge を作らせる** → `judge_reviews.md`
6. 人が決める → `atlas_notes.md`
7. commit（experiments/results/logs/ai_context をまとめて）

これで、あなたが期待する「思考の加速」が 5〜7 に集約されます。

------









＋＋＋＋＋＋＋＋＋＋＋＋＋＋＋

AIによる思考支援を最大化する

AIが強いのは、**比較・要約・論点抽出・反証・代替案生成**です。
弱いのは、入力が毎回バラバラで「何を見ればいいか」から推測させられる状況。

そのため、実験を記録し、データをAIが活用しやすい形にする。これは同時に実験の再現性やprovenanceのエビデンスとなる。

欲しい機能は：

- AIが思考を拡張してくれること
- 実験から仮説が自然に出てくること
- 論文化に近づくこと



##### 実験後にAIに渡すもの

- config.yaml
- logs/<exp>.json
- results/<exp>/metrics.json（あれば）

これだけ。

AIは：

- 要約
- 比較
- 仮説生成
- 追加実験提案

全部できる。

＋＋＋＋＋＋＋＋＋＋＋＋＋＋＋＋



実験結果としてconfig.yaml、logs/<exp>.json、results/<exp>/metrics.jsonがあれば、効果的な支援ができる。これに加えて必要な最低限のドキュメントは何？また、これらをベースに具体的な支援シナリオを示して

------



------

# 🧠 なぜこれで十分か？

AIは：

- metrics → 数値
- config → 条件
- logs → 実行状態
- intent → 意図

この4点があれば、

- 分析
- 比較
- 仮説生成
- 追加実験提案
- 論文化視点
- Reviewer視点

全部できます。



------

# 具体的支援シナリオ

## 🎬 シナリオ1：最適値探索

#### ① 自動比較表生成

#### ② 最適点候補提示

#### ③ 仮説生成

#### ④ 次実験提案

------

## 🎬 シナリオ2：論文化可能性判定

AIに聞く：

> この実験は論文化の核になるか？

AIが評価する観点：

- 既存理論との整合性
- 一般化可能性
- 比較不足点
- 新規性の強さ

------

## 🎬 シナリオ3：ディベートモード

あなた：

> 推進派として議論して

AI：

- MG + Taylor は新しい
- 数値的に明確な改善あり

あなた：

> Reviewerとして反論して

AI：

- 既存文献の拡張では？
- grid依存性未検証

これが思考拡張。

------

## 🎬 シナリオ4：異常検出

metrics.json が増えたとき：

AIは：

- 異常値
- 非単調挙動
- 収束不安定
- runtime不自然

を検出できます。

------

## シナリオ５：AIに比較させる

複数実験の logs を渡し：

> どれが最も論文価値が高い？

と評価させる。

------

## シナリオ６： 「新アイデア」を生むための最強パターン：二役ディベート

Judge を1体にせず、**同じ材料で2つの人格に議論させる**のが最も発想が出ます（追加の複雑さはほぼゼロ）。

- **Pro（推進派）**：「この結果は価値がある。次に何をすべきか」
- **Con（懐疑派/Reviewer #2）**：「それは本当に新規性？再現性？比較は？」

出力は `judge_reviews.md` に統合してOK。
この“対立構造”が、思考を強制的に立体化します。

------

## 3) 図表化と数値分析は「AIにやらせる」ではなく「AIが迷わない形で置く」

AIが図表や分析を“思考の燃料”にするために、results に最低限これだけ出すのが効きます（Poisson/MGやSOR/SSORにも共通）。

```
results/<exp>/
```

- `metrics.json`（最重要）
  - iterations, final_residual, runtime_sec, converged(true/false), key_params(ωなど)
- `residual_history.csv`（あれば）
- `figures/`（自動生成できるなら）

なぜ `metrics.json` が効くか：
AIが比較・要約・表化するとき、**CSVを逐一読ませなくても済む**からです。思考が速くなります。



AIはファイルが多い、意図が曖昧、情報が散乱、だと弱い。

AIが考えるための「最小の構造化」を用意するそれ以上は増やさない



- `ai_contentext/atlas_notes.md > project_notes.md`　プロジェクトに関する様々なこと
- `ai_contentext/codex_results.md > intent.md` 研究目的、ゴール、成功基準、比較

- `experiments/<exp>/notes.md`　各実験のメモ

- `config.yaml`　実験条件
- `logs/<exp>.json`　実験環境
- `results/<exp>/`　実験結果、データ、図表
- `results/<exp>/run_summary.json`

`ai_context`ディレクトリには、人が読むだけのメモなども可（この文書）



#### run_summary.json の設計

AIが「比較・仮説生成・異常検出・論文化評価」をしやすい最小構造を作成する

##### 🎯 設計原則

1. 人間が読める
2. AIが比較しやすい（フラット）
3. 実験間で必ず共通キーがある
4. solverに依存しすぎない
5. 将来拡張可能

- AI思考支援を最大化するための条件として、議論に必要なすべての実験でキーが揃っていること

- `results/<exp>/run_summary.json` を solver が出す。最初から出力する項目が決まればよいが、そうでない場合もある。あとで追加？



##### `notes.md` AIに結果を記述させる。セクションを分けて人の考えも書く

~~~
# Experiment Notes: sor_omega_1p25

---

## 1. AI Summary (Facts)
(自動生成)
metrics要約
収束状況
runtime
異常有無
解釈しない。
---

## 2. AI Analysis (Evaluation)
(自動生成)
intentとの整合性
次を必ず書かせる
### Hypotheses 仮説候補
- ...
### Counterarguments 他実験との比較
- ...
### Minimal Next Experiment 次実験提案
- exp_name: ...
  change: ...
---

## 3. Human Thoughts (Decision)
(手書き)
採用／却下
次exp名
投資判断
疑問
---
~~~



#### `history_stats.json`

ヒストリデータは生データの統計処理情報（特徴量）をjson化

~~~
{
  "history_stats": {
    "monotonic": false,
    "oscillation_detected": true,
    "diverged": true,
    "initial_residual": 1.0,
    "min_residual": 0.9640102,
    "final_residual": 3.605399,
    "convergence_rate_estimate": -0.12
  }
}

~~~



#### `run_summary.json`

~~~
{
  "timestamp": "2026-02-14T09:46:01",
  "config_path": "experiments/ssor_n128_omega1.8/config.yaml",
  "script": "scripts/run_solver.jl",

  "iterations": 162,
  "runtime_sec": 2.5548958778381348,

  "converged": false,
  "residual_l2": "inf",
  "error_l2": 4.391278058068776e146,
  "error_max": 2.0418659267440394e149,
  
  "artifacts": {
    "history": "history_ssor_nx128_ny128_nz128_steps162.txt"
  }
}
~~~

