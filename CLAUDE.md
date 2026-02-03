# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## プロジェクト概要

ADPoissonは、三次元ポアソン方程式をTaylor級数時間積分で解くJuliaソルバー。姉妹プロジェクト `ADburgers`（`/Users/Daily/Development/ADTM/ADburgers`）の設計を踏襲する。

## 開発規約

- 常に日本語で会話する
- TDD（テスト駆動開発）で進める：テスト作成 → 失敗確認 → コミット → 実装
- 2スペースインデント
- Kiro Spec-Driven Development（`.kiro/`）のワークフローに従う

## コマンド

```bash
# 依存関係インストール
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# テスト実行
julia --project=. -e 'using Pkg; Pkg.test()'

# 単一テストファイル実行
julia --project=. test/runtests_core.jl

# REPLでの開発
julia --project=.
```

## アーキテクチャ（ADburgers準拠の想定構成）

```
src/
  ADPoisson.jl    # モジュールエントリ（exports, includes）
  types.jl        # ProblemSpec, SolverSpec, Solution, TaylorArrays
  core.jl         # solve(), taylor_coeff!(), horner_update!()
  problems.jl     # テスト問題定義（解析解付き）
  factory.jl      # get_problem(id) ファクトリ
  analysis.jl     # 誤差計算（L2/L∞ノルム）、受容基準チェック
  visualization.jl # RecipesBaseプロットレシピ
test/
  runtests.jl          # テストハーネス
  runtests_core.jl     # Taylorエンジン単体テスト
  runtests_problems.jl # 問題検証テスト
scripts/               # 検証・可視化スクリプト
```

### 核心アルゴリズム

Taylor係数の再帰計算で時間発展を高次精度で実行。空間は中心差分で離散化し、`TaylorArrays`に事前確保したメモリ上で計算する。Horner法でTaylor級数を評価して時間ステップを進行。

### 主要な型

- **ProblemSpec**: 物理問題定義（領域、係数、初期条件、境界条件、厳密解）
- **SolverSpec**: ソルバー設定（N, Δt, K, 保存時刻）
- **Solution**: 計算結果（格子、時刻、解行列、メタデータ）
- **TaylorArrays**: Taylor係数の事前確保配列

## 参照実装

ADburgers（`/Users/Daily/Development/ADTM/ADburgers`）が完成済みリファレンス。設計パターン、型構造、テスト戦略はこれに従う。
