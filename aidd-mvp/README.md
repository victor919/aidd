# aidd-mvp · EGFR 端到端 AI 制药流水线（教学 MVP）

一条能真的跑通的 AIDD Pipeline：疾病 / 靶点 → 蛋白结构与口袋 → 公开活性数据 →
清洗 + 切分 → Baseline → Advanced Model → 对接 → ADMET → 干湿闭环。

> 目标不是发现新药，也不是刷 SOTA，而是**跑通一次**。
> 打通之后，以后每学一个模型，都知道它属于这条流水线上的哪个工位。

---

## 快速开始

```bash
# 环境：Python 3.11 + RDKit + scikit-learn + PyTorch + chemprop + AutoDock Vina
python scripts/01_target_rank.py
python scripts/02_download_structure.py
python scripts/03_extract_pocket.py
python scripts/04_fetch_activity.py        # 联网拉 ChEMBL；离线加 AIDD_OFFLINE=1
python scripts/05_standardize.py
python scripts/06_scaffold_split.py
python scripts/07_rf_baseline.py
python scripts/08_chemprop.py
python scripts/09_docking.py
python scripts/10_admet.py
python scripts/11_closed_loop.py
```

或一次跑完：

```bash
bash run_all.sh
```

所有参数集中在 `config.yaml`，脚本内**不允许**硬编码 `"P00533"` / `2048` / `42`。

---

## 可运行 / 可测试（三层保障，2026-10-17 实跑确认）

| 层 | 命令 | 检查什么 | 实跑结果 |
|---|---|---|---|
| 环境预检（左端守门） | `python scripts/00_preflight.py` | 必需依赖 / config 必填键 / 目录可写；缺必需 → 退出码 1，缺可选（vina、torch）→ WARN 并自动跳过对应步骤 | PASS |
| 产物自洽 | `python scripts/check_artifacts.py` | 特征行数 = 标签行数、切分无交集且覆盖全集、预测顺序对齐、对接分非空且已排序 | PASS |
| 不变量自检（右端守门） | `python scripts/99_selftest.py` | A 数据契约 / B 切分泄漏 / C 结构与口袋 / D 对接自证 / E 模型与闭环 | **19 PASS / 0 FAIL / 0 SKIP**（秒级） |
| 复现性重算 | `python scripts/99_selftest.py --full` | 追加 F 组：重算前 500 个分子 ECFP 逐位比对、同种子重算骨架切分 | **21 PASS / 0 FAIL**（15 秒） |
| 单元测试 | `python -m unittest discover -s tests -v` | 临时目录 + 8 棵树 + 400 条样本，真跑 05→06→07 | **OK**（14.3 秒） |

一键端到端：`AIDD_OFFLINE=1 bash run_all.sh` —— 预检 → 01~11 → 报告 → 产物自洽 → 不变量自检，实测 **EXIT=0，53 分 03 秒**。

> 设计原则：`99_selftest.py` 里每一条断言都对应一个真实踩过的坑（编号差 24、受体混进共晶配体、
> 随机切分"默写"、 chemprop 归一化写反、docking 与活性正交……）。文档会过时，注释会被忽略，只有断言会被执行。

---

## 目录约定

```
aidd-mvp/
├── config.yaml                 # 唯一的全局参数来源
├── run_all.sh
├── scripts/
│   ├── _common.py              # 配置、路径、vendor 引导、统一日志
│   ├── 01_target_rank.py       # 靶点多证据打分
│   ├── 02_download_structure.py# 下载 RCSB 结构
│   ├── 03_extract_pocket.py    # 口袋提取 + 编号偏移校验
│   ├── 04_fetch_activity.py    # ChEMBL 活性数据
│   ├── 05_standardize.py       # MolStandardize + pIC50 + 多测量聚合
│   ├── 06_scaffold_split.py    # Murcko 骨架切分
│   ├── 07_rf_baseline.py       # ECFP4 + RandomForest
│   ├── 08_chemprop.py          # D-MPNN 对照（同一 scaffold split）
│   ├── 09_docking.py           # AutoDock Vina（先自证，再筛分子）
│   ├── 10_admet.py             # 描述符面板 + 漏斗 + Pareto
│   └── 11_closed_loop.py       # 主动学习三策略对照
├── data/
│   ├── raw/                    # 原始下载物，永不手改
│   ├── interim/                # 中间产物（pdbqt / 对接 pose / 盒子）
│   └── processed/              # 可直接建模 / 交付的表格与结论
├── models/                     # 训练出来的权重
├── reports/
├── figures/
└── vendor/                     # 隔离安装的 biopython / meeko / gemmi
```

---

## 每一步的四问四答

| Stage | Input | 做什么 | Output | 下一步给谁 |
|---|---|---|---|---|
| 01 | 靶点证据表 | 五维加权打分 | `target_rank.csv` | 02 / 04 |
| 02 | UniProt `P00533` | 下载 PDB | `1M17.pdb` | 03 |
| 03 | PDB + 配体残基名 | 5 Å 口袋提取、编号偏移校验 | `pocket.json` | 09 |
| 04 | ChEMBL Target ID | 拉 IC50（relation `=`, nM, assay B） | `raw_activity.csv` | 05 |
| 05 | SMILES + IC50 | 标准化 → pIC50 → 中位数聚合 | `dataset.csv` | 06 |
| 06 | 数据集 | Murcko 骨架 + GroupShuffleSplit | `split_indices.json` | 07 / 08 |
| 07 | ECFP4 + 标签 | RandomForest 回归 | `predictions_*.csv` | 09 |
| 08 | 同一 split | D-MPNN（chemprop） | `08_chemprop.json` | 与 07 对照 |
| 09 | 结构 + Top-N 预测 | meeko/Vina 对接 + consensus | `docking.csv` | 10 |
| 10 | 候选 + 描述符 | Lipinski/Veber/QED + Pareto | `candidate_ranking.csv` | 湿实验 |
| 11 | 模型 + 实验回流 | 主动学习三策略 | `iteration_history.csv` | 下一轮设计 |

---

## 本机实跑结果（全部真实执行，非估算）

| 环节 | 关键数字 |
|---|---|
| 靶点排序 | EGFR Top1，且对权重扰动不敏感 |
| 结构 | PDB `1M17`，2.60 Å；5 Å 口袋 **20 个残基**，L718/G719/K745/T790/M793/C797/L844/T854/D855 全部命中 |
| **编号偏移** | 1M17 编号 = UniProt − **24**（结构里 C773 = 临床说的 C797），312/312 比对验证 |
| 活性数据 | ChEMBL `CHEMBL203`：18 381 条 → 10 468 个分子、1 910 个 assay |
| 清洗 | 0 条标准化失败；剔除 16 个离群；多测量
分子 3 710 个，其中位标准差 0.43 log |
| Baseline | ECFP4+RF：random MAE **0.510** → scaffold **0.602**（+18%） |
| Advanced | Chemprop D-MPNN 50 epoch：scaffold MAE **0.647**（略输 baseline） |
| 对接自检 | 重对接 erlotinib：晶体自打分 −6.02，最优 pose −7.38，**RMSD 1.41 Å**（<2 Å = 可信） |
| 虚拟筛选 | 30/30 成功；最优 −9.40 kcal/mol；预测 vs 对接 Pearson **r = −0.059**（两条证据正交） |
| ADMET 漏斗 | 10 468 → 活性 5 444 → +Lipinski/Veber 1 901 → +QED/剔除可疑活性 **1 798** |
| 通过率 | Lipinski 45.7% · Veber 86.7% · QED>0.3 63.5% · pIC50≥7 52.0% |
| Pareto 前沿 | **12 个**分子（0.11%）。最极端两组：pIC50 10.48 / QED 0.56，与 pIC50 4.92 / QED 0.94 —— 活性与成药性直接冲突 |
| 主动学习 | uncertainty 0.7470 最优，random 0.7777，greedy 0.8837（几乎没学到东西） |
| **转折点在很后面** | uncertainty 直到 **1 800** 个样本才反超 random；在此之前反而落后（1 200 样本时 0.8303 vs 0.7976） |

---

## 踩过的坑（每一条都实跑确认）

1. **`MolStandardize` 参数名**：是 `skipStandardization`，不是 `skipStandardize`。写错会静默返回 `None`，导致整列清洗失败。
2. **1M17 编号 ≠ 临床编号**：差 24。直接拿 pdb 里的 797 去找文献会一无所获。
3. **ChEMBL activity 端点偶发 500**：连 `limit=1` 都失败时是服务端故障，不是查询写错；脚本内置重试 + `AIDD_OFFLINE=1` 离线兜底。
4. **受体制备时务必删掉共晶配体**：留着它会导致配体-配体自碰撞，晶体 pose 自打分飙到 **+438 kcal/mol**；删掉后回到 −6.02。
5. **别用 RDKit 从 SMILES 生成的构象去做「重对接」**：它的坐标在原点附近，不在盒子里，Vina 会报 `ligand is outside the grid box`。要用 `AssignBondOrdersFromTemplate` 从晶体坐标恢复配体。
6. **chemprop 的 `ScaleTransform` 是正向标准化** `(X-mean)/scale`，不是反变换。误挂在 `FFN.output_transform` 上会让预测被二次标准化；本仓库实跑见证 MAE 从 0.60 直接飙到 12.16。
7. **不要把 pIC50 和 docking score 直接相加**：量纲不同、分布不同、可靠性不同。先 `StandardScaler` 再加权合并。
8. **公开数据里有物理上不可能的值**：pIC50 > 11 意味着 IC50 < 0.01 nM，通常源于单位/数值填写错误；不剔除的话，一个脏数据就能独占整条 Pareto 前沿。

---

## 环境说明

主依赖装在 `aidd` conda 环境；`biopython` / `meeko` / `gemmi` 用 `--target ./vendor` 隔离安装在工程内，
由 `scripts/_common.py` 在最前面插入 `sys.path`，**不污染原环境**。

```bash
python -m pip install --no-deps --target ./vendor biopython meeko gemmi
```

> 本文所有数字均来自本机实跑（MacBook Pro / M1 Pro），用于教学演示，不构成任何投资或用药建议。

## 文章与自动验证（2026-10-17 补充）

优化后的公众号文章在 `articles/公众号文章_优化版.html`，可编辑版本为同名 `.md`，正文包含核心脚本。环境版本见 `requirements.txt`。

```bash
python -m unittest discover -s tests -v
AIDD_TEST_FULL=1 python -m unittest discover -s tests -v
python scripts/check_artifacts.py
```

测试在临时目录运行，避免覆盖全量结果。全阶段冒烟使用低预算；它确认执行路径，不能替代模型质量与重对接质量验证。本次验证范围与限制见 `reports/validation.md`。原文“全部真实执行”是历史结果说明，本次未重跑全预算。`run_all.sh` 已修正为执行 01—11，并生成报告和验收产物。
