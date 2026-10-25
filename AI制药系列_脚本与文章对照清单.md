# AI 制药系列 · 脚本 ↔ 文章对照清单

> 生成时间：2026-10-25　｜　工作区：`/Users/dujiepeng/WorkBuddy/Work`
> 用途：回答三件事——**哪篇用了哪个脚本**、**哪张图是谁画的**、**哪些资产已经断链**。
> 统计口径：脚本说明取自各文件头部 docstring，配图归属用「脚本内出现的输出文件名」+「文件 mtime」双向比对确认。

---

## 0. 先看清全局：三层结构

```
   探索层  aidd_practice/        ← 试跑、验证素材、画配图（一次性，允许乱）
      │      产出：*.json 结论 + *配图/*.png
      ▼
   工程层  aidd-mvp/             ← 11 步流水线，config.yaml 驱动，可复跑、可自检
      │      产出：data/processed/* + reports/summary.md
      ▼
   文章层  *.html                ← 只消费上面两层的产物，本身不计算
```

**核心原则（系列一贯）**：文章里每个数字都必须能追溯到某个脚本的实跑输出。所以这份清单本质上是文章的**溯源表**——改任何一个数字，先回到这里找到它的来源脚本。

---

## 1. 系列文章索引

| # | 篇名 | 文章文件 | 主要脚本 | 配图目录 |
|---|---|---|---|---|
| ① | 起始篇·26 周学习计划 | `AI制药_公众号文章.html` | 无（规划类，无实跑） | — |
| ② | 第 1 课·读懂一个分子 | `AI制药第一课_公众号文章.html` | `aidd_ecfp_demo.py`（顶层） | — |
| ②附 | 30 行代码预测溶解度 | `DeepChem入门_公众号文章.html` | 无本地脚本留存 | — |
| ⑤ | 工具篇·RDKit→DeepChem | `AI制药工具篇_公众号文章.html` | **脚本已丢失** | `AI制药工具篇_配图/` |
| ⑥上 | 实践篇·七站流水线导论 | `AI制药实践篇_公众号文章.html` | `make_figures*.py`（1–5）、`pipeline_demo.py` | `AI制药实践篇_配图/` |
| ⑥下 | 实践篇 B·跑通 MVP | `AI制药实践篇B_实操篇_公众号文章.html` | `egfr_mvp_step1-5`、`mvp_trap_demo`、`make_figures_mvp`、`redraw_*` | `AI制药实操篇_配图/` |
| ⑦ | 行业调研篇 | `AI制药行业调研篇_公众号文章.html` | `make_figures_market`、`make_figures_source`、`make_figures_forecast` | `AI制药行业调研篇_配图/` |
| ⑧ | 工程篇·11 个脚本变流水线 | `AI制药工程篇_公众号文章.html` | `aidd-mvp/scripts/*`（01–11 + 00 + 99）、`make_figures_eng` | `AI制药工程篇_配图/` |
| ⑨ | 细节篇①·靶点与口袋 | `AI制药细节篇1_靶点与口袋_公众号文章.html` | `af_demo/e1,e3,e4,e5,e6,e7` + `make_figures_af` | `AI制药细节篇1_配图/` |
| ⑩ | 细节篇②·DiffDock（**待写**） | 无 | `diffdock_demo/d1,d2,d3` + `make_figures_lesson4` | `aidd_practice/AI制药第4课_配图/` |

> ③ 第 2 课（四层生化语义）、④ 政策篇（十五五）已发布但**本地未留存 HTML**，本清单不覆盖。

---

## 2. 逐篇对照（含注解）

### ⑤ 工具篇　`AI制药工具篇_公众号文章.html`

| 脚本 | 作用 | 产出 | 注解 |
|---|---|---|---|
| `aidd_ecfp_demo.py`（顶层） | ECFP4 指纹相似度演示 | `ecfp_similarity_matrix.png` | 顶层遗留脚本，② 与 ⑤ 共用 |
| ❓ **已丢失** | 画 fig1–fig5 | `fig1_tool_map.png` 等 5 张 | **全工作区搜索无生成脚本**（配图 09-27 12:01，脚本未留存）。要改这 5 张图只能重画 |

**断链提醒**：本篇 5 张正文配图无法复现。若日后需要改图，需按文章数据重新写绘图脚本。

---

### ⑥上 实践篇导论　`AI制药实践篇_公众号文章.html`

| 脚本 | 作用 | 产出（→ 文章图号） |
|---|---|---|
| `pipeline_demo.py` | PubChem 拉真实分子 + RDKit 计算，给正文提供可引用数字 | 无图，供数字 |
| `make_figures.py` | 首版配图 | fig2 `fig2_hit_to_lead`、fig3 `fig3_dry_wet_loop`、封面 `cover_500/900` |
| `make_figures2.py` | v2 新增 | fig6 `fig6_roadmap` |
| `make_figures3.py` | v3 **修** fig1/fig4/fig5 错乱 | fig4 `fig4_seams`、fig5 `fig5_two_wars` |
| `make_figures4.py` | v4 重画 | fig1 `fig1_seven_stations`、fig7 `fig7_three_layers` |
| `make_figures5.py` | v5，数据来自 `lora_loop.py` | fig8 `fig8_dmta_lora`、fig9 `fig9_lora_tradeoff` |
| `lora_loop.py` | 干湿闭环对照：不更新 / 全参微调 / LoRA | 供 fig8/fig9 数据 |
| `al_loop.py` | 主动学习 vs 随机选样（TDC Caco-2） | 供正文数字 |
| `covalent_warhead.py` | 共价弹头 SMARTS 过滤（KRAS 站） | 供正文数字 |

**注解**：本篇是「导论」，脚本特征是**多版迭代**（make_figures 1→5），每版修上一版的排版/口径问题。看数字出处时务必认准最终版（见 §4）。

---

### ⑥下 实践篇 B·实操　`AI制药实践篇B_实操篇_公众号文章.html`

| 脚本 | 作用 | 产出 | 注解 |
|---|---|---|---|
| `egfr_mvp_numbering.py` | 核实 1M17 PDB 编号 ↔ UniProt P00533 偏移 | 结论：偏移 **+24** | 后面所有脚本的地基 |
| `egfr_mvp_step1.py` | 靶点排序 + 结构 + 口袋提取 | `pocket.json`、`target_rank.csv` | |
| `egfr_mvp_step2.py` | ChEMBL 拉真实 IC50（带重试/断点续传） | `raw_activity.csv` | |
| `egfr_mvp_step3.py` | 清洗 → pIC50 → 骨架切分 → ECFP+RF | `dataset.csv`、`baseline_results.json` | |
| `egfr_mvp_step4.py` | ADMET → Pareto → 主动学习 | `al_history.json`、`dataset_with_props.csv` | |
| `egfr_mvp_step5.py` | 手写轻量 GCN 作对照 | `gcn_results.json` | Chemprop 需独立环境，本机未装 |
| `make_figures_mvp.py` | 实操篇主配图（334 行） | m1 `m1_mvp_map`、m3 `m3_pocket` | **m2/m4/m5/m6 已被重画覆盖** |
| `redraw_m_series.py` | 用 aidd-mvp 工程口径重画 | m2/m4/m5/m6 | 注解：旧版基于早期试跑，口径与最终工程不一致 |
| `mvp_trap_demo.py` | **验证流传伪代码的三处硬伤** | t1 `t1_box`（写死口袋坐标差 31.4 埃）、t2/t3 初版 | 437 行，本篇最重要脚本 |
| `redraw_t2_t3.py` | 重绘 t2/t3（修文字截顶） | t2、t3 | |
| `mvp_generation_sa.py` | BRICS 分解重组，验证 SAScore | t4 初版 | 枚举爆炸极慢 |
| `mvp_generation_sa_fast.py` | **SMILES 字符串重组快速版** | `t4_generation_sa.png` | ✅ 最终采用版 |
| `make_figures_traps.py` | 六坑地图 | `trap_map_A.png`、`trap_map_B.png` | |
| `mvp_target_sens.py` | 靶点排序权重敏感性（4 套权重） | `mvp_target_sensitivity.json` | EGFR 均第 1、KRAS 2↔4 洗牌 |

---

### ⑦ 行业调研篇　`AI制药行业调研篇_公众号文章.html`

| 脚本 | 作用 | 产出 |
|---|---|---|
| `make_figures_market.py` | 主配图（282 行，数据来自报告正文引用） | fig1–fig5 |
| `make_figures_source.py` | 来源图核对后的**高清重绘** | fig6 `fig6_traditional_vs_ai`、fig7 `fig7_tech_roadmap`、fig8 `fig8_value_chain` |
| `make_figures_forecast.py` | 五机构市场规模口径校准 | fig9 `fig9_global_forecast` |
| `make_cover_market.py` | 封面 1920×817 | `AI制药行业调研篇_公众号封面_1920x817.jpg` |

**注解**：
- `fig0_overview.jpg`（09-26）**无生成脚本**，是外部来源原图直接引用（已按「原图能否直接插入三判据」核过）。
- 配图目录里还有 `来源核对/` 子目录，是原始素材留存。
- 本篇数字来自**报告引用**而非本机计算，与其余各篇的「实跑」性质不同——引用时需注意标注来源。

---

### ⑧ 工程篇　`AI制药工程篇_公众号文章.html`

工程本体是 `aidd-mvp/`，不是 `aidd_practice/`。

| 文件 | 作用 | 注解 |
|---|---|---|
| `aidd-mvp/config.yaml` | 单一参数源 | 所有脚本从这里读参数，禁止硬编码 |
| `scripts/_common.py` | 公共工具：配置加载、路径、日志 | 每个脚本开头 `from _common import ROOT, cfg, log` |
| `scripts/00_preflight.py` | 环境预检，必需/可选分级 | 缺必需 FAIL（退出 1），缺可选只 WARN |
| `scripts/01_target_rank.py` | 靶点排序 | → `target_rank.csv` |
| `scripts/02_download_structure.py` | 下载实验结构 + 质检 | → `1M17.pdb` |
| `scripts/03_extract_pocket.py` | 口袋残基 + **编号偏移校验** | 本项目最大的坑都在这 |
| `scripts/04_fetch_activity.py` | ChEMBL 拉真实活性 | → `raw_activity.csv` |
| `scripts/05_standardize.py` | 清洗 → 去重 → pIC50 → 聚合 | MolStandardize 参数名是 `skipStandardization` |
| `scripts/06_scaffold_split.py` | 骨架切分 + ECFP 矩阵 | → `X_ecfp.npy`、`split_indices.json` |
| `scripts/07_rf_baseline.py` | ECFP4 + RF 基线 | **random MAE 0.5095 / scaffold 0.6015** |
| `scripts/08_chemprop.py` | D-MPNN 对照 | scaffold 0.6474 |
| `scripts/09_docking.py` | Vina 对接（280 行） | 重对接 RMSD **1.411 埃** |
| `scripts/10_admet.py` | ADMET + Pareto | 漏斗 10468→5444→1901→1798，Pareto 12 |
| `scripts/11_closed_loop.py` | 干湿闭环 | uncertainty 0.7470 < random 0.7777 |
| `scripts/99_selftest.py` | **把踩过的坑固化成断言**（19 项 / `--full` 21 项） | 与本清单 §2 的「注解」互补：清单给人看，断言给机器看 |
| `scripts/make_report.py` | 汇总所有 json → `reports/summary.md` | |
| `scripts/check_artifacts.py` | 只读检查现有产物 | 明确「检查 ≠ 重训」 |
| `tests/test_pipeline.py` | unittest 真跑 | 环境无 pytest，用 `python -m unittest discover -s tests` |
| `run_all.sh` | 一键跑通（实测 53m03s，EXIT=0） | |
| `make_figures_eng.py` | 工程篇配图 2 张 | `e1_structure.png`（11 脚本+三层数据+两端守门）、`e2_guard.png`（坑→断言映射） |

---

### ⑨ 细节篇①·靶点与口袋　`AI制药细节篇1_靶点与口袋_公众号文章.html`

目录 `aidd_practice/af_demo/`，是**目前组织得最整齐的一组**（e 系列正式脚本 + diag 系列诊断脚本分离）。

#### 正式脚本（e 系列，全部有结论产出）

| 脚本 | 作用 | 产出 | 对应文章位置 |
|---|---|---|---|
| `e1_af_vs_crystal.py` | AF vs 1M17 叠合分层 + pLDDT 能否当口袋探测器 | `results/e1_af_vs_crystal.json` | 第四节核心数据 |
| `e3_pocket_detection.py` | 手写几何口袋检测（LIGSITE 式埋藏度 + 6 连通聚类） | `results/e3_pocket_detection.json` | 第四节表：7983 / 26869 / 10191 立方埃 |
| `e4_visualize.py` | 输出「能看」的产物：PDB + PyMOL/ChimeraX 脚本 + 单文件交互 HTML | `viz/*.pdb`、`EGFR_AFPocket_3D.html`、`a6_3d_view.png` | **未入正文**（微信不支持 JS） |
| `e5_py3dmol.py` | py3Dmol 三段递进出图 | a7 `a7_py3dmol_hello`、a8 `a8_py3dmol_plddt`、a9 `a9_py3dmol_pocket` | 第六节第 4 步 |
| `e6_end2end.py` | **端到端闭环**：下载→裁剪→找口袋→出图，单文件可跑 | `results/e6_end2end.json` | 第六节四段代码即其核心片段 |
| `e7_rigor_check.py` | 口径审查补算（RMSD 四组对照 + 覆盖率链条） | `results/e7_rigor_check.json` | 因审稿意见新增 |
| `make_figures_af.py` | 本篇正文配图 | a1 `a1_folding`、a2 `a2_plddt_profile`、a3 `a3_rmsd_layers`、a4 `a4_plddt_rank`、a5 `a5_pocket_detection` | 图 1–5 |
| `make_cover_af.py` | 公众号封面 + 小红书封面 | `AI制药细节篇1_公众号封面_1920x817.jpg`、`..._小红书封面_1080x1440.jpg` | |
| `verify_article.py` | 交付前校验：标签闭合、图片存在、图号连续、表格宽度、禁用写法 | 无产出（校验器） | 每次改稿必跑 |

#### 诊断脚本（diag 系列，**任务完成，可归档**）

| 脚本 | 当时在查什么 | 结论 |
|---|---|---|
| `diag.py` | 序列比对实测 PDB↔UniProt 偏移 | offset = **+24**，312/312 全中 |
| `diag2.py` | 为何 common 残基为空 + 扫口袋检测参数 | 定位到 v4 残骸文件 + 参数区间 |
| `diag3.py` | 合成数据验证 Kabsch | 发现行向量要用 `P @ R.T + t` |
| `diag4.py` | scipy 独立实现交叉验证 | 确认全长 RMSD 11.675 埃非算法 bug |

**注解**：
- `e2` 编号已并入 `e1`（E1 文件里同时做两件事），所以没有独立 e2 文件——不是缺失。
- 本篇最重要的教训：**引用自己旧稿的数字要回 json 核对，别抄正文**。曾把「交集体积 453」抄成「配体体积 455」，导致覆盖率分母写错。

---

### ⑩ 细节篇②·DiffDock（素材已备，尚未成篇）

目录 `aidd_practice/diffdock_demo/`。

| 脚本 | 作用 | 产出 / 关键数字 |
|---|---|---|
| `d1_blindspots.py` | ECFP4 三个盲区：手性 / 构象 / 2D-3D 排序 | 对映体 ECFP4 相似度 **1.0000**；厄洛替尼 14 构象 RMSD 3.257 埃；Spearman ρ=0.0303 |
| `d2_vina_real_activity.py` | Vina 分数 vs 真实 pIC50 全面体检 | **ρ=−0.4984**（p=4.29e−07）；按 Vina 选 30 只命中 14 |
| `d3_box_sensitivity.py` | 三盒敏感性：正确口袋 / 偏移 8 埃 / 全蛋白盲搜 | 证明 Vina 分数是「搜索区域的函数」 |
| `make_figures_lesson4.py` | 配图 | c1 `c1_blindspots`、c2 `c2_vina_activity`、c4 `c4_positioning` |

**注解**：
- 本机**跑不通 DiffDock 全栈**（macOS arm64 无 torch_scatter/torch_cluster wheel + fair-esm 2.5GB + PyPI 53KB/s），故本篇顺延。
- 目录里还有 `venv/`（10-24 10:40）是当时的尝试残留，**可以安全删除**。
- 图号 c3 缺失（只有 c1/c2/c4），成篇时需补或重编号。

---

## 3. 断链与孤儿资产（需要留意）

| 类型 | 对象 | 说明 | 建议 |
|---|---|---|---|
| 🔴 断链 | 工具篇 5 张配图 | 全工作区无生成脚本（09-27） | 如需改图只能重画；建议补一个 `make_figures_tool.py` 留档 |
| 🟡 无脚本 | `AI制药行业调研篇_配图/fig0_overview.jpg` | 外部来源原图直接引用 | 已在 `来源核对/` 留素材，保持现状 |
| 🟡 无脚本 | `AI制药实操篇_配图/trap_overview_card.jpg`（10-17 20:34） | 晚于所有脚本，疑为后期手工合成 | 若需复用，记录合成步骤 |
| 🟢 未入正文 | `AI制药细节篇1_配图/a6_3d_view.png` | 由 `e4_visualize.py` 生成，作 3D 交互页的静态预览 | 微信不支持 JS，交互页只能外部链接 |
| 🟢 未入正文 | `EGFR_AFPocket_3D.html`（345 KB） | 单文件 4 视图交互页 | 可作为「阅读原文」链接 |
| ⚪ 可归档 | `af_demo/diag*.py`（4 个） | 诊断任务已完成，结论已沉淀进 e1/e3 与记忆 | 保留作溯源，不必再跑 |
| ⚪ 可清理 | `diffdock_demo/venv/` | DiffDock 环境尝试失败残留 | 可删 |
| ⚪ 未成篇 | `AI制药系列_第3课_分子性质预测_内容规划.md` | 早期规划，已被 ⑨ 取代 | 留作选题参考 |

---

## 4. 版本陷阱：同名图被多次重绘，哪个才是最终版？

绘图脚本迭代时**覆盖同名文件**，容易搞错生效版本。按 mtime 比对确认如下：

| 配图 | 候选脚本 | ✅ 最终生效 | 判据 |
|---|---|---|---|
| 实操篇 m2/m4/m5/m6 | `make_figures_mvp.py` vs `redraw_m_series.py` | **`redraw_m_series.py`** | 图 mtime 10-17 14:38 = 该脚本 mtime；脚本自述「旧版口径与工程不一致」 |
| 实操篇 m1/m3 | 仅 `make_figures_mvp.py` | `make_figures_mvp.py` | 无竞争者 |
| 实操篇 t2/t3 | `mvp_trap_demo.py` vs `redraw_t2_t3.py` | **`redraw_t2_t3.py`** | 图 mtime 10-17 14:28 |
| 实操篇 t1 | 仅 `mvp_trap_demo.py` | `mvp_trap_demo.py` | 图 mtime 10-17 12:55 |
| 实操篇 t4 | `mvp_generation_sa.py` vs `_fast.py` | **`mvp_generation_sa_fast.py`** | BRICSBuild 太慢，fast 版为最终采用；但图 mtime 18:18 晚于脚本，疑有末次微调 |
| 实践篇 fig1 | v1/v3/v4 | **`make_figures4.py`** | 图 mtime 21:09 |
| 实践篇 fig4/fig5 | v2/v3 | **`make_figures3.py`** | 图 mtime 20:45，v3 专为修这两张 |
| 实践篇 fig6 | 仅 v2 | `make_figures2.py` | 图 mtime 20:27 |
| 实践篇 fig2/fig3/封面 | 仅 v1 | `make_figures.py` | 图 mtime 20:01 |

**教训**：以后绘图脚本带版本迭代时，在文件头写清「本版覆盖哪几张图」，否则三个月后只能靠 mtime 考古。

---

## 5. 复现顺序（从零到全部产物）

### ⑨ 细节篇①（组织最清晰，推荐先看这套）

```bash
cd aidd_practice/af_demo
python e6_end2end.py        # 端到端：下载 → 裁剪 → 找口袋 → 出图（单文件可跑）
python e1_af_vs_crystal.py  # 分层 RMSD + pLDDT 探测口袋
python e3_pocket_detection.py
python e7_rigor_check.py    # 口径审查
python make_figures_af.py   # 正文图 1-5
python e5_py3dmol.py        # 正文图 7-9（需 Chrome headless 截图）
python e4_visualize.py      # 3D 交互页 + 图 6
python make_cover_af.py     # 封面
python verify_article.py    # 交付前校验
```

### ⑧ 工程篇（aidd-mvp）

```bash
cd aidd-mvp
bash run_all.sh                        # 一键跑通，约 53 分钟
python scripts/99_selftest.py --full   # 21 项断言
python -m unittest discover -s tests   # 注意：环境无 pytest
```

### ⑥下 实操篇

```bash
cd aidd_practice
python egfr_mvp_numbering.py   # 先确认编号偏移
python egfr_mvp_step1.py ... egfr_mvp_step5.py
python mvp_trap_demo.py        # 三处硬伤验证
python mvp_generation_sa_fast.py
python make_figures_mvp.py     # m1/m3
python redraw_m_series.py      # m2/m4/m5/m6
python redraw_t2_t3.py
python make_figures_traps.py
```

---

## 6. 维护建议

1. **改数字先查表**：文章里任何带单位的数字，先从本清单找到来源脚本，改脚本 → 重跑 → 更新 json → 再改正文，不要直接改 HTML。
2. **新增绘图脚本必写头部 docstring**：说明本版覆盖哪几张图、数据来源是哪个 json。§4 的考古工作本可以避免。
3. **诊断脚本与正式脚本分目录**：`af_demo/` 的 `e*`（正式）+ `diag*`（诊断）命名是目前最好的实践，建议推广到 `aidd_practice/` 其他篇。
4. **交付前必跑 `verify_article.py`**：校验标签闭合、图号连续、图片存在、表格宽度合计 100%、微信禁用写法（`<figure>`、负 margin）。
5. **补一个 `make_figures_tool.py`**：把工具篇 5 张图的绘制逻辑补录进仓库，消除唯一的硬断链。
