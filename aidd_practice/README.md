# AI 制药脚本与公众号文章索引

本目录是文章的实验与制图工作区，公众号 HTML 和多数配图目录位于上一级 `Work/`。正式可复跑工程在 [`../aidd-mvp/`](../aidd-mvp/README.md)。本索引根据本地文章的章节、图片引用和脚本的输入输出整理；“素材关联”不表示正文直接引用了脚本名，也不表示文章已发布。

## 按文章查脚本

### 实践篇：七站地图与闭环

正文：[AI制药实践篇_公众号文章.html](../AI制药实践篇_公众号文章.html)。

| 脚本 | 文章中的用途 | 输入、产物与注解 |
| --- | --- | --- |
| [pipeline_demo.py](pipeline_demo.py) | 实验：苗头到先导的分子性质、相似性与过滤演示 | 结果输出到终端；真实分子案例用于解释流水线，不构成药效验证。 |
| [al_loop.py](al_loop.py) | 实验：第七站，主动学习与随机选样对照 | TDC Caco-2 真实标签用于离线模拟回传；没有执行湿实验。 |
| [lora_loop.py](lora_loop.py) | 实验：第七站，DMTA 增量更新与 LoRA 对照 | TDC Lipophilicity 骨架切分；配图 make_figures5.py 使用手工录入的实测摘要，不自动读取本次输出。 |
| [covalent_warhead.py](covalent_warhead.py) | 实验：第四站及 KRAS 案例，共价弹头 SMARTS 过滤 | SMARTS 命中表示结构片段匹配，不能据此确认共价结合或真实活性。 |
| [make_figures.py](make_figures.py) | 配图：早期七站图、苗头到先导图、干湿闭环图 | 历史增量版本，不是五套完整图；同名图会被后运行的版本覆盖。 |
| [make_figures2.py](make_figures2.py) | 配图：v2：七站回环、接缝、双案例、学习路径 | 历史增量版本，不是五套完整图；同名图会被后运行的版本覆盖。 |
| [make_figures3.py](make_figures3.py) | 配图：v3：修订 fig1 / fig4 / fig5 | 历史增量版本，不是五套完整图；同名图会被后运行的版本覆盖。 |
| [make_figures4.py](make_figures4.py) | 配图：v4：修订 fig1，新增 fig7 三层支撑 | 历史增量版本，不是五套完整图；同名图会被后运行的版本覆盖。 |
| [make_figures5.py](make_figures5.py) | 配图：v5：fig8 DMTA / LoRA 与 fig9 更新代价 | 历史增量版本，不是五套完整图；同名图会被后运行的版本覆盖。数值为代码内实测摘要，重跑实验后须人工同步。 |

### 实践篇 B：EGFR 实操

正文：[AI制药实践篇B_实操篇_公众号文章.html](../AI制药实践篇B_实操篇_公众号文章.html)。

| 脚本 | 文章中的用途 | 输入、产物与注解 |
| --- | --- | --- |
| [egfr_mvp_step1.py](egfr_mvp_step1.py) | 早期实验 1：靶点排序、1M17 与口袋提取 | 写入 data/egfr_mvp/；属于正式 aidd-mvp 工程之前的试跑口径。 |
| [egfr_mvp_step2.py](egfr_mvp_step2.py) | 早期实验 2：ChEMBL EGFR IC50 数据抓取 | 输出 raw_activity.csv，带缓存和断点；接早期 step3，不能与正式工程数据目录混用。 |
| [egfr_mvp_step3.py](egfr_mvp_step3.py) | 早期实验 3：标准化、pIC50、骨架切分及 RF 基线 | 读取 step2，输出数据表、指纹及基线摘要；与正式工程的最终指标需区分。 |
| [egfr_mvp_step4.py](egfr_mvp_step4.py) | 早期实验 4：性质面板、Pareto 与主动学习 | 依赖 step3 的数据表及指纹；性质计算与离线主动学习不等于实验验证。 |
| [egfr_mvp_step5.py](egfr_mvp_step5.py) | 早期扩展：GCN 与 ECFP / RF 表示对照 | 依赖 step3；不是本篇正式流水线的第五阶段，文章当前章节不能按文件编号直接对应。 |
| [egfr_mvp_numbering.py](egfr_mvp_numbering.py) | 事实核查：PDB 与 UniProt 残基编号偏移 | 依赖 data/egfr_mvp/1M17.pdb；为口袋编号说明提供校验依据。 |
| [make_figures_mvp.py](make_figures_mvp.py) | 早期配图：m1–m6 | 读取 data/egfr_mvp/；最终 m2 / m4 / m5 / m6 由 redraw_m_series.py 按正式工程口径重画。 |
| [redraw_m_series.py](redraw_m_series.py) | 正式口径配图：m2 / m4 / m5 / m6 | 读取 ../aidd-mvp/data/raw 与 processed；应在 make_figures_mvp.py 之后运行以保留最终口径。 |
| [mvp_trap_demo.py](mvp_trap_demo.py) | 反例实验与配图：硬编码盒子、演示库、SA 与 docking | 读取正式 aidd-mvp 的配置和产物，写 mvp_trap_results.json 及 t1–t3。 |
| [mvp_target_sens.py](mvp_target_sens.py) | 补充实验：靶点排序权重敏感性 | 读取正式工程 target_rank.csv，输出 mvp_target_sensitivity.json。 |
| [mvp_generation_sa.py](mvp_generation_sa.py) | 补充实验：BRICS 玩具生成与 SAScore 对照 | 写真实分子属性缓存、mvp_generation_sa.json 及 t4；与 fast 版共享输出名但生成方法不同。 |
| [mvp_generation_sa_fast.py](mvp_generation_sa_fast.py) | 补充实验：SMILES 拼接玩具生成对照 | 先由 mvp_generation_sa.py 准备 mvp_real_props.csv；会覆盖同名 JSON 和 t4，不可混称 BRICS 结果。 |
| [redraw_t2_t3.py](redraw_t2_t3.py) | 版式修订：t2 / t3 | 依赖 mvp_trap_results.json 及正式工程数据；覆盖 mvp_trap_demo.py 生成的同名图。 |
| [make_figures_traps.py](make_figures_traps.py) | 配图：六坑地图 A / B | 读取正式工程数据及 mvp_trap_results.json；不重算实验结论。 |
| [make_cover_practice_b.py](make_cover_practice_b.py) | 封面：公众号初版 | 与 v2 写入同一 JPG；历史初版应先运行，避免覆盖 v2。 |
| [make_cover_practice_b_v2.py](make_cover_practice_b_v2.py) | 封面：公众号第二版 | 与初版共享输出名；文件名 v2 表示生成器版本，产物名不带 v2。 |
| [make_xhs_cover_practice.py](make_xhs_cover_practice.py) | 封面：小红书竖版 | 本篇跨平台封面，不是实验或正文配图。 |

### 工程篇：可复跑与自检

正文：[AI制药工程篇_公众号文章.html](../AI制药工程篇_公众号文章.html)。

| 脚本 | 文章中的用途 | 输入、产物与注解 |
| --- | --- | --- |
| [make_figures_eng.py](make_figures_eng.py) | 配图：e1 工程结构 / e2 坑与断言 | 真正的流水线、预检和自检在 ../aidd-mvp/；本脚本只绘制说明图。 |

### 行业调研篇

正文：[AI制药行业调研篇_公众号文章.html](../AI制药行业调研篇_公众号文章.html)。

| 脚本 | 文章中的用途 | 输入、产物与注解 |
| --- | --- | --- |
| [make_figures_market.py](make_figures_market.py) | 配图：fig1–fig5 市场、阶段、命中率、技术与巨头 | 图表数据写在代码内，来自文章研究素材；重绘不代表重新获取最新行业数据。 |
| [make_figures_source.py](make_figures_source.py) | 配图：来源重绘 fig6 / fig7 / fig8 | 当前文章引用 fig7；fig6 / fig8 为保留素材，不能把所有生成图都算作当前正文插图。 |
| [make_figures_forecast.py](make_figures_forecast.py) | 配图：fig9 市场预测口径对照 | 数据为代码内的历史研究摘要；当前网页引用 fig9_global_forecast.png。 |
| [make_cover_market.py](make_cover_market.py) | 封面：行业调研公众号横版 | 封面独立于正文 img 引用。 |

### 细节篇第 1 课：靶点与口袋

正文：[AI制药细节篇1_靶点与口袋_公众号文章.html](../AI制药细节篇1_靶点与口袋_公众号文章.html)。

| 脚本 | 文章中的用途 | 输入、产物与注解 |
| --- | --- | --- |
| [af_demo/e1_af_vs_crystal.py](af_demo/e1_af_vs_crystal.py) | 实验 E1 / E2：AF 与晶体偏差、pLDDT 与口袋 | 读取 struct/egfr_af_v6.pdb 和正式工程 1M17；输出 results/e1_af_vs_crystal.json。 |
| [af_demo/e3_pocket_detection.py](af_demo/e3_pocket_detection.py) | 实验 E3：几何口袋检测对照 | 复用 E1 解析与编号常量；输出 e3_pocket_detection.json；pLDDT 不是口袋标签。 |
| [af_demo/e4_visualize.py](af_demo/e4_visualize.py) | 实验 E4：结构、PyMOL / ChimeraX 与交互网页 | 复用 E1 / E3 算法；输出 viz/ 与上一级 EGFR_AFPocket_3D.html，提供补充可视化。 |
| [af_demo/e5_py3dmol.py](af_demo/e5_py3dmol.py) | 实验 E5：三段可视化示例与 a7–a9 配图 | 依赖 py3Dmol、Chrome 截图及已有结构文件；Chrome 路径写死为 macOS 安装位置。 |
| [af_demo/e6_end2end.py](af_demo/e6_end2end.py) | 正文第六节入口：下载、裁剪、找口袋、交互展示 | 独立小闭环，输出 struct/、results/e6_end2end.json 和 viz/e6_end2end_view.html；不自动产出全部正文图。 |
| [af_demo/e7_rigor_check.py](af_demo/e7_rigor_check.py) | 补算：RMSD 残基集合、叠合方式与覆盖率口径 | 复用 E1；输出 e7_rigor_check.json，供 make_figures_af.py 绘制偏差图。 |
| [af_demo/make_figures_af.py](af_demo/make_figures_af.py) | 正文配图：a1–a5 | 读取 E1 / E3 / E7 结果及结构；先完成实验再绘图，a7–a9 另由 E5 生成。 |
| [af_demo/make_cover_af.py](af_demo/make_cover_af.py) | 封面：公众号横版与小红书竖版 | 使用真实 AF Cα 坐标投影；依赖本地结构，不是正文实验。 |
| [af_demo/verify_article.py](af_demo/verify_article.py) | 交付校验：文章 HTML | 检查文章结构与引用等，不替代数值科学审查；执行前确认目标路径。 |
| [af_demo/diag.py](af_demo/diag.py) | 开发诊断：编号、口袋参数、Kabsch 与独立交叉验证 | 保留排错过程；不是正文推荐入口，不需要为了阅读文章依次执行。 |
| [af_demo/diag2.py](af_demo/diag2.py) | 开发诊断：编号、口袋参数、Kabsch 与独立交叉验证 | 保留排错过程；不是正文推荐入口，不需要为了阅读文章依次执行。 |
| [af_demo/diag3.py](af_demo/diag3.py) | 开发诊断：编号、口袋参数、Kabsch 与独立交叉验证 | 保留排错过程；不是正文推荐入口，不需要为了阅读文章依次执行。 |
| [af_demo/diag4.py](af_demo/diag4.py) | 开发诊断：编号、口袋参数、Kabsch 与独立交叉验证 | 保留排错过程；不是正文推荐入口，不需要为了阅读文章依次执行。 |

### 第 4 课储备素材：2D / 3D 与对接

未在相邻目录找到对应第 4 课正文，以下按脚本说明与产物归组，暂不绑定到某篇已存在文章。目录名 `diffdock_demo` 不代表已经运行 DiffDock；现有实验主要是 RDKit 与 Vina。

| 脚本 | 文章中的用途 | 输入、产物与注解 |
| --- | --- | --- |
| [diffdock_demo/d1_blindspots.py](diffdock_demo/d1_blindspots.py) | 第4课素材：2D 指纹的手性、构象与形状盲区 | 输出 results/d1_blindspots.json 等；不运行 DiffDock。 |
| [diffdock_demo/d2_vina_real_activity.py](diffdock_demo/d2_vina_real_activity.py) | 第4课素材：Vina 与真实活性的对照 | 读取正式工程结构和活性数据，执行 Vina；输出 d2 摘要与逐分子 CSV，不运行 DiffDock。 |
| [diffdock_demo/d3_box_sensitivity.py](diffdock_demo/d3_box_sensitivity.py) | 第4课素材：正确、偏移与盲搜盒子对照 | 依赖 D2 的 CSV；当前目录未发现 d3_box_sensitivity.json，绘图前需完成 D3。 |
| [diffdock_demo/make_figures_lesson4.py](diffdock_demo/make_figures_lesson4.py) | 第4课配图：盲区、排序、盒子与工具定位 | 读取 D1 / D2 / D3 JSON，输出本项目 AI制药第4课_配图/；相邻目录未找到对应第4课正文。 |

## 运行关系与图片覆盖顺序

以下是依赖关系，执行前请先阅读各脚本的路径与环境设置。许多脚本在导入时直接执行计算或写文件，不宜通过 `import` 来查看用途。

- **实践篇实验**：`pipeline_demo.py`、`covalent_warhead.py`、`al_loop.py`、`lora_loop.py` 可分别阅读和运行；后两者使用真实数据做离线回放，不涉及真实合成与湿实验。`run_output.txt` 是已有运行记录，不是所有脚本的统一结果源。
- **实践篇配图**：按 `make_figures.py → make_figures2.py → make_figures3.py → make_figures4.py → make_figures5.py` 的版本顺序可保留最新修订。最终图源：fig1 / fig7 来自第4版，fig2 / fig3 来自初版，fig4 / fig5 来自第3版，fig6 来自第2版，fig8 / fig9 来自第5版。图中部分数据直接写在代码中，重新训练不会自动更新图片。
- **早期 EGFR 实验**：`step1 → numbering` 核查结构与编号；`step2 → step3 → step4` 负责活性、建模、性质与闭环，`step3 → step5` 是额外 GCN 对照。编号脚本不是其他步骤的程序前置调用。所有这些试跑数据保存在 `data/egfr_mvp/`。
- **实操篇最终图**：`make_figures_mvp.py` 先生成 m1–m6；`redraw_m_series.py` 再以正式工程数据覆盖 m2 / m4 / m5 / m6。m1 / m3 仍来自早期绘图脚本；修改结构或口袋后需要同步检查 m3。`mvp_trap_demo.py → redraw_t2_t3.py` 生成并修订 t1–t3；`mvp_trap_demo.py → make_figures_traps.py` 生成六坑地图。
- **生成与 SA 对照**：BRICS 版和 SMILES 拼接版会写同一个 `mvp_generation_sa.json` 和 `t4_generation_sa.png`。快速版读取真实属性缓存，不能把两种方法的结果混为同一个实验；要比较须分别保留产物。现有同名文件的生产方法应结合 JSON 内容确认，不能只凭文件名判断。
- **细节篇数值与插图**：准备 `af_demo/struct/egfr_af_v6.pdb` 及正式工程 1M17 后，运行 E1、E3、E7，再运行 `make_figures_af.py` 得到 a1–a5。E3 / E7 复用 E1 的函数，不要求先读取 E1 JSON；制图脚本则需要这些结果文件。E4 提供额外结构展示；E5 生成正文 a7–a9。E6 是第六节的独立入口，不代替其余实验与全部配图。最后用 `verify_article.py` 做 HTML 交付校验。
- **第4课素材**：`D1` 独立；`D2 → D3` 后再运行 `make_figures_lesson4.py`。整理时发现 D1 / D2 摘要已存在，D3 有明细 CSV，但缺最终 JSON，因此不能认定当前素材已完整可制图。

## 数据口径与复现注意事项

1. **试跑与正式工程分开**：`data/egfr_mvp/` 是早期试跑；`../aidd-mvp/data/raw/` 和 `processed/` 是工程化数据。实操篇、工程篇的最终结论应对照正式工程配置、报告与断言；不能把早期 `stepN` 当成正式 `scripts/NN_*.py` 的一一映射。正式运行方式见 [工程 README](../aidd-mvp/README.md)。
2. **重绘不等于重新实验**：行业图、概念图及部分实验图使用代码内固定数据；另一些图读取 JSON / CSV。变更文章数字时，要同时更新实验来源、绘图输入、图注和正文。
3. **输出路径和覆盖**：多处写死 `/Users/dujiepeng/WorkBuddy/Work`、`~/WorkBuddy/Work` 或依赖相邻工程；运行会覆盖现有 JSON、CSV、图片和 HTML。迁移机器时先检查 `BASE / ROOT / WORK / OUT / FIG`。当前整理只增加文档与注解，未重跑会覆盖文章素材的实验。
4. **环境按任务准备**：常见实验依赖 RDKit、NumPy、pandas、scikit-learn；GCN / LoRA 依赖 PyTorch；结构与口袋实验依赖 SciPy；部分数据通过 TDC / requests 获取；对接依赖 Vina / Meeko 等工具；绘图多使用 matplotlib 与 macOS 中文字体。`vendor/`、`diffdock_demo/venv/` 是依赖目录，不属于自有文章脚本。
5. **结构版本记录待统一**：E1 读取文件名为 `egfr_af_v6.pdb`，其输出元信息仍含 “AlphaFold DB v4 model” 字样。本文档仅标记这一差异，未改动既有实验结果；后续更新文章数字时应先确认实际输入文件与版本。

## 当前文章引用与保留素材的边界

- 行业调研正文当前引用 fig1–fig5、fig7、fig9，以及 `fig0_overview.jpg`。本目录未找到 `fig0_overview.jpg` 的明确生成脚本；fig6 / fig8 虽有生成器，当前正文未引用。
- 实操篇还引用 `trap_overview_card.jpg`；本目录未找到该图片的明确生成脚本。
- 工程篇正文提及本目录的若干脚本，用于说明“草稿间”；实际 11 步流水线与自检代码归属相邻 `aidd-mvp`，本目录的 `make_figures_eng.py` 只负责配图。
- 工具篇正文引用 [`../aidd_ecfp_demo.py`](../aidd_ecfp_demo.py)，该脚本在本项目之外。入门总览、第一课、DeepChem 等文章也在相邻目录，本索引未将它们强行分配给本目录实验。

## 维护方式

新增脚本时，在文件顶部保留“文章关联 / 脚本职责 / 使用注解”，并在此表登记输入输出。新增配图版本时说明覆盖哪些图片；更换数据源或生成方法时，避免复用无法区分方法的结果文件名。文章未找到时明确标记为素材，不推断已发布标题或公众号链接。
