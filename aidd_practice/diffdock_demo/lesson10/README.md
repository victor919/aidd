# ⑩ 细节篇②：分子对接、评分与三维空间语义

本文代码与实验资料：区分“姿势摆得对不对”与“活性排得准不准”，复核 Vina 保存结果，并用 NumPy 理解等变性。

## 快速运行

从仓库根目录执行（Python 3.11，依赖 NumPy）：

```bash
python -m pip install -r aidd_practice/diffdock_demo/lesson10/requirements.txt
python aidd_practice/diffdock_demo/lesson10/scripts/recheck_d2.py
python aidd_practice/diffdock_demo/lesson10/scripts/geometry_checks.py
```

脚本按自身位置定位数据，运行后覆盖本目录 `audit/` 中的复核结果。这里重算的是已有 CSV 的统计量，没有重新运行 Vina，也没有在本机复现 DiffDock。

## 文件与来源

| 文件 | 内容与用途 |
|---|---|
| `scripts/recheck_d2.py` | 从保存的 92 行结果重新计算相关性、按顺序截取的预算命中数及分箱统计 |
| `scripts/geometry_checks.py` | 12 个随机点、24 次旋转：验证等变关系，并用有向体积说明镜像差异 |
| `audit/source/d2_vina_vs_activity.csv` | 原始 92 分子结果，保留 SMILES、实验 pIC50、Vina score 及历史派生列 |
| `audit/source/d2_vina_real_activity.json` | 历史汇总；包含下述已修正问题，不能当最新结论使用 |
| `audit/source/d2_original_for_audit.py` | 历史 d2 脚本存档，仅供溯源；不是推荐运行入口 |
| `audit/source/09_docking.json` | MVP 早期重对接与预测值相关性记录 |
| `audit/source/e8_equivariance.json` | 历史 28 点演示记录，区别于本次 12 点 NumPy 演示 |
| `audit/recalculated.json` | 当前 d2 复核汇总 |
| `audit/selected_molecules.csv` | 按 Vina score 排好的完整结果及实验 Top-30 标记 |
| `audit/geometry_checks.json` | 当前几何演示的数值结果 |
| `audit/sampling_check.json` | 上游全量数据哈希、分箱人数与抽样计数；不包含全量训练数据 |

CSV SHA-256：`52b6978bf4c5c721dc1550350958deddffa5fd6de7ac5b4b44c51c40eb9ba0fe`。

## 当前文章采用的结果

- 92 个分子：Vina score 与实验 pIC50 的 Pearson r = −0.5427068336，Spearman ρ = −0.4983816775。Vina 越低越有利，pIC50 越高越有利，符号应结合这一方向解释。
- 以这 92 个分子中的实验 Top-30 为目标，按 Vina 从低到高挑选 10 / 20 / 30 个，分别命中 **3 / 6 / 14** 个。随机选取的期望命中数分别约为 3.26 / 6.52 / 9.78；这是本样本描述，未给出显著性或外部泛化验证。
- 92 来自等宽活性分箱：第一箱只有 8 个，其余 7 箱各抽 12 个。不能写成“96 个中有 4 个对接失败”。
- 几何演示的旋转等变误差约 1.11e−15；镜像保留距离矩阵，却使有向体积从 +0.12 变成 −0.12。浮点误差的末位可能随环境变化。

## 必须纠正的历史解读

1. 原 d2 先将候选转为集合再截取，丢失排序，因此旧的预算 10/20 命中数 7/9 无效。应先截取有序序列，再求集合交集。
2. 历史 `_cons = 0.7 * pIC50_median - 0.3 * vina_score` 直接用了实验真值，存在标签泄漏。按正确顺序计算得到 10/20/28，但它**不是 QSAR 与对接联合预测的效果**，也不是严格性能上界。脚本仅为审计保留此计算，并明确标记 `NOT_predictive`。
3. MVP 的 |r| = 0.059 比较的是预测 pIC50 与 Vina score，不是 Vina 与实验活性。30 个样本又来自预测 Top-30，不能与本次分层样本直接比较，也不能据弱相关宣称二者独立或互补。
4. 重对接 RMSD 1.411 Å 来自旧脚本按原子顺序比较的 PDBQT 坐标，包含保留氢，未做图匹配或对称性校正。它不等价于标准的对称性校正重原子 RMSD，更不能证明活性排序准确。
5. 距离加权相对向量层对反射也等变；镜像后向量变化不能单独证明模型能识别手性。新演示用固定顶点顺序的有向体积说明这一点。

## 解释边界

相关性与命中数仅针对该 92 分子样本。质子化、互变异构、受体构象及活性测量条件等因素尚未逐一控制，不能直接归因。当前脚本不重新计算文章沿用的历史 p 值，不训练 QSAR，也不验证 DiffDock 的性能。

若要重新生成对接数据，请查上一级 `d2_vina_real_activity.py` 及工程层 `aidd-mvp/`，先阅读以上历史问题并修正统计逻辑；其 RDKit、Vina、受体与上游数据依赖不属于此轻量复核包。公众号 HTML、封面和配图遵循根目录约定，不随本目录入库。
