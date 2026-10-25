# 零基础入门 AI 制药 · 系列仓库

一个持续更新的公众号系列（「云溪杂谈」）及其配套代码的仓库。

**核心原则**：文章里每一个数字都来自本机实跑，不抄、不编。所以脚本、中间结论（json）与文章三者能互相对上。

**收录范围**：本仓库只收**脚本、数据与工程**。公众号文章 HTML、各篇配图与封面属于发布渠道资产，不入库（脚本可重跑生成，正文只发布于公众号）——见 `.gitignore` 第 9 节。

---

## 目录结构

```
.
├── AI制药系列_脚本与文章对照清单.md   ← 先读这个：哪篇用了哪个脚本、哪张图是谁画的
│
├── （不入库）AI制药*.html              公众号文章成稿 —— 只在公众号发布
├── （不入库）AI制药*_配图/              各篇配图 —— 由 make_figures*.py 重跑生成
├── （不入库）*_封面_*.jpg              公众号 / 小红书封面
│
├── aidd-mvp/                          ← 工程层：11 步流水线，可复跑、可自检
│   ├── config.yaml                    单一参数源，禁止硬编码
│   ├── scripts/00_preflight.py        环境预检（必需 FAIL / 可选 WARN）
│   ├── scripts/01..11_*.py            流水线主体
│   ├── scripts/99_selftest.py         把踩过的坑固化成断言
│   ├── scripts/make_report.py         汇总 json → reports/summary.md
│   ├── tests/                         unittest，真跑
│   └── run_all.sh                     一键跑通（实测 53m03s）
│
├── aidd_practice/                     ← 探索层：试跑、验证素材、画配图
│   ├── af_demo/                       ⑨ 细节篇①（e1–e7 正式脚本 + diag* 诊断脚本）
│   ├── diffdock_demo/                ⑩ 细节篇② 素材（d1–d3）
│   └── make_figures*.py               各篇绘图脚本
│
└── .gitignore                         排除 vendor / 训练模型 / 中间数据 / 文章配图（约 1.7GB）
```

## 三层资产关系

| 层 | 目录 | 职责 | 是否可重建 |
|---|---|---|---|
| 探索层 | `aidd_practice/` | 试跑、验证素材、出图 | 是（一次性，允许乱） |
| 工程层 | `aidd-mvp/` | 流水线、配置、自检 | 是（`run_all.sh`） |
| 文章层 | `*.html` | 只消费上面两层，本身不计算 | 否（人写的，**不入库**） |

**改任何一个正文数字的正确顺序**：改脚本 → 重跑 → 更新 json → 再改 HTML。不要直接改 HTML。

> 文章与配图留在本地工作区即可；仓库里保留的是**能重新产出它们的脚本**。需要文章原文请看公众号，需要配图请跑 `make_figures*.py`。

## 快速开始

```bash
# 工程层：一键跑通全部 11 步
cd aidd-mvp && bash run_all.sh
python scripts/99_selftest.py --full     # 21 项断言
python -m unittest discover -s tests     # 注意：环境无 pytest

# 探索层：⑨ 细节篇①（组织最完整，推荐先看这套）
cd aidd_practice/af_demo
python e6_end2end.py        # 端到端：下载 → 裁剪 → 找口袋 → 出图
python make_figures_af.py   # 正文图 1–5
python verify_article.py    # 交付前校验
```

## 环境

- Python 3.11（conda env `aidd`）
- RDKit 2026.03.5 / DeepChem 2.8.0 / scikit-learn 1.9.0 / numpy 2.4.6
- 可选：AutoDock Vina（缺了只 WARN，`00_preflight.py` 会跳过对接近一步）

## 两份索引，别看混了

| 文件 | 范围 | 定位 |
|---|---|---|
| `AI制药系列_脚本与文章对照清单.md`（根目录） | 全系列 | **溯源表**：哪篇用了哪个脚本、哪张图是谁画的、哪些资产断链、哪版才生效 |
| `aidd_practice/README.md` | 仅 `aidd_practice/` | 按文章查脚本的速查表，不含断链与版本判定 |

两者互补，不重复：前者用于「改数字前查来源」，后者用于「找脚本」。

## 已知断链

- **工具篇 5 张配图无生成脚本**（2026-09-27 产出，脚本未留存），如需改图只能重画。
- ③ 第 2 课、④ 政策篇已发布，但本地未留存 HTML。

详见 `AI制药系列_脚本与文章对照清单.md` 第 3 节。
