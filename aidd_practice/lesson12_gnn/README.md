# 第12篇：图神经网络与分子表征升级

公众号配套实验：在固定EGFR骨架测试集上比较2048位ECFP4＋随机森林与Chemprop D-MPNN。主实验4个拟合规模×3个种子×2种模型，共24次拟合。

## 内容

- `scripts/e18_prepare.py`：输入、骨架隔离、图特征与环境审计。
- `scripts/e19_benchmark.py`：RF/GNN训练、验证集选轮次、学习曲线。
- `scripts/e20_sensitivity.py`：节点特征遮蔽；保存差值为完整预测减遮蔽预测。
- `scripts/e21_figures.py`、`make_course_map.py`：科学配图与课程定位图。
- `scripts/e22_article.py`：从已有结果生成文章与微信发布版；产物不入库。
- `scripts/e23_attachment_audit.py`、`e24_supplement_figures.py`：附件协议、成对重抽样与缩放纠错。
- `scripts/e25_case_readout.py`：主案例读图、RF Top20位映射；差值为遮蔽预测减完整预测。
- `scripts/verify_article.py`：复核24次指标、模型重载与图文引用。
- `results/`：已完成实验的汇总与审计；逐分子预测、权重和输入由快照导入。

## 安装与复跑

使用Python 3.11；`requirements.txt`记录原实验精确版本，并非跨平台可安装保证。中文绘图需要中文字体，可调整绘图脚本的字体列表。

```bash
cd aidd_practice/lesson12_gnn
python -m pip install -r requirements.txt
# 从本地“AI制药第12篇_文章与实验配套.zip”解压后的目录导入：
python import_snapshot.py /path/to/AI制药第12篇_图神经网络
bash run_all.sh
```

本目录遵循系列仓库规则，不提交文章HTML、PNG、原始附件、训练权重与大型输入快照。因此仅克隆仓库不能完成全部复跑，必须先取得本篇完整配套包并导入。`run_all.sh`会跳过配置哈希一致的已有训练；更改配置应使用独立结果目录。仅重新生成图文时，依次运行e21、e24、e25、make_course_map、e22和verify_article，无须重训。

## 结果与边界

协议A：7148拟合、1550验证、1770测试；3种子平均测试MAE为RF 0.6334、GNN 0.7948。当前四个规模未见均值排名反转，不能外推为所有图模型的排名。

协议B：附件既有300维模型的复算，不与A混作同一训练。厄洛替尼9.8224→7.7220是同一权重的目标缩放纠错；A案例约7.2996。该分子在拟合池内，不是外部分子验证。节点置零敏感性不是注意力或化学因果贡献。

`results/verification.json`是本地完整包已完成的检查记录，不表示只克隆本目录即可通过。固定输入未记录ChEMBL发布号，不承诺在线重新抓取可获得同一数据集。
