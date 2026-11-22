"""从第12篇本地完整配套目录导入输入、权重与逐次结果；不联网，不训练。"""
from pathlib import Path
import argparse,shutil
p=argparse.ArgumentParser();p.add_argument("bundle",type=Path,help="解压后的AI制药第12篇_图神经网络目录");a=p.parse_args()
root=Path(__file__).resolve().parent
assert (a.bundle/"data/dataset_with_scaffold.csv").is_file(),"请提供解压后的配套目录"
assert a.bundle.resolve()!=root,"源目录与目标目录不能相同"
for folder in ["data","models","results","attachment_evidence"]:
    source=a.bundle/folder
    if source.is_dir():shutil.copytree(source,root/folder,dirs_exist_ok=True)
print("快照已导入；运行python scripts/verify_article.py前需先生成配图与文章。")
