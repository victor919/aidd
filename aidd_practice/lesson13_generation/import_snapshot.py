"""导入本地第13篇完整配套目录；不联网、不训练。"""
from pathlib import Path
import argparse,shutil
p=argparse.ArgumentParser();p.add_argument("bundle",type=Path,help="解压后的AI制药第13篇_分子生成目录");a=p.parse_args()
root=Path(__file__).resolve().parent
assert (a.bundle/"data/dataset_with_scaffold.csv").is_file(),"请提供配套包解压目录"
assert a.bundle.resolve()!=root,"源目录与目标目录不能相同"
for folder in ["data","models","results"]:
    if (a.bundle/folder).is_dir():shutil.copytree(a.bundle/folder,root/folder,dirs_exist_ok=True)
print("快照已导入；可复跑或从analyze开始重生成图文。")
