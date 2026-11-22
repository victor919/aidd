from pathlib import Path
import os,json,hashlib
ROOT=Path(__file__).resolve().parents[1]
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/mpl'))
os.environ.setdefault('XDG_CACHE_HOME',str(ROOT/'work/cache'))
import numpy as np,pandas as pd
from sklearn.model_selection import GroupShuffleSplit
from sklearn.metrics import mean_absolute_error,r2_score
from scipy.stats import spearmanr
CONFIG=json.loads((ROOT/'config.json').read_text())
def dump(p,obj):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(obj,ensure_ascii=False,indent=2))
def load_data():
 df=pd.read_csv(ROOT/'data/dataset_with_scaffold.csv',keep_default_na=False);X=np.load(ROOT/'data/X_ecfp.npy');y=df.pIC50_median.to_numpy()
 sp=json.loads((ROOT/'data/split_indices.json').read_text())['scaffold'];outer=np.array(sp['train']);test=np.array(sp['test'])
 a,b=next(GroupShuffleSplit(n_splits=1,test_size=CONFIG['validation_fraction_groups'],random_state=CONFIG['validation_seed']).split(outer,groups=df.iloc[outer].scaffold))
 fit,val=outer[a],outer[b]
 assert set(fit).isdisjoint(val) and set(outer).isdisjoint(test)
 assert len(set(fit)|set(val)|set(test))==len(df)
 for u,v in [(fit,val),(fit,test),(val,test)]:assert set(df.iloc[u].scaffold).isdisjoint(df.iloc[v].scaffold)
 return df,X,y,fit,val,test

def metrics(y,p):
 return {'MAE':float(mean_absolute_error(y,p)),'RMSE':float(np.sqrt(np.mean((y-p)**2))),'R2':float(r2_score(y,p)),'spearman':float(spearmanr(y,p).statistic)}
def make_model():
 from chemprop.nn import BondMessagePassing,MeanAggregation,RegressionFFN
 from chemprop.models import MPNN
 g=CONFIG['gnn'];h=g['hidden_dim']
 return MPNN(BondMessagePassing(d_h=h,depth=g['depth'],dropout=g['dropout']),MeanAggregation(),RegressionFFN(n_tasks=1,input_dim=h,hidden_dim=h,n_layers=g['ffn_layers'],dropout=g['dropout']))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
