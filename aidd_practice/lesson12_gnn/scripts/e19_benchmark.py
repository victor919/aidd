"""RF vs Chemprop D-MPNN，4样本量×3种子×2模型，共24次拟合。
python scripts/e19_benchmark.py        # 完整实验，可从结果断点续跑
python scripts/e19_benchmark.py --sizes 500 --seeds 42  # 小规模通路检查

手写 PyTorch 训练循环调用真实 Chemprop 模型，不使用伪API。
训练批次目标在loss内标准化，验证/测试输出手动反变换一次；不变更数据集标签。
固定预设配置，无超参搜索；GNN仅用验证集选择epoch。RF计时仅fit；GNN计时含fit与逐轮模型选择验证。不含共同预处理和最终测试推断。
"""
from common import *
import argparse,time,copy,random,torch
from torch.utils.data import DataLoader,Subset
from chemprop.data import MoleculeDataset,MoleculeDatapoint
from chemprop.data.collate import collate_batch
from chemprop.featurizers import SimpleMoleculeMolGraphFeaturizer
from rdkit import Chem,RDLogger
from sklearn.ensemble import RandomForestRegressor
RDLogger.DisableLog('rdApp.*')

def loader(ds,idx,shuffle=False,seed=42):
 return DataLoader(Subset(ds,[int(i) for i in idx]),batch_size=CONFIG['gnn']['batch_size'],shuffle=shuffle,num_workers=0,collate_fn=collate_batch,generator=torch.Generator().manual_seed(seed))
@torch.no_grad()
def predict(model,dl,mean,scale):
 model.eval();out=[]
 for b in dl:out.append(model(b.bmg,b.V_d,b.X_d).detach().cpu().numpy().reshape(-1)*scale+mean)
 return np.concatenate(out)

def fit_gnn(ds,tr,val,y,seed):
 random.seed(seed);np.random.seed(seed);torch.manual_seed(seed);torch.set_num_threads(CONFIG['gnn']['threads']);torch.use_deterministic_algorithms(True)
 g=CONFIG['gnn'];model=make_model();opt=torch.optim.Adam(model.parameters(),lr=g['init_lr']);mean=float(y[tr].mean());scale=float(y[tr].std(ddof=0));train=loader(ds,tr,True,seed);validation=loader(ds,val)
 best=float('inf');state=None;wait=0;history=[];t=time.perf_counter()
 for epoch in range(1,g['max_epochs']+1):
  lr=g['init_lr']+(g['max_lr']-g['init_lr'])*epoch/2 if epoch<=2 else g['max_lr']*(g['final_lr']/g['max_lr'])**((epoch-2)/(g['max_epochs']-2))
  for group in opt.param_groups:group['lr']=lr
  model.train();losses=[]
  for b in train:
   opt.zero_grad(set_to_none=True);p=model(b.bmg,b.V_d,b.X_d);target=(b.Y-mean)/scale;loss=((p-target)**2).mean();assert torch.isfinite(loss);loss.backward();torch.nn.utils.clip_grad_norm_(model.parameters(),5.);opt.step();losses.append(float(loss.detach()))
  pv=predict(model,validation,mean,scale);mae=float(np.abs(pv-y[val]).mean());history.append({'epoch':epoch,'train_normalized_MSE':float(np.mean(losses)),'validation_MAE':mae,'lr':lr})
  if mae<best-1e-8:best=mae;state=copy.deepcopy(model.state_dict());best_epoch=epoch;wait=0
  else:wait+=1
  if epoch%10==0:print(f'  GNN n={len(tr)} seed={seed} epoch={epoch} val_MAE={mae:.4f}',flush=True)
  if epoch>=g['min_epochs'] and wait>=g['patience']:break
 elapsed=time.perf_counter()-t;model.load_state_dict(state);model.eval()
 return model,mean,scale,history,{'fit_seconds':elapsed,'best_epoch':best_epoch,'epochs_executed':epoch,'validation_MAE':best,'parameters':sum(p.numel() for p in model.parameters() if p.requires_grad)}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--sizes',nargs='+');ap.add_argument('--seeds',nargs='+',type=int);args=ap.parse_args()
 sizes=args.sizes or [str(x) for x in CONFIG['sizes']];seeds=args.seeds or CONFIG['seeds'];df,X,y,fit,val,test=load_data()
 torch.set_num_threads(CONFIG['gnn']['threads'])
 dps=[MoleculeDatapoint(Chem.MolFromSmiles(s),y=np.array([float(t)],dtype=np.float32)) for s,t in zip(df.smiles,y)];ds=MoleculeDataset(dps,featurizer=SimpleMoleculeMolGraphFeaturizer());ds.cache=True
 test_loader=loader(ds,test);validation_loader=loader(ds,val)
 for seed in seeds:
  order=np.random.default_rng(seed).permutation(fit)
  for size in sizes:
   n=len(fit) if size=='full' else int(size);assert n<=len(fit);tr=order[:n];tag=f'n{n}_seed{seed}';record=ROOT/f'results/{tag}.json'
   if record.exists():
    old=json.loads(record.read_text());assert old['config_sha256']==sha(ROOT/'config.json');print('skip',tag,flush=True);continue
   print('START',tag,flush=True);dump(ROOT/f'results/indices_{tag}.json',tr.tolist());rows=[]
   t=time.perf_counter();rf=RandomForestRegressor(random_state=seed,**CONFIG['rf']);rf.fit(X[tr],y[tr]);seconds=time.perf_counter()-t
   preds={'RF':rf.predict(X[test])};rows.append({'model':'RF','n_train':n,'seed':seed,'fit_seconds':seconds,'validation_MAE':float(np.abs(rf.predict(X[val])-y[val]).mean()),**metrics(y[test],preds['RF'])})
   model,mean,scale,hist,stats=fit_gnn(ds,tr,val,y,seed);preds['GNN']=predict(model,test_loader,mean,scale);rows.append({'model':'GNN','n_train':n,'seed':seed,**stats,**metrics(y[test],preds['GNN'])})
   pd.DataFrame({'index':test,'smiles':df.iloc[test].smiles.to_numpy(),'scaffold':df.iloc[test].scaffold.to_numpy(),'true':y[test],**preds}).to_csv(ROOT/f'results/predictions_{tag}.csv',index=False)
   dump(ROOT/f'results/history_{tag}.json',hist)
   if n==len(fit) and seed==42:
    torch.save({'state_dict':model.state_dict(),'mean':mean,'scale':scale,'config':CONFIG},ROOT/'models/gnn_full_seed42.pt');dump(ROOT/'results/rf_full_seed42_importance.json',rf.feature_importances_.tolist())
   dump(record,{'config_sha256':sha(ROOT/'config.json'),'train_indices_sha256':sha(ROOT/f'results/indices_{tag}.json'),'rows':rows});print('DONE',tag,json.dumps(rows),flush=True)
 # 只聚合已完成的结果，输出状态不伪装成完整实验。
 rows=[]
 for p in (ROOT/'results').glob('n*_seed*.json'):rows+=json.loads(p.read_text())['rows']
 pd.DataFrame(rows).to_csv(ROOT/'results/benchmark.csv',index=False);print('TOTAL completed fits',len(rows),flush=True)
if __name__=='__main__':main()
