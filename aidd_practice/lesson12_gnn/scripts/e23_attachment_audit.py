"""附件快照复核与补充诊断，不重新训练：
1 NPZ复算vs JSON；2主实验/附件协议分别作骨架簇配对bootstrap；
3正确训练集缩放统计量重载300维模型，修正前200分子缩放和非规范SMILES成员判断。
"""
from common import *
import copy,torch
from rdkit import Chem
from chemprop.nn import BondMessagePassing,MeanAggregation,RegressionFFN
from chemprop.models import MPNN
from chemprop.data import MoleculeDatapoint,MoleculeDataset,build_dataloader
from chemprop.data.collate import BatchMolGraph
from chemprop.featurizers import SimpleMoleculeMolGraphFeaturizer
EVIDENCE=ROOT/'attachment_evidence'
def cluster_interval(delta,scaffolds,n_boot=2000):
 rng=np.random.default_rng(42);groups=np.unique(scaffolds);sums=np.array([delta[scaffolds==g].sum() for g in groups]);ns=np.array([(scaffolds==g).sum() for g in groups]);draw=rng.integers(0,len(groups),size=(n_boot,len(groups)));boot=sums[draw].sum(1)/ns[draw].sum(1)
 return {'mean_delta_GNN_minus_RF':float(np.mean(delta)),'cluster_bootstrap95':np.quantile(boot,[.025,.975]).tolist(),'n_test':int(len(delta)),'n_scaffold_clusters':len(groups),'n_bootstrap':n_boot,'scope':'条件于当前测试集和已训练模型；先按种子平均成对绝对误差，再按骨架整组重抽样。不含重训/换划分波动，不是普遍优劣证明。'}
def main():
 torch.set_num_threads(1);df,X,y,fit,val,test=load_data();splits=json.loads((ROOT/'data/split_indices.json').read_text());outer=np.array(splits['scaffold']['train']);canonical=[Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in df.smiles]
 records=json.loads((EVIDENCE/'results/g1_multiseed.json').read_text())['runs'];checked=[];mismatches=[];sources={}
 for p in sorted((EVIDENCE/'results').glob('pred_*.npz')):
  _,split,model,seedtext=p.stem.split('_');seed=int(seedtext[4:]);z=np.load(p);idx=z['idx'].astype(int);assert np.array_equal(idx,np.array(splits[split]['test']));assert np.allclose(z['true'],y[idx]);m=metrics(y[idx],z['pred']);row={'split':split,'model':model,'seed':seed,**m};checked.append(row);sources[split,model,seed]=z['pred'].astype(float)
  expected=next(r for r in records if r['split']==split and r['seed']==seed and r['model']==('D-MPNN' if model=='gnn' else 'ECFP4+RF'))
  for key,target in [('MAE','MAE'),('RMSE','RMSE'),('R2','R2'),('spearman','rho')]:
   if abs(m[key]-expected[target])>5.1e-5:mismatches.append({'file':p.name,'metric':key,'json_value':expected[target],'npz_recomputed':m[key]})
 # 分别核对原始g1快照；g4是另外一批重复记录，不与g1混合。
 scaf=[r for r in checked if r['split']=='scaffold'];assert len(scaf)==10
 assert not any(x for x in mismatches if 'scaffold' in x['file'])
 pd.DataFrame(checked).to_csv(ROOT/'results/attachment_snapshot_metrics.csv',index=False)
 main_deltas=[]
 for seed in CONFIG['seeds']:
  p=pd.read_csv(ROOT/f'results/predictions_n{len(fit)}_seed{seed}.csv');main_deltas.append(np.abs(p.GNN.to_numpy()-p.true.to_numpy())-np.abs(p.RF.to_numpy()-p.true.to_numpy()))
 extra_seeds=[42,7,13,21,33];extra_deltas=[np.abs(sources['scaffold','gnn',seed]-y[test])-np.abs(sources['scaffold','rf',seed]-y[test]) for seed in extra_seeds]
 intervals={'main_128dim_3seeds':cluster_interval(np.mean(main_deltas,axis=0),df.iloc[test].scaffold.to_numpy()),'attachment_300dim_5seeds':cluster_interval(np.mean(extra_deltas,axis=0),df.iloc[test].scaffold.to_numpy())}
 # Reconstruct EXACT original architecture, normalize from all8698 train labels (training source confirmed).
 model=MPNN(BondMessagePassing(d_h=300,depth=3,dropout=.05),MeanAggregation(),RegressionFFN(n_tasks=1,input_dim=300,hidden_dim=300,n_layers=2,dropout=.05))
 model.load_state_dict(torch.load(EVIDENCE/'results/model_scaffold_gnn_seed42.pt',map_location='cpu',weights_only=True));model.eval();mean=float(y[outer].mean());scale=float(y[outer].std(ddof=0));wrong_mean=float(y[outer[:200]].mean());wrong_scale=float(y[outer[:200]].std(ddof=0));feat=SimpleMoleculeMolGraphFeaturizer()
 ds=MoleculeDataset([MoleculeDatapoint(Chem.MolFromSmiles(df.iloc[i].smiles),y=np.array([y[i]])) for i in test],featurizer=feat);dl=build_dataloader(ds,batch_size=128,shuffle=False,num_workers=0);norm=[]
 with torch.no_grad():
  for b in dl:norm.extend(model(b.bmg).numpy().reshape(-1).tolist())
 torch.save({'state_dict':model.state_dict(),'mean':mean,'scale':scale,'train_indices':outer.tolist(),'config':{'depth':3,'hidden_dim':300,'ffn_layers':2,'dropout':.05,'aggregation':'mean','epochs':50,'seed':42},'source_checkpoint_sha256':sha(EVIDENCE/'results/model_scaffold_gnn_seed42.pt')},ROOT/'models/attachment_300dim_seed42_with_scaler.pt')
 corrected=np.array(norm)*scale+mean;maxdiff=float(np.max(np.abs(corrected-sources['scaffold','gnn',42])));assert maxdiff<5e-4,maxdiff
 pd.DataFrame({'index':test,'true':y[test],'prediction_reloaded_correct_scale':corrected,'prediction_source_npz':sources['scaffold','gnn',42]}).to_csv(ROOT/'results/attachment_model_reloaded.csv',index=False)
 raw=json.loads((EVIDENCE/'results/g5_attribution.json').read_text());mols=[]
 for original in raw['molecules']:
  mol=Chem.MolFromSmiles(original['smiles']);can=Chem.MolToSmiles(mol);ids=[i for i,s in enumerate(canonical) if s==can];bg=BatchMolGraph([feat(mol)])
  with torch.no_grad():
   base_norm=float(model(bg).item());deltas=[]
   for i in range(mol.GetNumAtoms()):
    masked=copy.deepcopy(bg);masked.V[i]=0;pn=float(model(masked).item());deltas.append((pn-base_norm)*scale)
  wrong=base_norm*wrong_scale+wrong_mean;fixed=base_norm*scale+mean
  assert abs(wrong-original['pred'])<5e-4,(original['tag'],wrong,original['pred'])
  assert np.max(np.abs(np.array(original['occlusion'])-(np.array(deltas)/scale*wrong_scale)))<5e-4
  mols.append({'tag':original['tag'],'canonical_smiles':can,'raw_string_found':bool((df.smiles==original['smiles']).any()),'canonical_dataset_indices':ids,'in_outer_training':any(i in set(outer) for i in ids),'original_reported_prediction':original['pred'],'wrong_scale_prediction_reproduced':wrong,'correct_scale_prediction':fixed,'correct_scale_masked_minus_full_deltas':deltas,'scale_correction_factor':scale/wrong_scale})
 summaries={}
 for model_name in ['rf','gnn']:
  values=np.array([r['MAE'] for r in scaf if r['model']==model_name]);summaries[model_name]={'mean':float(values.mean()),'std_ddof1':float(values.std(ddof=1)),'min':float(values.min()),'max':float(values.max()),'n_seeds':len(values)}
 out={'source_snapshot_sha256':{str(p.relative_to(EVIDENCE)):sha(p) for p in (EVIDENCE/'results').glob('*')},'scaffold_summary':summaries,'random_snapshot_discrepancies':mismatches,'scaffold_npz_records_verified':10,'random_npz_records_examined':6,'separate_g4_repeat_records':json.loads((EVIDENCE/'results/g4_random.json').read_text())['runs'],'model_reload_max_abs_difference':maxdiff,'paired_cluster_bootstrap':intervals,'scaling_correction':{'original_used_first_n':200,'correct_train_n':len(outer),'wrong_mean':wrong_mean,'wrong_scale':wrong_scale,'correct_mean':mean,'correct_scale':scale,'molecules':mols},'source_learning_curve_status':'Only JSON metrics and code supplied; no per-molecule predictions for g2. Do not mix single-seed attachment curve with main three-seed learning curve.'}
 dump(ROOT/'results/e23_attachment_audit.json',out)
 print(json.dumps({'scaffold':summaries,'mismatches':mismatches,'paired_bootstrap':intervals,'reload_difference':maxdiff,'scaling':out['scaling_correction']},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
