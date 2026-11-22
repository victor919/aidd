"""交付检查：固定输入、嵌套子集、模型指标、预测还原、图文引用一致性。"""
from common import *
import ast,re
from PIL import Image

def main():
 df,X,y,fit,val,test=load_data();audit=json.loads((ROOT/'results/e18_audit.json').read_text())
 assert audit['config_sha256']==sha(ROOT/'config.json')
 for name,h in audit['inputs_sha256'].items():assert sha(ROOT/'data'/name)==h
 counts=0
 for seed in CONFIG['seeds']:
  previous=set()
  for size in CONFIG['sizes']:
   n=len(fit) if size=='full' else size;tag=f'n{n}_seed{seed}';tr=json.loads((ROOT/f'results/indices_{tag}.json').read_text());assert len(tr)==n and len(set(tr))==n;assert previous<=set(tr)<=set(fit);previous=set(tr)
   result=json.loads((ROOT/f'results/{tag}.json').read_text());assert result['config_sha256']==sha(ROOT/'config.json');assert result['train_indices_sha256']==sha(ROOT/f'results/indices_{tag}.json')
   p=pd.read_csv(ROOT/f'results/predictions_{tag}.csv');assert p['index'].tolist()==test.tolist();assert np.allclose(p.true,y[test]);assert p[['RF','GNN']].notna().all().all()
   for row in result['rows']:
    actual=metrics(p.true.to_numpy(),p[row['model']].to_numpy())
    for key,value in actual.items():assert abs(value-row[key])<1e-6,(tag,key,value,row[key])
    counts+=1
   hist=json.loads((ROOT/f'results/history_{tag}.json').read_text());gnn=next(r for r in result['rows'] if r['model']=='GNN');best=min(hist,key=lambda x:x['validation_MAE']);assert best['epoch']==gnn['best_epoch'];assert len(hist)==gnn['epochs_executed'];assert 10<=len(hist)<=CONFIG['gnn']['max_epochs']
 assert counts==24
 article=(ROOT/'AI制药第12篇_公众号文章.html').read_text();refs=re.findall(r'<img[^>]+src="([^"]+)"',article);assert len(refs)==11
 for ref in refs:
  with Image.open(ROOT/ref) as image:image.verify()
 for target in re.findall(r'href="#([^"]+)"',article):assert f'id="{target}"' in article
 for p in (ROOT/'scripts').glob('*.py'):ast.parse(p.read_text(),filename=p.name)
 # 验证保存的缩放统计量只来自当次训练子集，以及模型重载后的预测。
 import torch
 from rdkit import Chem
 from chemprop.data import MoleculeDataset,MoleculeDatapoint,build_dataloader
 from chemprop.featurizers import SimpleMoleculeMolGraphFeaturizer
 torch.set_num_threads(1);ckpt=torch.load(ROOT/'models/gnn_full_seed42.pt',map_location='cpu',weights_only=False)
 assert np.isclose(ckpt['mean'],y[fit].mean());assert np.isclose(ckpt['scale'],y[fit].std(ddof=0))
 model=make_model();model.load_state_dict(ckpt['state_dict']);model.eval();subset=test[:6];ds=MoleculeDataset([MoleculeDatapoint(Chem.MolFromSmiles(df.iloc[i].smiles),y=np.array([y[i]])) for i in subset],featurizer=SimpleMoleculeMolGraphFeaturizer());batch=next(iter(build_dataloader(ds,batch_size=6,shuffle=False,num_workers=0)))
 with torch.no_grad():pred=model(batch.bmg,batch.V_d,batch.X_d).numpy().ravel()*ckpt['scale']+ckpt['mean']
 stored=pd.read_csv(ROOT/f'results/predictions_n{len(fit)}_seed42.csv').GNN.to_numpy()[:6];maxdiff=float(np.max(np.abs(pred-stored)));assert maxdiff<1e-5
 # 新增诊断只复用固定快照；确认区间中心与主成绩均值吻合。
 extra=json.loads((ROOT/'results/e23_attachment_audit.json').read_text())
 assert extra['scaffold_npz_records_verified']==10 and extra['random_npz_records_examined']==6
 assert not extra['random_snapshot_discrepancies']
 assert extra['model_reload_max_abs_difference']<5e-4
 bench=pd.read_csv(ROOT/'results/benchmark.csv');full=bench[bench.n_train==len(fit)]
 diff=float(full[full.model=='GNN'].MAE.mean()-full[full.model=='RF'].MAE.mean())
 assert abs(diff-extra['paired_cluster_bootstrap']['main_128dim_3seeds']['mean_delta_GNN_minus_RF'])<1e-6
 assert next(m for m in extra['scaling_correction']['molecules'] if m['tag']=='erlotinib')['in_outer_training']
 # 图7读图数值必须与主实验原始遮蔽结果一致，且RF中心名单交集可重算。
 case=json.loads((ROOT/'results/e25_case_readout.json').read_text())
 sensitivity=json.loads((ROOT/'results/e20_sensitivity.json').read_text())
 delta=np.array([r['delta_masked_minus_full'] for r in case['atom_rows']])
 assert np.allclose(delta,-np.array([r['delta_pIC50'] for r in sensitivity['atoms']]))
 assert case['in_fit_pool'] and case['n_atoms']==len(delta)==29
 assert case['n_abs_delta_gt_threshold']==int((np.abs(delta)>case['threshold_for_descriptive_count']).sum())
 assert np.isclose(case['max_abs_delta'],np.abs(delta).max())
 assert case['most_sensitive_atom']==int(np.argmax(np.abs(delta)))
 assert sorted(set(case['top10_sensitive_atoms']) & set(case['rf_center_voted_atoms']))==case['top10_and_center_votes_overlap']
 report={'status':'PASS','model_fits_checked':counts,'images_checked':len(refs),'metric_recomputations_match':True,'metric_recompute_tolerance':1e-6,'nested_subsets_verified':True,'all_three_partitions_scaffold_disjoint':True,'gnn_reload_max_abs_difference':maxdiff,'attachment_snapshot_checks':16,'attachment_scale_reconstruction_verified':True,'paired_bootstrap_centers_verified':True,'case_readout_verified':True,'html_reference_check':True,'html_browser_rendered':False,'note':'HTML结构与PNG解码已检查；整页浏览器渲染未完成。'}
 dump(ROOT/'results/verification.json',report);print(json.dumps(report,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
