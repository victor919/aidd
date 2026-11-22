"""验证输入隔离、预算、选择规则、评分口径、片段来源及发布链接。"""
from experiment import *
import ast,re
from PIL import Image
from urllib.parse import urlsplit

def main():
 audit=json.loads((ROOT/'results/audit.json').read_text());assert digest(ROOT/'config.json')==audit['config_sha256']
 for name,sha in audit['input_sha256'].items():assert digest(ROOT/'data'/name)==sha
 assert digest(ROOT/'models/gnn_full_seed42.pt')==audit['gnn_sha256']
 idx=json.loads((ROOT/'results/indices.json').read_text());df=pd.read_csv(ROOT/'data/dataset_with_scaffold.csv',keep_default_na=False)
 for u,v in [('fit','validation'),('fit','test'),('validation','test')]:assert set(df.iloc[idx[u]].scaffold).isdisjoint(df.iloc[idx[v]].scaffold)
 starts=json.loads((ROOT/'results/starts.json').read_text());assert all(x['index'] in idx['fit'] for x in starts)
 c=pd.read_csv(ROOT/'results/candidates.csv');p=pd.read_csv(ROOT/'results/proposals.csv');s=pd.read_csv(ROOT/'results/states.csv');attempt=pd.read_csv(ROOT/'results/attempts.csv');lookup=c.set_index('smiles');pool=json.loads((ROOT/'results/fragment_pool.json').read_text());fitmols=[Chem.MolFromSmiles(x) for x in df.iloc[idx['fit']].smiles];fitfp=[bit(m) for m in fitmols]
 assert len(p)==2880 and len(s)==396 and len(c)==audit['unique_rf_evaluations'];assert c.smiles.is_unique
 for (seed,start,arm),g in s.groupby(['seed','start_index','arm']):
  g=g.sort_values('round');assert g['round'].tolist()==list(range(11));assert g.budget.tolist()==list(range(0,81,8))
  assert g.iloc[0].smiles==next(x['smiles'] for x in starts if x['index']==start)
  for i in range(1,11):
   pp=p[(p.seed==seed)&(p.start_index==start)&(p.arm==arm)&(p['round']==i)];assert len(pp)==8
   assert set(pp.parent)=={g.iloc[i-1].smiles};assert g.iloc[i].smiles in [g.iloc[i-1].smiles]+pp.smiles.tolist()
   if arm=='B':assert np.isclose(g.iloc[i].RF,max([g.iloc[i-1].RF]+pp.RF.tolist()))
 # Sample all recorded successful transformations: donors belong to fitting-pool fragment library.
 successful=attempt[attempt.status=='accepted_proposal']
 for _,x in successful.iterrows():
  assert x.donor_fragment in pool[str(int(x.attachment_label))]
  parent=Chem.MolFromSmiles(x.parent);valid_cuts=[]
  for bond,labels in BRICS.FindBRICSBonds(parent):
   cut=BRICS.BreakBRICSBonds(parent,bonds=[(bond,labels)]);sides=Chem.GetMolFrags(cut,asMols=True)
   if len(sides)==2 and x.kept_fragment in [canon(z) for z in sides] and x.removed_fragment in [canon(z) for z in sides]:
    bt=parent.GetBondBetweenAtoms(*bond).GetBondType();valid_cuts.append(bt)
  assert valid_cuts
  rebuilt=stitch(Chem.MolFromSmiles(x.kept_fragment),Chem.MolFromSmiles(x.donor_fragment),valid_cuts[0]);assert canon(rebuilt)==x.candidate
 # Recompute reference similarity for deterministic sample; verify no label column added to generated candidates.
 sample=c.sample(min(40,len(c)),random_state=42)
 for _,row in sample.iterrows():
  sims=DataStructs.BulkTanimotoSimilarity(bit(Chem.MolFromSmiles(row.smiles)),fitfp);assert abs(max(sims)-row.nn_similarity)<1e-12;assert row.nn_index in idx['fit']
 assert not any(x in c.columns for x in ['true','pIC50_median','measured_activity'])
 extra=json.loads((ROOT/'results/attachment_readout.json').read_text())
 diag=pd.read_csv(ROOT/'results/historical_test_diagnostic.csv')
 mask=(diag.disagreement<.5)&(diag.RF_abs_error>2)&(diag.GNN_abs_error>2)
 assert extra['historical_both_wrong_count']==int(mask.sum())
 assert extra['triangle_inequality_max_violation']<1e-10
 assert len(extra['seed_pairs'])==3
 text=(ROOT/'AI制药第13篇_公众号文章.html').read_text();refs=re.findall(r'<img[^>]*src="([^"]+)"',text);assert len(refs)==9
 for ref in refs:
  with Image.open(ROOT/ref) as im:im.verify()
 for target in re.findall(r'href="#([^"]+)"',text):assert f'id="{target}"' in text
 pub=(ROOT.parent/'AI制药第13篇_公众号发布版.html').read_text();links=re.findall(r'<a[^>]*href="([^"]+)"',pub);assert links and all(urlsplit(x).hostname=='mp.weixin.qq.com' for x in links)
 for x in (ROOT/'scripts').glob('*.py'):ast.parse(x.read_text())
 report={'status':'PASS','trajectories':36,'proposal_slots':2880,'states':396,'unique_evaluated_molecules':len(c),'successful_actions_reconstructed':len(successful),'reference_similarity_spotchecks':len(sample),'scaffold_partitions_disjoint':True,'oracle_matches_lesson12':audit['rf_reload_prediction_maxdiff'],'greedy_choices_verified':True,'generated_labels_not_used':True,'images':9,'publication_links_wechat_only':True,'html_browser_rendered':False}
 dump('verification.json',report);print(json.dumps(report,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
