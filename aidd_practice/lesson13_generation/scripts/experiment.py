"""第13篇固定协议：匹配动作的随机游走/贪心搜索；无生成候选标签。"""
from pathlib import Path
import json,hashlib,time,random,sys,os,platform
import numpy as np,pandas as pd,joblib
from rdkit import Chem,DataStructs,RDLogger,RDConfig
from rdkit.Chem import BRICS,rdFingerprintGenerator,Descriptors,QED
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupShuffleSplit
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[1];os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/mpl'));os.environ.setdefault('XDG_CACHE_HOME',str(ROOT/'work/cache'));C=json.loads((ROOT/'config.json').read_text());RDLogger.DisableLog('rdApp.*')
sys.path.append(str(Path(RDConfig.RDContribDir)/'SA_Score'));import sascorer
GEN=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)
def dump(name,obj):(ROOT/'results'/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2))
def canon(m):return Chem.MolToSmiles(m,isomericSmiles=True)
def bit(m):return GEN.GetFingerprint(m)
def arr(fp):
 a=np.zeros(2048,dtype=np.uint8);DataStructs.ConvertToNumpyArray(fp,a);return a
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def stitch(a,b,bondtype):
 m=Chem.CombineMols(a,b);d=[x.GetIdx() for x in m.GetAtoms() if x.GetAtomicNum()==0]
 if len(d)!=2:raise ValueError('dummy count')
 n=[m.GetAtomWithIdx(i).GetNeighbors()[0].GetIdx() for i in d];rw=Chem.RWMol(m);rw.AddBond(n[0],n[1],bondtype)
 for i in sorted(d,reverse=True):rw.RemoveAtom(i)
 out=rw.GetMol();Chem.SanitizeMol(out);return out

def main():
 t=time.time();(ROOT/'work').mkdir(exist_ok=True)
 df=pd.read_csv(ROOT/'data/dataset_with_scaffold.csv',keep_default_na=False);X=np.load(ROOT/'data/X_ecfp.npy');y=df.pIC50_median.to_numpy();sp=json.loads((ROOT/'data/split_indices.json').read_text())['scaffold'];outer=np.array(sp['train']);test=np.array(sp['test']);a,b=next(GroupShuffleSplit(n_splits=1,test_size=.2,random_state=42).split(outer,groups=df.iloc[outer].scaffold));fit=outer[a];val=outer[b]
 assert len(fit)==7148 and len(val)==1550 and len(test)==1770
 for u,v in [(fit,val),(fit,test),(val,test)]:assert set(df.iloc[u].scaffold).isdisjoint(df.iloc[v].scaffold)
 dump('indices.json',{'fit':fit.tolist(),'validation':val.tolist(),'test':test.tolist()})
 model_path=ROOT/'work/rf.joblib'
 if model_path.exists():rf=joblib.load(model_path)
 else:
  print('Fitting fixed RF',flush=True);rf=RandomForestRegressor(n_estimators=500,max_features=1.0,min_samples_leaf=1,random_state=42,n_jobs=4);tr=np.random.default_rng(42).permutation(fit);rf.fit(X[tr],y[tr]);joblib.dump(rf,model_path)
 old=pd.read_csv(ROOT/'data/lesson12_predictions.csv');test_pred=rf.predict(X[test]);maxdiff=float(np.abs(test_pred-old.RF.to_numpy()).max());assert maxdiff<1e-6
 print('RF reproduced',maxdiff,flush=True)
 mols=[Chem.MolFromSmiles(s) for s in df.smiles];known={canon(m) for m in mols};fit_smiles={canon(mols[i]) for i in fit};refs=[bit(mols[i]) for i in fit];fit_scaffolds=set(df.iloc[fit].scaffold)
 # Materialize terminal single-dummy BRICS fragments, sorted before random selection.
 pool={}
 for i in fit:
  for s in BRICS.BRICSDecompose(mols[i]):
   f=Chem.MolFromSmiles(s);d=[a for a in f.GetAtoms() if a.GetAtomicNum()==0]
   if len(d)==1 and f.GetNumHeavyAtoms()>=1:pool.setdefault(d[0].GetIsotope(),set()).add(canon(f))
 pool={int(k):sorted(v) for k,v in sorted(pool.items())};dump('fragment_pool.json',pool)
 print('fragment pool',sum(map(len,pool.values())),flush=True)
 # Starts sampled without looking at measured activity or model predictions.
 eligible=[int(i) for i in sorted(fit) if 15<=mols[i].GetNumHeavyAtoms()<=50 and list(BRICS.FindBRICSBonds(mols[i]))]
 starts=random.Random(20261122).sample(eligible,C['start_count']);dump('starts.json',[{'index':i,'smiles':canon(mols[i])} for i in starts])
 cache={};logical=0
 def score_many(smis):
  nonlocal logical
  logical+=len(smis);missing=list(dict.fromkeys(s for s in smis if s not in cache))
  if missing:
   mm=[Chem.MolFromSmiles(s) for s in missing];fp=[bit(m) for m in mm];pred=rf.predict(np.array([arr(f) for f in fp]))
   for s,m,f,p in zip(missing,mm,fp,pred):
    sims=DataStructs.BulkTanimotoSimilarity(f,refs);j=int(np.argmax(sims));scaf=MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    cache[s]={'smiles':s,'RF':float(p),'nn_similarity':float(sims[j]),'nn_index':int(fit[j]),'MW':Descriptors.MolWt(m),'QED':QED.qed(m),'SA':sascorer.calculateScore(m),'heavy_atoms':m.GetNumHeavyAtoms(),'scaffold':scaf,'seen_fit':s in fit_smiles,'seen_snapshot':s in known,'scaffold_seen_fit':scaf in fit_scaffolds}
  return [cache[s] for s in smis]
 proposals=[];attempts=[];states=[]
 for seed in C['seeds']:
  for start_id in starts:
   initial=canon(mols[start_id]);score_many([initial])
   for arm in ['A','B']:
    rng=random.Random(seed*100000+start_id);current=initial
    states.append({'seed':seed,'start_index':start_id,'arm':arm,'round':0,'budget':0,**cache[current]})
    for step in range(1,C['rounds']+1):
     parent=current;pm=Chem.MolFromSmiles(parent);bonds=sorted(BRICS.FindBRICSBonds(pm));batch=[];actions=[];seen=set();ntry=0
     while len(batch)<C['proposals_per_round'] and ntry<C['max_attempts_per_round']:
      ntry+=1;row={'seed':seed,'start_index':start_id,'arm':arm,'round':step,'attempt':ntry,'parent':parent};reason=''
      try:
       if not bonds:raise ValueError('no_brics_bond')
       (i,j),labels=rng.choice(bonds);bt=pm.GetBondBetweenAtoms(i,j).GetBondType();broken=BRICS.BreakBRICSBonds(pm,bonds=[((i,j),labels)]);sides=Chem.GetMolFrags(broken,asMols=True,sanitizeFrags=True)
       if len(sides)!=2:raise ValueError('not_two_sides')
       discard=rng.randrange(2);keep=sides[1-discard];removed=sides[discard];ds=[a for a in removed.GetAtoms() if a.GetAtomicNum()==0]
       if len(ds)!=1:raise ValueError('not_single_dummy')
       label=ds[0].GetIsotope();donor_s=rng.choice(pool[label]);donor=Chem.MolFromSmiles(donor_s);product=stitch(keep,donor,bt);s=canon(product)
       row.update({'kept_fragment':canon(keep),'removed_fragment':canon(removed),'donor_fragment':donor_s,'attachment_label':label,'candidate':s})
       if not C['min_heavy_atoms']<=product.GetNumHeavyAtoms()<=C['max_heavy_atoms']:reason='size_cap'
       elif s==parent:reason='unchanged'
       elif s in seen:reason='duplicate_in_round'
       else:seen.add(s);batch.append(s);actions.append(row.copy());reason='accepted_proposal'
      except Exception as e:reason=str(e) if isinstance(e,ValueError) else 'chemistry_or_assembly_failure'
      row['status']=reason;attempts.append(row)
     # Fixed request budget: exhausted neighborhoods are explicit self-loops, not fabricated candidates.
     while len(batch)<C['proposals_per_round']:
      batch.append(parent);actions.append({'seed':seed,'start_index':start_id,'arm':arm,'round':step,'parent':parent,'candidate':parent,'status':'budget_padding_self_loop'})
     values=score_many(batch)
     options=[parent]+batch
     chosen=rng.randrange(len(options)) if arm=='A' else max(range(len(options)),key=lambda k:cache[options[k]]['RF'])
     current=options[chosen]
     for k,(s,act) in enumerate(zip(batch,actions),1):proposals.append({**act,'slot':k,'selected':k==chosen,**cache[s]})
     states.append({'seed':seed,'start_index':start_id,'arm':arm,'round':step,'budget':step*C['proposals_per_round'],'attempts':ntry,'selected_current':chosen==0,**cache[current]})
    print('trajectory',seed,start_id,arm,round(cache[current]['RF'],3),flush=True)
 pd.DataFrame(states).to_csv(ROOT/'results/states.csv',index=False);pd.DataFrame(proposals).to_csv(ROOT/'results/proposals.csv',index=False);pd.DataFrame(attempts).to_csv(ROOT/'results/attempts.csv',index=False)
 print('GNN posthoc scoring',len(cache),flush=True)
 import torch
 from chemprop.nn import BondMessagePassing,MeanAggregation,RegressionFFN
 from chemprop.models import MPNN
 from chemprop.featurizers import SimpleMoleculeMolGraphFeaturizer
 from chemprop.data.collate import BatchMolGraph
 torch.set_num_threads(1);ck=torch.load(ROOT/'models/gnn_full_seed42.pt',map_location='cpu',weights_only=False);g=ck['config']['gnn'];h=g['hidden_dim'];model=MPNN(BondMessagePassing(d_h=h,depth=g['depth'],dropout=g['dropout']),MeanAggregation(),RegressionFFN(n_tasks=1,input_dim=h,hidden_dim=h,n_layers=g['ffn_layers'],dropout=g['dropout']));model.load_state_dict(ck['state_dict']);model.eval();feat=SimpleMoleculeMolGraphFeaturizer();keys=list(cache)
 with torch.no_grad():
  for k in range(0,len(keys),128):
   ss=keys[k:k+128];bg=BatchMolGraph([feat(Chem.MolFromSmiles(s)) for s in ss]);p=model(bg).numpy().ravel()*ck['scale']+ck['mean']
   for s,v in zip(ss,p):cache[s]['GNN']=float(v);cache[s]['disagreement']=abs(float(v)-cache[s]['RF'])
 pd.DataFrame(cache.values()).to_csv(ROOT/'results/candidates.csv',index=False)
 # Historical test diagnostics stay outside search and do not choose a threshold.
 test_nn=[]
 for i in test:test_nn.append(max(DataStructs.BulkTanimotoSimilarity(bit(mols[i]),refs)))
 old['nn_similarity']=test_nn;old['disagreement']=abs(old.RF-old.GNN);old['RF_abs_error']=abs(old.RF-old.true);old['GNN_abs_error']=abs(old.GNN-old.true);old.to_csv(ROOT/'results/historical_test_diagnostic.csv',index=False)
 meta={'config':C,'config_sha256':digest(ROOT/'config.json'),'input_sha256':{p.name:digest(p) for p in (ROOT/'data').iterdir()},'gnn_sha256':digest(ROOT/'models/gnn_full_seed42.pt'),'n_fit':len(fit),'n_validation':len(val),'n_test':len(test),'rf_reload_prediction_maxdiff':maxdiff,'fragment_pool_size':sum(map(len,pool.values())),'fragment_pool_label_counts':{k:len(v) for k,v in pool.items()},'starts':starts,'trajectories':len(starts)*len(C['seeds'])*2,'logical_prediction_requests_including_starts':logical,'unique_rf_evaluations':len(cache),'search_proposal_slots':len(proposals),'attempts':len(attempts),'seconds':time.time()-t,'environment':{'python':platform.python_version(),'platform':platform.platform()},'limitations':['No experimental labels for generated candidates','Only one oracle, one local action family and six fitting-pool starts','Historical test is diagnostic only, not external prospective validation','No similarity threshold, docking or ADMET model used in selection','SA and QED are descriptors, not synthesis routes or ADMET evidence']};dump('audit.json',meta);print('DONE',meta['seconds'],flush=True)
if __name__=='__main__':main()
