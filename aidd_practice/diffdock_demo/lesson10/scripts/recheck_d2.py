"""Re-evaluate EXISTING CSV; do not redock or fit a model.
Run from bundle root: python scripts/recheck_d2.py
"""
from pathlib import Path
import csv,json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'audit/source/d2_vina_vs_activity.csv'
rows=list(csv.DictReader(p.open()))
y=np.array([float(r['pIC50_median']) for r in rows]);v=np.array([float(r['vina_score']) for r in rows])
if not np.isfinite(y).all() or not np.isfinite(v).all():raise ValueError('Non-finite input')
if len({r['smiles'] for r in rows})!=len(rows):raise ValueError('Duplicate SMILES; define deduplication first')
def rank(x):
 o=np.argsort(x,kind='stable');r=np.empty(len(x),float);i=0
 while i<len(x):
  j=i+1
  while j<len(x) and x[o[j]]==x[o[i]]:j+=1
  r[o[i:j]]=(i+j-1)/2;i=j
 return r
n=len(rows);true_order=sorted(range(n),key=lambda i:(-y[i],rows[i]['smiles']))
target=set(true_order[:30]);vo=sorted(range(n),key=lambda i:(v[i],rows[i]['smiles']))
# Demonstration ONLY: this score directly consumes evaluation labels.
leak=.7*y-.3*v;lo=sorted(range(n),key=lambda i:(-leak[i],rows[i]['smiles']))
res={'n':n,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'pearson_vina_vs_experimental':float(np.corrcoef(v,y)[0,1]),'spearman_vina_vs_experimental':float(np.corrcoef(rank(v),rank(y))[0,1]),'vina_hits':{str(k):len(target&set(vo[:k])) for k in (10,20,30)},'label_informed_mix_hits_NOT_predictive':{str(k):len(target&set(lo[:k])) for k in (10,20,30)},'random_expected_hits':{str(k):k*30/n for k in (10,20,30)},'top30_target_boundary_tie':bool(y[true_order[29]]==y[true_order[30]]),'vina_boundary_ties':{str(k):bool(v[vo[k-1]]==v[vo[k]]) for k in (10,20,30)},'experimental_range':[float(y.min()),float(y.max())],'legacy_cons_matches_true_label_formula':bool(np.allclose([float(r['_cons']) for r in rows],leak)), 'bins':[],'note':'Recalculated from saved scores; no Vina rerun. p-values in article are inherited, not recalculated.'}
for b in sorted({r['_vbin'] for r in rows},key=int):
 ix=[i for i,r in enumerate(rows) if r['_vbin']==b]
 res['bins'].append({'bin':int(b),'n':len(ix),'vina_min':float(v[ix].min()),'vina_max':float(v[ix].max()),'mean_experimental_pIC50':float(y[ix].mean()),'median_experimental_pIC50':float(np.median(y[ix]))})
(ROOT/'audit/recalculated.json').write_text(json.dumps(res,indent=2,ensure_ascii=False))
with (ROOT/'audit/selected_molecules.csv').open('w',newline='') as f:
 w=csv.writer(f);w.writerow(['vina_rank','smiles','vina_score','experimental_pIC50','in_experimental_top30'])
 for j,i in enumerate(vo):w.writerow([j+1,rows[i]['smiles'],v[i],y[i],i in target])
print(json.dumps(res,indent=2,ensure_ascii=False))
