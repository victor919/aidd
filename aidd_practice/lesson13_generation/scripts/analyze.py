from experiment import *
def main():
 audit=json.loads((ROOT/'results/audit.json').read_text());s=pd.read_csv(ROOT/'results/states.csv');p=pd.read_csv(ROOT/'results/proposals.csv');a=pd.read_csv(ROOT/'results/attempts.csv');c=pd.read_csv(ROOT/'results/candidates.csv');d=pd.read_csv(ROOT/'results/historical_test_diagnostic.csv')
 for seed in C['seeds']:
  for start in audit['starts']:
   for arm in ['A','B']:
    t=s[(s.seed==seed)&(s.start_index==start)&(s.arm==arm)].sort_values('round');assert len(t)==11;assert t.budget.tolist()==list(range(0,81,8))
    if arm=='B':assert np.all(np.diff(t.RF)>=-1e-10)
 assert len(p)==36*80
 init=s[s['round']==0].set_index(['seed','start_index','arm']);end=s[s['round']==10].set_index(['seed','start_index','arm']);changes=end[['RF','nn_similarity','SA','MW']]-init[['RF','nn_similarity','SA','MW']];changes.to_csv(ROOT/'results/paired_changes.csv')
 rows=[]
 for arm in ['A','B']:
  e=end.xs(arm,level='arm');v=changes.xs(arm,level='arm');pp=p[p.arm==arm];aa=a[a.arm==arm];unique=pp.drop_duplicates('smiles');valid=aa[aa.status.isin(['accepted_proposal','duplicate_in_round','unchanged'])];parsevalid=aa[aa.candidate.notna()]
  rows.append({'arm':arm,'trajectories':len(e),'slots':len(pp),'attempts':len(aa),'assembled_sanitized_attempts':len(parsevalid),'sanitized_attempt_rate':len(parsevalid)/len(aa),'admissible_attempts':len(valid),'unique_scored_candidates':len(unique),'novel_to_fit_fraction':float((~unique.seen_fit).mean()),'novel_to_snapshot_fraction':float((~unique.seen_snapshot).mean()),'scored_scaffolds':unique.scaffold.nunique(),'end_RF_mean':float(e.RF.mean()),'end_RF_std':float(e.RF.std()),'end_nn_mean':float(e.nn_similarity.mean()),'delta_RF_mean':float(v.RF.mean()),'delta_nn_mean':float(v.nn_similarity.mean()),'delta_SA_mean':float(v.SA.mean()),'delta_MW_mean':float(v.MW.mean()),'end_known_fit_states':int(e.seen_fit.sum()),'end_known_fit_unique':int(e.drop_duplicates('smiles').seen_fit.sum()),'end_unique':e.smiles.nunique(),'end_scaffolds':e.scaffold.nunique(),'selected_current_fraction':float(s[(s.arm==arm)&(s['round']>0)].selected_current.mean()),'padding_slots':int((pp.get('status',pd.Series(index=pp.index,dtype=str))=='budget_padding_self_loop').sum())})
 diag=[]
 for lo,hi in [(0,.4),(.4,.5),(.5,.6),(.6,.8),(.8,1.00001)]:
  q=d[(d.nn_similarity>=lo)&(d.nn_similarity<hi)];diag.append({'range':f'[{lo:.1f},{min(hi,1):.1f}{"]" if hi>1 else ")"}','n':len(q),'RF_MAE':float(q.RF_abs_error.mean()) if len(q) else None})
 dg=[]
 # Equal-sized descriptive disagreement bins; no threshold selected from labels.
 d['bin']=pd.qcut(d.disagreement,4,labels=False,duplicates='drop')
 for k,q in d.groupby('bin'):dg.append({'quartile':int(k)+1,'n':len(q),'mean_disagreement':float(q.disagreement.mean()),'RF_MAE':float(q.RF_abs_error.mean()),'GNN_MAE':float(q.GNN_abs_error.mean())})
 # Predeclared rank rules for illustrations. Cases may overlap; never hand-swap.
 final=end.reset_index().merge(c[['smiles','GNN','disagreement']],on='smiles',how='left');final.to_csv(ROOT/'results/final_candidates.csv',index=False)
 cases=[]
 for label,key,asc in [('RF最高','RF',False),('相似性最低','nn_similarity',True),('两模型分歧最大','disagreement',False),('SA最高','SA',False)]:
  row=final.sort_values([key,'smiles'],ascending=[asc,True]).iloc[0].to_dict();row['rule']=label;cases.append(row)
 funnel=[]
 for arm in ['A','B']:
  q=final[final.arm==arm].drop_duplicates('smiles');m=(q.MW>=150)&(q.MW<=600);sa=q.SA<=4.5;qed=q.QED>=.4
  funnel.append({'arm':arm,'unique_final':len(q),'MW_150_600':int(m.sum()),'plus_SA_le_4_5':int((m&sa).sum()),'plus_QED_ge_0_4':int((m&sa&qed).sum())})
 paired=changes.xs('B',level='arm')-changes.xs('A',level='arm');paired.to_csv(ROOT/'results/B_minus_A_changes.csv')
 summary={'audit':audit,'arms':rows,'initial_RF_mean':float(init.RF.mean()),'initial_nn_mean':float(init.nn_similarity.mean()),'B_minus_A_delta_RF_mean':float(paired.RF.mean()),'B_minus_A_delta_nn_mean':float(paired.nn_similarity.mean()),'paired_runs_RF_B_higher':int((paired.RF>0).sum()),'paired_runs_nn_B_lower':int((paired.nn_similarity<0).sum()),'historical_similarity_bins':diag,'historical_disagreement_quartiles':dg,'disagreement_RF_abs_error_spearman':float(spearmanr(d.disagreement,d.RF_abs_error).statistic),'cases':cases,'posthoc_funnel':funnel,'candidate_count':len(c),'verification':'PASS: 36 trajectories, matched budget, monotone greedy RF, no generated measured labels'}
 dump('summary.json',summary);print(json.dumps({k:v for k,v in summary.items() if k not in ['audit','cases']},ensure_ascii=False,indent=2))
if __name__=='__main__':main()
