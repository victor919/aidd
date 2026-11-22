"""借鉴附件的读图组织，仅复用本篇固定结果补充诊断，不重训不搜索。"""
from experiment import *
def main():
 s=pd.read_csv(ROOT/'results/states.csv');d=pd.read_csv(ROOT/'results/historical_test_diagnostic.csv');start=s[s['round']==0].set_index(['seed','start_index','arm']);end=s[s['round']==10].set_index(['seed','start_index','arm']);delta=end[['RF','nn_similarity','MW','SA','QED']]-start[['RF','nn_similarity','MW','SA','QED']];pairs=delta.xs('B',level='arm')-delta.xs('A',level='arm')
 seeds=[]
 for seed,p in pairs.groupby(level='seed'):seeds.append({'seed':int(seed),'pairs':len(p),'RF_delta_B_minus_A_mean':float(p.RF.mean()),'NN_delta_B_minus_A_mean':float(p.nn_similarity.mean())})
 mask=(d.disagreement<.5)&(d.RF_abs_error>2)&(d.GNN_abs_error>2);examples=d.loc[mask].sort_values(['true','index'],ascending=[False,True]).head(3)
 maxerr=np.maximum(d.RF_abs_error,d.GNN_abs_error);minerr=np.minimum(d.RF_abs_error,d.GNN_abs_error)
 dg=[]
 for k,idx in d.groupby(pd.qcut(d.disagreement,4,labels=False)).groups.items():
  q=d.loc[idx];dg.append({'quartile':int(k)+1,'n':len(q),'max_error_median':float(np.median(maxerr.loc[idx])),'min_error_median':float(np.median(minerr.loc[idx])),'max_error_gt1_fraction':float((maxerr.loc[idx]>1).mean())})
 stats=[]
 for arm in ['A','B']:
  q=delta.xs(arm,level='arm');stats.append({'arm':arm,**{col+'_change_mean':float(q[col].mean()) for col in q}})
 out={'source':'Own lesson13 frozen results; attachment supplies questions, not numbers','seed_pairs':seeds,'changes':stats,'historical_both_wrong_rule':'|RF-GNN| < 0.5 and both absolute errors > 2','historical_both_wrong_count':int(mask.sum()),'historical_both_wrong_examples':examples[['index','true','RF','GNN']].to_dict('records'),'historical_disagreement_bins':dg,'spearman_d_vs_max_error':float(spearmanr(d.disagreement,maxerr).statistic),'spearman_d_vs_min_error':float(spearmanr(d.disagreement,minerr).statistic),'triangle_inequality_max_violation':float(np.maximum(d.disagreement-2*maxerr,0).max()),'no_new_model_fits':True}
 dump('attachment_readout.json',out);print(json.dumps(out,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
