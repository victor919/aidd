"""汇总24次拟合并从真实结果绘制8张图。图1/2部分为方法示意；其余来源见图注。"""
from common import *
import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch
from rdkit import Chem,DataStructs
from rdkit.Chem import Draw,rdFingerprintGenerator
COL={'RF':'#3478A5','GNN':'#D58A29'};DARK='#173E35';GREEN='#14856A'
plt.rcParams.update({'font.sans-serif':['Hiragino Sans GB','Arial Unicode MS','DejaVu Sans'],'axes.unicode_minus':False,'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
def save(fig,name):
 fig.savefig(FIG/name,dpi=180,bbox_inches='tight',facecolor='white');plt.close(fig)
def box(ax,x,y,w,h,txt,fc='#EFF7F4',ec=GREEN,fs=11):
 ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.014',facecolor=fc,edgecolor=ec,transform=ax.transAxes));ax.text(x+w/2,y+h/2,txt,ha='center',va='center',transform=ax.transAxes,fontsize=fs,color=DARK,linespacing=1.7)
def arrow(ax,a,b):ax.add_patch(FancyArrowPatch(a,b,arrowstyle='->',mutation_scale=14,color=DARK,transform=ax.transAxes))
def molecular_image(mol,highlight=None,colors=None,legend=''):
 from rdkit.Chem.Draw import rdMolDraw2D
 from PIL import Image
 drawer=rdMolDraw2D.MolDraw2DCairo(600,330)
 drawer.DrawMolecule(mol,legend=legend,highlightAtoms=highlight or [],highlightBonds=[],highlightAtomColors=colors or {})
 drawer.FinishDrawing()
 return Image.open(io.BytesIO(drawer.GetDrawingText())).convert('RGB')
def draw_neighborhood_figure(mol):
 """用真实拓扑距离统一绘制0/1/2/3跳；两行布局便于公众号阅读。"""
 from rdkit.Chem import rdDepictor
 mol=Chem.Mol(mol);rdDepictor.Compute2DCoords(mol)
 # 选择连接环系的NH氮原子，避免仅依赖SMILES中的原子顺序。
 candidates=[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum()==7 and a.GetTotalNumHs()==1 and a.GetDegree()==2]
 assert len(candidates)==1
 center=candidates[0];dist=Chem.GetDistanceMatrix(mol);total=mol.GetNumHeavyAtoms()
 assert total==mol.GetNumAtoms()==29
 counts=[int((dist[center]<=k).sum()) for k in range(4)]
 assert counts==[1,3,7,12]
 fig,axes=plt.subplots(2,2,figsize=(10,10.2))
 fig.subplots_adjust(top=.84,bottom=.16,left=.035,right=.965,wspace=.09,hspace=.20)
 titles=['0跳｜只标出中心','1跳｜加入直接邻居','2跳｜继续向外扩展','3跳｜覆盖更远的局部环境']
 for k,ax in enumerate(axes.flat):
  atoms=[i for i in range(total) if dist[center,i]<=k]
  colors={i:((.98,.67,.25) if i==center else (.48,.78,.67)) for i in atoms}
  ax.imshow(molecular_image(mol,atoms,colors));ax.axis('off')
  ax.set_title(titles[k],fontsize=14,color=DARK,pad=30,fontweight='bold')
  ax.text(.5,1.06,f'覆盖 {counts[k]} / {total} 个重原子',transform=ax.transAxes,ha='center',va='bottom',fontsize=12,color=GREEN)
 fig.text(.5,.963,'图 1｜同一分子，逐步扩大的拓扑邻域',ha='center',fontsize=19,color=DARK,fontweight='bold')
 fig.text(.5,.920,'厄洛替尼 · 固定NH氮原子为中心 · 0跳 → 1跳 → 2跳 → 3跳',ha='center',fontsize=12,color='#5B7471')
 fig.text(.5,.105,'橙色：固定中心原子    绿色：该邻域内其余原子    未高亮：邻域外',ha='center',fontsize=11,color=DARK)
 fig.text(.5,.070,'“跳”按跨越的化学键数计算；高亮表示距中心的最短路径不超过当前跳数。',ha='center',fontsize=10.5,color='#647571')
 fig.text(.5,.036,'拓扑示意，不是注意力或活性贡献图；0/1/2/3跳不直接等同于D-MPNN的depth参数。',ha='center',fontsize=10,color='#647571')
 save(fig,'图01_分子图与邻域.png')
 dump(ROOT/'results/neighborhood_audit.json',{'molecule':'erlotinib','canonical_smiles':Chem.MolToSmiles(mol),'center_atom_index':center,'center_element':'N','center_total_hydrogens':1,'n_heavy_atoms':total,'radii':[0,1,2,3],'covered_atom_counts':counts,'atom_indices':{str(k):[int(i) for i in range(total) if dist[center,i]<=k] for k in range(4)},'meaning':'Topological shortest-path neighborhoods, not D-MPNN depth or attention.'})

def main():
 df,X,y,fit,val,test=load_data();audit=json.loads((ROOT/'results/e18_audit.json').read_text());smi=audit['erlotinib_smiles'];mol=Chem.MolFromSmiles(smi)
 # 1 真实拓扑距离生成的四阶段整合图
 draw_neighborhood_figure(mol)
 # 2 protocol
 fig,ax=plt.subplots(figsize=(11,4.6));ax.axis('off')
 box(ax,.015,.62,.2,.22,'10468 个分子\n固定外层骨架划分');arrow(ax,(.23,.73),(.29,.73))
 box(ax,.31,.62,.27,.22,f'外层训练 8698\n再按骨架分验证集');arrow(ax,(.59,.73),(.66,.73))
 box(ax,.68,.62,.29,.22,f'拟合池 {len(fit)}｜验证 {len(val)}\n测试 {len(test)} 独立保留')
 box(ax,.04,.12,.4,.25,'每个种子：嵌套抽样\n500 → 2000 → 5000 → 7148',fc='#F1F6FB',ec=COL['RF']);arrow(ax,(.43,.245),(.52,.245))
 box(ax,.55,.12,.4,.25,'相同分子子集 → RF / GNN\n3种子 × 4样本量 × 2模型 = 24次拟合',fc='#FFF7EB',ec=COL['GNN'])
 ax.set_title('图 2｜比较模型之前，先把分子与验证规则对齐',fontsize=15,color=DARK,pad=16)
 fig.text(.5,.015,'所有学习曲线共用验证/测试集；未搜索超参数。GNN用验证集选轮次，RF固定500树。',ha='center',fontsize=10,color='#666')
 save(fig,'图02_对照协议.png')
 b=pd.read_csv(ROOT/'results/benchmark.csv');assert len(b)==24,'完整学习曲线需完成24次拟合'
 stats=[]
 for (n,model),sub in b.groupby(['n_train','model']):
  assert len(sub)==3
  stats.append({'n_train':int(n),'model':model,**{f'{v}_{a}':float(getattr(sub[v],a)()) for v in ['MAE','RMSE','R2','spearman','fit_seconds'] for a in ['mean','std']}})
 full=b[b.n_train==len(fit)];pivot=full.pivot(index='seed',columns='model',values='MAE');delta=(pivot.GNN-pivot.RF)
 p=pd.read_csv(ROOT/f'results/predictions_n{len(fit)}_seed42.csv',keep_default_na=False)
 # compare scatter same molecular sample
 fig,axes=plt.subplots(1,2,figsize=(10,5),sharex=True,sharey=True,layout='constrained')
 for a,model in zip(axes,['RF','GNN']):
  a.scatter(p.true,p[model],s=8,alpha=.4,color=COL[model],rasterized=True);a.plot([1,11.5],[1,11.5],color='#888',ls='--');a.set(xlim=(1,11.5),ylim=(1,11.5),xlabel='聚合标签 pIC50',ylabel='预测 pIC50');a.set_aspect('equal');m=metrics(p.true.to_numpy(),p[model].to_numpy());a.set_title(f"{model}｜MAE {m['MAE']:.3f}｜R2 {m['R2']:.3f}",fontsize=12)
 fig.suptitle('图 4｜同一1770分子测试集、同一7148分子拟合集（seed 42）',fontsize=14,color=DARK);save(fig,'图04_同测试集散点.png')
 # 4 curves sd not CI
 fig,ax=plt.subplots(figsize=(9,5),layout='constrained')
 for model in ['RF','GNN']:
  rows=sorted([r for r in stats if r['model']==model],key=lambda r:r['n_train']);xx=np.array([r['n_train'] for r in rows]);mean=np.array([r['MAE_mean'] for r in rows]);sd=np.array([r['MAE_std'] for r in rows]);ax.plot(xx,mean,'o-',color=COL[model],label=model,lw=2);ax.fill_between(xx,mean-sd,mean+sd,color=COL[model],alpha=.13)
 ax.set(xlabel='实际拟合分子数',ylabel='测试 MAE / pIC50',xticks=[500,2000,5000,len(fit)]);ax.legend();ax.grid(alpha=.2);ax.set_title('图 5｜学习曲线：3个种子的均值 ± 标准差',fontsize=14,color=DARK)
 fig.text(.5,-.02,'阴影是3次运行的标准差，不是95%置信区间；不包含重新划分骨架测试集的波动。',ha='center',fontsize=10,color='#666');save(fig,'图05_学习曲线.png')
 #5 validation trajectories
 fig,ax=plt.subplots(figsize=(9,4.7),layout='constrained')
 for seed in CONFIG['seeds']:
  hist=json.loads((ROOT/f'results/history_n{len(fit)}_seed{seed}.json').read_text());values=[r['validation_MAE'] for r in hist];epochs=[r['epoch'] for r in hist];line=ax.plot(epochs,values,label=f'seed {seed}')[0];k=int(np.argmin(values));ax.scatter(epochs[k],values[k],s=70,marker='*',color=line.get_color())
 ax.set(xlabel='训练轮次',ylabel='验证 MAE / pIC50');ax.legend();ax.grid(alpha=.2);ax.set_title('图 3｜只根据验证集选择轮次；星号为该次最佳轮次',fontsize=14,color=DARK);save(fig,'图03_验证轨迹.png')
 #6 resource local wallclock
 fig,ax=plt.subplots(figsize=(9,4.7),layout='constrained')
 for model in ['RF','GNN']:
  rows=sorted([r for r in stats if r['model']==model],key=lambda r:r['n_train']);ax.plot([r['n_train'] for r in rows],[r['fit_seconds_mean'] for r in rows],'o-',color=COL[model],label=model)
 ax.set(xlabel='实际拟合分子数',ylabel='记录的训练流程平均墙钟秒数',xticks=[500,2000,5000,len(fit)]);ax.legend();ax.grid(alpha=.2);ax.set_title('图 6｜这台CPU上的成本，不是跨硬件速度排名',fontsize=14,color=DARK)
 fig.text(.5,-.02,'RF仅fit，4个并行作业；GNN含逐轮验证，1个Torch线程。共同预处理与最终测试不计时。',ha='center',fontsize=10,color='#666');save(fig,'图06_训练成本.png')
 #7 主实验A的具体读图：统一复用独立可运行函数
 from e25_case_readout import render_case
 render_case()
 #8 paired errors + cases chosen without hand picking
 e_rf=np.abs(p.RF-p.true);e_gnn=np.abs(p.GNN-p.true);difference=e_gnn-e_rf
 fig,ax=plt.subplots(1,2,figsize=(10,4.5),layout='constrained');ax[0].scatter(e_rf,e_gnn,s=9,alpha=.35,color=GREEN);lim=max(float(e_rf.max()),float(e_gnn.max()));ax[0].plot([0,lim],[0,lim],ls='--',color='#777');ax[0].set(xlabel='RF绝对误差',ylabel='GNN绝对误差',title='每点为同一个测试分子')
 ax[1].hist(difference,bins=45,color=COL['GNN'],alpha=.85);ax[1].axvline(0,color='#555',ls='--');ax[1].set(xlabel='GNN绝对误差 - RF绝对误差',ylabel='分子数',title='负值：这一个分子上GNN误差更小')
 fig.suptitle('图 8｜平均分之外，两个模型会在不同分子上犯错',fontsize=14,color=DARK);save(fig,'图08_成对误差.png')
 # Predeclared quantile-rank representative cases, no mechanism inference.
 cases=[]
 for label,index in [('GNN误差优势最大',int(np.argmin(difference))),('RF误差优势最大',int(np.argmax(difference))),('两者都较好',int(np.argmin(e_rf+e_gnn))),('两者共同困难',int(np.argmax(e_rf+e_gnn)))]:
  row=p.iloc[index];cases.append({'rule':label,'index':int(row['index']),'smiles':row.smiles,'true':float(row.true),'RF':float(row.RF),'GNN':float(row.GNN),'RF_abs_error':float(e_rf.iloc[index]),'GNN_abs_error':float(e_gnn.iloc[index])})
 out={'statistics':stats,'full_seed_paired_delta_GNN_minus_RF':{str(k):float(v) for k,v in delta.items()},'test_n':len(test),'fit_pool_n':len(fit),'val_n':len(val),'completed_fits':len(b),'seed42_fraction_GNN_lower_abs_error':float((e_gnn<e_rf).mean()),'cases':cases,'uncertainty_scope':'Sample SD over 3 seeds; fixed scaffold split. Includes subset sampling and initialization, not scaffold-split variability.'}
 dump(ROOT/'results/summary.json',out);print(json.dumps(out,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
