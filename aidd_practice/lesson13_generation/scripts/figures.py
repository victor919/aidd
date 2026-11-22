from experiment import *
import os,io,textwrap
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'work/mpl'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
from rdkit.Chem.Draw import rdMolDraw2D
plt.rcParams.update({'font.sans-serif':['Hiragino Sans GB','Arial Unicode MS','DejaVu Sans'],'axes.unicode_minus':False,'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
COL={'A':'#3E7E9D','B':'#D7892A'};LABEL={'A':'A 随机保留','B':'B 按RF贪心保留'}
def save(fig,n):fig.savefig(ROOT/f'figures/{n}.png',bbox_inches='tight',facecolor='white');plt.close(fig)
def molecule(s,legend='',size=(650,370)):
 m=Chem.MolFromSmiles(s);draw=rdMolDraw2D.MolDraw2DCairo(*size);draw.DrawMolecule(m);draw.FinishDrawing();return Image.open(io.BytesIO(draw.GetDrawingText()))
def boxes(title,texts,name):
 fig,ax=plt.subplots(figsize=(10,4.5));ax.axis('off');ax.set_title(title,color='#15594E',fontsize=17,pad=18)
 for i,(h,b) in enumerate(texts):
  x=.02+i*.33;ax.text(x,.62,h,transform=ax.transAxes,fontsize=14,color='#15594E',weight='bold',bbox=dict(facecolor='#EAF4F0',edgecolor='#B8D8CC',pad=13));ax.text(x,.24,b,transform=ax.transAxes,fontsize=11,linespacing=1.8)
  if i<2:ax.text(x+.265,.62,'→',transform=ax.transAxes,fontsize=25,color='#D7892A')
 save(fig,name)
def main():
 s=pd.read_csv(ROOT/'results/states.csv');d=pd.read_csv(ROOT/'results/historical_test_diagnostic.csv');c=pd.read_csv(ROOT/'results/candidates.csv');p=pd.read_csv(ROOT/'results/proposals.csv');r=json.loads((ROOT/'results/summary.json').read_text())
 # 对照开环与闭环，反馈箭头从保留状态返回下一轮提议。
 fig,ax=plt.subplots(figsize=(10,5.6));ax.set(xlim=(0,10),ylim=(0,6));ax.axis('off');ax.set_title('图1｜从一次评分，到反复按评分选择',color='#15594E',fontsize=17)
 for y,label,texts,color in [(4.6,'开环：给定分子库',['给定分子','冻结评分器','一次排序'],'#3E7E9D'),(2,'闭环：更新下一轮状态',['当前结构 → 提议','冻结评分器','选择并保留'],'#D7892A')]:
  ax.text(.2,y+.75,label,color=color,fontsize=13,weight='bold')
  for x,t in zip([1.7,5,8.3],texts):ax.text(x,y,t,ha='center',va='center',fontsize=13,color='white',bbox=dict(boxstyle='round,pad=.65',facecolor=color,edgecolor='none'))
  for x1,x2 in [(2.8,3.9),(6.1,7.2)]:ax.annotate('',xy=(x2,y),xytext=(x1,y),arrowprops=dict(arrowstyle='->',color=color,lw=2))
 ax.plot([8.3,8.3,1.7],[1.55,.65,.65],color='#D7892A',lw=2)
 ax.annotate('',xy=(1.7,1.55),xytext=(1.7,.65),arrowprops=dict(arrowstyle='->',color='#D7892A',lw=2))
 ax.text(5,.85,'保留结构成为下一轮起点；RF不更新',ha='center',fontsize=11,color='#9A651E')
 save(fig,'图01_反馈环')
 a=pd.read_csv(ROOT/'results/attempts.csv');q=a[a.status=='accepted_proposal'].iloc[0];fig,ax=plt.subplots(1,3,figsize=(11,4.3));
 for z,key,label in zip(ax,['parent','donor_fragment','candidate'],['父分子','替入片段（带连接点）','替换后的候选']):z.imshow(molecule(q[key]));z.axis('off');z.set_title(label)
 fig.suptitle('图2｜一次真实动作：在BRICS断键处替换单连接点片段',color='#15594E',fontsize=15);save(fig,'图02_替换动作')
 for metric,title,name,ylabel in [('RF','图3｜追求高分是否有效：看状态分数随查询预算变化','图03_分数轨迹','RF预测 pIC50'),('nn_similarity','图4｜分数之外：搜索与拟合池的最近邻相似性','图04_相似性轨迹','最近邻Tanimoto（参照仅7148拟合分子）')]:
  fig,ax=plt.subplots(figsize=(9,5))
  for arm in ['A','B']:
   q=s[s.arm==arm];g=q.groupby('budget')[metric].agg(['mean','std']);ax.plot(g.index,g['mean'],'o-',color=COL[arm],label=LABEL[arm]);ax.fill_between(g.index,g['mean']-g['std'],g['mean']+g['std'],color=COL[arm],alpha=.13)
  ax.set(xlabel='累计候选评分请求（含缓存命中）',ylabel=ylabel,title=title);ax.legend();ax.grid(alpha=.18);fig.text(.5,-.015,'均值及18条轨迹的标准差：6个共享起点 × 3个种子；不是独立项目置信区间。',ha='center',fontsize=10,color='#666');save(fig,name)
 fig,axs=plt.subplots(1,2,figsize=(11,4.7));bins=r['historical_disagreement_quartiles'];axs[0].bar([x['quartile'] for x in bins],[x['RF_MAE'] for x in bins],color='#3E7E9D');axs[0].set(xlabel='RF/GNN分歧四分位（由低到高）',ylabel='历史测试集RF平均绝对误差',title='有标签数据：分歧能否提示错误？')
 for arm in ['A','B']:
  q=s[(s.arm==arm)&(s['round']==10)].merge(c[['smiles','disagreement']],on='smiles');axs[1].scatter(q.RF,q.disagreement,color=COL[arm],alpha=.65,label=LABEL[arm])
 axs[1].set(xlabel='终点RF预测pIC50',ylabel='|RF预测 - GNN预测|',title='生成终点：只报告分歧，没有真值');axs[1].legend();fig.suptitle('图5｜有标签的诊断与无标签的候选，分开报告',color='#15594E',fontsize=15);save(fig,'图05_分歧诊断')
 fig,axs=plt.subplots(2,2,figsize=(11,8))
 for ax,x in zip(axs.flat,r['cases']):
  ax.imshow(molecule(x['smiles']));ax.axis('off');ax.set_title(f"{x['rule']}｜{x['arm']} seed{x['seed']}｜{'拟合池已有' if x['seen_fit'] else '拟合池未见'}\nRF {x['RF']:.2f} / GNN {x['GNN']:.2f} / NN {x['nn_similarity']:.2f}\nMW {x['MW']:.0f} / SA {x['SA']:.2f}",fontsize=11)
 fig.suptitle('图6｜按固定排名规则选例：重复出现也不替换',color='#15594E',fontsize=15);save(fig,'图06_候选卡片')
 fig,axs=plt.subplots(1,2,figsize=(11,4.7));bins=r['historical_similarity_bins'];axs[0].bar(range(len(bins)),[x['RF_MAE'] or 0 for x in bins],color='#3E7E9D');axs[0].set(xticks=range(len(bins)),xticklabels=[x['range'] for x in bins],xlabel='最近邻相似性分箱',ylabel='历史测试集RF MAE',title='有标签：误差与相似性')
 for i,x in enumerate(bins):axs[0].text(i,(x['RF_MAE'] or 0)+.02,f"n={x['n']}",ha='center',fontsize=9)
 for arm in ['A','B']:
  q=s[(s.arm==arm)&(s['round']==10)];axs[1].scatter(q.nn_similarity,q.RF,color=COL[arm],alpha=.65,label=LABEL[arm])
 axs[1].set(xlabel='最近邻相似性',ylabel='RF预测pIC50',title='无标签：终点的预测与结构位置');axs[1].legend();fig.suptitle('图7｜相似性是诊断线索，不是候选准确性的证明',color='#15594E',fontsize=15);save(fig,'图07_相似性与分数')
 boxes('图8｜从第9篇到第13篇：把候选和证据边界一起交付',[('9—10 结构与姿势','对象、口袋、摆放\n结构假设分别审查'),('11—12 活性预测','划分、基线、表示\n先检查评分器'),('13 生成与筛选','提出结构、记录搜索\n候选进入下一轮验证')],'图08_系列位置')
 # 补图A：用本篇数据展示动作与选择共同伴随的属性变化。
 fig,axs=plt.subplots(1,3,figsize=(12,4.5))
 for ax,key,title in zip(axs,['MW','SA','QED'],['分子量','SAScore（启发式）','QED（描述符组合）']):
  for arm in ['A','B']:
   q=s[s.arm==arm].groupby('budget')[key].mean();ax.plot(q.index,q.values,'o-',color=COL[arm],label=LABEL[arm])
  ax.set(xlabel='候选评分请求',ylabel=title);ax.grid(alpha=.18)
 axs[0].legend();fig.suptitle('补图A｜动作与选择的属性轨迹：不只盯着预测分数',color='#15594E',fontsize=15)
 fig.text(.5,-.015,'每臂18条轨迹均值。属性变化不能直接解释为活性变化或确定合成成本。',ha='center',fontsize=10,color='#666');save(fig,'图09_动作与属性变化')
 print('9 figures generated')
if __name__=='__main__':main()
