"""依据e23已复算的附件证据绘制两张补充图；不训练新模型。"""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
COL={'RF':'#3478A5','GNN':'#D58A29'};DARK='#173E35'
plt.rcParams.update({'font.sans-serif':['Hiragino Sans GB','Arial Unicode MS','DejaVu Sans'],'axes.unicode_minus':False,'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
def save(fig,name):fig.savefig(ROOT/'figures'/name,dpi=190,bbox_inches='tight',facecolor='white');plt.close(fig)
def main():
 audit=json.loads((ROOT/'results/e23_attachment_audit.json').read_text());main=pd.read_csv(ROOT/'results/benchmark.csv');fit_n=load_data()[3].size;full=main[main.n_train==fit_n];extra=pd.read_csv(ROOT/'results/attachment_snapshot_metrics.csv');extra=extra[extra.split=='scaffold']
 fig=plt.figure(figsize=(11,7.3),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1.2,1]);axes=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1])]
 for ax,source,title in [(axes[0],full,'协议A：7148拟合 / 128维 / 3种子'),(axes[1],extra,'协议B：附件8698拟合 / 300维 / 5种子')]:
  for i,model in enumerate(['RF','GNN']):
   sub=source[source.model.str.upper()==model];vals=sub.MAE.to_numpy();yy=i+np.linspace(-.13,.13,len(vals));ax.scatter(vals,yy,s=55,color=COL[model],zorder=3);ax.plot([vals.mean(),vals.mean()],[i-.22,i+.22],color=COL[model],lw=2)
  ax.set(xlim=(.57,.85),ylim=(-.55,1.55),yticks=[0,1],yticklabels=['RF','GNN'],xlabel='同一骨架测试集 MAE（越低越好）');ax.set_title(title,fontsize=12,color=DARK);ax.grid(axis='x',alpha=.2)
 ax=fig.add_subplot(gs[1,:]);keys=['main_128dim_3seeds','attachment_300dim_5seeds'];labels=['协议A','协议B']
 for i,key in enumerate(keys):
  r=audit['paired_cluster_bootstrap'][key];lo,hi=r['cluster_bootstrap95'];mean=r['mean_delta_GNN_minus_RF'];ax.errorbar(mean,i,xerr=[[mean-lo],[hi-mean]],fmt='o',capsize=5,color=[COL['GNN'],COL['RF']][i]);ax.text(hi+.004,i,f'{mean:+.3f} [{lo:+.3f}, {hi:+.3f}]',va='center',fontsize=10)
 ax.axvline(0,color='#888',ls='--');ax.set(xlim=(-.015,.30),ylim=(-.6,1.6),yticks=[0,1],yticklabels=labels,xlabel='平均成对绝对误差差值：GNN - RF（正值表示GNN误差更大）');ax.set_title('先对每个分子的种子误差取均值，再按771组测试骨架重抽样2000次',fontsize=11)
 fig.suptitle('补图 A｜两套训练协议分别报告；差值区间针对同一测试集合',fontsize=15,color=DARK)
 fig.supxlabel('上：每点为一次运行，竖线为均值。下：95%条件区间，不含重新训练或更换骨架划分的波动。',fontsize=10,color='#666')
 save(fig,'图09_双协议与配对差值.png')
 c=audit['scaling_correction'];mols=c['molecules'];old=[r['wrong_scale_prediction_reproduced'] for r in mols];new=[r['correct_scale_prediction'] for r in mols];labels=['厄洛替尼','测试集高标签案例','测试集低标签案例'];fig,ax=plt.subplots(figsize=(10,5.4),layout='constrained');x=np.arange(3)
 for offset,values,label,color in [(-.18,old,'错误：由前200个训练分子重算尺度','#A8B1B5'),(.18,new,'修正：使用全部8698训练分子的尺度','#16856B')]:
  bars=ax.bar(x+offset,values,width=.34,label=label,color=color)
  for bar,v in zip(bars,values):ax.text(bar.get_x()+bar.get_width()/2,v+.12,f'{v:.3f}',ha='center',fontsize=11)
 ax.set(xticks=x,xticklabels=labels,ylabel='还原后的预测 pIC50',ylim=(0,12.5));ax.legend(loc='upper center',frameon=False,fontsize=10);ax.set_title('补图 B｜同一300维模型权重，错误的缩放统计量会改变预测值',fontsize=14,color=DARK)
 fig.supxlabel(f'错误：mean {c["wrong_mean"]:.3f}, scale {c["wrong_scale"]:.3f}；训练时实际使用：mean {c["correct_mean"]:.3f}, scale {c["correct_scale"]:.3f}\n此处只修正预测还原，不是重训或模型精度提升试验；厄洛替尼位于训练集中。',fontsize=10,color='#666')
 save(fig,'图10_缩放还原审计.png');print('补充图9、10已生成。')
if __name__=='__main__':main()
