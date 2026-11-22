"""主实验A的具体读图：128维GNN节点遮蔽＋同次RF全局Top20位的环境中心计数。
不训练新模型。GNN采用delta=masked-full，符号与e20保存的full-masked相反。
中心投票与环境覆盖都存档，不能把中心投票当作完整子结构局部归因。
"""
from common import *
import io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator,rdmolops
from rdkit.Chem.Draw import rdMolDraw2D
plt.rcParams.update({'font.sans-serif':['Hiragino Sans GB','Arial Unicode MS','DejaVu Sans'],'axes.unicode_minus':False,'font.size':11,'axes.spines.top':False,'axes.spines.right':False})

def render_case():
 s=json.loads((ROOT/'results/e20_sensitivity.json').read_text());mol=Chem.MolFromSmiles(s['smiles']);delta=-np.array([r['delta_pIC50'] for r in s['atoms']]);imp=np.array(json.loads((ROOT/'results/rf_full_seed42_importance.json').read_text()))
 gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False);ao=rdFingerprintGenerator.AdditionalOutput();ao.AllocateBitInfoMap();gen.GetFingerprint(mol,additionalOutput=ao);info=ao.GetBitInfoMap();bits=np.argsort(-imp,kind='stable')[:20];center_votes=np.zeros(mol.GetNumAtoms(),int);coverage=np.zeros_like(center_votes);environments=[]
 for bit in bits:
  for center,radius in info.get(int(bit),[]):
   center_votes[center]+=1;atoms={int(center)}
   for bond in Chem.FindAtomEnvironmentOfRadiusN(mol,int(radius),int(center)):
    b=mol.GetBondWithIdx(int(bond));atoms.update([b.GetBeginAtomIdx(),b.GetEndAtomIdx()])
   for i in atoms:coverage[i]+=1
   environments.append({'bit':int(bit),'center':int(center),'radius':int(radius),'covered_atom_indices':sorted(atoms)})
 top=np.argsort(-np.abs(delta),kind='stable')[:10];voted=set(np.where(center_votes>0)[0].tolist());overlap=sorted(voted&set(top.tolist()))
 out={'protocol':'A: 7148 fit, 128 hidden, seed42; two models share training subset','smiles':s['smiles'],'base_prediction':s['base_prediction'],'in_fit_pool':bool(s['dataset_membership']['fit']),'atom_numbering':'zero-based canonical molecule atom indices','delta_definition':'prediction_masked - prediction_full','threshold_for_descriptive_count':.02,'n_abs_delta_gt_threshold':int((np.abs(delta)>.02).sum()),'n_atoms':mol.GetNumAtoms(),'max_abs_delta':float(np.max(np.abs(delta))),'most_sensitive_atom':int(top[0]),'rf_top_global_bits':bits.tolist(),'rf_center_votes':center_votes.tolist(),'rf_environment_coverage_votes':coverage.tolist(),'rf_center_voted_atoms':sorted(voted),'top10_sensitive_atoms':top.tolist(),'top10_and_center_votes_overlap':overlap,'atom_rows':[{'index':i,'element':mol.GetAtomWithIdx(i).GetSymbol(),'delta_masked_minus_full':float(delta[i]),'rf_center_votes':int(center_votes[i]),'rf_environment_coverage_votes':int(coverage[i])} for i in range(mol.GetNumAtoms())],'bit_environments':environments,'limitations':['zeroing features may be out of distribution','center votes are not local attribution or causal importance','top10 overlap is a descriptive mapping, not a test of chemical mechanism','erlotinib is in fit pool, not external validation']}
 dump(ROOT/'results/e25_case_readout.json',out)
 # Draw atom annotations and a signed color scale.
 for a in mol.GetAtoms():a.SetProp('atomNote',str(a.GetIdx()))
 limit=max(float(np.abs(delta).max()),1e-8);cm=plt.get_cmap('coolwarm');colors={i:tuple(cm(.5+.5*v/limit)[:3]) for i,v in enumerate(delta)};drawer=rdMolDraw2D.MolDraw2DCairo(760,460);drawer.DrawMolecule(mol,highlightAtoms=list(range(mol.GetNumAtoms())),highlightBonds=[],highlightAtomColors=colors);drawer.FinishDrawing();im=Image.open(io.BytesIO(drawer.GetDrawingText()))
 fig,ax=plt.subplots(1,2,figsize=(12,6),gridspec_kw={'width_ratios':[1.15,1]},layout='constrained');ax[0].imshow(im);ax[0].axis('off');ax[0].set_title(f'厄洛替尼｜主实验A预测 pIC50 {s["base_prediction"]:.3f}',fontsize=12)
 bar=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-limit,limit),cmap=cm),ax=ax[0],orientation='horizontal',fraction=.07,pad=.02);bar.set_label('遮蔽后 - 完整预测 / pIC50',fontsize=10)
 selected=top[::-1];vals=np.abs(delta[selected]);labels=[f'{mol.GetAtomWithIdx(int(i)).GetSymbol()}{i}｜RF中心票 {center_votes[i]}' for i in selected];bars=ax[1].barh(labels,vals,color='#D58A29')
 for b,v in zip(bars,vals):ax[1].text(v+.008,b.get_y()+b.get_height()/2,f'{v:.3f}',va='center',fontsize=10)
 ax[1].set(xlim=(0,limit*1.25),xlabel='GNN绝对预测变化 / pIC50');ax[1].set_title('GNN最敏感的10个原子，旁注RF映射参考',fontsize=12)
 fig.suptitle('图 7｜打开黑盒：输入敏感性与指纹位映射，是两种读图方式',fontsize=15,color='#173E35')
 fig.supxlabel('原子编号从0开始。RF中心票来自全局重要性Top20位在本分子中的环境中心，不是原子活性贡献。\n红色：遮蔽后预测升高；蓝色：遮蔽后预测降低。零向量可能超出训练分布；两种指标不能直接作机制比较。',fontsize=9.5,color='#666')
 fig.savefig(ROOT/'figures/图07_敏感性与指纹位.png',dpi=190,bbox_inches='tight',facecolor='white');plt.close(fig)
 print(json.dumps({k:out[k] for k in ['base_prediction','n_abs_delta_gt_threshold','max_abs_delta','most_sensitive_atom','rf_center_voted_atoms','top10_and_center_votes_overlap']},ensure_ascii=False))
 return out
if __name__=='__main__':render_case()
