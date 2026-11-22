"""输入特征遮蔽，不是注意力：保留所有原子/键，仅将一个原子的节点向量置零。
输出差值为完整预测-遮蔽预测。零向量可能偏离训练分布，不能解释成删原子的真实化学效应。
"""
from common import *
import copy,torch
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from chemprop.featurizers import SimpleMoleculeMolGraphFeaturizer
from chemprop.data.collate import BatchMolGraph

def main():
 torch.set_num_threads(CONFIG['gnn']['threads']);ckpt=torch.load(ROOT/'models/gnn_full_seed42.pt',map_location='cpu',weights_only=False);assert ckpt['config']==CONFIG
 model=make_model();model.load_state_dict(ckpt['state_dict']);model.eval();audit=json.loads((ROOT/'results/e18_audit.json').read_text());smi=audit['erlotinib_smiles'];mol=Chem.MolFromSmiles(smi)
 graph=SimpleMoleculeMolGraphFeaturizer()(mol);b=BatchMolGraph([graph]);mean,scale=ckpt['mean'],ckpt['scale']
 with torch.no_grad():
  original=float(model(b).item()*scale+mean);rows=[]
  for i,atom in enumerate(mol.GetAtoms()):
   masked=copy.deepcopy(b);masked.V[i]=0
   pred=float(model(masked).item()*scale+mean);rows.append({'atom':i,'symbol':atom.GetSymbol(),'prediction_masked':pred,'delta_pIC50':original-pred})
 imp=np.array(json.loads((ROOT/'results/rf_full_seed42_importance.json').read_text()));gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)
 from rdkit.Chem import rdFingerprintGenerator as fg
 extra=fg.AdditionalOutput();extra.AllocateBitInfoMap();fp=gen.GetFingerprint(mol,additionalOutput=extra);info=extra.GetBitInfoMap();bits=sorted(fp.GetOnBits(),key=lambda b:imp[b],reverse=True)[:6]
 df,_,_,fit,val,test=load_data();identities=[Chem.MolToSmiles(Chem.MolFromSmiles(s)) for s in df.smiles];idx=[i for i,s in enumerate(identities) if s==smi]
 source={k:[i for i in idx if i in set(v)] for k,v in [('fit',fit),('validation',val),('test',test)]}
 out={'smiles':smi,'base_prediction':original,'dataset_membership':source,'atoms':rows,'rf_present_top_bits':[{'bit':int(b),'global_impurity_importance':float(imp[b]),'environments':[list(e) for e in info[b]]} for b in bits], 'limitations':['node-feature zeroing is out-of-distribution sensitivity, not chemical deletion','no attention layer in this model','RF global impurity importance is neither local attribution nor causal importance','same hashed bit may map to multiple environments']}
 dump(ROOT/'results/e20_sensitivity.json',out);print(json.dumps(out,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
