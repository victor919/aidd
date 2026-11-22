"""固定拆分、版本和图特征审计；不从测试指标选择配置。"""
from common import *
import importlib.metadata,platform,torch
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
from chemprop.data import MoleculeDataset,MoleculeDatapoint
from chemprop.featurizers import SimpleMoleculeMolGraphFeaturizer

def main():
 df,X,y,fit,val,test=load_data();splits={'fit':fit.tolist(),'validation':val.tolist(),'test':test.tolist()};dump(ROOT/'results/protocol_indices.json',splits)
 feat=SimpleMoleculeMolGraphFeaturizer();m=Chem.MolFromSmiles('COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC');graph=feat(m)
 isomers=['CC(C)Cc1ccc([C@H](C)C(=O)O)cc1','CC(C)Cc1ccc([C@@H](C)C(=O)O)cc1'];graphs=[feat(Chem.MolFromSmiles(s)) for s in isomers]
 rows={}
 for chiral in [False,True]:
  gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=chiral)
  rows[str(chiral)]=DataStructs.TanimotoSimilarity(*[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in isomers])
 out={'n_total':len(df),'n_fit':len(fit),'n_validation':len(val),'n_test':len(test),'graph_shapes':{'atoms':list(graph.V.shape),'directed_bonds':list(graph.E.shape)},'erlotinib_smiles':Chem.MolToSmiles(m),'chirality':{'fingerprint_tanimoto':rows,'chemprop_atom_features_equal':bool(np.array_equal(graphs[0].V,graphs[1].V))},'environment':{'python':platform.python_version(),'platform':platform.platform(),'device':'cpu','mps_available':bool(torch.backends.mps.is_available()),'cuda_available':bool(torch.cuda.is_available())},'versions':{p:importlib.metadata.version(p) for p in ['chemprop','torch','lightning','rdkit','numpy','pandas','scipy','scikit-learn','matplotlib']},'inputs_sha256':{p.name:sha(p) for p in (ROOT/'data').glob('*')},'config_sha256':sha(ROOT/'config.json'),'protocol':CONFIG['protocol']}
 assert X.shape==(len(df),2048)
 # 检查固定快照的特征确实与所用生成规则一致，而非只相信文件名。
 gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=False)
 for i in np.random.default_rng(42).choice(len(df),20,replace=False):
  a=np.zeros(2048,dtype=np.uint8);DataStructs.ConvertToNumpyArray(gen.GetFingerprint(Chem.MolFromSmiles(df.iloc[i].smiles)),a);assert np.array_equal(a,X[i])
 dump(ROOT/'results/e18_audit.json',out);print(json.dumps(out,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
