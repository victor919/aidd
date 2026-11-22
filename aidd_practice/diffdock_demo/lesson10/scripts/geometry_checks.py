"""NumPy-only 12-point toy: rigid-motion equivariance, reflection and handedness.
No RDKit fallback or identity claim; not a DiffDock reproduction.
"""
import json
from pathlib import Path
import numpy as np
rng=np.random.default_rng(20261031);X=rng.normal(size=(12,3));errors=[]
def f(x):
 diff=x[:,None,:]-x[None,:,:];w=np.exp(-np.sum(diff*diff,axis=-1))
 return (w[:,:,None]*diff).sum(axis=1)
for _ in range(24):
 R,_=np.linalg.qr(rng.normal(size=(3,3)))
 if np.linalg.det(R)<0:R[:,0]*=-1
 errors.append(float(np.max(np.abs(f(X@R.T)-f(X)@R.T))))
M=np.diag([-1.,1.,1.]);P=np.array([[0.,0.,0.],[1.,0.,0.],[.5,.9,0.],[.5,.3,.8]])
def D(x):return np.linalg.norm(x[:,None,:]-x[None,:,:],axis=-1)
def signed_volume(x):return float(np.linalg.det(np.stack([x[1]-x[0],x[2]-x[0],x[3]-x[0]]))/6)
r={'input':'12 random points, seed=20261031','n_rotations':24,'max_rotation_equivariance_error':max(errors),'reflection_equivariance_error':float(np.max(np.abs(f(P@M.T)-f(P)@M.T))),'mirror_distance_matrix_error':float(np.max(np.abs(D(P)-D(P@M.T)))),'signed_volume_original':signed_volume(P),'signed_volume_mirror':signed_volume(P@M.T),'note':'Fixed ordered vertices. Signed volume is rotation invariant and changes sign under reflection; not reflection-invariant. Raw vector change alone does not prove chirality-sensitive scalar prediction.'}
assert max(errors)<1e-12 and r['reflection_equivariance_error']<1e-12
assert np.isclose(r['signed_volume_original'],-r['signed_volume_mirror'])
root=Path(__file__).resolve().parents[1];(root/'audit/geometry_checks.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
