from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

_FP = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def select(candidates, count, max_similarity=1.0):
    eligible = sorted((c for c in candidates if c.score.get("complete") and c.score.get("feasible")),
                      key=lambda c:(-c.score["fitness"],c.smiles))
    chosen, fps, seen = [], [], set()
    for c in eligible:
        if c.smiles in seen: continue
        fp = _FP.GetFingerprint(Chem.MolFromSmiles(c.smiles))
        if any(DataStructs.TanimotoSimilarity(fp,other)>max_similarity for other in fps): continue
        chosen.append(c); fps.append(fp); seen.add(c.smiles)
        if len(chosen)>=count: break
    return chosen if count > 0 else []
