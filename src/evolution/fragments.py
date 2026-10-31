"""Persistent BRICS bank. Association is separate from replacement evidence."""
import math
from rdkit import Chem
from rdkit.Chem import BRICS

def fragments(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return sorted(BRICS.BRICSDecompose(mol)) if mol is not None else []

class FragmentBank:
    def __init__(self, config, records=None, observations=None, comparisons=None):
        self.config = config
        self.records = records or {}
        self.observations = set(observations or [])
        self.comparisons = set(comparisons or [])

    def add(self, fragment, source, generation):
        mol = Chem.MolFromSmiles(fragment)
        if mol is None: return
        key = Chem.MolToSmiles(mol)
        if key not in self.records:
            self.records[key] = {"smiles":key, "attachments":[a.GetIsotope() for a in mol.GetAtoms() if a.GetAtomicNum()==0],
                "first_generation":generation, "sources":[], "association":{},
                "attempted":0, "valid":0, "duplicate":0, "evaluated":0, "improved":0,
                "reward_sum":0.0, "energy_delta_sum":{}, "removed_evaluated":0,
                "removed_reward_sum":0.0, "weight":1.0}
        r = self.records[key]
        if source not in r["sources"]: r["sources"].append(source)
        return key

    def ingest(self, smiles, generation):
        return [self.add(f, smiles, generation) for f in fragments(smiles)]

    def observe(self, smiles, energies, generation):
        keys = self.ingest(smiles, generation)
        for protein, energy in energies.items():
            token = protein+"\t"+smiles
            if token in self.observations: continue
            self.observations.add(token)
            for key in keys:
                r = self.records[key]["association"].setdefault(protein, {"n":0,"energy_sum":0.0,"best_energy":energy})
                r["n"] += 1; r["energy_sum"] += energy
                r["best_energy"] = min(r["best_energy"], energy)

    def status(self, event, status):
        for key in set(event.get("added", [])):
            if key in self.records: self.records[key][status] = self.records[key].get(status, 0) + 1

    def credit(self, parent, child, events, reward, old_energies, new_energies):
        token = parent+"\t"+child
        if token in self.comparisons or reward is None: return False
        self.comparisons.add(token)
        # Multi-step changes retain a combined reward in the event history only.
        if len(events) != 1 or events[0]["operator"] != "replace_fragment": return True
        event = events[0]
        for key in set(event.get("added", [])):
            r = self.records[key]; r["evaluated"] += 1
            r["improved"] += int(reward > 0); r["reward_sum"] += reward
            for p in new_energies:
                r["energy_delta_sum"][p] = r["energy_delta_sum"].get(p, 0)+new_energies[p]-old_energies[p]
        for key in set(event.get("removed", [])):
            if key in self.records:
                r = self.records[key]; r["removed_evaluated"] += 1; r["removed_reward_sum"] += reward
        return True

    def update_weights(self):
        for r in self.records.values():
            r["weight"] = math.exp(r["reward_sum"]/(r["evaluated"]+self.config.prior_strength)) if self.config.adaptive_fragments else 1.0

    def weights(self, keys=None):
        keys = sorted(self.records if keys is None else keys)
        if not keys: return {}
        total = sum(self.records[k]["weight"] for k in keys)
        e = self.config.exploration
        return {k:(1-e)*self.records[k]["weight"]/total+e/len(keys) for k in keys}

    def to_dict(self):
        return {"records":self.records, "observations":sorted(self.observations), "comparisons":sorted(self.comparisons)}
