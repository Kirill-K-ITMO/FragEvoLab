"""Bounded, seeded graph edits and BRICS attachment-aware assembly."""
from rdkit import Chem
from rdkit.Chem import BRICS
from .adaptation import choose
from .chemistry_rules import enumerate_mutations, apply_mutation, enumerate_brics_joins, apply_brics_join

def canonical(mol):
    try:
        Chem.SanitizeMol(mol)
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    except (ValueError, RuntimeError):
        return None

def attachment(mol):
    atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum()==0]
    return atoms

def join(left, right):
    """Enumerate compatible single attachments, retaining each input once."""
    results = []
    for plan in enumerate_brics_joins(left, right).plans:
        result = apply_brics_join(left, right, plan)
        if result.ok: results.append(result.smiles)
    return sorted(set(results))

class Operators:
    def __init__(self, rng, config, bank, adaptation):
        self.rng, self.config, self.bank, self.adaptation = rng, config, bank, adaptation
        self.freeze()

    def freeze(self):
        self.levels = self.adaptation.weights(self.config.level_probabilities)
        self.atoms = self.adaptation.weights(self.config.atom_probabilities)
        self.fragment_weights = self.bank.weights()

    def mutate(self, smiles, generation):
        current, events = smiles, []
        for _ in range(self.config.mutation_steps):
            level = choose(self.rng, self.levels)
            op = "replace_fragment" if level == "fragment" else choose(self.rng, self.atoms)
            event = {"operator":op, "level":level, "before":current, "added":[], "removed":[]}
            self.adaptation.record(op, "attempted")
            self.adaptation.record(level, "attempted")
            if level == "fragment":
                product = self.replace(current, event, generation)
            else: product = self.atom_edit(current, op, event)
            event["after"] = product
            event["status"] = "product" if product and product != current else "no_product"
            events.append(event)
            if not product or product == current: return None, events
            current = product
        return current, events

    def atom_edit(self, smiles, op, event):
        mol = Chem.MolFromSmiles(smiles)
        proposed = enumerate_mutations(mol, op)
        event["excluded_sites"] = proposed.rejected
        if not proposed.plans:
            event["reason"] = "no_eligible_site"
            return None
        # First sample a site, then its allowed edit. More element options do not
        # give a site greater selection probability.
        sites = sorted(set(plan.site for plan in proposed.plans))
        site = self.rng.choice(sites)
        plan = self.rng.choice([p for p in proposed.plans if p.site == site])
        result = apply_mutation(mol, plan)
        event.update(site=list(site), atomic_number=plan.atomic_number, bond_order=plan.bond_order)
        if not result.ok: event["reason"] = result.reason
        return result.smiles if result.ok else None

    def replace(self, smiles, event, generation):
        mol = Chem.MolFromSmiles(smiles)
        cuts = sorted(BRICS.FindBRICSBonds(mol))
        if not cuts:
            event["reason"] = "no_brics_cut"
            return None
        cut = self.rng.choice(cuts)
        parts = Chem.GetMolFrags(BRICS.BreakBRICSBonds(mol, bonds=[cut]), asMols=True)
        if len(parts)!=2: return None
        parts = sorted(parts,key=lambda m:(-m.GetNumHeavyAtoms(),Chem.MolToSmiles(m)))
        retained, removed = parts
        removed_s = Chem.MolToSmiles(removed)
        self.bank.add(removed_s, smiles, generation)
        compatible = {}
        products = {}
        for key, weight in self.fragment_weights.items():
            if len(self.bank.records[key]["attachments"]) != 1 or key == removed_s: continue
            candidates = [s for s in join(retained,Chem.MolFromSmiles(key)) if "*" not in s and s != smiles]
            if candidates: compatible[key] = weight; products[key] = candidates
        if not compatible:
            event["reason"] = "no_compatible_fragment"
            return None
        key = choose(self.rng, compatible)
        event.update(added=[key],removed=[removed_s],site=list(cut[0]),retained=Chem.MolToSmiles(retained))
        self.bank.status(event,"attempted")
        return self.rng.choice(products[key])

    def assemble(self):
        keys = [k for k in self.fragment_weights if self.bank.records[k]["attachments"]]
        if not keys: return None, [{"operator":"assembly","added":[],"removed":[],"status":"empty_bank"}]
        # Start from a weighted fragment and join at most four partners; no unbounded BRICSBuild.
        key = choose(self.rng,{k:self.fragment_weights[k] for k in keys})
        current = key; used = [key]
        for _ in range(4):
            mol = Chem.MolFromSmiles(current)
            if not attachment(mol): break
            partner = choose(self.rng,{k:self.fragment_weights[k] for k in keys})
            options = join(mol,Chem.MolFromSmiles(partner))
            if not options: continue
            current = self.rng.choice(options); used.append(partner)
        valid = not attachment(Chem.MolFromSmiles(current))
        event = {"operator":"assembly","added":used,"removed":[],"status":"product" if valid else "unclosed_attachments"}
        self.adaptation.record("assembly","attempted"); self.bank.status(event,"attempted")
        return (current if valid else None), [event]
