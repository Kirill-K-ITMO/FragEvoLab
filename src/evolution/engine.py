"""Evolutionary de novo orchestration. Single scalar fitness, no Pareto ranking."""
from pathlib import Path
import csv
import json
import logging
import random
import uuid
from collections import Counter
from src.ui.progress import evolution_progress
from src.core.calculation_cache import atomic_json
from src.evolution.config import EvolutionConfig, ProteinSpec
from src.evolution.validation import MolecularConstraints
from src.evolution.models import Candidate
from src.evolution.fitness import FitnessEvaluator, Score
from src.evolution.fragments import FragmentBank
from src.evolution.adaptation import Adaptation
from src.evolution.operators import Operators
from src.evolution.selection import select
from src.evolution import state
from src.evolution.docking_adapter import DockingAdapter

class EvolutionaryGenerator(DockingAdapter):
    def __init__(self, project_root=None):
        self.project_root = Path(project_root or Path.cwd()).resolve()
        self.output_dir = self.project_root/"output"/"evolutionary"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(__name__)
        self.evolution_history = []
        self.cache_stats = {"hits":0,"misses":0,"executed":0}
        self.exhaustiveness = 4
        self.force_preparation = self.force_docking = False
        self.preparation_seed = 42
        self.last_generation_stats = {}
        self.constraints = None

    @evolution_progress
    def run_evolution(self, initial_csv, docker, calculator, proteins=None,
                      top_n=10, num_generations=5, population_size=50, num_workers=3,
                      max_mw=500, max_logp=5.0, min_nitrogen=1, *,
                      exhaustiveness=4, force_preparation=False, force_docking=False,
                      preparation_seed=42, config=None, resume_from=None, stop_after_generation=None):
        """config is authoritative when supplied; generations excludes seed generation 0.

        proteins supplies target names when config is absent. EvolutionConfig
        supports both target and off-target settings.
        Resume requires identical settings and receptor/box/engine bytes, except that
        the generation limit may be extended. Checkpoints commit completed generations.
        """
        if config is None:
            config = EvolutionConfig(targets=[ProteinSpec(p) for p in (proteins or [])],
                parent_count=top_n, generations=num_generations,population_size=population_size,
                num_workers=num_workers,max_mw=max_mw,max_logp=max_logp,min_nitrogen=min_nitrogen,
                exhaustiveness=exhaustiveness,force_preparation=force_preparation,
                force_docking=force_docking,preparation_seed=preparation_seed)
        if isinstance(config,dict): config = EvolutionConfig.from_dict(config)
        self.config = config.validate()
        c = self.config
        self._progress.update(total=c.generations+1, completed=0)
        self.evaluator = FitnessEvaluator(c)
        self.constraints = MolecularConstraints(c.max_mw,c.max_logp,c.min_nitrogen)
        self.rng = random.Random(c.seed)
        from src.core.preflight import check_run
        check_run(self.project_root, config=c, seeds=initial_csv, resume=resume_from,
                  docker=docker, dependencies=False).require()
        self.fingerprint = state.fingerprint(c,docker)
        self.evolution_history, self.archive, self.attempts = [], {}, []
        self.cache_stats = {"hits":0,"misses":0,"executed":0}
        self.exhaustiveness, self.preparation_seed = c.exhaustiveness,c.preparation_seed
        self.force_preparation,self.force_docking = c.force_preparation,c.force_docking
        self.operators = None
        self.bank = FragmentBank(c)
        self.adaptation = Adaptation(c)
        self.run_id = uuid.uuid4().hex
        if resume_from:
            data = state.load(resume_from,self.fingerprint)
            self.run_id = data["run_id"]
            self.output_dir = Path(resume_from).resolve().parent
            self.rng.setstate(state.tuples(data["rng"]))
            self.archive = {s:Candidate(**v) for s,v in data["archive"].items()}
            self.population = [self.archive[s] for s in data["population"]]
            self.bank = FragmentBank(c,**data["bank"])
            self.adaptation = Adaptation(c,data["operators"])
            self.evolution_history = data["history"]
            self.cache_stats = data["cache_stats"]
            self.constraints.stats = Counter(data["validation"])
            self.attempts = data["attempts"]
            start = data["generation"]+1
            self._progress.update(completed=start)
            if c.generations < start-1: raise ValueError("Generation limit precedes checkpoint")
        else:
            self.output_dir = self.project_root/"output"/"evolutionary"/self.run_id
            self.output_dir.mkdir(parents=True,exist_ok=True)
            self.population = self._load_seeds(initial_csv)
            self._evaluate(self.population,docker,calculator,0)
            self._commit(0, {"seed":len(self.population)})
            start = 1
        self.operators = Operators(self.rng,c,self.bank,self.adaptation)
        self.stop_reason = "generation_limit"
        if stop_after_generation is not None and start > stop_after_generation:
            self.stop_reason = "paused"
            return self._report()
        for generation in range(start,c.generations+1):
            parents = select(self.population,c.parent_count,c.max_similarity)
            if not parents:
                self.stop_reason = "no_feasible_candidates"
                break
            self.operators.freeze()
            population, counts = self._offspring(parents,generation)
            if not population:
                self.stop_reason = "no_valid_offspring"
                break
            self.population = population
            self._evaluate(population,docker,calculator,generation)
            self._commit(generation,counts)
            if counts.get("new",0) == 0:
                self.stop_reason = "search_exhausted"
                break
            if stop_after_generation is not None and generation >= stop_after_generation:
                self.stop_reason = "paused"
                break
        if not select(self.population,1): self.stop_reason = "no_feasible_candidates"
        return self._report()

    def _load_seeds(self, path):
        found, seen = [], set()
        from src.preparation.inputs import read_ligand_rows
        for name,smiles in read_ligand_rows(path):
            canonical,reason = self.constraints.validate(smiles)
            if not canonical or canonical in seen:
                self.attempts.append({"generation":0,"input":smiles,"status":reason or "duplicate"})
                continue
            seen.add(canonical); found.append(Candidate(canonical,0))
        if not found: raise ValueError("No valid seeds satisfy the molecular constraints")
        # Explicit, seeded budget: do not silently exceed the configured population.
        if len(found)>self.config.population_size:
            found = self.rng.sample(sorted(found,key=lambda x:x.smiles),self.config.population_size)
        for cand in found: self.bank.ingest(cand.smiles,0)
        return found

    def _offspring(self, parents, generation):
        if getattr(self, "_progress", None): self._progress.update(description=f"Evolution G{generation}: mutation and assembly")
        c = self.config
        elite_n = min(len(parents), int(c.population_size*c.elite_fraction))
        # Validate elites through the same gate. Immutable parent records preserve origin.
        result = [p for p in parents[:elite_n] if self.constraints.validate(p.smiles)[0]]
        slots = c.population_size-len(result)
        mutation_n = int(slots*c.mutation_fraction+0.5)
        intended = ["mutation"]*mutation_n+["assembly"]*(slots-mutation_n)
        self.rng.shuffle(intended)
        counts = {"elite":len(result),"intended_mutations":mutation_n,"intended_assembly":slots-mutation_n,"new":0,"mutation":0,"assembly":0,"unfilled":0}
        seen = set(self.archive)
        fps = []
        from rdkit import Chem, DataStructs
        from rdkit.Chem import rdFingerprintGenerator
        maker = rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
        fps = [maker.GetFingerprint(Chem.MolFromSmiles(p.smiles)) for p in result]
        for kind in intended:
            accepted = False
            for _ in range(c.max_attempts_per_slot):
                parent = self.rng.choice(parents) if kind == "mutation" else None
                product,events = self.operators.mutate(parent.smiles,generation) if parent else self.operators.assemble()
                canonical,reason = self.constraints.validate(product) if product else (None,"no_product")
                if canonical and canonical in seen: reason = "duplicate"
                if canonical and not reason:
                    fp = maker.GetFingerprint(Chem.MolFromSmiles(canonical))
                    if any(DataStructs.TanimotoSimilarity(fp,other)>c.max_similarity for other in fps): reason = "similarity"
                status = reason or "valid"
                for event in events:
                    event["candidate_status"] = status
                    if status:
                        self.adaptation.record(event["operator"],status)
                        if "level" in event: self.adaptation.record(event["level"],status)
                        self.bank.status(event,status)
                self.attempts.append({"generation":generation,"origin":kind,"parent":parent.smiles if parent else None,
                                      "product":canonical,"status":status,"events":events})
                if reason: continue
                candidate = Candidate(canonical,generation,kind,parent.smiles if parent else None,events)
                result.append(candidate); seen.add(canonical); fps.append(fp)
                counts[kind] += 1; counts["new"] += 1; accepted = True
                break
            if not accepted: counts["unfilled"] += 1
        return result,counts

    def _evaluate(self, population, docker, calculator, generation):
        if getattr(self, "_progress", None): self._progress.update(description=f"Evolution G{generation}: preparation and docking")
        raw = self._dock_generation([p.smiles for p in population],docker,calculator,self.config.proteins,generation,self.config.num_workers)
        for i,candidate in enumerate(population,1):
            rows = raw.get(f"gen{generation}_mol{i}",{})
            candidate.energies = {p:r["affinity"] for p,r in rows.items()}
            score = self.evaluator.evaluate(candidate.energies)
            candidate.score = score.to_dict()
            self.bank.observe(candidate.smiles,candidate.energies,generation)
            if candidate.smiles not in self.archive and candidate.parent and score.complete:
                parent = self.archive[candidate.parent]
                reward = self.evaluator.improvement(Score(**parent.score),score)
                fresh = self.bank.credit(parent.smiles,candidate.smiles,candidate.events,reward,parent.energies,candidate.energies)
                if fresh and reward is not None:
                    # Chain has one combined operator category, never per-step credits.
                    op = candidate.events[0]["operator"] if len(candidate.events)==1 else "mutation_chain"
                    self.adaptation.record(op,"reward",reward)
                    if len(candidate.events)==1: self.adaptation.record(candidate.events[0]["level"],"reward",reward)
                    for event in candidate.events: event["chain_reward"] = reward
            if candidate.smiles not in self.archive and not score.complete:
                for event in candidate.events:
                    self.adaptation.record(event["operator"], "evaluation_failed")
                    self.bank.status(event, "evaluation_failed")
            self.archive[candidate.smiles] = candidate
        self.bank.update_weights()

    def _commit(self,generation,counts):
        self.evolution_history.append({"generation":generation,"molecules":[p.smiles for p in self.population],
            "results":{p.smiles:{"energies":dict(p.energies),"score":dict(p.score)} for p in self.population},
            "counts":counts,"cache_hits":self.last_generation_stats.get("hits",0),
            "new_dockings":self.last_generation_stats.get("executed",0),
            "exploration_count":counts.get("assembly",0),"exploitation_count":counts.get("mutation",0)+counts.get("elite",0),
            "fragment_weights":{k:v["weight"] for k,v in self.bank.records.items()},
            "operator_statistics":json.loads(json.dumps(self.adaptation.stats)),
            "used_level_probabilities":dict(self.operators.levels) if self.operators else {},
            "used_atom_probabilities":dict(self.operators.atoms) if self.operators else {}})
        payload = {"fingerprint":self.fingerprint,"run_id":self.run_id,"generation":generation,
            "config":self.config.to_dict(),"rng":self.rng.getstate(),
            "population":[p.smiles for p in self.population],"archive":{s:c.to_dict() for s,c in self.archive.items()},
            "bank":self.bank.to_dict(),"operators":self.adaptation.stats,"history":self.evolution_history,
            "attempts":self.attempts,"cache_stats":self.cache_stats,"validation":dict(self.constraints.stats)}
        state.save(self.output_dir/"checkpoint.json",payload)
        if getattr(self, "_progress", None): self._progress.update(completed=generation+1, description=f"Evolution: G{generation} saved")

    def _report(self):
        ranked = select(list(self.archive.values()),self.config.parent_count,self.config.max_similarity)
        report = {"run_id":self.run_id,"stop_reason":self.stop_reason,"config":self.config.to_dict(),
            "total_generations":len(self.evolution_history),"unique_molecules":len(self.archive),
            "total_dockings":self.cache_stats["executed"],"cache_hits":self.cache_stats["hits"],
            "cache_misses":self.cache_stats["misses"],"generations":self.evolution_history,
            "top_molecules":[c.to_dict() for c in ranked],
            "candidates":[c.to_dict() for c in self.archive.values()],
            "attempts":self.attempts,"operators":self.adaptation.stats,
            "fragments":self.bank.to_dict(),"validation":dict(self.constraints.stats),
            "checkpoint":str(self.output_dir/"checkpoint.json"),
            "json_report":str(self.output_dir/"report.json"),"html_report":str(self.output_dir/"report.html")}
        atomic_json(Path(report["json_report"]),report)
        with (self.output_dir/"candidates.csv").open("w",newline="",encoding="utf-8") as f:
            fields = ["smiles","generation","origin","parent","fitness","feasible"]+self.config.proteins
            writer = csv.DictWriter(f,fields); writer.writeheader()
            for cand in self.archive.values():
                writer.writerow({"smiles":cand.smiles,"generation":cand.generation,"origin":cand.origin,
                    "parent":cand.parent,"fitness":cand.score.get("fitness"),"feasible":cand.score.get("feasible"),**cand.energies})
        from src.reporting.evolution import write_report
        write_report(report, Path(report["html_report"]))
        return report

