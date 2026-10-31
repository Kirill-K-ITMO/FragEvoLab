"""Re-score a completed checkpoint without docking or reusing adaptive rewards."""
import argparse
import json
from pathlib import Path
from src.core.calculation_cache import atomic_json
from .config import EvolutionConfig
from .fitness import FitnessEvaluator
from .models import Candidate
from .selection import select

def reweight(checkpoint, config, output):
    checkpoint, output = Path(checkpoint).resolve(), Path(output).resolve()
    if checkpoint == output: raise ValueError("Output must not overwrite the source checkpoint")
    data = json.loads(checkpoint.read_text(encoding="utf-8-sig"))
    evaluator = FitnessEvaluator(config)
    candidates = []
    for row in data["archive"].values():
        candidate = Candidate(**row)
        candidate.score = evaluator.evaluate(candidate.energies).to_dict()
        candidates.append(candidate)
    report = {"source_checkpoint":str(checkpoint),"config":config.to_dict(),
              "purpose":"offline_reweighting_not_a_resumable_run",
              "candidates":[c.to_dict() for c in candidates],
              "top_molecules":[c.to_dict() for c in select(candidates,config.parent_count,config.max_similarity)],
              "incomplete":sum(not c.score["complete"] for c in candidates)}
    atomic_json(output,report)
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--checkpoint",type=Path,required=True)
    p.add_argument("--config",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    result=reweight(args.checkpoint,EvolutionConfig.from_dict(json.loads(args.config.read_text(encoding="utf-8-sig"))),args.output)
    print(f"Saved {len(result['candidates'])} candidates; incomplete: {result['incomplete']}")

if __name__ == "__main__": main()
