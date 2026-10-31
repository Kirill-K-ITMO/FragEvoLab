"""Evolution CLI: python -m src.evolution.cli --help"""
import argparse
import json
from pathlib import Path
from src.evolution.config import EvolutionConfig
from src.evolution.engine import EvolutionaryGenerator
from src.docking.vina import AutoDockVina
from src.docking.affinity import AffinityCalculator

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--project",type=Path,default=Path.cwd())
    p.add_argument("--config",type=Path,required=True)
    p.add_argument("--seeds",type=Path)
    p.add_argument("--resume",type=Path)
    p.add_argument("--vina",default=None)
    p.add_argument("--cpu",type=int,default=None,help="Maximum Vina threads per job")
    p.add_argument("--cpu-budget",type=int,default=None,help="Total Vina thread budget for this application/project")
    p.add_argument("--timeout",type=float,default=3600)
    args = p.parse_args()
    if not args.resume and not args.seeds: p.error("--seeds is required for a new run")
    config = EvolutionConfig.from_dict(json.loads(args.config.read_text(encoding="utf-8-sig")))
    docker = AutoDockVina(args.project,vina_executable=args.vina,cpu=args.cpu,cpu_budget=args.cpu_budget,seed=config.seed,timeout=args.timeout)
    from src.core.preflight import check_run, display
    readiness = check_run(args.project, config=config, seeds=args.seeds, resume=args.resume, docker=docker)
    display(readiness)
    if not readiness.ready:
        raise SystemExit(1)
    report = EvolutionaryGenerator(args.project).run_evolution(args.seeds,docker,AffinityCalculator(),config=config,resume_from=args.resume)
    print(json.dumps({k:report[k] for k in ("stop_reason","json_report","checkpoint")},indent=2))

if __name__ == "__main__": main()
