"""One evaluator for selection and matched parent-child rewards."""
from dataclasses import dataclass, asdict
import math

@dataclass
class Score:
    complete: bool
    feasible: bool
    fitness: float | None = None
    target_utility: float | None = None
    off_penalty: float | None = None
    reason: str = ""

    def to_dict(self): return asdict(self)

class FitnessEvaluator:
    def __init__(self, config):
        self.config = config.validate()

    def evaluate(self, energies):
        c = self.config
        if any(p not in energies or not isinstance(energies[p], (int,float)) or not math.isfinite(energies[p]) for p in c.proteins):
            return Score(False, False, reason="incomplete_docking")
        feasible = all(p.threshold is None or energies[p.name] <= p.threshold for p in c.targets)
        terms = [(-(energies[p.name]-p.offset)/p.scale, p.weight) for p in c.targets if p.weight > 0]
        utility = min(v for v,w in terms) if c.target_aggregation == "worst" else sum(v*w for v,w in terms)/sum(w for v,w in terms)
        penalty = max((p.weight * max(0.0, (p.threshold-energies[p.name])/p.scale) for p in c.off_targets), default=0.0)
        return Score(True, feasible, utility-c.penalty_strength*penalty, utility, penalty,
                     "" if feasible else "target_threshold")

    def improvement(self, parent, child):
        if not parent.complete or not child.complete: return None
        # Feasibility and fitness must agree with selection. Infeasible parents cannot
        # reproduce in the current policy, but this also protects imported comparisons.
        if parent.feasible != child.feasible: return 1.0 if child.feasible else -1.0
        return math.tanh((child.fitness-parent.fitness)/self.config.reward_scale)
