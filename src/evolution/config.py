"""Validated, JSON-serializable settings. Energies are in kcal/mol."""
from dataclasses import dataclass, field, asdict
import math

def finite(value, name, minimum=None):
    if not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")

def probabilities(values, name):
    if not values:
        raise ValueError(f"{name} is empty")
    for key, value in values.items():
        finite(value, key, 0)
    if not math.isclose(sum(values.values()), 1.0, abs_tol=1e-8):
        raise ValueError(f"{name} must sum to 1")

@dataclass
class ProteinSpec:
    name: str
    weight: float = 1.0
    threshold: float | None = None
    offset: float = 0.0
    scale: float = 1.0

    def validate(self, off=False):
        if not self.name or any(c in self.name for c in "/\\") or self.name in (".", ".."):
            raise ValueError("Protein name must be a file stem")
        finite(self.weight, "weight", 0)
        finite(self.offset, "offset")
        finite(self.scale, "scale")
        if self.scale <= 0:
            raise ValueError("scale must be positive")
        if off and self.threshold is None:
            raise ValueError("Every off-target needs an explicit energy threshold")
        if self.threshold is not None:
            finite(self.threshold, "threshold")

@dataclass
class EvolutionConfig:
    targets: list[ProteinSpec]
    off_targets: list[ProteinSpec] = field(default_factory=list)
    penalty_strength: float = 1.0
    target_aggregation: str = "mean"
    population_size: int = 50
    generations: int = 5
    parent_count: int = 10
    elite_fraction: float = 0.1
    mutation_fraction: float = 0.65  # fraction of NON-elite slots
    level_probabilities: dict = field(default_factory=lambda: {"atom": 0.4, "fragment": 0.6})
    atom_probabilities: dict = field(default_factory=lambda: {
        "add_atom": 0.3, "remove_atom": 0.2, "substitute_atom": 0.3, "change_bond": 0.2})
    mutation_steps: int = 1
    adaptive_fragments: bool = True
    adaptive_operators: bool = False
    exploration: float = 0.2
    prior_strength: float = 5.0
    reward_scale: float = 1.0
    min_operator_observations: int = 5
    max_attempts_per_slot: int = 30
    max_similarity: float = 1.0
    max_mw: float = 500.0
    max_logp: float = 5.0
    min_nitrogen: int = 1
    seed: int = 42
    preparation_seed: int = 42
    exhaustiveness: int = 4
    num_workers: int = 3
    force_preparation: bool = False
    force_docking: bool = False

    def validate(self):
        for name in ("adaptive_fragments", "adaptive_operators", "force_preparation", "force_docking"):
            if not isinstance(getattr(self,name),bool): raise ValueError(f"{name} must be boolean")
        if not self.targets:
            raise ValueError("At least one target is required")
        for p in self.targets: p.validate()
        for p in self.off_targets: p.validate(off=True)
        names = self.proteins
        if len(set(names)) != len(names):
            raise ValueError("Protein lists must be unique and disjoint")
        if sum(p.weight for p in self.targets) <= 0:
            raise ValueError("At least one target weight must be positive")
        if self.target_aggregation not in ("mean", "worst"):
            raise ValueError("target_aggregation must be mean or worst")
        for name in ("population_size", "parent_count", "mutation_steps", "max_attempts_per_slot", "num_workers", "exhaustiveness", "min_operator_observations"):
            v = getattr(self, name)
            if isinstance(v, bool) or not isinstance(v, int) or v < 1: raise ValueError(name)
        for name in ("generations", "min_nitrogen", "seed", "preparation_seed"):
            v = getattr(self, name)
            if isinstance(v, bool) or not isinstance(v, int) or v < 0: raise ValueError(name)
        if self.preparation_seed > 2147483647: raise ValueError("preparation_seed must be in 0..2147483647")
        if self.seed > 2147483647: raise ValueError("seed must be in 1..2147483647")
        if self.seed == 0: raise ValueError("seed must be positive: Vina seed=0 is not reproducible")
        for name in ("elite_fraction", "mutation_fraction", "exploration", "max_similarity"):
            v = getattr(self, name); finite(v, name, 0)
            if v > 1: raise ValueError(name)
        if self.elite_fraction >= 1: raise ValueError("elite_fraction must be < 1")
        for name in ("penalty_strength", "max_mw", "prior_strength", "reward_scale"):
            finite(getattr(self, name), name, 0)
        if self.prior_strength <= 0 or self.reward_scale <= 0: raise ValueError("positive smoothing required")
        finite(self.max_logp, "max_logp")
        if set(self.level_probabilities) != {"atom", "fragment"}: raise ValueError("Unknown mutation level")
        if set(self.atom_probabilities) != {"add_atom", "remove_atom", "substitute_atom", "change_bond"}: raise ValueError("Unknown atom operator")
        probabilities(self.level_probabilities, "level_probabilities")
        probabilities(self.atom_probabilities, "atom_probabilities")
        return self

    @property
    def proteins(self):
        return [p.name for p in self.targets + self.off_targets]

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, values):
        values = dict(values)
        for key in ("targets", "off_targets"):
            values[key] = [ProteinSpec(**p) if isinstance(p, dict) else ProteinSpec(p) for p in values.get(key, [])]
        return cls(**values).validate()
