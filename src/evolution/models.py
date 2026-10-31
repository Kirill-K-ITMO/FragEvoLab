from dataclasses import dataclass, field, asdict

@dataclass
class Candidate:
    smiles: str
    generation: int
    origin: str = "seed"
    parent: str | None = None
    events: list = field(default_factory=list)
    energies: dict = field(default_factory=dict)
    score: dict = field(default_factory=dict)

    def to_dict(self): return asdict(self)
