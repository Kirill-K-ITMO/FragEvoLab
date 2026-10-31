"""Smoothed bounded weights; disabled operators never become enabled."""
import math

class Adaptation:
    def __init__(self, config, stats=None):
        self.config = config
        self.stats = stats or {}

    def record(self, operator, status, reward=None):
        s = self.stats.setdefault(operator, {"attempted":0,"valid":0,"duplicate":0,"evaluated":0,"improved":0,"reward_sum":0.0})
        if status != "reward": s[status] = s.get(status, 0) + 1
        if reward is not None:
            s["evaluated"] += 1
            s["improved"] += int(reward > 0)
            s["reward_sum"] += reward

    def weights(self, base):
        c = self.config
        active = {k:v for k,v in base.items() if v > 0}
        if not c.adaptive_operators: return active
        adjusted = {}
        for k,v in active.items():
            s = self.stats.get(k, {})
            n = s.get("evaluated",0)
            effect = s.get("reward_sum",0)/(n+c.prior_strength) if n >= c.min_operator_observations else 0
            adjusted[k] = v * math.exp(effect)
        total = sum(adjusted.values())
        return {k:(1-c.exploration)*v/total + c.exploration*active[k]/sum(active.values()) for k,v in adjusted.items()}

def choose(rng, weights):
    keys = sorted(k for k,v in weights.items() if v > 0)
    if not keys: raise ValueError("No enabled choices")
    return rng.choices(keys, weights=[weights[k] for k in keys], k=1)[0]
