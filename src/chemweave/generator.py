"""Incremental, rule-based network enumeration."""

from collections import Counter, deque
from dataclasses import asdict, dataclass
from itertools import product

from rdkit import rdBase

from .models import Network, Reaction, Species


@dataclass(frozen=True, slots=True)
class GenerationConfig:
    max_carbon: int = 8
    max_species: int = 1000
    max_reactions: int = 5000


class Generator:
    def __init__(self, rules, config=None):
        self.rules = list(rules)
        self.config = config or GenerationConfig()
        if len({rule.id for rule in self.rules}) != len(self.rules):
            raise ValueError("Rule IDs must be unique")
        if min(asdict(self.config).values()) < 1:
            raise ValueError("Generation limits must be positive")

    def generate(self, seeds):
        """Enumerate the closure within max_carbon, or stop at a size budget.

        Each species is matched against each template slot once. Only newly
        available reactant tuples are executed. Ordered tuples are retained:
        template slots need not be interchangeable.
        """
        network = Network(
            rules=[rule.to_dict() for rule in self.rules], config=asdict(self.config),
            stats={"applications": 0, "rejected": {}, "rdkit_version": rdBase.rdkitVersion},
        )
        identities = {}
        pending = deque()
        for species in seeds:
            if species.carbon_count > self.config.max_carbon:
                raise ValueError(f"Seed exceeds max_carbon: {species.smiles}")
            if species.key not in identities:
                if len(network.species) >= self.config.max_species:
                    raise ValueError("Seeds exceed max_species")
                sid = f"S{len(network.species) + 1:05d}"
                identities[species.key] = sid
                network.species[sid] = species
                pending.append(sid)
                network.seeds.append(sid)

        pools = [[[] for _ in range(rule.arity)] for rule in self.rules]
        applications = set()
        reaction_index = {}
        rejected = Counter()
        while pending:
            sid = pending.popleft()
            species = network.species[sid]
            for ri, rule in enumerate(self.rules):
                for slot in range(rule.arity):
                    if not rule.matches(species, slot):
                        continue
                    pools[ri][slot].append(sid)
                    choices = [([sid] if j == slot else pool) for j, pool in enumerate(pools[ri])]
                    for reactant_ids in product(*choices):
                        application = (rule.id, reactant_ids)
                        if application in applications:
                            continue
                        applications.add(application)
                        network.stats["applications"] += 1
                        reactants = tuple(network.species[rid] for rid in reactant_ids)
                        for products, provenance in rule.apply(reactants, rejected):
                            if any(s.carbon_count > self.config.max_carbon for s in products):
                                rejected["max_carbon"] += 1
                                continue
                            lhs = tuple(sorted(Counter(reactant_ids).items()))
                            # Assign provisional IDs so budget checks do not leave
                            # orphan species behind when a reaction cannot be added.
                            new_species = {}
                            product_ids = []
                            for s in products:
                                pid = identities.get(s.key) or new_species.get(s.key, (None,))[0]
                                if pid is None:
                                    pid = f"S{len(network.species) + len(new_species) + 1:05d}"
                                    new_species[s.key] = (pid, s)
                                product_ids.append(pid)
                            rhs = tuple(sorted(Counter(product_ids).items()))
                            if lhs == rhs:
                                rejected["identity"] += 1
                                continue
                            key = (rule.id, lhs, rhs)
                            match = {"reactants": list(reactant_ids),
                                     "products": product_ids, **provenance}
                            if key in reaction_index:
                                reaction_index[key].matches.append(match)
                                continue
                            if len(network.species) + len(new_species) > self.config.max_species:
                                network.termination = "max_species"
                            elif len(network.reactions) >= self.config.max_reactions:
                                network.termination = "max_reactions"
                            if network.termination != "exhausted":
                                network.stats["rejected"] = dict(rejected)
                                return network
                            for identity, (pid, s) in new_species.items():
                                identities[identity] = pid
                                network.species[pid] = s
                                pending.append(pid)
                            reaction = Reaction(
                                id=f"R{len(network.reactions) + 1:05d}",
                                reactants=lhs, products=rhs, rule_id=rule.id,
                                family=rule.family, kind=rule.kind, matches=[match],
                            )
                            network.reactions.append(reaction)
                            reaction_index[key] = reaction
        network.stats["rejected"] = dict(rejected)
        return network
