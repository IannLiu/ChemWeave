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
    max_depth: int | None = None


class Generator:
    def __init__(self, rules, config=None):
        self.rules = list(rules)
        self.config = config or GenerationConfig()
        if len({rule.id for rule in self.rules}) != len(self.rules):
            raise ValueError("Rule IDs must be unique")
        if min(self.config.max_carbon, self.config.max_species, self.config.max_reactions) < 1:
            raise ValueError("Generation limits must be positive")
        if self.config.max_depth is not None and self.config.max_depth < 0:
            raise ValueError("max_depth must be nonnegative")

    def generate(self, seeds, *, targets=(), progress=None, stop_requested=None):
        """Grow from seeds until exhausted, a limit, or all target species appear.

        Each species is matched against each template slot once. Only newly
        available reactant tuples are executed. Ordered tuples are retained:
        template slots need not be interchangeable. Seeds have depth zero;
        products have depth one greater than the deepest reactant. FIFO order
        makes the first discovery of each species its minimum reachable depth.
        Targets are Species objects; identity includes phase and site.
        Optional callbacks report completed applications and request a stop
        between candidates, preserving complete accepted reaction records.
        """
        network = Network(
            rules=[rule.to_dict() for rule in self.rules], config=asdict(self.config),
            stats={"applications": 0, "rejected": {}, "rdkit_version": rdBase.rdkitVersion},
        )
        identities = {}
        pending = deque()
        target_records = {}
        for target in targets:
            target_records.setdefault(target.key, {
                "smiles": target.smiles, "phase": target.phase, "site": target.site,
                "name": target.name, "species_id": None,
            })
        network.targets = list(target_records.values())
        remaining_targets = set(target_records)
        for species in seeds:
            if species.carbon_count > self.config.max_carbon:
                raise ValueError(f"Seed exceeds max_carbon: {species.smiles}")
            if species.key not in identities:
                if len(network.species) >= self.config.max_species:
                    raise ValueError("Seeds exceed max_species")
                sid = f"S{len(network.species) + 1:05d}"
                identities[species.key] = sid
                network.species[sid] = species
                network.species_depths[sid] = 0
                pending.append(sid)
                network.seeds.append(sid)
                if species.key in remaining_targets:
                    target_records[species.key]["species_id"] = sid
                    remaining_targets.remove(species.key)
            network.seed_inputs.append({
                "species_id": identities[species.key], "mapped_smiles": species.mapped_smiles,
                "input_maps": species.input_maps, "name": species.name,
            })
        if progress is not None:
            progress(network)
        if target_records and not remaining_targets:
            network.termination = "targets_reached"
            return network

        pools = [[[] for _ in range(rule.arity)] for rule in self.rules]
        applications = set()
        reaction_index = {}
        channel_index = {}
        rejected = Counter()
        while pending:
            if stop_requested is not None and stop_requested():
                network.termination = "user_stopped"
                break
            sid = pending.popleft()
            depth = network.species_depths[sid]
            # All remaining species are at least this deep. Boundary species
            # are retained, but cannot participate in another allowed step.
            if self.config.max_depth is not None and depth >= self.config.max_depth:
                network.termination = "max_depth"
                break
            species = network.species[sid]
            for ri, rule in enumerate(self.rules):
                for slot in range(rule.arity):
                    if not rule.matches(species, slot):
                        continue
                    pools[ri][slot].append(sid)
                    choices = [([sid] if j == slot else pool) for j, pool in enumerate(pools[ri])]
                    for reactant_ids in product(*choices):
                        if stop_requested is not None and stop_requested():
                            network.termination = "user_stopped"
                            network.stats["rejected"] = dict(rejected)
                            return network
                        application = (rule.id, reactant_ids)
                        if application in applications:
                            continue
                        applications.add(application)
                        network.stats["applications"] += 1
                        reactants = tuple(network.species[rid] for rid in reactant_ids)
                        for products, provenance in rule.apply(reactants, rejected):
                            if stop_requested is not None and stop_requested():
                                network.termination = "user_stopped"
                                network.stats["rejected"] = dict(rejected)
                                return network
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
                            reaction = reaction_index.get(key)
                            if reaction is None:
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
                                    network.species_depths[pid] = depth + 1
                                    pending.append(pid)
                                    if identity in remaining_targets:
                                        target_records[identity]["species_id"] = pid
                                        remaining_targets.remove(identity)
                                reaction = Reaction(
                                    id=f"R{len(network.reactions) + 1:05d}",
                                    reactants=lhs, products=rhs, rule_id=rule.id,
                                    family=rule.family, kind=rule.kind,
                                    depth=depth + 1,
                                )
                                network.reactions.append(reaction)
                                reaction_index[key] = reaction
                            channel_key = (reaction.id, provenance["symmetry_key"])
                            if channel_key in channel_index:
                                channel_index[channel_key]["match_count"] += 1
                            else:
                                channel = {
                                    "reactants": list(reactant_ids),
                                    "products": product_ids, **provenance,
                                    "match_count": 1,
                                }
                                reaction.channels.append(channel)
                                channel_index[channel_key] = channel
                            if target_records and not remaining_targets:
                                network.termination = "targets_reached"
                                network.stats["rejected"] = dict(rejected)
                                return network
                        if progress is not None:
                            network.stats["rejected"] = dict(rejected)
                            progress(network)
        network.stats["rejected"] = dict(rejected)
        return network
