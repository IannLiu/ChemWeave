"""Recipe-independent orbits of complete mapped reaction events.

Reactant graph automorphisms act on the atom/bond changes and product partition.
A stabilizer chain finds generators without enumerating spectator permutations.
Only the event orbit is enumerated; identical reactant factorials apply once.
"""

from collections import Counter, deque
from functools import lru_cache
from importlib.resources import files
import json
from math import factorial, prod

from rdkit import Chem

from .mapping import reaction_symmetry_key
from .rules import load_rules
from .sites import set_site_atoms, site_identity, site_label


RMG_PY_COMMIT = "032e9d20e953158e801f2de897dc5677d65bbccd"
CONVENTION = (
    "Fixed, independently labelled reactant occurrences. Enumerate the orbit of the complete "
    "mapped event under reactant graph automorphisms, retaining product species grouping and "
    "phase/site states. Deduplicate exact labelled events, not reactant-exchange orbits. "
    "n_e = unique_labeled_paths / product(factorial(identical reactant multiplicity)); apply this "
    "statistical correction exactly once. Atom-map numbers are identities, not symmetry colors. "
    "Constitutional paths in the supplied electronic graph; no resonance-ensemble, facial, "
    "rotational, lattice, orientation or spin multiplicity. No match_count or recipe input."
)


@lru_cache(maxsize=1)
def path_family_rules():
    """Example/validation templates only; never consulted by the factor engine."""
    return tuple(load_rules(files("chemweave").joinpath("data/path_families.json")))


def _automorphism_generators(labels, edges, atom_count, max_checks):
    """Find orbit transversals along an atom-point stabilizer chain.

    At level i, fix atoms < i and find one automorphism i -> j for each
    reachable j. Their union generates the full action on chemical atoms.
    Auxiliary component/bond nodes cannot introduce additional atom actions.
    """
    vocabulary = sorted(set(labels))
    graph = Chem.RWMol()
    for label in labels:
        atom = Chem.Atom(0)
        atom.SetNoImplicit(True)
        atom.SetIsotope(vocabulary.index(label) + 1)
        graph.AddAtom(atom)
    for a, b in edges:
        graph.AddBond(a, b, Chem.BondType.SINGLE)
    graph = graph.GetMol()
    graph.UpdatePropertyCache(strict=False)
    generators, checks, group_order = [], 0, 1
    marker = len(vocabulary) + 1
    for base in range(atom_count):
        ranks = Chem.CanonicalRankAtoms(graph, breakTies=False, includeChirality=False)
        query = Chem.Mol(graph)
        query.GetAtomWithIdx(base).SetIsotope(marker)
        canonical = Chem.MolToSmiles(query, canonical=True)
        query_order = json.loads(query.GetProp("_smilesAtomOutputOrder"))
        orbit_size = 1
        for candidate in range(base + 1, atom_count):
            if ranks[candidate] != ranks[base]:
                continue
            checks += 1
            if checks > max_checks:
                return None, checks, None
            target = Chem.Mol(graph)
            target.GetAtomWithIdx(candidate).SetIsotope(marker)
            if Chem.MolToSmiles(target, canonical=True) == canonical:
                target_order = json.loads(target.GetProp("_smilesAtomOutputOrder"))
                correspondence = dict(zip(query_order, target_order))
                generators.append(tuple(correspondence[i] for i in range(atom_count)))
                orbit_size += 1
        group_order *= orbit_size
        graph.GetAtomWithIdx(base).SetIsotope(marker)
        marker += 1
    return generators, checks, group_order


def calculate_path_factor(mapped_reactants, mapped_products, *, reactant_states=None,
                          product_states=None, site_atoms=None, max_paths=10000, max_checks=50000):
    """Return n_e and evidence from complete explicit-H mapped species lists.

    Each list entry is ONE species, including a disconnected associated ion pair.
    States are (phase, legacy_site_label) pairs; omitted phases infer surface for
    explicit dummy-site species, otherwise gas. Unspecified sites mean abstract X.
    No network, reaction family, recipe, parameter class or kinetic units required.
    Non-success results have status/reason/reason_code and never contain a factor.
    """
    if max_paths < 1 or max_checks < 1:
        raise ValueError("max_paths and max_checks must be positive")
    result = {"method": "reaction_path_degeneracy", "algorithm": "mapped_event_orbit_v1",
              "convention": CONVENTION, "source": "Automorphisms of the supplied mapped reactant graph"}
    params = Chem.SmilesParserParams()
    params.removeHs = False
    metadata = {int(k): v for k, v in (site_atoms or {}).items()}
    sides, states = [], []
    for strings, supplied_states in ((mapped_reactants, reactant_states), (mapped_products, product_states)):
        if isinstance(strings, str) or not strings or (supplied_states is not None and len(supplied_states) != len(strings)):
            return {**result, "status": "invalid_input", "reason_code": "species_partition",
                    "reason": "Provide nonempty lists of mapped species and one state per species"}
        mols, side_states = [], []
        for index, smiles in enumerate(strings):
            mol = Chem.MolFromSmiles(smiles, params)
            if mol is None or not mol.GetNumAtoms():
                return {**result, "status": "invalid_input", "reason_code": "invalid_structure", "reason": "Invalid mapped species SMILES"}
            mol = Chem.AddHs(mol)
            if any(not a.GetAtomMapNum() for a in mol.GetAtoms()):
                return {**result, "status": "missing_information", "reason_code": "incomplete_atom_mapping",
                        "reason": "Complete atom maps, including transferred and explicit spectator H, are required; generated channels supply them"}
            state = tuple(supplied_states[index]) if supplied_states is not None else (
                "surface" if any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()) else "gas", None)
            if len(state) != 2:
                return {**result, "status": "invalid_input", "reason_code": "species_state", "reason": "Each state must contain phase and optional site label"}
            if (state[0] == "surface" or state[1] is not None) and not any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
                return {**result, "status": "missing_information", "reason_code": "missing_site_nodes",
                        "reason": "Catalytic species need explicit mapped abstract site nodes, not only a species-level label"}
            try:
                set_site_atoms(mol, {a.GetAtomMapNum(): metadata.get(a.GetAtomMapNum(), {})
                                    for a in mol.GetAtoms() if a.GetAtomicNum() == 0})
                site_identity(mol, state[0])
            except ValueError as error:
                return {**result, "status": "invalid_input", "reason_code": "site_state", "reason": str(error)}
            mols.append(mol)
            side_states.append(state)
        sides.append(mols)
        states.append(side_states)
    inventories = [[a.GetAtomMapNum() for m in mols for a in m.GetAtoms()] for mols in sides]
    if any(len(nums) != len(set(nums)) for nums in inventories):
        return {**result, "status": "invalid_input", "reason_code": "duplicate_atom_maps", "reason": "Atom maps must be unique on each reaction side, including repeated reactant occurrences"}
    if set(inventories[0]) != set(inventories[1]):
        return {**result, "status": "invalid_input", "reason_code": "atom_inventory", "reason": "Both reaction sides must contain the same mapped atom inventory"}
    atoms = [{a.GetAtomMapNum(): a for m in mols for a in m.GetAtoms()} for mols in sides]
    if set(metadata) - {n for n, a in atoms[0].items() if a.GetAtomicNum() == 0}:
        return {**result, "status": "invalid_input", "reason_code": "site_map_reference", "reason": "Site metadata must reference mapped catalyst nodes"}
    if any((a.GetAtomicNum(), a.GetIsotope(), site_label(a)) !=
           (atoms[1][n].GetAtomicNum(), atoms[1][n].GetIsotope(), site_label(atoms[1][n])) for n, a in atoms[0].items()):
        return {**result, "status": "invalid_input", "reason_code": "atom_conservation", "reason": "Elements, isotopes and catalyst-node attributes must be conserved by the mapping"}
    if sum(a.GetFormalCharge() for a in atoms[0].values()) != sum(a.GetFormalCharge() for a in atoms[1].values()):
        return {**result, "status": "invalid_input", "reason_code": "charge_balance", "reason": "Mapped reaction does not conserve total charge"}
    if any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for side in atoms for a in side.values()) or any(
            b.GetStereo() != Chem.BondStereo.STEREONONE for mols in sides for m in mols for b in m.GetBonds()):
        return {**result, "status": "unsupported_algorithm", "reason_code": "stereochemical_action",
                "reason": "Assigned stereochemistry is present, but parity-preserving event orbits are not implemented"}

    numbers = sorted(atoms[0])
    indices = {n: i for i, n in enumerate(numbers)}
    atom_states, bonds, partitions = [], [], []
    for mols, side_states in zip(sides, states):
        atom_states.append({indices[a.GetAtomMapNum()]: (a.GetFormalCharge(), a.GetNumRadicalElectrons(), a.GetIsAromatic())
                            for m in mols for a in m.GetAtoms()})
        bonds.append({tuple(sorted((indices[b.GetBeginAtom().GetAtomMapNum()], indices[b.GetEndAtom().GetAtomMapNum()]))): str(b.GetBondType())
                      for m in mols for b in m.GetBonds()})
        partitions.append(tuple(sorted((json.dumps(state), tuple(sorted(indices[a.GetAtomMapNum()] for a in m.GetAtoms())))
                                       for m, state in zip(mols, side_states))))
    # The event contains changes AND the resulting species partition. The latter
    # distinguishes associated ions, phase transfers and disconnected adsorbates.
    changed_atoms = tuple(sorted((i, atom_states[0][i], value) for i, value in atom_states[1].items() if value != atom_states[0][i]))
    changed_bonds = tuple(sorted((a, b, bonds[0].get((a, b), ""), bonds[1].get((a, b), ""))
                                for a, b in bonds[0].keys() | bonds[1].keys() if bonds[0].get((a, b)) != bonds[1].get((a, b))))
    event = (changed_atoms, changed_bonds, partitions[1])
    if not changed_atoms and not changed_bonds and partitions[0] == partitions[1]:
        return {**result, "status": "invalid_input", "reason_code": "empty_event", "reason": "The mapped reaction contains no atom, bond or species-state change"}
    labels = [json.dumps(("atom", atoms[0][n].GetAtomicNum(), atoms[0][n].GetIsotope(),
                          atom_states[0][indices[n]], site_label(atoms[0][n]))) for n in numbers]
    edges = []
    for state, members in partitions[0]:
        node = len(labels)
        labels.append(json.dumps(("species", state)))
        edges.extend((node, i) for i in members)
    for (a, b), order in bonds[0].items():
        node = len(labels)
        labels.append(json.dumps(("bond", order)))
        edges.extend(((node, a), (node, b)))
    generators, checks, group_order = _automorphism_generators(labels, edges, len(numbers), max_checks)
    if generators is None:
        return {**result, "status": "incomplete", "reason_code": "automorphism_limit",
                "reason": "Automorphism search limit reached; no partial factor assigned", "counting": {"isomorphism_checks": checks}}
    paths, pending, actions = {event: tuple(range(len(numbers)))}, deque([event]), 0
    while pending:
        current = pending.popleft()
        for permutation in generators:
            moved = (tuple(sorted((permutation[i], before, after) for i, before, after in current[0])),
                     tuple(sorted((*sorted((permutation[a], permutation[b])), before, after) for a, b, before, after in current[1])),
                     tuple(sorted((state, tuple(sorted(permutation[i] for i in members))) for state, members in current[2])))
            actions += 1
            if moved not in paths:
                if len(paths) >= max_paths:
                    return {**result, "status": "incomplete", "reason_code": "event_limit",
                            "reason": "Event orbit limit reached; no partial factor assigned",
                            "counting": {"discovered_paths": len(paths), "max_paths": max_paths}}
                paths[moved] = tuple(permutation[i] for i in paths[current])
                pending.append(moved)
    identities = []
    for m, state in zip(sides[0], states[0]):
        identities.append(reaction_symmetry_key([m], [m], [state], [state]))
    multiplicities = sorted(Counter(identities).values(), reverse=True)
    divisor = prod(factorial(n) for n in multiplicities)
    count = len(paths)
    examples = []
    for permutation in list(paths.values())[:64]:
        products = []
        for original in sides[1]:
            m = Chem.Mol(original)
            for atom in m.GetAtoms():
                atom.SetAtomMapNum(numbers[permutation[indices[atom.GetAtomMapNum()]]])
            products.append(Chem.MolToSmiles(m))
        examples.append({"mapped_reactants": list(mapped_reactants), "mapped_products": products})
    key = reaction_symmetry_key(*sides, *states)
    def mapped_partition(partition):
        return [{"state": json.loads(state), "atoms": [numbers[i] for i in members]} for state, members in partition]
    return {**result, "status": "calculated", "factor": count / divisor, "n_e": count / divisor,
            "kinetic_channel": key,
            "event": {"atom_state_fields": ["charge", "radicals", "aromatic"],
                      "site_atoms": {str(n): json.loads(site_label(a)) for n, a in atoms[0].items() if site_label(a)},
                      "atom_changes": [{"atom": numbers[i], "before": list(before), "after": list(after)} for i, before, after in changed_atoms],
                      "bond_changes": [{"atoms": [numbers[a], numbers[b]], "before": before or None, "after": after or None} for a, b, before, after in changed_bonds],
                      "reactant_partition": mapped_partition(partitions[0]), "product_partition": mapped_partition(partitions[1])},
            "counting": {"unique_labeled_paths": count, "identical_reactant_divisor": divisor,
                         "identical_reactant_multiplicities": multiplicities, "correction": 1 / divisor,
                         "reactant_exchange_quotiented": False, "correction_applications": 1,
                         "automorphism_group_order": group_order if group_order <= 2**53 - 1 else str(group_order),
                         "event_stabilizer_order": group_order // count if group_order // count <= 2**53 - 1 else str(group_order // count),
                         "automorphism_generators": len(generators), "isomorphism_checks": checks,
                         "event_actions": actions, "duplicate_events_removed": actions + 1 - count,
                         "paths": examples, "paths_truncated": count > len(examples),
                         "representation": "supplied electronic graph; no resonance ensemble or unassigned stereo multiplicity",
                         "rmg_py_commit": RMG_PY_COMMIT, "rmg_runtime_cross_validated": False}}


def reaction_path_degeneracy(*, network, reaction, channel, center):
    """Network adapter; the pure calculator does not depend on rules or families."""
    result = calculate_path_factor(channel["mapped_reactant_smiles"], channel["mapped_product_smiles"],
        reactant_states=[(network.species[sid].phase, network.species[sid].site) for sid in channel["reactants"]],
        product_states=[(network.species[sid].phase, network.species[sid].site) for sid in channel["products"]],
        site_atoms=channel.get("site_atoms", {}))
    result["family"] = reaction.family  # provenance only, including "unspecified"
    return result
