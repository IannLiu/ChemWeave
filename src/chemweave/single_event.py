"""Optional, post-generation single-event annotations; no rate fitting or solver.

The built-in calculator covers forward gas-phase H + methane H/D abstraction.
It uses proper rotational symmetry numbers of a tetrahedral reactant and an
axial C--H--H transition state, never template matches or graph automorphisms.
"""

from collections import Counter
from copy import deepcopy
from math import isclose, isfinite, prod

from rdkit import Chem


METHANE_CONVENTION = (
    "Forward gas-phase H/D atom + methane(H/D): tetrahedral methane; axial "
    "C--H/D--H/D saddle with three spectator ligands. Proper rotations only: "
    "sigma_R = 12 (4 identical), 3 (3+1), 2 (2+2); sigma_TS = 3 when all "
    "spectators are identical, otherwise 1. Atomic reactant sigma = 1. "
    "n_SE = product(sigma_R)/sigma_TS; no product symmetry, spin, tunneling, "
    "or nuclear-spin weights included. Intrinsic parameters exclude n_SE."
)


def reaction_center(mapped_reaction_smiles):
    """Extract bond and atom-state changes using reaction-local atom maps."""
    params = Chem.SmilesParserParams()
    params.removeHs = False
    sides = [Chem.MolFromSmiles(s, params) for s in mapped_reaction_smiles.split(">>")]
    if len(sides) != 2 or any(m is None for m in sides):
        raise ValueError("A complete mapped reaction SMILES is required")
    inventories, bonds, states = [], [], []
    for mol in sides:
        atoms = {a.GetAtomMapNum(): a for a in mol.GetAtoms()}
        if 0 in atoms or len(atoms) != mol.GetNumAtoms():
            raise ValueError("Reaction maps must be nonzero and unique on each side")
        inventories.append({n: (a.GetAtomicNum(), a.GetIsotope()) for n, a in atoms.items()})
        states.append({n: {"charge": a.GetFormalCharge(), "radicals": a.GetNumRadicalElectrons(),
                           "stereo": a.GetProp("_CIPCode") if a.HasProp("_CIPCode") else ""}
                       for n, a in atoms.items()})
        bonds.append({tuple(sorted((b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum()))):
                      {"order": b.GetBondTypeAsDouble(), "stereo": str(b.GetStereo())}
                      for b in mol.GetBonds()})
    if inventories[0] != inventories[1]:
        raise ValueError("Reaction maps must conserve elements and isotopes bijectively")
    changes, atom_changes, center = [], [], set()
    for pair in sorted(bonds[0].keys() | bonds[1].keys()):
        before, after = bonds[0].get(pair), bonds[1].get(pair)
        if before != after:
            changes.append({"atoms": list(pair), "before": before, "after": after})
            center.update(pair)
    for number in sorted(states[0]):
        if states[0][number] != states[1][number]:
            atom_changes.append({"atom": number, "before": states[0][number], "after": states[1][number]})
            center.add(number)
    return {"atom_maps": sorted(center), "bond_changes": changes, "atom_changes": atom_changes}


def _methane_abstraction(channel, center):
    """Recognize the complete elementary transformation, independent of rule ID."""
    params = Chem.SmilesParserParams()
    params.removeHs = False
    left = Chem.MolFromSmiles(channel["mapped_reaction_smiles"].split(">>")[0], params)
    atoms = {a.GetAtomMapNum(): a for a in left.GetAtoms()}
    carbons = [a for a in left.GetAtoms() if a.GetAtomicNum() == 6]
    if left.GetNumAtoms() != 6 or len(carbons) != 1:
        return None
    carbon = carbons[0]
    hydrogens = list(carbon.GetNeighbors())
    abstractors = [a for a in left.GetAtoms() if a.GetDegree() == 0 and a.GetAtomicNum() == 1]
    if (carbon.GetDegree() != 4 or len(abstractors) != 1
            or any(a.GetAtomicNum() != 1 or a.GetIsotope() not in (0, 1, 2) for a in hydrogens)
            or any(a.GetFormalCharge() for a in left.GetAtoms())
            or any(a.GetNumRadicalElectrons() != (1 if a.GetIdx() == abstractors[0].GetIdx() else 0) for a in left.GetAtoms())
            or any(b.GetBondTypeAsDouble() != 1 for b in left.GetBonds())):
        return None
    abstractor = abstractors[0]
    if abstractor.GetIsotope() not in (0, 1, 2):
        return None
    # RDKit groups isotope 0 and explicit 1H into separate mapping channels.
    # Mixing them would assign the full physical factor to both channels.
    if {a.GetIsotope() for a in hydrogens} >= {0, 1}:
        return None
    broken = [c for c in center["bond_changes"] if c["before"] == {"order": 1.0, "stereo": "STEREONONE"} and c["after"] is None]
    formed = [c for c in center["bond_changes"] if c["before"] is None and c["after"] == {"order": 1.0, "stereo": "STEREONONE"}]
    if len(center["bond_changes"]) != 2 or len(broken) != 1 or len(formed) != 1:
        return None
    carbon_map, abstractor_map = carbon.GetAtomMapNum(), abstractor.GetAtomMapNum()
    transfer = set(broken[0]["atoms"]) - {carbon_map}
    if carbon_map not in broken[0]["atoms"] or len(transfer) != 1:
        return None
    transfer = transfer.pop()
    if set(formed[0]["atoms"]) != {transfer, abstractor_map}:
        return None
    # No extra charge/radical/stereo changes may hide behind the same bond edits.
    expected = {carbon_map: (0, 1), abstractor_map: (1, 0)}
    if (len(center["atom_changes"]) != 2 or any(
            c["atom"] not in expected
            or (c["before"]["radicals"], c["after"]["radicals"]) != expected[c["atom"]]
            or c["before"]["charge"] != 0 or c["after"]["charge"] != 0
            or c["before"]["stereo"] or c["after"]["stereo"]
            for c in center["atom_changes"])):
        return None
    isotope = lambda a: a.GetIsotope() or 1
    counts = Counter(isotope(a) for a in hydrogens)
    sigma_r = {(4,): 12, (1, 3): 3, (2, 2): 2}[tuple(sorted(counts.values()))]
    spectators = sorted(isotope(a) for a in hydrogens if a.GetAtomMapNum() != transfer)
    sigma_ts = 3 if len(set(spectators)) == 1 else 1
    subclass = (f"methane_abstraction:C{carbon.GetIsotope() or 12}:"
                f"transfer{isotope(atoms[transfer])}:atom{isotope(abstractor)}:"
                f"spectators{''.join(map(str, spectators))}")
    return {"status": "calculated", "kinetic_subclass": subclass, "parameter_class": subclass,
            "factor": sigma_r / sigma_ts, "sigma_reactants": [sigma_r, 1], "sigma_ts": sigma_ts,
            "convention": METHANE_CONVENTION, "source": "methane_h_abstraction_v1",
            "references": ["https://cccbdb.nist.gov/thermox.asp", "https://cccbdb.nist.gov/ts2x.asp?ri=8",
                           "https://web.stanford.edu/group/Zarelab/publinks/764.pdf"]}


def annotate_single_events(network, *, calculator="methane_h_abstraction_v1", factors=(), parameters=None):
    """Annotate channels in place and return the network; leave generation generic.

    Imports select a 1-based channel and reaction ID, and must provide an
    independent factor (or sigma_reactants/sigma_ts), convention, source,
    kinetic_subclass and parameter_class. All changes commit only on success.
    Overall rules remain unsupported. Parameter classes share one registry
    record, initialized with null A/n/Ea for later experimental fitting.
    """
    if calculator not in ("methane_h_abstraction_v1", None):
        raise ValueError(f"Unknown single-event calculator: {calculator}")
    imports = {}
    for record in factors:
        key = (record["reaction_id"], record["channel"])
        if key in imports:
            raise ValueError(f"Duplicate imported single-event factor: {key}")
        imports[key] = record
    registry, annotations = {}, []
    for reaction in network.reactions:
        for index, channel in enumerate(reaction.channels, 1):
            center = reaction_center(channel["mapped_reaction_smiles"])
            annotation = {"status": "unsupported", "reaction_center": center,
                          "reason": "No calculator for this elementary family and state"}
            elementary = reaction.kind == "elementary"
            states = [network.species[sid] for sid in channel["reactants"] + channel["products"]]
            if elementary and calculator and all(s.phase == "gas" and s.site is None for s in states):
                annotation.update(_methane_abstraction(channel, center) or {})
            record = imports.pop((reaction.id, index), None)
            if record is not None:
                if not elementary:
                    raise ValueError(f"{reaction.id}: single-event factors require kind='elementary'")
                required = ("kinetic_subclass", "parameter_class", "convention", "source")
                if any(not record.get(k) for k in required):
                    raise ValueError(f"{reaction.id}: imported factors need {', '.join(required)}")
                annotation = {"reaction_center": center, **{k: deepcopy(record[k]) for k in required},
                              "status": "imported"}
                sigma_r, sigma_ts = record.get("sigma_reactants"), record.get("sigma_ts")
                ratio = None
                if sigma_r is not None or sigma_ts is not None:
                    if not sigma_r or sigma_ts is None or any(not isfinite(v) or v <= 0 for v in [*sigma_r, sigma_ts]):
                        raise ValueError("Symmetry numbers must be finite and positive")
                    ratio = prod(sigma_r) / sigma_ts
                    annotation.update(sigma_reactants=list(sigma_r), sigma_ts=sigma_ts)
                factor = record.get("factor", ratio)
                if factor is None or not isfinite(factor) or factor <= 0:
                    raise ValueError("An independent finite positive single-event factor is required")
                if ratio is not None and not isclose(factor, ratio, rel_tol=1e-9):
                    raise ValueError("Imported factor disagrees with sigma_reactants / sigma_ts")
                annotation["factor"] = factor
            if annotation["status"] != "unsupported":
                annotation.pop("reason", None)
                class_id = annotation["parameter_class"]
                order = sum(n for _, n in reaction.reactants)
                units = "s^-1" if order == 1 else f"(m^3/mol)^{order - 1} s^-1"
                if class_id in registry and registry[class_id]["molecularity"] != order:
                    raise ValueError("A parameter class cannot mix molecularities")
                registry.setdefault(class_id, {"A": None, "n": None, "Ea": None, "T_ref": 1.0,
                                               "A_units": units, "Ea_units": "J/mol", "T_ref_units": "K",
                                               "molecularity": order, "status": "unfitted", "source": ""})
            elif not elementary:
                annotation["reason"] = "Overall transformation; elementary kinetics not assigned"
            annotations.append((channel, annotation))
    if imports:
        raise ValueError(f"Imported factors select unknown reaction/channel: {list(imports)}")
    for class_id, values in (parameters or {}).items():
        if class_id not in registry:
            raise ValueError(f"Unknown parameter class: {class_id}")
        allowed = {"A", "n", "Ea", "T_ref", "status", "source"}
        if values.keys() - allowed:
            raise ValueError(f"Unknown parameter fields: {values.keys() - allowed}")
        for name in ("A", "n", "Ea", "T_ref"):
            value = values.get(name)
            if value is not None and (not isfinite(value) or (name in ("A", "T_ref") and value <= 0)):
                raise ValueError(f"Invalid parameter {name}: {value}")
        if "T_ref" in values and values["T_ref"] is None:
            raise ValueError("T_ref must be positive")
        registry[class_id].update(deepcopy(values))
    for channel, annotation in annotations:
        channel["single_event"] = annotation
    network.single_event = {
        "options": deepcopy({"calculator": calculator, "factors": list(factors), "parameters": parameters or {}}),
        "rate_law": "k_channel(T) = factor * A * (T/T_ref)^n * exp(-Ea/(R*T))",
        "parameters": registry,
        "unmapped_reactions": [r.id for r in network.reactions if not r.channels],
    }
    return network
