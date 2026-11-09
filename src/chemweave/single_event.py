"""Family-independent reaction-center and single-event annotation framework."""

from copy import deepcopy
from math import isclose, isfinite, prod

from rdkit import Chem


_CALCULATORS = {}


def register_single_event_calculator(name, calculator):
    """Register a named callback for explicit selection and saved options.

    The callback receives keyword arguments network, reaction, channel, center.
    Return a factor with convention/source, a structured non-success status,
    or None when inapplicable. Kinetic classes and units are optional metadata.
    Register custom callbacks again before replaying a saved network's options.
    """
    _CALCULATORS[name] = calculator


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


def annotate_single_events(network, *, calculator=None, factors=(), parameters=None):
    """Annotate channels in place and return the network; leave generation generic.

    Imports select a 1-based channel and reaction ID, and must provide an
    independent factor (or sigma_reactants/sigma_ts), convention and source.
    Kinetic classes and units are optional. All changes commit only on success.
    No family calculator is selected by default. Overall rules remain
    unsupported. Parameter classes share one registry record, initialized with null A/n/Ea for later experimental fitting.
    """
    if calculator is not None and calculator not in _CALCULATORS:
        raise ValueError(f"Unknown single-event calculator: {calculator}")
    calculate = _CALCULATORS.get(calculator)
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
                          "family": reaction.family, "method": "not_calculated",
                          "reason": "No applicable calculator selected for this mapped reaction"}
            # Elementary character is an input assumption, not an eligibility
            # test. Legacy 'kind' metadata must not block saved networks either.
            imported = imports.pop((reaction.id, index), None)
            record = imported
            if record is None and calculate is not None:
                record = calculate(network=network, reaction=reaction, channel=channel, center=center)
            if imported is None and record is not None and record.get("status", "calculated") != "calculated":
                annotation.update(record)
                record = None
            if record is not None:
                required = ("convention", "source")
                if any(not record.get(k) for k in required):
                    raise ValueError(f"{reaction.id}: factor assignments need {', '.join(required)}")
                annotation = {"reaction_center": center, **{k: deepcopy(record[k]) for k in required},
                              "status": "imported" if imported is not None else "calculated"}
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
                    raise ValueError("Assigned factor disagrees with sigma_reactants / sigma_ts")
                annotation["factor"] = factor
                annotation["n_e"] = factor
                annotation["method"] = record.get("method", "symmetry_ratio" if ratio is not None else "imported_factor")
                annotation["family"] = record.get("family", reaction.family)
                for name in ("counting", "event", "algorithm", "kinetic_subclass", "parameter_class", "rate_basis", "A_units", "kinetic_channel"):
                    if name in record:
                        annotation[name] = deepcopy(record[name])
                if "references" in record:
                    annotation["references"] = deepcopy(record["references"])
            if annotation["status"] in {"calculated", "imported"}:
                annotation.pop("reason", None)
            if annotation["status"] in {"calculated", "imported"} and annotation.get("parameter_class"):
                class_id = annotation["parameter_class"]
                order = sum(n for _, n in reaction.reactants)
                gas = all(network.species[sid].phase == "gas" for sid, _ in (*reaction.reactants, *reaction.products))
                units = ("s^-1" if order == 1 else f"(m^3/mol)^{order - 1} s^-1") if gas else None
                units = annotation.get("A_units", units)
                basis = annotation.get("rate_basis", "gas_concentration" if gas else None)
                if class_id in registry and (registry[class_id]["A_units"] != units or registry[class_id]["rate_basis"] != basis):
                    raise ValueError("A parameter class cannot mix rate normalization conventions")
                if class_id in registry and registry[class_id]["molecularity"] != order:
                    raise ValueError("A parameter class cannot mix molecularities")
                registry.setdefault(class_id, {"A": None, "n": None, "Ea": None, "T_ref": 1.0,
                                               "A_units": units, "Ea_units": "J/mol", "T_ref_units": "K",
                                               "molecularity": order, "status": "unfitted", "source": ""})
                registry[class_id]["rate_basis"] = basis
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
    kinetic_channels = {}
    for channel, annotation in annotations:
        channel["single_event"] = annotation
        if "kinetic_channel" in annotation:
            key = annotation["kinetic_channel"]
            kinetic_channels.setdefault(key, {"factor": annotation["factor"], "parameter_class": annotation.get("parameter_class")})
    network.single_event = {
        "options": deepcopy({"calculator": calculator, "factors": list(factors), "parameters": parameters or {}}),
        **({"rate_law": "k_channel(T) = factor * A * (T/T_ref)^n * exp(-Ea/(R*T))"} if registry else {}),
        "parameters": registry,
        "kinetic_channels": kinetic_channels,
        "unmapped_reactions": [r.id for r in network.reactions if not r.channels],
    }
    return network
