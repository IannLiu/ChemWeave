"""Explicit abstract catalyst nodes, independent of kinetic assumptions.

Atomic number zero in a species means a catalyst node, never a query wildcard.
Metal nodes bind adsorbates; Brønsted anchors bind an explicit framework O.
All metadata is immutable through a reaction; protonation lives on O/H atoms.
"""

import json

from rdkit import Chem


SITE_PROPERTY = "_chemweave_site"


def site_label(atom):
    return atom.GetProp(SITE_PROPERTY) if atom.HasProp(SITE_PROPERTY) else ""


def set_site_atoms(mol, records):
    """Attach JSON attributes keyed by the supplied, unique atom-map numbers."""
    atoms = {a.GetAtomMapNum(): a for a in mol.GetAtoms() if a.GetAtomMapNum()}
    for number, attributes in records.items():
        number = int(number)
        if number not in atoms or atoms[number].GetAtomicNum() != 0:
            raise ValueError("site_atoms must select mapped atomic-number-zero catalyst nodes")
        if sum(a.GetAtomMapNum() == number for a in mol.GetAtoms()) != 1:
            raise ValueError("Explicit site map numbers must be unique")
        attributes = dict(attributes)
        attributes.setdefault("kind", "surface")
        attributes.setdefault("type", "X")
        if attributes["kind"] not in {"surface", "metal", "bronsted"} or not attributes["type"]:
            raise ValueError("Site kind must be surface/metal/bronsted and type must be nonempty")
        attributes.setdefault("capacity", 1)
        if attributes["capacity"] != 1:
            raise ValueError("Each explicit node represents one site; use multiple nodes for multi-site adsorption")
        atoms[number].SetProp(SITE_PROPERTY, json.dumps(attributes, sort_keys=True, separators=(",", ":")))


def site_identity(mol, phase):
    """Validate nodes and return a numbering-independent colored identity."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0 and not site_label(atom):
            atom.SetProp(SITE_PROPERTY, json.dumps({"capacity": 1, "kind": "surface", "type": "X"}, sort_keys=True, separators=(",", ":")))
    labels = sorted({site_label(a) for a in mol.GetAtoms() if site_label(a)})
    colored = Chem.Mol(mol)
    for atom in colored.GetAtoms():
        atom.SetAtomMapNum(0)
        if atom.GetAtomicNum() != 0:
            continue
        label = site_label(atom)
        if not label:
            raise ValueError("Every '*' species atom needs explicit site_atoms metadata")
        attributes = json.loads(label)
        if phase != "surface" or atom.GetDegree() > attributes["capacity"]:
            raise ValueError("Catalyst nodes require surface phase and cannot exceed site capacity")
        if atom.GetIsotope() or atom.GetFormalCharge() or atom.GetNumRadicalElectrons():
            raise ValueError("Catalyst nodes are neutral anchors; place charge on explicit chemical atoms")
        if any(b.GetBondTypeAsDouble() not in {1, 2, 3, 4} for b in atom.GetBonds()):
            raise ValueError("Abstract sites support single through quadruple attachment bonds")
        if attributes["kind"] == "bronsted" and (atom.GetDegree() != 1 or atom.GetNeighbors()[0].GetAtomicNum() != 8):
            raise ValueError("A Brønsted anchor must attach to an explicit framework oxygen")
        atom.SetIsotope(labels.index(label) + 1)
    return json.dumps([labels, Chem.MolToSmiles(colored)], separators=(",", ":")) if labels else ""
