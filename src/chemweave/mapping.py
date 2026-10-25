"""Reaction-instance numbering and whole-reaction symmetry signatures."""

import hashlib
import json

from rdkit import Chem


def prepare_reactants(species):
    """Give every atom a unique reaction-local map, retaining source labels.

    Species already have explicit, numbered hydrogens. Reserve all their labels
    before replacing collisions, so a replacement cannot steal a later label.
    Repeated molecules are separate occurrences with separate source records.
    """
    next_map = max((a.GetAtomMapNum() for s in species for a in s._mapped_mol.GetAtoms()), default=0) + 1
    used = set()
    mols, origins = [], []
    for slot, s in enumerate(species):
        mol = Chem.Mol(s._mapped_mol)
        input_maps = dict(s.input_maps)
        for atom in mol.GetAtoms():
            source_map = atom.GetAtomMapNum()
            reaction_map = source_map
            if reaction_map in used:
                reaction_map = next_map
                next_map += 1
            used.add(reaction_map)
            atom.SetAtomMapNum(reaction_map)
            origins.append({
                "reactant": slot, "atom_index": atom.GetIdx(),
                "species_map": source_map, "input_map": input_maps[source_map],
                "reaction_map": reaction_map,
            })
        mols.append(mol)
    return tuple(mols), origins


def reaction_symmetry_key(reactants, products, reactant_states, product_states):
    """Canonicalize a colored graph of the complete before/after reaction.

    The same atom node participates in both sides. Separate component nodes
    preserve molecular occurrences and phase/site labels. Bond nodes encode
    before/after bond orders and defined stereo. This tests joint symmetry,
    rather than independently replacing equivalent atoms by the same map.
    Numbering is ignored; isotope, charge, radicals and assigned stereo remain.
    No list of molecular automorphisms is enumerated.
    """
    before = {a.GetAtomMapNum(): a for mol in reactants for a in mol.GetAtoms()}
    after = {a.GetAtomMapNum(): a for mol in products for a in mol.GetAtoms()}
    labels, edges, atom_nodes = [], [], {}
    for number, atom in before.items():
        descriptors = []
        for a in (atom, after[number]):
            descriptors.append((a.GetAtomicNum(), a.GetIsotope(), a.GetFormalCharge(),
                                a.GetNumRadicalElectrons(), a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else ""))
        atom_nodes[number] = len(labels)
        labels.append(("atom", *descriptors))
    for side, mols, states in (("reactant", reactants, reactant_states), ("product", products, product_states)):
        for mol, state in zip(mols, states):
            component = len(labels)
            labels.append((side, *state))
            edges.extend((component, atom_nodes[a.GetAtomMapNum()]) for a in mol.GetAtoms())
    bonds = {}
    for side, mols in enumerate((reactants, products)):
        for mol in mols:
            for bond in mol.GetBonds():
                pair = tuple(sorted((bond.GetBeginAtom().GetAtomMapNum(), bond.GetEndAtom().GetAtomMapNum())))
                bonds.setdefault(pair, [None, None])[side] = (str(bond.GetBondType()), str(bond.GetStereo()))
    for (a, b), states in bonds.items():
        node = len(labels)
        labels.append(("bond", *states))
        edges.extend(((node, atom_nodes[a]), (node, atom_nodes[b])))
    # Dummy atoms with isotope colors encode a graph, not a chemical molecule.
    # Include the color vocabulary in the key to distinguish different labels.
    text_labels = [json.dumps(label, separators=(",", ":")) for label in labels]
    vocabulary = sorted(set(text_labels))
    colors = {label: i + 1 for i, label in enumerate(vocabulary)}
    graph = Chem.RWMol()
    for label in text_labels:
        atom = Chem.Atom(0)
        atom.SetIsotope(colors[label])
        graph.AddAtom(atom)
    for a, b in edges:
        graph.AddBond(a, b, Chem.BondType.SINGLE)
    canonical = Chem.MolToSmiles(graph.GetMol(), canonical=True)
    payload = json.dumps((vocabulary, canonical), separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()
