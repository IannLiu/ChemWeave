"""Opt-in family-specific calculators, separate from the annotation framework.

Each calculator owns its chemical applicability, subclassification and symmetry
convention. Methane H/D abstraction is the first validation implementation.
"""

from collections import Counter

from rdkit import Chem


METHANE_CONVENTION = (
    "Forward gas-phase H/D atom + methane(H/D): tetrahedral methane; axial "
    "C--H/D--H/D saddle with three spectator ligands. Proper rotations only: "
    "sigma_R = 12 (4 identical), 3 (3+1), 2 (2+2); sigma_TS = 3 when all "
    "spectators are identical, otherwise 1. Atomic reactant sigma = 1. "
    "n_SE = product(sigma_R)/sigma_TS; no product symmetry, spin, tunneling, "
    "or nuclear-spin weights included. Intrinsic parameters exclude n_SE."
)


def methane_h_abstraction(*, network, reaction, channel, center):
    """Recognize the complete elementary transformation, independent of rule ID."""
    states = [network.species[sid] for sid in channel["reactants"] + channel["products"]]
    if any(s.phase != "gas" or s.site is not None for s in states):
        return None
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
    return {"kinetic_subclass": subclass, "parameter_class": subclass,
            "factor": sigma_r / sigma_ts, "sigma_reactants": [sigma_r, 1], "sigma_ts": sigma_ts,
            "convention": METHANE_CONVENTION, "source": "methane_h_abstraction_v1",
            "references": ["https://cccbdb.nist.gov/thermox.asp", "https://cccbdb.nist.gov/ts2x.asp?ri=8",
                           "https://web.stanford.edu/group/Zarelab/publinks/764.pdf"]}

