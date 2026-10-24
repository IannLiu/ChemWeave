"""Species identity, reaction stoichiometry, and network persistence."""

from collections import Counter
from dataclasses import asdict, dataclass, field
import csv
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


@dataclass(frozen=True, slots=True)
class Species:
    """A molecular graph in a specified phase and optional single-site state.

    Use ``from_smiles`` or ``from_mol`` so the structure is standardized once.
    The name is a label; it does not contribute to chemical identity.
    """

    smiles: str
    formula: str
    elements: tuple[tuple[str, int], ...]
    charge: int
    carbon_count: int
    phase: str = "gas"
    site: str | None = None
    name: str = ""
    _mol: Chem.Mol = field(repr=False, compare=False, hash=False, default=None)

    @classmethod
    def from_smiles(cls, smiles, *, name="", phase="gas", site=None):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or mol.GetNumAtoms() == 0:
            raise ValueError(f"Invalid or empty molecular SMILES: {smiles!r}")
        return cls.from_mol(mol, name=name, phase=phase, site=site)

    @classmethod
    def from_mol(cls, mol, *, name="", phase="gas", site=None):
        mol = Chem.Mol(mol)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
            for prop in ("old_mapno", "react_idx", "react_atom_idx"):
                if atom.HasProp(prop):
                    atom.ClearProp(prop)
        # Remove ordinary explicit H atoms and sanitize the reaction product.
        # Molecular H2 and isotope-labelled H are retained by RDKit.
        mol = Chem.RemoveHs(mol)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        elements = Counter(a.GetSymbol() for a in Chem.AddHs(mol).GetAtoms())
        return cls(
            smiles=Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True),
            formula=rdMolDescriptors.CalcMolFormula(mol),
            elements=tuple(sorted(elements.items())),
            charge=Chem.GetFormalCharge(mol),
            carbon_count=elements["C"],
            phase=phase,
            site=site,
            name=name,
            _mol=mol,
        )

    @property
    def key(self):
        return self.smiles, self.phase, self.site


@dataclass(slots=True)
class Reaction:
    """Directed reaction with integer stoichiometry and template provenance.

    Matches record template executions, NOT a SEMK degeneracy factor.
    """

    id: str
    reactants: tuple[tuple[str, int], ...]
    products: tuple[tuple[str, int], ...]
    rule_id: str
    family: str
    kind: str
    matches: list[dict] = field(default_factory=list)


@dataclass(slots=True)
class Network:
    species: dict[str, Species] = field(default_factory=dict)
    reactions: list[Reaction] = field(default_factory=list)
    seeds: list[str] = field(default_factory=list)
    rules: list[dict] = field(default_factory=list)
    config: dict = field(default_factory=dict)
    stats: dict = field(default_factory=dict)
    termination: str = "exhausted"

    def summary(self):
        return {
            "species": len(self.species),
            "reactions": len(self.reactions),
            "termination": self.termination,
            **self.stats,
        }

    def export_json(self, path):
        """Save all chemical records and the exact rules/configuration used."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {
            "schema_version": 1,
            "species": [
                {"id": sid, "smiles": s.smiles, "formula": s.formula,
                 "phase": s.phase, "site": s.site, "name": s.name}
                for sid, s in self.species.items()
            ],
            "reactions": [asdict(r) for r in self.reactions],
            "seeds": self.seeds,
            "rules": self.rules,
            "config": self.config,
            "stats": self.stats,
            "termination": self.termination,
        }
        path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    @classmethod
    def from_json(cls, path):
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        if data["schema_version"] != 1:
            raise ValueError(f"Unsupported network schema: {data['schema_version']}")
        species = {
            record["id"]: Species.from_smiles(
                record["smiles"], name=record["name"],
                phase=record["phase"], site=record["site"],
            ) for record in data["species"]
        }
        reactions = [
            Reaction(**{
                **record,
                "reactants": tuple((sid, n) for sid, n in record["reactants"]),
                "products": tuple((sid, n) for sid, n in record["products"]),
            }) for record in data["reactions"]
        ]
        return cls(
            species=species, reactions=reactions, seeds=data["seeds"],
            rules=data["rules"], config=data["config"],
            stats=data["stats"], termination=data["termination"],
        )

    def export_csv(self, directory):
        """Write species, reaction, and long-form stoichiometry tables."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        with (directory / "species.csv").open("w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["id", "name", "smiles", "formula", "charge", "phase", "site"])
            for sid, s in self.species.items():
                writer.writerow([sid, s.name, s.smiles, s.formula, s.charge, s.phase, s.site])
        with (directory / "reactions.csv").open("w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["id", "rule_id", "family", "kind", "reactants", "products", "template_matches"])
            for r in self.reactions:
                left = " + ".join(f"{n} {sid}" for sid, n in r.reactants)
                right = " + ".join(f"{n} {sid}" for sid, n in r.products)
                writer.writerow([r.id, r.rule_id, r.family, r.kind, left, right, len(r.matches)])
        with (directory / "stoichiometry.csv").open("w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["reaction_id", "side", "species_id", "coefficient"])
            for r in self.reactions:
                for side, terms in (("reactant", r.reactants), ("product", r.products)):
                    writer.writerows((r.id, side, sid, n) for sid, n in terms)

    def export_dot(self, path):
        """Use reaction nodes to preserve multi-reactant/multi-product edges."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        lines = ["digraph reactions {", "  rankdir=LR;"]
        for sid, s in self.species.items():
            label = f"{sid}: {s.name or s.smiles}"
            if s.site is not None:
                label += f" [{s.phase}, {s.site}]"
            lines.append(f"  {sid} [label={json.dumps(label, ensure_ascii=False)}];")
        for r in self.reactions:
            label = json.dumps(f"{r.id}: {r.rule_id}")
            lines.append(f"  {r.id} [shape=box, label={label}];")
            for sid, n in r.reactants:
                lines.append(f'  {sid} -> {r.id} [label="{n}"];')
            for sid, n in r.products:
                lines.append(f'  {r.id} -> {sid} [label="{n}"];')
        lines.append("}")
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
