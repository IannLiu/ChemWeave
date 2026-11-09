"""Species identity, reaction stoichiometry, and network persistence."""

from collections import Counter
from dataclasses import asdict, dataclass, field, replace
import csv
import json
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdMolDescriptors

from .sites import set_site_atoms, site_identity, site_label


@dataclass(frozen=True, slots=True)
class Species:
    """A molecular graph with phase, explicit catalyst nodes and a legacy site label.

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
    mapped_smiles: str = ""
    input_maps: tuple[tuple[int, int], ...] = ()
    symmetry_classes: tuple[tuple[int, int], ...] = ()
    site_atoms: dict = field(default_factory=dict, compare=False)
    _site_identity: str = field(default="", repr=False)
    _mol: Chem.Mol = field(repr=False, compare=False, hash=False, default=None)
    _mapped_mol: Chem.Mol = field(repr=False, compare=False, hash=False, default=None)

    @classmethod
    def from_smiles(cls, smiles, *, name="", phase="gas", site=None, site_atoms=None):
        params = Chem.SmilesParserParams()
        params.removeHs = False
        mol = Chem.MolFromSmiles(smiles, params)
        if mol is None or mol.GetNumAtoms() == 0:
            raise ValueError(f"Invalid or empty molecular SMILES: {smiles!r}")
        set_site_atoms(mol, site_atoms or {})
        return cls.from_mol(mol, name=name, phase=phase, site=site)

    @classmethod
    def from_mol(cls, mol, *, name="", phase="gas", site=None):
        mol = Chem.Mol(mol)
        Chem.SanitizeMol(mol)
        mol = Chem.AddHs(mol)
        next_map = max((a.GetAtomMapNum() for a in mol.GetAtoms()), default=0) + 1
        used, input_maps = set(), []
        for atom in mol.GetAtoms():
            original = atom.GetAtomMapNum()
            number = original
            if not number or number in used:
                number = next_map
                next_map += 1
            used.add(number)
            atom.SetAtomMapNum(number)
            input_maps.append((number, original))
            for prop in ("old_mapno", "react_idx", "react_atom_idx"):
                if atom.HasProp(prop):
                    atom.ClearProp(prop)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        identity = site_identity(mol, phase)
        mapped_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        ranked = Chem.Mol(mol)
        labels = sorted({site_label(a) for a in mol.GetAtoms() if site_label(a)})
        for a in ranked.GetAtoms():
            if site_label(a):
                a.SetIsotope(labels.index(site_label(a)) + 1)
        ranks = Chem.CanonicalRankAtoms(ranked, breakTies=False, includeAtomMaps=False)
        symmetry_classes = tuple(sorted((a.GetAtomMapNum(), rank) for a, rank in zip(mol.GetAtoms(), ranks)))
        elements = Counter(a.GetSymbol() for a in mol.GetAtoms())
        structure = Chem.Mol(mol)
        for atom in structure.GetAtoms():
            atom.SetAtomMapNum(0)
        with rdBase.BlockLogs():
            structure = Chem.RemoveHs(structure)
        smiles = Chem.MolToSmiles(structure, canonical=True, isomericSmiles=True)
        # Reparse the canonical compact structure to give registry atoms a
        # consistent order and normalize bracketed versus implicit H counts.
        with rdBase.BlockLogs():
            structure = Chem.MolFromSmiles(smiles)
        return cls(
            smiles=smiles,
            formula=rdMolDescriptors.CalcMolFormula(mol),
            elements=tuple(sorted(elements.items())),
            charge=Chem.GetFormalCharge(mol),
            carbon_count=elements["C"],
            phase=phase,
            site=site,
            name=name,
            mapped_smiles=mapped_smiles,
            input_maps=tuple(sorted(input_maps)),
            symmetry_classes=symmetry_classes,
            site_atoms={a.GetAtomMapNum(): json.loads(site_label(a)) for a in mol.GetAtoms() if site_label(a)},
            _site_identity=identity,
            _mol=structure,
            _mapped_mol=mol,
        )

    @property
    def key(self):
        return self.smiles, self.phase, self.site, self._site_identity


@dataclass(slots=True)
class Reaction:
    """Directed reaction with integer stoichiometry and template provenance.

    Each channel stores one complete mapping and an equivalent-match count.
    The count is a template execution statistic, not a SEMK degeneracy factor.
    """

    id: str
    reactants: tuple[tuple[str, int], ...]
    products: tuple[tuple[str, int], ...]
    rule_id: str
    family: str
    kind: str
    channels: list[dict] = field(default_factory=list)
    depth: int | None = None


@dataclass(slots=True)
class Network:
    species: dict[str, Species] = field(default_factory=dict)
    reactions: list[Reaction] = field(default_factory=list)
    seeds: list[str] = field(default_factory=list)
    seed_inputs: list[dict] = field(default_factory=list)
    species_depths: dict[str, int] = field(default_factory=dict)
    targets: list[dict] = field(default_factory=list)
    rules: list[dict] = field(default_factory=list)
    config: dict = field(default_factory=dict)
    stats: dict = field(default_factory=dict)
    termination: str = "exhausted"
    single_event: dict | None = None

    def summary(self):
        return {
            "species": len(self.species),
            "reactions": len(self.reactions),
            "channels": sum(len(r.channels) for r in self.reactions),
            "template_matches": sum(c["match_count"] for r in self.reactions for c in r.channels),
            "termination": self.termination,
            "depth_reached": max(self.species_depths.values(), default=None),
            "targets_total": len(self.targets),
            "targets_found": sum(t["species_id"] is not None for t in self.targets),
            "missing_targets": [t for t in self.targets if t["species_id"] is None],
            **self.stats,
        }

    def to_dict(self):
        """Return the persisted representation, shared by files and the web UI."""
        return {
            "schema_version": 6,
            "species": [
                {"id": sid, "smiles": s.smiles, "formula": s.formula,
                 "phase": s.phase, "site": s.site, "name": s.name,
                 "mapped_smiles": s.mapped_smiles, "input_maps": s.input_maps,
                 "symmetry_classes": s.symmetry_classes,
                 "site_atoms": s.site_atoms,
                 "depth": self.species_depths.get(sid)}
                for sid, s in self.species.items()
            ],
            "reactions": [asdict(r) for r in self.reactions],
            "seeds": self.seeds,
            "seed_inputs": self.seed_inputs,
            "targets": self.targets,
            "rules": self.rules,
            "config": self.config,
            "stats": self.stats,
            "termination": self.termination,
            "single_event": self.single_event,
        }

    def export_json(self, path):
        """Save all chemical records and the exact rules/configuration used."""
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    @classmethod
    def from_json(cls, path):
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))

    @classmethod
    def from_dict(cls, data):
        """Restore a saved network without requiring an intermediate file."""
        if data["schema_version"] not in (1, 2, 3, 4, 5, 6):
            raise ValueError(f"Unsupported network schema: {data['schema_version']}")
        species = {}
        for record in data["species"]:
            s = Species.from_smiles(
                record.get("mapped_smiles", record["smiles"]), name=record["name"],
                phase=record["phase"], site=record["site"],
                site_atoms=record.get("site_atoms"),
            )
            if "input_maps" in record:
                s = replace(s, input_maps=tuple(tuple(pair) for pair in record["input_maps"]))
            species[record["id"]] = s
        reactions = []
        for record in data["reactions"]:
            record = dict(record)
            if data["schema_version"] == 2:
                matches = record.pop("matches", [])
                record["channels"] = [
                    {**matches[c["match_indices"][0]], "match_count": len(c["match_indices"])}
                    for c in record.get("channels", [])
                ]
            elif data["schema_version"] == 1:
                # Version 1 has no complete mapped reaction representatives.
                record.pop("matches", None)
                record["channels"] = []
            reactions.append(Reaction(**{
                **record,
                "reactants": tuple((sid, n) for sid, n in record["reactants"]),
                "products": tuple((sid, n) for sid, n in record["products"]),
            }))
        return cls(
            species=species, reactions=reactions, seeds=data["seeds"],
            seed_inputs=data.get("seed_inputs", []),
            species_depths={r["id"]: r["depth"] for r in data["species"] if r.get("depth") is not None},
            targets=data.get("targets", []),
            rules=data["rules"], config=data["config"],
            stats=data["stats"], termination=data["termination"],
            single_event=data.get("single_event"),
        )

    def export_csv(self, directory):
        """Write species, reaction, and long-form stoichiometry tables."""
        directory = Path(directory)
        directory.mkdir(parents=True, exist_ok=True)
        with (directory / "species.csv").open("w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["id", "name", "smiles", "formula", "charge", "phase", "site", "mapped_smiles", "depth", "site_atoms"])
            for sid, s in self.species.items():
                writer.writerow([sid, s.name, s.smiles, s.formula, s.charge, s.phase, s.site, s.mapped_smiles,
                                 self.species_depths.get(sid), json.dumps(s.site_atoms)])
        with (directory / "reactions.csv").open("w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["id", "rule_id", "family", "kind", "reactants", "products", "template_matches", "channels", "depth"])
            for r in self.reactions:
                left = " + ".join(f"{n} {sid}" for sid, n in r.reactants)
                right = " + ".join(f"{n} {sid}" for sid, n in r.products)
                count = sum(c["match_count"] for c in r.channels)
                writer.writerow([r.id, r.rule_id, r.family, r.kind, left, right, count, len(r.channels), r.depth])
        with (directory / "channels.csv").open("w", encoding="utf-8-sig", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["reaction_id", "channel", "mapped_reaction_smiles", "template_matches"])
            for r in self.reactions:
                for i, channel in enumerate(r.channels, 1):
                    writer.writerow([r.id, i, channel["mapped_reaction_smiles"], channel["match_count"]])
        if self.single_event is not None:
            with (directory / "single_event.csv").open("w", encoding="utf-8-sig", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["reaction_id", "channel", "status", "kinetic_subclass", "parameter_class",
                                 "factor", "sigma_reactants", "sigma_ts", "convention", "source", "reason", "reaction_center",
                                 "method", "family", "rate_basis", "kinetic_channel", "counting",
                                 "n_e", "algorithm", "reason_code", "event", "unique_labeled_paths", "identical_reactant_divisor"])
                for r in self.reactions:
                    for i, channel in enumerate(r.channels, 1):
                        se = channel.get("single_event", {})
                        writer.writerow([r.id, i, *[se.get(k) for k in ("status", "kinetic_subclass", "parameter_class", "factor")],
                                         json.dumps(se.get("sigma_reactants")), se.get("sigma_ts"),
                                         se.get("convention"), se.get("source"), se.get("reason"),
                                         json.dumps(se.get("reaction_center")),
                                         *[se.get(k) for k in ("method", "family", "rate_basis", "kinetic_channel")],
                                         json.dumps(se.get("counting")), se.get("n_e", se.get("factor")),
                                         se.get("algorithm"), se.get("reason_code"), json.dumps(se.get("event")),
                                         (se.get("counting") or {}).get("unique_labeled_paths"),
                                         (se.get("counting") or {}).get("identical_reactant_divisor")])
            with (directory / "kinetic_parameters.csv").open("w", encoding="utf-8-sig", newline="") as f:
                writer = csv.writer(f)
                fields = ["A", "n", "Ea", "T_ref", "A_units", "Ea_units", "T_ref_units", "molecularity", "status", "source", "rate_basis"]
                writer.writerow(["parameter_class", *fields])
                for class_id, record in self.single_event["parameters"].items():
                    writer.writerow([class_id, *[record.get(k) for k in fields]])
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
