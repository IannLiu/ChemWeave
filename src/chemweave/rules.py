"""Declarative reaction rules; RDKit templates are compiled once per rule."""

from collections import Counter
from dataclasses import dataclass, field
from importlib.resources import files
import json
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdChemReactions

from .models import Species


@dataclass(slots=True)
class ReactionRule:
    id: str
    smarts: str
    family: str = "unspecified"
    kind: str = "overall"
    version: str = "1"
    source: str = ""
    notes: str = ""
    exclude_smarts: tuple[str, ...] = ()
    reactant_phases: tuple[str, ...] = ()
    product_phases: tuple[str, ...] = ()
    reactant_sites: tuple[str | None, ...] = ()
    product_sites: tuple[str | None, ...] = ()
    _reaction: object = field(init=False, repr=False)
    _exclusions: tuple = field(init=False, repr=False)

    def __post_init__(self):
        self._reaction = rdChemReactions.ReactionFromSmarts(self.smarts)
        if self._reaction is None:
            raise ValueError(f"Invalid reaction SMARTS for {self.id}")
        self._reaction.Initialize(silent=True)
        nr = self._reaction.GetNumReactantTemplates()
        np = self._reaction.GetNumProductTemplates()
        for attr, size, default in (
            ("reactant_phases", nr, "gas"), ("product_phases", np, "gas"),
            ("reactant_sites", nr, None), ("product_sites", np, None),
        ):
            values = tuple(getattr(self, attr)) or (default,) * size
            if len(values) != size:
                raise ValueError(f"{self.id}: {attr} needs {size} entries")
            setattr(self, attr, values)
        self.exclude_smarts = tuple(self.exclude_smarts)
        self._exclusions = tuple(Chem.MolFromSmarts(s) for s in self.exclude_smarts)
        if any(pattern is None for pattern in self._exclusions):
            raise ValueError(f"Invalid exclusion SMARTS for {self.id}")

    @property
    def arity(self):
        return self._reaction.GetNumReactantTemplates()

    def matches(self, species, slot):
        return (
            species.phase == self.reactant_phases[slot]
            and species.site == self.reactant_sites[slot]
            and species._mol.HasSubstructMatch(self._reaction.GetReactantTemplate(slot))
            and not any(species._mol.HasSubstructMatch(p) for p in self._exclusions)
        )

    def apply(self, reactants, rejected):
        """Yield valid product tuples and their raw template-match provenance.

        ``reactants`` follow template slot order. The generator supplies only
        matching tuples; this operation sanitizes products and checks balances.
        ``rejected`` collects the reasons for discarded template outcomes.
        """
        left_elements = Counter()
        for s in reactants:
            left_elements.update(dict(s.elements))
        left_sites = Counter(s.site for s in reactants if s.site is not None)
        left_charge = sum(s.charge for s in reactants)
        # Explicit index labels make provenance interpretable after JSON export,
        # even when canonical SMILES parsing would change the atom order.
        indexed_reactants = []
        for s in reactants:
            mol = Chem.Mol(s._mol)
            for i, atom in enumerate(mol.GetAtoms()):
                atom.SetAtomMapNum(i + 1)
            indexed_reactants.append(Chem.MolToSmiles(mol))
        # maxProducts=0 avoids RDKit's silent default product cap.
        outcomes = self._reaction.RunReactants(tuple(s._mol for s in reactants), maxProducts=0)
        for raw_products in outcomes:
            atom_maps = [
                [
                    {key: atom.GetIntProp(key) for key in ("react_idx", "react_atom_idx", "old_mapno")
                     if atom.HasProp(key)}
                    for atom in mol.GetAtoms()
                ] for mol in raw_products
            ]
            try:
                with rdBase.BlockLogs():
                    products = tuple(
                        Species.from_mol(mol, phase=self.product_phases[i], site=self.product_sites[i])
                        for i, mol in enumerate(raw_products)
                    )
            except Chem.rdchem.MolSanitizeException:
                rejected["invalid_structure"] += 1
                continue
            right_elements = Counter()
            for s in products:
                right_elements.update(dict(s.elements))
            right_sites = Counter(s.site for s in products if s.site is not None)
            if (left_elements != right_elements
                or left_charge != sum(s.charge for s in products)
                or left_sites != right_sites):
                rejected["unbalanced"] += 1
                continue
            yield products, {
                "indexed_reactant_smiles": indexed_reactants,
                "product_atom_properties": atom_maps,
            }

    def to_dict(self):
        return {
            name: getattr(self, name) for name in (
                "id", "smarts", "family", "kind", "version", "source", "notes",
                "exclude_smarts", "reactant_phases", "product_phases",
                "reactant_sites", "product_sites",
            )
        }


def load_rules(path=None):
    """Load a JSON rule list; no path selects the bundled demo rules."""
    if path is None:
        content = files("chemweave").joinpath("data/demo_rules.json").read_text(encoding="utf-8")
    else:
        content = Path(path).read_text(encoding="utf-8")
    return [ReactionRule(**record) for record in json.loads(content)]
