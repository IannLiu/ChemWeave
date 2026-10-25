"""Declarative reaction rules; RDKit templates are compiled once per rule."""

from collections import Counter
from dataclasses import dataclass, field
from importlib.resources import files
import json
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import rdChemReactions

from .models import Species
from .mapping import prepare_reactants, reaction_symmetry_key


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
    _template_sources: dict = field(init=False, repr=False)

    def __post_init__(self):
        self._reaction = rdChemReactions.ReactionFromSmarts(self.smarts)
        if self._reaction is None:
            raise ValueError(f"Invalid reaction SMARTS for {self.id}")
        self._reaction.Initialize(silent=True)
        self._reaction.GetSubstructParams().useChirality = True
        self._template_sources = {
            atom.GetAtomMapNum(): slot
            for slot, mol in enumerate(self._reaction.GetReactants()) for atom in mol.GetAtoms()
            if atom.GetAtomMapNum()
        }
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
            and species._mapped_mol.HasSubstructMatch(self._reaction.GetReactantTemplate(slot), useChirality=True)
            and not any(species._mapped_mol.HasSubstructMatch(p) for p in self._exclusions)
        )

    def apply(self, reactants, rejected):
        """Execute explicit-H templates, restore maps and describe each match.

        Matching tuples are supplied in template order. Each source atom must
        appear exactly once in the products, with element and isotope retained.
        Symmetry keys describe the full reaction, including both molecular sides.
        """
        mols, origins = prepare_reactants(reactants)
        sources = {a.GetAtomMapNum(): a for mol in mols for a in mol.GetAtoms()}
        mapped_reactants = [Chem.MolToSmiles(mol) for mol in mols]
        reactant_states = [(s.phase, s.site) for s in reactants]
        product_states = list(zip(self.product_phases, self.product_sites))
        left_sites = Counter(s.site for s in reactants if s.site is not None)
        left_charge = sum(s.charge for s in reactants)
        right_sites = Counter(site for site in self.product_sites if site is not None)
        symmetry_keys = {}
        for raw_products in self._reaction.RunReactants(mols, maxProducts=0):
            try:
                with rdBase.BlockLogs():
                    for mol in raw_products:
                        Chem.SanitizeMol(mol)
                    product_mols = tuple(Chem.AddHs(mol) for mol in raw_products)
                    for mol in product_mols:
                        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
            except Chem.rdchem.MolSanitizeException:
                rejected["invalid_structure"] += 1
                continue
            # Restore source labels on both template atoms and untouched context.
            atom_properties, numbers, valid = [], [], True
            for mol in product_mols:
                properties = []
                for atom in mol.GetAtoms():
                    props = {key: atom.GetIntProp(key)
                             for key in ("react_idx", "react_atom_idx", "old_mapno") if atom.HasProp(key)}
                    slot = props.get("react_idx", self._template_sources.get(props.get("old_mapno")))
                    if slot is not None and "react_atom_idx" in props:
                        source = mols[slot].GetAtomWithIdx(props["react_atom_idx"])
                        number = source.GetAtomMapNum()
                    else:
                        number = atom.GetAtomMapNum()
                        source = sources.get(number)
                    if source is None or number not in sources:
                        valid = False
                        break
                    if (atom.GetAtomicNum(), atom.GetIsotope()) != (source.GetAtomicNum(), source.GetIsotope()):
                        valid = False
                        break
                    atom.SetAtomMapNum(number)
                    numbers.append(number)
                    properties.append({**props, "map_num": number})
                atom_properties.append(properties)
                if not valid:
                    break
            if not valid or len(numbers) != len(sources) or set(numbers) != set(sources):
                rejected["atom_mapping"] += 1
                continue
            if left_charge != sum(Chem.GetFormalCharge(m) for m in product_mols) or left_sites != right_sites:
                rejected["unbalanced"] += 1
                continue
            products = tuple(
                Species.from_mol(mol, phase=self.product_phases[i], site=self.product_sites[i])
                for i, mol in enumerate(product_mols)
            )
            mapped_products = [Chem.MolToSmiles(mol) for mol in product_mols]
            mapped_reaction = ".".join(mapped_reactants) + ">>" + ".".join(mapped_products)
            if mapped_reaction not in symmetry_keys:
                symmetry_keys[mapped_reaction] = reaction_symmetry_key(
                    mols, product_mols, reactant_states, product_states,
                )
            yield products, {
                "mapped_reaction_smiles": mapped_reaction,
                "mapped_reactant_smiles": mapped_reactants,
                "mapped_product_smiles": mapped_products,
                "atom_mapping": origins,
                "product_atom_properties": atom_properties,
                "symmetry_key": symmetry_keys[mapped_reaction],
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
