**English** | [简体中文](README.zh-CN.md)

# ChemWeave

A Python package for generating chemical reaction networks from RDKit reaction templates. Provide seed species, reaction rules, and a carbon-number limit to produce species tables, reaction tables, stoichiometry, and rule provenance. Rules are maintained separately from the enumeration algorithm so the chemistry can evolve as the mechanism is developed.

The current release provides network generation. Four bundled rules demonstrate overall structural transformations: alkane dehydrogenation, olefin double-bond migration, restricted terminal-olefin coupling, and olefin scission. **These examples do not cover a complete aromatization network and are not a validated elementary mechanism for Ga/H-ZSM-5.**

## Installation and testing

Requires Python 3.11 or later and RDKit. You can use an existing environment, such as Conda's `base` or `ml` environment. From the directory containing the project:

```bash
conda activate base
cd chemweave
python -m pip install -e .
python -m pytest -q
```

pytest is declared as an optional development dependency in `.[test]`. To install it with the package, use `python -m pip install -e ".[test]"`. An environment that already has pytest can use the commands above.

You can also run directly from the project directory without installing the package:

```bash
PYTHONPATH=src python -m chemweave generate \
  --config examples/propane.toml --output results/propane
```

## Python API

```python
from chemweave import Species, Generator, GenerationConfig, load_rules, Network

network = Generator(
    rules=load_rules(),
    config=GenerationConfig(max_carbon=6, max_species=200, max_reactions=1000),
).generate(seeds=[Species.from_smiles("CCC", name="propane")])

print(network.summary())
network.export_json("results/network.json")
network.export_csv("results")
network.export_dot("results/network.dot")

restored = Network.from_json("results/network.json")
```

`load_rules()` loads the bundled examples; `load_rules("my_rules.json")` loads your own rule file. Create `Species` objects with `from_smiles` or `from_mol` to standardize their structures, molecular formulas, and identities.

The C6 example generated 22 species, 45 reactions, and 45 symmetry-grouped channels in the tested `base` environment (Python 3.12, RDKit 2026.03.6). These counts illustrate a run; they do not establish chemical validity.

## Custom rules

Rules are stored as a JSON list. This example represents the overall transformation of propane dehydrogenation, including H2:

```json
[
  {
    "id": "propane_dehydrogenation",
    "smarts": "[C;H3:1](-[H:4])-[C;H2:2](-[H:5])-[C;H3:3]>>[C:1]=[C:2]-[C:3].[H:4]-[H:5]",
    "family": "dehydrogenation",
    "kind": "overall",
    "version": "1",
    "source": "Add the literature reference or rule rationale here",
    "notes": "Overall transformation; the catalytic cycle is not expanded"
  }
]
```

You can also define rules directly in Python with `ReactionRule(...)`. The main fields are:

| Field | Meaning |
| --- | --- |
| `id` | Unique rule name. Different mechanisms can use different rule IDs even when their overall stoichiometry is identical. |
| `smarts` | RDKit Reaction SMARTS; atom-map numbers specify the structural transformation. |
| `family` | Reaction-family label. |
| `kind` | Step classification, defaulting to `overall`. Use `elementary` after reviewing the mechanism. |
| `version`, `source`, `notes` | Rule version, supporting evidence, and assumptions. |
| `exclude_smarts` | Skip the rule if any reactant matches any of these patterns. |
| `reactant_phases`, `product_phases` | Species phases in template order; all default to `gas`. |
| `reactant_sites`, `product_sites` | Single-site labels in template order; all default to `null`. |

Rules are not automatically reversed. Supply a template for the reverse reaction when needed; some double-bond migration templates can themselves generate both directions. Setting `kind` to `elementary` does not mean the program has established the mechanism.

## Network expansion

1. Standardize seed structures and identify species by structure, phase, and site label.
2. Match each new species against each reactant slot of each rule once.
3. Combine newly added candidates with existing candidates and execute ordered reactant tuples that have not been processed. A species may occupy two slots to represent self-coupling.
4. Execute templates with explicit hydrogens. Restore source atom maps and check a one-to-one atom inventory (including elements and isotopes), charge balance, and site-label counts.
5. Apply `max_carbon`, then merge duplicate reactions by rule ID and stoichiometry. Store one complete representative per symmetry-distinct mapping in `channels`; equivalent executions increment its `match_count` without retaining duplicate records.
6. Queue new species until the queue is exhausted or the species/reaction budget is reached.

Checks occur when inputs are read and candidate reactions are generated; the full network is not repeatedly revalidated during enumeration. Molecular descriptors and reaction templates are prepared when their objects are created. Enumeration does not estimate barriers or fluxes, and species are not removed merely because they have not been detected experimentally.

`max_carbon` applies to each species, including the seeds. It defines the scope of a run; it does not establish that pathways outside that scope are unimportant in the actual system.

## Surface states and single-event modeling

You can create `Species.from_smiles("C=CC", phase="surface", site="BAS")`. Site labels contribute to species identity. Rules can constrain input and output sites, and the generator checks that the number of each site label is conserved.

This is an **abstract single-site representation**, suitable for transformations between explicitly defined molecular structures on the same site. Dedicated representations and bundled rules for vacant sites, proton exchange, two-site species, and Ga catalytic cycles are not yet included. Surface structures, hydrogen, and charge must be defined consistently when extending the mechanism; attaching a label to a gas-phase molecule alone does not establish an adsorption mechanism.

The abstract site labels do not change atom numbering or symmetry equivalence: phases and site labels are included in reaction-channel comparison.

## Atom mapping and symmetry

Species identity uses unmapped canonical SMILES, phase, and site. `species.mapped_smiles` is a separate explicit-H representation. Existing non-conflicting map numbers are retained; missing numbers, including those of added hydrogens, are allocated above the largest supplied number. Repeated labels within an input are disambiguated. `species.input_maps` records `(assigned_map, supplied_map)` pairs; a supplied value of zero means no label was provided.

```python
seed = Species.from_smiles("[CH3:101][CH2:205][CH3:309]")
print(seed.smiles)             # CCC: used for species identity
print(seed.mapped_smiles)      # Retains carbon maps 101, 205, 309; includes mapped H
print(dict(seed.input_maps))   # Assigned labels mapped back to supplied labels
print(seed.symmetry_classes)  # (map_number, symmetry_class) pairs
```

For a multi-reactant step, each occurrence is numbered separately. Labels that conflict with a previous occurrence are reassigned above all reserved labels. Each match's `atom_mapping` records the reactant slot, source atom index, species-level map, supplied map, and reaction-local map. Thus, two copies of a species can react without sharing atom IDs. Structural seed aliases are stored in `network.seed_inputs`; enumeration uses one representative per species identity.

Templates run on explicit-H molecules. Rules that transfer hydrogen must move mapped H atoms and their bonds explicitly. Creating untraceable atoms, dropping source atoms, or producing duplicate product labels rejects that candidate with reason `atom_mapping`. Charge and site balance are checked separately. Existing implicit-H templates may need rewriting; this release does not silently infer their hydrogen correspondence. The bundled version-2 rules specify overall atom correspondence, not validated elementary hydrogen-transfer mechanisms.

`reaction.channels` stores one complete representative per symmetry-distinct mapped reaction. A record includes `mapped_reaction_smiles`, the mapped reactant/product SMILES lists, source numbering in `atom_mapping`, product-atom provenance, `symmetry_key`, and `match_count`. Reactant occurrence order is retained even when stoichiometry combines two identical species. A duplicate increments the existing channel's `match_count`; its full mapping is discarded. There is no separate `reaction.matches` list or list of duplicate-match indices.

For propane dehydrogenation, 24 template matches yield one stored channel with `match_count=24`. The summary field `template_matches` and CSV columns with that name sum the counters; they do not count stored records.

Channel grouping uses a canonical graph of the complete reaction before and after transformation. Component nodes preserve molecular grouping, phases, and sites. Atomic labels preserve elements, isotopes, charges, radicals, and assigned atom CIP labels; bond labels preserve both sides' bond orders and defined stereo. This avoids treating independent membership in atom symmetry classes as proof of joint reaction equivalence. It also permits exchange of chemically identical reactant occurrences.

Basic atom classes use `CanonicalRankAtoms(breakTies=False, includeAtomMaps=False)`, without enumerating self-matches. Explicit-H neopentane has class sizes 1, 4, and 12. True atom maps remain unique and are never replaced by symmetry-class labels. Defined stereochemistry is retained; new stereoisomers are not exhaustively generated.

These are single-step mappings. Numbers are local to a species representative or a reaction occurrence, not globally unique atom identities throughout the network. Multi-step ancestry requires composing the individual mappings and reconciling product atoms with the species registry representatives; map integers alone do not establish that correspondence. Neither channel counts nor `match_count` values are SEMK degeneracy factors. Physical pathway counting and kinetic parameter sharing remain separate work. This version does not calculate rates, conversions, or selectivities.

New JSON exports use schema version 3 and preserve compact channel representatives, counters, numbering, symmetry classes, and seed inputs. Loading a version-2 file converts its raw matches and channel indices to the compact representation. Version-1 files remain readable as structural networks; their incomplete legacy provenance is not retained as mapped channels. Source atom mappings cannot be recovered from an old unmapped structure.

## Command line and outputs

After installation, run:

```bash
chemweave generate --config examples/propane.toml --output results/propane
```

Use `examples/mapped_propane.toml` to run the same network from a partially mapped seed.

A custom `rules = "my_rules.json"` entry is resolved relative to the configuration file's directory. The output directory is resolved relative to the working directory of the command.

| File | Contents |
| --- | --- |
| `network.json` | Complete species, reactions, representative channel mappings and counts, rules, configuration, and RDKit version; can be reloaded. |
| `species.csv` | Species ID, name, unmapped structure, formula, charge, phase, site, and mapped SMILES. |
| `reactions.csv` | Reaction ID, rule, family, step classification, stoichiometry on both sides, template-match count, and channel count. |
| `channels.csv` | One representative mapped reaction and its `match_count` for each symmetry-grouped channel. |
| `stoichiometry.csv` | Species and positive integer coefficients for each side of each reaction. Assign a negative sign to reactant coefficients when constructing a stoichiometric matrix. |
| `network.dot` | Graphviz graph with a separate node for each reaction, preserving multiple-reactant and multiple-product relationships. |

DOT export does not require Graphviz. If Graphviz is installed, render the graph with `dot -Tsvg results/propane/network.dot -o results/propane/network.svg`.

`termination="exhausted"` means there are no pending combinations within the selected rules and carbon-number range. `max_species` or `max_reactions` indicates a budget stop; the saved result is a partial network. The command exits with code 0 on completion and 2 on a budget stop. Input errors are reported by Python/RDKit.

## Further development

Extend the rule set with sources and applicability conditions for cracking, branching, cyclization, and aromatic formation, then compare results at higher carbon-number limits. Add vacant-site and catalyst-state representations based on the selected surface mechanism, followed by single-event parameter organization and kinetic solving.

## References

- [RDKit chemical reaction API](https://www.rdkit.org/docs/GettingStartedInPython.html#chemical-reactions)
- [RDKit Reaction SMARTS](https://www.rdkit.org/docs/RDKit_Book.html#reaction-smarts)
- [Python packaging guide](https://packaging.python.org/en/latest/tutorials/packaging-projects/)
