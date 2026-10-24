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

The C6 example generated 22 species and 45 reactions in the tested `base` environment (Python 3.12, RDKit 2026.03.6). These counts illustrate a run; they do not establish chemical validity.

## Custom rules

Rules are stored as a JSON list. This example represents the overall transformation of propane dehydrogenation, including H2:

```json
[
  {
    "id": "propane_dehydrogenation",
    "smarts": "[CH3:1][CH2:2][CH3:3]>>[CH2:1]=[CH:2][CH3:3].[H][H]",
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
4. Standardize products and check valence, element balance (including hydrogen), charge balance, and site-label counts.
5. Apply `max_carbon`, then merge duplicate reactions by rule ID, reactant stoichiometry, and product stoichiometry. Equivalent template matches are retained in `matches`.
6. Queue new species until the queue is exhausted or the species/reaction budget is reached.

Checks occur when inputs are read and candidate reactions are generated; the full network is not repeatedly revalidated during enumeration. Molecular descriptors and reaction templates are prepared when their objects are created. Enumeration does not estimate barriers or fluxes, and species are not removed merely because they have not been detected experimentally.

`max_carbon` applies to each species, including the seeds. It defines the scope of a run; it does not establish that pathways outside that scope are unimportant in the actual system.

## Surface states and single-event modeling

You can create `Species.from_smiles("C=CC", phase="surface", site="BAS")`. Site labels contribute to species identity. Rules can constrain input and output sites, and the generator checks that the number of each site label is conserved.

This is an **abstract single-site representation**, suitable for transformations between explicitly defined molecular structures on the same site. Dedicated representations and bundled rules for vacant sites, proton exchange, two-site species, and Ga catalytic cycles are not yet included. Surface structures, hydrogen, and charge must be defined consistently when extending the mechanism; attaching a label to a gas-phase molecule alone does not establish an adsorption mechanism.

Each reaction retains its rule and template-match records. These include reactant order, product order, atom-index-labeled reactant SMILES, and RDKit's raw product-atom provenance. Map numbers in the indexed SMILES equal the input atom index plus one; `react_atom_idx` in the provenance is zero-based. Atoms introduced by the template, such as those in H2, may have no input-atom provenance.

These records support structural traceability. **`len(reaction.matches)` is not a single-event degeneracy factor.** Rotational symmetry numbers, physical counting of equivalent pathways, transition states, and kinetic parameter sharing require separate implementation. This version does not calculate rates, conversions, or selectivities.

## Command line and outputs

After installation, run:

```bash
chemweave generate --config examples/propane.toml --output results/propane
```

A custom `rules = "my_rules.json"` entry is resolved relative to the configuration file's directory. The output directory is resolved relative to the working directory of the command.

| File | Contents |
| --- | --- |
| `network.json` | Complete species, reactions, template matches, rules, configuration, and RDKit version; can be reloaded. |
| `species.csv` | Species ID, name, structure, molecular formula, charge, phase, and site. |
| `reactions.csv` | Reaction ID, rule, family, step classification, stoichiometry on both sides, and template-match count. |
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
