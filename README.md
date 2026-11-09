**English** | [简体中文](README.zh-CN.md)

# ChemWeave

ChemWeave helps you explore how a reaction network grows from your starting species and reaction templates. Generate and inspect networks, follow product formation pathways, and save the results through a local web interface or Python. It also calculates **reaction-path factors (`n_e`) for single-event kinetics**, supporting gas-phase reactions and abstract catalytic sites.

## Contents

- [Installation and testing](#installation)
- [Browser interface](#web-interface)
- [Python API](#python-api)
- [Custom rules](#custom-templates)
- [Network expansion](#network-expansion)
- [Surface states and single-event modeling](#catalyst-sites)
- [Atom mapping and symmetry](#atom-mapping)
- [Depth and target stopping](#search-limits)
- [Formation pathways](#formation-pathways)
- [Command line and outputs](#cli-and-outputs)
- [Single-event kinetics: calculating n_e](#single-event-kinetics)
  - [Automatic n_e calculation](#automatic-ne)
  - [Combined gas and catalytic example](#combined-example)
  - [Standalone catalytic example: H₂ dissociative adsorption](#catalytic-example)
  - [Earlier symmetry-ratio validation example](#symmetry-example)
- [Progress and user cancellation in Python](#progress-and-stop)
- [Further development](#development)
- [References](#references)

<a id="installation"></a>

## Installation and testing

Requires Python 3.11 or later and RDKit. You can use an existing environment, such as Conda's `base` or `ml` environment. From the directory containing the project:

```bash
conda activate base
cd chemweave
python -m pip install -e .
python -m pytest -q
```

pytest is declared as an optional development dependency in `.[test]`. To install it with the package, use `python -m pip install -e ".[test]"`. An environment that already has pytest can use the commands above.

On macOS/Linux, with Python and RDKit already available in your environment, you can run directly from the project source without installing ChemWeave:

```bash
PYTHONPATH=src python -m chemweave generate \
  --config examples/propane.toml --output results/propane
```

<a id="web-interface"></a>

## Browser interface

Windows, macOS, and Linux use the same Python entry point. Open a terminal in the project directory containing `pyproject.toml`, using a Python 3.11+ environment. Install the optional web dependencies once, then launch the local workspace:

```bash
python -m pip install -e ".[web]"
python -m chemweave web
```

After installation, `chemweave web` is an equivalent command when that environment is active. Keep the terminal open while using the interface. The browser opens at `http://127.0.0.1:8765`; if it does not open automatically, enter this address manually. The interface uses English labels. Add seed SMILES and optional names, select rules, set the carbon/depth limits, and click **Generate**. Optional target rows track experimental products. **Stop when all targets are found** is enabled by default; turn it off to keep searching for alternative formation routes. Expand each species row to set its phase and site.

The center panel displays species structures and separate reaction nodes with stoichiometric coefficients. Click a species or reaction to inspect its details, switch mapping channels, or display the atom-map numbers actually stored in that reaction. Search, depth filtering, and neighborhood views help explore the graph. Networks above 80 species initially display up to depth 2; the depth selector can show the remaining layers. Filtering retains complete visible reactions, including every reactant and coproduct. Statistics update while generation runs; the graph is drawn after completion.

**Stop** requests a cooperative stop between candidates. An active native RDKit matching call must return before that request can take effect. Accepted reactions remain complete; a stopped run has `termination="user_stopped"` and can be saved. Its match counters and alternative pathways may be incomplete.

Use **Save project** to download a JSON file containing inputs, selected rules, settings, and the network. **Open** accepts either this project file or an existing `network.json`, including supported older schemas. **Export** provides network JSON, a ZIP of CSV tables plus DOT, or a PNG of the current graph view. Browser downloads use the browser's normal download location.

The local server keeps one workspace in memory, shared by browser tabs. Refreshing restores that workspace; saving a project preserves it across server restarts. **New** creates a blank workspace after confirmation: no seed species or reaction templates are loaded. **Load example** replaces the current workspace with the gas + catalytic example inputs without running them. **Generate** runs only the current inputs. The server listens only on the local loopback interface. Shut it down with Ctrl+C in its terminal.

To keep dependencies in a separate environment, create it once from the project directory:

```text
python -m venv .venv-web
```

Then install and start using that environment's Python directly; activation is not required.

**Windows (PowerShell or Command Prompt):**

```text
.\.venv-web\Scripts\python.exe -m pip install -e ".[web]"
.\.venv-web\Scripts\python.exe -m chemweave web
```

**macOS / Linux:**

```bash
.venv-web/bin/python -m pip install -e ".[web]"
.venv-web/bin/python -m chemweave web
```

For later sessions, run only the launch command for your environment. If `python` is named `python3` on your system, use `python3` when creating the environment. An existing Conda environment can also run the common installation and launch commands above. To reuse packages such as RDKit from the Python environment that creates the virtual environment, add `--system-site-packages` to the `venv` command.

Use `python -m chemweave web --port 8766` to choose another port, or add `--no-browser` to start without opening a browser automatically. When using the separate environment, replace `python` with its platform-specific path shown above. Close the service with Ctrl+C in the terminal.

The backend uses [FastAPI](https://fastapi.tiangolo.com/) and Uvicorn. RDKit draws molecule SVGs. [Cytoscape.js](https://js.cytoscape.org/) 3.34.4 and its MIT license are bundled in the package; the UI does not use a CDN at runtime. No Node.js or frontend build step is needed. Library and command-line generation remain usable without installing the web extra.

<a id="python-api"></a>

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

<a id="custom-templates"></a>

## Custom rules

Rules are stored as a JSON list. This example represents a mapped propane dehydrogenation step, including H2. Supplied steps are treated as elementary inputs; the program does not establish elementary character:

```json
[
  {
    "id": "propane_dehydrogenation",
    "smarts": "[C;H3:1](-[H:4])-[C;H2:2](-[H:5])-[C;H3:3]>>[C:1]=[C:2]-[C:3].[H:4]-[H:5]",
    "family": "dehydrogenation",
    "kind": "elementary",
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
| `kind` | Optional legacy metadata, defaulting to `elementary`. All supplied steps are treated as elementary inputs; this field does not gate factor calculation. |
| `version`, `source`, `notes` | Rule version, supporting evidence, and assumptions. |
| `exclude_smarts` | Skip the rule if any reactant matches any of these patterns. |
| `reactant_phases`, `product_phases` | Species phases in template order; all default to `gas`. |
| `reactant_sites`, `product_sites` | Legacy species-level labels in template order; all default to `null`. |
| `site_constraints` | Required attributes keyed by template atom-map number; explicit catalyst queries use `[#0]`. |

Rules are not automatically reversed. Supply a template for the reverse reaction when needed; some double-bond migration templates can themselves generate both directions. Setting `kind` to `elementary` does not mean the program has established the mechanism.

<a id="network-expansion"></a>

## Network expansion

1. Standardize seed structures and identify species by structure, phase, legacy site label, and explicit catalyst-node attributes.
2. Match each new species against each reactant slot of each rule once.
3. Combine newly added candidates with existing candidates and execute ordered reactant tuples that have not been processed. A species may occupy two slots to represent self-coupling.
4. Execute templates with explicit hydrogens. Restore source atom maps and check a one-to-one atom inventory (including elements and isotopes), charge balance, legacy site-label counts, and explicit catalyst-node conservation/capacity.
5. Apply `max_carbon`, then merge duplicate reactions by rule ID and stoichiometry. Store one complete representative per symmetry-distinct mapping in `channels`; equivalent executions increment its `match_count` without retaining duplicate records.
6. Queue new species until the queue is exhausted or the species/reaction budget is reached.

Checks occur when inputs are read and candidate reactions are generated; the full network is not repeatedly revalidated during enumeration. Molecular descriptors and reaction templates are prepared when their objects are created. Enumeration does not estimate barriers or fluxes, and species are not removed merely because they have not been detected experimentally.

`max_carbon` applies to each species, including the seeds. It defines the scope of a run; it does not establish that pathways outside that scope are unimportant in the actual system.

<a id="catalyst-sites"></a>

## Surface states and single-event modeling

Species can contain explicit mapped catalyst nodes. In species SMILES, `[*:1]` is an abstract node; in reaction SMARTS use `[#0:1]` to avoid the any-atom wildcard. Supply its type and meaning separately:

```python
metal = Species.from_smiles("[*:1]", phase="surface",
    site_atoms={1: {"kind": "metal", "type": "Fe_model", "capacity": 1}})
acid = Species.from_smiles("[*:1]O", phase="surface",
    site_atoms={1: {"kind": "bronsted", "type": "ZSM5_model", "capacity": 1}})
```

An abstract node can be vacant or attached to one atom with bond order 1–4. Missing attributes default to kind `surface`, type `X`, capacity 1. Multiple nodes represent multi-site occupancy. The Brønsted model contains a real O–H bond and transferable mapped H; protonation forms a surface-associated framework O−/carbenium ion pair, retained as one species. Attributes, attachments, isotope and charge participate in identity, matching, mapping, conservation and persistence. Lattice adjacency and adsorption orientation are outside this abstract graph model. See [the counting and catalyst-model report](docs/path_degeneracy.md) for the precise assumptions and source-linked tests.

Legacy `Species.from_smiles("C=CC", phase="surface", site="BAS")` remains usable for generic labelled networks, but a species-level label is insufficient for automatic catalytic path counting. In the browser, expand **Phase and catalyst sites** to enter the node metadata JSON. Save JSON or the accompanying `site_atoms` data: compact SMILES alone cannot preserve catalyst types.

<a id="atom-mapping"></a>

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

Templates run on explicit-H molecules. Rules that transfer hydrogen must move mapped H atoms and their bonds explicitly. Creating untraceable atoms, dropping source atoms, or producing duplicate product labels rejects that candidate with reason `atom_mapping`. Charge and site balance are checked separately. Existing implicit-H templates may need rewriting; this release does not silently infer their hydrogen correspondence. The bundled version-2 rules provide complete atom correspondence. As with user templates, elementary character is an input assumption, not a software validation result.

`reaction.channels` stores one complete representative per symmetry-distinct mapped reaction. A record includes `mapped_reaction_smiles`, the mapped reactant/product SMILES lists, source numbering in `atom_mapping`, product-atom provenance, `symmetry_key`, and `match_count`. Reactant occurrence order is retained even when stoichiometry combines two identical species. A duplicate increments the existing channel's `match_count`; its full mapping is discarded. There is no separate `reaction.matches` list or list of duplicate-match indices.

For propane dehydrogenation, 24 template matches yield one stored channel with `match_count=24`. The summary field `template_matches` and CSV columns with that name sum the counters; they do not count stored records.

Channel grouping uses a canonical graph of the complete reaction before and after transformation. Component nodes preserve molecular grouping, phases, and sites. Atomic labels preserve elements, isotopes, charges, radicals, and assigned atom CIP labels; bond labels preserve both sides' bond orders and defined stereo. This avoids treating independent membership in atom symmetry classes as proof of joint reaction equivalence. It also permits exchange of chemically identical reactant occurrences.

Basic atom classes use `CanonicalRankAtoms(breakTies=False, includeAtomMaps=False)`, without enumerating self-matches. Explicit-H neopentane has class sizes 1, 4, and 12. True atom maps remain unique and are never replaced by symmetry-class labels. Defined stereochemistry is retained; new stereoisomers are not exhaustively generated.

These are single-step mappings. Numbers are local to a species representative or a reaction occurrence, not globally unique atom identities throughout the network. Multi-step ancestry requires composing the individual mappings and reconciling product atoms with the species registry representatives; map integers alone do not establish that correspondence. Neither channel counts nor `match_count` values are SEMK degeneracy factors. Independent factors and kinetic parameter sharing are handled by the optional single-event annotation module. The package does not fit parameters or solve rates, conversions, or selectivities.

New JSON exports use schema version 6 and preserve compact channel representatives, counters, numbering, symmetry classes, seed inputs, species/reaction depths, and target status. Version-3, version-4 and version-5 files remain readable; depths absent from older files are left unknown. Loading a version-2 file converts its raw matches and channel indices to the compact representation. Version-1 files remain readable as structural networks; their incomplete legacy provenance is not retained as mapped channels. Source atom mappings cannot be recovered from an old unmapped structure.

<a id="search-limits"></a>

## Depth and target stopping

Set `GenerationConfig(max_depth=2)` to allow at most two successive reaction steps from the seeds. Seeds have depth 0. For each reaction:

```text
reaction depth = 1 + max(depth of every reactant)
new product depth = reaction depth
```

A species keeps its minimum reachable depth. FIFO processing discovers species in increasing depth order, including when several seeds or bimolecular reactions are present. Cycles do not increase a registered species's depth. For a bimolecular step, this measures the longest dependency chain; parallel preparation branches are not added together. A reaction returning to an earlier species still has its own reaction depth.

Products at the depth limit are retained, along with the reaction and mappings that generated them. They are not expanded further. `max_depth=0` keeps only the seeds; omitting it leaves depth unrestricted. An exhausted network may stop before reaching the configured depth. With the bundled rules, propane is at depth 0, propene and H2 at depth 1, and propene self-coupling produces a C6 olefin at depth 2. Try `examples/depth_limited.toml`.

Pass specific target species to stop once **all** have appeared:

```python
network = Generator(
    load_rules(),
    GenerationConfig(max_carbon=6, max_depth=4),
).generate(
    [Species.from_smiles("CCC")],
    targets=[Species.from_smiles("C=CC"), Species.from_smiles("C=C")],
)
print(network.summary())
print(network.species_depths)  # species ID -> minimum depth
print(network.targets)        # requested identities and discovered species IDs
```

Targets use the same identity as the registry: canonical structure, phase, and site, ignoring atom-map numbers and names. Duplicate targets count once; targets already present among seeds count immediately. Specifying targets does not add them to the network. This is a structural reachability condition, independent of experimental concentrations or detection limits. `max_species` continues to limit the total number of stored species.

In TOML, use `[[targets]]` entries with the same fields as `[[seeds]]`:

```toml
[[seeds]]
smiles = "CCC"

[[targets]]
smiles = "C=CC"

[[targets]]
smiles = "C=C"

[generation]
max_carbon = 6
max_depth = 4
```

`examples/target_species.toml` contains this configuration with labels and size budgets. The two conditions can be combined: stop at the depth boundary, or earlier when all targets appear. The summary reports `depth_reached`, `targets_total`, `targets_found`, and `missing_targets`. Each target in `network.targets` has its discovered `species_id`, or `None` when missing.

Target stopping saves the complete accepted reaction, all its coproducts, and the channel mapping that satisfied the condition, then returns immediately. Unvisited matches, channels, and alternative pathways are not enumerated, so match counters can be incomplete. The resulting partial network depends on seed/rule enumeration order and does not establish that every shortest pathway has been collected.

<a id="formation-pathways"></a>

## Formation pathways

After generating or opening a network, use **Formation pathways** in the left panel. Select a product, set **Max reactions per route** and **Max pathways**, then click **Analyze formation**. You can also start from a species inspector's **Analyze formation** button. Turn off **Stop when all targets are found** before generating the network if you want to collect alternative routes. Analysis uses the stored network; it does not expand a previously truncated run.

Choose one pathway to see its reactions in an executable order, or **Union of found pathways** to compare routes. Each reaction keeps all reactants, coproducts, coefficients and mapping channels. The union is a comparison view, not a new executable route. Pathway views bypass the network depth display filter so complete reactions remain visible; **Reset view** returns to the full network. **Stop** also works during analysis and retains routes already found. **Save project** includes the latest analysis, and **Export pathways JSON** saves the analysis with the underlying network.

```python
from chemweave import GenerationConfig, Generator, Species, find_pathways, load_rules

network = Generator(load_rules(), GenerationConfig(
    max_carbon=6, max_depth=4, stop_when_targets_found=False,
)).generate([Species.from_smiles("CCC")], targets=[Species.from_smiles("C=C")])

target_id = network.targets[0]["species_id"]
if target_id is not None:
    result = find_pathways(network, target_id, max_reactions=8, max_pathways=100)
    for route in result["pathways"]:
        print(route["reactions"])  # Reaction IDs in one valid execution order
    print(result["termination"], result["network_search_complete"])
```

This example yields two four-reaction routes to ethene with the bundled rules. `examples/formation_pathways.toml` generates the same network from the CLI. To analyze a saved network, replace `S00008` below with the target's ID in that file:

```bash
python -m chemweave generate --config examples/formation_pathways.toml --output results
python -m chemweave pathways --network results/network.json --target S00008 --max-reactions 8 --max-pathways 100 --output results/pathways.json
```

`stop_when_targets_found = false` belongs in the TOML `[generation]` section. Omitting it preserves the existing early-stop behavior. The `pathways` command exits with code 2 if enumeration hits its pathway cap or the stored network was stopped early; it still writes the result. An exhaustive analysis with no matching routes exits 0 and reports an empty list.

The analyzer first traces candidate producers backward, then enumerates seed-supported reaction sets forward. All reactants of a reaction must be available. A route is retained only if removing any one reaction makes the target unreachable within that route. This removes unnecessary preparation branches and product-recycling detours, while retaining longer alternative routes. Shared preparation reactions count once, and coproducts may support other branches. A structural cycle is not automatically discarded when its steps produce necessary coproducts; what is excluded is circular dependence on material that has no seed-supported source. Independent execution orders of the same reaction set count as one route. A target already among the seeds has one zero-reaction route.

Limits and interpretation:

- `max_reactions` counts distinct reactions across **all** branches. For `A → B`, `A → C`, `B + C → P`, the route contains three reactions although its dependency depth is two. The stored network's minimum depths do not bound every alternative route's length.
- Results describe qualitative species availability, not finite inventories, repeated reaction events, net yields or kinetics. For example, self-coupling requires a source of that species but does not check that two molecules remain in an inventory. Stoichiometric coefficients are retained for inspection.
- Mapping channels stay on the original reaction records. Routes are not expanded into channel combinations or assigned a continuous atom map across steps.
- `termination="exhausted"` means enumeration finished within the **stored network and the chosen route-size bound**. `max_pathways` means the cap was reached and completeness was not confirmed; `user_stopped` means cancellation. Traversal follows stored reaction order, not a shortest-path or kinetic ranking.
- `network_search_complete` is true only for network termination `exhausted` or `max_depth`, relative to its rules and configured bounds. Target stopping, species/reaction budgets, and user cancellation can omit alternative routes even when pathway enumeration itself finishes. The web interface reports these two sources of incompleteness separately.

<a id="cli-and-outputs"></a>

## Command line and outputs

After installation, run:

```bash
chemweave generate --config examples/propane.toml --output results/propane
```

Use `examples/mapped_propane.toml` to run the same network from a partially mapped seed.

A custom `rules = "my_rules.json"` entry is resolved relative to the configuration file's directory. The output directory is resolved relative to the working directory of the command.

| File | Contents |
| --- | --- |
| `network.json` | Complete species, reactions, depths, target status, representative channel mappings and counts, optional single-event annotations/parameters, rules, configuration, and RDKit version; can be reloaded. |
| `species.csv` | Species ID, name, unmapped structure, formula, charge, phase, site, mapped SMILES, and minimum depth. |
| `reactions.csv` | Reaction ID, rule, family, step classification, stoichiometry on both sides, template-match count, channel count, and reaction depth. |
| `channels.csv` | One representative mapped reaction and its `match_count` for each symmetry-grouped channel. |
| `stoichiometry.csv` | Species and positive integer coefficients for each side of each reaction. Assign a negative sign to reactant coefficients when constructing a stoichiometric matrix. |
| `network.dot` | Graphviz graph with a separate node for each reaction, preserving multiple-reactant and multiple-product relationships. |

DOT export does not require Graphviz. If Graphviz is installed, render the graph with `dot -Tsvg results/propane/network.dot -o results/propane/network.svg`.

The saved `termination` states why generation stopped:

| Value | Meaning |
| --- | --- |
| `exhausted` | No pending combinations remain within the selected rules and carbon-number range. Targets may still be missing. |
| `max_depth` | The configured depth boundary was reached; products at that depth are included. |
| `targets_reached` | All requested target species are present, possibly already among the seeds. |
| `user_stopped` | A caller or the browser requested a cooperative stop. |
| `max_species`, `max_reactions` | A size budget prevented adding the next complete reaction; a consistent partial network is saved. |

The command exits with code 0 for exhaustion, depth completion, or target success **provided no requested targets are missing**. A size-budget stop, a user stop, or missing targets returns 2; outputs are still saved. Input errors are reported by Python/RDKit.

<a id="single-event-kinetics"></a>

## Single-event kinetics: calculating n_e

SEMK is one use of ChemWeave. `Generator` continues to enumerate general rule-based networks without kinetic assumptions. Call `annotate_single_events(network, ...)` afterwards, add a `[single_event]` TOML section for CLI generation, or enable **Single-event annotation → Annotate after generation** in the English browser interface. Omitting the option leaves networks unannotated. All supplied mapped steps, including the default demo rules, are eligible for path-factor calculation. Historical `kind="overall"` labels in saved projects do not block calculation. Reapply annotation to replace old skipped results; regeneration is not required.

`single_event.py` is family-independent: it extracts reaction centers, dispatches an explicitly selected calculator, accepts imported factor assignments, and manages shared parameter records. No calculator is selected by default: `annotate_single_events(network)` extracts centers without assuming a reaction family or assigning a factor. Automatic constitutional-path counting lives in `path_degeneracy.py`, with optional example templates in `data/path_families.json` that are not used by the counting algorithm. The earlier methane-specific symmetry example remains separate in `single_event_calculators.py`.

Alternative calculation methods can use `register_single_event_calculator(name, callback)` and then `annotate_single_events(network, calculator=name)`, without changing the core module. A callback receives keyword arguments `network`, `reaction`, `channel`, `center`; it returns `None` when inapplicable, or a dictionary containing `convention`, `source`, and either `factor` or `sigma_reactants`/`sigma_ts`. Optional `references`, `kinetic_subclass` and `parameter_class` are retained; these are not factor prerequisites. The callback owns the chemical and phase/site applicability checks; the core checks and stores the common factor/parameter records. Imported records override calculator output for the selected channels. Custom callbacks must be registered again before replaying their saved options in a new process; executable Python code is not serialized.

<a id="automatic-ne"></a>

### Automatic n_e calculation

Select `calculator="mapped_event_v1"` (`path_degeneracy_v1` is a compatibility alias). The calculator extracts bond changes, atom-state changes and product species partitions directly from a complete mapped elementary reaction. It enumerates the event's orbit under reactant-graph symmetry. New templates need no second recipe and no built-in family name. The saved method is **`reaction_path_degeneracy`**, algorithm **`mapped_event_orbit_v1`**; this is a structural path factor under the documented counting convention, not a three-dimensional rotational-symmetry calculation.

```python
from chemweave import GenerationConfig, Generator, ReactionRule, Species, annotate_single_events

rule = ReactionRule(id="user_CH_homolysis",
                    smarts="[C;X4:1]-[#1:2]>>[C;H0:1].[#1;H0:2]")
network = Generator([rule], GenerationConfig(max_depth=1)).generate([Species.from_smiles("C")])
annotate_single_events(network, calculator="mapped_event_v1")
assert network.reactions[0].channels[0]["single_event"]["n_e"] == 4
```

The pure `calculate_path_factor(mapped_reactants, mapped_products, ...)` interface accepts lists of fully mapped, explicit-H species SMILES. Each list entry denotes one species occurrence, including associated disconnected species. The generator supplies these full mappings automatically from a mapped template. Optional phase/site states and explicit site attributes preserve catalyst distinctions. Parameter classes, shared parameters, `rate_basis` and `A_units` are not required; pure calculation creates no parameter classes or rate law.

**Counting:** keep reactant occurrences distinguishable during event enumeration. Deduplicate identical events on those labelled atoms to obtain `N = unique_labeled_paths`; this does not quotient identical-reactant exchange. Apply `D = product(m_i!)` once, where `m_i` counts identical reactant species occurrences: **`n_e = N / D`**. H2 + 2 X → 2 XH gives `2 / 2 = 1`; the reverse gives `1 / 2 = 0.5`. `match_count`, product identity alone and the number of saved channels never define the factor. The saved counting evidence makes the single correction explicit. Users remain responsible for using this factor consistently in their kinetic expressions, without applying a second identical-reactant divisor.

Run the new unlisted-template example:

```bash
PYTHONPATH=src python -m chemweave generate \
  --config examples/custom_ne.toml --output results/custom-ne
```

In the browser, expand **Custom templates (JSON)**, enter a list of templates, and click **Use these templates**. Set seed SMILES (and phase/site attributes when relevant), enable **Annotate after generation**, leave `{"calculator":"mapped_event_v1"}` in the options, and click **Generate**. Select a reaction to see its factor, status, event changes and count/divisor. **Apply to current network** also annotates saved mapped networks. Gas/metal/acid examples are conveniences, not the only supported inputs.

There is no C/H, acyclic or family whitelist. Validation covers homolysis, abstraction, disproportionation/recombination, metal-site adsorption/desorption/hydrogenation and acid-site proton transfer, hydride shift and β-scission, plus oxygen and ring inputs. The algorithm counts the supplied electronic graph, not a resonance ensemble. Assigned stereochemistry currently returns `unsupported_algorithm`; incomplete maps or absent catalyst nodes return `missing_information`; inconsistent input returns `invalid_input`; exhausted computational limits return `incomplete`. These states carry a reason and no factor. Manual factor imports remain available. Abstract sites do not imply lattice, adsorption-orientation or three-dimensional multiplicities.

[The implementation and validation report](docs/path_degeneracy.md) records the counting proof, 11 document-linked or explicitly derived cases, manual enumeration, independent brute-force checks and fixed RMG source references. RMG itself was **not run** locally. This route does not use the 2023 hydrocracking free-state cancellation assumption or claim equality with all Vandewiele rotational factors: methyl/ethyl disproportionation yields path factor 3, versus 6 under that paper's symmetry convention.

<a id="combined-example"></a>

### Combined gas and catalytic example

Click **Load example** to load all five templates together, with propane and a vacant metal-site seed; optional factor annotation is enabled in this demo. A new workspace starts blank. Propane dehydrogenation produces H₂ at depth 1, which participates in H₂ + 2 M → 2 MH at depth 2. The adsorption template is selected alongside the gas templates and can be unchecked independently. Existing saved projects retain their own templates and seeds; choose **Load example** to load the combined example, then **Generate** to run it.

`examples/combined_gas_catalytic.toml` runs the same combined setup:

```bash
PYTHONPATH=src python -m chemweave generate \
  --config examples/combined_gas_catalytic.toml --output results/combined-gas-catalytic
```

The coupling between these demo steps is the generated H₂ species; this does not introduce a metal-mediated propane dehydrogenation mechanism.

<a id="catalytic-example"></a>

### Standalone catalytic example: H₂ dissociative adsorption

Alongside the custom gas template, `examples/custom_catalytic_rules.json` supplies a user-defined catalytic template with no family registration or additional recipe. `examples/custom_catalytic_ne.toml` supplies hydrogen and an explicit vacant metal site species:

```bash
PYTHONPATH=src python -m chemweave generate \
  --config examples/custom_catalytic_ne.toml --output results/custom-catalytic-ne
```

Source: Lozano-Blanco et al. (2011), p.427 Table 2 (1), H₂ + 2 M → 2 MH. The node `[*:1]` has `phase="surface"` and attributes `kind="metal", type="M", capacity=1`. One seed site species supplies two independent identical reactant occurrences. The template matches sites with `[#0]`, breaks H1–H2 and forms H1–M3 and H2–M4. The distinct labelled assignments are `(H1–M3, H2–M4)` and `(H1–M4, H2–M3)`; simultaneous swapping generates a duplicate assignment, not a third path. Thus `N=2`, `D=2!`, **n_e=1**. Sites and their properties persist into the adsorbates. No lattice-neighbor or orientation multiplicity is assumed.

In the English web interface, **Single-event annotation → H₂ adsorption on metal → Generate** loads the same reaction as a built-in convenience. To use the standalone custom example, import its rules JSON and enter the seeds above, or **Open** the CLI-generated `results/custom-catalytic-ne/network.json`. The custom template has its own ID and omits a family label; it uses the same generic calculation. The reaction inspector and CSV record both the two labelled paths and the divisor 2.

<a id="symmetry-example"></a>

### Earlier symmetry-ratio validation example

The earlier symmetry validation family is forward gas-phase **H/D atom + methane(H/D) → dihydrogen + methyl radical**. Its elementary rule is distributed separately in `src/chemweave/data/single_event_rules.json`; **Load methane validation example** loads it in the browser without changing the default rules. Run the same example from the CLI:

```bash
PYTHONPATH=src python -m chemweave generate \
  --config examples/single_event.toml --output results/single-event
```

The calculator `methane_h_abstraction_v1` adopts a tetrahedral methane reactant and an axial C–H/D–H/D saddle with three spectator ligands. **Proper rotational symmetry numbers** define `n_SE = product(sigma_reactants) / sigma_TS`. They are not full graph automorphism counts. The atomic reactant has symmetry number 1. This convention excludes product symmetry, electronic/nuclear spin weights, tunneling, and recrossing; these are not inferred from atom mapping. There are no internal rotors in this validation family.

| Methane | Transferred isotope | Methane sigma | TS sigma | n_SE |
| --- | --- | --- | --- | --- |
| CH4 | H | 12 | 3 | 4 |
| CH3D | H / D | 3 | 1 / 3 | 3 / 1 |
| CH2D2 | H / D | 2 | 1 / 1 | 2 / 2 |
| CHD3 | H / D | 3 | 3 / 1 | 1 / 3 |
| CD4 | D | 12 | 3 | 4 |

These H/D results follow from the stated idealized geometry and isotope-preserving proper rotations: Td has 12 rotations, C3v has 3, C2v has 2, and Cs has 1. [NIST's symmetry-number explanation](https://cccbdb.nist.gov/thermox.asp) documents the convention; [NIST CCCBDB reaction 8](https://cccbdb.nist.gov/ts2x.asp?ri=8) identifies the H + CH4 elementary reaction. [Hu et al., J. Phys. Chem. A 2006](https://web.stanford.edu/group/Zarelab/publinks/764.pdf) describes its C3v saddle. The isotope variants validate this symmetry convention, not measured kinetic isotope effects. No fitted or ab initio rate constants are bundled.

All calculators retain `single_event.reaction_center`: changed bonds (including order/stereo), changed atom charge/radical/stereo states, and the participating **reaction-local atom maps**. The methane symmetry calculator recognizes the complete bond/radical transformation and gas/no-site state, rather than trusting a rule name. Use consistent protium labels: methane mixing unlabeled H with explicit 1H is unsupported because those labels split RDKit channels. Its kinetic subclass includes carbon isotope, transferred isotope, abstractor isotope, and spectator isotope composition. Distinct isotope environments therefore do not automatically share fitted parameters. Inputs outside this older methane calculator’s scope retain a center and an explicit `unsupported` status with no factor. Historical step-classification labels do not restrict the generic path calculator. Reactions without mapping channels are listed in `network.single_event["unmapped_reactions"]`.

**`match_count` and channel count are never single-event factors.** Validation tests deliberately change the template to produce 4 versus 12 matches while the independently calculated factor remains 4. Early-stop counters likewise do not change the symmetry ratio. This older symmetry calculator covers only the forward methane family. The separate path calculator has the scope and normalization described above; it does not extend the symmetry-ratio method.

```python
from chemweave import Network, annotate_single_events

network = Network.from_json("results/single-event/network.json")
annotate_single_events(network, calculator="methane_h_abstraction_v1")
channel = network.reactions[0].channels[0]
print(channel["single_event"]["factor"])  # 4.0, from 12 / 3
class_id = channel["single_event"]["parameter_class"]
print(network.single_event["parameters"][class_id])
# A, n, Ea start as None; status is "unfitted".
network.export_json("results/single-event/annotated.json")
network.export_csv("results/single-event")
```

The earlier calculator and explicitly supplied parameter records retain optional shared-parameter support. A channel may reference `parameter_class`, with values stored once in `network.single_event["parameters"]`. The legacy modified Arrhenius expression is retained only when parameter classes are present. The generic path calculator neither creates these records nor requires a rate expression, units or fitted values; no fitting or rate evaluation is performed.

Factors can also be imported as positive finite numbers, or as reactant/TS symmetry numbers whose ratio defines the factor. An import requires `convention` and `source`; subclass, shared class, `rate_basis` and `A_units` are optional. If factor and symmetry numbers are both supplied, they must agree. The following legacy methane input illustrates optional shared records, **not experimental fit data**:

```json
{
  "calculator": null,
  "factors": [{
    "reaction_id": "R00001", "channel": 1,
    "kinetic_subclass": "methane-H-abstraction",
    "parameter_class": "methane-H-fit",
    "sigma_reactants": [12, 1], "sigma_ts": 3,
    "convention": "Forward axial saddle; proper rotations; n_SE = 12/3",
    "source": "Independent symmetry assignment under the documented methane convention"
  }],
  "parameters": {
    "methane-H-fit": {"A": null, "n": null, "Ea": null, "status": "unfitted", "source": "Awaiting experimental fitting"}
  }
}
```

Save these options as `single-event.json`, then annotate a saved network:

```bash
python -m chemweave annotate --network results/single-event/network.json \
  --annotations single-event.json --output results/annotated
```

The same keys are accepted by `annotate_single_events(network, **options)` and the browser's **Calculation options and imported factors (JSON)** editor. **Apply to current network** updates the loaded network without regeneration. Channel indices are 1-based and identify records in that specific saved network; recheck selectors after regenerating. A shared class cannot mix molecularities, rate bases or units. The older methane calculator assigns parameter classes; the generic path calculator does not. Applying options replaces previous annotations and parameter records, so include any imported factors and parameter values you want to retain. `calculator: null` (also the default when omitted) disables automatic factor calculation while preserving center extraction. The browser validation button explicitly fills in the methane calculator; an empty options object selects no family. Running `annotate` without an options file likewise extracts centers only.

Network JSON schema 6 and saved browser projects retain centers, subclass/class references, independent factors, conventions, provenance, shared parameter records, and reproducible annotation options. Schema 1–5 remains readable; old records without complete mapping remain unannotated until usable mappings are available. CSV export adds **`single_event.csv`** (one annotation per mapped channel, including unsupported status and center) and **`kinetic_parameters.csv`** (one record per shared class) only when annotation has been applied. `species.csv` appends `site_atoms`; annotation exports append `n_e`, method, algorithm, family, status/reason, extracted event, kinetic-channel identity, labelled path count and identical-reactant divisor. Existing columns are retained. In the English reaction inspector, select a channel to view its center, independent factor, symmetry/TS convention, subclass, shared class, and fitted/unfitted values. Save/Open and JSON/CSV downloads preserve the same data.

<a id="progress-and-stop"></a>

## Progress and user cancellation in Python

`Generator.generate` accepts optional `progress(network)` and `stop_requested()` callbacks. Progress is reported after seeds are registered and after a complete template/reactant application; the callback sees the current network and should return promptly. When the stop callback returns true, generation returns a consistent partial network with `termination="user_stopped"`. Existing calls need no changes. `Network.to_dict()` and `Network.from_dict()` share the JSON representation used by file persistence.

<a id="development"></a>

## Further development

Further development can improve template authoring, comparison of alternative rule sets, formation-pathway exploration, and performance on larger networks. The optional path-factor algorithm can extend its stereochemistry and resonance handling. Experimental fitting and kinetic solving can consume exported networks as separate applications.

<a id="references"></a>

## References

- [RDKit chemical reaction API](https://www.rdkit.org/docs/GettingStartedInPython.html#chemical-reactions)
- [RDKit Reaction SMARTS](https://www.rdkit.org/docs/RDKit_Book.html#reaction-smarts)
- [Python packaging guide](https://packaging.python.org/en/latest/tutorials/packaging-projects/)
