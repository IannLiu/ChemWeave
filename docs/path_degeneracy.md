# Generic mapped-event path factors

Updated 2026-11-09. The current calculator is `mapped_event_v1`, with output method `reaction_path_degeneracy` and algorithm `mapped_event_orbit_v1`. `path_degeneracy_v1` is retained as a compatibility selection name for this new calculator. Previously saved annotations remain unchanged until explicitly recalculated.

**The input is a complete mapped elementary reaction, not a second family recipe.** The user supplies the elementary mechanism; the software does not establish elementary character, estimate barriers, solve a reactor model or fit experimental data. Family names are optional provenance. The bundled templates are examples and validation inputs only and are never consulted by the counting engine.

## Interfaces

```python
from chemweave import calculate_path_factor

result = calculate_path_factor(
    ["[H:1][H:2]", "[*:3]", "[*:4]"],
    ["[H:1][*:3]", "[H:2][*:4]"],
)
assert result["counting"]["unique_labeled_paths"] == 2
assert result["counting"]["identical_reactant_divisor"] == 2
assert result["n_e"] == 1
```

Each list entry is **one species occurrence**, possibly an associated disconnected species. This preserves the distinction between two separate vacant sites and one species containing two sites. All atoms including H have unique nonzero maps on each side; maps match bijectively across sides. The normal generator already supplies these complete explicit-H channels, so users writing templates do not need to manually number all spectator H atoms in each molecule.

Optional arguments:

- `reactant_states`, `product_states`: lists of `(phase, legacy_site_label)` pairs. Omission infers `surface` for species containing an explicit site node, otherwise `gas`. Supply states when different phase assignments are intended.
- `site_atoms`: catalyst attributes keyed by reaction-local map number; omitted attributes mean abstract `X`.
- `max_paths=10000`, `max_checks=50000`: limits on distinct events and individualized graph comparisons. Exceeding a limit returns `incomplete`, never a partial factor.

No network, family, recipe, parameter class, shared parameter, `rate_basis` or `A_units` is required by this pure interface.

For a network:

```python
from chemweave import GenerationConfig, Generator, ReactionRule, Species, annotate_single_events

rule = ReactionRule(id="user_CH_homolysis",
                    smarts="[C;X4:1]-[#1:2]>>[C;H0:1].[#1;H0:2]")
network = Generator([rule], GenerationConfig(max_depth=1)).generate([Species.from_smiles("C")])
annotate_single_events(network, calculator="mapped_event_v1")
assert network.reactions[0].channels[0]["single_event"]["n_e"] == 4
assert network.single_event["parameters"] == {}
```

All supplied mapped steps are treated as elementary inputs; the software does not prove elementary character. New `ReactionRule` inputs and bundled demos use `kind="elementary"`. The `kind` field is retained as compatibility metadata only: legacy `overall` labels never block calculation or manual imports. Reapply annotation to an old saved project to replace previously skipped results without regeneration. Omitting annotation leaves the generic generator unchanged. The earlier methane rotational-symmetry calculator and optional imported parameter records remain supported separately.

CLI: `PYTHONPATH=src python -m chemweave generate --config examples/custom_ne.toml --output results/custom-ne`. This uses the standalone user template in `examples/custom_ne_rules.json`; no corresponding factor recipe exists.

## Event extraction, counting and normalization

Let `R` be the full reactant graph and `Delta` the mapped event. Graph colors preserve element, isotope, formal charge, radical count, aromaticity, bond order, molecular occurrence membership, phase, legacy site label and explicit catalyst attributes. Map numbers identify atoms but do not act as symmetry colors.

The event contains atom-state changes, bond changes **and the complete product species partition/state**. Therefore association of disconnected ions and adsorption-state changes are not silently lost just because a covalent bond is absent. The atom-state array field names are recorded in `event.atom_state_fields`.

```
P = {g(Delta) : g is an automorphism of R}
N = |P| = |Aut(R)| / |Stab_R(Delta)|
D = product_i factorial(m_i)
n_e = N / D
```

`m_i` is the multiplicity of an identical reactant species occurrence, using its full graph and state. These are graph group orders, not rotational symmetry numbers. The implementation applies the correction exactly once and saves `correction_applications: 1` and `reactant_exchange_quotiented: false`.

The engine constructs generators using a point-stabilizer chain: fix preceding atom positions, individually mark a source/target atom, and compare their canonical colored graphs. For matching canonical SMILES, RDKit's `_smilesAtomOutputOrder` gives a concrete node correspondence. The transversals at each stabilizer level generate the full action on chemical atoms. Only the orbit of the **event** is then traversed. Interchangeable spectators do not require enumerating their factorially many full atom permutations. Molecular/bond auxiliary nodes encode structure and grouping; counts use the action on chemical atoms.

Every orbit step canonicalizes the exact labelled event, so duplicate actions add no path. The visited set is complete before the final factor is assigned. Up to 64 representative mapped paths are included in the evidence; `paths_truncated` reports display-sample truncation, which does not truncate the count. Group orders exceeding JavaScript's safe-integer range are decimal strings to preserve audit precision in browser JSON.

The statistical correction follows the previously selected RMG-style identical-reactant convention. It is separate from deleting duplicate representations. For example:

- H2 + X + X: two H-to-site assignments remain; `N=2`, `D=2!`, `n_e=1`.
- 2 ethyl radicals disproportionating: either labelled radical donates one of three beta-H atoms; `N=6`, `D=2!`, `n_e=3`.
- 2 methyl radicals recombining: one labelled new bond; `N=1`, `D=2!`, `n_e=0.5`.
- The derived three-H encounter forming H2 + H: choose one of three pairs; `N=3`, `D=3!`, `n_e=0.5`.

These fractional factors are intentional. The software does not also quotient the event orbit by identical-reactant exchanges and then apply another factorial. It does not incorporate `match_count`. Consumers are responsible for using these factors consistently in their chosen kinetic expressions; this module does not require or prescribe a rate equation, units or fitted parameters.

A full mapped event defines one symmetry orbit. Non-equivalent sites discovered by a general template remain separate channels. Different events that happen to have the same product species are not merged. Duplicate templates producing the same full mapped event share an event identifier irrespective of their family label.

## Abstract catalytic sites

In species SMILES, `*` is an abstract catalyst node; in SMARTS, `[#0]` selects it, while `*` is an any-atom query. `Species.from_smiles("[*:1]", phase="surface")` defaults to `{kind: "surface", type: "X", capacity: 1}`. Optional type/element/facet and other attributes distinguish nodes as graph colors; they do not trigger a geometry calculation. Constraints can be supplied in a template's `site_constraints` if needed for generation.

A node has at most one attached atom, with single through quadruple adsorption bond order. Multiple nodes describe multi-site species. A species can also contain a disconnected adsorbate/site association; the species partition preserves it. This is an abstract mean-field representation, not a real lattice. There are no required 3D coordinates, adjacency multiplicities or adsorption orientations, and none are added to `n_e`.

For the source-linked Brønsted examples, `[*:1]O` represents anchor–O–H and `[*:1][O-].C[CH+]C` is one associated framework-anion/carbenium species. Actual O/H and charge changes participate in mapping. This avoids treating acid proton transfer as metal adsorption. Multiple anchors in a declared associated species are handled as that supplied graph; no additional geometric association is guessed.

Site metadata and node inventory are conserved during template execution, and participate in species/channel identity and persistence. A surface species with only a legacy `site` label and no explicit graph node lacks the site representation needed for this catalytic interface.

## Statuses and limits

| Status | Meaning |
|---|---|
| `calculated` | Complete event orbit and correction computed; `factor` and `n_e` are present. |
| `unsupported_algorithm` | Input contains assigned tetrahedral or bond stereochemistry, but parity-preserving event actions are not implemented. This is an algorithm limit, not missing coordinates. |
| `missing_information` | Required atom correspondence or explicit catalytic nodes are absent. Generated valid channels normally already supply these. |
| `invalid_input` | Invalid structure, repeated maps on one side, nonconserved atoms/isotopes/charge/site attributes, invalid partition/site state, or no event. |
| `incomplete` | Search limit reached; no factor assigned. |
| `imported` / legacy `unsupported` | Annotation-layer imported factor / calculator-free or calculator-specific legacy annotation. |

Failures retain a machine-readable `reason_code` and readable `reason`. There is no default `n_e=1` on failure.

Elements, rings, aromatic molecules, localized radicals/ions, unknown families and extra site labels are not rejected by a family whitelist. Counting uses the supplied RDKit-normalized electronic graph, including aromatic bond representation. It does not enumerate a resonance ensemble or infer stereoisomer/facial multiplicities for unspecified stereo. This is the definition of the present graph calculation, not a request for additional TS data. Assigned stereo is explicitly refused until its permutation action is implemented correctly. The algorithm has finite search limits and does not claim polynomial performance for arbitrary graphs.

## Document-linked cases and results

The original documents were located and inspected; printed page numbers are used below. Input templates and independent hand expectations are in `tests/data/path_cases.json`; new templates are embedded in those test cases rather than added to the bundled family library.

| Reaction / source | Event and manual enumeration | N | D | Expected and actual n_e |
|---|---|---:|---:|---:|
| H + CH4 → CH3 + H2; Vandewiele 2015 p.187 Scheme 8 | Break one of four equivalent C–H, form H–H; move radical | 4 | 1 | 4 |
| CH3 + C2H5 disproportionation; Vandewiele p.189 Table 4 | Transfer one of three beta-H; C–C becomes C=C, both radicals removed | 3 | 1 | 3 |
| 2 C2H5 disproportionation; **derived** identical-reactant check | Two donor occurrences × three H | 6 | 2 | 3 |
| H2 + 2 M → 2 MH; Lozano-Blanco 2011 p.427 Table 2 (1) | Two distinct H-to-site assignments | 2 | 2 | 1 |
| 2 MH → H2 + 2 M; same equation, reverse | One labelled bond-forming event | 1 | 2 | 0.5 |
| M–C2H5 + M–H → C2H6 + 2 M; Table 2 (8), n=2 | Break C–M/H–M, form C–H; reactants differ | 1 | 1 | 1 |
| Acid + propene protonation; aromatization summary MTO Fig.1, Wang 2023 p.125 Table 2 | Two inequivalent alkene ends; one path in each channel; explicit O−/C+ | 1 each | 1 | 1 each |
| Associated isopropyl cation deprotonation; summary MTO Fig.1 reverse step | Two methyl groups × three H; framework O− accepts H | 6 | 1 | 6 |
| 3,4-dimethyl-2-pentyl carbenium beta-scission; Wang pp.128–129 Figs.7–8 | C2=C3 formation/C3–C4 cleavage in the selected source example; association retained | 1 | 1 | 1 |
| **New user template:** CH4 → CH3 + H; Vandewiele p.187 Scheme 7 | Four equivalent C–H homolysis positions; absent from bundled recipes | 4 | 1 | 4 |
| **New user template:** 2-methylheptan-3-yl to 2-methylheptan-2-yl carbenium hydride shift; Wang p.125 Table 2 | One H on adjacent tertiary donor; transfer H and positive charge in the explicit framework model | 1 | 1 | 1 |

The hydride-shift template restricts the donor to a tertiary carbon, corresponding to the illustrated branch. Other donor positions would be separate generated events; the calculator does not impose this restriction. Similarly, beta-scission's source template selects its secondary/tertiary leaving cation; the calculator has no such built-in chemistry filter. Unspecified stereo in these case inputs is counted constitutionally.

Vandewiele's Table 4 gives 6 for methyl/ethyl disproportionation under its molecular/TS symmetry convention. The present supplied-graph path count is 3. This is a documented convention difference, not a claim to reproduce that 6. The catalytic results do not use Wang's free-state symmetry-cancellation assumption. No TS or restricted-rotation factor is inferred from any of these graph counts.

Additional derived tests cover propane's equivalent/inequivalent positions (6 and 2), H/D distinctions (3 and 1), reordered and renumbered reactants/products, duplicate and early-stopped template matches, different mapped events with the same species products, generic oxygen chemistry, ring chemistry, explicit multi-bond surface attachments, same versus different site types, product association, missing/invalid inputs and search limits. Six small source-linked cases are cross-checked against an **independent exhaustive RDKit atom-self-match enumeration using complete labelled product graphs**, without the engine's event signatures or stabilizer-chain algorithm. A long-chain test checks that millions of spectator automorphisms do not need full enumeration.

## RMG comparison and verification

The counting convention was checked against fixed RMG-Py commit `032e9d20e953158e801f2de897dc5677d65bbccd`, particularly [find_degenerate_reactions / reduce_same_reactant_degeneracy](https://github.com/ReactionMechanismGenerator/RMG-Py/blob/032e9d20e953158e801f2de897dc5677d65bbccd/rmgpy/data/kinetics/common.py) and [calculate_degeneracy / _generate_reactions](https://github.com/ReactionMechanismGenerator/RMG-Py/blob/032e9d20e953158e801f2de897dc5677d65bbccd/rmgpy/data/kinetics/family.py). The earlier recipe implementation also inspected database commit `6391cc4f75341803605d9ad17417e14cf2db1f59`. The current engine does not load either database.

RMG is not installed in the local base/ml environments. **RMG runtime cross-validation has not been performed**, and is saved as false. Independent exhaustive graph enumeration and manual cases validate the implemented convention; neither implies agreement with every RMG family/resonance convention.

Verification: 125 tests passed in the web-enabled environment; 34 path tests passed in ml (Python 3.13.1 / RDKit 2025.03.5), as well as base RDKit 2026.03.6. The custom-template CLI example, Node syntax check, whitespace check and wheel build passed. A browser UI run entered the new CH homolysis template and methane seed manually, generated the network without family/recipe/parameter inputs, and displayed n_e=4 with 4/1 counting evidence; screenshot: `results/generic-ne-ui.png`. Web API tests submit new user gas/acid templates without family/kind/parameter/unit inputs, calculate factors, export CSV, save/open projects and regenerate. The only full-suite warning is an existing Starlette/httpx deprecation.

## Integration and changed files

- `path_degeneracy.py`: pure `calculate_path_factor`, generic event extraction, automorphism generators/orbits, correction, audit data and network adapter. Recipe lookup is gone from calculation.
- `single_event.py`, `__init__.py`: expose/register the new interface, preserve structured statuses, and make kinetic classes/units optional. No rate-law record is created for pure factor annotations.
- `sites.py`, `rules.py`: generic default X, multiple adsorption bond orders, and elementary-by-default user templates. The generator itself remains unchanged in this revision.
- `models.py`: retain new result fields in schema 6 JSON and append `n_e`, algorithm, reason code, event, raw path count and divisor columns to CSV. Existing schema 1–5 data remains readable; old annotation fields remain exportable.
- `web.py`, `static/index.html`, `static/app.js`: template validation and an editable custom-template JSON panel, generic default calculator when enabling annotation, factor-first inspection, and status-specific explanations. Imported JSON templates use the same validation and generation path.
- Tests, `examples/custom_ne_rules.json`, `examples/custom_ne.toml`, README translations and packaging manifest: reproducible custom-template use and regression coverage. Existing methane symmetry tests and optional parameter-import functionality are retained.

The browser workflow is: edit seed species → **Custom templates (JSON)** / Import → **Use these templates** → enable **Annotate after generation** → **Generate** → click a reaction diamond. **Apply to current network** annotates existing mapped channels. Save/Open and JSON/CSV exports retain the same factors and evidence. No fitting module was added.

## Literature

- Vandewiele et al., *Symmetry calculation for molecules and transition states*, J. Comput. Chem. 2015, 36, 181–192, [10.1002/jcc.23788](https://doi.org/10.1002/jcc.23788). Supplied local PDF; the previous SI investigation is not a runtime dependency.
- Lozano-Blanco et al., *Extension of the Single-Event Methodology to Metal Catalysis: Application to Fischer–Tropsch Synthesis*, OGST 2011, 66, 423–435, [10.2516/ogst/2009075](https://doi.org/10.2516/ogst/2009075). Supplied local filename begins `xtension…`.
- Wang et al., *Single event kinetic modeling for paraffin hydrocracking over an industrial Ni–W silica–aluminum catalyst*, RCE 2023, 8, 123–136, [10.1039/D2RE00286H](https://doi.org/10.1039/D2RE00286H).
- `烷烃芳构化总结.docx`, MTO Figure 1. Figure identifier is used because Word layout page numbering is not stable.

### Legacy classification regression (2026-11-09)

Removed the historical `overall` eligibility gate from both automatic and imported factor annotation. The browser displays “Elementary input” as the input assumption. A regression generates all four default templates from propane, checks calculated factors for every channel, opens a simulated old project containing `overall` labels and skipped results, reapplies annotation without regeneration, verifies Save/Open, and then regenerates the historical templates with identical factors.

Validation for this correction: 28 annotation/web tests passed, including the new default-template and legacy-project regression; JavaScript syntax and whitespace checks passed.

### Combined default network (2026-11-09)

The default rule library now includes the four gas templates plus `user_metal_H2_adsorption`. A new browser workspace seeds propane and one vacant metal-site species and enables optional path annotation for the demo. H2 is generated by dehydrogenation at depth 1 and feeds the catalytic adsorption at depth 2. Existing projects retain their own inputs. `examples/combined_gas_catalytic.toml` reproduces the combined network: 11 species, 11 reactions, propane dehydrogenation n_e=12, H2 adsorption n_e=1. This couples the templates through H2; it does not add a metal-mediated propane dehydrogenation mechanism. Full suite: 127 passed (one existing Starlette/httpx deprecation). Browser verification confirmed all five templates selected and the catalytic count 2/2=1 in the same network; screenshot `results/combined-gas-catalytic-ui.png`.

### Blank projects and explicit example loading (2026-11-09)

This supersedes the automatic browser default loading described above: a fresh workspace and **New** now contain no seeds or templates. **Load example** explicitly loads the five-template gas/catalytic setup, clears prior results, and does not generate. **Generate** uses only current inputs. The frontend no longer falls back to bundled rules when a request has no rules. README translations describe the new workflow.

Verification: 15 web API tests passed; JavaScript syntax and whitespace checks passed. Browser checks confirmed blank Generate requests seed input, Load example remains “Not generated”, changing example depth from 4 to 1 produces 4 species and 1 reaction, and New after that run clears inputs, graph, annotations and export availability. Screenshot: `results/new-blank-workspace.png`.

### Frontend cache consistency (2026-11-09)

The reported screenshot combined the new HTML controls with old JavaScript behavior (propane-only initialization and the raw catalyst template ID). The active server on port 8765 served the current on-disk JavaScript and had an idle empty workspace. Added versioned application-script/style URLs; the entry route derives versions from content hashes and returns `Cache-Control: no-store`. The static HTML also has a changed fallback version for already-running processes that still serve it directly. Tests check content-derived asset URLs, served asset bytes and the entry cache header. Actual port-8765 browser verification with a fresh URL confirmed New leaves zero seeds/templates and no generated network. Existing server data was already empty; the user's service was not restarted.
