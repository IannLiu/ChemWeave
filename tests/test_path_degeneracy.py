"""Source-linked elementary mechanisms plus explicitly labelled derived checks."""

from collections import Counter
from copy import deepcopy
import csv
import json
from pathlib import Path

import pytest
from rdkit import Chem

from chemweave import GenerationConfig, Generator, Network, ReactionRule, Species, annotate_single_events, path_family_rules, calculate_path_factor
from chemweave.mapping import reaction_symmetry_key
from chemweave.sites import site_label


CASES = json.loads((Path(__file__).parent / "data/path_cases.json").read_text())
METAL = {"kind": "metal", "type": "Fe_model", "capacity": 1}
ACID = {"kind": "bronsted", "type": "ZSM5_model", "capacity": 1}


def case_network(case):
    rule = ReactionRule(**case["template"]) if "template" in case else next(r for r in path_family_rules() if r.id == case["rule"])
    reactants = [Species.from_smiles(**s) for s in case["reactants"]]
    network = Generator([rule], GenerationConfig(max_carbon=10, max_depth=1)).generate(reactants)
    annotate_single_events(network, calculator="path_degeneracy_v1")
    wanted = Counter(s.key for s in reactants)
    annotations = [c["single_event"] for r in network.reactions for c in r.channels
                   if Counter(network.species[sid].key for sid in c["reactants"]) == wanted]
    return network, annotations


@pytest.mark.parametrize("case", CASES, ids=lambda x: x["id"])
def test_document_cases(case):
    network, annotations = case_network(case)
    observed = []
    for se in annotations:
        assert se["status"] == "calculated", se
        assert se["method"] == "reaction_path_degeneracy"
        counts = se["counting"]
        observed.append({"paths": counts["unique_labeled_paths"],
                         "divisor": counts["identical_reactant_divisor"], "factor": se["factor"]})
        assert not counts["reactant_exchange_quotiented"]
        assert not counts["rmg_runtime_cross_validated"]
        assert len(counts["paths"]) == counts["unique_labeled_paths"]
        assert "sigma_ts" not in se
        assert network.single_event["parameters"] == {}
        assert not {"parameter_class", "rate_basis", "A_units"} & se.keys()
        assert counts["automorphism_group_order"] == counts["event_stabilizer_order"] * counts["unique_labeled_paths"]
    assert sorted(observed, key=str) == sorted(case["expected"], key=str)


def test_equivalent_and_inequivalent_sites_without_molecule_lookup():
    rule = path_family_rules()[0]
    n = Generator([rule], GenerationConfig(max_depth=1)).generate([Species.from_smiles("CCC"), Species.from_smiles("[H]")])
    annotate_single_events(n, calculator="path_degeneracy_v1")
    annotations = [c["single_event"] for r in n.reactions for c in r.channels]
    assert sorted(a["factor"] for a in annotations) == [2, 6]
    assert len({a["kinetic_channel"] for a in annotations}) == 2
    acid_case = next(c for c in CASES if c["id"] == "propene_acid_protonation")
    _, annotations = case_network(acid_case)
    assert len({a["kinetic_channel"] for a in annotations}) == 2


def test_independent_enumeration_ignores_template_and_early_stop():
    redundant = ReactionRule(id="redundant", kind="elementary", smarts=
        "[C;X4:1](-[H:2])(-[H:4]).[#1;D0:3]>>[C;H0:1]-[H:4].[H:2]-[H:3]")
    seeds = [Species.from_smiles("C"), Species.from_smiles("[H]")]
    full = Generator([redundant], GenerationConfig(max_depth=1)).generate(seeds)
    stopped = Generator([redundant]).generate(seeds, targets=[Species.from_smiles("[CH3]")])
    assert full.reactions[0].channels[0]["match_count"] == 12
    assert stopped.reactions[0].channels[0]["match_count"] == 1
    for n in (full, stopped):
        annotate_single_events(n, calculator="path_degeneracy_v1")
        assert n.reactions[0].channels[0]["single_event"]["factor"] == 4


def test_occurrence_swaps_and_atom_relabeling_do_not_change_factor_or_class():
    case = next(c for c in CASES if c["id"] == "H2_adsorption_two_sites")
    first, a = case_network(case)
    changed = deepcopy(case)
    changed["reactants"] = [{"smiles": "[*:90]", "phase": "surface", "site_atoms": {90: METAL}},
                            {"smiles": "[H:90][H:1]"},
                            {"smiles": "[*:4]", "phase": "surface", "site_atoms": {4: METAL}}]
    second, b = case_network(changed)
    for field in ("factor", "kinetic_channel"):
        assert a[0][field] == b[0][field]
    assert a[0]["counting"]["unique_labeled_paths"] == b[0]["counting"]["unique_labeled_paths"] == 2
    # Deduplication leaves both fixed-ID H-to-site assignments. No second /2.
    assert len(first.reactions[0].channels) == len(second.reactions[0].channels) == 1


def test_identical_radical_recombination_event_rate_convention():
    rule = next(r for r in path_family_rules() if r.family == "R_Recombination")
    n = Generator([rule], GenerationConfig(max_depth=1)).generate([Species.from_smiles("[CH3]")])
    annotate_single_events(n, calculator="path_degeneracy_v1")
    se = n.reactions[0].channels[0]["single_event"]
    assert se["counting"]["unique_labeled_paths"] == 1
    assert se["factor"] == 0.5
    # The loss of two molecules per event is separate from the event coefficient.
    assert 2 * se["factor"] == 1


def test_isotope_paths_and_site_types_are_not_collapsed():
    n = Generator([path_family_rules()[0]], GenerationConfig(max_depth=1)).generate(
        [Species.from_smiles("[2H]C"), Species.from_smiles("[H]")])
    annotate_single_events(n, calculator="path_degeneracy_v1")
    assert sorted(c["single_event"]["factor"] for r in n.reactions for c in r.channels) == [1, 3]
    case = deepcopy(next(c for c in CASES if c["id"] == "H2_adsorption_two_sites"))
    case["reactants"][2]["site_atoms"]["1"]["type"] = "Fe_step"
    _, annotations = case_network(case)
    assert annotations[0]["counting"]["identical_reactant_divisor"] == 1
    assert annotations[0]["factor"] == 2


def test_site_identity_mapping_conservation_and_roundtrip(tmp_path):
    case = next(c for c in CASES if c["id"] == "propene_acid_protonation")
    n, _ = case_network(case)
    for r in n.reactions:
        for c in r.channels:
            assert len(c["site_atoms"]) == 1
            assert sum(n.species[sid].charge for sid in c["reactants"]) == sum(n.species[sid].charge for sid in c["products"]) == 0
            product = n.species[c["products"][0]]
            assert sorted(a.GetFormalCharge() for a in product._mapped_mol.GetAtoms() if a.GetFormalCharge()) == [-1, 1]
    data = json.loads(json.dumps(n.to_dict()))
    restored = Network.from_dict(data)
    assert [s.key for s in n.species.values()] == [s.key for s in restored.species.values()]
    annotate_single_events(restored, **restored.single_event["options"])
    assert json.loads(json.dumps(restored.to_dict())) == data
    restored.export_csv(tmp_path)
    with (tmp_path / "single_event.csv").open(encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    assert rows[0]["method"] == "reaction_path_degeneracy" and json.loads(rows[0]["counting"])["paths"]
    a = Species.from_smiles("[*:1]CC[*:2]", phase="surface", site_atoms={1: METAL, 2: METAL})
    b = Species.from_smiles("[*:99]CC[*:88]", phase="surface", site_atoms={99: METAL, 88: METAL})
    assert len(a.site_atoms) == 2 and a.key == b.key
    changed = {**METAL, "type": "other"}
    c = Species.from_smiles("[*:1]CC[*:2]", phase="surface", site_atoms={1: METAL, 2: changed})
    assert a.key != c.key
    assert reaction_symmetry_key([a._mapped_mol], [a._mapped_mol], [("surface", None)], [("surface", None)]) != reaction_symmetry_key(
        [c._mapped_mol], [c._mapped_mol], [("surface", None)], [("surface", None)])


def test_site_constraints_capacity_and_conservation():
    assert Species.from_smiles("*", phase="surface").site_atoms
    with pytest.raises(ValueError, match="capacity"):
        Species.from_smiles("C[*:1]C", phase="surface", site_atoms={1: METAL})
    with pytest.raises(ValueError, match="wildcard"):
        ReactionRule(id="wildcard", smarts="[*:1]>>[*:1]", site_constraints={1: {"kind": "metal"}})
    rule = next(r for r in path_family_rules() if r.id == "path_h2_adsorption")
    acid = Species.from_smiles("[*:1]O", phase="surface", site_atoms={1: ACID})
    assert not rule.matches(acid, 1)
    wrong = rule.to_dict()
    wrong["site_constraints"] = {3: {"kind": "metal", "type": "other"}, 4: {"kind": "metal"}}
    m = Species.from_smiles("[*:1]", phase="surface", site_atoms={1: METAL})
    n = Generator([ReactionRule(**wrong)]).generate([Species.from_smiles("[H][H]"), m])
    assert not n.reactions
    # A source node cannot disappear even if a template forgets to retain it.
    drop = ReactionRule(id="drop", smarts="[#0:1][H:2]>>[H:2]", reactant_phases=["surface"])
    n = Generator([drop]).generate([Species.from_smiles("[*:1][H]", phase="surface", site_atoms={1: METAL})])
    assert not n.reactions and n.stats["rejected"]["atom_mapping"]


def test_duplicate_rules_share_kinetic_channel_not_double_rate():
    rule = path_family_rules()[0]
    twin = ReactionRule(**{**rule.to_dict(), "id": "another_template"})
    n = Generator([rule, twin], GenerationConfig(max_depth=1)).generate([Species.from_smiles("C"), Species.from_smiles("[H]")])
    annotate_single_events(n, calculator="path_degeneracy_v1")
    assert len(n.reactions) == 2
    assert len(n.single_event["kinetic_channels"]) == 1


def test_same_products_do_not_merge_different_mapped_mechanisms():
    # Derived counterexample: the incoming H attaches to C while two original
    # methane H atoms form H2. Same species products, different bond-change path.
    exchange = ReactionRule(id="exchange", kind="elementary", family="H_Abstraction", smarts=
        "[C;X4:1](-[#1:2])(-[#1:4]).[#1;D0:3]>>[C;H0:1]-[#1:3].[#1:2]-[#1:4]")
    n = Generator([path_family_rules()[0], exchange], GenerationConfig(max_depth=1)).generate(
        [Species.from_smiles("C"), Species.from_smiles("[H]")])
    annotate_single_events(n, calculator="path_degeneracy_v1")
    assert len(n.reactions) == 2 and n.reactions[0].products == n.reactions[1].products
    first, other = [r.channels[0]["single_event"] for r in n.reactions]
    assert first["factor"] == 4
    assert other["status"] == "calculated" and other["factor"] == 6
    assert first["kinetic_channel"] != other["kinetic_channel"]


def test_assigned_stereo_is_an_algorithm_limit_not_missing_input():
    n = Generator([path_family_rules()[0]], GenerationConfig(max_depth=1)).generate(
        [Species.from_smiles("C[C@H](CC)CCC"), Species.from_smiles("[H]")])
    annotate_single_events(n, calculator="mapped_event_v1")
    assert n.reactions
    assert all(c["single_event"]["status"] == "unsupported_algorithm" and
               c["single_event"]["reason_code"] == "stereochemical_action"
               for r in n.reactions for c in r.channels)


def test_no_calculator_is_not_one_and_import_still_works():
    from chemweave import load_rules
    n = Generator(load_rules(), GenerationConfig(max_depth=1)).generate([Species.from_smiles("CCC")])
    annotate_single_events(n)
    assert all(c["single_event"]["status"] == "unsupported" and "factor" not in c["single_event"] for r in n.reactions for c in r.channels)
    case = deepcopy(next(c for c in CASES if c["id"] == "H2_adsorption_two_sites"))
    for s in case["reactants"][1:]:
        s["site_atoms"]["1"]["orientation"] = "unresolved"
    n, annotations = case_network(case)
    # Additional site attributes are symmetry colors, not a new geometry model.
    assert annotations[0]["status"] == "calculated" and annotations[0]["factor"] == 1
    annotate_single_events(n, calculator="path_degeneracy_v1", factors=[{
        "reaction_id": n.reactions[0].id, "channel": 1, "factor": 2.5,
        "convention": "Independent surface calculation",
        "source": "Test import"}])
    assert n.reactions[0].channels[0]["single_event"]["status"] == "imported"


def test_new_templates_never_consult_recipe_library(monkeypatch):
    import chemweave.path_degeneracy as module
    rule = ReactionRule(id="user_oxygen_homolysis", smarts="[O:1]-[#1:2]>>[O;H0:1].[#1;H0:2]")
    assert rule.kind == "elementary"
    n = Generator([rule], GenerationConfig(max_depth=1)).generate([Species.from_smiles("O")])
    def no_recipes():
        raise AssertionError("The pure factor calculation must not load family recipes")
    monkeypatch.setattr(module, "path_family_rules", no_recipes)
    annotate_single_events(n, calculator="mapped_event_v1")
    se = n.reactions[0].channels[0]["single_event"]
    assert se["factor"] == 2 and se["family"] == "unspecified"
    assert n.single_event["parameters"] == {}
    assert not {"rate_law", "rate_definition"} & n.single_event.keys()


def test_three_identical_occurrences_and_incomplete_search():
    # Derived counting check: a three-H encounter, two form H2, one is unchanged.
    left = ["[H:1]", "[H:2]", "[H:3]"]
    right = ["[H:1][H:2]", "[H:3]"]
    result = calculate_path_factor(left, right)
    assert result["factor"] == 0.5
    assert result["counting"]["unique_labeled_paths"] == 3
    assert result["counting"]["identical_reactant_divisor"] == 6
    assert result["counting"]["correction_applications"] == 1
    for kwargs in ({"max_paths": 1}, {"max_checks": 1}):
        limited = calculate_path_factor(left, right, **kwargs)
        assert limited["status"] == "incomplete" and "factor" not in limited


def test_missing_information_and_invalid_input_are_not_algorithm_limits():
    incomplete = calculate_path_factor(["[CH4:1]"], ["[CH3:1]", "[H:2]"])
    assert incomplete["status"] == "missing_information" and incomplete["reason_code"] == "incomplete_atom_mapping"
    duplicate = calculate_path_factor(["[H:1]", "[H:1]"], ["[H:1][H:2]"])
    assert duplicate["status"] == "invalid_input" and duplicate["reason_code"] == "duplicate_atom_maps"
    loss = calculate_path_factor(["[H:1][H:2]"], ["[H:1]"])
    assert loss["status"] == "invalid_input" and "factor" not in loss
    legacy_site = calculate_path_factor(["[H:1]", "[H:2]"], ["[H:1][H:2]"],
                                        reactant_states=[("surface", "M"), ("gas", None)])
    assert legacy_site["status"] == "missing_information" and legacy_site["reason_code"] == "missing_site_nodes"
    empty = calculate_path_factor(["[H:1][H:2]"], ["[H:1][H:2]"])
    assert empty["status"] == "invalid_input" and empty["reason_code"] == "empty_event"


def test_rings_and_localized_radical_graphs_have_no_family_scope_gate():
    rule = ReactionRule(id="user_homolysis", smarts="[C;X4:1]-[#1:2]>>[C;H0:1].[#1;H0:2]")
    for smiles, expected in [("C1CCCCC1", 12), ("[CH2]C=CC", 3)]:
        n = Generator([rule], GenerationConfig(max_depth=1)).generate([Species.from_smiles(smiles)])
        annotate_single_events(n, calculator="mapped_event_v1")
        assert n.reactions and n.reactions[0].channels[0]["single_event"]["factor"] == expected


def test_terminal_hydrogen_spectators_do_not_require_full_automorphism_enumeration():
    rule = ReactionRule(id="long_chain", smarts="[C;X4:1]-[#1:2]>>[C;H0:1].[#1;H0:2]")
    n = Generator([rule], GenerationConfig(max_carbon=20)).generate([Species.from_smiles("C" * 20)], targets=[Species.from_smiles("[H]")])
    annotate_single_events(n, calculator="mapped_event_v1")
    se = n.reactions[0].channels[0]["single_event"]
    assert se["factor"] == 6
    assert se["counting"]["automorphism_group_order"] > 1_000_000
    assert se["counting"]["isomorphism_checks"] < 200


def test_generic_surface_nodes_multibonds_and_associated_product_partition():
    # One actual C atom binds X with double/triple/quadruple bond. Occupancy is
    # one site per node, independent of adsorption bond order.
    for bond in ("=", "#", "$"):
        s = Species.from_smiles("[*:1]" + bond + "C", phase="surface")
        assert len(s.site_atoms) == 1 and next(iter(s.site_atoms.values()))["kind"] == "surface"
    first = calculate_path_factor(["[H:1][H:2]", "[*:3]", "[*:4]"], ["[H:1][*:3]", "[H:2][*:4]"])
    assert first["factor"] == 1
    # Grouping both free sites in ONE species changes the statistical divisor.
    paired = calculate_path_factor(["[H:1][H:2]", "[*:3].[*:4]"], ["[H:1][*:3]", "[H:2][*:4]"])
    assert paired["factor"] == 2 and paired["counting"]["identical_reactant_divisor"] == 1
    assert first["kinetic_channel"] != paired["kinetic_channel"]
    # Association alone is an event; no fabricated forming bond is required.
    association = calculate_path_factor(["[Na+:1]", "[Cl-:2]"], ["[Na+:1].[Cl-:2]"])
    assert association["factor"] == 1 and not association["event"]["bond_changes"]


@pytest.mark.parametrize("case_id", ["methane_H", "methyl_ethyl_disproportionation", "ethyl_self_disproportionation",
                                    "H2_adsorption_two_sites", "H2_desorption_two_adsorbates", "isopropyl_acid_deprotonation"])
def test_against_independent_full_atom_permutation_enumeration(case_id):
    """Small graphs only: enumerate ALL chemical graph self-matches directly.

    This oracle uses labelled full products, not event signatures or the engine's
    colored auxiliary graph / stabilizer-chain generators.
    """
    case = next(c for c in CASES if c["id"] == case_id)
    n, annotations = case_network(case)
    wanted = Counter(Species.from_smiles(**s).key for s in case["reactants"])
    ch = next(c for r in n.reactions for c in r.channels if Counter(n.species[sid].key for sid in c["reactants"]) == wanted)
    params = Chem.SmilesParserParams(); params.removeHs = False
    left = [Chem.MolFromSmiles(s, params) for s in ch["mapped_reactant_smiles"]]
    combined = left[0]
    for mol in left[1:]:
        combined = Chem.CombineMols(combined, mol)
    numbers = [a.GetAtomMapNum() for a in combined.GetAtoms()]
    labels = [(a.GetAtomicNum(), a.GetIsotope(), a.GetFormalCharge(), a.GetNumRadicalElectrons(),
               json.dumps(ch.get("site_atoms", {}).get(a.GetAtomMapNum(), {}), sort_keys=True)) for a in combined.GetAtoms()]
    parts = sorted((n.species[sid].phase, tuple(sorted(a.GetAtomMapNum() for a in mol.GetAtoms()))) for sid, mol in zip(ch["reactants"], left))
    matches = combined.GetSubstructMatches(combined, uniquify=False, maxMatches=100000)
    assert len(matches) < 100000
    outcomes = set()
    for match in matches:
        if any(labels[i] != labels[j] for i, j in enumerate(match)):
            continue
        move = dict(zip(numbers, (numbers[j] for j in match)))
        if sorted((phase, tuple(sorted(move[i] for i in members))) for phase, members in parts) != parts:
            continue
        output = []
        for sid, smi in zip(ch["products"], ch["mapped_product_smiles"]):
            mol = Chem.MolFromSmiles(smi, params)
            for atom in mol.GetAtoms():
                atom.SetAtomMapNum(move[atom.GetAtomMapNum()])
            output.append((n.species[sid].phase, Chem.MolToSmiles(mol)))
        outcomes.add(tuple(sorted(output)))
    assert len(outcomes) == annotations[0]["counting"]["unique_labeled_paths"]
