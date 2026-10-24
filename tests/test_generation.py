"""Chemical and enumeration regressions, independent of the demo graph size."""

from collections import Counter
import csv
import json

import pytest

from chemweave import GenerationConfig, Generator, Network, ReactionRule, Species, load_rules


@pytest.fixture
def demo_rules():
    return {r.id: r for r in load_rules()}


def test_species_identity_and_explicit_hydrogen():
    propane = Species.from_smiles("CCC")
    assert propane.key == Species.from_smiles("[CH3:5][CH2:2][CH3:8]").key
    assert Species.from_smiles("[H]C([H])([H])[H]").key == Species.from_smiles("C").key
    assert propane.formula == "C3H8"
    assert Species.from_smiles("[H][H]").elements == (("H", 2),)
    assert propane.key != Species.from_smiles("CCC", phase="surface", site="BAS").key
    assert Species.from_smiles("CC=CC").key != Species.from_smiles("C=CCC").key


def test_dehydrogenation_keeps_hydrogen_and_match_provenance(demo_rules):
    net = Generator([demo_rules["alkane_dehydrogenation"]]).generate([Species.from_smiles("CCC")])
    assert {s.formula for s in net.species.values()} == {"C3H8", "C3H6", "H2"}
    assert len(net.reactions) == 1
    reaction = net.reactions[0]
    assert len(reaction.matches) > 1  # equivalent template mappings, not SEMK degeneracy
    assert len(reaction.products) == 2
    match = reaction.matches[0]
    indexed = Species.from_smiles(match["indexed_reactant_smiles"][0])
    assert indexed.smiles == "CCC"
    assert any(atom.get("old_mapno") == 1
               for product in match["product_atom_properties"] for atom in product)


def test_new_olefin_can_react_with_itself(demo_rules):
    rules = [demo_rules["alkane_dehydrogenation"], demo_rules["terminal_olefin_coupling"]]
    net = Generator(rules, GenerationConfig(max_carbon=6)).generate([Species.from_smiles("CCC")])
    propene = next(sid for sid, s in net.species.items() if s.smiles == "C=CC")
    coupling = next(r for r in net.reactions if r.family == "oligomerization")
    assert coupling.reactants == ((propene, 2),)
    assert net.species[coupling.products[0][0]].formula == "C6H12"
    assert net.termination == "exhausted"


def test_binary_slots_and_seed_order(demo_rules):
    rule = demo_rules["terminal_olefin_coupling"]
    # Asymmetric slots require ethene first and propene second.
    rule = ReactionRule(
        id="asymmetric", family=rule.family,
        smarts="[CH2:1]=[CH2:2].[CH2:3]=[CH:4][CH3:5]>>[C:1]-[C:2]-[C:3]=[C:4]-[C:5]",
    )
    for inputs in (("C=C", "C=CC"), ("C=CC", "C=C")):
        net = Generator([rule], GenerationConfig(max_carbon=5)).generate(
            [Species.from_smiles(s) for s in inputs]
        )
        assert len(net.reactions) == 1
        assert any(s.formula == "C5H10" for s in net.species.values())


def test_isomerization_cycles_terminate_and_keep_directions(demo_rules):
    net = Generator([demo_rules["olefin_double_bond_shift"]]).generate([Species.from_smiles("C=CCC")])
    assert {s.smiles for s in net.species.values()} == {"C=CCC", "CC=CC"}
    assert len(net.reactions) == 2
    assert net.reactions[0].reactants == net.reactions[1].products
    assert net.reactions[0].products == net.reactions[1].reactants
    assert net.termination == "exhausted"


def test_distinct_rules_with_same_stoichiometry_are_retained(demo_rules):
    first = demo_rules["alkane_dehydrogenation"]
    second = ReactionRule(**{**first.to_dict(), "id": "alternative_mechanism"})
    net = Generator([first, second]).generate([Species.from_smiles("CCC")])
    assert len(net.reactions) == 2
    assert len(net.species) == 3


@pytest.mark.parametrize("smarts,reason", [
    ("[CH3:1][CH3:2]>>[CH2:1]=[CH2:2]", "unbalanced"),
    ("[CH3:1][CH3:2]>>[CH2+:1][CH3:2].[H+]", "unbalanced"),
    ("[C:1]>>[C:1](C)(C)(C)(C)C", "invalid_structure"),
])
def test_bad_candidate_products_are_discarded(smarts, reason):
    net = Generator([ReactionRule(id="bad", smarts=smarts)]).generate([Species.from_smiles("CC")])
    assert not net.reactions
    assert len(net.species) == 1
    assert net.stats["rejected"][reason] > 0


def test_surface_state_matching_and_site_conservation(demo_rules):
    spec = demo_rules["alkane_dehydrogenation"].to_dict()
    spec.update(reactant_phases=["surface"], reactant_sites=["BAS"],
                product_phases=["surface", "gas"], product_sites=["BAS", None])
    surface_rule = ReactionRule(**spec)
    gas = Species.from_smiles("CCC")
    surface = Species.from_smiles("CCC", phase="surface", site="BAS")
    net = Generator([surface_rule]).generate([gas, surface])
    assert len(net.reactions) == 1
    assert net.reactions[0].reactants == ((net.seeds[1], 1),)
    assert any(s.site == "BAS" and s.formula == "C3H6" for s in net.species.values())
    spec.update(product_phases=["gas", "gas"], product_sites=[None, None])
    net = Generator([ReactionRule(**spec)]).generate([surface])
    assert not net.reactions
    assert net.stats["rejected"]["unbalanced"] > 0


def test_demo_network_conserves_elements_charge_sites_and_references():
    net = Generator(load_rules(), GenerationConfig(max_carbon=6)).generate([Species.from_smiles("CCC")])
    assert net.termination == "exhausted"
    for reaction in net.reactions:
        totals = []
        for terms in (reaction.reactants, reaction.products):
            elements, sites, charge = Counter(), Counter(), 0
            for sid, coefficient in terms:
                species = net.species[sid]
                elements.update({e: n * coefficient for e, n in species.elements})
                charge += species.charge * coefficient
                if species.site is not None:
                    sites[species.site] += coefficient
                assert species.carbon_count <= 6
            totals.append((elements, sites, charge))
        assert totals[0] == totals[1]
    # A completed run should be deterministic on repeated use of the same rules.
    again = Generator(load_rules(), GenerationConfig(max_carbon=6)).generate([Species.from_smiles("CCC")])
    assert [(r.rule_id, r.reactants, r.products) for r in net.reactions] == [
        (r.rule_id, r.reactants, r.products) for r in again.reactions
    ]


@pytest.mark.parametrize("limits,termination", [
    ({"max_species": 2}, "max_species"),
    ({"max_reactions": 1}, "max_reactions"),
])
def test_budget_stop_preserves_consistent_partial_network(limits, termination):
    net = Generator(load_rules(), GenerationConfig(max_carbon=6, **limits)).generate([Species.from_smiles("CCC")])
    assert net.termination == termination
    referenced = set(net.seeds)
    for r in net.reactions:
        referenced.update(sid for sid, _ in r.reactants + r.products)
    assert referenced == set(net.species)


def test_roundtrip_and_exports_preserve_full_network(tmp_path):
    net = Generator(load_rules(), GenerationConfig(max_carbon=6)).generate([Species.from_smiles("CCC")])
    path = tmp_path / "network.json"
    net.export_json(path)
    restored = Network.from_json(path)
    second_path = tmp_path / "restored.json"
    restored.export_json(second_path)
    assert json.loads(path.read_text()) == json.loads(second_path.read_text())
    restored.export_csv(tmp_path)
    restored.export_dot(tmp_path / "network.dot")
    with (tmp_path / "species.csv").open(encoding="utf-8-sig") as f:
        assert len(list(csv.DictReader(f))) == len(net.species)
    with (tmp_path / "stoichiometry.csv").open(encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
        assert any(row["coefficient"] == "2" for row in rows)
    dot = (tmp_path / "network.dot").read_text()
    assert "shape=box" in dot
    assert len(net.rules) == 4

