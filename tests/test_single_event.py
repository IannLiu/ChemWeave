"""Physical symmetry examples and independence from template enumeration."""

from copy import deepcopy
import csv
import json
from pathlib import Path

import pytest

from chemweave import (GenerationConfig, Generator, Network, ReactionRule, Species,
                       annotate_single_events, load_rules, reaction_center, register_single_event_calculator)


CALCULATOR = "methane_h_abstraction_v1"
RULES = Path(__file__).resolve().parents[1] / "src/chemweave/data/single_event_rules.json"


def methane_network(smiles="C", rules=None, atom="[H]"):
    return Generator(rules or load_rules(RULES), GenerationConfig(max_depth=1)).generate(
        [Species.from_smiles(smiles), Species.from_smiles(atom)])


def test_default_annotation_does_not_choose_a_family():
    network = annotate_single_events(methane_network())
    se = network.reactions[0].channels[0]["single_event"]
    assert se["reaction_center"]["bond_changes"]
    assert se["status"] == "unsupported" and "factor" not in se
    assert network.single_event["options"]["calculator"] is None
    assert network.single_event["parameters"] == {}


def test_register_another_family_and_replay_saved_options():
    rule = ReactionRule(id="chlorine_abstraction", kind="elementary", family="halogen_abstraction",
                        smarts="[#1;D0:1].[Cl:2]-[Cl:3]>>[#1:1]-[Cl:2].[Cl;H0:3]")
    network = Generator([rule], GenerationConfig(max_depth=1)).generate(
        [Species.from_smiles("[H]"), Species.from_smiles("ClCl")])
    assert len(network.reactions) == 1

    def external_calculator(*, network, reaction, channel, center):
        if reaction.family != "halogen_abstraction":
            return None
        assert len(center["bond_changes"]) == 2
        assert all(sid in network.species for sid in channel["reactants"])
        return {"kinetic_subclass": "H+Cl2", "parameter_class": "H+Cl2-fit",
                "sigma_reactants": [1, 2], "sigma_ts": 1,
                "convention": "Gas-phase collinear H-Cl-Cl saddle; proper rotations",
                "source": "Test-only external family adapter"}

    register_single_event_calculator("test_halogen_abstraction", external_calculator)
    annotate_single_events(network, calculator="test_halogen_abstraction")
    annotation = network.reactions[0].channels[0]["single_event"]
    assert annotation["status"] == "calculated" and annotation["factor"] == 2
    assert annotation["parameter_class"] == "H+Cl2-fit"
    restored = Network.from_dict(json.loads(json.dumps(network.to_dict())))
    annotate_single_events(restored, **restored.single_event["options"])
    assert restored.reactions[0].channels[0]["single_event"] == annotation
    # Applicability is owned by the registered adapter, not the core framework.
    unrelated = annotate_single_events(methane_network(), calculator="test_halogen_abstraction")
    assert unrelated.reactions[0].channels[0]["single_event"]["status"] == "unsupported"


@pytest.mark.parametrize("smiles,expected", [
    ("C", [(12, 3, 4)]),
    ("[2H]C", [(3, 1, 3), (3, 3, 1)]),
    ("[2H]C[2H]", [(2, 1, 2), (2, 1, 2)]),
    ("[2H]C([2H])[2H]", [(3, 1, 3), (3, 3, 1)]),
    ("[2H]C([2H])([2H])[2H]", [(12, 3, 4)]),
])
def test_methane_hd_rotational_symmetry(smiles, expected):
    network = annotate_single_events(methane_network(smiles), calculator=CALCULATOR)
    observed = []
    for reaction in network.reactions:
        for channel in reaction.channels:
            se = channel["single_event"]
            assert se["status"] == "calculated"
            assert len(se["reaction_center"]["atom_maps"]) == 3
            assert len(se["reaction_center"]["bond_changes"]) == 2
            assert len(se["reaction_center"]["atom_changes"]) == 2
            observed.append((se["sigma_reactants"][0], se["sigma_ts"], se["factor"]))
            parameter = network.single_event["parameters"][se["parameter_class"]]
            assert parameter["A"] is parameter["n"] is parameter["Ea"] is None
            assert parameter["status"] == "unfitted"
    assert sorted(observed) == sorted(expected)
    assert len(network.single_event["parameters"]) == len(expected)


def test_factor_independent_of_template_multiplicity_and_partial_counters():
    original = load_rules(RULES)[0]
    redundant = ReactionRule(id="extra_spectator", kind="elementary", smarts=
        "[C;X4;H4:1](-[H:2])(-[H:4]).[#1;D0:3]>>[C;H0:1]-[H:4].[H:2]-[H:3]")
    counts = []
    subclasses = []
    for rule in (original, redundant):
        network = methane_network(rules=[rule])
        channel = network.reactions[0].channels[0]
        counts.append(channel["match_count"])
        annotate_single_events(network, calculator=CALCULATOR)
        se = channel["single_event"]
        assert se["factor"] == 4
        subclasses.append(se["kinetic_subclass"])
        channel["match_count"] = 1  # An early-stop counter has no physical role.
        annotate_single_events(network, calculator=CALCULATOR)
        assert channel["single_event"]["factor"] == 4
    assert counts == [4, 12]
    assert subclasses[0] == subclasses[1]
    stopped = Generator([original]).generate(
        [Species.from_smiles("C"), Species.from_smiles("[H]")], targets=[Species.from_smiles("[CH3]")])
    assert stopped.termination == "targets_reached"
    assert stopped.reactions[0].channels[0]["match_count"] == 1
    annotate_single_events(stopped, calculator=CALCULATOR)
    assert stopped.reactions[0].channels[0]["single_event"]["factor"] == 4


def test_map_numbers_do_not_change_class_or_factor():
    classes = []
    for smiles, atom in (("C", "[H]"), ("[CH4:98]", "[H:98]")):
        network = annotate_single_events(methane_network(smiles, atom=atom), calculator=CALCULATOR)
        se = network.reactions[0].channels[0]["single_event"]
        assert se["factor"] == 4
        classes.append(se["kinetic_subclass"])
    assert classes[0] == classes[1]


def test_imported_factors_share_parameters_and_roundtrip(tmp_path):
    network = methane_network("[2H]C")
    factors = [{"reaction_id": r.id, "channel": 1, "factor": 2.5,
                "kinetic_subclass": "externally_classified", "parameter_class": "shared_fit",
                "convention": "User-defined path convention", "source": "Independent external calculation"}
               for r in network.reactions]
    parameters = {"shared_fit": {"A": 1e6, "n": 0, "Ea": 30000,
                                 "status": "fitted", "source": "Synthetic test data; not an experimental result"}}
    annotate_single_events(network, calculator=None, factors=factors, parameters=parameters)
    assert len(network.single_event["parameters"]) == 1
    assert all(r.channels[0]["single_event"]["parameter_class"] == "shared_fit" for r in network.reactions)
    path = tmp_path / "network.json"
    network.export_json(path)
    restored = Network.from_json(path)
    assert json.loads(json.dumps(restored.to_dict())) == json.loads(json.dumps(network.to_dict()))
    assert restored.single_event["parameters"]["shared_fit"]["A"] == 1e6
    restored.export_csv(tmp_path)
    with (tmp_path / "single_event.csv").open(encoding="utf-8-sig") as f:
        records = list(csv.DictReader(f))
    assert len(records) == 2 and all(float(r["factor"]) == 2.5 for r in records)
    with (tmp_path / "kinetic_parameters.csv").open(encoding="utf-8-sig") as f:
        assert len(list(csv.DictReader(f))) == 1
    # Saved options can reproduce the exact annotation after opening in the UI.
    annotate_single_events(restored, **restored.single_event["options"])
    assert json.loads(json.dumps(restored.to_dict())) == json.loads(json.dumps(network.to_dict()))


def test_imported_symmetry_ratio_and_transactional_errors():
    network = annotate_single_events(methane_network(), calculator=CALCULATOR)
    record = {"reaction_id": "R00001", "channel": 1, "sigma_reactants": [12, 1], "sigma_ts": 3,
              "kinetic_subclass": "imported", "parameter_class": "fit", "source": "Reference",
              "convention": "Forward rotational symmetry ratio"}
    annotate_single_events(network, factors=[record])
    assert network.reactions[0].channels[0]["single_event"]["factor"] == 4
    before = deepcopy(network.to_dict())
    for change in ({"factor": 8}, {"channel": 0}, {"factor": float("nan")}, {"source": ""}, {"sigma_ts": 0}):
        with pytest.raises(ValueError):
            annotate_single_events(network, factors=[{**record, **change}])
        assert network.to_dict() == before
    with pytest.raises(ValueError, match="Duplicate"):
        annotate_single_events(network, factors=[record, record])
    with pytest.raises(ValueError, match="Unknown parameter"):
        annotate_single_events(network, parameters={"nonexistent": {"A": 3}})
    assert network.to_dict() == before


def test_opt_in_overall_unsupported_other_elementary_and_old_schema():
    network = Generator(load_rules(), GenerationConfig(max_depth=1)).generate([Species.from_smiles("CCC")])
    assert network.single_event is None
    assert all("single_event" not in c for r in network.reactions for c in r.channels)
    annotate_single_events(network, calculator=CALCULATOR)
    se = network.reactions[0].channels[0]["single_event"]
    assert se["status"] == "unsupported" and "factor" not in se
    assert se["reaction_center"]["bond_changes"]
    with pytest.raises(ValueError, match="elementary"):
        annotate_single_events(network, factors=[{"reaction_id": network.reactions[0].id, "channel": 1}])
    # Even a methane-like mapping must satisfy the calculator's phase convention.
    other = methane_network()
    other.reactions[0].kind = "overall"
    annotate_single_events(other, calculator=CALCULATOR)
    assert other.reactions[0].channels[0]["single_event"]["status"] == "unsupported"
    data = methane_network().to_dict()
    data["schema_version"] = 4
    data.pop("single_event")
    restored = Network.from_dict(data)
    assert restored.single_event is None
    annotate_single_events(restored, calculator=CALCULATOR)
    assert restored.reactions[0].channels[0]["single_event"]["factor"] == 4
    restored.reactions[0].channels.clear()
    annotate_single_events(restored, calculator=CALCULATOR)
    assert restored.single_event["unmapped_reactions"] == ["R00001"]


def test_calculator_scope_excludes_surface_reverse_and_mixed_protium_labels():
    surface_rule = load_rules(RULES)[0].to_dict()
    surface_rule.update(reactant_phases=["surface", "gas"], product_phases=["surface", "gas"])
    network = Generator([ReactionRule(**surface_rule)], GenerationConfig(max_depth=1)).generate(
        [Species.from_smiles("C", phase="surface"), Species.from_smiles("[H]")])
    annotate_single_events(network, calculator=CALCULATOR)
    assert network.reactions[0].channels[0]["single_event"]["status"] == "unsupported"
    network = methane_network()
    reaction = network.reactions[0]
    channel = reaction.channels[0]
    reaction.reactants, reaction.products = reaction.products, reaction.reactants
    channel["reactants"], channel["products"] = channel["products"], channel["reactants"]
    channel["mapped_reaction_smiles"] = ">>".join(reversed(channel["mapped_reaction_smiles"].split(">>")))
    annotate_single_events(network, calculator=CALCULATOR)
    assert channel["single_event"]["status"] == "unsupported"
    mixed = annotate_single_events(methane_network("[1H]C"), calculator=CALCULATOR)
    assert mixed.reactions
    assert all(c["single_event"]["status"] == "unsupported" for r in mixed.reactions for c in r.channels)


def test_reaction_center_bond_order_charge_and_invalid_mapping():
    center = reaction_center("[CH3:1][CH3:2]>>[CH2:1]=[CH2:2]")
    assert center["bond_changes"][0]["before"]["order"] == 1
    assert center["bond_changes"][0]["after"]["order"] == 2
    center = reaction_center("[NH4+:1].[OH-:2]>>[NH3:1].[OH2:2]")
    assert center["atom_maps"] == [1, 2]
    assert len(center["atom_changes"]) == 2
    with pytest.raises(ValueError, match="unique"):
        reaction_center("[CH3:1][CH3:1]>>[CH3:1][CH3:2]")
    with pytest.raises(ValueError, match="bijectively"):
        reaction_center("[H:1][H:2]>>[H:1][2H:2]")
