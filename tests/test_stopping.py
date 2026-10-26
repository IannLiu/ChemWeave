"""Depth bounds and identity-based target stopping."""

import csv
import json

import pytest

from chemweave import GenerationConfig, Generator, Network, Species, load_rules


@pytest.mark.parametrize("depth,formulas", [
    (0, {"C3H8"}),
    (1, {"C3H8", "C3H6", "H2"}),
    (2, {"C3H8", "C3H6", "H2", "C6H12"}),
])
def test_depth_bounds_retain_boundary_products(depth, formulas):
    net = Generator(load_rules(), GenerationConfig(max_carbon=6, max_depth=depth)).generate(
        [Species.from_smiles("CCC")]
    )
    assert net.termination == "max_depth"
    assert {s.formula for s in net.species.values()} == formulas
    assert max(net.species_depths.values()) == depth
    for reaction in net.reactions:
        expected = 1 + max(net.species_depths[sid] for sid, _ in reaction.reactants)
        assert reaction.depth == expected <= depth
    if depth == 0:
        assert net.stats["applications"] == 0


def test_first_discovery_depths_are_minimal_and_ignore_seed_order():
    seeds = [Species.from_smiles(s) for s in ("CCC", "C=C")]
    results = []
    for inputs in (seeds, seeds[::-1]):
        net = Generator(load_rules(), GenerationConfig(max_carbon=6)).generate(inputs)
        # Independently relax the complete reaction graph to shortest depths.
        shortest = {sid: (0 if sid in net.seeds else float("inf")) for sid in net.species}
        changed = True
        while changed:
            changed = False
            for reaction in net.reactions:
                depth = 1 + max(shortest[sid] for sid, _ in reaction.reactants)
                for sid, _ in reaction.products:
                    if depth < shortest[sid]:
                        shortest[sid] = depth
                        changed = True
        assert net.species_depths == shortest
        for reaction in net.reactions:
            assert reaction.depth == 1 + max(shortest[sid] for sid, _ in reaction.reactants)
        results.append({s.key: shortest[sid] for sid, s in net.species.items()})
    assert results[0] == results[1]


def test_isomerization_cycle_does_not_inflate_species_depth():
    rule = next(r for r in load_rules() if r.id == "olefin_double_bond_shift")
    net = Generator([rule], GenerationConfig(max_depth=5)).generate([Species.from_smiles("C=CCC")])
    assert net.termination == "exhausted"
    assert {s.smiles: net.species_depths[sid] for sid, s in net.species.items()} == {
        "C=CCC": 0, "CC=CC": 1,
    }
    assert [r.depth for r in net.reactions] == [1, 2]


def test_all_targets_must_appear_and_last_reaction_is_complete():
    targets = [Species.from_smiles(s) for s in ("C=CC", "C=C")]
    net = Generator(load_rules(), GenerationConfig(max_carbon=6)).generate(
        [Species.from_smiles("CCC")], targets=targets,
    )
    assert net.termination == "targets_reached"
    assert net.summary()["targets_found"] == net.summary()["targets_total"] == 2
    assert net.summary()["missing_targets"] == []
    assert [net.species_depths[t["species_id"]] for t in net.targets] == [1, 4]
    last = net.reactions[-1]
    assert len(last.products) == 2  # Keep the coproduct of the target-forming step.
    assert last.channels[0]["mapped_reaction_smiles"]
    assert all(sid in net.species for sid, _ in last.reactants + last.products)
    referenced = set(net.seeds)
    for reaction in net.reactions:
        referenced.update(sid for sid, _ in reaction.reactants + reaction.products)
    assert referenced == set(net.species)


def test_seed_targets_ignore_map_labels_and_duplicate_requests():
    net = Generator(load_rules(), GenerationConfig(max_depth=0)).generate(
        [Species.from_smiles("CCC"), Species.from_smiles("[H][H]")],
        targets=[Species.from_smiles("[CH3:9][CH2:7][CH3:3]"), Species.from_smiles("CCC")],
    )
    assert net.termination == "targets_reached"
    assert net.summary()["targets_total"] == 1
    assert net.targets[0]["species_id"] == net.seeds[0]
    assert net.stats["applications"] == 0
    assert len(net.species) == 2  # All requested seeds are retained.


def test_targets_match_phase_and_site():
    target = Species.from_smiles("CCC", phase="surface", site="BAS")
    net = Generator([]).generate([Species.from_smiles("CCC")], targets=[target])
    assert net.termination == "exhausted"
    assert net.summary()["targets_found"] == 0
    assert net.summary()["missing_targets"][0]["site"] == "BAS"


@pytest.mark.parametrize("depth,termination", [(0, "max_depth"), (1, "max_depth"), (None, "exhausted")])
def test_unreachable_targets_are_reported(depth, termination):
    net = Generator(load_rules(), GenerationConfig(max_carbon=3, max_depth=depth)).generate(
        [Species.from_smiles("CCC")], targets=[Species.from_smiles("C=C")],
    )
    assert net.termination == termination
    assert net.targets[0]["species_id"] is None
    assert net.summary()["targets_found"] == 0
    assert net.summary()["missing_targets"][0]["smiles"] == "C=C"


def test_species_budget_prevents_partial_target_reaction():
    net = Generator(load_rules(), GenerationConfig(max_species=2)).generate(
        [Species.from_smiles("CCC")], targets=[Species.from_smiles("C=CC")],
    )
    assert net.termination == "max_species"
    assert len(net.species) == 1
    assert not net.reactions
    assert net.targets[0]["species_id"] is None


def test_depth_and_target_records_roundtrip_and_csv(tmp_path):
    net = Generator(load_rules(), GenerationConfig(max_carbon=6, max_depth=4)).generate(
        [Species.from_smiles("CCC")], targets=[Species.from_smiles("C=C", name="ethene")],
    )
    path = tmp_path / "network.json"
    net.export_json(path)
    restored = Network.from_json(path)
    assert restored.species_depths == net.species_depths
    assert restored.targets == net.targets
    assert restored.summary() == net.summary()
    assert [r.depth for r in restored.reactions] == [r.depth for r in net.reactions]
    restored.export_csv(tmp_path)
    with (tmp_path / "species.csv").open(encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    assert {r["id"]: int(r["depth"]) for r in rows} == net.species_depths
    with (tmp_path / "reactions.csv").open(encoding="utf-8-sig") as f:
        assert [int(r["depth"]) for r in csv.DictReader(f)] == [r.depth for r in net.reactions]


def test_version_three_network_loads_without_inventing_depths(tmp_path):
    net = Generator(load_rules(), GenerationConfig(max_depth=1)).generate([Species.from_smiles("CCC")])
    path = tmp_path / "old.json"
    net.export_json(path)
    data = json.loads(path.read_text())
    data["schema_version"] = 3
    data.pop("targets")
    for record in data["species"] + data["reactions"]:
        record.pop("depth")
    path.write_text(json.dumps(data))
    restored = Network.from_json(path)
    assert restored.species_depths == {}
    assert restored.targets == []
    assert all(r.depth is None for r in restored.reactions)


def test_progress_and_cooperative_stop_preserve_completed_reaction():
    from threading import Event

    cancel = Event()
    snapshots = []

    def report(network):
        snapshots.append(network.summary())
        if network.reactions:
            cancel.set()

    net = Generator(load_rules(), GenerationConfig(max_carbon=6)).generate(
        [Species.from_smiles('CCC')], progress=report, stop_requested=cancel.is_set,
    )
    assert net.termination == 'user_stopped'
    assert snapshots[0]['species'] == 1
    assert snapshots[-1]['reactions'] == 1
    assert len(net.species) == 3
    assert net.reactions[0].channels[0]['match_count'] == 24
    assert Network.from_dict(net.to_dict()).to_dict() == net.to_dict()
