"""Atom inventories, external numbering, and whole-reaction equivalence."""

from collections import Counter
import json

import pytest
from rdkit import Chem

from chemweave import GenerationConfig, Generator, Network, ReactionRule, Species, load_rules
from chemweave.mapping import prepare_reactants, reaction_symmetry_key


@pytest.fixture
def dehydrogenation():
    return next(rule for rule in load_rules() if rule.id == "alkane_dehydrogenation")


@pytest.mark.parametrize("smiles", [
    "CCC", "[CH3:101][CH2:205][CH3:309]", "[CH3:309][CH2:205][CH3:101]",
    "[H:90][C:101]([H])([H])[C:205]([H])([H])[C:309]([H])([H])[H]",
])
def test_plain_and_mapped_inputs_generate_same_structures(smiles, dehydrogenation):
    net = Generator([dehydrogenation]).generate([Species.from_smiles(smiles)])
    assert {s.smiles for s in net.species.values()} == {"CCC", "C=CC", "[H][H]"}
    assert len(net.reactions) == 1
    assert len(net.reactions[0].channels) == 1
    assert net.stats["rejected"] == {}


def test_partial_input_maps_preserve_existing_labels(dehydrogenation):
    s = Species.from_smiles("[CH3:101][CH2:205][CH3:309]")
    assert {number for number, original in s.input_maps if original} == {101, 205, 309}
    assert all(number > 309 for number, original in s.input_maps if not original)
    net = Generator([dehydrogenation]).generate([s])
    for match in net.reactions[0].channels:
        original = {row["input_map"]: row["reaction_map"] for row in match["atom_mapping"] if row["input_map"]}
        assert original == {101: 101, 205: 205, 309: 309}
        # All hydrogens in H2 originate at the input, rather than being invented.
        h2 = next(smi for smi in match["mapped_product_smiles"]
                  if Species.from_smiles(smi).formula == "H2")
        params = Chem.SmilesParserParams()
        params.removeHs = False
        assert all(a.GetAtomMapNum() > 309 for a in Chem.MolFromSmiles(h2, params).GetAtoms())


def test_duplicate_input_labels_are_recorded_and_disambiguated():
    s = Species.from_smiles("[CH3:7][CH2:7][CH3:9]")
    assert s.key == Species.from_smiles("CCC").key
    assert sum(original == 7 for _, original in s.input_maps) == 2
    assert len({number for number, _ in s.input_maps}) == 11
    assert (7, 7) in s.input_maps and (9, 9) in s.input_maps


def test_self_coupling_numbers_each_occurrence_separately():
    rule = next(r for r in load_rules() if r.id == "terminal_olefin_coupling")
    s = Species.from_smiles("[CH2:11]=[CH:12][CH3:13]")
    net = Generator([rule], GenerationConfig(max_carbon=6)).generate([s])
    assert len(net.reactions) == 1
    assert net.reactions[0].reactants == ((net.seeds[0], 2),)
    for match in net.reactions[0].channels:
        origins = match["atom_mapping"]
        left = {o["reaction_map"] for o in origins if o["reactant"] == 0}
        right = {o["reaction_map"] for o in origins if o["reactant"] == 1}
        assert left.isdisjoint(right)
        assert len(left) == len(right) == 9
        assert all(o["reaction_map"] == o["species_map"] for o in origins if o["reactant"] == 0)


def test_complete_reaction_maps_are_bijective_in_demo():
    net = Generator(load_rules(), GenerationConfig(max_carbon=6)).generate([Species.from_smiles("CCC")])
    params = Chem.SmilesParserParams()
    params.removeHs = False
    for reaction in net.reactions:
        assert len({c["symmetry_key"] for c in reaction.channels}) == len(reaction.channels)
        assert all(c["match_count"] >= 1 for c in reaction.channels)
        for match in reaction.channels:
            before, after = match["mapped_reaction_smiles"].split(">>")
            atoms = []
            for side in (before, after):
                mol = Chem.MolFromSmiles(side, params)
                inventory = {a.GetAtomMapNum(): (a.GetAtomicNum(), a.GetIsotope()) for a in mol.GetAtoms()}
                assert len(inventory) == mol.GetNumAtoms()
                assert 0 not in inventory
                atoms.append(inventory)
            assert atoms[0] == atoms[1]


def test_symmetry_key_ignores_numbering_but_keeps_stereo(dehydrogenation):
    keys = []
    for smiles in ("CCC", "[CH3:101][CH2:205][CH3:309]", "[CH3:309][CH2:205][CH3:101]"):
        net = Generator([dehydrogenation]).generate([Species.from_smiles(smiles)])
        keys.append(net.reactions[0].channels[0]["symmetry_key"])
    assert len(set(keys)) == 1
    trans = Species.from_smiles("C/C=C/C")
    cis = Species.from_smiles("C/C=C\\C")
    assert trans.key != cis.key
    mols, _ = prepare_reactants([trans])
    same, _ = prepare_reactants([trans])
    changed, _ = prepare_reactants([cis])
    states = [("gas", None)]
    assert reaction_symmetry_key(mols, same, states, states) != reaction_symmetry_key(mols, changed, states, states)


def test_neopentane_symmetry_is_a_partition_without_self_match_enumeration(monkeypatch):
    def forbidden(*args, **kwargs):
        raise AssertionError("Symmetry classification must not enumerate self matches")
    monkeypatch.setattr(Chem.Mol, "GetSubstructMatches", forbidden)
    s = Species.from_smiles("CC(C)(C)C")
    groups = Counter(rank for _, rank in s.symmetry_classes)
    assert sorted(groups.values()) == [1, 4, 12]
    assert len(s.symmetry_classes) == 17
    assert len({number for number, _ in s.symmetry_classes}) == 17


def test_joint_reaction_symmetry_is_not_independent_atom_class_relabeling():
    s = Species.from_smiles("c1ccccc1")
    mols, _ = prepare_reactants([s])
    mol = mols[0]
    carbons = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 6]
    keys = []
    for distance in (1, 2):
        product = Chem.Mol(mol)
        a, b = carbons[0], carbons[distance]
        pairs = [(a.GetIdx(), b.GetIdx())]
        pairs.extend(zip([n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() == 1],
                         [n.GetIdx() for n in b.GetNeighbors() if n.GetAtomicNum() == 1]))
        for ai, bi in pairs:
            amap = mol.GetAtomWithIdx(ai).GetAtomMapNum()
            bmap = mol.GetAtomWithIdx(bi).GetAtomMapNum()
            product.GetAtomWithIdx(ai).SetAtomMapNum(bmap)
            product.GetAtomWithIdx(bi).SetAtomMapNum(amap)
        keys.append(reaction_symmetry_key(mols, [product], [("gas", None)], [("gas", None)]))
    assert keys[0] != keys[1]


def test_isotope_hydrogen_is_preserved(dehydrogenation):
    net = Generator([dehydrogenation]).generate([Species.from_smiles("[2H:99]CC")])
    assert len(net.reactions) >= 1
    params = Chem.SmilesParserParams()
    params.removeHs = False
    for r in net.reactions:
        for match in r.channels:
            mol = Chem.MolFromSmiles(".".join(match["mapped_product_smiles"]), params)
            deuterium = [a for a in mol.GetAtoms() if a.GetIsotope() == 2]
            assert len(deuterium) == 1 and deuterium[0].GetAtomMapNum() == 99


def test_seed_aliases_and_mapping_survive_json(tmp_path):
    seeds = [Species.from_smiles("[CH3:10][CH2:20][CH3:30]"), Species.from_smiles("[CH3:1][CH2:2][CH3:3]")]
    net = Generator([]).generate(seeds)
    assert len(net.species) == len(net.seeds) == 1
    assert len(net.seed_inputs) == 2
    path = tmp_path / "network.json"
    net.export_json(path)
    restored = Network.from_json(path)
    second = tmp_path / "restored.json"
    restored.export_json(second)
    assert json.loads(path.read_text()) == json.loads(second.read_text())


def test_identical_scission_products_keep_separate_atom_occurrences():
    rule = next(r for r in load_rules() if r.id == "olefin_scission")
    net = Generator([rule]).generate([Species.from_smiles("CCC=C")])
    assert len(net.reactions) == 1
    product_id = next(sid for sid, s in net.species.items() if s.smiles == "C=C")
    assert net.reactions[0].products == ((product_id, 2),)
    params = Chem.SmilesParserParams()
    params.removeHs = False
    for match in net.reactions[0].channels:
        labels = [{a.GetAtomMapNum() for a in Chem.MolFromSmiles(smi, params).GetAtoms()}
                  for smi in match["mapped_product_smiles"]]
        assert labels[0].isdisjoint(labels[1])


def test_version_one_json_remains_readable(tmp_path):
    data = {
        "schema_version": 1,
        "species": [{"id": "S00001", "smiles": "CCC", "formula": "C3H8",
                     "phase": "gas", "site": None, "name": "propane"}],
        "reactions": [], "seeds": ["S00001"], "rules": [], "config": {},
        "stats": {}, "termination": "exhausted",
    }
    path = tmp_path / "old_network.json"
    path.write_text(json.dumps(data))
    net = Network.from_json(path)
    assert net.species["S00001"].smiles == "CCC"
    assert net.seed_inputs == []
    assert all(original == 0 for _, original in net.species["S00001"].input_maps)


@pytest.mark.parametrize("smiles,products,match_count", [
    ("CCC", {"C=CC"}, 24),
    ("CCCC", {"C=CCC", "CC=CC"}, 32),
])
def test_compact_channels_keep_one_mapping_and_only_increment_count(smiles, products, match_count, dehydrogenation, tmp_path):
    net = Generator([dehydrogenation]).generate([Species.from_smiles(smiles)])
    assert len(net.reactions) == len(products)
    observed = set()
    for reaction in net.reactions:
        assert not hasattr(reaction, "matches")
        assert len(reaction.channels) == 1
        channel = reaction.channels[0]
        assert "match_indices" not in channel
        assert channel["atom_mapping"] and channel["product_atom_properties"]
        observed.update(net.species[sid].smiles for sid, _ in reaction.products if net.species[sid].carbon_count)
    assert observed == products
    assert net.summary()["template_matches"] == match_count
    path = tmp_path / "compact.json"
    net.export_json(path)
    data = json.loads(path.read_text())
    assert data["schema_version"] == 3
    assert all("matches" not in r for r in data["reactions"])
    assert sum(len(r["channels"]) for r in data["reactions"]) == len(products)


def test_version_two_json_compacts_redundant_records_on_load(dehydrogenation, tmp_path):
    net = Generator([dehydrogenation]).generate([Species.from_smiles("CCC")])
    path = tmp_path / "network.json"
    net.export_json(path)
    data = json.loads(path.read_text())
    data["schema_version"] = 2
    reaction = data["reactions"][0]
    representative = reaction["channels"][0]
    count = representative.pop("match_count")
    reaction["matches"] = [dict(representative) for _ in range(count)]
    reaction["channels"] = [{
        "symmetry_key": representative["symmetry_key"],
        "mapped_reaction_smiles": representative["mapped_reaction_smiles"],
        "match_indices": list(range(count)),
    }]
    path.write_text(json.dumps(data))
    restored = Network.from_json(path)
    assert not hasattr(restored.reactions[0], "matches")
    assert len(restored.reactions[0].channels) == 1
    assert restored.reactions[0].channels[0]["match_count"] == 24
    assert restored.summary() == net.summary()
    restored.export_json(path)
    compact = json.loads(path.read_text())
    assert compact["schema_version"] == 3
    assert "matches" not in compact["reactions"][0]
