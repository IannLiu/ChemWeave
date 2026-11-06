import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


def test_cli_single_event_generation_and_saved_annotation(tmp_path):
    root = Path(__file__).resolve().parents[1]
    env = {**os.environ, "PYTHONPATH": str(root / "src")}
    output = tmp_path / "generated"
    command = [sys.executable, "-m", "chemweave"]
    result = subprocess.run(command + ["generate", "--config", str(root / "examples/single_event.toml"),
                                      "--output", str(output)], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    data = json.loads((output / "network.json").read_text())
    assert data["reactions"][0]["channels"][0]["single_event"]["factor"] == 4
    options = tmp_path / "options.json"
    options.write_text(json.dumps({"calculator": None, "factors": [{
        "reaction_id": "R00001", "channel": 1, "factor": 4,
        "kinetic_subclass": "methane", "parameter_class": "methane-fit",
        "source": "Independent symmetry calculation", "convention": "sigma_R / sigma_TS = 12 / 3",
    }]}))
    result = subprocess.run(command + ["annotate", "--network", str(output / "network.json"),
                                      "--annotations", str(options), "--output", str(tmp_path / "annotated")],
                            env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    saved = json.loads((tmp_path / "annotated/network.json").read_text())
    assert saved["reactions"][0]["channels"][0]["single_event"]["status"] == "imported"
    assert (tmp_path / "annotated/kinetic_parameters.csv").exists()


def test_cli_custom_rules_are_relative_to_config(tmp_path):
    rules = [{
        "id": "dehydrogenation",
        "smarts": "[C;H3:1](-[H:4])-[C;H2:2](-[H:5])-[C;H3:3]>>[C:1]=[C:2]-[C:3].[H:4]-[H:5]",
    }]
    (tmp_path / "rules.json").write_text(json.dumps(rules))
    config = tmp_path / "input.toml"
    config.write_text('rules = "rules.json"\n[[seeds]]\nsmiles = "CCC"\n')
    source = str(Path(__file__).resolve().parents[1] / "src")
    completed = subprocess.run(
        [sys.executable, "-m", "chemweave", "generate", "--config", str(config),
         "--output", str(tmp_path / "output")],
        cwd=tmp_path, env={**os.environ, "PYTHONPATH": source}, text=True, capture_output=True,
    )
    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout)["reactions"] == 1
    assert (tmp_path / "output/network.json").exists()


def test_cli_budget_stop_is_reported(tmp_path):
    config = tmp_path / "input.toml"
    config.write_text('[[seeds]]\nsmiles = "CCC"\n[generation]\nmax_species = 1\n')
    source = str(Path(__file__).resolve().parents[1] / "src")
    completed = subprocess.run(
        [sys.executable, "-m", "chemweave", "generate", "--config", str(config),
         "--output", str(tmp_path / "output")],
        cwd=tmp_path, env={**os.environ, "PYTHONPATH": source}, text=True, capture_output=True,
    )
    assert completed.returncode == 2
    assert json.loads(completed.stdout)["termination"] == "max_species"
    assert (tmp_path / "output/network.json").exists()



@pytest.mark.parametrize("settings,targets,termination,exit_code", [
    ('max_depth = 1', '', 'max_depth', 0),
    ('max_depth = 1', '[[targets]]\nsmiles = "C=CC"\n[[targets]]\nsmiles = "[H][H]"\n', 'targets_reached', 0),
    ('max_depth = 1', '[[targets]]\nsmiles = "C=C"\n', 'max_depth', 2),
    ('max_carbon = 3', '[[targets]]\nsmiles = "C=C"\n', 'exhausted', 2),
])
def test_cli_depth_and_target_conditions(tmp_path, settings, targets, termination, exit_code):
    config = tmp_path / "input.toml"
    config.write_text('[[seeds]]\nsmiles = "CCC"\n[generation]\n' + settings + '\n' + targets)
    source = str(Path(__file__).resolve().parents[1] / "src")
    completed = subprocess.run(
        [sys.executable, "-m", "chemweave", "generate", "--config", str(config),
         "--output", str(tmp_path / "output")],
        cwd=tmp_path, env={**os.environ, "PYTHONPATH": source}, text=True, capture_output=True,
    )
    assert completed.returncode == exit_code, completed.stderr
    assert json.loads(completed.stdout)["termination"] == termination
    data = json.loads((tmp_path / "output/network.json").read_text())
    assert all(record["depth"] is not None for record in data["species"])
    assert len(data["targets"]) == targets.count('[[targets]]')


def test_cli_pathways_use_saved_network_and_report_limits(tmp_path):
    from chemweave import GenerationConfig, Generator, Species, load_rules
    net = Generator(load_rules(), GenerationConfig(max_carbon=6, max_depth=4,
                    stop_when_targets_found=False)).generate(
                        [Species.from_smiles('CCC')], targets=[Species.from_smiles('C=C')])
    net.export_json(tmp_path / 'network.json')
    source = str(Path(__file__).resolve().parents[1] / 'src')
    for count, code in [(100, 0), (1, 2)]:
        completed = subprocess.run(
            [sys.executable, '-m', 'chemweave', 'pathways', '--network', str(tmp_path/'network.json'),
             '--target', net.targets[0]['species_id'], '--max-pathways', str(count),
             '--output', str(tmp_path/'pathways.json')],
            cwd=tmp_path, env={**os.environ, 'PYTHONPATH': source}, text=True, capture_output=True)
        assert completed.returncode == code, completed.stderr
        report = json.loads((tmp_path/'pathways.json').read_text())
        assert len(report['analysis']['pathways']) == min(count, 2)
        assert report['network']['reactions']
