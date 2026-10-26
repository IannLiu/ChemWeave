import json
import os
from pathlib import Path
import subprocess
import sys

import pytest


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
