import json
import os
from pathlib import Path
import subprocess
import sys


def test_cli_custom_rules_are_relative_to_config(tmp_path):
    rules = [{
        "id": "dehydrogenation",
        "smarts": "[CH3:1][CH2:2][CH3:3]>>[CH2:1]=[CH:2][CH3:3].[H][H]",
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

