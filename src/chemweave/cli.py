"""Command-line entry point for reproducible TOML-configured runs."""

import argparse
import json
from pathlib import Path
import tomllib

from .generator import GenerationConfig, Generator
from .models import Species
from .rules import load_rules


def main():
    parser = argparse.ArgumentParser(description="Generate a rule-based chemical reaction network")
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser("generate", help="Generate a network from a TOML configuration")
    generate.add_argument("--config", type=Path, required=True)
    generate.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()

    data = tomllib.loads(args.config.read_text(encoding="utf-8"))
    rule_path = data.get("rules")
    if rule_path is not None:
        rule_path = args.config.parent / rule_path
    network = Generator(
        rules=load_rules(rule_path), config=GenerationConfig(**data.get("generation", {})),
    ).generate([Species.from_smiles(**record) for record in data["seeds"]])
    network.export_json(args.output / "network.json")
    network.export_csv(args.output)
    network.export_dot(args.output / "network.dot")
    print(json.dumps(network.summary(), ensure_ascii=False, indent=2))
    return 0 if network.termination == "exhausted" else 2

