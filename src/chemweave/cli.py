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
    web = subparsers.add_parser("web", help="Open the local browser interface")
    web.add_argument("--port", type=int, default=8765)
    web.add_argument("--no-browser", action="store_true", help="Do not open a browser automatically")
    args = parser.parse_args()

    if args.command == "web":
        try:
            from .web import serve
        except ModuleNotFoundError as error:
            if error.name not in {"fastapi", "uvicorn", "starlette", "pydantic"}:
                raise
            parser.error('Install the web interface first: pip install "chemweave[web]"')
        serve(port=args.port, open_browser=not args.no_browser)
        return 0

    data = tomllib.loads(args.config.read_text(encoding="utf-8"))
    rule_path = data.get("rules")
    if rule_path is not None:
        rule_path = args.config.parent / rule_path
    network = Generator(
        rules=load_rules(rule_path), config=GenerationConfig(**data.get("generation", {})),
    ).generate(
        [Species.from_smiles(**record) for record in data["seeds"]],
        targets=[Species.from_smiles(**record) for record in data.get("targets", [])],
    )
    network.export_json(args.output / "network.json")
    network.export_csv(args.output)
    network.export_dot(args.output / "network.dot")
    print(json.dumps(network.summary(), ensure_ascii=False, indent=2))
    complete = network.termination in ("exhausted", "max_depth", "targets_reached")
    missing_targets = any(t["species_id"] is None for t in network.targets)
    return 0 if complete and not missing_targets else 2
