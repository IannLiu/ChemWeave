"""Local browser interface; web dependencies are optional for library users."""

from dataclasses import dataclass, field
from contextlib import asynccontextmanager
from functools import lru_cache
from io import BytesIO
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import Event, Lock, Thread, Timer
from time import monotonic
import webbrowser
from zipfile import ZipFile

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import uvicorn

from . import __version__
from .generator import GenerationConfig, Generator
from .models import Network, Species
from .pathways import find_pathways
from .rules import ReactionRule, load_rules


class SpeciesInput(BaseModel):
    smiles: str
    name: str = ""
    phase: str = "gas"
    site: str | None = None


class RunInput(BaseModel):
    seeds: list[SpeciesInput] = Field(min_length=1)
    targets: list[SpeciesInput] = Field(default_factory=list)
    rules: list[dict] | None = Field(default=None, min_length=1)
    config: GenerationConfig = Field(default_factory=lambda: GenerationConfig(max_carbon=6, max_depth=4))


class PathwayInput(BaseModel):
    target: str
    max_reactions: int = Field(default=8, ge=0)
    max_pathways: int = Field(default=100, ge=1)


@dataclass
class Workspace:
    lock: Lock = field(default_factory=Lock)
    cancel: Event = field(default_factory=Event)
    status: str = "idle"
    summary: dict = field(default_factory=dict)
    network: Network | None = None
    request: dict | None = None
    error: str | None = None
    analysis: dict | None = None
    started: float = 0
    elapsed: float = 0
    revision: int = 0
    worker: Thread | None = None


@lru_cache(maxsize=512)
def molecule_svg(smiles, mapped=False):
    """Draw actual input/reaction maps without inventing new display numbers."""
    params = Chem.SmilesParserParams()
    params.removeHs = not mapped
    mol = Chem.MolFromSmiles(smiles, params)
    if mol is None or not mol.GetNumAtoms():
        raise ValueError(f"Invalid or empty molecular SMILES: {smiles!r}")
    if not mapped:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        mol = Chem.RemoveHs(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(400, 210)
    drawer.drawOptions().clearBackground = False
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def create_app(browser_url=None):
    workspace = Workspace()

    @asynccontextmanager
    async def lifespan(app):
        if browser_url:
            Timer(0.5, webbrowser.open, args=[browser_url]).start()
        yield
        workspace.cancel.set()

    app = FastAPI(title="ChemWeave", version=__version__, docs_url=None,
                  redoc_url=None, lifespan=lifespan)
    app.state.workspace = workspace
    static = Path(__file__).with_name("static")
    app.mount("/static", StaticFiles(directory=static), name="static")

    @app.get("/")
    def index():
        return FileResponse(static / "index.html")

    @app.get("/api/rules")
    def rules():
        return {"version": __version__, "rules": [r.to_dict() for r in load_rules()]}

    @app.post("/api/preview")
    def preview(record: SpeciesInput):
        try:
            species = Species.from_smiles(**record.model_dump())
            return {"smiles": species.smiles, "formula": species.formula,
                    "svg": molecule_svg(species.smiles)}
        except ValueError as error:
            raise HTTPException(422, str(error)) from error

    @app.get("/api/molecule")
    def molecule(smiles: str, mapped: bool = False):
        try:
            return Response(molecule_svg(smiles, mapped), media_type="image/svg+xml")
        except ValueError as error:
            raise HTTPException(422, str(error)) from error

    @app.post("/api/run", status_code=202)
    def run(request: RunInput):
        try:
            seeds = [Species.from_smiles(**s.model_dump()) for s in request.seeds]
            targets = [Species.from_smiles(**s.model_dump()) for s in request.targets]
            selected_rules = load_rules() if request.rules is None else [ReactionRule(**r) for r in request.rules]
            generator = Generator(selected_rules, request.config)
        except (ValueError, TypeError) as error:
            raise HTTPException(422, str(error)) from error
        with workspace.lock:
            if workspace.status in {"running", "analyzing", "stopping"}:
                raise HTTPException(409, "A calculation is already running")
            workspace.cancel.clear()
            workspace.status = "running"
            workspace.network = None
            workspace.analysis = None
            workspace.summary = {}
            workspace.error = None
            workspace.started = monotonic()
            workspace.elapsed = 0
            workspace.request = {**request.model_dump(mode="json"), "rules": [r.to_dict() for r in selected_rules]}
            workspace.revision += 1

        def work():
            last_report = 0

            def report(network):
                nonlocal last_report
                now = monotonic()
                if now - last_report < 0.25:
                    return
                snapshot = network.summary()
                with workspace.lock:
                    workspace.summary = snapshot
                last_report = now

            try:
                network = generator.generate(seeds, targets=targets, progress=report,
                                             stop_requested=workspace.cancel.is_set)
                with workspace.lock:
                    workspace.network = network
                    workspace.summary = network.summary()
                    workspace.status = "finished"
                    workspace.elapsed = monotonic() - workspace.started
                    workspace.revision += 1
            except Exception as error:
                # Keep a background failure visible to the browser instead of
                # leaving it indefinitely in the running state.
                with workspace.lock:
                    workspace.error = str(error)
                    workspace.status = "error"
                    workspace.elapsed = monotonic() - workspace.started
                    workspace.revision += 1

        workspace.worker = Thread(target=work, daemon=True, name="chemweave-generation")
        workspace.worker.start()
        return {"status": "running"}

    @app.get("/api/status")
    def status():
        with workspace.lock:
            elapsed = monotonic() - workspace.started if workspace.status in {"running", "analyzing", "stopping"} else workspace.elapsed
            return {"status": workspace.status, "summary": workspace.summary,
                    "error": workspace.error, "elapsed": round(elapsed, 2),
                    "revision": workspace.revision,
                    "analysis_termination": workspace.analysis["termination"] if workspace.analysis else None}

    @app.post("/api/stop")
    def stop():
        with workspace.lock:
            if workspace.status in {"running", "analyzing"}:
                workspace.status = "stopping"
                workspace.cancel.set()
            return {"status": workspace.status}

    @app.get("/api/network")
    def network():
        with workspace.lock:
            if workspace.network is None:
                raise HTTPException(409, "No completed network is available")
            return {"network": workspace.network.to_dict(), "request": workspace.request,
                    "analysis": workspace.analysis}

    @app.post("/api/pathways", status_code=202)
    def analyze(request: PathwayInput):
        with workspace.lock:
            if workspace.status in {"running", "analyzing", "stopping"}:
                raise HTTPException(409, "A calculation is already running")
            network = workspace.network
            if network is None:
                raise HTTPException(409, "Generate or open a network first")
            if request.target not in network.species:
                raise HTTPException(422, "Select a target species present in the network")
            workspace.cancel.clear()
            workspace.analysis = None
            workspace.status = "analyzing"
            workspace.error = None
            workspace.started = monotonic()
            workspace.elapsed = 0
            workspace.revision += 1

        def work():
            try:
                result = find_pathways(network, **request.model_dump(),
                                       stop_requested=workspace.cancel.is_set)
                with workspace.lock:
                    workspace.analysis = result
                    workspace.status = "finished"
                    workspace.elapsed = monotonic() - workspace.started
                    workspace.revision += 1
            except Exception as error:
                with workspace.lock:
                    workspace.error = str(error)
                    workspace.status = "error"
                    workspace.elapsed = monotonic() - workspace.started
                    workspace.revision += 1

        workspace.worker = Thread(target=work, daemon=True, name="chemweave-pathways")
        workspace.worker.start()
        return {"status": "analyzing"}

    @app.post("/api/open")
    def open_project(data: dict):
        try:
            network = Network.from_dict(data.get("network", data))
            request = data.get("request")
            if request is None:
                request = {
                    "seeds": [
                        {"smiles": r["mapped_smiles"], "name": r.get("name", ""),
                         "phase": network.species[r["species_id"]].phase,
                         "site": network.species[r["species_id"]].site}
                        for r in network.seed_inputs
                    ] or [
                        {"smiles": network.species[sid].mapped_smiles,
                         "name": network.species[sid].name, "phase": network.species[sid].phase,
                         "site": network.species[sid].site} for sid in network.seeds
                    ],
                    "targets": [{k: t.get(k) for k in ("smiles", "name", "phase", "site")} for t in network.targets],
                    "rules": network.rules, "config": network.config,
                }
            request = RunInput.model_validate(request).model_dump(mode="json")
        except (ValueError, TypeError, KeyError) as error:
            raise HTTPException(422, str(error)) from error
        with workspace.lock:
            if workspace.status in {"running", "analyzing", "stopping"}:
                raise HTTPException(409, "Stop the current calculation before opening a file")
            workspace.network = network
            workspace.analysis = data.get("analysis")
            workspace.request = request
            workspace.status = "finished"
            workspace.error = None
            workspace.summary = network.summary()
            workspace.elapsed = 0
            workspace.revision += 1
        return {"network": network.to_dict(), "request": request, "analysis": workspace.analysis}

    @app.post("/api/new")
    def new_project():
        with workspace.lock:
            if workspace.status in {"running", "analyzing", "stopping"}:
                raise HTTPException(409, "Stop the current calculation before starting a new project")
            workspace.network = None
            workspace.analysis = None
            workspace.request = None
            workspace.summary = {}
            workspace.error = None
            workspace.status = "idle"
            workspace.elapsed = 0
            workspace.revision += 1
        return {"status": "idle"}

    @app.get("/api/export/{format}")
    def export(format: str):
        with workspace.lock:
            network, request, analysis = workspace.network, workspace.request, workspace.analysis
        if network is None:
            raise HTTPException(409, "No completed network is available")
        if format in {"project", "json", "pathways"}:
            data = network.to_dict()
            if format == "project":
                data = {"chemweave_project": 1, "request": request, "network": data, "analysis": analysis}
            elif format == "pathways":
                if analysis is None:
                    raise HTTPException(409, "Analyze formation pathways first")
                data = {"analysis": analysis, "network": data}
            body = json.dumps(data, ensure_ascii=False, indent=2).encode()
            return Response(body, media_type="application/json", headers={
                "Content-Disposition": f'attachment; filename="chemweave-{format}.json"',
            })
        if format == "csv":
            with TemporaryDirectory() as directory:
                network.export_csv(directory)
                network.export_dot(Path(directory) / "network.dot")
                buffer = BytesIO()
                with ZipFile(buffer, "w") as archive:
                    for path in Path(directory).iterdir():
                        archive.write(path, path.name)
            return Response(buffer.getvalue(), media_type="application/zip", headers={
                "Content-Disposition": 'attachment; filename="chemweave-tables.zip"',
            })
        raise HTTPException(404, "Unknown export format")

    return app


def serve(port=8765, open_browser=True):
    """Bind locally and open the interface once the server has started."""
    app = create_app(browser_url=f"http://127.0.0.1:{port}" if open_browser else None)
    uvicorn.run(app, host="127.0.0.1", port=port, log_level="info")
