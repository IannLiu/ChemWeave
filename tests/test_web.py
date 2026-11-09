"""Real chemistry through the local UI API, including cancellation and persistence."""

from io import BytesIO
import json
from threading import Event
from time import monotonic, sleep
from zipfile import ZipFile

import pytest

pytest.importorskip('fastapi')
pytest.importorskip('httpx')
from fastapi.testclient import TestClient

from chemweave import GenerationConfig, Generator, Network, Species, load_rules
from chemweave.web import create_app


@pytest.fixture
def client():
    with TestClient(create_app()) as instance:
        yield instance


def finish(client):
    deadline = monotonic() + 10
    while monotonic() < deadline:
        state = client.get('/api/status').json()
        if state['status'] not in {'running', 'analyzing', 'stopping'}:
            return state
        sleep(.02)
    pytest.fail('Background generation did not finish')


def test_assets_and_molecule_previews(client):
    assert client.get('/').status_code == 200
    assert 'ChemWeave' in client.get('/').text
    for path in ['app.js', 'style.css', 'vendor/cytoscape.min.js']:
        assert client.get('/static/' + path).status_code == 200
    assert len(client.get('/api/rules').json()['rules']) == 4
    preview = client.post('/api/preview', json={'smiles': '[CH3:15][CH2:8][CH3:7]'})
    assert preview.status_code == 200
    assert preview.json()['smiles'] == 'CCC'
    assert preview.json()['formula'] == 'C3H8'
    assert '<svg' in preview.json()['svg']
    assert client.post('/api/preview', json={'smiles': 'not a molecule'}).status_code == 422
    mapped = client.get('/api/molecule', params={'smiles': '[H:4][H:9]', 'mapped': True})
    assert mapped.status_code == 200
    assert 'image/svg+xml' in mapped.headers['content-type']
    assert 'atom-0' in mapped.text and 'atom-1' in mapped.text


def test_single_event_example_annotation_export_and_restore(client):
    request = client.get('/api/single-event/example').json()
    assert request['rules'][0]['kind'] == 'elementary'
    assert client.post('/api/single-event', json={}).status_code == 409
    assert client.post('/api/run', json=request).status_code == 202
    assert finish(client)['status'] == 'finished'
    payload = client.get('/api/network').json()
    se = payload['network']['reactions'][0]['channels'][0]['single_event']
    assert se['status'] == 'calculated' and se['factor'] == 4
    assert se['sigma_reactants'] == [12, 1] and se['sigma_ts'] == 3
    assert 'Independent single-event factor' in client.get('/static/app.js').text
    class_id = se['parameter_class']
    options = {**request['single_event'], 'parameters': {class_id: {'Ea': 42000, 'source': 'Illustrative fit input'}}}
    updated = client.post('/api/single-event', json=options)
    assert updated.status_code == 200
    assert updated.json()['network']['single_event']['parameters'][class_id]['Ea'] == 42000
    before = client.get('/api/network').json()
    assert client.post('/api/single-event', json={'calculator': 'unknown'}).status_code == 422
    assert client.get('/api/network').json() == before
    with ZipFile(BytesIO(client.get('/api/export/csv').content)) as archive:
        assert {'single_event.csv', 'kinetic_parameters.csv'} <= set(archive.namelist())
        assert '42000' in archive.read('kinetic_parameters.csv').decode('utf-8-sig')
    project = client.get('/api/export/project').json()
    assert client.post('/api/new').status_code == 200
    assert client.post('/api/open', json=project).json() == before
    raw = client.post('/api/open', json=before['network']).json()
    assert raw['request']['single_event'] == raw['network']['single_event']['options']
    assert raw['network'] == before['network']


def test_generation_downloads_and_project_restore(client):
    request = {'seeds': [{'smiles': 'CCC', 'name': 'propane'}],
               'config': {'max_carbon': 6, 'max_depth': 2}}
    assert client.post('/api/run', json=request).status_code == 202
    state = finish(client)
    assert state['status'] == 'finished'
    assert state['summary']['termination'] == 'max_depth'
    assert state['summary']['species'] == 4
    payload = client.get('/api/network').json()
    network = Network.from_dict(payload['network'])
    assert network.summary() == state['summary']
    assert network.reactions[0].channels[0]['match_count'] == 24
    project = client.get('/api/export/project').json()
    assert project['request']['seeds'][0]['smiles'] == 'CCC'
    assert len(project['request']['rules']) == 4
    with ZipFile(BytesIO(client.get('/api/export/csv').content)) as archive:
        assert set(archive.namelist()) == {'species.csv', 'reactions.csv', 'channels.csv', 'stoichiometry.csv', 'network.dot'}
        assert 'depth' in archive.read('species.csv').decode('utf-8-sig')
    assert client.get('/api/export/json').json() == payload['network']
    assert client.post('/api/new').status_code == 200
    assert client.get('/api/network').status_code == 409
    restored = client.post('/api/open', json=project)
    assert restored.status_code == 200
    assert restored.json() == payload
    raw = client.post('/api/open', json=payload['network'])
    assert raw.status_code == 200
    assert raw.json()['request']['seeds'][0]['name'] == 'propane'
    assert raw.json()['network'] == payload['network']


@pytest.mark.parametrize('kind, expected', [('gas', [2, 6]), ('metal', [1]), ('acid', [1, 1])])
def test_path_examples_sites_export_and_regeneration(client, kind, expected):
    request = client.get('/api/single-event/path-example/' + kind).json()
    for seed in request['seeds']:
        assert client.post('/api/preview', json=seed).status_code == 200
    assert client.post('/api/run', json=request).status_code == 202
    assert finish(client)['status'] == 'finished'
    payload = client.get('/api/network').json()
    annotations = [c['single_event'] for r in payload['network']['reactions'] for c in r['channels']]
    assert sorted(a['factor'] for a in annotations) == expected
    assert all(a['method'] == 'reaction_path_degeneracy' for a in annotations)
    if kind != 'gas':
        assert any(s['site_atoms'] for s in payload['network']['species'])
    with ZipFile(BytesIO(client.get('/api/export/csv').content)) as archive:
        assert 'reaction_path_degeneracy' in archive.read('single_event.csv').decode('utf-8-sig')
        assert 'site_atoms' in archive.read('species.csv').decode('utf-8-sig')
    project = client.get('/api/export/project').json()
    assert client.post('/api/new').status_code == 200
    assert client.post('/api/open', json=project).json() == payload
    restored = client.post('/api/open', json=payload['network']).json()
    assert restored['network'] == payload['network']
    assert client.post('/api/run', json=restored['request']).status_code == 202
    assert finish(client)['status'] == 'finished'
    regenerated = client.get('/api/network').json()['network']
    assert regenerated['single_event'] == payload['network']['single_event']


def test_targets_and_input_errors(client):
    request = {'seeds': [{'smiles': 'CCC'}], 'targets': [{'smiles': 'C=C'}],
               'config': {'max_carbon': 6, 'max_depth': 4}}
    assert client.post('/api/run', json=request).status_code == 202
    summary = finish(client)['summary']
    assert summary['termination'] == 'targets_reached'
    assert summary['targets_found'] == 1
    assert summary['depth_reached'] == 4
    assert client.post('/api/run', json={'seeds': []}).status_code == 422
    assert client.post('/api/run', json={'seeds': [{'smiles': 'bad!'}]}).status_code == 422
    assert client.post('/api/run', json={'seeds': [{'smiles': 'CCC'}], 'rules': []}).status_code == 422
    assert client.post('/api/run', json={'seeds': [{'smiles': 'CCC'}], 'config': {'max_depth': -1}}).status_code == 422
    assert client.post('/api/open', json={'wrong': 'format'}).status_code == 422


def test_background_error_is_reported(client):
    assert client.post('/api/run', json={'seeds': [{'smiles': 'CCC'}], 'config': {'max_carbon': 1}}).status_code == 202
    state = finish(client)
    assert state['status'] == 'error'
    assert 'max_carbon' in state['error']
    assert client.get('/api/network').status_code == 409


def test_stop_keeps_a_consistent_network_and_prevents_overlapping_jobs(client, monkeypatch):
    entered, release = Event(), Event()
    original = Generator.generate

    def wait_for_stop(self, seeds, **kwargs):
        entered.set()
        assert release.wait(5)
        return original(self, seeds, **kwargs)

    monkeypatch.setattr(Generator, 'generate', wait_for_stop)
    request = {'seeds': [{'smiles': 'CCC'}]}
    assert client.post('/api/run', json=request).status_code == 202
    assert entered.wait(5)
    try:
        assert client.post('/api/run', json=request).status_code == 409
        assert client.post('/api/new').status_code == 409
        assert client.post('/api/stop').json()['status'] == 'stopping'
        assert client.get('/api/status').json()['status'] == 'stopping'
    finally:
        release.set()
    state = finish(client)
    assert state['summary']['termination'] == 'user_stopped'
    saved = client.get('/api/network').json()['network']
    assert len(saved['species']) == 1
    assert saved['reactions'] == []
    assert client.get('/api/export/project').status_code == 200


def test_pathways_continue_generation_export_and_restore(client):
    request = {'seeds': [{'smiles': 'CCC'}], 'targets': [{'smiles': 'C=C'}],
               'config': {'max_carbon': 6, 'max_depth': 4, 'stop_when_targets_found': False}}
    assert client.post('/api/run', json=request).status_code == 202
    assert finish(client)['summary']['termination'] == 'max_depth'
    payload = client.get('/api/network').json()
    target = payload['network']['targets'][0]['species_id']
    assert client.post('/api/pathways', json={'target': target}).status_code == 202
    assert finish(client)['analysis_termination'] == 'exhausted'
    result = client.get('/api/network').json()['analysis']
    assert len(result['pathways']) == 2
    assert all(len(p['reactions']) == 4 for p in result['pathways'])
    report = client.get('/api/export/pathways').json()
    assert report['analysis'] == result
    assert report['network'] == payload['network']
    project = client.get('/api/export/project').json()
    assert project['request']['config']['stop_when_targets_found'] is False
    assert client.post('/api/new').status_code == 200
    assert client.post('/api/open', json=project).json()['analysis'] == result
    assert client.post('/api/pathways', json={'target': target, 'max_pathways': 1}).status_code == 202
    assert finish(client)['analysis_termination'] == 'max_pathways'
    assert client.post('/api/open', json=payload['network']).json()['analysis'] is None
    assert client.post('/api/pathways', json={'target': 'not-present'}).status_code == 422
    assert client.post('/api/pathways', json={'target': target, 'max_reactions': -1}).status_code == 422


def test_pathway_analysis_can_be_stopped_without_discarding_network(client, monkeypatch):
    from chemweave import web
    request = {'seeds': [{'smiles': 'CCC'}], 'config': {'max_depth': 1}}
    assert client.post('/api/run', json=request).status_code == 202
    finish(client)
    before = client.get('/api/network').json()['network']
    target = before['species'][-1]['id']
    entered, release = Event(), Event()
    original = web.find_pathways

    def wait_for_stop(*args, **kwargs):
        entered.set()
        assert release.wait(5)
        return original(*args, **kwargs)

    monkeypatch.setattr(web, 'find_pathways', wait_for_stop)
    assert client.post('/api/pathways', json={'target': target}).status_code == 202
    assert entered.wait(5)
    try:
        assert client.get('/api/status').json()['status'] == 'analyzing'
        assert client.post('/api/pathways', json={'target': target}).status_code == 409
        assert client.post('/api/run', json=request).status_code == 409
        assert client.post('/api/open', json=before).status_code == 409
        assert client.post('/api/new').status_code == 409
        assert client.post('/api/stop').json()['status'] == 'stopping'
    finally:
        release.set()
    assert finish(client)['analysis_termination'] == 'user_stopped'
    after = client.get('/api/network').json()
    assert after['network'] == before
    assert after['analysis']['pathways'] == []


@pytest.mark.parametrize('case_id', ['document_CH_homolysis', 'document_hydride_shift'])
def test_user_templates_compute_without_family_or_recipe(client, case_id):
    from pathlib import Path
    cases = json.loads((Path(__file__).parent / 'data/path_cases.json').read_text())
    case = next(c for c in cases if c['id'] == case_id)
    template = {k: v for k, v in case['template'].items() if k not in {'family', 'kind'}}
    checked = client.post('/api/rules/validate', json=[template])
    assert checked.status_code == 200
    assert checked.json()['rules'][0]['kind'] == 'elementary'
    assert checked.json()['rules'][0]['family'] == 'unspecified'
    request = {'seeds': case['reactants'], 'rules': [template],
               'config': {'max_carbon': 10, 'max_depth': 1},
               'single_event': {'calculator': 'mapped_event_v1'}}
    assert client.post('/api/run', json=request).status_code == 202
    assert finish(client)['status'] == 'finished'
    payload = client.get('/api/network').json()
    se = payload['network']['reactions'][0]['channels'][0]['single_event']
    assert se['n_e'] == case['expected'][0]['factor']
    assert se['algorithm'] == 'mapped_event_orbit_v1'
    assert not {'parameter_class', 'rate_basis', 'A_units'} & se.keys()
    assert payload['network']['single_event']['parameters'] == {}
    with ZipFile(BytesIO(client.get('/api/export/csv').content)) as archive:
        evidence = archive.read('single_event.csv').decode('utf-8-sig')
        assert 'mapped_event_orbit_v1' in evidence and 'identical_reactant_divisor' in evidence
    project = client.get('/api/export/project').json()
    assert client.post('/api/new').status_code == 200
    assert client.post('/api/open', json=project).json() == payload
    raw = client.post('/api/open', json=payload['network']).json()
    assert client.post('/api/run', json=raw['request']).status_code == 202
    assert finish(client)['status'] == 'finished'
    second = client.get('/api/network').json()['network']['reactions'][0]['channels'][0]['single_event']
    assert second == se
    assert client.post('/api/rules/validate', json=[]).status_code == 422
    assert client.post('/api/rules/validate', json=[{'id': 'bad', 'smarts': 'bad'}]).status_code == 422
