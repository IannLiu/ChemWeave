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
