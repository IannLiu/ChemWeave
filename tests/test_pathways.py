"""Hypergraph routes: complete branches, shared coproducts and no redundant loops."""

from itertools import combinations
import random

import pytest

from chemweave import Network, Reaction, Species, find_pathways


def graph(edges, seeds=('A',)):
    # Symbolic chemistry isolates the graph algorithm from SMARTS execution.
    ids = dict.fromkeys([*seeds, *(s for sides in edges for side in sides for s in side)])
    return Network(
        species={sid: Species.from_smiles('C', name=sid) for sid in ids},
        seeds=list(seeds),
        reactions=[Reaction(f'R{i}', tuple((s, left.count(s)) for s in dict.fromkeys(left)),
                            tuple((s, right.count(s)) for s in dict.fromkeys(right)),
                            'test', 'test', 'overall')
                   for i, (left, right) in enumerate(edges, 1)],
    )


def routes(result):
    return {frozenset(p['reactions']) for p in result['pathways']}


def test_alternatives_require_both_branches_and_keep_coproducts():
    net = graph([('A', 'B'), ('A', 'C'), ('BC', 'PQ'), ('A', 'D'), ('D', 'P'),
                 ('P', 'X'), ('X', 'P'), ('A', 'Z')])
    result = find_pathways(net, 'P')
    assert routes(result) == {frozenset(['R1','R2','R3']), frozenset(['R4','R5'])}
    assert result['termination'] == 'exhausted'
    assert result['network_search_complete']
    branch = next(p for p in result['pathways'] if 'R3' in p['reactions'])
    assert set(branch['species']) == set('ABCPQ')
    assert branch['reactions'][-1] == 'R3'
    assert routes(find_pathways(net, 'P', max_reactions=2)) == {frozenset(['R4','R5'])}
    assert not find_pathways(net, 'P', max_reactions=1)['pathways']


def test_unseeded_cycles_and_missing_coreactants_are_not_routes():
    net = graph([('A', 'B'), ('BC', 'P'), ('D', 'C'), ('C', 'D')])
    assert find_pathways(net, 'P')['pathways'] == []
    net = graph([('A', 'B'), ('BC', 'PQ'), ('Q', 'C')])
    assert find_pathways(net, 'P')['pathways'] == []


def test_shared_precursor_coproducts_and_self_coupling():
    net = graph([('A', 'BC'), ('B', 'D'), ('C', 'E'), ('DE', 'PQ')])
    result = find_pathways(net, 'P', max_reactions=4)
    assert routes(result) == {frozenset(['R1','R2','R3','R4'])}
    assert len(result['pathways']) == 1  # Independent branch orders are identical routes.
    assert result['pathways'][0]['reactions'][0] == 'R1'
    assert result['pathways'][0]['reactions'][-1] == 'R4'
    net = graph([('A', 'B'), ('BB', 'P')])
    assert len(find_pathways(net, 'P')['pathways']) == 1
    assert net.reactions[1].reactants == (('B', 2),)


def test_useful_coproducts_are_not_lost_by_blanket_species_cycle_removal():
    net = graph([('A', 'B'), ('B', 'CY'), ('C', 'BX'), ('XY', 'P')])
    result = find_pathways(net, 'P')
    assert result['pathways'][0]['reactions'] == ['R1', 'R2', 'R3', 'R4']
    # B reappears, but X and Y are necessary. The causal execution has no cycle.
    assert len(result['pathways']) == 1


def test_redundant_branch_is_removed_even_if_it_ran_before_the_target():
    net = graph([('A', 'B'), ('A', 'BC'), ('BC', 'P')])
    assert routes(find_pathways(net, 'P')) == {frozenset(['R2','R3'])}


def test_seed_target_limits_cancellation_and_partial_network():
    net = graph([('A', 'B'), ('B', 'P'), ('A', 'P')])
    result = find_pathways(net, 'A', max_reactions=0)
    assert result['pathways'] == [{'reactions': [], 'species': ['A']}]
    net.termination = 'targets_reached'
    result = find_pathways(net, 'P', max_pathways=1)
    assert not result['network_search_complete']
    assert result['termination'] == 'max_pathways'
    assert len(result['pathways']) == 1
    result = find_pathways(net, 'P', stop_requested=lambda: True)
    assert result['termination'] == 'user_stopped'
    assert result['pathways'] == []
    with pytest.raises(ValueError, match='Unknown target'):
        find_pathways(net, 'missing')


def test_small_networks_agree_with_exhaustive_subset_oracle():
    rng = random.Random(19)
    for _ in range(30):
        edges = [( ''.join(rng.sample('ABCDP', rng.randint(1, 2))),
                   ''.join(rng.sample('ABCDP', rng.randint(1, 2)))) for _ in range(7)]
        net = graph(edges, seeds=('A', 'B'))
        if 'P' not in net.species:
            continue
        successful = set()
        for size in range(5):
            for subset in combinations(net.reactions, size):
                ids = frozenset(r.id for r in subset)
                if any(route <= ids for route in successful):
                    continue
                available = set(net.seeds)
                while True:
                    previous = set(available)
                    for reaction in subset:
                        if all(s in available for s, _ in reaction.reactants):
                            available.update(s for s, _ in reaction.products)
                    if available == previous:
                        break
                if 'P' in available:
                    successful.add(ids)
        assert routes(find_pathways(net, 'P', max_reactions=4, max_pathways=1000)) == successful


def test_cancellation_retains_completed_routes():
    net = graph([('A', 'P'), ('A', 'B'), ('B', 'P')])
    checks = 0

    def stop_after_first_route():
        nonlocal checks
        checks += 1
        return checks > 2

    result = find_pathways(net, 'P', stop_requested=stop_after_first_route)
    assert result['termination'] == 'user_stopped'
    assert routes(result) == {frozenset(['R1'])}
