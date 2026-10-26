"""Seed-supported formation routes in a directed reaction hypergraph."""

from collections import defaultdict


def find_pathways(network, target, *, max_reactions=8, max_pathways=100,
                  stop_requested=None):
    """Enumerate minimal reaction sets that form a target species ID.

    Every reactant must be available before a reaction runs. Availability is
    qualitative: quantities are not consumed, and each reaction is used at
    most once. Coproducts can support other branches. A route is minimal when
    deleting any one of its reactions makes the target unreachable within
    that route. Alternative routes need not have the same length.

    Results contain reaction IDs in one executable order, not every ordering
    of independent branches or every combination of atom-mapping channels.
    Enumeration follows stored reaction order; a capped result is not a
    ranking of the shortest or most likely routes.
    """
    if target not in network.species:
        raise ValueError(f"Unknown target species ID: {target}")
    if max_reactions < 0 or max_pathways < 1:
        raise ValueError("max_reactions must be nonnegative and max_pathways positive")

    seeds = set(network.seeds)
    reactions = {r.id: r for r in network.reactions}
    inputs = {r.id: {sid for sid, _ in r.reactants} for r in network.reactions}
    outputs = {r.id: {sid for sid, _ in r.products} for r in network.reactions}
    producers = defaultdict(list)
    for rid in reactions:
        for sid in outputs[rid]:
            producers[sid].append(rid)

    # Reverse traversal only narrows the search; it does not certify a route.
    candidates, visited_species = set(), set()
    pending = [target]
    while pending:
        sid = pending.pop()
        if sid in seeds or sid in visited_species:
            continue
        visited_species.add(sid)
        for rid in producers[sid]:
            candidates.add(rid)
            pending.extend(inputs[rid])
    candidates = [rid for rid in reactions if rid in candidates]

    def reaches_target(route):
        """Test a reduced route without assuming its original execution order."""
        available = set(seeds)
        pending = list(route)
        while pending and target not in available:
            enabled = [rid for rid in pending if inputs[rid] <= available]
            if not enabled:
                break
            for rid in enabled:
                available.update(outputs[rid])
                pending.remove(rid)
        return target in available

    result = {
        "target": target,
        "max_reactions": max_reactions,
        "max_pathways": max_pathways,
        "network_termination": network.termination,
        "network_search_complete": network.termination in {"exhausted", "max_depth"},
        "termination": "exhausted",
        "pathways": [],
        "states_examined": 0,
    }
    # Explicit DFS avoids recursion limits. Sets identify routes, so different
    # execution orders of independent branches cannot produce duplicates.
    stack = [((), seeds)]
    seen = {frozenset()}
    while stack:
        if stop_requested is not None and stop_requested():
            result["termination"] = "user_stopped"
            break
        route, available = stack.pop()
        result["states_examined"] += 1
        if target in available:
            if any(reaches_target(route[:i] + route[i + 1:]) for i in range(len(route))):
                continue
            species = seeds & {target}
            for rid in route:
                species.update(inputs[rid] | outputs[rid])
            result["pathways"].append({
                "reactions": list(route),
                "species": [sid for sid in network.species if sid in species],
            })
            if len(result["pathways"]) >= max_pathways:
                result["termination"] = "max_pathways"
                break
            continue
        if len(route) >= max_reactions:
            continue
        for rid in reversed(candidates):
            if not inputs[rid] <= available or outputs[rid] <= available:
                continue
            key = frozenset((*route, rid))
            if key in seen:
                continue
            seen.add(key)
            stack.append(((*route, rid), available | outputs[rid]))
    return result
