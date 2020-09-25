# Walter Fontana at 15.09.2020

import networkx as nx


def make_nxgraph(komplex):
    """
    Generate a networkx graph. For now this assumes one connected component (such as a molecular species).
    (Note: update this to handle patterns.)

    :param komplex:
    :return:
    """
    if not komplex.bonds:  # we are dealing with a singleton node
        nxG = nx.Graph()
        name = next(iter(komplex.agents))
        nxG.add_node(name)
    else:
        if komplex.is_multigraph():
            nxG = nx.MultiGraph()
        else:
            nxG = nx.Graph()
        for (a1, s1), (a2, s2) in komplex.bonds:
            nxG.add_edge(a1, a2)

    # set node attributes
    attr = {}
    for node, nodedata in nxG.nodes.items():
        attr[node] = {'type': node.split(komplex.idsep[0])[0], 'id': komplex.extract_identifier(node)[1]}
    nx.set_node_attributes(nxG, attr)

    # if komplex.bonds:   # there's a problem with multigraphs of size 2, because (a1, a2) is sorted...
    #     # set edge attributes
    #     attr = {}
    #     for ((a1, s1), (a2, s2)) in komplex.bonds:
    #         txt = komplex.info[a1]['type'] + '.' + s1 + '<-->' + komplex.info[a2]['type'] + '.' + s2
    #         attr[(a1, a2)] = {'sites': txt}
    #     nx.set_edge_attributes(nxG, attr)

    return nxG


def get_cycle(nxG):
    """
    A wrapper for networkx find_cycle()
    :param nxG: networkX graph
    :return: list of tuples, such as [(0, 1), (1, 2), (0, 2)]
             or []
    """
    try:
        cycle = nx.find_cycle(nxG, orientation='ignore')
        edge_list = []
        if cycle:
            if len(cycle[0]) == 3:
                edge_list = [(tail, head) for tail, head, discard in cycle]
            elif len(cycle[0]) == 4:
                edge_list = [(tail, head) for tail, head, discard, discard in cycle]
        return edge_list
    except nx.NetworkXNoCycle:
        return []


def get_minimum_cycle_basis(nxG):
    """
    A wrapper for networkx cycle basis finder
    :param nxG: networkX graph
    :return: list of edge lists, such as [[(0, 1), (1, 2), (2, 3)]...]
             or []
             and the number of self-loops discarded in the conversion from multi-graph to graph
    """
    if nxG.is_multigraph():
        G, N = convert_multigraph_to_graph(nxG)
    else:
        G = nxG
        N = 0
    # A list of cycle lists. Each cycle list is a list of nodes which forms a cycle (loop) in G.
    # The nodes are not necessarily returned in a order by which they appear in the cycle...
    basis = nx.minimum_cycle_basis(G)

    # convert into edge lists...
    edge_lists = []
    for node_list in basis:
        edges = []
        start = node_list[0]
        p = start
        node_list.remove(p)
        while node_list:
            p_neighbors = {nbor for a, nbor in G.edges(p)}
            for q in p_neighbors:
                if q in node_list:
                    # edge between p and q
                    edges.append((p, q))
                    node_list.remove(q)
                    p = q
                    break
        edges.append((p, start))  # return to origin
        edge_lists.append(edges)
    return edge_lists, N


def convert_multigraph_to_graph(nxG):
    if not nxG.is_multigraph():
        return nxG, 0
    G = nx.Graph()
    n = 0
    for u, v, discard in nxG.edges:
        if not G.has_edge(u, v):
            G.add_edge(u, v)
        else:
            n += 1
    return G, n


def delete_edgelists(nxG, edge_list=[]):
    # to unify handling, convert to a list of lists (such as coming from a cycle basis)
    if not isinstance(edge_list[0], list):
        edge_list = [edge_list]

    newG = nxG.copy()
    for list_of_edges in edge_list:
        for e in list_of_edges:
            newG.remove_edge(e[0], e[1])
    return newG


def delete_nodelist(nxG, node_list=[]):
    newG = nxG.copy()
    for node in node_list:
        newG.remove_node(node)
    return newG


if __name__ == '__main__':
    import kappathings as kt
    import kappaviz as viz
    import kappasnap as ks

    # usage scenarios

    # kapparing2 = 'A(r[7] l[1]),A(r[1] l[2]),A(r[2] l[3]),A(r[3] l[4]),A(r[4] l[5]),A(r[5] l[6]),A(r[6] l[7])'
    # kappanoring = 'A(r[.] l[1]),A(r[1] l[2]),A(r[2] l[3]),A(r[3] l[4]),A(r[4] l[5]),A(r[5] l[6]),A(r[6] l[.])'
    kapparing = 'A(r[.] l[1]),A(r[1] l[2] m[7]),A(r[2] l[3]),A(r[3] l[4]),A(r[4] l[5] m[7]),A(r[5] l[6]),A(r[6] l[.])'
    # (that's the normalized=False flag).
    c = kt.KappaComplex(kapparing)
    c.show()
    g = make_nxgraph(c)
    cycle = get_cycle(g)
    basis, n = get_minimum_cycle_basis(g)
    print(cycle)
    print(basis)
    r = viz.Renderer(c)
    # r.html_render()
    r.render()
    r.color_edgelists(edge_list=cycle)
    viz.show()