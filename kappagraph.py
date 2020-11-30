# Walter Fontana, 2020

import networkx as nx


class KappaGraph:
    def __init__(self, komplex):
        """
        Kappa complex to networkx graph (plus multi-graph flag)
        :param komplex:
        :return:
        """
        self.komplex = komplex
        self.multigraph = False
        if not komplex.bonds:  # a singleton node
            self.nxGraph = nx.Graph()
            name = next(iter(komplex.agents))
            self.nxGraph.add_node(name)
        else:
            if komplex.is_multigraph():
                self.nxGraph = nx.MultiGraph()
                self.multigraph = True
            else:
                self.nxGraph = nx.Graph()
            for (a1, s1), (a2, s2) in komplex.bonds:
                txt = a1 + '@' + s1 + '--' + a2 + '@' + s2
                self.nxGraph.add_edge(a1, a2, bond=txt)

        # set node attributes
        nattr = {}
        for node, nodedata in self.nxGraph.nodes.items():
            nattr[node] = {'type': node.split(komplex.idsep[0])[0], 'id': komplex.extract_identifier(node)[1]}
        nx.set_node_attributes(self.nxGraph, nattr)

    def get_cycle(self):
        """
        A wrapper for networkx find_cycle()
        :param nxG: networkX graph
        :return: list of tuples, such as [(0, 1), (1, 2), (0, 2)]
                 or []
        """
        try:
            cycle = nx.find_cycle(self.nxGraph, orientation='ignore')
            edge_list = []
            if cycle:
                if len(cycle[0]) == 3:
                    edge_list = [(tail, head) for tail, head, discard in cycle]
                elif len(cycle[0]) == 4:
                    edge_list = [(tail, head) for tail, head, discard, discard in cycle]
            return edge_list
        except nx.NetworkXNoCycle:
            return []

    def get_cycle_basis(self, minimum=True):
        """
        A wrapper for networkx cycle basis finder
        :param minimum: True (default), if minimum cycle basis is sought
        :param nxG: networkX graph
        :return: list of edge lists, such as [[(0, 1), (1, 2), (2, 3)]...]
                 or []
                 and the number of self-loops discarded in the conversion from multi-graph to graph
        """
        if self.multigraph:
            G, N = self.convert_multigraph_to_graph()
        else:
            G = self.nxGraph
            N = 0
        # A list of cycle lists. Each cycle list is a list of nodes which forms a cycle (loop) in G.
        # The nodes are not necessarily returned in a order by which they appear in the cycle...
        if minimum:
            basis = nx.minimum_cycle_basis(G)
        else:
            basis = nx.cycle_basis(G)

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

    def convert_multigraph_to_graph(self):
        G = nx.Graph()
        n = 0
        for u, v, discard in self.nxGraph.edges:
            if not G.has_edge(u, v):
                G.add_edge(u, v)
            else:
                n += 1
        return G, n

    def delete_edgelists(self, edge_list=[]):
        # to unify handling, convert to a list of lists (such as coming from a cycle basis)
        if not isinstance(edge_list[0], list):
            edge_list = [edge_list]

        for list_of_edges in edge_list:
            for e in list_of_edges:
                self.nxGraph.remove_edge(e[0], e[1])

    def delete_nodelist(self, node_list=[]):
        for node in node_list:
            self.nxGraph.remove_node(node)

    def write_dot(self, filename='complex.dot', uniform=True, shape='oval'):
        palette = ('blue', 'red', 'green', 'cyan', 'magenta', 'yellow', 'khaki', 'silver')
        shapelette = ('circle', 'triangle', 'polygon', 'oval', 'diamond', 'house', 'hexagon',
                      'parallelogram', 'pentagon', 'rectangle')
        # assign colors to nodes
        color = {}
        shapes = {}
        i = 0
        # fill palette index in order of (descending) frequency
        for type in self.komplex.composition.keys():
            color[type] = i % len(palette)
            shapes[type] = i % len(shapelette)
            i += 1
        for node in self.nxGraph.nodes():
            self.nxGraph.nodes[node]['style'] = 'filled'
            self.nxGraph.nodes[node]['fillcolor'] = palette[color[self.nxGraph.nodes[node]['type']]]
            if not uniform:
                self.nxGraph.nodes[node]['shape'] = shapelette[shapes[self.nxGraph.nodes[node]['type']]]
            else:
                self.nxGraph.nodes[node]['shape'] = shape
        nx.nx_agraph.write_dot(self.nxGraph, filename)


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
    g = KappaGraph(c)
    cycle = g.get_cycle()
    basis, n = g.get_cycle_basis()
    print(cycle)
    print(basis)
    r = viz.Renderer(c)
    # r.html_render()
    r.render()
    r.color_edgelists(edge_list=cycle)
    r.show()
