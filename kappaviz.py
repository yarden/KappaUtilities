import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


def render(komplex, prog='neato', node_size=400, font_size=11, labels=True, pdf=None):
    """
    render a networkx graph

    prog can be any of
        neato, dot, twopi, circo, fdp, nop, wc, acyclic, gvpr,
        gvcolor, ccomps, sccmap, tred, sfdp, unflatten

    colors can be one of {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'}
    """

    if komplex.nxgraph is None:
        komplex.get_nxgraph_from_structure()
    nxgraph = komplex.nxgraph

    options = {'font_size': font_size,
               'font_color': 'white',
               'font_weight': 'bold'}

    palette = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

    if labels:
        # extract node labels for graphic display (these would be typically shorter than the node names)
        label = {}
        for node in nxgraph.nodes:
            label[node] = nxgraph.nodes[node]['id']
        options['labels'] = label
    else:
        options['with_labels'] = False

    color = []
    clr = {}
    i = 0
    # fill palette in order of (descending) frequency
    for node in komplex.composition.keys():
        clr[node] = palette[i % 8]
        i += 1
    # fill color list
    for node in nxgraph.nodes:
        t = node.split(komplex.idsep[0])[0]
        color += [clr[t]]

    options['node_color'] = color
    options['node_size'] = node_size

    # positions = nx.nx_pydot.graphviz_layout(nxgraph, prog=prog)
    positions = nx.nx_agraph.graphviz_layout(nxgraph, prog=prog)
    nx.draw_networkx(nxgraph, pos=positions, **options)

    # edge_labels = nx.get_edge_attributes(nxgraph, 'sites')
    # nx.draw_networkx_edge_labels(nxgraph, positions, edge_labels=edge_labels, font_size=4)

    colors = [f'{clr[n]}' for n in clr.keys()]
    items = [Line2D([0, 1], [0, 1],
                    color='white',
                    marker='o',
                    markersize=7,
                    markerfacecolor=clr,
                    linewidth=0) for clr in colors]
    labels = [f'{node}' for node in clr.keys()]
    plt.legend(items, labels)

    if pdf:
        plt.savefig(pdf)
    else:
        plt.show()


def show_ranked_complexes(snapshot, sort='size', cutoff=3, cols=3, rows=1, prog='neato', node_size=20, font_size=4):
    """
    display the ranked complexes of a snapshot

    :param snapshot: a snapshot object
    :param sort: 'size' (default) or 'count'
    :param cols: # of columns of plots
    :param rows: # of rows of plots
    :param prog: layout program
    :param cutoff: size or count cutoff
    :param node_size: size of disk representing a node
    :param font_size: size of labels
    """
    ranking = []
    if sort == 'size':
        ranking = sorted(snapshot.complex, key=lambda x: x.size, reverse=True)
    elif sort == 'count':
        ranking = sorted(snapshot.complex, key=lambda x: x.count, reverse=True)

    plt.figure(frameon=False, figsize=(5., 5.))
    i = 1
    for c in ranking[0:cutoff]:
        plt.subplot(rows, cols, i)
        plt.title(sort + ' ' + str(c.size))
        render(c, prog=prog, node_size=node_size, font_size=font_size, labels=False)
        i += 1
    plt.show()


if __name__ == '__main__':
    import kappathings as kt
    import re

    # usage scenarios

    line = open('TestData/bigly.ka', 'r').read()
    # remove newlines that might occur in the file
    line = re.sub(r'\n+', ' ', line)
    # create a KappaComplex with whatever assignment of node identifiers arises
    # (that's the normalized=False flag).
    c1 = kt.KappaComplex(line, normalize=True)
    print(c1)
    print(f'is multi-graph: {c1.is_multigraph()}')
    print("visualize...")
    # these seem to work; other layout programs fail (a networkx issue)
    render(c1, prog='sfdp', node_size=80, font_size=4)
    # render(c1, prog='fdp', node_size=80, font_size=4)
    # render(c1, prog='neato', node_size=80, font_size=4)
    # render(c1, prog='twopi', node_size=80, font_size=4)
    # render(c1, prog='dot', node_size=80, font_size=4)
    # render(c1, prog='circo', node_size=80, font_size=4)
