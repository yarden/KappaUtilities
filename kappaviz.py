import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import plotly.graph_objects as go


def get_nxgraph(komplex):
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

    # set edge attributes
    attr = {}
    for ((a1, s1), (a2, s2)) in komplex.bonds:
        txt = komplex.info[a1]['type'] + '.' + s1 + '<-->' + komplex.info[a2]['type'] + '.' + s2
        attr[(a1, a2)] = {'sites': txt}
    nx.set_edge_attributes(nxG, attr)

    return nxG


def write_dot(komplex, file_name='complex.dot',
              palette=('blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'khaki', 'silver'),
              shapalette=('polygon', 'oval', 'circle', 'triangle', 'diamond', 'house', 'hexagon',
                          'parallelogram', 'pentagon', 'rectangle'),
              uniform=True, shape='oval'):
    nxG = get_nxgraph(komplex)

    # assign colors to nodes
    color = {}
    shapes = {}
    i = 0
    # fill palette index in order of (descending) frequency
    for type in komplex.composition.keys():
        color[type] = i % len(palette)
        shapes[type] = i % len(shapalette)
        i += 1
    for node in nxG.nodes():
        nxG.nodes[node]['style'] = 'filled'
        nxG.nodes[node]['fillcolor'] = palette[color[nxG.nodes[node]['type']]]
        if not uniform:
            nxG.nodes[node]['shape'] = shapalette[shapes[nxG.nodes[node]['type']]]
        else:
            nxG.nodes[node]['shape'] = shape

    nx.nx_agraph.write_dot(nxG, file_name)


class Renderer:
    def __init__(self, komplex, prog='neato', node_info=True):
        """
        In establishing a Renderer object, a layout of nodes is triggered.
        Subsequently, various display methods can be invoked.

        :param komplex:
        :param prog: any of neato, dot, twopi, circo, fdp, nop, wc, acyclic, gvpr, gvcolor, ccomps,
                            sccmap, tred, sfdp, unflatten
        :param node_info:
        """
        # (komplex is never assigned, so it should be effectively 'pass by reference')

        # create the networkx version of 'komplex'
        self.nxG = get_nxgraph(komplex)
        self.composition = komplex.composition

        # positions = nx.nx_pydot.graphviz_layout(self.nxG, prog=prog)
        self.positions = nx.nx_agraph.graphviz_layout(self.nxG, prog=prog)

        # self.node_hover_text = []
        # if node_info:
        #     for node in self.nxG.nodes:
        #         iface = komplex.agents[node]
        #         iface = dict(sorted(iface.items()))
        #         info = f"<b>{node}</b><br>"
        #         for site in iface.keys():
        #             info += f"{site:>10}  ->  state: {iface[site]['state']:<5} bond: {iface[site]['bond']}<br>"
        #         self.node_hover_text += [info[:-4]]

        self.node_hover_text = {}
        if node_info:
            for type in self.composition.keys():
                self.node_hover_text[type] = []
            for node in self.nxG.nodes:
                iface = komplex.agents[node]
                iface = dict(sorted(iface.items()))
                info = f"<b>{node}</b><br>"
                for site in iface.keys():
                    info += f"{site:>10}  ->  state: {iface[site]['state']:<5} bond: {iface[site]['bond']}<br>"
                self.node_hover_text[self.nxG.nodes[node]['type']] += [info[:-4]]

    def nx_render(self,
                  file_name='',
                  labels='short',
                  node_size=400,
                  font_size=11,
                  palette=('b', 'g', 'r', 'c', 'm', 'y', 'k', 'w')):
        """
        Render a networkx graph with matplotlib.

        :param file_name:
        :param labels:
        :param node_size:
        :param font_size:
        :param palette:
        :return:
        """
        nx_options = {'font_size': font_size,
                      'font_color': 'white',
                      'font_weight': 'bold',
                      'node_size': node_size}
        # set labels
        if labels == 'type':
            nx_options['labels'] = {node: self.nxG.nodes[node]['type'] for node in self.nxG.nodes}
        elif labels == 'short':
            nx_options['labels'] = {node: self.nxG.nodes[node]['id'] for node in self.nxG.nodes}
        elif labels == 'full':
            nx_options['labels'] = {node: self.nxG.nodes[node]['type'] +
                                          self.nxG.nodes[node]['id'] for node in self.nxG.nodes}
        else:
            nx_options['with_labels'] = False

        # assign colors to nodes
        clr = {}
        i = 0
        # fill palette index in order of (descending) frequency
        for type in self.composition.keys():
            clr[type] = i % len(palette)
            i += 1
        nx_options['node_color'] = []
        for node in self.nxG.nodes:
            nx_options['node_color'] += [palette[clr[self.nxG.nodes[node]['type']]]]

        nx.draw_networkx(self.nxG, pos=self.positions, **nx_options)

        # edge_labels = nx.get_edge_attributes(self.nxG, 'sites')
        # nx.draw_networkx_edge_labels(self.nxG, self.positions, edge_labels=edge_labels, font_size=4)

        colors = [f'{palette[clr[n]]}' for n in clr.keys()]
        items = [Line2D([0, 1], [0, 1],
                        color='white',
                        marker='o',
                        markersize=7,
                        markerfacecolor=clr,
                        linewidth=0) for clr in colors]
        labels = [f'{node}' for node in clr.keys()]
        plt.legend(items, labels)

        if file_name:
            plt.savefig(file_name)
        else:
            plt.show()

    def html_render(self, file_name='', node_size=15,
                    palette=('blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'khaki', 'silver')):
        """
        Render a networkx graph with plotly.

        :param file_name:
        :param node_size:
        :param palette:
        :return:
        """
        edge_x = []
        edge_y = []
        for edge in self.nxG.edges():
            x0, y0 = self.positions[edge[0]]
            x1, y1 = self.positions[edge[1]]
            edge_x.append(x0)
            edge_x.append(x1)
            edge_x.append(None)
            edge_y.append(y0)
            edge_y.append(y1)
            edge_y.append(None)

        # edges as lines in a scatter trace
        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=1, color='gray'),
            hoverinfo='none',
            showlegend=False,
            mode='lines')

        # stratify by type (so we can have a legend for each type)
        clr = {}
        i = 0
        node_x = {}
        node_y = {}
        for type in self.composition.keys():
            node_x[type] = []
            node_y[type] = []
            # fill palette index in order of (descending) frequency
            clr[type] = i % len(palette)
            i += 1

        for node in self.nxG.nodes():
            x, y = self.positions[node]
            node_x[self.nxG.nodes[node]['type']].append(x)
            node_y[self.nxG.nodes[node]['type']].append(y)

        # nodes as a scatter traces
        node_trace = {}
        for type in self.composition.keys():
            node_trace[type] = go.Scatter(
                x=node_x[type],
                y=node_y[type],
                name=type,
                legendgroup=type,
                showlegend=True,
                mode='markers',
                text=self.node_hover_text[type],
                hoverinfo='text',
                marker=dict(
                    color=palette[clr[type]],
                    size=node_size,
                    line_width=2
                )
            )

        # # info for edges
        # # can't do that in plotly right now.
        # edge_trace.text = []
        # for n1, n2, data in self.nxG.edges.data():
        #     print(data['sites'])

        layout = go.Layout(
            plot_bgcolor='white',
            hovermode='closest',
            margin=dict(b=20, l=5, r=5, t=40),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            hoverlabel=dict(
                # bgcolor="white",
                font_size=11,
                font_family="Consolas"
            )
        )

        fig = go.Figure(layout=layout)

        fig.add_trace(edge_trace)
        for type in self.composition.keys():
            fig.add_trace(node_trace[type])

        config = dict({'scrollZoom': True,
                       'displayModeBar': True,
                       'displaylogo': False,
                       'modeBarButtonsToAdd': ['drawline', 'drawopenpath', 'drawclosedpath',
                                               'drawcircle', 'drawrect', 'eraseshape']
                       }
                      )

        if file_name:
            if file_name.endswith('.html'):
                fig.write_html(file_name)
            else:
                fig.write_image(file_name)
        else:
            fig.show(config=config)


def show_ranked_complexes(snapshot, sort='size', cutoff=3, cols=3, rows=1, prog='neato', node_size=20, font_size=4):
    """
    Display the ranked complexes of a snapshot.

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
        r = Renderer(c, prog=prog)
        r.nx_render()
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
    # write_dot(c1, 'complex.dot')
    r = Renderer(c1)
    # r.nx_render(node_size=100, labels='none')
    # r.html_render()
    r.html_render(file_name='complex.html')
    # r.nx_render(c1, prog='sfdp', node_size=80, font_size=4)
    # r.nx_render(c1, prog='fdp', node_size=80, font_size=4)
    # r.nx_render(c1, prog='neato', node_size=80, font_size=4)
    # r.nx_render(c1, prog='twopi', node_size=80, font_size=4)
    # r.nx_render(c1, prog='dot', node_size=80, font_size=4)
    # r.nx_render(c1, prog='circo', node_size=80, font_size=4)
