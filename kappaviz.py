# Walter Fontana, 2020

import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.collections as artcoll
import plotly.graph_objects as go
import kappagraph as kg
import traceback

import sys
from PyQt5.QtWidgets import QDialog, QApplication, QVBoxLayout, QLineEdit, QTextBrowser

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar


def show():
    plt.show(block=False)


class Window(QDialog):
    def __init__(self, fig, parent=None, rows=1, cols=1, size=(10, 10)):
        super(Window, self).__init__(parent)

        # Canvas Widget that displays the figure
        self.canvas = FigureCanvas(fig)

        # Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)

        # self.lineEdit = QLineEdit()
        # self.textBrowser = QTextBrowser()

        # set the layout
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        # layout.addWidget(self.textBrowser)
        # layout.addWidget(self.lineEdit)
        self.setLayout(layout)
        # self.lineEdit.returnPressed.connect(self.process)

    # def process(self):
    #     userInput = self.lineEdit.text()
    #     if userInput == "Comment vas tu?":
    #         self.textBrowser.append("Très bien.")
    #     else:
    #         self.textBrowser.append("What?")
    #     self.lineEdit.clear()


class Canvas:
    def __init__(self, rows=1, cols=1, size=(10, 10), qt5=False, name=None):
        if name is None:
            (filename, line_number, function_name, text) = traceback.extract_stack()[-2]
            name = text[:text.find('=')].strip()
        self.name = name
        self.figure = None
        self.qt = None
        self.rows = rows
        self.cols = cols

        # axs is a numpy ndarray
        self.figure, self.axes = plt.subplots(nrows=rows, ncols=cols, figsize=size, squeeze=False, frameon=False)
        for i in range(0, self.rows):
            for j in range(0, self.cols):
                self.axes[i, j].axis("off")
                self.axes[i, j].spines["top"].set_visible(False)
                self.axes[i, j].spines["right"].set_visible(False)
                self.axes[i, j].spines["left"].set_visible(False)
                self.axes[i, j].spines["bottom"].set_visible(False)

        if qt5:
            self.qt = Window(self.figure, rows=rows, cols=cols, size=size)

    def number_of_plot_areas(self):
        return self.rows * self.cols

    def __del__(self):
        if self.figure is not None:
            plt.close(self.figure)
        # print(f'{self.name} deleted')

    def clear(self, i, j):  # assuming i and j start with 1...
        self.axes[i-1, j-1].clear()

    def panel2index(self, panel):  # starting with 1,1
        return (panel[0] - 1) * self.cols + panel[1]

    def index2panel(self, n):  # starting with 1
        i = -(- n // self.cols)
        j = n - self.cols * (i - 1)
        return i, j


class Renderer:
    def __init__(self, komplex, prog='neato', node_info=False):
        """
        In establishing a Renderer object, a layout of nodes is triggered.
        Subsequently, various display methods can be invoked.

        :param komplex:
        :param prog: any of neato, dot, twopi, circo, fdp, nop, wc, acyclic, gvpr, gvcolor, ccomps,
                            sccmap, tred, sfdp, unflatten
        :param node_info:
        """

        self.ax = None
        self.figure = None

        self.Graph = kg.KappaGraph(komplex)
        self.nxGraph = self.Graph.nxGraph

        self.node_hover_text = {}
        if node_info:
            for type in komplex.composition.keys():
                self.node_hover_text[type] = []
            for node in self.nxGraph.nodes:
                iface = komplex.agents[node]
                iface = dict(sorted(iface.items()))
                info = f"<b>{node}</b><br>"
                for site in iface.keys():
                    info += f"{site:>10}  ->  state: {iface[site]['state']:<5} bond: {iface[site]['bond']}<br>"
                self.node_hover_text[self.nxGraph.nodes[node]['type']] += [info[:-4]]

        self.nx_palette = ('c', 'r', 'b', 'g', 'm', 'y', 'k', 'w')
        self.html_palette = ('blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'khaki', 'silver')

        self.nx_options = {'font_size': 11,
                           'font_color': 'white',
                           'font_weight': 'bold',
                           'node_size': 400,
                           'labels': {},
                           'edge_color': 'black',
                           'width': 1
                           }

        # assign colors to node types
        self.type_color = {}
        i = 0
        # fill palette index in order of (descending) frequency
        for typ in komplex.composition.keys():
            self.type_color[typ] = i % len(self.nx_palette)
            i += 1

        self.nx_options['node_color'] = []
        for node in self.nxGraph.nodes:
            self.nx_options['node_color'] += [self.nx_palette[self.type_color[self.nxGraph.nodes[node]['type']]]]
        self.legend_colors = [f'{self.nx_palette[self.type_color[n]]}' for n in self.type_color.keys()]

        # layout
        self.positions = nx.nx_agraph.graphviz_layout(self.nxGraph, prog=prog)

    def __del__(self):
        if self.figure is not None:
            plt.close(self.figure)

    def layout(self):
        self.positions = nx.nx_agraph.graphviz_layout(self.nxGraph, prog='neato')
        self.nx_options['node_color'] = []
        for node in self.nxGraph.nodes:
            self.nx_options['node_color'] += [self.nx_palette[self.type_color[self.nxGraph.nodes[node]['type']]]]
        self.legend_colors = [f'{self.nx_palette[self.type_color[n]]}' for n in self.type_color.keys()]

    def set_palette(self, palette):
        self.nx_palette = palette

    def set_html_palette(self, palette):
        self.html_palette = palette

    def render(self,
               canvas=None,
               panel=(1, 1),
               labels='short',
               node_size=40,
               font_size=9,
               line_width=1,
               edge_color='gray',
               legend=False,
               title="",
               title_font_size=10,
               figure_size=(6, 6)):
        """
        Render a networkx graph with matplotlib.

        :param figure_size:
        :param canvas:
        :param title_font_size:
        :param title:
        :param panel:
        :param legend:
        :param edge_color:
        :param line_width:
        :param labels:
        :param node_size:
        :param font_size:
        :return:
        """
        self.nx_options['font_size'] = font_size
        self.nx_options['node_size'] = node_size
        # set labels
        self.nx_options['with_labels'] = True
        if labels == 'type':
            self.nx_options['labels'] = {node: self.nxGraph.nodes[node]['type'] for node in self.nxGraph.nodes}
        elif labels == 'short':
            self.nx_options['labels'] = {node: self.nxGraph.nodes[node]['id'] for node in self.nxGraph.nodes}
        elif labels == 'full':
            self.nx_options['labels'] = {node: self.nxGraph.nodes[node]['type'] + self.nxGraph.nodes[node]['id']
                                         for node in self.nxGraph.nodes}
        else:
            self.nx_options['with_labels'] = False

        self.nx_options['edge_color'] = edge_color
        self.nx_options['width'] = line_width  # edge width

        if canvas is not None:
            self.ax = canvas.axes[panel[0] - 1, panel[1] - 1]
            self.figure = canvas.figure
            # we clear the panel since we are drawing the whole network
            self.ax.clear()
        else:
            if self.figure is not None:
                plt.close(self.figure)
            self.figure, self.ax = plt.subplots(figsize=figure_size, frameon=False)
            self.ax.axis("off")
            self.ax.spines["top"].set_visible(False)
            self.ax.spines["right"].set_visible(False)
            self.ax.spines["left"].set_visible(False)
            self.ax.spines["bottom"].set_visible(False)

        nx.draw_networkx(self.nxGraph, pos=self.positions, ax=self.ax, **self.nx_options)

        # the legend
        if legend:
            items = [Line2D([0, 1], [0, 1], color='white', marker='o', markersize=7, markerfacecolor=clr, linewidth=0)
                     for clr in self.legend_colors]
            labels = [f'{node}' for node in self.type_color.keys()]
            self.ax.legend(items, labels)

        if title:
            self.ax.set_title(title, fontsize=title_font_size)

    def color_edgelists(self, edge_list=[], line_width=1, edge_color='r'):
        # to unify handling, convert to a list of lists (such as coming from a cycle basis)
        if edge_list:
            if not isinstance(edge_list[0], list):
                edge_list = [edge_list]

        self.delete_edgelists(edge_list=edge_list)
        # draw requested edges in new style
        for list_of_edges in edge_list:
            nx.draw_networkx_edges(self.nxGraph, self.positions, ax=self.ax, edgelist=list_of_edges, width=line_width,
                                   edge_color=edge_color)

    def color_nodelist(self, node_list=[], color='b', line_width=2):
        node_color = []
        for node in node_list:
            node_color += [self.nx_palette[self.type_color[self.nxGraph.nodes[node]['type']]]]
        nx.draw_networkx_nodes(self.nxGraph, nodelist=node_list, pos=self.positions, ax=self.ax,
                               node_size=self.nx_options['node_size'], node_color=node_color,
                               linewidths=line_width, edgecolors=color)

    def delete_edgelists(self, edge_list=[]):
        # to unify handling, convert to a list of lists (such as coming from a cycle basis)
        if edge_list:
            if not isinstance(edge_list[0], list):
                edge_list = [edge_list]

        untouched_edges = set([frozenset(e) for e in self.nxGraph.edges()])
        for list_of_edges in edge_list:
            untouched_edges = untouched_edges - set([frozenset(e) for e in list_of_edges])
        remaining_edges = [tuple(x) for x in untouched_edges]

        self.remove_all_edges()
        # redraw what is left in old style
        nx.draw_networkx_edges(self.nxGraph, self.positions, ax=self.ax, edgelist=remaining_edges, **self.nx_options)

    def delete_nodelist(self, node_list=[]):
        untouched_nodes = set(n for n in self.nxGraph.nodes()) - set(n for n in node_list)
        remaining_nodes = [x for x in untouched_nodes]

        self.ax.cla()  # clear the whole figure
        node_color = []
        for node in remaining_nodes:
            node_color += [self.nx_palette[self.type_color[self.nxGraph.nodes[node]['type']]]]
        nx.draw_networkx_nodes(self.nxGraph, nodelist=remaining_nodes, pos=self.positions,
                               ax=self.ax, node_size=self.nx_options['node_size'], node_color=node_color)

        # remove the edges incident on the removed nodes
        e_to_delete = []
        for node in node_list:
            e_to_delete += list(self.nxGraph.edges(node))
        edges = set([frozenset(e) for e in self.nxGraph.edges()]) - set([frozenset(e) for e in e_to_delete])
        remaining_edges = [tuple(x) for x in edges]
        nx.draw_networkx_edges(self.nxGraph, self.positions, ax=self.ax, edgelist=remaining_edges, **self.nx_options)

    def remove_all_edges(self):
        for artist in self.ax.get_children():
            if isinstance(artist, artcoll.LineCollection):
                artist.remove()

    # outdated at this time...
    def html_render(self, filename='', node_size=15, line_width=1, cycle=[]):
        """
        Render a networkx graph with plotly.

        :param line_width:
        :param cycle:
        :param filename:
        :param node_size:
        :param palette:
        :return:
        """
        edge_x = []
        edge_y = []
        for edge in self.nxGraph.edges():
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
            line=dict(width=line_width, color='gray'),
            hoverinfo='none',
            showlegend=False,
            mode='lines')

        cyc_edge_trace = None
        if cycle:
            cyc_edge_x = []
            cyc_edge_y = []
            for item in cycle:
                x0, y0 = self.positions[item[0]]
                x1, y1 = self.positions[item[1]]
                cyc_edge_x.append(x0)
                cyc_edge_x.append(x1)
                cyc_edge_x.append(None)
                cyc_edge_y.append(y0)
                cyc_edge_y.append(y1)
                cyc_edge_y.append(None)

            # edges as lines in a scatter trace
            cyc_edge_trace = go.Scatter(
                x=cyc_edge_x, y=cyc_edge_y,
                line=dict(width=(line_width * 1.5), color='red'),
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
            clr[type] = i % len(self.html_palette)
            i += 1

        for node in self.nxGraph.nodes():
            x, y = self.positions[node]
            node_x[self.nxGraph.nodes[node]['type']].append(x)
            node_y[self.nxGraph.nodes[node]['type']].append(y)

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
                    color=self.html_palette[clr[type]],
                    size=node_size,
                    line_width=0
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
        if cyc_edge_trace:
            fig.add_trace(cyc_edge_trace)
        for type in self.composition.keys():
            fig.add_trace(node_trace[type])

        config = dict({'scrollZoom': True,
                       'displayModeBar': True,
                       'displaylogo': False,
                       'modeBarButtonsToAdd': ['drawline', 'drawopenpath', 'drawclosedpath',
                                               'drawcircle', 'drawrect', 'eraseshape']
                       }
                      )

        if filename:
            if filename.endswith('.html'):
                fig.write_html(filename)
            else:
                fig.write_image(filename)
        else:
            fig.show(config=config)


# to be overhauled...
def show_ranked_complexes(snapshot, sort='size', cutoff=3, cols=3, rows=1, prog='neato'):
    """
    Display the ranked complexes of a snapshot.

    :param snapshot: a snapshot object
    :param sort: 'size' (default) or 'count'
    :param cols: # of columns of plots
    :param rows: # of rows of plots
    :param prog: layout program
    :param cutoff: size or count cutoff
    """
    ranking = []
    if sort == 'size':
        ranking = sorted(snapshot.complexes, key=lambda x: x.size, reverse=True)
    elif sort == 'count':
        ranking = sorted(snapshot.complexes, key=lambda x: x.count, reverse=True)

    plt.figure(frameon=False, figsize=(5., 5.))
    i = 1
    for c in ranking[0:cutoff]:
        plt.subplot(rows, cols, i)
        plt.title(sort + ' ' + str(c.size))
        r = Renderer(c, prog=prog)
        r.render()
        i += 1
    plt.show()


if __name__ == '__main__':
    import kappathings as kt
    import re

    # usage scenarios

    # line = open('TestData/bigly.ka', 'r').read()
    # remove newlines that might occur in the file
    # line = re.sub(r'\n+', ' ', line)
    # create a KappaComplex with whatever assignment of node identifiers arises
    # (that's the normalized=False flag).
    line = "A(l[19] r[.] p[2]), A(l[53] r[19] p[42]), A(l[37] r[53] p[45]), A(l[.] r[37] p[29]), P(a1[3] a2[51] a3[" \
           "29] d[.]), A(l[20] r[.] p[3]), A(l[.] r[20] p[27]), P(a1[27] a2[.] a3[.] d[44]), P(a1[.] a2[.] a3[.] d[" \
           "44]), A(l[.] r[.] p[51]), P(a1[14] a2[.] a3[45] d[22]), A(l[24] r[.] p[14]), A(l[.] r[24] p[50]), " \
           "P(a1[50] a2[.] a3[30] d[.]), A(l[13] r[16] p[30]), A(l[.] r[13] p[21]), P(a1[21] a2[35] a3[1] d[9]), " \
           "A(l[.] r[.] p[35]), A(l[.] r[26] p[1]), A(l[26] r[.] p[32]), P(a1[32] a2[.] a3[.] d[9]), A(l[16] r[54] p[" \
           ".]), A(l[54] r[38] p[40]), A(l[38] r[.] p[7]), P(a1[.] a2[.] a3[7] d[52]), P(a1[18] a2[.] a3[4] d[52]), " \
           "A(l[.] r[.] p[18]), A(l[23] r[.] p[4]), A(l[49] r[23] p[47]), A(l[28] r[49] p[11]), A(l[6] r[28] p[31]), " \
           "A(l[.] r[6] p[12]), P(a1[.] a2[12] a3[.] d[.]), P(a1[.] a2[31] a3[.] d[15]), P(a1[.] a2[5] a3[.] d[15]), " \
           "A(l[.] r[.] p[5]), P(a1[8] a2[.] a3[11] d[36]), A(l[.] r[39] p[8]), A(l[39] r[.] p[25]), P(a1[41] a2[.] " \
           "a3[25] d[33]), A(l[.] r[.] p[41]), P(a1[.] a2[43] a3[.] d[33]), A(l[.] r[46] p[43]), A(l[46] r[.] p[10]), " \
           "P(a1[10] a2[.] a3[.] d[22]), P(a1[34] a2[.] a3[.] d[36]), A(l[.] r[.] p[34]), P(a1[.] a2[47] a3[.] d[.]), " \
           "P(a1[.] a2[40] a3[42] d[17]), P(a1[48] a2[2] a3[.] d[17]), A(l[.] r[.] p[48]) "
    line2 = "A(l[.] r[4] p[1]), A(l[4] r[.] p[3]), P(a1[3] a2[1] a3[.] d[2]), P(a1[.] a2[.] a3[.] d[2])"
    c1 = kt.KappaComplex(line)

    # r = Renderer(c1)
    # r.render(labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')
    # show()

    c2 = kt.KappaComplex(line2)

    plt.ion()
    canvas = Canvas(2, 2)
    r = Renderer(c1, canvas)
    r2 = Renderer(c2, canvas)
    r.render((1, 2), labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray', title="t=0.37")
    # input()
    r.render((2, 1), labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')
    # input()
    r2.render((2, 1), labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')
    # input()
    r.render((2, 2), labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')
    # input()

    # r.render(node_size=10, labels='none')
    # r.show()
    # r.html_render()
    # r.html_render(filename='complex.html')
    # r.nx_render(c1, prog='sfdp', node_size=80, font_size=4)
    # r.nx_render(c1, prog='fdp', node_size=80, font_size=4)
    # r.nx_render(c1, prog='neato', node_size=80, font_size=4)
    # r.nx_render(c1, prog='twopi', node_size=80, font_size=4)
    # r.nx_render(c1, prog='dot', node_size=80, font_size=4)
    # r.nx_render(c1, prog='circo', node_size=80, font_size=4)

    # app = QApplication(sys.argv)
    # canvas = Canvas(2, 2, qt5=True)
    # r = Renderer(c1, canvas)
    # r2 = Renderer(c2, canvas)
    # r.render((1, 2), labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')
    # r.render((2, 1), labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')
    # r2.render((2, 1), labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')
    # r.render((2, 2), labels='no', node_size=20, font_size=9, line_width=1, edge_color='gray')
    # canvas.qt.show()
    # sys.exit(app.exec_())
