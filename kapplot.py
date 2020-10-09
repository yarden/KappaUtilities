# Walter Fontana at 20.09.2020

import plotly.graph_objects as go
import matplotlib.pyplot as plt


def show(pdf=''):
    """
    :param pdf: if not '', write to pdf file
    """
    if pdf == '':
        plt.show(block=False)
    else:
        plt.savefig(pdf)


class XY_plot:
    def __init__(self, df, x='', y='', title='', xmajor=0, ymajor=0, params={}):
        """
        :param ymajor: multiple for major y-tick marks (0 for auto)
        :param xmajor: multiple for major x-tick marks (0 for auto)
        :param title: plot title
        :param params: parameter dict to be passed to plot
        :param df: pandas dataframe
        :param y: name of x column
        :param x: name of y column
        """
        self.default_x = 0
        self.default_y = 1
        self.title = ''
        self.parameters = {'linestyle': '',
                           'linewidth': 0.5,
                           'marker': 'o',
                           'label': '',
                           'markersize': 0
                           }

        self.fig, self.ax = plt.subplots()
        self.add(df, x=x, y=y, title=title, xmajor=xmajor, ymajor=ymajor, params=params)

    def add(self, df, x='', y='', title='', xmajor=0, ymajor=0, params={}):
        """
        :param ymajor: multiple for major y-tick marks (0 for auto)
        :param xmajor: multiple for major x-tick marks (0 for auto)
        :param title: plot title
        :param params: parameter dict to be passed to plot
        :param df: pandas dataframe
        :param y: name of x column
        :param x: name of y column
        """
        self.parameters = {**self.parameters, **params}

        if x == '':
            for idx, c in enumerate(df.columns):
                if idx == self.default_x:
                    x = c
                    break
        if y == '':
            for idx, c in enumerate(df.columns):
                if idx == self.default_y:
                    y = c
                    break
        if title == '':
            self.title = f'{x} vs {y}'
        else:
            self.title = title

        self.ax.plot(df[x], df[y], 'o-', **self.parameters)
        if xmajor != 0:
            self.ax.xaxis.set_major_locator(plt.MultipleLocator(xmajor))
        if ymajor != 0:
            self.ax.yaxis.set_major_locator(plt.MultipleLocator(ymajor))
        plt.grid(color='lightgrey')
        # plt.vlines(df[x], 0., df[y], color='b', linestyles='solid', label='')
        # plt.hlines(0, 0., df[x][len(df)-1], color='orange', linestyles='dashed')
        self.ax.set_xlabel(x)
        self.ax.set_ylabel(y)
        self.ax.set_title(title)


if __name__ == '__main__':
    import pandas as pd
    import kappasnap as ks

    snap1 = ks.SnapShot('TestData/snap19.ka')
    sd_df1 = pd.DataFrame(snap1.get_size_distribution(dictionary=True))
    plot = XY_plot(sd_df1, xmajor=2, ymajor=2000, params={'linestyle': '-',
                                                          'label': 'snap19',
                                                          'color': 'r',
                                                          'markerfacecolor': 'r'})
    # show()
    snap2 = ks.SnapShot('TestData/snap98.ka')
    sd_df2 = pd.DataFrame(snap2.get_size_distribution(dictionary=True))
    plot.add(sd_df2, xmajor=2, ymajor=2000, params={'linestyle': '-',
                                                    'label': 'snap98',
                                                    'color': 'g',
                                                    'markerfacecolor': 'g'})
    plot.ax.legend()
    show()
