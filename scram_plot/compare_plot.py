from pylab import *  # @UnusedWildImport
import matplotlib.pyplot as plt  # @Reimport
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.io import output_notebook
from bokeh.models import HoverTool
from collections import OrderedDict
import csv
import profile_plot as pp
import math
import os.path

def compare_plot(file_prefix, nt_list, seq1, seq2, plot_type, browser, save_plot, pub, fig_size, xylim):
    """
    Compare plot
    :param file_prefix: path/to/file prefix
    :param nt_list: read length list to plot
    :param seq1: x seq name
    :param seq2: y seq name
    :param plot_type: plot type
    :param browser: output to html
    :param save_plot: bool to save plot
    :param pub: bool to remove axes and legend
    :param fig_size: figure dimensions
    :param xylim: x/y axi limit
    """
    #try:
    if nt_list[0]=="mir":
        fname = file_prefix + "_miR.csv"
        if os.path.isfile(fname):
            try:
                _format_compare_data(fname, "mir", browser, plot_type, pub, save_plot,
                                     seq1, seq2, fig_size, xylim)
            except:
                print("\nCannot load and process {}".format(fname))
                sys.exit()
        else:
            print("\n{} does not exist at this location".format(fname))
            sys.exit()

    else:
        for nt in nt_list:
            fname = "{0}_{1}.csv".format(file_prefix, nt)
            if os.path.isfile(fname):
                try:
                    _format_compare_data(fname, int(nt), browser, plot_type, pub, save_plot,
                                         seq1, seq2, fig_size, xylim)
                except:
                    print("\nCannot load and process {}".format(fname))
                    sys.exit()
            else:
                print("\n{} does not exist at this location".format(fname))
                sys.exit()


def _format_compare_data(file_name, nt, browser, plot_type, pub, save_plot, seq1, seq2, fig_size, xylim):
    if "/" in file_name:
        file_path = file_name.rsplit('/', 1)[0]
    else:
        file_path="."

    if browser:
        output_file(file_path + '/{0}_{1}_{2}.html'.format(seq1, seq2, nt))
    else:
        output_notebook(hide_banner=True)
    first_line = True
    x_vals_line = []
    x_vals_point = []
    xerr = []
    max_x = 0.0
    y_vals_line = []
    y_vals_point = []
    header = []
    yerr = []
    max_y = 0.0
    log_max, max_x = _extract_alignment_data(file_name, first_line, header, max_x, max_y, x_vals_line, x_vals_point,
                                             xerr, xylim, y_vals_line, y_vals_point, yerr)
    # Interactive
    #Hack as bokeh tooltip text wrapping for large x values not working properly
    for plot_point in range(len(header)):
        if math.log10(x_vals_point[plot_point]+0.0001)> 0.35*math.log10(max_x):
            header[plot_point]=header[plot_point][:40]

    _plot_type(fig_size, file_path, header, log_max, nt, plot_type, pub, save_plot, seq1, seq2, x_vals_line,
               x_vals_point, xerr, y_vals_line, y_vals_point, yerr)


def _plot_type(fig_size, file_path, header, log_max, nt, plot_type, pub, save_plot, seq1, seq2, x_vals_line,
               x_vals_point, xerr, y_vals_line, y_vals_point, yerr):
    if plot_type == "log" or plot_type == "all":
        _plot_compare_plot(file_path, header, log_max, nt, seq1, seq2, [], x_vals_point, [], y_vals_point, [], [],
                           save_plot,
                           pub, fig_size)
    if plot_type == "log_error" or plot_type == "all":
        _plot_compare_plot(file_path, header, log_max, nt, seq1, seq2, x_vals_line, x_vals_point, y_vals_line,
                           y_vals_point, xerr, yerr, save_plot, pub, fig_size)


def _extract_alignment_data(file_name, first_line, header, max_x, max_y, x_vals_line, x_vals_point, xerr, xylim,
                            y_vals_line, y_vals_point, yerr):
    with open(file_name) as csvfile:
        line_reader = csv.reader(csvfile)
        for line in line_reader:
            if first_line:
                first_line = False
            else:
                # calc max value
                if float(line[-4]) > max_x:
                    max_x = float(line[-4])
                if float(line[-2]) > max_y:
                    max_y = float(line[-2])
                # line
                line[0] = line[0].strip()
                x_se = [float(line[-4]) - float(line[-3]), float(line[-4]) + float(line[-3])]
                y_se = [float(line[-2]) - float(line[-1]), float(line[-2]) + float(line[-1])]
                xerr.append(float(line[-3]))
                x_vals_line.append(x_se)
                y_vals_line.append([float(line[-2]), float(line[-2])])
                x_vals_line.append([float(line[-4]), float(line[-4])])
                y_vals_line.append(y_se)
                yerr.append(float(line[-1]))
                # point
                x_vals_point.append(float(line[-4]))
                y_vals_point.append(float(line[-2]))

                header.append(line[0])
    if xylim == "auto":
        _max = max([max_x, max_y])  # sets up max x and y scale values
        log_max = _max + _max / 2
    else:
        log_max = int(xylim)
    csvfile.close()
    return log_max, max_x


def _plot_compare_plot(file_path, header, log_max, nt, seq1, seq2, x_vals_line, x_vals_point, y_vals_line, y_vals_point,
                       xerr, yerr, save_plot, pub_plot, fig_size = 4):
    # Std Error bars
    hover = HoverTool(
        tooltips=[
            ("(x,y)", "($x, $y)"),
            ("header", "@Desc")
        ],
        names=["circle", ]
    )
    p = figure(plot_width=600, plot_height=600,
               x_axis_type="log", y_axis_type="log",
               x_range=(0.1, log_max), y_range=(0.1, log_max),
               toolbar_location="above", tools=[hover, 'save', 'box_zoom', 'reset'])
    source_point = ColumnDataSource(data=OrderedDict(x=x_vals_point,
                                                     y=y_vals_point, Desc=header, )
                                    )
    if nt == "mir":
        p.circle('x', 'y', name="circle", source=source_point, size=3, color=pp._nt_colour(nt),
                 legend="microRNA")
    else:
        p.circle('x', 'y', name="circle", source=source_point, size=3, color=pp._nt_colour(nt), legend="{0} nt".format(
        nt))
    p.legend.location = "top_left"
    p.line([0.1, log_max], [0.1, log_max])
    if xerr != []:
        p.multi_line(xs=x_vals_line, ys=y_vals_line, color=pp._nt_colour(nt), alpha=0.5, )
    p.xaxis.axis_label = seq1
    p.yaxis.axis_label = seq2
    show(p)
    if save_plot:
        _compare_plot_to_file(fig_size, file_path, log_max, nt, pub_plot, seq1, seq2, x_vals_point, xerr, y_vals_point, yerr)


def _compare_plot_to_file(fig_size, file_path, log_max, nt, pub_plot, seq1, seq2, x_vals_point, xerr, y_vals_point, yerr):
    fig = plt.figure(figsize=(fig_size, fig_size))
    if xerr != []:
        plt.errorbar(x_vals_point, y_vals_point, xerr=xerr, yerr=yerr, capsize=0, ls='none', color=pp._nt_colour(
            nt),
                     elinewidth=0.5)
    plt.plot([0.1, log_max], [0.1, log_max], alpha=0.9, linewidth=1)
    if nt == "mir":
        plt.scatter(x_vals_point, y_vals_point, color=pp._nt_colour(nt), s=3, label="microRNA".format(nt))
    else:
        plt.scatter(x_vals_point, y_vals_point, color=pp._nt_colour(nt), s=3, label="{0} nt".format(nt))
    plt.xlim([0.1, log_max])
    plt.ylim([0.1, log_max])
    plt.xscale('log')
    plt.yscale('log')
    if pub_plot:
        pp._pub_plot()
        plt.minorticks_off()
    else:
        plt.grid(linestyle='-', alpha=0.2)
        plt.xlabel(seq1)
        plt.ylabel(seq2)
        plt.legend()
    plt.savefig(file_path + '/{0}_{1}_{2}.png'.format(seq1, seq2, nt), dpi=300)
