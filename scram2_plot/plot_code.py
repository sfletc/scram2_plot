
import time
import csv
import sys
import numpy
from pylab import *  # @UnusedWildImport
import matplotlib.pyplot as plt  # @Reimport
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.io import output_notebook
from bokeh.models import HoverTool
from collections import OrderedDict


class DNA(object):
    dna_alphabet = set("AGCTN")

    def __init__(self, sequence):
        self.sequence = sequence.upper()

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, key):
        return self.sequence[key]

    def __hash__(self):
        return hash(self.sequence)

    def __repr__(self):
        return self.sequence

    def __eq__(self, other):
        return self.sequence == other.sequence


def import_scram2_den(in_file):
    """
    Import a SCRAM2 csv file to a dictionary
    """
    alignments = {}
    srna_len=0
    with open(in_file, 'r') as f:
        first_line = True
        for line in f:
            if first_line:
                first_line = False
            else:
                line = line.strip().rsplit(',', 5)
                srna_len = len(line[2])
                if line[0] not in alignments:
                    alignments[line[0]] = [(int(line[1]), DNA(line[2]), int(line[3]), float(line[4]),float(line[5]))]
                else:
                    alignments[line[0]].append((int(line[1]), DNA(line[2]), int(line[3]), float(line[4]),float(line[5])))
    return alignments, srna_len


def extract_header_alignment(header, alignments):
    """
    With a provided complete header, extract the alignment and process to correct format for fill in zeros
    """
    try:
        extracted_alignments = alignments[header]
        sorted_fwd_alignment = []
        sorted_rvs_alignment = []
        aln_count = 0.0
        ref_len = 0
        for alignment in extracted_alignments:
            ref_len = alignment[0]
            if alignment[3] > 0:
                sorted_fwd_alignment.append((alignment[2], alignment[3], alignment[4]))
                aln_count += alignment[3]
            elif alignment[3] < 0:
                sorted_rvs_alignment.append((alignment[2], alignment[3], alignment[4]))
                aln_count -= alignment[3]
    except:
        sorted_fwd_alignment = []
        sorted_rvs_alignment = []
        aln_count = 0.0
        ref_len = 0
    return [sorted_fwd_alignment, sorted_rvs_alignment, aln_count], ref_len




def fill_in_zeros_se(fwd_rvs_align_list, ref_len, nt):
    """
    Generate alignment counts for every nucleotide in the reference
    :param fwd_rvs_align_list:  list of sorted forwards and reverse alignments
    :param ref_len: number of nucleotides in the reference sequence (int)
    :param nt: length of the aligned reads (int)
    :return: reference_x_axis ([0,0,...] (list(int)) - length of refseq seq,
             fwd_alignment_y_axis [2,4,5.2,6,....] (list(float)) - sense strand alignment count (positive),
             fwd_rvs_align_list [-3,-4,-5.6,...] (list(float)) - antisense strand alignment count (negative)
    """
    sorted_fwd_alignment = fwd_rvs_align_list[0]
    sorted_rvs_alignment = fwd_rvs_align_list[1]

    fwd_alignment_y_axis_upper = [0] * ref_len
    fwd_alignment_y_axis_lower = [0] * ref_len
    revs_alignment_y_axis_upper = [0] * ref_len
    revs_alignment_y_axis_lower = [0] * ref_len

    reference_x_axis = list(range(0, ref_len))
    # Note alignment position for graphing is in the centre of the read (and not the 5' end)
    try:
        for i in sorted_fwd_alignment:
            fwd_alignment_y_axis_upper[i[0]] = i[1] + i[2]
            fwd_alignment_y_axis_lower[i[0]] = i[1] - i[2]
        for i in sorted_rvs_alignment:
            revs_alignment_y_axis_upper[i[0]] = i[1] + i[2]
            revs_alignment_y_axis_lower[i[0]] = i[1] - i[2]
    except:
        pass
    # #Coverage per nucleotide instead - maybe use?
    #     for i in sorted_fwd_alignment:
    #         for j in range(nt):
    #             fwd_alignment_y_axis[i[0]+j]+=i[1]
    #     for i in sorted_rvs_alignment:
    #         for j in range(nt):
    #             revs_alignment_y_axis[i[0]-j]+=i[1]

    return reference_x_axis, fwd_alignment_y_axis_upper, fwd_alignment_y_axis_lower, \
           revs_alignment_y_axis_upper, revs_alignment_y_axis_lower

def _smoothed_for_plot_se(graph_processed, smooth_win_size):
    """
    Return fwd and rvs smoothed profiles
    :param graph_processed:
    :param smooth_win_size:
    :return:
    """
    y_fwd_smoothed_upper = smooth(numpy.array(graph_processed[1]),
                            smooth_win_size, window='blackman')
    y_fwd_smoothed_lower = smooth(numpy.array(graph_processed[2]),
                            smooth_win_size, window='blackman')
    y_rvs_smoothed_upper = smooth(numpy.array(graph_processed[3]),
                            smooth_win_size, window='blackman')
    y_rvs_smoothed_lower = smooth(numpy.array(graph_processed[4]),
                            smooth_win_size, window='blackman')
    return y_fwd_smoothed_upper, y_fwd_smoothed_lower, y_rvs_smoothed_upper,  y_rvs_smoothed_lower


def multi_header_plot(search_terms, in_files, cutoff, plot_y_lim, win, pub):
    """

    :param search_terms:
    :param cutoff:
    :param f:
    :param no_display:
    :param search_terms:
    :param in_21:
    :param in_22:
    :param in_24:
    :param cutoff:
    :param plot_y_lim:
    :param pub:
    :param f:
    :param no_display:
    :return:
    """

    all_files_present = 0
# for file_name in in_files:
    try:
# if file_name.strip().split(".")[-2][-2:] == "21":
        print("Loading 21 nt Alignment File\n")
        in_21, _ = import_scram2_den(in_files+"_21.csv")
        # all_files_present +=1
        # elif file_name.strip().split(".")[-2][-2:] == "22":
        print("Loading 22 nt Alignment File\n")
        in_22, _ = import_scram2_den(in_files+"_22.csv")
        # all_files_present +=1
        # elif file_name.strip().split(".")[-2][-2:] == "24":
        print("Loading 24 nt Alignment File\n")
        in_24, _ = import_scram2_den(in_files+"_24.csv")
        # all_files_present +=1
    except:
        print("\n21nt, 22nt and 24nt alignment files are required to proceed")
        sys.exit()
    substring = " ".join(search_terms)

    print("Extracting headers\n")
    all_keys=set()
    for nt in [in_21, in_22, in_24]:
        for header in nt.keys():
            all_keys.add(header)
    for header in all_keys:
        if substring.lower() in header.lower():
            alignment_21, a = extract_header_alignment(header, in_21)
            alignment_22, b = extract_header_alignment(header, in_22)
            alignment_24, c = extract_header_alignment(header, in_24)
            if alignment_21[2] >= cutoff or alignment_22[2] >= cutoff or alignment_24[2] >= cutoff:
                print (header)
                ref_len = max(a, b, c)
                if win ==0:
                    win = int(ref_len / 30)
                if win % 2 != 0 or win == 0: win += 1
                graph_processed_21 = fill_in_zeros_se(alignment_21, ref_len, 21)
                graph_processed_22 = fill_in_zeros_se(alignment_22, ref_len, 22)
                graph_processed_24 = fill_in_zeros_se(alignment_24, ref_len, 24)
                x_ref = graph_processed_21[0]
                y_fwd_smoothed_upper_21, y_fwd_smoothed_lower_21, y_rvs_smoothed_upper_21, y_rvs_smoothed_lower_21 = _smoothed_for_plot_se(graph_processed_21, win)
                y_fwd_smoothed_upper_22, y_fwd_smoothed_lower_22, y_rvs_smoothed_upper_22, y_rvs_smoothed_lower_22 = _smoothed_for_plot_se(graph_processed_22, win)
                y_fwd_smoothed_upper_24, y_fwd_smoothed_lower_24, y_rvs_smoothed_upper_24, y_rvs_smoothed_lower_24 = _smoothed_for_plot_se(graph_processed_24, win)
                den_multi_plot_21_22_24_se(x_ref, y_fwd_smoothed_upper_21, y_fwd_smoothed_lower_21,
                                           y_rvs_smoothed_upper_21, y_rvs_smoothed_lower_21,
                                           y_fwd_smoothed_upper_22, y_fwd_smoothed_lower_22,
                                           y_rvs_smoothed_upper_22, y_rvs_smoothed_lower_22,
                                           y_fwd_smoothed_upper_24, y_fwd_smoothed_lower_24,
                                           y_rvs_smoothed_upper_24, y_rvs_smoothed_lower_24,
                                           header, plot_y_lim, pub)

def single_header_plot(search_terms, in_file, cutoff, plot_y_lim, win, pub):
    """

    :param f:
    :param no_display:
    :param pub:
    :param search_terms:
    :param in_file:
    :param cutoff:
    :param plot_y_lim:
    :return:
    """
    substring = " ".join(search_terms)
    print("Loading Alignment Files\n")
    in_x, nt = import_scram2_den(in_file)
    print("Extracting headers\n")
    for header in in_x.keys():
        if substring.lower() in header.lower():
            alignment_x, ref_len = extract_header_alignment(header, in_x)
            if alignment_x[2] >= cutoff:
                print(header)
                if win ==0:
                    win = int(ref_len / 30)
                if win % 2 != 0 or win == 0: win += 1
                graph_processed_x = fill_in_zeros_se(alignment_x, ref_len, nt)
                x_ref = graph_processed_x[0]
                y_fwd_smoothed_upper, y_fwd_smoothed_lower, y_rvs_smoothed_upper, y_rvs_smoothed_lower = _smoothed_for_plot_se(graph_processed_x, win)
                den_plot_se(x_ref, y_fwd_smoothed_upper, y_fwd_smoothed_lower, y_rvs_smoothed_upper, y_rvs_smoothed_lower,
                            nt, header, plot_y_lim, pub)



def smooth(x, window_len, window='hamming'):
    """
    Smoothing function from scipy cookbook
    :param x:
    :param window_len:
    :param window:
    :return:
    """

    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 6:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = numpy.r_[x[window_len - 1:0:-1], x, x[-1:-window_len:-1]]
    if window == 'flat':  # moving average
        w = numpy.ones(window_len, 'd')
    else:
        w = eval('numpy.' + window + '(window_len)')

    y = numpy.convolve(w / w.sum(), s, mode='valid')
    return y[int(window_len / 2 - 1):-int(window_len / 2)]



def den_plot_se(x_ref, y_fwd_smoothed_upper, y_fwd_smoothed_lower, y_rvs_smoothed_upper, y_rvs_smoothed_lower,
                nt, xlab, plot_y_lim, pub):
    """

    :param f:
    :param no_display:
    :param plot_y_lim:
    :param pub:
    :param x_ref:
    :param y_fwd_smoothed_upper:
    :param y_fwd_smoothed_lower:
    :param y_rvs_smoothed_upper:
    :param y_rvs_smoothed_lower:
    :param nt:
    :param xlab:
    :return:
    """
    fig = plt.figure(figsize=(10, 5))
    plt.plot(x_ref, y_fwd_smoothed_upper, color=_nt_colour(nt), label='{0} nt'.format(nt), lw=1, alpha=0.2)
    plt.plot(x_ref, y_fwd_smoothed_lower, color=_nt_colour(nt), lw=1, alpha=0.2)
    plt.fill_between(x_ref, y_fwd_smoothed_upper,y_fwd_smoothed_lower, color=_nt_colour(nt), alpha=0.5)
    plt.plot(x_ref, y_rvs_smoothed_upper, color=_nt_colour(nt), lw=1, alpha=0.2)
    plt.plot(x_ref, y_rvs_smoothed_lower, color=_nt_colour(nt), lw=1, alpha=0.2)
    plt.fill_between(x_ref, y_rvs_smoothed_upper,y_rvs_smoothed_lower, color=_nt_colour(nt), alpha=0.5)
    axhline(y=0)
    if pub:
        _pub_plot()
    else:
        xlabel(xlab)
        ylabel('Reads per million reads')
        plt.legend(loc='best', fancybox=True, framealpha=0.5)
    _generate_profile(plot_y_lim)
    _generate_profile(plot_y_lim)


def den_multi_plot_21_22_24_se(x_ref, y_fwd_smoothed_upper_21, y_fwd_smoothed_lower_21,
                           y_rvs_smoothed_upper_21, y_rvs_smoothed_lower_21,
                           y_fwd_smoothed_upper_22, y_fwd_smoothed_lower_22,
                           y_rvs_smoothed_upper_22, y_rvs_smoothed_lower_22,
                           y_fwd_smoothed_upper_24, y_fwd_smoothed_lower_24,
                           y_rvs_smoothed_upper_24, y_rvs_smoothed_lower_24,
                           header, plot_y_lim, pub):
    fig = plt.figure(figsize=(10, 5))
    #21
    plt.plot(x_ref, y_fwd_smoothed_upper_21, color=_nt_colour(21), label='{0} nt'.format(21), lw=1, alpha=0.2)
    plt.plot(x_ref, y_fwd_smoothed_lower_21, color=_nt_colour(21), lw=1, alpha=0.2)
    plt.fill_between(x_ref, y_fwd_smoothed_upper_21,y_fwd_smoothed_lower_21, color=_nt_colour(21), alpha=0.5)
    plt.plot(x_ref, y_rvs_smoothed_upper_21, color=_nt_colour(21), lw=1, alpha=0.2)
    plt.plot(x_ref, y_rvs_smoothed_lower_21, color=_nt_colour(21), lw=1, alpha=0.2)
    plt.fill_between(x_ref, y_rvs_smoothed_upper_21,y_rvs_smoothed_lower_21, color=_nt_colour(21), alpha=0.5)
    #22
    plt.plot(x_ref, y_fwd_smoothed_upper_22, color=_nt_colour(22), label='{0} nt'.format(22), lw=1, alpha=0.2)
    plt.plot(x_ref, y_fwd_smoothed_lower_22, color=_nt_colour(22), lw=1, alpha=0.2)
    plt.fill_between(x_ref, y_fwd_smoothed_upper_22,y_fwd_smoothed_lower_22, color=_nt_colour(22), alpha=0.5)
    plt.plot(x_ref, y_rvs_smoothed_upper_22, color=_nt_colour(22), lw=1, alpha=0.2)
    plt.plot(x_ref, y_rvs_smoothed_lower_22, color=_nt_colour(22), lw=1, alpha=0.2)
    plt.fill_between(x_ref, y_rvs_smoothed_upper_22,y_rvs_smoothed_lower_22, color=_nt_colour(22), alpha=0.5)
    #24
    plt.plot(x_ref, y_fwd_smoothed_upper_24, color=_nt_colour(24), label='{0} nt'.format(24), lw=1, alpha=0.2)
    plt.plot(x_ref, y_fwd_smoothed_lower_24, color=_nt_colour(24), lw=1, alpha=0.2)
    plt.fill_between(x_ref, y_fwd_smoothed_upper_24,y_fwd_smoothed_lower_24, color=_nt_colour(24), alpha=0.5)
    plt.plot(x_ref, y_rvs_smoothed_upper_24, color=_nt_colour(24), lw=1, alpha=0.2)
    plt.plot(x_ref, y_rvs_smoothed_lower_24, color=_nt_colour(24), lw=1, alpha=0.2)
    plt.fill_between(x_ref, y_rvs_smoothed_upper_24,y_rvs_smoothed_lower_24, color=_nt_colour(24), alpha=0.5)
    
    
    axhline(y=0)
    if pub:
        _pub_plot()
    else:
        xlabel(header)
        ylabel('Reads per million reads')
        plt.legend(loc='best', fancybox=True, framealpha=0.5)
    _generate_profile(plot_y_lim)
    
def _generate_profile(plot_y_lim):
    """
    Generate profile
    :param file_fig: output plot to pdf (bool)
    :param file_name: output filename (str)
    :param onscreen: show plot on screen (bool)
    :param plot_y_lim: + / - y-axis limit (int)
    """
    if plot_y_lim != 0:
        ylim(-plot_y_lim, plot_y_lim)
    #fig1 = plt.gcf()
    plt.show()
    #plt.close(fig1)

def _pub_plot():
    """
    Remove axis, labels, legend from plot
    """
    plt.tick_params(
        axis='both',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom='on',  # ticks along the bottom edge are off
        top='off',
        right='on',
        left='on',  # ticks along the top edge are off
        labelbottom='off',
        labelleft='off',
        labelright='off',
        labelsize=15)  # labels along the bottom edge are off
    _clear_frame()


def _clear_frame(ax=None):
    """
    Removes frame for publishing plots
    """
    if ax is None:
        ax = plt.gca()
    ax.xaxis.set_visible(True)
    ax.yaxis.set_visible(True)
    for spine in ax.spines.values():
        spine.set_visible(False)


def _nt_colour(nt):
    """
    Set default colours for 21, 22 and 24 nt sRNAs
    :param nt: aligned read length (int)
    :return: colour code (str)
    """
    if nt == 21:
        col = '#00CC00'
    elif nt == 22:
        col = '#FF3399'
    elif nt == 24:
        col = '#3333FF'
    else:
        col = 'black'
    return col



def cdp_plot_bokeh(file_name, seq1, seq2, plot_type, browser):
    output_notebook()
    if browser:
        output_file("plot.html")
    else:
        output_notebook()
    first_line=True
    x_vals_line = []
    x_vals_point = []
    ellipse_width=[]
    max_x=0.0
    y_vals_line = []
    y_vals_point = []
    ellipse_height=[]
    header=[]
    line_header=[]
    max_y=0.0
    try:
        nt = int(file_name.strip().split('.')[-2][-2:])
    except:
        nt = 20 # this is just for plot colour, so makes black if nt can't be parsed from filename
    with open(file_name) as csvfile:
        line_reader=csv.reader(csvfile)
        for line in line_reader:
            if first_line:
                first_line=False
            else:
                #calc max value
                if float(line[-4]) > max_x:
                    max_x = float(line[-4])
                if float(line[-2]) > max_y:
                    max_y = float(line[-2])
                #line
                line[0]=line[0].strip()
                x_se = [float(line[-4])-float(line[-3]), float(line[-4])+float(line[-3])]
                y_se = [float(line[-2])-float(line[-1]), float(line[-2])+float(line[-1])]

                x_vals_line.append(x_se)
                y_vals_line.append([float(line[-2]),float(line[-2])])
                x_vals_line.append([float(line[-4]),float(line[-4])])
                y_vals_line.append(y_se)
                #point
                x_vals_point.append(float(line[-4]))
                y_vals_point.append(float(line[-2]))
                #height and width
                ellipse_width.append(2 * float(line[-3]))
                ellipse_height.append(2 * float(line[-1]))
                header.append(line[0])
                line_header=line_header+[line[0]]+[line[0]]

    _max = max([max_x,max_y])  # sets up max x and y scale values
    log_max = _max + _max / 2
    csvfile.close()

    #Interactive
    if plot_type=="log" or plot_type=="all":
        log_dot_plot(header, log_max, nt, seq1, seq2, x_vals_point, y_vals_point)
    if plot_type == "log_error" or plot_type == "all":
        log_se_plot(header, line_header, log_max, nt, seq1, seq2, x_vals_line, x_vals_point, y_vals_line, y_vals_point)

    if plot_type == "linear" or plot_type == "all":
        linear_ellipse_plot(_max, ellipse_height, ellipse_width, header, nt, seq1, seq2, x_vals_point, y_vals_point)


def linear_ellipse_plot(_max, ellipse_height, ellipse_width, header, nt, seq1, seq2, x_vals_point, y_vals_point):
    # elipse
    print("Interactive linear plot")
    hover = HoverTool(
        tooltips=[
            ("(x,y)", "($x, $y)"),
            ("header", "@Desc")
        ]
    )
    linear_max = _max + _max / 5
    p = figure(plot_width=600, plot_height=600,
               x_range=(1, linear_max), y_range=(1, linear_max),
               toolbar_location="above", tools=[hover, 'save', 'box_zoom', 'reset'])
    source_point = ColumnDataSource(data=OrderedDict(x=x_vals_point,
                                                     y=y_vals_point,
                                                     Desc=header, )
                                    )
    source_ellipse = ColumnDataSource(data=OrderedDict(x=x_vals_point,
                                                       y=y_vals_point, width=ellipse_width, height=ellipse_height,
                                                       Desc=header, )
                                      )
    p.ellipse('x', 'y', 'width', 'height', source=source_ellipse, color=_nt_colour(nt), alpha=0.2)
    p.circle('x', 'y', source=source_point, size=4, color=_nt_colour(nt))
    p.line([1, linear_max], [1, linear_max])
    p.xaxis.axis_label = seq1
    p.yaxis.axis_label = seq2
    show(p)


def log_se_plot(header, line_header, log_max, nt, seq1, seq2, x_vals_line, x_vals_point, y_vals_line, y_vals_point):
    # Std Error bars
    print("Interactive log plot with se bars")
    hover = HoverTool(
        tooltips=[
            ("(x,y)", "($x, $y)"),
            ("header", "@Desc")
        ]
    )
    p = figure(plot_width=600, plot_height=600,
               x_axis_type="log", y_axis_type="log",
               x_range=(0.1, log_max), y_range=(0.1, log_max),
               toolbar_location="above", tools=[hover, 'save', 'box_zoom', 'reset'])
    source_point = ColumnDataSource(data=OrderedDict(x=x_vals_point,
                                                     y=y_vals_point, Desc=header, )
                                    )
    source_line = ColumnDataSource(data=OrderedDict(x=x_vals_line,
                                                    y=y_vals_line, Desc=line_header, )
                                   )
    p.multi_line('x', 'y', source=source_line, color=_nt_colour(nt), alpha=0.5, )
    p.circle('x', 'y', source=source_point, size=2, color=_nt_colour(nt))
    p.line([0.1, log_max], [0.1, log_max])
    p.xaxis.axis_label = seq1
    p.yaxis.axis_label = seq2
    show(p)


def log_dot_plot(header, log_max, nt, seq1, seq2, x_vals_point, y_vals_point):
    print("Interactive log plot")
    hover = HoverTool(
        tooltips=[
            ("(x,y)", "($x, $y)"),
            ("header", "@Desc")
        ]
    )
    p = figure(plot_width=600, plot_height=600,
               x_axis_type="log", y_axis_type="log",
               x_range=(0.1, log_max), y_range=(0.1, log_max),
               toolbar_location="above", tools=[hover, 'save', 'box_zoom', 'reset'])
    source = ColumnDataSource(data=OrderedDict(x=x_vals_point,
                                               y=y_vals_point,
                                               Desc=header, )
                              )
    p.circle('x', 'y', source=source, size=5, color=_nt_colour(nt))
    p.line([0.1, log_max], [0.1, log_max])
    p.xaxis.axis_label = seq1
    p.yaxis.axis_label = seq2
    show(p)