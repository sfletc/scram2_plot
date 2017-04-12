import csv
import numpy
from pylab import *  # @UnusedWildImport
import matplotlib.pyplot as plt  # @Reimport
from bokeh.plotting import figure, output_file, show, ColumnDataSource
from bokeh.io import output_notebook
from bokeh.models import HoverTool
from collections import OrderedDict


class DNA(object):
    """
    DNA class
    """
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
    :param in_file: path/to/profile string
    :return: alignments dictionary and snra length in the alignment
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
    :param header: reference sequence header string 
    :param alignments: alignments dictionary
    :return: sorted_fwd_alignment, sorted_rvs_alignment, aln_count list
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


def fill_in_zeros_se(fwd_rvs_align_list, ref_len):
    """
    Generate alignment counts for every nucleotide in the reference
    :param fwd_rvs_align_list:  list of sorted forwards and reverse alignments
    :param ref_len: number of nucleotides in the reference sequence (int)
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
    :param graph_processed: list of fwd and rvs upper and lower se bounds
    :param smooth_win_size: smoothing window size
    :return: list of smoothed fwd and rvs upper and lower se bound
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


def multi_header_plot(nt_list, search_terms, in_files, cutoff, plot_y_lim, win, pub, save_plot):
    """
    21,22,24nt profile plot
    :param nt_list: 
    :param search_terms: header search terms list
    :param in_files: alignment files prefix
    :param cutoff: highest count of the most abundant alignment of 21,22,24 nt profiles
    :param plot_y_lim: set y limits on plot 
    :param win: smoothing window size
    :param pub: remove box and axis labels
    """
    select_win = False
    try:
        alignment_file_list=[]
        for nt in nt_list:
            print("Loading {0} nt Alignment File\n".format(nt))
            in_file, _ = import_scram2_den(in_files+"_"+nt+".csv")
            alignment_file_list.append(in_file)
    except:
        print("\nProblem loading alignment files.  Possible a missing file for the sRNA lengths provided\n")
        sys.exit()
    substring = " ".join(search_terms)

    print("Extracting headers\n")
    all_keys=set()
    for nt in alignment_file_list:
        for header in nt.keys():
            all_keys.add(header)
    for header in all_keys:
        if substring.lower() in header.lower():
            header_alignment_tuple=()
            ref_len_tuple=()
            for alignment_file in alignment_file_list:
                alignment, ref_len = extract_header_alignment(header, alignment_file)
                header_alignment_tuple=header_alignment_tuple+(alignment,)
                ref_len_tuple=ref_len_tuple+(ref_len,)
            above_cutoff=False
            for alignment in header_alignment_tuple:
                if alignment[2]>=cutoff:
                    above_cutoff=True
            if above_cutoff:
                plot_name = save_file_name(header[1:-2])
                print (header)
                max_ref_len = max(ref_len_tuple)
                if win ==0 or select_win:
                    win = int(max_ref_len / 30)
                    select_win=True
                if win % 2 != 0 or win == 0: win += 1
                graph_processed_list=[]
                nt_pos=0
                for alignment in header_alignment_tuple:
                    graph_processed_list.append(fill_in_zeros_se(alignment, max_ref_len))
                    nt_pos+=1

                x_ref = graph_processed_list[0][0]
                smoothed_for_plot_tuple=()
                for graph_processed in graph_processed_list:
                    y_fwd_smoothed_upper, y_fwd_smoothed_lower, y_rvs_smoothed_upper, \
                    y_rvs_smoothed_lower = _smoothed_for_plot_se(graph_processed, win)
                    smoothed_for_plot_tuple=smoothed_for_plot_tuple+((y_fwd_smoothed_upper, y_fwd_smoothed_lower,
                                                                      y_rvs_smoothed_upper, y_rvs_smoothed_lower),)

                den_multi_plot(nt_list, x_ref, smoothed_for_plot_tuple, header, plot_y_lim, pub, save_plot, plot_name)


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


def den_multi_plot(nt_list, x_ref, smoothed_for_plot_tuple, header, plot_y_lim, pub, save_plot, plot_name):

    fig = plt.figure(figsize=(10, 5))
    nt_pos=0
    for smoothed_for_plot in smoothed_for_plot_tuple:
        plt.plot(x_ref, smoothed_for_plot[0], color=_nt_colour(int(nt_list[nt_pos])), label='{0} nt'.format(nt_list[
                                                                                                            nt_pos]),
                 lw=1, alpha=0.2)
        plt.plot(x_ref, smoothed_for_plot[1], color=_nt_colour(int(nt_list[nt_pos])), lw=1, alpha=0.2)
        plt.fill_between(x_ref, smoothed_for_plot[0], smoothed_for_plot[1], color=_nt_colour(int(nt_list[nt_pos])),
                         alpha=0.5)
        plt.plot(x_ref, smoothed_for_plot[2], color=_nt_colour(int(nt_list[nt_pos])), lw=1, alpha=0.2)
        plt.plot(x_ref, smoothed_for_plot[3], color=_nt_colour(int(nt_list[nt_pos])), lw=1, alpha=0.2)
        plt.fill_between(x_ref, smoothed_for_plot[2], smoothed_for_plot[3], color=_nt_colour(int(nt_list[nt_pos])),
                         alpha=0.5)
        nt_pos+=1
    axhline(y=0)
    if pub:
        _pub_plot()
    else:
        xlabel(header)
        ylabel('Reads per million reads')
        plt.legend(loc='best', fancybox=True, framealpha=0.5)
    if plot_y_lim != 0:
        ylim(-plot_y_lim, plot_y_lim)
    if save_plot:
        plt.savefig('{0}.png'.format(plot_name),dpi=300)
    plt.show()


def _pub_plot():
    """
    Remove axis, labels, legend from plot
    """
    plt.tick_params(
        axis='both',  # changes apply to the x-axis
        which='both',  # both major and minor ticks are affected
        bottom='on',  # ticks along the bottom edge are off
        top='on',
        right='on',
        left='on',  # ticks along the top edge are off
        labelbottom='off',
        labelleft='off',
        labelright='off',
        labelsize=15)  # labels along the bottom edge are off
    _clear_frame()

def save_file_name(header):
    out_file_name=""
    for i in header:
        if len(out_file_name)>100:
            break
        else:
            if i==" " or not i.isalnum():
                out_file_name+="_"
            else:
                out_file_name+=i
    return out_file_name

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



def cdp_plot_bokeh(file_name, seq1, seq2, plot_type, browser, save_plot, pub):
    output_notebook()
    if browser:
        output_file("plot.html")
    else:
        output_notebook()
    first_line=True
    x_vals_line = []
    x_vals_point = []
    # ellipse_width=[]
    xerr=[]
    max_x=0.0
    y_vals_line = []
    y_vals_point = []
    # ellipse_height=[]
    header=[]
    yerr=[]
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
                xerr.append(float(line[-3]))
                x_vals_line.append(x_se)
                y_vals_line.append([float(line[-2]),float(line[-2])])
                x_vals_line.append([float(line[-4]),float(line[-4])])
                y_vals_line.append(y_se)
                yerr.append(float(line[-1]))
                #point
                x_vals_point.append(float(line[-4]))
                y_vals_point.append(float(line[-2]))

                header.append(line[0])

    _max = max([max_x,max_y])  # sets up max x and y scale values
    log_max = _max + _max / 2
    csvfile.close()

    #Interactive
    if plot_type=="log" or plot_type=="all":
        compare_plot(header, log_max, nt, seq1, seq2, [], x_vals_point, [], y_vals_point, [], [], save_plot, pub)
    if plot_type == "log_error" or plot_type == "all":
        compare_plot(header, log_max, nt, seq1, seq2, x_vals_line, x_vals_point, y_vals_line,
                    y_vals_point, xerr, yerr, save_plot, pub)


def compare_plot(header, log_max, nt, seq1, seq2, x_vals_line, x_vals_point, y_vals_line, y_vals_point,
                xerr, yerr, save_plot, pub_plot):
    # Std Error bars
    print("Interactive log plot with se bars")
    hover = HoverTool(
        tooltips=[
            ("(x,y)", "($x, $y)"),
            ("header", "@Desc")
        ],
        names=["circle",]
    )
    p = figure(plot_width=600, plot_height=600,
               x_axis_type="log", y_axis_type="log",
               x_range=(0.1, log_max), y_range=(0.1, log_max),
               toolbar_location="above", tools=[hover, 'save', 'box_zoom', 'reset'])
    source_point = ColumnDataSource(data=OrderedDict(x=x_vals_point,
                                                     y=y_vals_point, Desc=header, )
                                    )
    p.circle('x', 'y', name="circle", source=source_point, size=3, color=_nt_colour(nt), legend="{0} nt".format(nt))
    p.legend.location = "top_left"
    p.line([0.1, log_max], [0.1, log_max])
    if xerr !=[]:
        p.multi_line(xs=x_vals_line, ys=y_vals_line, color=_nt_colour(nt), alpha=0.5, )
    p.xaxis.axis_label = seq1
    p.yaxis.axis_label = seq2
    show(p)
    if save_plot:

        fig = plt.figure(figsize=(8, 8))
        if xerr!=[]:
            plt.errorbar(x_vals_point, y_vals_point, xerr=xerr, yerr=yerr, capsize=0, ls='none', color=_nt_colour(nt),
                         elinewidth=0.5)
        plt.grid(linestyle='-',alpha=0.2)
        plt.plot([0.1, log_max], [0.1, log_max], alpha=0.3)
        plt.scatter(x_vals_point, y_vals_point, color=_nt_colour(nt), s=3, label="{0} nt".format(nt))
        plt.xlim([0.1,log_max])
        plt.ylim([0.1,log_max])
        plt.xscale('log')
        plt.yscale('log')
        if pub_plot:
            _pub_plot()
        else:
            plt.xlabel(seq1)
            plt.ylabel(seq2)
            plt.legend()

        plt.savefig('{0}_{1}_{2}.png'.format(seq1,seq2,nt), dpi=300)
        #plt.show()