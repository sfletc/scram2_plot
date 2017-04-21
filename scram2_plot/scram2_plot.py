#!/usr/bin/env python3

import plot_code as pc



from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter

__all__ = []
__version__ = 0.65
__date__ = '2016-02-08'
__updated__ = '2016-03-08'


def main(argv=None):
    """Command line options."""


    try:
        # Setup argument parser
        parser = ArgumentParser()

        subparsers = parser.add_subparsers(help="Select profile or scatter plot", dest="command")
        parser_profile = subparsers.add_parser("profile",
                                               help="Generates alignment profile/s for 1 or more reference sequences")

        #profile plot
        parser_profile.add_argument('-a', '--alignment',
                            type=str, help="sRNA alignment file prefix used by SCRAM2 profile (i.e. exclude _21.csv, _22.csv, "
                                           "_24.csv)")
        parser_profile.add_argument('-cutoff','--cutoff', type = int, default=1,
                            help = "Min. alignment RPMR from the most abundant profile (if multi) to generate plot")

        parser_profile.add_argument('-s','--search', type=str, help="Full header or substring of header", nargs='*')

        parser_profile.add_argument('-nt','--nt', type=str,help="Comma-seperated list of sRNA lengths to plot.  "
                                                                "SCRAM2 alignment files must be available for each "
                                                                "sRNA "
                                                                "length")


        parser_profile.add_argument('-ylim', '--ylim',
                            type=float, help='+/- y axis limit',
                            default=0)
        parser_profile.add_argument('-win','--win', type = int, help = 'Smoothing window size (default=auto)',
                                    default=0)

        parser_profile.add_argument('-pub', '--publish', action='store_true',
                            default=False,
                            help='Remove all labels from profiles for editing for publication')
        parser_profile.add_argument('-png', '--png', action='store_true',
                            default=False,
                            help='Export plot/s as 300 dpi .png file/s')


        #Compare plot
        parser_cdp = subparsers.add_parser("compare", help = "Generates a scatter plot for a SCRAM2 cpd alignment")
        parser_cdp.add_argument('-plot_type', '--plot_type', default="log_error",
                                    help='Bokeh plot type to display (log, log_error or all)')
        parser_cdp.add_argument('-a', '--alignment',
                            type=str, help="sRNA alignment file prefix used by SCRAM2 profile (i.e. exclude _21.csv, _22.csv, "
                                           "_24.csv)")
        parser_cdp.add_argument('-l','--length', type=str,help="Comma-seperated list of sRNA lengths to plot.  "
                                                                "SCRAM2 alignment files must be available for each "
                                                                "sRNA "
                                                                "length")
        parser_cdp.add_argument('-xlab', '--x_label', default=["Treatment 1"],
                                help='x label - corresponds to -s1 treatment in SCRAM2 arguments', nargs='*')
        parser_cdp.add_argument('-ylab', '--y_label', default=["Treatment 2"],
                                help='y label - corresponds to -s2 treatment in SCRAM2 arguments', nargs='*')
        parser_cdp.add_argument('-html', '--html', default=False, action='store_true',
                                help='If not using Jupyter Notebook, output interactive plot to browser as save to .html')
        parser_cdp.add_argument('-pub', '--publish', action='store_true',
                            default=False,
                            help='Remove all labels from profiles for editing for publication')
        parser_cdp.add_argument('-png', '--png', action='store_true',
                            default=False,
                            help='Export plot/s as 300 dpi .png file/s')


        # Process arguments
        args = parser.parse_args()

        if args.command == "profile":
            search_term = args.search
            alignment_prefix = args.alignment
            cutoff = args.cutoff
            ylim = args.ylim
            pub = args.publish
            win = args.win
            nt_list=args.length.split(',')
            save_plot=args.png
            pc.multi_header_plot(nt_list,search_term, alignment_prefix, cutoff, ylim, win, pub, save_plot)

        if args.command == "compare":
            alignment_prefix = args.alignment
            nt_list=args.nt.split(',')
            xlab= " ".join(args.x_label)
            ylab= " ".join(args.y_label)
            plot_type=args.plot_type
            browser=args.html
            save_plot=args.png
            pub = args.publish
            pc.cdp_plot_bokeh(alignment_prefix, nt_list, xlab, ylab, plot_type, browser, save_plot, pub)

    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0


if __name__ == "__main__":

    main()