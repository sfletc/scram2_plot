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
        profile_subparsers=parser_profile.add_subparsers(help="Select single (x nt plot) or multi (21,22 and 24 nt "
                                                              "plot)", dest="subcommand")

        profile_single = profile_subparsers.add_parser("single", help = "single (x nt) profile plot")
        profile_multi = profile_subparsers.add_parser("multi", help="multi (21, 22, 24 nt) profile plot")

        #parser_profile.add_argument('-plot', '--plot', type=str, help="single or multi")

        #x nt plot
        profile_single.add_argument('-a', '--alignment',
                            type=str, help="siRNA alignment file generated by SCRAM2 profile")
        profile_single.add_argument('-cutoff','--cutoff', type = int,
                            help = "Min. alignment RPMR from the most abundant profile (if multi) to generate plot")

        profile_single.add_argument('-s','--search', type=str, help="Full header or substring of header", nargs='*')

        profile_single.add_argument('-ylim', '--ylim',
                            type=float, help='+/- y axis limit',
                            default=0)
        profile_single.add_argument('-win','--win', type = int, help = 'Smoothing window size (default=auto)',
                                    default=0)

        profile_single.add_argument('-pub', '--publish', action='store_true',
                            default=False,
                            help='Remove all labels from profiles for editing for publication')


        #21,22,24nt plot
        profile_multi.add_argument('-a', '--alignments',
                            type=str, help="21nt, 22nt and 24nt profile alignment csvs must be included and must have the sRNA length at the terminus of the file name (e.g. _21.csv", nargs='*')
        profile_multi.add_argument('-cutoff','--cutoff', type = int,
                            help = "Min. alignment RPMR from the most abundant profile (if multi) to generate plot")

        profile_multi.add_argument('-s','--search', type=str, help="Full header or substring of header", nargs='*')

        profile_multi.add_argument('-ylim', '--ylim',
                            type=float, help='+/- y axis limit',
                            default=0)
        profile_multi.add_argument('-win','--win', type = int, help = 'Smoothing window size (default=auto)',
                                    default=0)
        profile_multi.add_argument('-pub', '--publish', action='store_true',
                            default=False,
                            help='Remove all labels from (non-bokeh) profiles for \
                            publication')

        #Generate a CDP scatter plot
        # cdp_plot_bokeh(file_name, seq1, seq2, nt)
        parser_cdp = subparsers.add_parser("scatter", help = "Generates a scatter plot for a SCRAM2 cpd alignment")
        parser_cdp.add_argument('-plot_type', '--plot_type', default="log_error",
                                    help='Bokeh plot type to display (log, log_error, linear or all)')
        parser_cdp.add_argument('-a','--alignment', help="SCRAM2 cdp alignment file")
        parser_cdp.add_argument('-xlab', '--x_label', default=["Treatment 1"],
                                help='x label - corresponds to -s1 treatment in SCRAM2 arguments', nargs='*')
        parser_cdp.add_argument('-ylab', '--y_label', default=["Treatment 2"],
                                help='y label - corresponds to -s2 treatment in SCRAM2 arguments', nargs='*')
        parser_cdp.add_argument('-browser', '--browser', default=False, action='store_true',
                                help='If not using Jupyter Notebook, output plot to browser')

        # Process arguments
        args = parser.parse_args()

        if args.command == "profile":
            if args.subcommand == "multi":
                search_term = args.search
                a = args.alignments
                cutoff = args.cutoff
                ylim = args.ylim
                pub = args.publish
                win = args.win
                pc.multi_header_plot(search_term, a, cutoff, ylim, win, pub)

            elif args.subcommand =="single":
                search_term = args.search
                a = args.alignment
                cutoff = args.cutoff
                ylim = args.ylim
                pub = args.publish
                win = args.win
                pc.single_header_plot(search_term, a, cutoff, ylim, win, pub)

        if args.command == "scatter":
            alignment_file = args.alignment
            xlab= " ".join(args.x_label)
            ylab= " ".join(args.y_label)
            plot_type=args.plot_type
            browser=args.browser
            pc.cdp_plot_bokeh(alignment_file, xlab, ylab, plot_type, browser)

    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0


if __name__ == "__main__":

    main()