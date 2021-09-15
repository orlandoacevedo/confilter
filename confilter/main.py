import argparse
import time
import sys
from confilter.__init__ import FEATURES, VERSION, FILEFORMAT
from confilter.confilter import BulkProcess, PlotSamples


def parsecmd():
    """Parse command line input"""
    def parse_remove_chars(line):
        line = line.replace('[',' ').replace(']',' ').replace(';',',')
        return line

    def parse_in_line(line,bobcon=False,boacon=False):
        line = parse_remove_chars(line)
        line = line.replace('-',' ').replace('_',' ')
        lt = line.split(',')
        ls = []
        for t in lt:
            lp = t.split()
            if not len(lp): continue
            if bobcon:
                if len(lp) == 2:
                    ls.append(lp)
                else:
                    print('Warning: wrong bond connection < {:} >'.format(line))
                    raise ValueError('wrongly defined')
            elif boacon:
                if len(lp) == 3:
                    ls.append(lp)
                else:
                    print('Warning: wrong angle connection < {:} >'.format(line))
                    raise ValueError('wrongly defined')
            else:
                ls.append(lp)

        if not len(ls): return []

        lm = []
        for t in ls:
            lx = []
            for v in t:
                try:
                    lx.append(int(v))
                except ValueError:
                    print('Warning: wrong input < {:} >'.format(line))
                    raise ValueError('wrongly defined')
            lm.append(lx)

        return lm


    parser = argparse.ArgumentParser(
        description='Conformation Filtration',
        allow_abbrev=False,
    )
    parser.add_argument(
        '-v','--version',
        action='version',
        version=VERSION,
    )
    parser.add_argument(
        '-f','--datafilelist',
        help='Data files, separate by space or comma',
        nargs='+',
        metavar='file',
    )
    parser.add_argument(
        '-d','--indexfilelist',
        help='Index files, as the filtration reference',
        nargs='+',
        metavar='file',
    )
    parser.add_argument(
        '--static',
        help='turn on static mode calculation, default is in dynamic/all mode',
        action='store_true',
    )
    parser.add_argument(
        '--separate',
        help='valid in dynamic mode, change to dynamic/separate mode',
        action='store_true',
    )
    parser.add_argument(
        '--vndx',
        help='valid in static mode, set filtration index value',
        type=float,
    )
    parser.add_argument(
        '--borandom',
        help='valid in static mode, rather than by lowest-bit, filtering out molecules randomly',
        action='store_true',
    )
    parser.add_argument(
        '--seed',
        help='random seed, (optional, highest priority)',
        type=int,
    )
    parser.add_argument(
        '-btol','--btol',
        help='bonds tolerance, any changes in bcon smaller than it will be filtered out',
        type=float,
    )
    parser.add_argument(
        '-atol','--atol',
        help='angles tolerance, any changes in acon smaller than it will be filtered out',
        type=float,
    )
    parser.add_argument(
        '-bcon',
        help='Bond connections, in pairs, separate by comma',
        nargs='+',
        metavar='B1 B2, B1-B2',
    )
    parser.add_argument(
        '-acon',
        help='Angle connections, in pairs, separate by comma',
        nargs='+',
        metavar='A1 A2 A3, A1-A2-A3',
    )
    parser.add_argument(
        '-g','--fragments',
        help='fragments, separate by comma',
        nargs='+',
        metavar='G',
    )
    parser.add_argument(
        '--no-oball',
        help='turn off bonds probability overall calculation, Boolean',
        action='store_true',
    )
    parser.add_argument(
        '-obpar','--obpar',
        help='turn on bonds probability parameters calculation, Boolean',
        action='store_true',
    )
    parser.add_argument(
        '--no-oaall',
        help='turn off angles probability overall calculation, Boolean',
        action='store_true',
    )
    parser.add_argument(
        '-oapar','--oapar',
        help='turn on angles probability parameters calculation, Boolean',
        action='store_true',
    )
    parser.add_argument(
        '-nc','--no-force-double-check',
        help='turn off double check prompt info before execution',
        action='store_true',
    )
    parser.add_argument(
        '--features',
        help='show development features',
        action='store_true',
    )
    parser.add_argument(
        '-p','--file-format-explanations',
        help='show input system file format explanations',
        action='store_true',
    )
    parser.add_argument(
        '-nu','--no-userinputs',
        help='specify the indexes of input connections start at 0',
        action='store_true',
    )
    parser.add_argument(
        '-o','--fname',
        help='output system file name',
    )
    parser.add_argument(
        '-ft','--ftype',
        help='output system file type, [txt, xsf, xyz]',
    )
    subparser = parser.add_subparsers(title='continuous subcommand')
    sub = subparser.add_parser(
        'plot',
        help='plot cross comparsions based on different chosen samples',
        allow_abbrev=False
    )
    sub.set_defaults(command='plot')
    sub.add_argument(
        '-bf','--probdatafilelist',
        help='bulk process probability files',
        metavar='B',
        nargs='+',
    )
    sub.add_argument(
        '-nl','--nmlist',
        help='highest priority, select number of samples for plot, (optional)',
        metavar='n',
        nargs='+',
        type=int,
    )
    sub.add_argument(
        '-ns','--nmsamples',
        help='choose number of samples for plot, default is 5',
        metavar='n',
        type=int,
    )
    sub.add_argument(
        '-sn','--startndx',
        help='start index for the inputs, (optional)',
        metavar='n',
        type=int,
    )
    sub.add_argument(
        '-en','--endndx',
        help='end index for the inputs, (optional)',
        metavar='n',
        type=int,
    )
    sub.add_argument(
        '-inc','--incndx',
        help='increments for choose, (optional)',
        metavar='n',
        type=int,
    )
    sub.add_argument(
        '-nr','--nmranges',
        help='random ranges for increments, (optional)',
        metavar='n',
        type=int,
    )

    # annoying part
    bo = False
    if len(sys.argv) == 1: bo = True
    if not bo and 'plot' not in sys.argv:
        if '-h' in sys.argv or '--help' in sys.argv: bo = True
    if bo:
        parser.print_help()
        exit()

    # due to the potential bug, args may not have enough namespace
    args,left = sub.parse_known_args()
    if 'plot' in left:
        left.remove('plot')
        opts = parser.parse_args(left)
        # conclude them
        args = argparse.Namespace(**vars(opts),**vars(args))
    else:
        args = parser.parse_args(sys.argv[1:])

    if 'features' in args and args.features:
        for i in FEATURES:
            print(i)
        exit()

    if 'file_format_explanations' in args and args.file_format_explanations:
        print(FILEFORMAT)
        exit()

    # default settings
    fdict = {
        'datafilelist'              :   None,
        'indexfilelist'             :   None,
        'bcon'                      :   None,
        'btol'                      :   None,
        'oball'                     :   True,
        'obpar'                     :   False,
        'acon'                      :   None,
        'atol'                      :   None,
        'oaall'                     :   True,
        'oapar'                     :   False,
        'mode'                      :   None,
        'vndx'                      :   None,
        'borandom'                  :   None,
        'boall'                     :   True,
        'fragments'                 :   None,
        'bool_force_double_check'   :   True,
        'userinputs'                :   True,
        'fname'                     :   None,
        'ftype'                     :   None,
        'probdatafilelist'          :   None,
        'nmlist'                    :   None,
        'nmsamples'                 :   None,
        'startndx'                  :   None,
        'endndx'                    :   None,
        'incndx'                    :   None,
        'nmranges'                  :   None,
        'seed'                      :   None,
    }

    bod = False
    if 'datafilelist' in args and args.datafilelist is not None:
        if len(args.datafilelist): bod = True
    bog = False
    if 'probdatafilelist' in args and args.probdatafilelist is not None:
        if len(args.probdatafilelist): bog = True
    if 'command' in args:
        if not bod and not bog:
            print('Fatal: no input: missing:  -f/--datafilelist or -bf/--probdatafilelist')
            exit()
    else:
        if not bod:
            print('Warning: -f/--datafilelist is missing')
            exit()
    
    if bod:
        stmp = parse_remove_chars(' '.join(args.datafilelist)).replace(',',' ')
        fdict['datafilelist'] = stmp.split()
    
    if bog:
        stmp = parse_remove_chars(' '.join(args.probdatafilelist)).replace(',',' ')
        fdict['probdatafilelist'] = stmp.split()

    if 'indexfilelist' in args and args.indexfilelist:
        stmp = parse_remove_chars(' '.join(args.indexfilelist)).replace(',',' ')
        fdict['indexfilelist'] = stmp.split()

    if 'static' in args and args.static: fdict['mode'] = 'static'
    if 'separate' in args and args.separate: fdict['boall'] = False
    if 'vndx' in args: fdict['vndx'] = args.vndx     # be aware of 0.0
    if 'borandom' in args and args.borandom: fdict['borandom'] = True

    if 'bcon' in args and args.bcon:
        fdict['bcon'] = parse_in_line(' '.join(args.bcon),bobcon=True)
    if 'acon' in args and args.acon:
        fdict['acon'] = parse_in_line(' '.join(args.acon),boacon=True)
    if 'fragments' in args and args.fragments:
        fdict['fragments'] = parse_in_line(' '.join(args.fragments))

    if 'btol' in args and args.btol: fdict['btol'] = args.btol
    if 'atol' in args and args.atol: fdict['atol'] = args.atol
    if 'oball' in args and args.no_oball: fdict['oball'] = False
    if 'obpar' in args and args.obpar: fdict['obpar'] = True
    if 'oaall' in args and args.no_oaall: fdict['oaall'] = False
    if 'oapar' in args and args.oapar: fdict['oapar'] = True
    if 'no_force_double_check' in args and args.no_force_double_check:
        fdict['bool_force_double_check'] = False
    if 'no_userinputs' in args and args.no_userinputs: fdict['userinputs'] = False
    if 'fname' in args and args.fname: fdict['fname'] = args.fname
    if 'ftype' in args and args.ftype: fdict['ftype'] = args.ftype
    if 'nmlist' in args and args.nmlist: fdict['nmlist'] = args.nmlist
    if 'nmsamples' in args and args.nmsamples: fdict['nmsamples'] = args.nmsamples
    if 'startndx' in args and args.startndx: fdict['startndx'] = args.startndx
    if 'endndx' in args and args.endndx: fdict['endndx'] = args.endndx
    if 'incndx' in args and args.incndx: fdict['incndx'] = args.incndx
    if 'nmranges' in args and args.nmranges: fdict['nmranges'] = args.nmranges
    if 'seed' in args and args.seed: fdict['seed'] = args.seed

    print('Note: time: {:}'.format(time.ctime()))
    if 'command' in args:
        print('Note: processing plot ...')
        PS = PlotSamples(**fdict)
    else:
        print('Note: processing data files ...')
        PS = BulkProcess(**fdict)
    if not PS.nice:
        print(PS.info)
        return
    PS.run()
    if not PS.nice: print(PS.info)


if __name__ == '__main__':
    parsecmd()


