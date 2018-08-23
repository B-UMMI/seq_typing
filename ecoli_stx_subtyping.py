#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
ecoli_stx_subtyping.py - Gets E. coli subtypes
<https://github.com/B-UMMI/seq_typing/modules/>

Copyright (C) 2018 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: August 22, 2018

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import os.path
import time

from modules import utils
from modules import parse_results


version = '1.0'


def stx_subtype_parser(report_types, stx1_reference_file, stx2_reference_file, stx2_alternative_sequenced_covered,
                       stx2_alternative_sequence_identity):
    """
    Parse the stx subtypes, specially regarding stx2 that might have multiple sequences

    Parameters
    ----------
    report_types : str
        Path to seq_typing.report_types.tab file
    stx1_reference_file : str
        Path to stx1 reference sequences file
    stx2_reference_file : str
        Path to stx2 reference sequences file
    stx2_alternative_sequenced_covered : float
        Minimal percentage of sequence covered to consider extra stx2 subtypes
    stx2_alternative_sequence_identity : float
        Minimal sequence identity to consider extra stx2 subtypes

    Returns
    -------
    stx1_result : str
        String containing stx1 subtype, NT if seq_typing ran but no subtype could be determined, or NA if it didn't run
    stx2_result_main : str
        String containing stx2 subtype, NT if seq_typing ran but no subtype could be determined, or NA if it didn't run.
        If extra subtypes were found, they will be present inside brackets separated by ;
    """

    stx1_result = 'NA'
    stx2_result_main = None
    stx2_result_other = []
    with open(report_types, 'rt') as reader:
        stx1_run = False
        stx2_run = False
        for line in reader:
            if not line.startswith('#'):
                line = line.split('\t')
                subtype = line[2]
                if os.path.basename(line[1]) == stx1_reference_file:
                    if line[0] == 'selected':
                        stx1_result = subtype
                    stx1_run = True
                if os.path.basename(line[1]) == stx2_reference_file:
                    if line[0] == 'selected':
                        stx2_result_main = subtype
                    elif line[0] == 'other_probable_type':
                        if float(line[4]) >= stx2_alternative_sequenced_covered and \
                                float(line[6]) >= stx2_alternative_sequence_identity and \
                                subtype != stx2_result_main and \
                                subtype not in stx2_result_other:
                            stx2_result_other.append(subtype)
                    stx2_run = True
        if stx1_run and stx1_result == 'NA':
            stx1_result = 'NT'
        if stx2_result_main is not None and len(stx2_result_other) > 0:
            stx2_result_main = '{main}({other})'.format(main=stx2_result_main,
                                                        other=';'.join(sorted(stx2_result_other)))
        elif stx2_result_main is None:
            if stx2_run:
                stx2_result_main = 'NT'
            else:
                stx2_result_main = 'NA'

    return stx1_result, stx2_result_main


def main():
    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling "python3 seq_typing.py"')

    sys.path.append('..')
    from seq_typing import python_arguments

    parser, parser_reads, parser_assembly, parser_blast = python_arguments()
    parser.prog = 'ecoli_stx_subtyping.py'
    parser.description = 'Gets E. coli subtypes'

    # Add specific arguments
    parser_reads.add_argument('--stx2covered', type=float,
                              metavar='95',
                              help='Minimal percentage of sequence covered to consider extra stx2'
                                   ' subtypes',
                              required=False, default=100)
    parser_reads.add_argument('--stx2identity', type=float,
                              metavar='95',
                              help='Minimal sequence identity to consider extra stx2 subtypes',
                              required=False, default=100)

    parser_assembly.add_argument('--stx2covered', type=float,
                                 metavar='95',
                                 help='Minimal percentage of sequence covered to consider extra stx2'
                                      ' subtypes',
                                 required=False, default=100)
    parser_assembly.add_argument('--stx2identity', type=float,
                                 metavar='95',
                                 help='Minimal sequence identity to consider extra stx2 subtypes',
                                 required=False, default=100)

    args = parser.parse_args()

    start_time = time.time()

    args.outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    # Start logger
    logfile, time_str = utils.start_logger(args.outdir)

    script_path = utils.general_information(logfile, version, args.outdir, time_str)
    del script_path
    print('\n')

    folders_2_remove = []

    # Create modules pickles folder
    pickles_folder = os.path.join(args.outdir, 'pickles', '')
    if not os.path.isdir(pickles_folder):
        os.makedirs(pickles_folder)
    folders_2_remove.append(pickles_folder)

    # Run functions
    folders_2_remove_func, references_results, reference, references_headers = args.func(args)
    folders_2_remove.extend(folders_2_remove_func)

    # Parse results
    _, _, _, _, _ = parse_results.parse_results(references_results, reference, references_headers, args.outdir,
                                                args.minGeneCoverage, args.minDepthCoverage, args.typeSeparator)

    stx1_result, stx2_result = stx_subtype_parser(
        os.path.join(args.outdir, 'seq_typing.report_types.tab'),
        [ref_file for ref_file in reference if 'stx1' in os.path.basename(ref_file).lower()][0],
        [ref_file for ref_file in reference if 'stx2' in os.path.basename(ref_file).lower()][0],
        args.stx2covered, args.stx2identity)

    print('\n'
          'E. coli stx_subtyping - {stx1_result}:{stx2_result}\n'
          '\n'.format(stx1_result=stx1_result, stx2_result=stx2_result))
    with open(os.path.join(args.outdir, 'seq_typing.ecoli_stx_subtyping.txt'), 'wt') as writer:
        writer.write(':'.join([stx1_result, stx2_result]))

    if not args.debug:
        for folder in folders_2_remove:
            utils.removeDirectory(folder)

    _ = utils.runTime(start_time)


if __name__ == "__main__":
    main()
