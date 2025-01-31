#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
ecoli_stx_subtyping.py - Gets E. coli stx subtypes
<https://github.com/B-UMMI/seq_typing/modules/>

Copyright (C) 2019 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: January 10, 2019

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
import argparse

try:
    from __init__ import __version__

    import modules.utils as utils
    import modules.parse_results as parse_results

    from seq_typing import python_arguments as python_arguments
except ImportError:
    from seqtyping.__init__ import __version__

    from seqtyping.modules import utils as utils
    from seqtyping.modules import parse_results as parse_results

    from seqtyping.seq_typing import python_arguments


def stx_subtype_parser(
        report_types,
        stx1_reference_file, stx2_reference_file,
        stx2_alternative_sequenced_covered, stx2_alternative_sequence_identity
    ):
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
    stx1_result = "NT"
    stx2_result_main = None
    stx2_result_other = []
    with open(report_types, "rt") as reader:
        for line in reader:
            if not line.startswith("#"):
                line = line.split("\t")
                subtype = line[2]
                if stx1_reference_file is not None and line[1] == stx1_reference_file:
                    if line[0] == "selected":
                        stx1_result = subtype
                if stx2_reference_file is not None and line[1] == stx2_reference_file:
                    if line[0] == "selected":
                        stx2_result_main = subtype
                    elif line[0] == "other_probable_type":
                        if float(line[4]) >= stx2_alternative_sequenced_covered and \
                                float(line[6]) >= stx2_alternative_sequence_identity and \
                                subtype != stx2_result_main and \
                                subtype not in stx2_result_other:
                            stx2_result_other.append(subtype)

        if stx2_result_main is not None and len(stx2_result_other) > 0:
            stx2_result_main = "{main}({other})".format(main=stx2_result_main,
                                                        other=";".join(sorted(stx2_result_other)))
        elif stx2_result_main is None:
            stx2_result_main = "NT"

    return stx1_result, stx2_result_main


def get_ref_file(refereces_list, locus_name):
    ref_file = None
    for ref_file_ in refereces_list:
        if locus_name.lower() in os.path.basename(ref_file_).lower():
            ref_file = ref_file_
            break
    if ref_file is None:
        print(f"No reference file found for {locus_name}")
    else:
        print(f"The following reference file will be used for {locus_name} : {ref_file}")
    return ref_file


def main():
    program_name = "ecoli_stx_subtyping.py"

    if sys.version_info[0] < 3:
        sys.exit("Must be using Python 3. Try calling \"python3 {}\"".format(program_name))

    parser, parser_reads, _, parser_assembly, _ = python_arguments(program_name=program_name, version=__version__)
    parser.description = "Gets E. coli stx subtypes"

    # Add specific arguments
    parser_reads.add_argument("--stx2Acovered", type=float,
                              metavar="N",
                              help="Minimal percentage of sequence covered to consider extra stx2A"
                                   " subtypes (value between [0, 100])",
                              required=False, default=100)
    parser_reads.add_argument("--stx2Aidentity", type=float,
                              metavar="N",
                              help="Minimal sequence identity to consider extra stx2A"
                                   " subtypes (value between [0, 100])",
                              required=False, default=99.5)

    parser_assembly.add_argument("--stx2Acovered", type=float,
                                 metavar="N",
                                 help="Minimal percentage of sequence covered to consider extra stx2A"
                                      " subtypes (value between [0, 100])",
                                 required=False, default=100)
    parser_assembly.add_argument("--stx2Aidentity", type=float,
                                 metavar="N",
                                 help="Minimal sequence identity to consider extra stx2A"
                                      " subtypes (value between [0, 100])",
                                 required=False, default=99.5)

    parser_reads.add_argument("--stx2Bcovered", type=float,
                              metavar="N",
                              help="Minimal percentage of sequence covered to consider extra stx2B"
                                   " subtypes (value between [0, 100])",
                              required=False, default=95)
    parser_reads.add_argument("--stx2Bidentity", type=float,
                              metavar="N",
                              help="Minimal sequence identity to consider extra stx2B"
                                   " subtypes (value between [0, 100])",
                              required=False, default=95)

    parser_assembly.add_argument("--stx2Bcovered", type=float,
                                 metavar="N",
                                 help="Minimal percentage of sequence covered to consider extra stx2B"
                                      " subtypes (value between [0, 100])",
                                 required=False, default=95)
    parser_assembly.add_argument("--stx2Bidentity", type=float,
                                 metavar="N",
                                 help="Minimal sequence identity to consider extra stx2B"
                                      " subtypes (value between [0, 100])",
                                 required=False, default=95)

    args = parser.parse_args()

    msg = []
    if args.minGeneCoverage < 0 or args.minGeneCoverage > 100:
        msg.append("--minGeneCoverage should be a value between [0, 100]")
    # if args.minGeneIdentity < 0 or args.minGeneIdentity > 100:
    #     msg.append("--minGeneIdentity should be a value between [0, 100]")
    if args.stx2Acovered < 0 or args.stx2Acovered > 100:
        msg.append("--stx2Acovered should be a value between [0, 100]")
    if args.stx2Aidentity < 0 or args.stx2Aidentity > 100:
        msg.append("--stx2Aidentity should be a value between [0, 100]")
    if args.stx2Bcovered < 0 or args.stx2Bcovered > 100:
        msg.append("--stx2Bcovered should be a value between [0, 100]")
    if args.stx2Bidentity < 0 or args.stx2Bidentity > 100:
        msg.append("--stx2Bidentity should be a value between [0, 100]")
    if args.blast is None and args.org != ["stx", "subtyping"]:
        msg.append("Use \"--org stx subtyping\" with {}".format(program_name))

    if len(msg) > 0:
        argparse.ArgumentParser(prog="{} options".format(program_name)).error("\n".join(msg))

    start_time = time.time()

    args.outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    # Start logger
    logfile, _ = utils.start_logger(args.outdir)

    _ = utils.general_information(script_name=program_name, logfile=logfile, version=__version__)
    print("\n")

    folders_2_remove = []

    # Create modules pickles folder
    pickles_folder = os.path.join(args.outdir, "pickles", "")
    if not os.path.isdir(pickles_folder):
        os.makedirs(pickles_folder)
    folders_2_remove.append(pickles_folder)

    # Run functions
    folders_2_remove_func, references_results, reference, references_headers, assembly, args_config_dict = \
        args.func(args)
    folders_2_remove.extend(folders_2_remove_func)

    min_gene_identity = args_config_dict["minimum_gene_identity"]
    min_gene_coverage = args_config_dict["minimum_gene_coverage"]
    type_separator = args_config_dict["type_separator"]
    length_extra_seq = args_config_dict["length_extra_seq"]

    type_in_new = not args.typeNotInNew

    # Parse results
    _, _, _, _, _ = parse_results.parse_results(references_results, reference, references_headers, args.outdir,
                                                min_gene_coverage, args.minDepthCoverage, type_separator,
                                                sample=args.sample, save_new_allele=args.saveNewAllele,
                                                assembly=assembly, extra_seq=length_extra_seq,
                                                min_identity=min_gene_identity,
                                                type_in_new=type_in_new)

    stx1A_result, stx2A_result = stx_subtype_parser(
        os.path.join(args.outdir, "seq_typing.report_types.tab"),
        get_ref_file(reference, "stx1A"), get_ref_file(reference, "stx2A"),
        args.stx2Acovered, args.stx2Aidentity
    )

    stx1B_result, stx2B_result = stx_subtype_parser(
        os.path.join(args.outdir, "seq_typing.report_types.tab"),
        get_ref_file(reference, "stx1B"), get_ref_file(reference, "stx2B"),
        args.stx2Bcovered, args.stx2Bidentity
    )

    # Rename the file to keep ecoli_stx_subtyping stamp
    if os.path.isfile(os.path.join(args.outdir, "seq_typing.report_types.tab")):
        os.rename(os.path.join(args.outdir, "seq_typing.report_types.tab"),
                  os.path.join(args.outdir, "seq_typing.ecoli_stx_subtyping.report_types.tab"))

    # Remove the file to only keep the ecoli_stx_subtyping one
    if os.path.isfile(os.path.join(args.outdir, "seq_typing.report.txt")):
        os.remove(os.path.join(args.outdir, "seq_typing.report.txt"))

    stx_subtype_str = ":".join([stx1A_result, stx1B_result, stx2A_result, stx2B_result])
    print(
        "\n"
        f"E. coli stx_subtyping - {stx_subtype_str}\n"
        "\n"
    )
    with open(os.path.join(args.outdir, "seq_typing.ecoli_stx_subtyping.txt"), "wt", encoding="utf-8") as writer:
        print(stx_subtype_str, file=writer)

    if not args.debug:
        for folder in folders_2_remove:
            utils.removeDirectory(folder)

    _ = utils.run_time(program_name, start_time)


if __name__ == "__main__":
    main()
