#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
seq_typing.py - Determine which reference sequence is more likely to be present
in a given sample
<https://github.com/B-UMMI/seq_typing/>

Copyright (C) 2018 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: October 15, 2018

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

import argparse
import os
import time
import sys

import modules.utils as utils
import modules.run_rematch as run_rematch
import modules.parse_results as parse_results
import modules.run_blast as run_blast

version = '2.1'


def parse_config(config_file):
    config = {'length_extra_seq': None, 'minimum_depth_presence': None, 'minimum_depth_call': None,
              'minimum_gene_coverage': None}

    with open(config_file, 'rtU') as reader:
        field = None
        for line in reader:
            line = line.rstrip('\r\n')
            if len(line) > 0:
                line = line.split(' ')[0]
                if line.startswith('#'):
                    line = line[1:].split(' ')[0]
                    field = line
                else:
                    if field is not None:
                        if field in ['length_extra_seq', 'minimum_depth_presence', 'minimum_depth_call',
                                     'minimum_gene_coverage']:
                            line = int(line)
                            if field == 'minimum_gene_coverage':
                                if line < 0 or line > 100:
                                    sys.exit('minimum_gene_coverage in config file must be an integer between 0 and'
                                             ' 100')
                        config[field] = line
                        field = None

    for field in config:
        if config[field] is None:
            sys.exit('{} in config file is missing'.format(field))

    return config


def get_fasta_config(species):
    """
    Get the reference fasta file and config file for the species provided

    Parameters
    ----------
    species : list
        List with strings that correspond to the species name, e.g. ('escherichia', 'coli')

    Returns
    -------
    fasta : list
        Sorted list of fasta files to be used as references
    config : str
        File path to config file
    """

    fasta = []
    config = None

    file_path = os.path.abspath(__file__)
    species_folder = os.path.join(os.path.dirname(file_path), 'reference_sequences', '_'.join(species), '')

    files = [f for f in os.listdir(species_folder) if not f.startswith('.') and
             os.path.isfile(os.path.join(species_folder, f))]
    for file_found in files:
        file_found = os.path.join(species_folder, file_found)
        if file_found.endswith('.fasta'):
            fasta.append(file_found)
        elif file_found.endswith('.config'):
            config = str(file_found)

    return sorted(fasta), config


def get_species_allowed():
    """
    Get the name of the species with reference sequences provided for serotyping

    Returns
    -------
    species : list
        List with species names, e.g. ['escherichia coli', 'streptococcus agalactiae']
    """

    file_path = os.path.abspath(__file__)
    serotyping_folder = os.path.join(os.path.dirname(file_path), 'reference_sequences', '')
    species = [d.replace('_', ' ') for d in os.listdir(serotyping_folder) if not d.startswith('.') and
               os.path.isdir(os.path.join(serotyping_folder, d))]
    return species


def include_rematch_dependencies_path():
    command = ['which', 'rematch.py']
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
    if run_successfully:
        rematch_script = stdout.splitlines()[0]
        utils.setPATHvariable(False, rematch_script)
        return rematch_script
    else:
        sys.exit('ReMatCh not found in the PATH')


def clean_header(header, problematic_characters):
    new_header = header
    if any(x in header for x in problematic_characters):
        for x in problematic_characters:
            new_header = new_header.replace(x, '_')
    return header, new_header


def parse_reference(reference, problematic_characters):
    reference_dict = {}
    headers_correspondence = {}
    with open(reference, 'rtU') as reader:
        header = None
        sequence = ''
        for line in reader:
            line = line.rstrip('\r\n')
            if len(line) > 0:
                if line.startswith('>'):
                    if header is not None:
                        reference_dict[header] = sequence
                    original_header, new_header = clean_header(line[1:], problematic_characters)
                    if new_header in headers_correspondence:
                        sys.exit('Possible conflicting sequence header in {reference} file:\n'
                                 '{original_header} header might be the same as {first_header} header after problematic'
                                 ' characters ({problematic_characters}) replacement (new header: {new_header})'
                                 ''.format(reference=reference, original_header=original_header,
                                           first_header=headers_correspondence[new_header],
                                           problematic_characters=problematic_characters, new_header=new_header))
                    header = str(new_header)
                    headers_correspondence[header] = str(original_header)
                    sequence = ''
                else:
                    sequence += line.replace(' ', '').upper()
        if len(sequence) > 0:
            reference_dict[header] = sequence
    return reference_dict, headers_correspondence


def rename_duplicated_headers(references_headers, reference, reference_dict, headers_correspondence,
                              problematic_characters):
    renamed_reference_dict, renamed_headers_correspondence, headers_changed = {}, {}, []
    for ref in references_headers:
        if any(x in references_headers[ref].keys() for x in reference_dict):
            for header in reference_dict:
                if header in references_headers[ref]:
                    original_header, new_header = clean_header('_'.join([header,
                                                                         os.path.basename(reference),
                                                                         'SeqTyping']),
                                                               problematic_characters)
                    renamed_reference_dict[new_header] = reference_dict[header]
                    renamed_headers_correspondence[new_header] = headers_correspondence[header]
                    headers_changed.append(header)
                    utils.Bcolors_print('{header} sequence header from {reference} file was found already'
                                        ' in {ref} file. Therefore it was renamed'
                                        ' to {new_header}'.format(header=header, reference=os.path.basename(reference),
                                                                  ref=os.path.basename(ref),
                                                                  new_header=new_header), 'WARNING')
    for header in reference_dict:
        if header not in headers_changed:
            renamed_reference_dict[header] = reference_dict[header]
            renamed_headers_correspondence[header] = headers_correspondence[header]
    return renamed_reference_dict, renamed_headers_correspondence


def write_sequence(reference_dict, writer):
    for header, sequence in reference_dict.items():
        writer.write('>' + header + '\n')
        writer.write('\n'.join(utils.chunkstring(sequence, 80)) + '\n')


def prepare_references(references, map_ref_together, references_dir):
    problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/", ":"]
    references_headers = {}
    references_files = []

    if map_ref_together:
        ref_file = os.path.join(references_dir, 'references_together.fasta')
        writer = open(ref_file, 'wt')
        references_files.append(ref_file)

    for reference in references:
        reference_dict, headers_correspondence = parse_reference(reference, problematic_characters)
        reference_dict, headers_correspondence = rename_duplicated_headers(references_headers, reference,
                                                                           reference_dict, headers_correspondence,
                                                                           problematic_characters)

        references_headers[reference] = dict(headers_correspondence)

        if not map_ref_together:
            ref_file = os.path.join(references_dir, os.path.basename(reference))
            writer = open(ref_file, 'wt')
            references_files.append(ref_file)

        write_sequence(reference_dict, writer)
        writer.flush()

        if not map_ref_together:
            writer.close()

    if map_ref_together:
        writer.close()

    return references_files, references_headers


def assembly_subcommand(args):
    msg = []
    # if args.blast is None and args.org is None:
    #     msg.append('--blast or --org must be provided')
    if args.blast is not None and args.type is None:
        msg.append('With --blast option you must provide the --type')

    if len(msg) > 0:
        argparse.ArgumentParser.error('\n'.join(msg))

    if args.type == 'nucl':
        utils.required_programs({'blastn': ['-version', '>=', '2.6.0']})
    elif args.type == 'prot':
        utils.required_programs({'blastp': ['-version', '>=', '2.6.0']})

    args.fasta = os.path.abspath(args.fasta[0].name)
    if args.blast is not None:
        args.blast = [blast.name for blast in args.blast]
    else:
        args.blast, config = get_fasta_config(args.org)
        config = parse_config(config)
        if args.type != 'nucl':
            print('\n'
                  'ATTENTION: Blast DB type provided was not "nucl"'
                  'It was changed to "nucl"'
                  '\n')
        args.type = 'nucl'
        args.minGeneCoverage = config['minimum_gene_coverage']
        args.typeSeparator = '_'

        print('\n'
              'Settings that will be used:\n'
              '    DB reference file: {reference}\n'
              '    Blast DB type: nucl\n'
              '    minGeneCoverage: {minGeneCoverage}\n'
              '    Type separator character: {typeSeparator}'
              '\n'.format(reference=args.blast, minGeneCoverage=args.minGeneCoverage, typeSeparator=args.typeSeparator))

    folders_2_remove_all = []
    references_results_all = {}
    references_headers_all = {}
    blast_files = []
    for blast in args.blast:
        _, folders_2_remove, blast_results, blast, headers_correspondence = run_blast.run_blast(blast, args.outdir,
                                                                                                args.type, args.fasta)
        folders_2_remove_all.extend(folders_2_remove)
        references_results_all[blast] = blast_results
        references_headers_all[blast] = headers_correspondence
        blast_files.append(blast)

    return folders_2_remove_all, references_results_all, blast_files, references_headers_all


def blast_subcommand(args):
    msg = []
    if args.fasta is not None and args.type is None:
        msg.append('With --fasta option you must provide the --type')
    # if args.fasta is None and args.org is None:
    #     msg.append('--fasta or --org must be provided')

    if len(msg) > 0:
        argparse.ArgumentParser.error('\n'.join(msg))

    utils.required_programs({'makeblastdb': ['-version', '>=', '2.6.0']})

    if args.fasta is not None:
        args.fasta = [os.path.abspath(fasta.name) for fasta in args.fasta]
    else:
        args.fasta, _ = get_fasta_config(args.org)
        if args.type != 'nucl':
            print('\n'
                  'ATTENTION: Blast DB type provided was not "nucl"'
                  'It was changed to "nucl"'
                  '\n')
        args.type = 'nucl'

        print('\n'
              'Settings that will be used:\n'
              '    fasta: {reference}\n'
              '    Blast DB type: nucl\n'
              '\n'.format(reference=args.fasta))

    utils.removeDirectory(os.path.join(args.outdir, 'pickles', ''))

    error_msg = []
    for fasta in args.fasta:
        # Create DB
        blast_db = os.path.join(args.outdir, '{blast_DB}'.format(blast_DB=os.path.basename(fasta)))
        db_exists, original_file = run_blast.check_db_exists(blast_db)
        if not db_exists and not original_file:
            db_exists = run_blast.create_blast_db(fasta, blast_db, args.type)
            if db_exists:
                print('Blast DB created for {file} in {outdir}'.format(file=fasta, outdir=args.outdir))
                # sys.exit(0)
            else:
                error_msg.append('It was not possible to create Blast DB or {}'.format(fasta))
        elif db_exists and original_file:
            error_msg.append('Blast DB already found for {file} in {outdir} as {blast_db}'.format(file=fasta,
                                                                                                  outdir=args.outdir,
                                                                                                  blast_db=blast_db))
        else:
            error_msg.append('It was found only Blast DB files or the original fasta file from which the Blast DB'
                             ' should be produced ({file}). Either include the missing files or remove the ones present'
                             ' (usually the original fasta file)'.format(file=fasta))

    if len(error_msg) == 0:
        sys.exit(0)
    else:
        sys.exit('\n'.join(error_msg))


def reads_subcommand(args):
    # if args.reference is None and args.org is None:
    #     argparse.ArgumentParser.error('--reference or --org must be provided')

    rematch_script = include_rematch_dependencies_path()

    utils.required_programs({'rematch.py': ['--version', '>=', '4.0']})

    args.fastq = [os.path.abspath(fastq.name) for fastq in args.fastq]

    if args.reference is not None:
        args.reference = [reference.name for reference in args.reference]
    else:
        args.reference, config = get_fasta_config(args.org)
        config = parse_config(config)
        args.extraSeq = config['length_extra_seq']
        args.minCovPresence = config['minimum_depth_presence']
        args.minCovCall = config['minimum_depth_call']
        args.minGeneCoverage = config['minimum_gene_coverage']
        args.typeSeparator = '_'

        print('\n'
              'Settings that will be used:\n'
              '    reference: {reference}\n'
              '    extraSeq: {extraSeq}\n'
              '    minCovPresence: {minCovPresence}\n'
              '    minCovCall: {minCovCall}\n'
              '    minGeneCoverage: {minGeneCoverage}\n'
              '    Type separator character: {typeSeparator}'
              '\n'.format(reference=args.reference, extraSeq=args.extraSeq, minCovPresence=args.minCovPresence,
                          minCovCall=args.minCovCall, minGeneCoverage=args.minGeneCoverage,
                          typeSeparator=args.typeSeparator))

    folders_2_remove = []

    # Prepare references
    references_dir = os.path.join(args.outdir, 'references', '')
    if not os.path.isdir(references_dir):
        os.makedirs(references_dir)
    folders_2_remove.append(references_dir)
    references_files, references_headers = prepare_references(args.reference, args.mapRefTogether, references_dir)

    pickles_folder = os.path.join(args.outdir, 'pickles', '')

    # Run ReMatCh
    pickle_file = os.path.join(pickles_folder, 'rematch_module.pkl')
    if args.resume and os.path.isfile(pickle_file):
        print('ReMatCh module already run')
        references_results, module_dir = utils.extractVariableFromPickle(pickle_file)
        folders_2_remove.append(module_dir)
    else:
        _, references_results, module_dir = run_rematch.run_rematch(rematch_script, args.outdir, references_files,
                                                                    args.fastq, args.threads, args.extraSeq,
                                                                    args.minCovPresence, args.minCovCall,
                                                                    args.minFrequencyDominantAllele,
                                                                    args.minGeneCoverage, args.minGeneIdentity,
                                                                    args.debug, args.doNotRemoveConsensus)
        folders_2_remove.append(module_dir)
        utils.saveVariableToPickle([references_results, module_dir], pickle_file)

    return folders_2_remove, references_results, args.reference, references_headers


def python_arguments(program_name, version):
    """
    Sets pythons arguments

    Parameters
    ----------
    program_name : str
        String with the name of the program
    version : str
        String with version

    Returns
    -------
    parser : argparse.ArgumentParser
        General argparse
    parser_reads : argparse.ArgumentParser.add_subparsers.add_parser
        Reads subparser
        For running the program using fastq files
    parser_assembly : argparse.ArgumentParser.add_subparsers.add_parser
        Assembly subparser
        For running the program using a fasta file
    parser_blast : argparse.ArgumentParser.add_subparsers.add_parser
        Blast subparser
        For creating Blast DB
    """
    parser = argparse.ArgumentParser(prog=program_name,
                                     description='Determines which reference sequence is more likely to be present in a'
                                                 ' given sample',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version',
                        version='{prog} v{version}'.format(prog=parser.prog, version=version))

    subparsers = parser.add_subparsers(title='Subcommands', description='Valid subcommands', help='Additional help')

    parser_reads = subparsers.add_parser('reads', description='Run {} using fastq files.'.format(parser.prog),
                                         help='reads --help')
    parser_assembly = subparsers.add_parser('assembly', description='Run {prog} using a fasta file. If running'
                                                                    ' multiple samples using the same DB sequence file,'
                                                                    ' consider use first "{prog} blast"'
                                                                    ' subcommand.'.format(prog=parser.prog),
                                            help='assembly --help')
    parser_blast = subparsers.add_parser('blast', description='Creates Blast DB. This is useful when running the same'
                                                              ' DB sequence file for different samples.',
                                         help='blast --help')

    parser_reads_required = parser_reads.add_argument_group('Required options')
    parser_reads_required.add_argument('-f', '--fastq', nargs='+', action=utils.required_length((1, 2), '--fastq'),
                                       type=argparse.FileType('r'), metavar='/path/to/input/file.fq.gz',
                                       help='Path to single OR paired-end fastq files. If two files are passed, they'
                                            ' will be assumed as being the paired fastq files',
                                       required=True)

    parser_reads_reference = parser_reads.add_mutually_exclusive_group(required=True)
    parser_reads_reference.add_argument('-r', '--reference', nargs='+', type=argparse.FileType('r'),
                                        metavar='/path/to/reference_sequence.fasta',
                                        help='Fasta file containing reference sequences. If more than one file is'
                                             ' passed, a type for each file will be determined. Give the files name in'
                                             ' the same order that the type must be determined.')
    parser_reads_reference.add_argument('--org', nargs=2, type=str.lower, metavar=('escherichia', 'coli'),
                                        help='Name of the organism with reference sequences provided together'
                                             ' with {} for typing ("reference_sequences" folder)'.format(parser.prog),
                                        action=utils.arguments_choices_words(get_species_allowed(), '--org'))

    parser_reads_optional_general = parser_reads.add_argument_group('General facultative options')
    parser_reads_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                               help='Path to the directory where the information will be stored'
                                                    ' (default: ./',
                                               required=False, default='.')
    parser_reads_optional_general.add_argument('-j', '--threads', type=int, metavar='N',
                                               help='Number of threads to use (default: 1)', required=False, default=1)
    parser_reads_optional_general.add_argument('--mapRefTogether', action='store_true',
                                               help='Map the reads against all references together')
    parser_reads_optional_general.add_argument('--typeSeparator', type=str, metavar='_',
                                               help='Last single character separating the general sequence header from'
                                                    ' the last part containing the type (default: _)',
                                               required=False, default='_')
    parser_reads_optional_general.add_argument('--extraSeq', type=int, metavar='N',
                                               help='Sequence length added to both ends of target sequences (usefull to'
                                                    ' improve reads mapping to the target one) that will be trimmed in'
                                                    ' ReMatCh outputs (default when not using --org: 0)',
                                               required=False, default=0)
    parser_reads_optional_general.add_argument('--minCovPresence', type=int, metavar='N',
                                               help='Reference position minimum coverage depth to consider the position'
                                                    ' to be present in the sample (default when not using --org: 5)',
                                               required=False, default=5)
    parser_reads_optional_general.add_argument('--minCovCall', type=int, metavar='N',
                                               help='Reference position minimum coverage depth to perform a base call'
                                                    ' (default when not using --org: 10)',
                                               required=False, default=10)
    parser_reads_optional_general.add_argument('--minFrequencyDominantAllele', type=float, metavar='0.6',
                                               help=argparse.SUPPRESS, required=False, default=0.6)
    # parser_optional_general.add_argument('--minFrequencyDominantAllele', type=float, metavar='0.6',
    #                                      help='Minimum relative frequency of the dominant allele coverage depth'
    #                                           ' (value between [0, 1]). Positions with lower values will be'
    #                                           ' considered as having multiple alleles (and will be coded as N)',
    #                                      required=False, default=0.6)
    parser_reads_optional_general.add_argument('--minGeneCoverage', type=int, metavar='N',
                                               help='Minimum percentage of target reference sequence covered to'
                                                    ' consider a sequence to be present (value between [0, 100])'
                                                    ' (default when not using --org: 60)',
                                               required=False, default=60)
    parser_reads_optional_general.add_argument('--minDepthCoverage', type=int, metavar='N',
                                               help='Minimum depth of coverage of target reference sequence to'
                                                    ' consider a sequence to be present (default: 2)',
                                               required=False, default=2)
    # parser_reads_optional_general.add_argument('--minGeneIdentity', type=int, metavar='N', help=argparse.SUPPRESS,
    #                                            required=False, default=80)
    parser_reads_optional_general.add_argument('--minGeneIdentity', type=int, metavar='N',
                                               help='Minimum percentage of identity of reference sequence covered to'
                                                    ' consider a gene to be present (value between [0, 100]). One INDEL'
                                                    ' will be considered as one difference (default: 80)',
                                               required=False, default=80)
    parser_reads_optional_general.add_argument('--doNotRemoveConsensus', action='store_true',
                                               help='Do not remove ReMatCh consensus sequences')
    parser_reads_optional_general.add_argument('--debug', action='store_true',
                                               help='Debug mode: do not remove temporary files')
    parser_reads_optional_general.add_argument('--resume', action='store_true',
                                               help='Resume %(prog)s')
    parser_reads_optional_general.add_argument('--notClean', action='store_true',
                                               help='Do not remove intermediate files')

    parser_assembly_required = parser_assembly.add_argument_group('Required options')
    parser_assembly_required.add_argument('-f', '--fasta', nargs=1, type=argparse.FileType('r'),
                                          metavar='/path/to/query/assembly_file.fasta',
                                          help='Path to fasta file containing the query sequences from which the'
                                               ' types should be assessed',
                                          required=True)

    parser_assembly_reference = parser_assembly.add_mutually_exclusive_group(required=True)
    parser_assembly_reference.add_argument('-b', '--blast', nargs='+', type=argparse.FileType('r'),
                                           metavar='/path/to/Blast/db.sequences.file',
                                           help='Path to DB sequence file. If Blast DB was already produced only'
                                                ' provide the file that do not end with ".n*" something (do not use for'
                                                ' example /blast_db.sequences.fasta.nhr). If no Blast DB is found for'
                                                ' the DB sequence file, one will be created in --outdir. If more than'
                                                ' one Blast DB file is passed, a type for each file will be determined.'
                                                ' Give the files in the same order that the type must be determined.')
    parser_assembly_reference.add_argument('--org', nargs=2, type=str.lower, metavar=('escherichia', 'coli'),
                                           help='Name of the organism with DB sequence file provided'
                                                ' ("reference_sequences" folder) together'
                                                ' with seq_typing.py for typing',
                                           action=utils.arguments_choices_words(get_species_allowed(), '--org'))

    parser_assembly_optional_reference = parser_assembly.add_argument_group('General facultative options')
    parser_assembly_optional_reference.add_argument('-t', '--type', choices=['nucl', 'prot'], type=str, metavar='nucl',
                                                    help='Blast DB type (available options: %(choices)s)')

    parser_assembly_optional_general = parser_assembly.add_argument_group('General facultative options')
    parser_assembly_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                                  help='Path to the directory where the information will be stored'
                                                       ' (default: ./)',
                                                  required=False, default='.')
    parser_assembly_optional_general.add_argument('-j', '--threads', type=int, metavar='N', required=False,
                                                  help='Number of threads to use (default: 1)', default=1)
    parser_assembly_optional_general.add_argument('--typeSeparator', type=str, metavar='_',
                                                  help='Last single character separating the general sequence header'
                                                       ' from the last part containing the type (default: _)',
                                                  required=False, default='_')
    parser_assembly_optional_general.add_argument('--minGeneCoverage', type=int, metavar='N',
                                                  help='Minimum percentage of target reference sequence covered to'
                                                       ' consider a sequence to be present (value between [0, 100])'
                                                       ' (default when not using --org: 60)',
                                                  required=False, default=60)
    parser_assembly_optional_general.add_argument('--minGeneIdentity', type=int, metavar='N',
                                                  help='Minimum percentage of identity of reference sequence covered'
                                                       ' to consider a gene to be present (value between [0, 100])'
                                                       ' (default: 80)',
                                                  required=False, default=80)
    parser_assembly_optional_general.add_argument('--minDepthCoverage', type=int, metavar='N', help=argparse.SUPPRESS,
                                                  required=False, default=1)
    parser_assembly_optional_general.add_argument('--debug', action='store_true',
                                                  help='Debug mode: do not remove temporary files')
    parser_assembly_optional_general.add_argument('--resume', action='store_true', help=argparse.SUPPRESS)

    parser_blast_required = parser_blast.add_argument_group('Required options')
    parser_blast_required.add_argument('-t', '--type', choices=['nucl', 'prot'], type=str, metavar='nucl',
                                       help='Blast DB type (available options: %(choices)s)', required=True)

    parser_blast_reference = parser_blast.add_mutually_exclusive_group(required=True)
    parser_blast_reference.add_argument('-f', '--fasta', nargs='+', type=argparse.FileType('r'),
                                        metavar='/path/to/db.sequences.fasta',
                                        help='Path to DB sequence file. If more than one file is passed, a Blast DB for'
                                             ' each file will be created.')
    parser_blast_reference.add_argument('--org', nargs=2, type=str.lower, metavar=('escherichia', 'coli'),
                                        help='Name of the organism with DB sequence file provided'
                                             ' ("reference_sequences" folder) together'
                                             ' with seq_typing.py for typing',
                                        action=utils.arguments_choices_words(get_species_allowed(), '--org'))

    parser_blast_optional_general = parser_blast.add_argument_group('General facultative options')
    parser_blast_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                               help='Path to the directory where the information will be stored'
                                                    ' (default: ./)',
                                               required=False, default='.')

    parser_reads.set_defaults(func=reads_subcommand)
    parser_assembly.set_defaults(func=assembly_subcommand)
    parser_blast.set_defaults(func=blast_subcommand)

    return parser, parser_reads, parser_assembly, parser_blast


def main():
    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling "python3 seq_typing.py"')

    parser, _, _, _ = python_arguments('seq_typing.py', version)
    args = parser.parse_args()

    msg = []
    if args.minGeneCoverage < 0 or args.minGeneCoverage > 100:
        msg.append('--minGeneCoverage should be a value between [0, 100]')
    if args.minGeneIdentity < 0 or args.minGeneIdentity > 100:
        msg.append('--minGeneIdentity should be a value between [0, 100]')

    if len(msg) > 0:
        argparse.ArgumentParser.error('\n'.join(msg))

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

    if not args.debug:
        for folder in folders_2_remove:
            utils.removeDirectory(folder)

    _ = utils.runTime(start_time)


if __name__ == "__main__":
    main()
