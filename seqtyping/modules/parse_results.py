import os.path
from functools import partial
from itertools import groupby as itertools_groupby

try:
    import modules.utils as utils
except ImportError:
    from seqtyping.modules import utils as utils


extra_blast_fields = ['query', 'q_start', 'q_end', 'q_length', 's_start', 's_end', 's_length', 'evalue']


def get_best_sequence(data_by_gene, min_gene_coverage, min_depth_coverage):
    sequence = {}
    probable_sequences = {}
    improbable_sequences = {}

    fields = ['gene_coverage', 'gene_mean_read_coverage', 'gene_identity'] + extra_blast_fields
    for gene, rematch_results in data_by_gene.items():
        if rematch_results['gene_mean_read_coverage'] < min_depth_coverage:
            improbable_sequences[rematch_results['gene_coverage']] = gene
        else:
            if rematch_results['gene_coverage'] >= min_gene_coverage:
                if rematch_results['gene_coverage'] not in sequence:
                    sequence[rematch_results['gene_coverage']] = gene
                else:
                    if data_by_gene[sequence[rematch_results['gene_coverage']]]['gene_mean_read_coverage'] < \
                            rematch_results['gene_mean_read_coverage']:
                        probable_sequences[sequence[rematch_results['gene_coverage']]] = \
                            [rematch_results['gene_coverage']] + \
                            [data_by_gene[sequence[rematch_results['gene_coverage']]][field] for field in fields[1:]]
                        sequence[rematch_results['gene_coverage']] = gene
                    else:
                        probable_sequences[gene] = [rematch_results[field] for field in fields]
            else:
                improbable_sequences[rematch_results['gene_coverage']] = gene

    if len(sequence) == 1:
        sequence = next(iter(sequence.values()))  # Get the first one
    elif len(sequence) == 0:
        sequence = None
    else:
        for gene in [sequence[i] for i in sorted(sequence.keys(), reverse=True)[1:]]:
            probable_sequences[gene] = [data_by_gene[gene][field] for field in fields]
        sequence = sequence[sorted(sequence.keys(), reverse=True)[0]]

    if len(improbable_sequences) > 0:
        improbable_sequences = improbable_sequences[sorted(improbable_sequences.keys(), reverse=True)[0]]
    else:
        improbable_sequences = None

    return sequence, probable_sequences, improbable_sequences


def get_results(references_results, min_gene_coverage, min_depth_coverage, type_separator, references_files,
                references_headers):
    intermediate_results = {}
    for reference, data_by_gene in references_results.items():
        intermediate_results[reference] = get_best_sequence(data_by_gene, min_gene_coverage, min_depth_coverage)

    results = {}
    results_info = {}
    probable_results = {}
    improbable_results = {}
    for original_reference in references_headers.keys():
        results[original_reference] = 'NT'  # For None Typeable
        probable_results[original_reference] = []
        improbable_results[original_reference] = []

    fields = ['gene_coverage', 'gene_mean_read_coverage', 'gene_identity'] + extra_blast_fields
    for reference, data in intermediate_results.items():
        sequence = data[0]
        probable_sequences = data[1]
        improbable_sequences = data[2]
        if sequence is not None:
            for original_reference, headers in references_headers.items():
                for new_header, original_header in headers.items():
                    if sequence == new_header:
                        results[original_reference] = original_header.rsplit(type_separator, 1)[1]
                        results_info[original_reference] = \
                            [original_header] + \
                            [references_results[reference][new_header][field] for field in fields]
                    if len(probable_sequences) > 0:
                        if new_header in probable_sequences:
                            probable_results[original_reference].append([original_header] +
                                                                        [probable_sequences[new_header][x]
                                                                         for x in range(0, len(fields))])
                    if sequence == new_header and \
                            (len(probable_sequences) == 0 or
                             len(probable_results[original_reference]) == len(probable_sequences)):
                        break
        else:
            if improbable_sequences is not None:
                for original_reference, headers in references_headers.items():
                    for new_header, original_header in headers.items():
                        if improbable_sequences == new_header:
                            print('\n' +
                                  '\n'.join(['NONE TYPEABLE REFERENCE',
                                             'Reference file: {}'.format(original_reference),
                                             'Most likely type: {}'.format(original_header.rsplit(type_separator, 1)[1]),
                                             'Most likely sequence: {}'.format(original_header),
                                             'Sequenced covered: {}'.format(
                                                 references_results[reference][new_header]['gene_coverage']),
                                             'Coverage depth: {}'.format(
                                                 references_results[reference][new_header]['gene_mean_read_coverage']),
                                             'Sequence identity: {}'.format(
                                                 references_results[reference][new_header]['gene_identity'])]) +
                                  '\n')
                            improbable_results[original_reference] = \
                                [original_header] + [references_results[reference][new_header][field]
                                                     for field in fields]

    return ':'.join([results[reference] for reference in references_files]), results_info, probable_results, \
           improbable_results


def add_empty_blast_info(data_dict):
    """
    Add empty extra Blast results to read mapping results

    Parameters
    ----------
    data_dict : dict
        ReMatCh results for a single sequence (already without the counter)

    Returns
    -------
    data_dict : dict
        Return the input dictionary already with ['query', 'q_start', 'q_end', 'q_length', 's_start', 's_end',
        's_length', 'evalue'] fields as 'NA'
    """
    for field in extra_blast_fields:
        data_dict[field] = 'NA'
    return data_dict


def split_references_results_by_references(references_results, references_headers):
    organized_references_results = {}
    if len(references_headers) > 1:
        for reference, pickleFile in references_results.items():
            data_by_gene = utils.extractVariableFromPickle(pickleFile)
            for counter, data in data_by_gene.items():
                for reference_file, headers in references_headers.items():
                    if reference_file not in organized_references_results:
                        organized_references_results[reference_file] = {}
                    if data['header'] in headers.keys():
                        organized_references_results[reference_file][data['header']] = add_empty_blast_info(data)
    else:
        organized_references_results[next(iter(references_results.keys()))] = {}
        for reference, pickleFile in references_results.items():
            data_by_gene = utils.extractVariableFromPickle(pickleFile)
            for counter, data in data_by_gene.items():
                organized_references_results[reference][data['header']] = add_empty_blast_info(data)
    return organized_references_results


def save_new_allele_reads(sample, new_allele_dir, rematch_consensus, sequence_selected, selected_type):
    """
    Save the new allele found using reads in a reference file folder, under the file name of the selected type and with
    sample's name as header

    Parameters
    ----------
    sample : str
        Sample's name
    new_allele_dir : str
        Path to the folder where the new allele will be saved
    rematch_consensus : str
        Path to ReMatCh consensus file for the reference file used and from which the new allele will be retreived
    sequence_selected : str
        Header of the selected sequence
    selected_type : str
        Type from the selected sequence

    Returns
    -------

    """
    reference_file = os.path.basename(rematch_consensus)
    os.makedirs(os.path.join(new_allele_dir, reference_file, ''))
    with open(rematch_consensus, mode='rt', newline=None) as reader:
        for line in reader:
            if line.startswith('>'):
                header = line.rstrip('\r\n')[1:]
                if header == sequence_selected:
                    seq = ''
                    for line in reader:
                        if not line.startswith('>'):
                            seq += line.rstrip('\r\n')
                        else:
                            break
                    with open(os.path.join(new_allele_dir, reference_file,
                                           '{selected_type}.fasta'.format(selected_type=selected_type)), 'wt') as writer:
                        writer.write('>{}\n'.format(sample))
                        writer.write('\n'.join(utils.chunkstring(seq, 80)) + '\n')
                    break


def save_new_allele_assembly(sample, new_allele_dir, reference_file, query, selected_type, assembly,
                             q_start, q_end, q_length, s_start, s_end, s_length):
    """
    Save the new allele found using the assembly in a reference file folder, under the file name of the selected type
    and with sample's name as header

    Parameters
    ----------
    sample : str
        Sample's name
    new_allele_dir : str
        Path to the folder where the new allele will be saved
    reference_file : str
        Path to the reference file used as subject in Blast search
    query : str
        Header of the assembly that give a hit
    selected_type : str
        Type from the selected sequence
    assembly : str
        Path to the assembly used as query in Blast search
    q_start : int
        Query hit starting position
    q_end : int
        Query hit ending position
    q_length : int
        Query length
    s_start : int
        Subject hit starting position
    s_end : int
        Subject hit ending position
    s_length : int
        Subject length

    Returns
    -------

    """

    # TODO: remove code redundancy when have time (maybe never!)

    found_query = False

    os.makedirs(os.path.join(new_allele_dir, os.path.basename(reference_file), ''))
    reader = open(assembly, mode='rt', newline=None)
    fasta_iter = (g for k, g in itertools_groupby(reader, lambda x: x.startswith('>')))
    for header in fasta_iter:
        header = header.__next__()[1:].rstrip('\r\n')
        seq = ''.join(s.rstrip('\r\n') for s in fasta_iter.__next__())
        if header == query:
            found_query = True
            if s_start > s_end:
                starting_position = q_start - (s_length - s_start)
                ending_position = q_end + s_end - 1
            else:
                starting_position = q_start - (s_start - 1)
                ending_position = q_end + (s_length - s_end)

            partial_seq = ''
            if starting_position < 1:
                starting_position = 1
                partial_seq = '_partial'
            if ending_position > q_length:
                ending_position = q_length
                partial_seq = '_partial'

            seq = seq[(starting_position - 1):ending_position]
            if s_start > s_end:
                seq = utils.reverse_complement(seq)

            with open(os.path.join(new_allele_dir, os.path.basename(reference_file),
                                   '{selected_type}.fasta'.format(selected_type=selected_type)), 'wt') as writer:
                writer.write('>{sample}{partial_seq}\n'.format(sample=sample, partial_seq=partial_seq))
                writer.write('\n'.join(utils.chunkstring(seq, 80)) + '\n')
            break
    reader.close()

    if not found_query:
        print('New allele found, but the Blast query name was not exactly the same as the header!\n'
              'Trying looking for a sequence starting with the Blast query name...')
        os.makedirs(os.path.join(new_allele_dir, os.path.basename(reference_file), ''))
        reader = open(assembly, mode='rt', newline=None)
        fasta_iter = (g for k, g in itertools_groupby(reader, lambda x: x.startswith('>')))
        for header in fasta_iter:
            header = header.__next__()[1:].rstrip('\r\n')
            seq = ''.join(s.rstrip('\r\n') for s in fasta_iter.__next__())
            if header.startswith(query):
                if s_start > s_end:
                    starting_position = q_start - (s_length - s_start)
                    ending_position = q_end + s_end - 1
                else:
                    starting_position = q_start - (s_start - 1)
                    ending_position = q_end + (s_length - s_end)

                partial_seq = ''
                if starting_position < 1:
                    starting_position = 1
                    partial_seq = '_partial'
                if ending_position > q_length:
                    ending_position = q_length
                    partial_seq = '_partial'

                seq = seq[(starting_position - 1):ending_position]
                if s_start > s_end:
                    seq = utils.reverse_complement(seq)

                with open(os.path.join(new_allele_dir, os.path.basename(reference_file),
                                       '{selected_type}.fasta'.format(selected_type=selected_type)), 'wt') as writer:
                    writer.write('>{sample}{partial_seq}\n'.format(sample=sample, partial_seq=partial_seq))
                    writer.write('\n'.join(utils.chunkstring(seq, 80)) + '\n')
                break
        reader.close()


def write_reports(outdir, seq_type, seq_type_info, probable_results, improbable_results, type_separator, assembly,
                  sample='sample', save_new_allele=False):
    new_allele_found = False
    new_allele_dir = os.path.join(outdir, 'new_allele', '')
    if save_new_allele:
        utils.removeDirectory(new_allele_dir)
        os.makedirs(new_allele_dir)

    with open(os.path.join(outdir, 'seq_typing.report.txt'), 'wt') as writer:
        writer.write(seq_type + '\n')

    with open(os.path.join(outdir, 'seq_typing.report_types.tab'), 'wt') as writer:
        writer.write('\t'.join(['#sequence_type', 'reference_file', 'type', 'sequence', 'sequenced_covered',
                                'coverage_depth', 'sequence_identity'] + extra_blast_fields) +
                     '\n')

        print('\n' + 'TYPEABLE REFERENCES')
        typeable_references = False
        for reference, data in seq_type_info.items():
            if len(data) > 0:
                print('\n' +
                      '\n'.join(['Reference file: {}'.format(reference),
                                 'Type: {}'.format(data[0].rsplit(type_separator, 1)[1]),
                                 'Sequence: {}'.format(data[0]),
                                 'Sequenced covered: {}'.format(data[1]),
                                 'Coverage depth: {}'.format(data[2]),
                                 'Sequence identity: {}'.format(data[3])]) +
                      '\n')
                writer.write('\t'.join(['selected', reference, data[0].rsplit(type_separator, 1)[1]] +
                                       list(map(str, data))) + '\n')

                if save_new_allele and (data[1] < 100 and data[3] < 100):
                    new_allele_found = True
                    if assembly is None:
                        rematch_consensus = os.path.join(outdir, 'rematch', 'new_allele',
                                                         os.path.basename(reference))
                        save_new_allele_reads(sample=sample, new_allele_dir=new_allele_dir,
                                              rematch_consensus=rematch_consensus, sequence_selected=data[0],
                                              selected_type=data[0].rsplit(type_separator, 1)[1])
                    else:
                        save_new_allele_assembly(sample=sample, new_allele_dir=new_allele_dir, reference_file=reference,
                                                 query=data[4], selected_type=data[0].rsplit(type_separator, 1)[1],
                                                 assembly=assembly, q_start=data[5], q_end=data[6], q_length=data[7],
                                                 s_start=data[8], s_end=data[9], s_length=data[10])

                typeable_references = True

        if typeable_references is False:
            print('No references return a type')

        print('\n' + 'Types found:' + '\n')
        print(seq_type + '\n')

        for reference, types in probable_results.items():
            if len(types) > 0:
                print('\n' +
                      'Other possible types found for {reference}!\n'
                      'Check seq_typing.report.other_probable_types.tab file.'.format(reference=reference) +
                      '\n')
                for probable_type in types:
                    writer.write('\t'.join(['other_probable_type',
                                            reference,
                                            probable_type[0].rsplit(type_separator, 1)[1]] +
                                           list(map(str, probable_type))) + '\n')

        for reference, data in improbable_results.items():
            if len(data) > 0:
                writer.write('\t'.join(['most_likely', reference, data[0].rsplit(type_separator, 1)[1]] +
                                       list(map(str, data))) + '\n')

    if not new_allele_found and os.path.isdir(new_allele_dir):
        utils.removeDirectory(new_allele_dir)

    if save_new_allele and assembly is None:
        utils.removeDirectory(os.path.join(outdir, 'rematch', 'new_allele', ''))


@utils.timer('Module Parse results')
def parse_results(references_results, references_files, references_headers, outdir, min_gene_coverage,
                  min_depth_coverage, type_separator, sample='sample', save_new_allele=False, assembly=None):
    if assembly is None:
        print('Parsing reads results')
        references_results = split_references_results_by_references(references_results, references_headers)
    else:
        print('Parsing assembly results')

    seq_type, seq_type_info, probable_results, improbable_results = get_results(references_results,
                                                                                min_gene_coverage,
                                                                                min_depth_coverage, type_separator,
                                                                                references_files,
                                                                                references_headers)
    write_reports(outdir, seq_type, seq_type_info, probable_results, improbable_results, type_separator, assembly,
                  sample=sample, save_new_allele=save_new_allele)

    return seq_type, seq_type_info, probable_results, improbable_results
