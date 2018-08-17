import os.path
from functools import partial

import modules.utils as utils


extra_blast_fields = ['query', 'q_start', 'q_end', 's_start', 's_end', 'evalue']


def get_best_sequence(data_by_gene, min_gene_coverage, min_depth_coverage):
    sequence = {}
    probable_sequences = {}
    improbable_sequences = {}

    fields = ['gene_coverage', 'gene_mean_read_coverage', 'gene_identity'] + extra_blast_fields
    for gene, rematch_results in data_by_gene.items():
        if min_depth_coverage is not None and rematch_results['gene_mean_read_coverage'] < min_depth_coverage:
            continue
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

    return ':'.join([results[reference] for reference in references_files]), results_info, probable_results, improbable_results


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
        Return the input dictionary already with ['query', 'q_start', 'q_end', 's_start', 's_end', 'evalue'] fields
        as 'NA'
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


def write_reports(outdir, seq_type, seq_type_info, probable_results, improbable_results, type_separator):
    with open(os.path.join(outdir, 'seq_typing.report.txt'), 'wt') as writer:
        writer.write(seq_type + '\n')

    with open(os.path.join(outdir, 'seq_typing.report_types.tab'), 'wt') as writer:
        writer.write('\t'.join(['#sequence_type', 'type', 'reference_file', 'sequence', 'sequenced_covered',
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


module_timer = partial(utils.timer, name='Module Parse results')


@module_timer
def parse_results(references_results, references_files, references_headers, outdir, min_gene_coverage,
                  min_depth_coverage, type_separator):
    try:
        references_results = split_references_results_by_references(references_results, references_headers)
    except TypeError:
        print('Parsing assembly results')
    else:
        print('Parsing reads results')
    finally:
        seq_type, seq_type_info, probable_results, improbable_results = get_results(references_results,
                                                                                    min_gene_coverage,
                                                                                    min_depth_coverage, type_separator,
                                                                                    references_files,
                                                                                    references_headers)
        write_reports(outdir, seq_type, seq_type_info, probable_results, improbable_results, type_separator)

    return seq_type, seq_type_info, probable_results, improbable_results
