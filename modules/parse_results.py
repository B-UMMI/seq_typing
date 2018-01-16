import os.path

import modules.utils as utils


def get_best_sequence(data_by_gene, minGeneCoverage):
    sequence = {}
    probable_sequences = {}

    for gene, rematch_results in data_by_gene.items():
        if rematch_results['gene_coverage'] >= minGeneCoverage:
            if rematch_results['gene_coverage'] not in sequence:
                sequence[rematch_results['gene_coverage']] = gene
            else:
                if data_by_gene[sequence[rematch_results['gene_coverage']]]['gene_mean_read_coverage'] < rematch_results['gene_mean_read_coverage']:
                    probable_sequences[sequence[rematch_results['gene_coverage']]] = (rematch_results['gene_coverage'], data_by_gene[sequence[rematch_results['gene_coverage']]]['gene_mean_read_coverage'], data_by_gene[sequence[rematch_results['gene_coverage']]]['gene_identity'])
                    sequence[rematch_results['gene_coverage']] = gene
                else:
                    probable_sequences[gene] = (rematch_results['gene_coverage'], rematch_results['gene_mean_read_coverage'], rematch_results['gene_identity'])

    if len(sequence) == 1:
        sequence = next(iter(sequence.values()))  # Get the first one
    elif len(sequence) == 0:
        sequence = None
    else:
        for gene in [sequence[i] for i in sorted(sequence.keys(), reverse=True)[1:]]:
            probable_sequences[gene] = (data_by_gene[gene]['gene_coverage'], data_by_gene[gene]['gene_mean_read_coverage'], data_by_gene[gene]['gene_identity'])
        sequence = sequence[sorted(sequence.keys(), reverse=True)[0]]

    return sequence, probable_sequences


def get_results(references_results, minGeneCoverage, typeSeparator, references_files, references_headers):
    intermediate_results = {}
    for reference, data_by_gene in references_results.items():
        intermediate_results[reference] = get_best_sequence(data_by_gene, minGeneCoverage)

    results = {}
    results_info = {}
    probable_results = {}
    for original_reference in references_headers.keys():
        results[original_reference] = 'NT'  # For None Typeable
        probable_results[original_reference] = []

    for reference, data in intermediate_results.items():
        sequence = data[0]
        probable_sequences = data[1]
        if sequence is not None:
            for original_reference, headers in references_headers.items():
                for new_header, original_header in headers.items():
                    if sequence == new_header:
                        results[original_reference] = original_header.rsplit('_', 1)[1]
                        results_info[original_reference] = (original_header, references_results[original_reference][new_header]['gene_coverage'], references_results[original_reference][new_header]['gene_mean_read_coverage'], references_results[original_reference][new_header]['gene_identity'])
                    if len(probable_sequences) > 0:
                        if new_header in probable_sequences:
                            probable_results[original_reference].append((original_header, probable_sequences[new_header][0], probable_sequences[new_header][1]), probable_sequences[new_header][2])
                    if sequence == new_header and (len(probable_sequences) == 0 or len(probable_results[original_reference]) == len(probable_sequences)):
                        break

    return ':'.join([results[reference] for reference in references_files]), results_info, probable_results


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
                        organized_references_results[reference_file][data['header']] = data
    else:
        organized_references_results[next(iter(references_results.keys()))] = {}
        for reference, pickleFile in references_results.items():
            data_by_gene = utils.extractVariableFromPickle(pickleFile)
            for counter, data in data_by_gene.items():
                organized_references_results[reference][data['header']] = data
    return organized_references_results


def write_reports(outdir, seq_type, seq_type_info, probable_results):
    with open(os.path.join(outdir, 'seq_typing.report.txt'), 'wt') as writer:
        writer.write(seq_type + '\n')
        print('\n' + 'Types found:' + '\n')
        print(seq_type + '\n')
        for reference, data in seq_type_info.items():
            print('\n' + '\n'.join(['Reference_file: {}'.format(reference), 'Sequence: {}'.format(data[0]), 'Sequenced covered: {}'.format(data[1]), 'Coverage depth: {}'.format(data[2]), 'Sequence identity: {}'.format(data[3])]) + '\n')
    with open(os.path.join(outdir, 'seq_typing.report.other_probable_types.tab'), 'wt') as writer:
        header_other_probable_types = False
        for reference, types in probable_results.items():
            if len(types) > 0:
                if not header_other_probable_types:
                    writer.write('\t'.join(['#reference_file', 'sequence', 'sequenced_covered', 'coverage_depth', 'sequence_identity']) + '\n')
                    header_other_probable_types = True
                    print('\n' + 'Other possible types found! Check seq_typing.report.other_probable_types.tab file)' + '\n')
                for probable_type in types:
                    writer.write('\t'.join([reference] + map(str, probable_type)) + '\n')


def parse_results(references_results, references_files, references_headers, outdir, minGeneCoverage, typeSeparator):
    references_results = split_references_results_by_references(references_results, references_headers)
    seq_type, seq_type_info, probable_results = get_results(references_results, minGeneCoverage, typeSeparator, references_files, references_headers)
    write_reports(outdir, seq_type, seq_type_info, probable_results)
    return seq_type, seq_type_info, probable_results


# references_headers = {reference_file: {new_header: original_header}}

# references_results = {reference_mapped: {counter: {'header': new_header, 'gene_coverage': 0.0, 'gene_low_coverage': 0, 'gene_number_positions_multiple_alleles': 0, 'gene_mean_read_coverage': 0, 'gene_identity': 0}}}
