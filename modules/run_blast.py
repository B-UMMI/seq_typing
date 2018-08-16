import os.path
import sys
from modules.utils import runCommandPopenCommunicate as RUN_subprocess


def check_db_exists(db_path):
    """
    Check if Blast DB exists

    Parameters
    ----------
    db_path : str
        Path to Blast DB files, e.g. /output/blast_db.nucl.sequences.fasta

    Returns
    -------
    db_exists : bool
        Tells if the Blast DB already exists
    original_file : bool
        Tells if the original fasta file from which the Blast DB was produced is present
    """
    db_exists = False
    files = [f for f in os.listdir(os.path.dirname(db_path))
             if not f.startswith('.')
             and os.path.isfile(os.path.join(os.path.dirname(db_path), f))]
    counter = 0
    original_file = False
    for file_found in files:
        if file_found == os.path.basename(db_path):
            original_file = True
        elif os.path.splitext(file_found)[0] == os.path.basename(db_path):
            counter += 1
    if counter > 0:
        db_exists = True
    return db_exists, original_file


def create_blast_db(db_sequences, db_output, db_type):
    """
    Creates a Blast DB

    Parameters
    ----------
    db_sequences : str
        Path to fasta file containing the sequences from which Blast DB will be created, e.g. /input/sequences.fasta
    db_output : str
        Path to Blast output DB files, e.g. /output/blast_db.nucl.sequences.fasta
    db_type : str
        Blast DB type. Can only be 'nucl' or 'prot'

    Returns
    -------
    run_successfully : bool
        Tells if the Blast DB was successfully created
    """

    run_successfully = False

    # Check Blast DB type
    if db_type not in ('nucl', 'prot'):
        exit('Wrong Blast DB type provided ({db_type}). Use one of the following: nucl, prot'.format(db_type=db_type))

    # Check if db_output dir exists
    if not os.path.isdir(os.path.dirname(db_output)):
        os.makedirs(os.path.dirname(db_output))
    else:
        db_exists, original_file = check_db_exists(db_output)
        if db_exists and original_file:
            print('Blast DB already found at {db_output}'
                  ' for {db_sequences}'.format(db_output=os.path.dirname(db_output),
                                               db_sequences=db_sequences,
                                               file_found=os.path.basename(db_sequences)))

            run_successfully = True
        elif db_exists and not original_file:
            print('The original fasta file ({db_sequences}) from which the Blast DB was produced is not present. Make'
                  ' sure it is found in {db_dir} and it is'
                  ' named {original_file_name}'.format(db_dir=os.path.dirname(db_output), db_sequences=db_sequences,
                                                       original_file_name=os.path.basename(db_output)))
        elif not db_exists and original_file:
            print('The original fasta file ({db_sequences}) from which the Blast DB was supposed to be produced is'
                  ' present in Blast DB directory ({db_dir}), but no Blast DB files were found'
                  ' there.'.format(db_sequences=db_sequences, db_dir=os.path.dirname(db_output)))
        else:
            run_successfully, _, _ = RUN_subprocess(['makeblastdb', '-parse_seqids', '-dbtype', db_type, '-in',
                                                     db_sequences, '-out', db_output], False, None, True)
            if run_successfully:
                from shutil import copyfile
                copyfile(db_sequences, db_output)

    return run_successfully


def run_blast_command(query_file, blast_db, db_type, blast_output, threads=1):
    """
    Run Blast: blastn or blastp

    Parameters
    ----------
    query_file : str
        Path to fasta file containing the query sequences, e.g. /input/queries.fasta
    blast_db : str
        Path to Blast DB files, e.g. /input/blast_db.nucl.sequences.fasta
    db_type : str
        Blast DB type. Can only be 'nucl' or 'prot'
    blast_output : str
        Path to Blast output tabular file
    threads : int
        Number of CPUs to use during Blast

    Returns
    -------
    run_successfully : bool
        Tells if the Blast DB was successfully created
    """

    # Check Blast DB type
    if db_type not in ('nucl', 'prot'):
        exit('Wrong Blast DB type provided ({db_type}).'
             ' Use one of the following: "nucl", "prot"'.format(db_type=db_type))

    # Check if db_output dir exists
    if not os.path.isdir(os.path.dirname(blast_output)):
        os.makedirs(os.path.dirname(blast_output))

    command = ['', '-query', query_file, '-db', blast_db, '-out', blast_output, '-outfmt', '', '-dust', 'no',
               '-culling_limit', '1', '-num_threads', str(threads)]

    if db_type == 'nucl':
        command[0] = 'blastn -task blastn'
    else:
        command[0] = 'blastp -task blastp'

    command[8] = "'7 qseqid qlen sseqid slen qstart qend sstart send evalue length pident nident mismatch gaps'"

    run_successfully, _, _ = RUN_subprocess(command, False, None, True)

    return run_successfully


def parse_blast_output(blast_output):
    """
    Parse Blast output

    Parameters
    ----------
    blast_output : str
        Path to Blast output tabular file

    Returns
    -------
    run_successfully : bool
        Tells if the parser worked
    """

    # TODO: update Retuns

    # Blast fields
    # 'query_id', 'query_length', 'subject_id', 'subject_length', 'q_start', 'q_end', 's_start', 's_end', 'evalue',
    # 'alignment_length', 'percentage_identity', 'identical', 'mismatches', 'gaps'

    output_blast = {}
    # To make it compatible with ReMatCh results
    subject_index = {}

    # Read Blast output
    with open(blast_output, 'rt') as reader:
        for line in reader:
            line = line.rstrip('\r\n')
            if len(line) > 0:
                if not line.startswith('#'):
                    line = line.split('\t')

                    # By subject
                    counter = len(subject_index) + 1
                    if line[2] not in subject_index:
                        subject_index[line[2]] = counter
                        output_blast[counter] = {'header': line[2],
                                                 'gene_coverage': int(line[9]) / float(line[3]),
                                                 'gene_low_coverage': 0,
                                                 'gene_number_positions_multiple_alleles': 0,
                                                 'gene_mean_read_coverage': 1,
                                                 'gene_identity': float(line[10]),
                                                 'q_start': int(line[4]),
                                                 'q_end': int(line[5]), 's_start': int(line[6]),
                                                 's_end': int(line[7]), 'evalue': float(line[8]),
                                                 'gaps': int(line[13]), 'query': line[0],
                                                 'alignment_length': int(line[9]), 'subject_length': int(line[3])}
                    else:
                        previous_blast_results = output_blast[subject_index[line[2]]]
                        present_blast_results = {'header': line[2], 'gene_coverage': int(line[9]) / float(line[3]),
                                                 'gene_low_coverage': 0, 'gene_number_positions_multiple_alleles': 0,
                                                 'gene_mean_read_coverage': 1, 'gene_identity': float(line[10]),
                                                 'q_start': int(line[4]), 'q_end': int(line[5]),
                                                 's_start': int(line[6]), 's_end': int(line[7]),
                                                 'evalue': float(line[8]),
                                                 'gaps': int(line[13]), 'query': line[0],
                                                 'alignment_length': int(line[9]), 'subject_length': int(line[3])}

                        to_change = False
                        if previous_blast_results['alignment_length'] < int(line[9]):
                            to_change = True
                        elif previous_blast_results['alignment_length'] == int(line[9]) and \
                                previous_blast_results['gene_coverage'] <= present_blast_results['gene_coverage']:
                            if previous_blast_results['gene_coverage'] < present_blast_results['gene_coverage']:
                                to_change = True
                            elif previous_blast_results['gene_coverage'] == present_blast_results['gene_coverage'] and \
                                    previous_blast_results['gene_identity'] <= float(line[10]):
                                if previous_blast_results['gene_identity'] < float(line[10]):
                                    to_change = True
                                elif previous_blast_results['gene_identity'] == float(line[10]) and \
                                        previous_blast_results['evalue'] >= float(line[8]):
                                    if previous_blast_results['evalue'] > float(line[8]):
                                        to_change = True
                                    elif previous_blast_results['evalue'] == float(line[8]) and \
                                            previous_blast_results['gaps'] > int(line[13]):
                                        to_change = True

                        if to_change:
                            output_blast[subject_index[line[2]]] = present_blast_results

    return output_blast


def run_blast(fasta_file, db_path, threads, outdir):
    # references_results = {'/mnt/beegfs/scratch/ONEIDA/mpmachado/test_seq_typing_dengue/test_out/references/references_together.fasta': '/mnt/beegfs/scratch/ONEIDA/mpmachado/test_seq_typing_dengue/test_out/rematch/references_together.fasta_0.pkl'}

    # references_headers = {'/mnt/beegfs/scratch/ONEIDA/mpmachado/test_seq_typing_dengue/subtype_2.fasta': {'gb_KY474330_Organism_Dengue_virus_2_Strain_ame_00104-P_Segment_null_Subtype_2-AsianAmerican_Host_Human_seqTyping_2-AsianAmerican': 'gb:KY474330|Organism:Dengue_virus_2|Strain_ame:00104-P|Segment:null|Subtype:2-AsianAmerican|Host:Human_seqTyping_2-AsianAmerican', 'gb_KY474334_Organism_Dengue_virus_2_Strain_ame_00121-P_Segment_null_Subtype_2-AsianAmerican_Host_Human_seqTyping_2-AsianAmerican': 'gb:KY474334|Organism:Dengue_virus_2|Strain_ame:00121-P|Segment:null|Subtype:2-AsianAmerican|Host:Human_seqTyping_2-AsianAmerican', 'gb_KY794785_Organism_Dengue_virus_2_Strain_ame_100511_Segment_null_Subtype_2-Cosmopolitan_Host_Human_seqTyping_2-Cosmopolitan': 'gb:KY794785|Organism:Dengue_virus_2|Strain_ame:100511|Segment:null|Subtype:2-Cosmopolitan|Host:Human_seqTyping_2-Cosmopolitan'}, '/mnt/beegfs/scratch/ONEIDA/mpmachado/test_seq_typing_dengue/subtype_3.fasta': {'gb_KY586703_Organism_Dengue_virus_Strain_ame_Ser3_Thailand_Bangkok_Seq1_Segment_null_Subtype_3-II_Host_Human_seqTyping_3-II': 'gb:KY586703|Organism:Dengue_virus|Strain_ame:Ser3_Thailand_Bangkok_Seq1|Segment:null|Subtype:3-II|Host:Human_seqTyping_3-II', 'gb_KY586715_Organism_Dengue_virus_Strain_Name_Ser3_Thailand_Bangkok_Seq10_Segment_null_Subtype_3-II_Host_Human_seqTyping_3-II': 'gb:KY586715|Organism:Dengue_virus|Strain_Name:Ser3_Thailand_Bangkok_Seq10|Segment:null|Subtype:3-II|Host:Human_seqTyping_3-II'}}

    # data_by_gene = {1: {'header': 'gb_KY474330_Organism_Dengue_virus_2_Strain_ame_00104-P_Segment_null_Subtype_2-AsianAmerican_Host_Human_seqTyping_2-AsianAmerican', 'gene_coverage': 10.763786493659694, 'gene_low_coverage': 19.086757990867582, 'gene_number_positions_multiple_alleles': 0, 'gene_mean_read_coverage': 77.5123287671233, 'gene_identity': 85.02283105022832}, 2: {'header': 'gb_KY474334_Organism_Dengue_virus_2_Strain_ame_00121-P_Segment_null_Subtype_2-AsianAmerican_Host_Human_seqTyping_2-AsianAmerican', 'gene_coverage': 10.999705101739892, 'gene_low_coverage': 17.247542448614837, 'gene_number_positions_multiple_alleles': 0, 'gene_mean_read_coverage': 76.30741733690796, 'gene_identity': 84.62913315460233}, 3: {'header': 'gb_KY794785_Organism_Dengue_virus_2_Strain_ame_100511_Segment_null_Subtype_2-Cosmopolitan_Host_Human_seqTyping_2-Cosmopolitan', 'gene_coverage': 28.87053966381599, 'gene_low_coverage': 15.253660197480423, 'gene_number_positions_multiple_alleles': 13, 'gene_mean_read_coverage': 801.5263874702077, 'gene_identity': 81.68198842356145}, 4: {'header': 'gb_KY586703_Organism_Dengue_virus_Strain_ame_Ser3_Thailand_Bangkok_Seq1_Segment_null_Subtype_3-II_Host_Human_seqTyping_3-II', 'gene_coverage': 99.65753424657534, 'gene_low_coverage': 0.0, 'gene_number_positions_multiple_alleles': 1, 'gene_mean_read_coverage': 457.106529209622, 'gene_identity': 96.56357388316151}, 5: {'header': 'gb_KY586715_Organism_Dengue_virus_Strain_Name_Ser3_Thailand_Bangkok_Seq10_Segment_null_Subtype_3-II_Host_Human_seqTyping_3-II', 'gene_coverage': 100.0, 'gene_low_coverage': 0.0, 'gene_number_positions_multiple_alleles': 0, 'gene_mean_read_coverage': 7656.738151425762, 'gene_identity': 93.31366764995083}}
    # gene_coverage, gene_identity, gene_mean_read_coverage
    # gene_mean_read_coverage = 1

    # references_files = ['/mnt/beegfs/scratch/ONEIDA/mpmachado/test_seq_typing_dengue/subtype_2.fasta', '/mnt/beegfs/scratch/ONEIDA/mpmachado/test_seq_typing_dengue/subtype_3.fasta']
    pass
