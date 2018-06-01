import os.path
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
    """
    db_exists = False
    files = [f for f in os.listdir(os.path.dirname(db_path))
             if not f.startswith('.')
             and os.path.isfile(os.path.join(os.path.dirname(db_path), f))]
    counter = 0
    for file_found in files:
        if os.path.splitext(file_found)[0] == os.path.basename(db_path):
            counter += 1
    if counter > 1:
        db_exists = True
    return db_exists


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

    # Check Blast DB type
    if db_type not in ('nucl', 'prot'):
        exit('Wrong Blast DB type provided ({db_type}). Use one of the following: nucl, prot'.format(db_type=db_type))

    # Check if db_output dir exists
    if not os.path.isdir(os.path.dirname(db_output)):
        os.makedirs(os.path.dirname(db_output))
    else:
        if check_db_exists(db_output):
            exit('Blast DB already found at {dirname} for'
                 ' {basename}: {file_found}'.format(dirname=os.path.dirname(db_output),
                                                    basename=os.path.basename(db_output), file_found=file_found))

    run_successfully, _, _ = RUN_subprocess(['makeblastdb', '-parse_seqids', '-dbtype', db_type, '-in', db_sequences,
                                             '-out', db_output], False, None, True)

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


def run_blast(fasta_file, db_path, threads, outdir):
    # references_results = {'/mnt/beegfs/scratch/ONEIDA/mpmachado/test_seq_typing_dengue/test_out/references/references_together.fasta': '/mnt/beegfs/scratch/ONEIDA/mpmachado/test_seq_typing_dengue/test_out/rematch/references_together.fasta_0.pkl'}

    # references_headers = {'/mnt/beegfs/scratch/ONEIDA/mpmachado/test_seq_typing_dengue/subtype_2.fasta': {'gb_KY474330_Organism_Dengue_virus_2_Strain_ame_00104-P_Segment_null_Subtype_2-AsianAmerican_Host_Human_seqTyping_2-AsianAmerican': 'gb:KY474330|Organism:Dengue_virus_2|Strain_ame:00104-P|Segment:null|Subtype:2-AsianAmerican|Host:Human_seqTyping_2-AsianAmerican', 'gb_KY474334_Organism_Dengue_virus_2_Strain_ame_00121-P_Segment_null_Subtype_2-AsianAmerican_Host_Human_seqTyping_2-AsianAmerican': 'gb:KY474334|Organism:Dengue_virus_2|Strain_ame:00121-P|Segment:null|Subtype:2-AsianAmerican|Host:Human_seqTyping_2-AsianAmerican', 'gb_KY794785_Organism_Dengue_virus_2_Strain_ame_100511_Segment_null_Subtype_2-Cosmopolitan_Host_Human_seqTyping_2-Cosmopolitan': 'gb:KY794785|Organism:Dengue_virus_2|Strain_ame:100511|Segment:null|Subtype:2-Cosmopolitan|Host:Human_seqTyping_2-Cosmopolitan'}, '/mnt/beegfs/scratch/ONEIDA/mpmachado/test_seq_typing_dengue/subtype_3.fasta': {'gb_KY586703_Organism_Dengue_virus_Strain_ame_Ser3_Thailand_Bangkok_Seq1_Segment_null_Subtype_3-II_Host_Human_seqTyping_3-II': 'gb:KY586703|Organism:Dengue_virus|Strain_ame:Ser3_Thailand_Bangkok_Seq1|Segment:null|Subtype:3-II|Host:Human_seqTyping_3-II', 'gb_KY586715_Organism_Dengue_virus_Strain_Name_Ser3_Thailand_Bangkok_Seq10_Segment_null_Subtype_3-II_Host_Human_seqTyping_3-II': 'gb:KY586715|Organism:Dengue_virus|Strain_Name:Ser3_Thailand_Bangkok_Seq10|Segment:null|Subtype:3-II|Host:Human_seqTyping_3-II'}}

    # data_by_gene = {1: {'header': 'gb_KY474330_Organism_Dengue_virus_2_Strain_ame_00104-P_Segment_null_Subtype_2-AsianAmerican_Host_Human_seqTyping_2-AsianAmerican', 'gene_coverage': 10.763786493659694, 'gene_low_coverage': 19.086757990867582, 'gene_number_positions_multiple_alleles': 0, 'gene_mean_read_coverage': 77.5123287671233, 'gene_identity': 85.02283105022832}, 2: {'header': 'gb_KY474334_Organism_Dengue_virus_2_Strain_ame_00121-P_Segment_null_Subtype_2-AsianAmerican_Host_Human_seqTyping_2-AsianAmerican', 'gene_coverage': 10.999705101739892, 'gene_low_coverage': 17.247542448614837, 'gene_number_positions_multiple_alleles': 0, 'gene_mean_read_coverage': 76.30741733690796, 'gene_identity': 84.62913315460233}, 3: {'header': 'gb_KY794785_Organism_Dengue_virus_2_Strain_ame_100511_Segment_null_Subtype_2-Cosmopolitan_Host_Human_seqTyping_2-Cosmopolitan', 'gene_coverage': 28.87053966381599, 'gene_low_coverage': 15.253660197480423, 'gene_number_positions_multiple_alleles': 13, 'gene_mean_read_coverage': 801.5263874702077, 'gene_identity': 81.68198842356145}, 4: {'header': 'gb_KY586703_Organism_Dengue_virus_Strain_ame_Ser3_Thailand_Bangkok_Seq1_Segment_null_Subtype_3-II_Host_Human_seqTyping_3-II', 'gene_coverage': 99.65753424657534, 'gene_low_coverage': 0.0, 'gene_number_positions_multiple_alleles': 1, 'gene_mean_read_coverage': 457.106529209622, 'gene_identity': 96.56357388316151}, 5: {'header': 'gb_KY586715_Organism_Dengue_virus_Strain_Name_Ser3_Thailand_Bangkok_Seq10_Segment_null_Subtype_3-II_Host_Human_seqTyping_3-II', 'gene_coverage': 100.0, 'gene_low_coverage': 0.0, 'gene_number_positions_multiple_alleles': 0, 'gene_mean_read_coverage': 7656.738151425762, 'gene_identity': 93.31366764995083}}
    # gene_coverage, gene_identity, gene_mean_read_coverage
    # gene_mean_read_coverage = 1

    # references_files = ['/mnt/beegfs/scratch/ONEIDA/mpmachado/test_seq_typing_dengue/subtype_2.fasta', '/mnt/beegfs/scratch/ONEIDA/mpmachado/test_seq_typing_dengue/subtype_3.fasta']
    pass
