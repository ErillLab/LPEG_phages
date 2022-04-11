# -*- coding: utf-8 -*-
"""
MOTIF-MATCHING SCORING (MMS)

[description]

"""



import json
from Bio import motifs
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO
import numpy as np
import time
import os
from Bio import SeqFeature as SF
import shutil
import warnings


config_json_file_name = "motif_matching_scoring_config.json"


def set_up():
    
    global ANALYSIS_NAME
    global TF_NAME
    global BS_LIST
    global PSEUDOCOUNT_VALUE
    global EMAIL
    global API_KEY
    global SLEEP_TIME
    global ATTEMPTS_LIMIT
    global BS_LIST_FILENAME
    global DRAW_GENOMES_FROM_PREVIOUS_RUN
    
    file_reader = json.load(open(config_json_file_name))
    
    # Name of the folder for the results from this analysis
    ANALYSIS_NAME = file_reader['analysis_name']
    
    # TF name
    TF_NAME = file_reader['TF_name']
    
    # Binding sites
    BS_LIST_FILENAME = file_reader['binding_sites_filename']
    
    # Blast parameters
    PSEUDOCOUNT_VALUE = file_reader['pseudocount_value']
    
    # Entrez parameters
    EMAIL = file_reader['entrez_parameters']['email']
    API_KEY = file_reader['entrez_parameters']['api_key']
    SLEEP_TIME = file_reader['entrez_parameters']['sleep_time']
    ATTEMPTS_LIMIT = file_reader['entrez_parameters']['attempts_limit']
    
    # Draw genomes from previous run
    DRAW_GENOMES_FROM_PREVIOUS_RUN = file_reader['draw_genomes_from_previous_run']


def check_dir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)


def gene_to_site_distance(g_left, g_right, g_strand, site_pos, site_length,
                          circular_genome = False, genome_length = None):
    '''
    

    Parameters
    ----------
    g_left : INT
        Gene left-end position.
    g_right : INT
        Gene right-end position.
    g_strand : INT or STR
        Gene strand.
    site_pos : INT
        Position of the TFBS (meant to be the left-end position of the TFBS).
    site_length : INT
        Length of the TFBS.
    circular_genome : BOOL, optional
        True when the genome is circular. Used to take into account the cases
        where a gene and a TFBS are at the two extremities of the genome. The
        default is False.
    genome_length : INT, optional
        Needed to calculate distances on circular genomes. The default is None.

    Returns
    -------
    result : DICT
        The function returns a dictionary with two keys: 'distance' and 'location'.
        The 'distance' key provides an integer. Its absolute value is the distance
        between the gene start position and the closest edge of the TFBS. It's
        going to be 0 if the gene start is contained into the TFBS; it's going
        to be negative for TFBS that are upstream of the gene start; it's going
        to be positive for TFBS that are downstream of the gene start.
        The 'location' key is redundant, and it provides a string indicating
        whether the TFBS is located upstream or downstream of the gene start.
        
        
        EXAMPLE 1:
            A 23 bp TFBS located at position 1000, will be reported to be +3 bp
            from a gene start at position 997, -5 bp from a gene start at
            position 1028, and 0 bp from a gene start at position 1016 (because
            the TFBS would contain the gene start).
        
        EXAMPLE 2:
            In a circular genome of 1,000,000 bp a TFBS located at position
            1000 would be reported to be at +1030 bp from a gene start located
            at position 999,970.
        
    '''
    
    # Define site center "position" (it can be non-integer)
    edge_to_center = (site_length - 1)/2
    site_center = site_pos + edge_to_center
    
    # List of possible pairs of coordinates
    # Only one pair is necessary for non circular genomes
    coordinates = [[g_left, g_right]]
    
    # Three pairs of coordinates (in total) for circular genomes
    if circular_genome == True:
        # second pair (firt + genome_length)
        coordinates.append([g_left + genome_length, g_right + genome_length])
        # third pair (firt - genome_length)
        coordinates.append([g_left - genome_length, g_right - genome_length])
    
    # In this list, a single distance will be appended for non circular genomes,
    # while for circular genomes three distances will be recorded (for the
    # three coordinate systems)
    tmp_distances = []
    
    for coord in coordinates:
        gene_left = coord[0]
        gene_right = coord[1]
        
        # If gene is on forward strand
        if g_strand == 1 or g_strand == '+' or g_strand == '1':
            gene_start = gene_left

            tmp_distance = site_center - gene_start
        
        # If gene is on reverse strand
        elif g_strand == -1 or g_strand == '-' or g_strand == '-1':
            gene_start = gene_right

            tmp_distance = gene_start - site_center
        
        else:
            print('problem ...')
        
        tmp_distances.append(tmp_distance)
    
    tmp_absolute_distances = [abs(x) for x in tmp_distances]
    # Choose the distance with the lowest absolute value
    gene_to_siteCenter = tmp_distances[np.argmin(tmp_absolute_distances)]
   
    # Define distance
    if abs(gene_to_siteCenter) < edge_to_center:
        distance = 0
    else:
        # Reduce the absoulte value of the distance by  edge_to_center
        if gene_to_siteCenter > 0:
            gene_to_siteEdge = gene_to_siteCenter - edge_to_center
        else:
            gene_to_siteEdge = gene_to_siteCenter + edge_to_center
        distance = gene_to_siteEdge
    
    if distance < 0:
        location = 'upstream'
    elif distance == 0:
        location = 'overlapping'
    else:
        location = 'downstream'
    
    result = {'distance': int(distance), 'location': location}
    return result


def is_intragenic(g_left, g_right, site_pos, site_length):
    '''
    

    Parameters
    ----------
    g_left : INT
        Gene left-end position.
    g_right : INT
        Gene right-end position.
    site_pos : INT
        Position of the TFBS (meant to be the left-end position of the TFBS).
    site_length : INT
        Length of the TFBS.

    Returns
    -------
    intragenic : BOOL

    '''
    
    TFBS_left_end = site_pos
    TFBS_right_end = site_pos + site_length
    
    intragenic = False
    
    # If the TFBS is intragenic at least one of the two ends (left or rigth or both)
    # needs to be inside the gene
    
    if g_left <= TFBS_left_end and TFBS_left_end <= g_right:
        intragenic = True
    
    if g_left <= TFBS_right_end and TFBS_right_end <= g_right:
        intragenic = True
    
    return intragenic


def set_threshold(pssm):
    '''
    The threshold is defined so that it satisfyies (roughly) the equality
    between the −log of the false-positive rate and the information content, as
    proposed in
    
    Hertz GZ, Stormo GD. Identifying DNA and protein patterns with
    statistically significant alignments of multiple sequences.
    Bioinformatics. 1999 Jul-Aug;15(7-8):563-77.
    doi: 10.1093/bioinformatics/15.7.563. PMID: 10487864.
    
    This method is used in patser software by G. Z. Hertz, and threshold_patser
    (used here) mimcs that.
    '''
    distribution = pssm.distribution(precision=10**4)
    threshold = distribution.threshold_patser()
    return threshold


def combine_forward_and_reverse_scores(f_scores, r_scores):
    '''
    Generates "effective binding scores" by integrating the binding scores for
    the forward strand with the scores for the reversed strand.
    
    f_scores : numpy array
        Binding scores on the forward strand.
    r_scores : numpy array
        Binding scores on the reverse strand.
    effective_scores : numpy array
        Integrated binding scores.
    
    Effective binding scores, take into account contributions from both
    orientations, following the method in:
    
    Hobbs ET, Pereira T, O'Neill PK, Erill I. A Bayesian inference method for
    the analysis of transcriptional regulatory networks in metagenomic data.
    Algorithms Mol Biol. 2016 Jul 8;11:19. doi: 10.1186/s13015-016-0082-8.
    PMID: 27398089; PMCID: PMC4938975.
    '''
    effective_scores = np.log2(2**f_scores + 2**r_scores)
    return effective_scores


def download_genome(genome_id, genomes_dir):
    
    downloading_failed = False
    reading_failed = False
    
    print('downloading genome from NCBI ...')
    
    for n in range(ATTEMPTS_LIMIT):
        try:
            handle = Entrez.efetch(db="nucleotide", id=genome_id,
                                   rettype="gbwithparts", retmode="text")
            time.sleep(SLEEP_TIME)
            break
        except:
            print('Attempt', n + 1, 'failed')
            if n == ATTEMPTS_LIMIT - 1:
                downloading_failed = True
    
    if downloading_failed:
        return None, downloading_failed, reading_failed
    
    # If download didn't fail, try to read handle:
    for n in range(ATTEMPTS_LIMIT):
        try:
            genome_rec = SeqIO.read(handle, "gb")
            time.sleep(SLEEP_TIME)
            break
        except:
            print('Reading attempt', n + 1, 'failed')
            time.sleep(SLEEP_TIME)
            if n == ATTEMPTS_LIMIT - 1:
                handle.close()
                reading_failed = True
    
    if reading_failed:
        return None, downloading_failed, reading_failed
    
    # If donloading and reading were successfull, save into 'genomes' directory
    # and return genome record
    out_file_name = os.path.join(genomes_dir, genome_id + ".gb")
    SeqIO.write(genome_rec, out_file_name, "gb")
    return genome_rec, downloading_failed, reading_failed


def is_DTR(feat):
    '''
    Returns True if the input GenBank feature is a Direct Terminal Repeat.
    Returns False otherwise.
    '''
    if feat.type == 'repeat_region':
        qual = feat.qualifiers
        if 'rpt_type' in qual.keys():
            rpt_type = qual['rpt_type']
            if 'inverted' not in rpt_type and 'terminal' in rpt_type:
                return True
        elif 'note' in qual.keys():
            for note in qual['note']:
                if (note.startswith('terminal repeat region predicted') or
                    note.startswith('termnal repeat region predicted')):
                    # This is true for MT708545, MT708546 (Pasto, Paso)
                    return True
        else:
            return False


def get_dtr_coordinates(genome_rec):
    '''
    [docstring here]  # !!!
    '''
    
    coord_list = []
    for x in genome_rec.features:
        if is_DTR(x):
            coord = [int(x.location.start), int(x.location.end)]
            coord_list.append(coord)
    
    # Warning if more than 2 DTR copies are in coord_list
    if len(coord_list) > 2:
        warnings.warn("More than 2 'direct terminal repeats' found on genome " +
                      str(genome_rec.id) + ", (" + str(genome_rec.description) + ")")
    
    if len(coord_list) == 0:
        return 'No DTR found'
    else:
        return coord_list


def save_scores(effective_scores, genome_id):
    '''
    Convert the binding scores into a list of float, and save it into 'scores'
    directory as a json file.
    
    effective_scores : numpy array
    '''
    effective_scores_list = [float(x) for x in effective_scores]
    scores_dirpath = os.path.join(
        ANALYSIS_NAME, 'MMS_output', 'scores')
    check_dir(scores_dirpath)
    scores_filepath = os.path.join(scores_dirpath, genome_id + "_scores.json")
    with open(scores_filepath, 'w') as outfile:
        json.dump(effective_scores_list, outfile)


def report_failed_downloads(failed_downloads):
    
    failed_downloads_filepath = os.path.join(
        ANALYSIS_NAME, 'MMS_output', 'failed_downloads.json')
    if len(failed_downloads) > 0:
        with open(failed_downloads_filepath, 'w') as outfile:
            json.dump(failed_downloads, outfile)
        print('\n{} failed downloads.'.format(len(failed_downloads)))
        print('Check "failed_downloads.json" in the "MMS_output" directory.')
    else:
        # Remove previous file if present. If no error was encountered no files
        # with this name will be present in the ANALYSIS_NAME directory
        if os.path.exists(failed_downloads_filepath):
            os.remove(failed_downloads_filepath)


def report_failed_reads(failed_reads):
    
    failed_reads_filepath = os.path.join(
        ANALYSIS_NAME, 'MMS_output', 'failed_reads.json')
    if len(failed_reads) > 0:
        with open(failed_reads_filepath, 'w') as outfile:
            json.dump(failed_reads, outfile)
        print('\n{} failed reads.'.format(len(failed_reads)))
        print('Check "failed_reads.json" in the "MMS_output" directory.')
    else:
        # Remove previous file if present. If no error was encountered no files
        # with this name will be present in the ANALYSIS_NAME directory
        if os.path.exists(failed_reads_filepath):
            os.remove(failed_reads_filepath)


def report_genomes_with_no_matches(no_matches_genomes_list):
    
    no_matches_filepath = os.path.join(
        ANALYSIS_NAME, 'MMS_output', 'no_matches_genomes.json')
    if len(no_matches_genomes_list) > 0:        
        with open(no_matches_filepath, 'w') as outfile:
            json.dump(no_matches_genomes_list, outfile)
        print('\n{} genomes with 0 matches.'.format(
            len(no_matches_genomes_list)))
        print('Check "no_matches_genomes.json" in the "MMS_output" directory.')
    else:
        # Remove previous file if present. If no error was encountered no files
        # with this name will be present in the ANALYSIS_NAME directory
        if os.path.exists(no_matches_filepath):
            os.remove(no_matches_filepath)


def report_genomes_without_genes(no_genes_genomes_dict):
    
    no_genes_filepath = os.path.join(
        ANALYSIS_NAME, 'MMS_output', 'no_genes_genomes.json')
    if len(no_genes_genomes_dict) > 0:
        with open(no_genes_filepath, 'w') as outfile:
            json.dump(no_genes_genomes_dict, outfile)
        print('\n{} genomes without annotated genes.'.format(
            len(no_genes_genomes_dict)))
        print('Check "no_genes_genomes.json" in the "MMS_output" directory.')
    else:
        # Remove previous file if present. If no error was encountered no files
        # with this name will be present in the ANALYSIS_NAME directory
        if os.path.exists(no_genes_filepath):
            os.remove(no_genes_filepath)


def get_genetic_context(match, genome_rec, motif_len, genome_len):
    
    
    # If the match is found to be within a gene, this variable will be set to
    # "intragenic". Otherwise, it remains "intergenic". However, if the genome
    # doesn't contain any annotated genes, the variable will be eventually set
    # to "no_gene".
    inter_or_intra_genic = 'intergenic'
    
    genes = []
    genes_idx = []
    distances = []
    
    for g_idx in range(len(genome_rec.features)):
        gene = genome_rec.features[g_idx]
        
        # Ignore the feature if it's not 'CDS' nor 'tRNA' nor 'rRNA'
        if gene.type != 'CDS' and gene.type != 'tRNA' and gene.type != 'rRNA':
            continue
        
        # Get gene location
        strand = gene.location.strand
        left = gene.location.start.position
        right = gene.location.end.position
        # Where is the match relative to the gene
        relative_pos = gene_to_site_distance(
            left, right, strand, match, motif_len,
            circular_genome = True, genome_length = genome_len
        )
        # Distance (from site to gene start)
        distance = relative_pos['distance']
        # As soon as one gene is found to contain the match, the variable
        # 'inter_or_intra_genic' is set to "intragenic"
        if is_intragenic(left, right, match, motif_len):
            inter_or_intra_genic = 'intragenic'
        
        genes.append(gene)
        genes_idx.append(g_idx)
        distances.append(distance)
    
    # Now get info about the closets gene. It will be stored in this dictionary
    closest_gene_info = {}
    
    # If the genome contains annotated genes
    if len(genes) > 0:
        no_annotated_genes = False
        
        # Get closest gene record
        abs_distances = [abs(x) for x in distances]
        j = np.argmin(abs_distances)  # j-th element was the closest gene
        closest_gene = genes[j]
        
        # Get its index in genome_rec.features
        closest_gene_info['rec_idx'] = genes_idx[j]
        
        # Get its gene id
        closest_gene_info['gene_id'] = closest_gene.id
                
        # Get locus tag
        if 'locus_tag' in closest_gene.qualifiers.keys():
            closest_gene_info['locus_tag'] = closest_gene.qualifiers['locus_tag']
        else:
            closest_gene_info['locus_tag'] = 'no_locus_tag'
        
        # Get db_xref
        if 'db_xref' in closest_gene.qualifiers.keys():
            closest_gene_info['db_xref'] = closest_gene.qualifiers['db_xref']
        else:
            closest_gene_info['db_xref'] = 'no_db_xref'
        
        # Get its distance from match
        closest_gene_info['distance'] = distances[j]
        
        # Get location
        closest_gene_info['strand'] = closest_gene.location.strand
        closest_gene_info['left_end'] = closest_gene.location.start.position
        closest_gene_info['right_end'] = closest_gene.location.end.position
        
        # Get sequence type
        closest_gene_info['type'] = closest_gene.type
        
        # Get product
        if 'product' in closest_gene.qualifiers.keys():
            closest_gene_info['product'] = closest_gene.qualifiers['product']
        else:
            closest_gene_info['product'] = 'no_product'
        
        # Get protein_id
        if 'protein_id' in closest_gene.qualifiers.keys():
            closest_gene_info['prot_id'] = closest_gene.qualifiers['protein_id']
        else:
            closest_gene_info['prot_id'] = 'no_protein_id'
    
    else:
        no_annotated_genes = True
                
        inter_or_intra_genic = 'no_gene'
        
        closest_gene_info['rec_idx'] = 'no_gene'
        closest_gene_info['gene_id'] = 'no_gene'
        closest_gene_info['locus_tag'] = 'no_gene'
        closest_gene_info['db_xref'] = 'no_gene'
        closest_gene_info['distance'] = 'no_gene'
        closest_gene_info['strand'] = 'no_gene'
        closest_gene_info['left_end'] = 'no_gene'
        closest_gene_info['right_end'] = 'no_gene'
        closest_gene_info['type'] = 'no_gene'
        closest_gene_info['product'] = 'no_gene'
        closest_gene_info['prot_id'] = 'no_gene'
    
    return inter_or_intra_genic, closest_gene_info, no_annotated_genes


def report_warnings(failed_downloads, failed_reads, no_matches, no_genes):
    
    if (len(failed_downloads) + len(failed_reads) +
        len(no_matches) + len(no_genes) > 0):
        print('\n\nWARNINGS:')
    
    report_failed_downloads(failed_downloads)
    report_failed_reads(failed_reads)
    report_genomes_with_no_matches(no_matches)
    report_genomes_without_genes(no_genes)


def add_TFBS_site_feature_to_gbk(genome_record, TFname, start, end, locus_tag, score):
    
    # Start location and end location for the feature
    my_start_pos = SF.ExactPosition(start)
    my_end_pos = SF.ExactPosition(end)
    
    # Use the locations do define a FeatureLocation
    my_feature_location = SF.FeatureLocation(my_start_pos,my_end_pos)
    
    # Define a feature type
    my_feature_type = TFname + "_site"
    
    # Create a SeqFeature
    my_feature = SF.SeqFeature(my_feature_location,type=my_feature_type)
    
    # Add qualifiers
    my_feature.qualifiers['locus_tag'] = locus_tag
    my_feature.qualifiers['score'] = score
    
    # Add your newly created SeqFeature to SeqRecord
    for i in range(len(genome_record.features)):
        feat = genome_record.features[i]
        if feat.location.start > start:
            break
    new_feature_list = (genome_record.features[:i] +
                        [my_feature] +
                        genome_record.features[i:])
    genome_record.features = new_feature_list


def get_genomes_from_previous_run(genomes_acc_list):
    
    for genome_id in genomes_acc_list:
        genomefilename = genome_id + ".gb"
        # 'genomes' directory of the previous run
        previous_run_genomes_dir = os.path.join(
            DRAW_GENOMES_FROM_PREVIOUS_RUN, 'genomes')
        # 'genomes' directory of the current run
        current_run_genomes_dir = os.path.join(
            ANALYSIS_NAME, 'genomes')
        # If the genome was already downloaded in the previous run
        prev_run_genomes = os.listdir(previous_run_genomes_dir)
        if genomefilename in prev_run_genomes:
            # Copy it into the 'genomes' directory of the current run
            genomefilepath = os.path.join(previous_run_genomes_dir,
                                          genomefilename)
            shutil.copy2(genomefilepath, current_run_genomes_dir)



def main():
    
    # Output directory
    check_dir(ANALYSIS_NAME)
    MMS_output_dir = os.path.join(ANALYSIS_NAME, 'MMS_output')
    check_dir(MMS_output_dir)
    
    # Set email and API key
    Entrez.email = EMAIL
    Entrez.apikey = API_KEY
    
    # Load list of genome accessions
    genomes_filepath = os.path.join(ANALYSIS_NAME, 'RRG_output', 'relevant_genomes.json')
    genomes_acc_list = json.load(open(genomes_filepath))
    
    print('Generating motif model and setting threshold ...')
    
    # Generate motif model
    bs_list_filepath = os.path.join('input_datasets', BS_LIST_FILENAME)
    bs_list = json.load(open(bs_list_filepath))
    bs_seq_list = [Seq(x) for x in bs_list]
    m = motifs.create(bs_seq_list)
    
    # Generate PWM (and reverse complement)
    pwm = m.counts.normalize(pseudocounts=PSEUDOCOUNT_VALUE)
    rpwm = pwm.reverse_complement()
    
    # Generate PSSM (and reverse complement)
    pssm = pwm.log_odds()
    rpssm = rpwm.log_odds()
    
    # Motif length
    motif_len = pssm.length
    
    # Set score threshold
    threshold = set_threshold(pssm)
    
    # Create CSV file where results are going to be stored
    output_filepath = os.path.join(MMS_output_dir, 'MMS_results.tsv')
    f = open(output_filepath, "w")
    # Write headers
    # First 7 columns
    f.write('Genome ID\tDescription\tGenome length\tDTR coordinates\t' +
            'Match position\tMatch score\tMatch sequence\tIntragenic\t')
    # More 11 columns refer to the closest gene to the match
    f.write('Gene index\tGene ID\tGene tag\tdb_xref\tGene distance\t' +
            'Gene strand\tGene left end\tGene right end\tType\tProduct\t' +
            'Protein ID')
    f.close()
    
    # Create directories to store the genomes if they don't exist yet
    genomes_dir = os.path.join(ANALYSIS_NAME, 'genomes')
    genomes_TFBS_annotation_dir = os.path.join(
        ANALYSIS_NAME, 'genomes_' + TF_NAME + '_sites_annotated')
    check_dir(genomes_dir)
    check_dir(genomes_TFBS_annotation_dir)
    
    
    if DRAW_GENOMES_FROM_PREVIOUS_RUN:
        get_genomes_from_previous_run(genomes_acc_list)
    
    
    failed_downloads = []
    failed_reads = []
    
    # If a genome record doesn't contain gene records, info about the genome
    # and the position of the matches found on that genome will be stored in
    # this dictionary
    no_genes = {}
    
    # Genomes on which no match with the PSSM is found are stored in this list
    no_matches = []
    
    i = 1  # counter to print progress
    
    for genome_id in genomes_acc_list:
        
        print('\n\t', i, '/', len(genomes_acc_list))
        i += 1
        print('\nGENOME', genome_id)
        
        # GET GENOME SEQUENCE
        
        # Genomes already downloaded into the genomes directory
        genomes_dir_content = os.listdir(genomes_dir)
        
        # If this genome has not been downloaded yet, it's downloaded and saved
        # in 'genomes' directory        
        if genome_id + ".gb" not in genomes_dir_content:
            genome_rec, download_failed, read_failed = download_genome(
                genome_id, genomes_dir)
            
            if download_failed:
                print('Download failed for', genome_id)
                print('Appending it to failed_downloads')
                failed_downloads.append(genome_id)
                continue
            
            if read_failed == True:
                print('Reading failed for', genome_id)
                print('Appending it to failed_reads')
                failed_reads.append(genome_id)
                continue
        
        # If this genome is already present, it's read from the 'genomes'
        # directory
        else:
            print('retrieving genome from genomes directory...')
            file_name = os.path.join(genomes_dir, genome_id + ".gb")
            genome_rec = SeqIO.read(file_name, "gb")
        
        # Genome description
        description = genome_rec.description
        
        # If DTRs are annotated, get their coordinates
        dtr_coord = get_dtr_coordinates(genome_rec)
        
        
        # SCAN THE GENOME
        print('scanning genome ...')
        f_scores = pssm.calculate(genome_rec.seq)  # Scan on forward strand
        r_scores = rpssm.calculate(genome_rec.seq)  # Scan on reverse strand
        effective_scores = combine_forward_and_reverse_scores(f_scores, r_scores)
        
        # Save scores into 'scores' folder
        save_scores(effective_scores, genome_id)
        
        # Find matches
        matches = np.argwhere(effective_scores > threshold).flatten()
        matches = [int(n) for n in matches]  # convert to list of int
        
        
        ######## Annotate TFBS on gbk ###########
        for k in range(len(matches)):
            match_start = matches[k]
            score = effective_scores[match_start]
            match_end = match_start + motif_len
            locus_tag = TF_NAME + "_site_" + str(k + 1)
            add_TFBS_site_feature_to_gbk(genome_rec, TF_NAME, match_start,
                                         match_end, locus_tag, score)
        filepath = os.path.join(genomes_TFBS_annotation_dir, genome_id + ".gb")
        with open(filepath, "w") as f:
            SeqIO.write(genome_rec, f, "genbank")
        #########################################
        
        
        # Find closest gene to each match
        print('finding genes close to the matches and appending results to ' +
              'the output file ...')
        
        genome_len = len(genome_rec.seq)
        
        if len(matches) == 0:
            no_matches.append(genome_id)
        
        for match in matches:
            
            inter_or_intra_genic, closest_gene_info, no_annotated_genes = (
                get_genetic_context(match, genome_rec, motif_len, genome_len)
            )
            
            # if the genome doesn't contain annotated genes
            if no_annotated_genes:
                # compile no_genes dictionary
                if genome_id in no_genes.keys():
                    no_genes[genome_id].append(match)
                else:
                    no_genes[genome_id] = [match]
            
            match_sequence = genome_rec.seq[match:match+motif_len]
            
            # Append results to the output file
            f = open(output_filepath, "a")
            f.write('\n')
            f.write(genome_id)
            f.write('\t')
            f.write(description)
            f.write('\t')
            f.write(str(genome_len))
            f.write('\t')
            f.write(str(dtr_coord))  # DTR coordinates
            f.write('\t')
            f.write(str(match))  # match position
            f.write('\t')
            f.write(str(effective_scores[match]))  # match score
            f.write('\t')
            f.write(str(match_sequence))  # match sequence
            f.write('\t')
            f.write(inter_or_intra_genic)  # match is intergenic/intragenic
            f.write('\t')
            f.write(str(closest_gene_info['rec_idx']))
            f.write('\t')
            f.write(str(closest_gene_info['gene_id']))
            f.write('\t')
            f.write(str(closest_gene_info['locus_tag']))
            f.write('\t')
            f.write(str(closest_gene_info['db_xref']))
            f.write('\t')
            f.write(str(closest_gene_info['distance']))
            f.write('\t')
            f.write(str(closest_gene_info['strand']))
            f.write('\t')
            f.write(str(closest_gene_info['left_end']))
            f.write('\t')
            f.write(str(closest_gene_info['right_end']))
            f.write('\t')
            f.write(str(closest_gene_info['type']))
            f.write('\t')
            f.write(str(closest_gene_info['product']))
            f.write('\t')
            f.write(str(closest_gene_info['prot_id']))
            f.close()
    
    # Final messages to the reader
    report_warnings(failed_downloads, failed_reads, no_matches, no_genes)
    print('\nMMS: done.\n')



if __name__ == "__main__":
    
    # Read config file and set all the variables
    set_up()
    
    # Main program
    main()








