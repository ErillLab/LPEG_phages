# -*- coding: utf-8 -*-
"""
RETRIEVE RELEVANT GENOMES (RRG)

[Description]

"""


from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
import json
import time
import pandas as pd
import os


config_json_file_name = "retrieve_relevant_genomes_config.json"


def set_up():
    
    global ANALYSIS_NAME
    global PROTEINS_FILENAME
    global DATABASE
    global TAXID
    global EVAL_THRESHOLD
    global QCOV_THRESHOLD
    global ATTEMPTS_LIMIT
    global EMAIL
    global API_KEY
    global SLEEP_TIME
    global MINIMUM_GENOME_LENGTH
    global ALPHA_PHAGES_FILENAME
    global PILITROPIC_PHAGES_FILENAME
    global FLAGELLOTROPIC_PHAGES_FILENAME
    global USE_TBLASTN_HITS_FROM_PREVIOUS_RUN
    
    file_reader = json.load(open(config_json_file_name))
    
    # Name of the folder for the results from this analysis
    ANALYSIS_NAME = file_reader['analysis_name']
    
    # Set file names
    PROTEINS_FILENAME = file_reader['proteins_filename']
    
    # Set BLAST parameters
    DATABASE = file_reader['blast_parameters']['database']
    TAXID = file_reader['blast_parameters']['taxid']
    EVAL_THRESHOLD = file_reader['blast_parameters']['eval_threshold']
    QCOV_THRESHOLD = file_reader['blast_parameters']['qcov_threshold']
    MINIMUM_GENOME_LENGTH = file_reader['blast_parameters']['minimum_genome_length']
    ATTEMPTS_LIMIT = file_reader['blast_parameters']['attempts_limit']
    
    # Set Entrez parameters
    EMAIL = file_reader['entrez_parameters']['email']
    API_KEY = file_reader['entrez_parameters']['api_key']
    SLEEP_TIME = file_reader['entrez_parameters']['sleep_time']
        
    # Files for alpha-infecting phages and pilitropic phages
    ALPHA_PHAGES_FILENAME = file_reader["alpha_phages_file_name"]
    PILITROPIC_PHAGES_FILENAME = file_reader["plitropic_phages_file_name"]
    FLAGELLOTROPIC_PHAGES_FILENAME = file_reader["flagellotropic_phages_file_name"]
    
    # Use the tBLASTn hits obtained in a previous run, instead of searching again
    USE_TBLASTN_HITS_FROM_PREVIOUS_RUN = file_reader["use_tblastn_hits_from_previous_run"]


def tblastn_search(proteins, tag='0'):
    
    # If download fails after as many attempts as ATTEMPTS_LIMIT, the protein will
    # be appended to this list
    failed_list = []
    
    info_table = None
    
    # taxon parameter is defined. It will be used to restrict the BLAST search
    if TAXID != None:
        taxon = "txid" + str(TAXID) + "[orgn]"
    else:
        taxon = None
    
    list_1_query_prot_acc = []
    list_2_hit_genome_acc = []
    list_3_alignment_start_pos = []
    list_4_alignment_end_pos = []
    list_5_hit_name = []
    list_6_hit_E_value = []
    list_7_hit_query_cover = []
    
    for protein in proteins:
        
        print('Blast search for protein', protein.id)
        
        failed = False
        
        for i in range(ATTEMPTS_LIMIT):
            
            try:
                # Blast search
                result_handle = NCBIWWW.qblast("tblastn", DATABASE, protein.seq, entrez_query=taxon, expect=EVAL_THRESHOLD)
                time.sleep(SLEEP_TIME)
                break
            
            except:
                print('Attempt', i + 1, 'failed')
                if i == ATTEMPTS_LIMIT - 1:
                    failed = True
        
        if failed == True:
            print('Search failed for', protein, '\nAppending it to failed_list')
            failed_list.append(protein)
            continue
        
        for n in range(ATTEMPTS_LIMIT):
            try:
                blast_record = NCBIXML.read(result_handle)
                time.sleep(SLEEP_TIME)
                break
            except:
                print('Attempt', n + 1, 'failed')
                if n == ATTEMPTS_LIMIT - 1:
                    failed = True
        
        if failed == True:
            print('Reading failed for', protein, '\nAppending it to failed_list')
            failed_list.append(protein)
            continue
        
        # Get Hits
        for alignment in blast_record.alignments:
            
            if iscompletegenome(alignment.hit_def):
                
                for hsp in alignment.hsps:
                                    
                    # Compute query cover:
                    # Length of the query sequence
                    query_length = len(protein.seq)
                    # Length of the part of the query sequence involved in the alignment
                    query_span = hsp.query_end - hsp.query_start + 1
                    # query cover
                    query_cover = query_span / query_length
                    
                    # Check whether query cover is acceptable
                    if query_cover > QCOV_THRESHOLD:               
                        # For info Table
                        list_1_query_prot_acc.append(protein.id)
                        list_2_hit_genome_acc.append(alignment.accession)
                        list_3_alignment_start_pos.append(hsp.sbjct_start)
                        list_4_alignment_end_pos.append(hsp.sbjct_end)
                        list_5_hit_name.append(alignment.hit_def)
                        list_6_hit_E_value.append(hsp.expect)
                        list_7_hit_query_cover.append(query_cover)
        
        # NEW TMP INFO TABLE
        # Better names
        list_5_hit_name_short = [s.split(', complete')[0] for s in list_5_hit_name]
        
        # Info Table with E-values and query cover values
        info_table = pd.DataFrame({'Delta phage protein': list_1_query_prot_acc,
                                   'Hit genome accession': list_2_hit_genome_acc,
                                   'Start position of the alignment': list_3_alignment_start_pos,
                                   'End position of the alignment': list_4_alignment_end_pos,
                                   'Phage name': list_5_hit_name_short,
                                   'E_value': list_6_hit_E_value,
                                   'Q_cov': list_7_hit_query_cover})
        
        # Save temporary table
        tmp_table_filename = 'tmp_' + tag + '_tBLASTn_hits_Info_Table.tsv'
        RRG_out_dir = os.path.join(ANALYSIS_NAME, 'RRG_output')
        tmp_table_filepath = os.path.join(RRG_out_dir, tmp_table_filename)
        info_table.to_csv(tmp_table_filepath, sep='\t', index=False)
    
    # Return info_table and list of failed searches
    return info_table, failed_list


def isrefseq(genome_accession):
    ''' Returns True if the input is a RefSeq genome accession '''
    return genome_accession[2] == '_'


def iscompletegenome(description):
    return 'complete genome' in description and 'partial' not in description


def check_dir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)


def organize_genomes_by_taxid(genomes):
    
    # Genome accessions found for each species
    taxid_to_accessions = {}  # key: species taxid -> value: list of accessions
    taxid_to_title = {}  # key: species taxid -> value: species title
    failed_downloads = []
    failed_reads = []
    short_genomes = []  # Genomes that are descarted because too short
    incomplete_genomes = []  # Genomes that are descarted because incomplete
    
    # Keep track of which species have a refseq
    #species_with_refseq = []
    
    for acc in genomes:
        
        failed = False  # Becomes True when either downloading or reading fails
        
        # Try to download genome info
        for n in range(ATTEMPTS_LIMIT):
            try:
                handle = Entrez.esummary(db="nucleotide", id=acc)
                time.sleep(SLEEP_TIME)                
                break
            except:
                print('Attempt', n + 1, 'failed')
                if n == ATTEMPTS_LIMIT - 1:
                    failed = True
          
        if failed:
            print('Download failed for', acc)
            print('Appending it to failed_downloads')
            failed_downloads.append(acc)
            continue
        
        # If download didn't fail, try to read handle:
        for n in range(ATTEMPTS_LIMIT):
            try:
                record = Entrez.read(handle)
                handle.close()
                break
            except:
                print('Reading attempt', n + 1, 'failed')
                time.sleep(SLEEP_TIME)
                if n == ATTEMPTS_LIMIT - 1:
                    handle.close()
                    failed = True
        
        if failed:
            print('Reading failed for', acc)
            print('Appending it to failed_reads')
            failed_reads.append(acc)
            continue
        
        # Process record
        genome_title = record[0]['Title']
        genome_taxid = int(record[0]['TaxId'])
        genome_length = int(record[0]['Length'])
        print(genome_title)
        
        # Check genome length
        if genome_length > MINIMUM_GENOME_LENGTH:
            
            # Check genome completeness
            if iscompletegenome(genome_title):
                
                # Compile species-specific list of genome accessions
                if genome_taxid in taxid_to_accessions.keys():
                    taxid_to_accessions[genome_taxid].append(acc)
                else:
                    taxid_to_accessions[genome_taxid] = [acc]
                    taxid_to_title[genome_taxid] = genome_title
            
            # Discard if incomplete
            else:
                incomplete_genomes.append(acc)
        
        # Discard if too short
        else:
            short_genomes.append(acc)
    
    # Check failed downloads
    if len(failed_downloads) > 0:
        print('\nCheck failed_downloads list: {} failed downloads'.format(
            len(failed_downloads))
        )
    
    # Check failed reads
    if len(failed_reads) > 0:
        print('\nCheck failed_reads list: {} failed reads'.format(
            len(failed_reads))
        )
    
    # Check short genomes
    if len(short_genomes) > 0:
        print('\n' + str(len(short_genomes)) + ' genomes were discarted ' +
              'because shorter than the reqired minimum length of ' +
              str(MINIMUM_GENOME_LENGTH) + ' base-pairs')
    
    # Check incomplete genomes
    if len(incomplete_genomes) > 0:
        print('\n{} genomes were discarted because incomplete'.format(
            len(incomplete_genomes))
        )
    
    return [taxid_to_accessions, taxid_to_title, failed_downloads,
            failed_reads, short_genomes, incomplete_genomes]


def loop_over_failed_searches(info_table, failed_list):
    ''' Re-runs the tblastn_search function until failed_list stops shrinking. ''' 
    
    # While loop to try again failed searches
    failed_list_old = []
    c = 0
    while failed_list != failed_list_old:
        print('tblastn_search iteration {}:\t{} failed tBLASTn queries'.format(
            c, len(failed_list))
        )
        failed_list_old = failed_list
        c += 1
        info_table_to_add, failed_list = tblastn_search(failed_list_old, str(c))
        # Add new hits to info_table
        if isinstance(info_table_to_add, pd.DataFrame):
            if len(info_table_to_add) != 0:
                info_table = pd.concat([info_table, info_table_to_add])
    
    # Clean up temporary tBLASTn tables
    RRG_out_dir = os.path.join(ANALYSIS_NAME, 'RRG_output')
    analysis_dir_content = os.listdir(RRG_out_dir)
    for f in analysis_dir_content:
        if f[:4] == 'tmp_' and f[-15:] == '_Info_Table.tsv':
            print('\nRemoving file ' + f)
            filepath = os.path.join(RRG_out_dir, f)
            os.remove(filepath)
    
    return info_table, failed_list


def select_one_genome_per_taxid(taxid_to_accessions):
    # Select accessions for the analysis
    selected_acc_dict = {}
    for taxid in taxid_to_accessions.keys():
        accessions = taxid_to_accessions[taxid]
        for acc in accessions:
            if isrefseq(acc):
                # if a refseq is present, it becomes the only element of the
                # list. Otherwise, 'accession' is going to be a list of one or
                # more genbank accesssions.
                accessions = [acc]
        # The first element (alphabetically) is going to be selected. (It's
        # going to be the refseq if 'accessions' only contained the refseq).
        accessions.sort()
        selected_acc = accessions[0]
        selected_acc_dict[taxid] = selected_acc
    return selected_acc_dict


def main():
    
    print('\nRUN NAME:', ANALYSIS_NAME, '\n')
    
    # Make directory for this run if it doesn't exist yet
    check_dir(ANALYSIS_NAME)
    # Output subfolder
    RRG_out_dir = os.path.join(ANALYSIS_NAME, 'RRG_output')
    check_dir(RRG_out_dir)
    
    # Set email and API key
    Entrez.email = EMAIL
    Entrez.apikey = API_KEY
    
    # tBLASTn
    
    if USE_TBLASTN_HITS_FROM_PREVIOUS_RUN:
        # Previous run RRG output directory
        prev_outdir = os.path.join(USE_TBLASTN_HITS_FROM_PREVIOUS_RUN, 'RRG_output')
        
        # Load tBLASTn hits from previous run
        filepath = os.path.join(prev_outdir, 'tBLASTn_hits.json')
        tblastn_accessions = json.load(open(filepath))
        
        # Load tBLASTn hits Info Table
        filepath = os.path.join(prev_outdir, 'tBLASTn_hits_Info_Table.tsv')
        info_table = pd.read_csv(filepath, sep='\t')
        
    else:
        # Perform tBLASTn search
        
        # Query proteins
        proteins_filepath = os.path.join('input_datasets', PROTEINS_FILENAME)
        proteins = SeqIO.parse(proteins_filepath, format="fasta")
        
        # Run tBLASTn
        info_table, failed_list = tblastn_search(proteins)
        
        # Try again with the failed searches until failed_list stops shrinking
        info_table, failed_list = loop_over_failed_searches(info_table, failed_list)
        
        # If there are still failed searches after the while loop, report them
        if len(failed_list) > 0:
            failed_filepath = os.path.join(RRG_out_dir, 'failed_tBLASTn_searches.json')            
            print("{} genomes couldn't be retrieved. Saving them into: {}".format(
                len(failed_list), failed_filepath))
            with open(failed_filepath, 'w') as outfile:
                json.dump(failed_list, outfile)
        
        # Hit genomes
        tblastn_accessions = list(info_table['Hit genome accession'])
        tblastn_accessions = list(set(tblastn_accessions))  # Remove duplicates
    
    # Save tBLASTn hit genomes
    filepath = os.path.join(RRG_out_dir, 'tBLASTn_hits.json')
    with open(filepath, 'w') as outfile:
        json.dump(tblastn_accessions, outfile)
    
    # Save a table about tBLASTn output with complete info
    filepath = os.path.join(RRG_out_dir, 'tBLASTn_hits_Info_Table.tsv')
    info_table.to_csv(filepath, sep='\t', index=False)
    
    
    # Load the other datasets
    
    # Alpha-infecting phage genomes
    alpha_filepath = os.path.join('input_datasets', ALPHA_PHAGES_FILENAME)
    alpha_accessions = json.load(open(alpha_filepath))
    
    # Pilitropic phage genomes
    pili_filepath = os.path.join('input_datasets', PILITROPIC_PHAGES_FILENAME)
    pilitropic_accessions = json.load(open(pili_filepath))
    
    # Flagellotropic phage genomes
    flagel_filepath = os.path.join('input_datasets', FLAGELLOTROPIC_PHAGES_FILENAME)
    flagellotropic_accessions = json.load(open(flagel_filepath))
    
    # All the genomes
    genomes = (tblastn_accessions +
               alpha_accessions +
               pilitropic_accessions +
               flagellotropic_accessions)
    genomes = list(set(genomes))  # Remove duplicates
    
    print('\n\tOrganizing genomes by taxid:')
    (taxid_to_accessions, taxid_to_title, failed_download, failed_reading,
     short_genomes, incomplete_genomes) = organize_genomes_by_taxid(genomes)
    
    # Select one representative genome per taxid
    selected_acc_dict = select_one_genome_per_taxid(taxid_to_accessions)
    
    # accession_to_taxid dictionary
    accession_to_taxid = {}
    for taxid in taxid_to_accessions.keys():
        accessions = taxid_to_accessions[taxid]
        for acc in accessions:
            accession_to_taxid[acc] = taxid
    
    # All the accessions collected
    collected_accessions = list(accession_to_taxid.keys())
    
    acc_list = []
    taxid_list = []
    title_list = []
    tblastn_list = []
    alpha_list = []
    pili_list = []
    flagel_list = []
    used_list = []
    
    for acc in genomes:
        
        # Get taxid and title
        if acc in collected_accessions:
            taxid = accession_to_taxid[acc]
            title = taxid_to_title[taxid]
        else:
            taxid = 'ignored accession'
            title = 'ignored accession'
        
        # Define source(s)
        # Is from tblastn set
        if acc in tblastn_accessions:
            tblastn_tag = 'yes'
        else:
            tblastn_tag = 'no'
        # Is from alpha-infecting set
        if acc in alpha_accessions:
            alpha_tag = 'yes'
        else:
            alpha_tag = 'no'
        # Is from pilitropic set
        if acc in pilitropic_accessions:
            pili_tag = 'yes'
        else:
            pili_tag = 'no'
        # Is from flagellotropic set
        if acc in flagellotropic_accessions:
            flagel_tag = 'yes'
        else:
            flagel_tag = 'no'
        
        # Define used_tag
        if acc in failed_download or acc in failed_reading:
            used_tag = 'No: failed request to NCBI'
        elif acc in short_genomes:
            used_tag = 'No: short genome'
        elif acc in incomplete_genomes:
            used_tag = 'No: incomplete genome'
        elif acc not in selected_acc_dict.values():
            used_tag = 'No: replaced with ' + selected_acc_dict[taxid]
        else:
            used_tag = 'Yes'
        
        acc_list.append(acc)
        taxid_list.append(taxid)
        title_list.append(title)
        tblastn_list.append(tblastn_tag)
        alpha_list.append(alpha_tag)
        pili_list.append(pili_tag)
        flagel_list.append(flagel_tag)
        used_list.append(used_tag)
    
    report_table = pd.DataFrame({'accession': acc_list,
                                 'taxid': taxid_list,
                                 'title': title_list,
                                 'from_tblastn_search': tblastn_list,
                                 'from_alpha-infecting_set': alpha_list,
                                 'from_pilitropic_set': pili_list,
                                 'from_flagellotropic_set': flagel_list,
                                 'used': used_list})
    # Save Report Table
    filepath = os.path.join(RRG_out_dir, 'relevant_genomes_Report_Table.tsv')
    report_table.to_csv(filepath, sep='\t', index=False)
    
    # All selected accessions
    selected_df = report_table.loc[report_table['used'] == 'Yes']
    all_selected_accessions = list(selected_df['accession'])
    
    # Save list of selected accessions
    filepath = os.path.join(RRG_out_dir, 'relevant_genomes.json')
    with open(filepath, 'w') as outfile:
        json.dump(all_selected_accessions, outfile)
    print('\nRRG: done.\n')



if __name__ == "__main__":
    
    # Read config file and set all the variables
    set_up()
    
    # Main program
    main()
    
    
    












