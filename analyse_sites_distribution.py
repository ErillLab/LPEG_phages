# -*- coding: utf-8 -*-
"""
ANALYZE SITES DISTRIBUTION (ASD)

[description]

"""


import json
import numpy as np
import pandas as pd
import os
import warnings


config_json_file_name = "analyse_sites_distribution_config.json"

def set_up():
    
    global ANALYSIS_NAME
    global N_BINS
    global N_SITES_THRESHOLD
    global REMOVE_TERM_REPEAT
    global USE_DOUBLE_BINNING
    global ENTROPY_THRESHOLD
    
    file_reader = json.load(open(config_json_file_name))
    
    # Name of the folder for the results from this analysis
    ANALYSIS_NAME = file_reader['analysis_name']
    
    # Number of bins
    N_BINS = file_reader["n_bins"]
    
    # Minimum number of sites on genome required for the analysis of the
    # distribution to be computed
    N_SITES_THRESHOLD = file_reader["n_sites_threshold"]
    
    # Remove terminal repeats or not
    REMOVE_TERM_REPEAT = file_reader["remove_term_repeat"]
    
    # Use double binning method or not
    USE_DOUBLE_BINNING = file_reader["use_double_binning"]
    
    # Entropy threshold: genomes with normalized positional entropy below this
    # value are considered as "LPEG" (Low Positional Entropy Genomes)
    ENTROPY_THRESHOLD = file_reader["lpeg_threshold"]


def entropy(counts):
    counts_vector = np.array(counts)
    frequencies = counts_vector / counts_vector.sum()
    H = 0
    for p in frequencies:
        if p != 0:
            H = H - p * np.log(p)
    return H


def norm_entropy(counts):
    '''
    
    Entropy divided by the maximum entropy possible with that number of counts
    and that number of bins.
    
    Parameters
    ----------
    counts : array-like object
        Counts associated to each class.

    Returns
    -------
    rel_possible_ent : float
        Ranges from 0, when entropy is 0, to 1, when entropy is the maximum
        possible entropy. The maximum possible entropy depends on the number of
        counts and bins, and it's achieved when the counts are distributed as
        evenly as possible among the bins. Example: with 10 bins and 12 counts,
        maximum possible entropy is the entropy of the distribution where 2
        bins contain 2 counts, and 8 bins contain 1 count.

    '''
    
    counts_vector = np.array(counts)
    n_obs = counts_vector.sum()
    n_bins = len(counts_vector)
    if n_obs == 1:
        rel_possible_ent = 1
    else:
        # Compute max entropy possible with that number of obs and bins
        quotient = n_obs // n_bins
        remainder = n_obs % n_bins
        chunk_1 = np.repeat(quotient, n_bins - remainder)
        chunk_2 = np.repeat(quotient + 1, remainder)
        values = np.hstack((chunk_1, chunk_2))  # values distr as evenly as possible
        max_possible_entropy = entropy(values)
        # Compute relative entropy
        rel_possible_ent = entropy(counts) / max_possible_entropy
    return rel_possible_ent
    

def gini_coeff(values_for_each_class):
    '''
    Gini coefficient measures distribution inequality.

    Parameters
    ----------
    values_for_each_class : array-like object
        Values associated to each class.
        They don't need to be already sorted and/or normalized.

    Returns
    -------
    gini_coeff : float
        Ranges from 0 (perfect equality) to 1 (maximal inequality).
    '''
    
    import numpy as np
    values = np.array(values_for_each_class)
    norm_values = values / values.sum()  # normalize
    
    # Generate Lorenz curve
    norm_values.sort()
    cum_distr = np.cumsum(norm_values)
    cum_distr = list(cum_distr)
    cum_distr.insert(0, 0)
    
    # Get area under Lorenz curve
    n_classes = len(cum_distr)-1
    under_lorenz = np.trapz(y = cum_distr, dx = 1/n_classes)
    
    # Area under Perfect Equality curve
    # It's the area of a triangle with base = 1 and height = 1
    under_PE = 0.5
    
    # Compute Gini coefficient
    gini_coeff = (under_PE - under_lorenz) / under_PE
    
    return gini_coeff


def norm_gini_coeff(values_for_each_class):
    '''
    Normalized Gini coefficient.
    The minimum and maximum possible Gini coefficient with that number of
    bins and observations are computed. Then, norm_Gini_coefficient is
    defined as
    norm_Gini_coefficient := (Gini - min_Gini) / (max_Gini - min_Gini)

    Parameters
    ----------
    values_for_each_class : array-like object
        Values associated to each class.
        They don't need to be already sorted and/or normalized.

    Returns
    -------
    norm_gini_coeff : float
        Ranges from 0 (minimal inequality possible) to 1 (maximal
        inequality possible).
    '''
    
    import numpy as np

    # Compute Gini coefficient
    nuber_of_bins = len(values_for_each_class)
    number_of_obs = np.array(values_for_each_class).sum()
    Gini = gini_coeff(values_for_each_class)
    
    # Compute minimum possible Gini coefficient
    quotient = number_of_obs // nuber_of_bins
    remainder = number_of_obs % nuber_of_bins
    chunk_1 = np.repeat(quotient, nuber_of_bins - remainder)
    chunk_2 = np.repeat(quotient + 1, remainder)
    vect = np.hstack((chunk_1, chunk_2))  # values distr as evenly as possible
    min_Gini = gini_coeff(vect)
    
    # Compute maximum possible Gini coefficient
    chunk_1 = np.repeat(0, nuber_of_bins - 1)
    chunk_2 = np.repeat(number_of_obs, 1)
    vect = np.hstack((chunk_1, chunk_2))  # values distr as unevenly as possible
    vect = [int(v) for v in vect]
    max_Gini = gini_coeff(vect)
    
    # Compute normalized Gini coefficient
    if max_Gini - min_Gini == 0:
        norm_gini = 0
    else:
        norm_gini = (Gini - min_Gini) / (max_Gini - min_Gini)
    
    return norm_gini


def study_sites_distribution(matches_list, genome_length):
    
    # Number of sites
    n_sites = len(matches_list)
    
    # Site density (sites per thousand bp)
    site_density = 1000 * n_sites / genome_length
    
    # Counts in each bin
    counts, bins = np.histogram(
        matches_list, bins=N_BINS, range = (0, genome_length)
    )
    
    entr = entropy(counts)  # Positional entropy
    norm_entr = norm_entropy(counts)  # Normalized positional entropy
    gini = gini_coeff(counts)  # Gini coefficient
    norm_gini = norm_gini_coeff(counts)  # Normalized Gini coefficient
    
    # Shifting bins by half the size of a bin -> changing the coordinate system
    half_bin_size = int((bins[1] - bins[0])/2)
    # Change coordinates (the start point moved from 0 to half_bin_size)
    shifted_matches = [m_pos - half_bin_size for m_pos in matches_list]
    for i in range(len(shifted_matches)):
        if shifted_matches[i] < 0:
            shifted_matches[i] = genome_length + shifted_matches[i]
    shifted_matches.sort()   
    
    # Counts in each shifted bin
    counts_sh, bins_sh = np.histogram(
        shifted_matches, bins=N_BINS, range = (0, genome_length)
    )
    
    if USE_DOUBLE_BINNING:
        entr_sh = entropy(counts_sh)
        norm_entr_sh = norm_entropy(counts_sh)
        gini_sh = gini_coeff(counts_sh)
        norm_gini_sh = norm_gini_coeff(counts_sh)
        
        entr = min(entr, entr_sh)
        norm_entr = min(norm_entr, norm_entr_sh)
        gini = max(gini, gini_sh)
        norm_gini = max(norm_gini, norm_gini_sh)
    
    distr_metrics = {
        "counts": counts,  # they will be used for the plot
        "site_density": site_density,
        "entropy": entr,
        "norm_entropy": norm_entr,
        "gini": gini,
        "norm_gini": norm_gini
    }
    
    return distr_metrics


def organize_info_from_MMS(input_filepath):
    
    length_dict =  {}
    description_dict = {}
    term_repeat_dict = {}
    matches_dict = {}
    scores_dict = {}
    inter_intra_genic_dict = {}
    
    chunksize = 10 ** 4  # To avoid loading into memory all the file at once
    
    start = 0
    
    for chunk in pd.read_csv(input_filepath, sep='\t', chunksize=chunksize):
        
        stop = start + chunk.shape[0]
            
        for i in range(start, stop):
            
            genome_id = chunk['Genome ID'][i]
            description = chunk['Description'][i]
            genome_length = chunk['Genome length'][i]
            dtr_coordinates = chunk['DTR coordinates'][i]
            match_pos = chunk['Match position'][i]
            match_score = chunk['Match score'][i]
            inter_intra = chunk['Intragenic'][i]
            
            if genome_id not in length_dict.keys():
                length_dict[genome_id] = int(genome_length)
            
            if genome_id not in description_dict.keys():
                description_dict[genome_id] = description
            
            if genome_id not in term_repeat_dict.keys():
                term_repeat_dict[genome_id] = dtr_coordinates
            
            if genome_id not in matches_dict.keys():
                matches_dict[genome_id] = [match_pos]
            else:
                matches_dict[genome_id].append(match_pos)
            
            if genome_id not in scores_dict.keys():
                scores_dict[genome_id] = [match_score]
            else:
                scores_dict[genome_id].append(match_score)         
            
            if genome_id not in inter_intra_genic_dict.keys():
                inter_intra_genic_dict[genome_id] = [inter_intra]
            else:
                inter_intra_genic_dict[genome_id].append(inter_intra)    
        
        start += chunksize
    
    genomes_list = list(matches_dict.keys())
    
    info = {
        'length_dict': length_dict,
        'description_dict': description_dict,
        'term_repeat_dict': term_repeat_dict,
        'matches_dict': matches_dict,
        'scores_dict': scores_dict,
        'inter_intra_genic_dict': inter_intra_genic_dict
    }
    
    return genomes_list, info


def get_avg_score(scores):
    
    n_sites = len(scores)
    if n_sites == 0:
        average_score = 'No sites'
    else:
        average_score = sum(scores)/n_sites
    
    return average_score


def get_intergenic_freq(locations):
    
    n_intergenic = locations.count('intergenic')
    n_intragenic = locations.count('intragenic')
    
    if n_intergenic + n_intragenic == 0:
        intergenic_freq = 'no_gene_was_annotated'
    else:
        intergenic_freq = n_intergenic / (n_intergenic + n_intragenic)
    
    return intergenic_freq


def dtr_freq_for_intergenic(gen_context_inside_DTR_or_GAP, gen_context_outside_DTR_or_GAP):
    '''
    Returns the frequency of intergenic sites that lie inside a DTR/gap region.
    
    gen_context_inside_DTR_or_GAP : list. Each element is a string representing
    a site inside the DTR/gap region. The string can be either "intergenic" or
    "intragenic".
    
    gen_context_outside_DTR_or_GAP : list. Each element is a string representing
    a site outside the DTR/gap region. The string can be either "intergenic" or
    "intragenic".
    '''
    n_intergenic_inside = gen_context_inside_DTR_or_GAP.count('intergenic')
    n_intergenic_outside = gen_context_outside_DTR_or_GAP.count('intergenic')
    if n_intergenic_inside + n_intergenic_outside == 0:
        frequency = 'no_intergenic_gene_found'
    else:
        frequency = n_intergenic_inside / (n_intergenic_inside + n_intergenic_outside)
    return frequency


def is_in_region(position, region_left_end, region_right_end):
    '''
    [docstring here]  # !!!
    '''
    if region_left_end < region_right_end:
        if region_left_end < position and position < region_right_end:
            return True
        else:
            return False
    
    elif region_left_end > region_right_end:
        if position < region_right_end or region_left_end < position:
            return True
        else:
            return False
    
    else:
        raise Exception('Left end and Right end of the region cannot coincide.')


def check_dir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)


def main():
    
    # Output directory
    ASD_output_dir = os.path.join(ANALYSIS_NAME, 'ASD_output')
    check_dir(ASD_output_dir)  # Make directory if it doesn't exist yet
    
    # Get data from the output of "motif_matching_scoring.py"
    input_filepath = os.path.join(ANALYSIS_NAME, 'MMS_output', 'MMS_results.tsv')
    genomes_list, info = organize_info_from_MMS(input_filepath)
    
    description_dict = info['description_dict']
    length_dict = info['length_dict']
    matches_dict = info['matches_dict']
    scores_dict = info['scores_dict']
    inter_intra_genic_dict = info['inter_intra_genic_dict']
    term_repeat_dict = info['term_repeat_dict']
    
    
    # Import datasets to define categories
    RRG_report_table_filepath = os.path.join(
        ANALYSIS_NAME, 'RRG_output', 'relevant_genomes_Report_Table.tsv')
    RRG_df = pd.read_csv(RRG_report_table_filepath, sep='\t')
    
    pilitropic_genomes = list(RRG_df.loc[RRG_df['from_pilitropic_set'] == 'yes', 'accession'])
    flagellotropic_genomes = list(RRG_df.loc[RRG_df['from_flagellotropic_set'] == 'yes', 'accession'])
    alpha_infecting_genomes = list(RRG_df.loc[RRG_df['from_alpha-infecting_set'] == 'yes', 'accession'])
    
    
    #---------------------------------------------------------------------
    # COMPUTE STATISTICS FOR THE PREDICTED TF-BINDING SITES IN EACH GENOME
    #---------------------------------------------------------------------
    
    # The following genomes have a GAP region to be studied. They are:
    # Delta, Lacusarx, MMDA13
    genomes_with_gap_region = ['MN862068', 'NC_041927', 'MN820898']
    
    
    genome_length_list = []
    description_list = []
    n_sites_list = []
    site_density_list = []
    entropy_list = []
    norm_entropy_list = []
    gini_coeff_list = []
    norm_gini_coeff_list = []
    
    average_score_list = []
    avg_score_inside_DTR_or_GAP_list = []
    avg_score_outside_DTR_list = []
    
    intergenic_freq_list = []
    intergenic_freq_inside_DTR_or_GAP_list = []
    intergenic_freq_outside_DTR_or_GAP_list = []
    
    dtr_or_gap_freq_for_intergenic_sites_list = []
    
    category_list = []
    pilus_list = []
    flagellum_list = []
    pe_category_list = []
    
    
    # For superimposed distributions plot
    for_lineplot_accession_list = []
    for_lineplot_short_name_list = []
    for_lineplot_category_list = []
    for_lineplot_y_values_list = []
    for_lineplot_selection_criterion_list = []
    
    for genome in genomes_list:
                
        description = description_dict[genome]        
        genome_length = length_dict[genome]
        matches_list = matches_dict[genome]  # Positions of hits
        scores_list = scores_dict[genome]  # Scores of the hits
        inter_intra_list = inter_intra_genic_dict[genome]  # Locations of the hits        
        dtr_list = term_repeat_dict[genome]  # Coordinates of all annotated DTR copies
        
        # Decode string representation to list of lists
        if dtr_list != 'No DTR found':
            dtr_list = json.loads(dtr_list)  # string to list of lists
        
        # Manually add DTR coordinates for 'Mesorhizobium phage Cp1R7A-A1' and
        # 'Caulobacter phage CcrColossus': the DTR was not annotated in the
        # genbank files for these two phages
        if description == 'Mesorhizobium phage Cp1R7A-A1, complete genome':
            dtr_list = [[0, 12081]]
        if description == 'Caulobacter phage CcrColossus, complete genome':
            dtr_list = [[0,16700]]
        
        # If asked, when both the copies of the DTR are reported, one of them
        # is removed (to avoid counting the sites in the DTR twice).
        # The copy at the end of the genome is the one that will be considered.
        if REMOVE_TERM_REPEAT:
            if dtr_list != 'No DTR found':
                if len(dtr_list) > 1:
                    # Remove first copy of the DTR
                    first_copy = dtr_list[0]
                    first_copy_start, first_copy_end = first_copy
                    
                    # Check that the first DTR copy is at the beginning of the
                    # genome as expected
                    if first_copy_start not in [0,1]:
                        warnings.warn("\nFirst DTR copy starts at position " +
                                      str(first_copy_start) + " for phage " +
                                      description + ".\n")
                    
                    # Remove first copy of the DTR
                    dtr_list = dtr_list[1:]
                    
                    # Update Genome length
                    genome_length = genome_length - first_copy_end
                    
                    # Update matches positions
                    matches_list = [p - first_copy_end for p in matches_list]
                    
                    # Filter out matches that were in the first copy of the DTR
                    new_matches_list = []
                    new_scores_list = []
                    new_inter_intra_list = []
                    
                    for i in range(len(matches_list)):
                        if matches_list[i] > 0:
                            new_matches_list.append(matches_list[i])
                            new_scores_list.append(scores_list[i])
                            new_inter_intra_list.append(inter_intra_list[i])
                    
                    matches_list = new_matches_list
                    scores_list = new_scores_list
                    inter_intra_list = new_inter_intra_list
                    
                    # Update DTRs coordinates
                    for dtr in dtr_list:
                        dtr[0] -= first_copy_end  # Update left_end
                        dtr[1] -= first_copy_end  # Update right_end
        
        
        """
        Study Direct Terminal Repeat (DTR) region or gap region:
          - Average score: inside VS outside DTR/GAP
          - Proportion of intergenic sites: inside VS outside DTR/GAP
        """
        
        # If this genome doesn't have DTR nor a GAP region
        if (dtr_list == 'No DTR found' and
            genome not in genomes_with_gap_region):
            
            avg_score_inside_DTR_or_GAP = 'No DTR/GAP region'
            avg_score_outside_DTR_or_GAP = 'No DTR/GAP region'
            intergenic_freq_inside_DTR_or_GAP = 'No DTR/GAP region'
            intergenic_freq_outside_DTR_or_GAP = 'No DTR/GAP region'
            dtr_or_gap_freq_for_intergenic_sites = 'No DTR/GAP region'
            
        else:
            # Compare inside VS outside the DTR/GAP region
            
            # Each "score" is a motif-matching score
            # Scores inside and outside DTR/GAP
            scores_inside_DTR_or_GAP = []
            scores_outside_DTR_or_GAP = []
            
            # Each "genetic context" can be either "intergenic" or "intragenic"
            # Genetic contexts inside and outside DTR/GAP
            gen_context_inside_DTR_or_GAP = []
            gen_context_outside_DTR_or_GAP = []
            
            # Define the genomic region to be studied (a DTR or a gap region)
            
            # DELTA
            if genome == 'MN862068':
                region_left_end, region_right_end = 86518, 606
            # LACUSARX
            elif genome == 'NC_041927':
                region_left_end, region_right_end = 119462, 130138
            # MMDA13
            elif genome == 'MN820898':
                region_left_end, region_right_end = 47759, 48946
            # ALL THE PHAGES WITH A DTR
            else:
                region_left_end, region_right_end = dtr_list[-1]
            
            # Perform DTR/GAP region analysis
            
            for i in range(len(matches_list)):
                if is_in_region(matches_list[i], region_left_end, region_right_end):
                    scores_inside_DTR_or_GAP.append(scores_list[i])
                    gen_context_inside_DTR_or_GAP.append(inter_intra_list[i])
                else:
                    scores_outside_DTR_or_GAP.append(scores_list[i])
                    gen_context_outside_DTR_or_GAP.append(inter_intra_list[i])
            
            #For new violins (alternative to line plot)
            #all_putative_sites[description] = [scores_inside_DTR_or_GAP, scores_outside_DTR_or_GAP]
            #all_putative_sites_context[description] = [gen_context_inside_DTR_or_GAP, gen_context_outside_DTR_or_GAP]
            
            # Average score inside and outside DTR/GAP region
            avg_score_inside_DTR_or_GAP = get_avg_score(scores_inside_DTR_or_GAP)
            avg_score_outside_DTR_or_GAP = get_avg_score(scores_outside_DTR_or_GAP)
            # Intergenic frequency inside and outside DTR/GAP region
            intergenic_freq_inside_DTR_or_GAP = get_intergenic_freq(gen_context_inside_DTR_or_GAP)
            intergenic_freq_outside_DTR_or_GAP = get_intergenic_freq(gen_context_outside_DTR_or_GAP)
            
            # DTR frequency for intergenic sites
            dtr_or_gap_freq_for_intergenic_sites = dtr_freq_for_intergenic(
                gen_context_inside_DTR_or_GAP, gen_context_outside_DTR_or_GAP)
        
        
        # Compute site density, entropy and Gini
        distr_metrics = study_sites_distribution(matches_list, genome_length)
        
        counts = distr_metrics['counts']
        site_density = distr_metrics['site_density']
        entr = distr_metrics['entropy']
        norm_entr = distr_metrics['norm_entropy']
        gini = distr_metrics['gini']
        norm_gini = distr_metrics['norm_gini']
        
        # Number of sites
        n_sites = len(matches_list)
        
        # Average score
        average_score = get_avg_score(scores_list)
        
        # Proportion of intergenic sites
        intergenic_freq = get_intergenic_freq(inter_intra_list)
        
        
        # Define category (later used for plotting with different colors)
        if genome == 'MN862068':
            category = 'vB_BsubS-Delta'
        elif 'Caulobacter' in description:
            category = 'Caulobacter phages'
        elif genome in alpha_infecting_genomes:
            category = 'Alphaproteobacteria phages'
        else:
            category = 'Other phages'
        
        # Define pilus_tag (later used for plotting with different shapes)
        if genome in pilitropic_genomes:
            pilus_tag = 'pilitropic'
        else:
            pilus_tag = 'unknown'
        
        # Define flagellum tag (later used for plotting with different shapes)
        if genome in flagellotropic_genomes:
            flagellum_tag = 'flagellotropic'
        else:
            flagellum_tag = 'unknown'
        
        # Define PE category (Positional Entropy category)
        if norm_entr < ENTROPY_THRESHOLD:
            pe_category = 'LPEG'
        else:
            pe_category = 'Others'
        
        # Select phages for superimposed distr plot
        selected_for_superimposed_distr_plot = False
        if n_sites >= N_SITES_THRESHOLD:
            if (category == 'Caulobacter phages' and
                (pilus_tag == 'pilitropic' or flagellum_tag == 'flagellotropic')):
                selected_for_superimposed_distr_plot = True
                selected_because = 'piltropic/flagellotropic Caulobacter phage'
            elif category == 'vB_BsubS-Delta':
                selected_for_superimposed_distr_plot = True
                selected_because = 'Delta phage'
            elif norm_entr < ENTROPY_THRESHOLD:
                selected_for_superimposed_distr_plot = True
                selected_because = 'LPEG'
        
        if selected_for_superimposed_distr_plot:
            for_lineplot_accession_list.append(genome)
            for_lineplot_short_name_list.append(description.split(',')[0])
            for_lineplot_category_list.append(category)
            for_lineplot_y_values_list.append([y / sum(counts) for y in list(counts)])
            for_lineplot_selection_criterion_list.append(selected_because)
        
        
        # APPEND METRICS TO THEIR LIST
        genome_length_list.append(genome_length)
        description_list.append(description)
        n_sites_list.append(n_sites)
        site_density_list.append(site_density)
        entropy_list.append(entr)
        norm_entropy_list.append(norm_entr)
        gini_coeff_list.append(gini)
        norm_gini_coeff_list.append(norm_gini)
        average_score_list.append(average_score)
        avg_score_inside_DTR_or_GAP_list.append(avg_score_inside_DTR_or_GAP)
        avg_score_outside_DTR_list.append(avg_score_outside_DTR_or_GAP)
        intergenic_freq_list.append(intergenic_freq)
        intergenic_freq_inside_DTR_or_GAP_list.append(intergenic_freq_inside_DTR_or_GAP)
        intergenic_freq_outside_DTR_or_GAP_list.append(intergenic_freq_outside_DTR_or_GAP)
        dtr_or_gap_freq_for_intergenic_sites_list.append(dtr_or_gap_freq_for_intergenic_sites)
        category_list.append(category)
        pilus_list.append(pilus_tag)
        flagellum_list.append(flagellum_tag)
        pe_category_list.append(pe_category)
    
    
    # Save data for the lineplots
    for_lineplot = pd.DataFrame({
        'accession': for_lineplot_accession_list,
        'short_name': for_lineplot_short_name_list,
        'category': for_lineplot_category_list,
        'y_values': for_lineplot_y_values_list,
        'selection_criterion': for_lineplot_selection_criterion_list})
    filepath = os.path.join(ASD_output_dir, 'table_for_lineplot.tsv')
    for_lineplot.to_csv(filepath, sep='\t', index=False)
    
    
    # Generate complete dataframe
    data = {'Genome accession': genomes_list,
            'Description': description_list,
            'Genome length' : genome_length_list,
            'Number of sites': n_sites_list,
            'Site density (sites per thousand bp)': site_density_list,
            'Positional entropy': entropy_list,
            'Normalized positional entropy': norm_entropy_list,
            'Gini coefficient': gini_coeff_list,
            'Normalized Gini coefficient': norm_gini_coeff_list,
            'Average score': average_score_list,
            'Average score inside DTR/gap': avg_score_inside_DTR_or_GAP_list,
            'Average score outside DTR/gap': avg_score_outside_DTR_list,
            'Intergenic freq': intergenic_freq_list,
            'Intergenic freq inside DTR/gap': intergenic_freq_inside_DTR_or_GAP_list,
            'Intergenic freq outside DTR/gap': intergenic_freq_outside_DTR_or_GAP_list,
            'DTR/gap frequency for intergenic sites': dtr_or_gap_freq_for_intergenic_sites_list,
            'Category': category_list,
            'Pilus tag': pilus_list,
            'Flagellum tag': flagellum_list,
            'PE category': pe_category_list}
    
    # Whole dataframe
    whole_df = pd.DataFrame(data)
    
    # Sliced dataframe (according to the threshold on the number of sites)
    sliced_df = whole_df.loc[whole_df['Number of sites'] >= N_SITES_THRESHOLD]
    
    # Save whole dataframe
    whole_df_filename = "results_table_unsliced.tsv"
    whole_df_filepath = os.path.join(ASD_output_dir, whole_df_filename)
    whole_df.to_csv(whole_df_filepath, sep='\t', index=False)
    
    # Save sliced dataframe
    sliced_df_filename = "results_table.tsv"
    sliced_df_filepath = os.path.join(ASD_output_dir, sliced_df_filename)
    sliced_df.to_csv(sliced_df_filepath, sep='\t', index=False)
    
    # Done
    print('\nASD: done.\n')
    


if __name__ == "__main__":
    
    # Read config file and set all the variables
    set_up()
    
    # Main program
    main()
    













