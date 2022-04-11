# -*- coding: utf-8 -*-
"""
GENERATE PLOTS



"""


import json
import os
import numpy as np
import pandas as pd
import random
import scipy.stats as stats
import matplotlib.pyplot as plt

from plotnine import (
    ggplot, geom_point, geom_violin, aes, theme, element_blank,
    scale_shape_manual, scale_fill_manual, ylim
)

from analyse_sites_distribution import norm_entropy, norm_gini_coeff


config_json_file_name = 'generate_plots_config.json'


def set_up():
    
    global ANALYSIS_NAME
    global LPEG_THRESHOLD
    global IMAGE_WIDTH
    global IMAGE_HEIGHT
    global IMAGE_RESOLUTION
    global POINT_TRANSPARENCY
    global N_RANDOM_CONTROLS_PER_N_SITES
    global N_BINS_USED
    global VIOL_HEIGHT
    global VIOL_WIDTH
    global VIOL_DPI
    
    file_reader = json.load(open(config_json_file_name))
    
    ANALYSIS_NAME = file_reader['analysis_name']
    LPEG_THRESHOLD = file_reader['lpeg_threshold']
    IMAGE_WIDTH = file_reader['image_width']
    IMAGE_HEIGHT = file_reader['image_height']
    IMAGE_RESOLUTION = file_reader['image_resolution']
    POINT_TRANSPARENCY = file_reader['point_transparency']
    
    # Number of random positional distributions to be generated for each
    # possible number of sites (to be used as statistical control).
    N_RANDOM_CONTROLS_PER_N_SITES = file_reader['n_random_controls_per_n_sites']
    #  For example,
    # if genomes have from 6 to 200 predicted binding sites, there are 195
    # possible s values (where s is the number of sites). For each s value,
    # the script generates r random positionbal distributions of s sites,
    # where r is equal to  N_RANDOM_CONTROLS_PER_N_SITES.
    
    # Number of bins used in ASD
    # ("n_bins" parameter in  analyse_sites_distribution_config.json)
    N_BINS_USED = file_reader['n_bins_used']
    
    # Violin plots parameters
    VIOL_HEIGHT = file_reader['violin_plots_parameters']['height']
    VIOL_WIDTH = file_reader['violin_plots_parameters']['width']
    VIOL_DPI = file_reader['violin_plots_parameters']['dpi']


def check_dir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)


def df_column_to_float(df, df_column):
    df[df_column] = df[df_column].astype(float)


def save_mwu(outfilename, title, mwu):
    filepath = os.path.join(ANALYSIS_NAME, 'plots', outfilename)
    f = open(filepath, "w")
    f.write(title + '\n')
    f.write(str(mwu))
    f.close()


def inside_vs_outside_prepare_df(res_df, inside_colname, outside_colname, feature_name):
    
    lpeg_df = get_lpeg_regions_df(res_df)
    
    inside_list = list(lpeg_df[inside_colname])
    outside_list = list(lpeg_df[outside_colname])
    total_list = inside_list + outside_list
    classes = ['Inside DTR/gap'] * len(inside_list) + ['Outside DTR/gap'] * len(outside_list)
    acc_names = list(lpeg_df['Genome accession']) * 2
    dict_for_violins = {'Genome accessions': acc_names,
                        feature_name: total_list,
                        'Class': classes}
    df_for_violins = pd.DataFrame(dict_for_violins)
        
    return df_for_violins, inside_list, outside_list


def make_superimposed_lineplot_top(lineplot_df):
    
    palette_1 = ['dodgerblue', 'red', 'limegreen', 'greenyellow', 'magenta',
                 'blueviolet', 'darkorange', 'pink', 'green', 'blue', 'sienna',
                 'silver']
    
    plt.figure(figsize=(15,10))
    x_values = np.arange(0, 1, 0.02)
    c = 0
    legend_for_pili_caulo_and_delta = []
    for i in range(len(lineplot_df)):
        y_values = lineplot_df.loc[i, 'y_values']
        y_values = json.loads(y_values)  # string to list of floats
        shortname = lineplot_df.loc[i, 'short_name']
        crit = lineplot_df.loc[i, 'selection_criterion']
        if crit == 'piltropic/flagellotropic Caulobacter phage':
            plt.plot(x_values, y_values, color = palette_1[c % len(palette_1)])
            c += 1
            legend_for_pili_caulo_and_delta.append(shortname)
        elif crit == 'Delta phage':
            plt.plot(x_values, y_values, color = 'black')
            legend_for_pili_caulo_and_delta.append(shortname)
    
    #plt.legend(legend_for_pili_caulo_and_delta, bbox_to_anchor=(1.04,1), loc="upper left")
    plt.legend(legend_for_pili_caulo_and_delta)
    plt.xlabel('Normalized genome position')
    plt.ylabel('Normalized counts')
    plots_dir = os.path.join(ANALYSIS_NAME, 'plots')
    filepath = os.path.join(plots_dir, 'lineplot_top.png')
    plt.savefig(filepath, bbox_inches="tight")
    plt.close()
    

def make_superimposed_lineplot_bottom(lineplot_df):
    
    plt.figure(figsize=(15,10))
    x_values = np.arange(0, 1, 0.02)
    for i in range(len(lineplot_df)):
        y_values = lineplot_df.loc[i, 'y_values']
        y_values = json.loads(y_values)  # string to list of floats
        name = lineplot_df.loc[i, 'short_name']
        ctgr = lineplot_df.loc[i, 'category']
        crit = lineplot_df.loc[i, 'selection_criterion']
        if crit == 'piltropic/flagellotropic Caulobacter phage':
            plt.plot(x_values, y_values, color = 'lightgrey')
        elif crit == 'Delta phage':
            plt.plot(x_values, y_values, color = 'lightgrey')
        elif name == 'Ruegeria phage vB_RpoS-V16':
            plt.plot(x_values, y_values, linestyle='dotted')
        elif ctgr == 'Caulobacter phages':
            plt.plot(x_values, y_values)
        elif ctgr == 'Alphaproteobacteria phages':
            plt.plot(x_values, y_values, linestyle='dashed')
        else:
            plt.plot(x_values, y_values, linestyle='dotted')
    
    legend_list = list(lineplot_df['short_name'])
    plt.legend(legend_list, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('Normalized genome position')
    plt.ylabel('Normalized counts')
    plots_dir = os.path.join(ANALYSIS_NAME, 'plots')
    filepath = os.path.join(plots_dir, "lineplot_bottom.png")
    plt.savefig(filepath, bbox_inches="tight", dpi=900)
    plt.close()


def make_superimposed_lineplot(lineplot_df):
    
    data = lineplot_df
    
    rainbow_idx = 0
    stiles = ['solid', 'dashed', 'dashdot']
    
    plt.figure(figsize=(15,10))
    x_values = np.arange(0, 1, 0.02)
    for idx in data.index:
        y_values = data.loc[idx, 'y_values']
        y_values = json.loads(y_values)  # string to list of floats
        name = data.loc[idx, 'short_name']
        if name == 'Brevundimonas phage vB_BsubS-Delta':
            plt.plot(x_values, y_values, color='black')
        else:
            plt.plot(x_values, y_values,
                     linestyle=stiles[rainbow_idx // 10], alpha=0.8)
            rainbow_idx += 1
    
    legend_list = list(data['short_name'])
    plt.legend(legend_list, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xlabel('Normalized genome position')
    plt.ylabel('Normalized counts')
    filepath = os.path.join(ANALYSIS_NAME, 'plots', "lineplot.png")
    plt.savefig(filepath, bbox_inches="tight", dpi=600)
    plt.close()


def get_lpeg_regions_df(res_df):
    
    # Slice dataframe keeping only genomes where DTR/gap region is present and
    # contains at least one site, so that the average score among the sites
    # that are inside/outside the region can be compared
    only_with_populated_DTR_df = res_df.loc[
        (res_df['Average score inside DTR/gap']!='No DTR/GAP region') &
        (res_df['Average score inside DTR/gap']!= 'No sites') &
        (res_df['Average score outside DTR/gap']!= 'No sites')
    ].copy()
    
    # Convert average scores and intergenic frequencies from strings to float
    columns_to_float = ['Average score inside DTR/gap',
                        'Average score outside DTR/gap',
                        'Intergenic freq inside DTR/gap',
                        'Intergenic freq outside DTR/gap']
    for column in columns_to_float:
        df_column_to_float(only_with_populated_DTR_df, column)
    
    # This study is performed on LPEG
    lpeg_df = only_with_populated_DTR_df.loc[only_with_populated_DTR_df['PE category']=='LPEG']
    return lpeg_df



def get_dtr_freq_for_interg_sites_df(res_df):
    # Slice dataframe keeping only genomes where DTR/gap region is present and
    # at least one intergenic site is present
    only_with_DTR_df = res_df.loc[
        (res_df['DTR/gap frequency for intergenic sites']!='No DTR/GAP region') &
        (res_df['DTR/gap frequency for intergenic sites']!= 'no_intergenic_gene_found')
    ].copy()
    
    # Convert average scores and intergenic frequencies from strings to float
    df_column_to_float(only_with_DTR_df, 'DTR/gap frequency for intergenic sites')
    
    # This study is performed on LPEG
    lpeg_df = only_with_DTR_df.loc[only_with_DTR_df['PE category']=='LPEG']
    return lpeg_df




def get_scatterplot(dataframe, y_axis_units, legend_pos, neg_control_df):
    
    # Generate a dataframe for pilitropic/flagellotropic phages
    # (it will be used to plot circles/squares)
    only_pili_df = dataframe.loc[dataframe['Pilus tag'] == 'pilitropic'].copy()
    only_flag_df = dataframe.loc[dataframe['Flagellum tag'] == 'flagellotropic'].copy()
    only_pili_df['pilus_flag_tag'] = only_pili_df['Pilus tag']
    only_flag_df['pilus_flag_tag'] = only_flag_df['Flagellum tag']
    pili_flag_df = pd.concat([only_pili_df, only_flag_df], ignore_index=True, sort=False)
    pili_flag_df = pili_flag_df.replace('flagellotropic', 'flagellotropic' + (' ' * 24))
    pili_flag_df = pili_flag_df.assign(
        pilus_flag_tag = pd.Categorical(
            pili_flag_df['pilus_flag_tag'], ['pilitropic', 'flagellotropic' + (' ' * 24)])
    )
    
    scatterplot = (
        
        ggplot()
        
        + geom_point(
            neg_control_df,
            aes(x='Number of sites', y=y_axis_units, fill='Category'),
            color='None', alpha=POINT_TRANSPARENCY)
        
        + geom_point(
            dataframe,
            aes(x='Number of sites', y=y_axis_units, fill='Category'),
            color='None', alpha=POINT_TRANSPARENCY)
        
        + scale_fill_manual(['#91db57', '#db5f57', '#57d3db', '#a157db', '#b3b3b3'])
        
        # Plot empty circles/squares for pilitropic/flagellotropic phages
        + geom_point(
            pili_flag_df,
            aes(x='Number of sites', y=y_axis_units, shape="pilus_flag_tag"),
            color='black', alpha=POINT_TRANSPARENCY, fill='None', size=2)
        
        + scale_shape_manual(['o', 's'])
        
        + theme(legend_position=legend_pos)
        + theme(legend_title=element_blank())
        + theme(legend_background = element_blank())
    
    )
    
    return scatterplot


def get_2D_plot(dataframe, x_axis_units, y_axis_units):
    
    # Generate a dataframe for pilitropic/flagellotropic phages
    # (it will be used to plot circles/squares)
    only_pili_df = dataframe.loc[dataframe['Pilus tag'] == 'pilitropic'].copy()
    only_flag_df = dataframe.loc[dataframe['Flagellum tag'] == 'flagellotropic'].copy()
    only_pili_df['pilus_flag_tag'] = only_pili_df['Pilus tag']
    only_flag_df['pilus_flag_tag'] = only_flag_df['Flagellum tag']
    pili_flag_df = pd.concat([only_pili_df, only_flag_df], ignore_index=True, sort=False)
    pili_flag_df = pili_flag_df.replace('flagellotropic', 'flagellotropic' + (' ' * 24))
    pili_flag_df = pili_flag_df.assign(
        pilus_flag_tag = pd.Categorical(
            pili_flag_df['pilus_flag_tag'], ['pilitropic', 'flagellotropic' + (' ' * 24)])
    )
    
    scatterplot = (
        
        ggplot()
        
        + geom_point(
            dataframe,
            aes(x=x_axis_units, y=y_axis_units, fill='Category'),
            color='None', alpha=POINT_TRANSPARENCY)
        
        + scale_fill_manual(['#db5f57', '#91db57', '#57d3db', '#a157db', '#b3b3b3'])
        
        # Plot empty circles/squares for pilitropic/flagellotropic phages
        + geom_point(
            pili_flag_df,
            aes(x=x_axis_units, y=y_axis_units, shape="pilus_flag_tag"),
            color='black', alpha=POINT_TRANSPARENCY, fill='None', size=2)
        
        + scale_shape_manual(['o', 's'])
                
        + theme(legend_title=element_blank())
        + theme(legend_background = element_blank())
    
    )
    
    return scatterplot


def generate_neg_control(min_n_sites, max_n_sites):
        
    n_sites_list = []
    entropy_list = []
    gini_list = []
    
    for rep in range(N_RANDOM_CONTROLS_PER_N_SITES):
        
        for n_sites in range(min_n_sites, max_n_sites+1):
            
            # New histogram
            counts = [0]*N_BINS_USED
            for i in range(n_sites):
                selected_bin = random.randint(0, N_BINS_USED-1)
                counts[selected_bin] = counts[selected_bin] + 1
            
            norm_ent = norm_entropy(counts)
            norm_gini = norm_gini_coeff(counts)
            
            # Shift bins by half bin size
            for i in range(len(counts)):
                bin_count = counts[i]
                for site in range(bin_count):
                    if random.random() < 0.5:
                        counts[i] -= 1
                        counts[i-1] += 1
            norm_ent_sh = norm_entropy(counts)
            norm_gini_sh = norm_gini_coeff(counts)
            
            norm_ent = min(norm_ent, norm_ent_sh)
            norm_gini = max(norm_gini, norm_gini_sh)
            
            n_sites_list.append(n_sites)
            entropy_list.append(norm_ent)
            gini_list.append(norm_gini)
    
    control_df = pd.DataFrame({'Number of sites': n_sites_list,
                               'Normalized positional entropy': entropy_list,
                               'Normalized Gini coefficient': gini_list})
    # Add 'Category' column
    control_df['Category'] = ['control'] * len(control_df)
    
    return control_df


def jitter(dataframe, columns, jitter_param):
    '''
    Add a random value (from a uniform distr) to generate jitter.
    '''
    
    
    jitter_df = dataframe.copy()
    
    for col in columns:
        val_range = jitter_df[col].max() - jitter_df[col].min()
        max_shift = abs(val_range/(jitter_param * IMAGE_WIDTH))
        
        for index in jitter_df.index:
            shift = random.uniform(-max_shift, max_shift)
            jitter_df.loc[index, col] = jitter_df.loc[index, col] + shift
    
    return jitter_df


def main():
    
    plots_dir = os.path.join(ANALYSIS_NAME, 'plots')
    check_dir(plots_dir)
    
    # Load data from ASD (Analyze Site Distribution)
    
    # Load main results table
    filepath = os.path.join(ANALYSIS_NAME, 'ASD_output', 'results_table.tsv')
    res_df = pd.read_csv(filepath, sep='\t')
    res_df.index = res_df['Genome accession']
    # Load data for lineplots
    filepath = os.path.join(ANALYSIS_NAME, 'ASD_output', 'table_for_lineplot.tsv')
    lineplot_df = pd.read_csv(filepath, sep='\t')
    
    
    # LINE PLOTS (to show positional distributions)
    # ---------------------------------------------
    
    print('Generating line plot ...')
    # Line plot
    make_superimposed_lineplot(lineplot_df)
    
    
    # VIOLINS (Average score and intergenic frequency)
    # ------------------------------------------------
    
    print('Generating violin plots + Mann Whitney U tests ...')
    # Prepare dataframe for DTR avg score comparison (inside VS outside)
    in_VS_out_dtr_avg_scr_df, inside_list, outside_list = inside_vs_outside_prepare_df(
        res_df, 'Average score inside DTR/gap',
        'Average score outside DTR/gap', 'Average score')
    
    # Define y-range for both the violin plots showing average scores
    min_avg_genome = min(res_df['Average score'])
    max_avg_genome = max(res_df['Average score'])
    min_avg_in_out = min(in_VS_out_dtr_avg_scr_df['Average score'])
    max_avg_in_out = max(in_VS_out_dtr_avg_scr_df['Average score'])
    avg_scr_min = min(min_avg_genome, min_avg_in_out)
    avg_scr_max = max(max_avg_genome, max_avg_in_out)
    
    # (A)  AVERAGE SCORES: LPEG VS Others
    
    # Make and save violin plots
    avg_scores_violin = (
        ggplot() +
        geom_violin(res_df, aes(x='PE category', y='Average score')) +
        ylim(avg_scr_min, avg_scr_max)
    )
    filepath = os.path.join(ANALYSIS_NAME, 'plots', 'average_scores_LPEG_VS_others.png')
    avg_scores_violin.save(filepath, height=VIOL_HEIGHT, width=VIOL_WIDTH, dpi=VIOL_DPI, verbose = False)
    
    # Mann-Whitney U test
    lpeg_avg_scr_list = list(res_df.loc[res_df['PE category']=='LPEG'].loc[:,'Average score'])
    others_avg_scr_list = list(res_df.loc[res_df['PE category']=='Others'].loc[:,'Average score'])
    mwu = stats.mannwhitneyu(lpeg_avg_scr_list, others_avg_scr_list, alternative='greater')
    save_mwu('average_scores_LPEG_VS_others_MWU_test.txt',
             'Average scores distribution: LPEG VS Others', mwu)
    
    
    # (B)  LPEG AVG SCORES: inside VS outside DTR/gap
    
    # Make and save violin plots
    LPEG_avg_scores_inside_outside = (
        ggplot() +
        geom_violin(in_VS_out_dtr_avg_scr_df, aes(x='Class', y='Average score')) +
        ylim(avg_scr_min, avg_scr_max)
    )
    filepath = os.path.join(ANALYSIS_NAME, 'plots', 'LPEG_average_scores_inside_VS_outside.png')
    LPEG_avg_scores_inside_outside.save(filepath, height=VIOL_HEIGHT, width=VIOL_WIDTH, dpi=VIOL_DPI, verbose = False)
    
    # Mann-Whitney U test
    mwu = stats.mannwhitneyu(inside_list, outside_list, alternative='greater')
    save_mwu('LPEG_average_scores_inside_VS_outside_MWU_test.txt',
             'LPEG average scores distribution: Inside VS Outside DTR/GAP region', mwu)
    
    
    # (C)  INTERGENIC FREQUENCIES: LPEG VS Others
    
    # Slice dataframe keeping only genomes with annotated genes, so that the
    # intergenic frequences can all be turned from string into float
    only_with_genes_df = res_df.loc[res_df['Intergenic freq'] != 'no_gene_was_annotated'].copy()
    df_column_to_float(only_with_genes_df, 'Intergenic freq')
    
    # Make and save violin plots
    intergenic_freq_violin = (
        ggplot() + geom_violin(only_with_genes_df, aes(x='PE category', y='Intergenic freq'))
    )
    filepath = os.path.join(ANALYSIS_NAME, 'plots', 'intergenic_frequency_LPEG_VS_others.png')
    intergenic_freq_violin.save(filepath, height=VIOL_HEIGHT, width=VIOL_WIDTH, dpi=VIOL_DPI, verbose = False)
    
    # Mann-Whitney U test
    lpeg_interg_freq_list = list(only_with_genes_df.loc[only_with_genes_df['PE category']=='LPEG'].loc[:,'Intergenic freq'])
    others_interg_freq_list = list(only_with_genes_df.loc[only_with_genes_df['PE category']=='Others'].loc[:,'Intergenic freq'])
    mwu = stats.mannwhitneyu(lpeg_interg_freq_list, others_interg_freq_list, alternative='greater')
    save_mwu('intergenic_frequency_LPEG_VS_others_MWU_test.txt',
             'Intergenic frequency distribution: LPEG VS Others', mwu)
    
    
    # - - - Supplementary plots - - -
    
    # (A)  AVERAGE SCORES (LPEG VS Other Alpha-infecting phages)
    
    # Keep 'Alphaproteobacteria phages', 'Caulobacter phages' and 'vB_BsubS-Delta'
    # because they infect alphaproteobacteria (i.e. exclude 'Other phages' category)
    only_alphas = res_df[res_df['Category'] != 'Other phages'].copy()
    only_alphas['PE category'] = only_alphas['PE category'].replace(
        'Others', 'Other Alpha-infecting phages')
    
    avg_score_lpeg_vs_other_alpha = (
        ggplot() +
        geom_violin(only_alphas, aes(x='PE category', y='Average score')) +
        ylim(avg_scr_min, avg_scr_max)
    )
    filepath = os.path.join(ANALYSIS_NAME, 'plots', 'average_scores_LPEG_VS_alphas.png')
    avg_score_lpeg_vs_other_alpha.save(filepath, height=VIOL_HEIGHT, width=VIOL_WIDTH, dpi=VIOL_DPI, verbose = False)
    
    # Mann-Whitney U test
    lpeg_avg_scr_list = list(only_alphas.loc[only_alphas['PE category']=='LPEG'].loc[:,'Average score'])
    other_alpha_avg_scr_list = list(only_alphas.loc[only_alphas['PE category']!='LPEG'].loc[:,'Average score'])
    mwu = stats.mannwhitneyu(lpeg_avg_scr_list, other_alpha_avg_scr_list, alternative='greater')
    save_mwu('average_scores_LPEG_VS_alphas_MWU_test.txt',
              'Average scores distribution: LPEG VS other Alphaproteobacteria-infecting phages', mwu)
    
    
    # (B)  INTERGENIC FREQUENCIES (LPEG VS Other Alpha-infecting phages)
    
    # Keep 'Alphaproteobacteria phages', 'Caulobacter phages' and 'vB_BsubS-Delta'
    # because they infect alphaproteobacteria (i.e. exclude 'Other phages' category)
    only_alphas_only_genes = only_with_genes_df[only_with_genes_df['Category'] != 'Other phages'].copy()
    only_alphas_only_genes['PE category'] = only_alphas_only_genes['PE category'].replace(
        'Others', 'Other Alpha-infecting phages')
    
    interg_f_lpeg_vs_other_alpha = (
        ggplot() +
        geom_violin(only_alphas_only_genes, aes(x='PE category', y='Intergenic freq'))
    )
    filepath = os.path.join(ANALYSIS_NAME, 'plots', 'intergenic_frequency_LPEG_VS_alphas.png')
    interg_f_lpeg_vs_other_alpha.save(filepath, height=VIOL_HEIGHT, width=VIOL_WIDTH, dpi=VIOL_DPI, verbose = False)
    
    # Mann-Whitney U test
    lpeg_interg_freq_list = list(only_alphas_only_genes.loc[only_alphas_only_genes['PE category']=='LPEG'].loc[:,'Intergenic freq'])
    other_alpha_interg_freq_list = list(only_alphas_only_genes.loc[only_alphas_only_genes['PE category']!='LPEG'].loc[:,'Intergenic freq'])
    mwu = stats.mannwhitneyu(lpeg_interg_freq_list, other_alpha_interg_freq_list, alternative='greater')
    save_mwu('intergenic_frequency_LPEG_VS_alphas_MWU_test.txt',
              'Intergenic frequency distribution: LPEG VS other Alphaproteobacteria-infecting phages', mwu)
    
    
    # New supplementary violin
    dtr_freq_for_interg_df = get_dtr_freq_for_interg_sites_df(res_df)
    
    interg_inside_dtr_freq = (
        ggplot() +
        geom_violin(dtr_freq_for_interg_df, aes(x='PE category', y='DTR/gap frequency for intergenic sites'))
    )
    filepath = os.path.join(ANALYSIS_NAME, 'plots', 'inside_DTR_frequency_for_intergenic_sites.png')
    interg_inside_dtr_freq.save(filepath, height=VIOL_HEIGHT, width=VIOL_WIDTH, dpi=VIOL_DPI, verbose = False)
    
    # - - - - - - - - - - - - - - - 
    
    
    # ENTROPY AND GINI PLOTS
    # ----------------------
    
    print('Generating normalized Entropy and normalized Gini plots ...')
    
    # Generate a negative control set
    min_n_sites, max_n_sites = min(res_df['Number of sites']), max(res_df['Number of sites'])
    neg_control = generate_neg_control(min_n_sites, max_n_sites)
    
    # Apply jitter
    res_df_j = jitter(
        res_df,
        ['Number of sites', 'Normalized positional entropy', 'Normalized Gini coefficient'],
        50)
    
    # Control order of the legend labels
    legend_labels_order = [
        'Caulobacter phages',
        'Alphaproteobacteria phages',
        'Other phages',
        'vB_BsubS-Delta',
        'control']
    res_df_j = res_df_j.assign(
        Category = pd.Categorical(res_df_j['Category'], legend_labels_order, ordered=True)
    )
    
    # Entropy plot
    entr_plot = get_scatterplot(res_df_j, 'Normalized positional entropy',
                                legend_pos=(.75, .34), neg_control_df=neg_control)
    
    filepath = os.path.join(plots_dir, 'norm_pos_entropy_VS_number_of_sites.png')
    entr_plot.save(filepath, height=IMAGE_HEIGHT, width=IMAGE_WIDTH,
                   dpi=IMAGE_RESOLUTION, verbose = False)
    
    # Gini plot
    gini_plot = get_scatterplot(res_df_j, 'Normalized Gini coefficient',
                                legend_pos=(.75, .64), neg_control_df=neg_control)
    
    filepath = os.path.join(plots_dir, 'norm_gini_coeff_VS_number_of_sites.png')
    gini_plot.save(filepath, height=IMAGE_HEIGHT, width=IMAGE_WIDTH,
                   dpi=IMAGE_RESOLUTION, verbose = False)
    
    evenn_plot = get_scatterplot(res_df_j, 'Normalized Evenness',
                                legend_pos=(.75, .34), neg_control_df=neg_control)
    filepath = os.path.join(plots_dir, 'norm_evenness_VS_number_of_sites.png')
    evenn_plot.save(filepath, height=IMAGE_HEIGHT, width=IMAGE_WIDTH,
                   dpi=IMAGE_RESOLUTION, verbose = False)
    
    plt.scatter(res_df_j['Normalized Evenness'], res_df_j['Number of sites'])
    
    # ENTROPY VS INTERGENIC FREQUENCY
    # -------------------------------
    
    print('Generating normalized Entropy VS Intergenic frequency plot ...')
    
    # Prepare dataset
    res_df_E_VS_I = res_df[res_df['Intergenic freq'] != 'no_gene_was_annotated'].copy()
    res_df_E_VS_I['Intergenic freq'] = res_df_E_VS_I['Intergenic freq'].astype(float)
    res_df_E_VS_I = res_df_E_VS_I.sort_values(by=['Intergenic freq'], ascending=False)
    
    # Apply jitter
    res_df_E_VS_I_j = jitter(
        res_df_E_VS_I,
        ['Normalized positional entropy', 'Intergenic freq'],
        50)
    
    nnn_plot = get_2D_plot(res_df_E_VS_I_j, 'Intergenic freq', 'Normalized positional entropy')
    filepath = os.path.join(plots_dir, 'norm_pos_entropy_VS_intergenic_freq.png')
    nnn_plot.save(filepath, height=IMAGE_HEIGHT, width=IMAGE_WIDTH,
                   dpi=IMAGE_RESOLUTION, verbose = False)
    
    print('\nDone. Plots are saved in', plots_dir)



if __name__=="__main__":
    
    set_up()
    
    main()
    
    





