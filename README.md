# LPEG phages
Code and data for the analysis of TF-binding sites in bacteriophage genomes.  
Phages that show a skewed positional distribution of predicted TF-binding sites over their genome sequence are named LPEG (Low Positional Entropy Genome) phages.

## Pipeline used to identify LPEGs in our work on Alphaproteobacteria-infecting phages

### 1. Retrieve genomes
Script used to retrieve the accessions of the genomes used in the analysis.  
``retrieve_relevant_genomes.py``  

Settings file in JSON format.  
``retrieve_relevant_genomes_config.json``  

Remember to provide your email (and your API key if you have one) in the "entrez_parameters". Also check that the "sleep_time" parameter (wait time in seconds between queries to the database) conforms with the NCBI's Entrez query policies, to avoid being 'blacklisted' by NCBI.

### 2. PSSM-scan
Script used to look for putative TF-binding sites in the genome sequences using a PSSM.  
``motif_matching_scoring.py``  

Settings file in JSON format.  
``motif_matching_scoring_config.json``  

Remember to provide your email (and your API key if you have one) in the "entrez_parameters". Also check that the "sleep_time" parameter (wait time in seconds between queries to the database) conforms with the NCBI's Entrez query policies, to avoid being 'blacklisted' by NCBI.

### 3. Statistical analysis
Script used to perform the statistical analysis on the results from the PSSM-scan.  
``analyze_sites_distribution.py``  
It was used to study the positional distribution over the genome sequences of the predicted binding sites, as well as the frequency of intergenic predicted binding sites and their average score.  

Settings file in JSON format.  
``analyze_sites_distribution_config.json``

### 4. Generate plots
Script used to generate the Figures for our work on Alphaproteobacteria-infecting phages.  
``generate_plots.py``  

Settings file in JSON format.  
``generate_plots_config.json``


