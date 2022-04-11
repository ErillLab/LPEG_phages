# LPEG_phages
Code and data for analysis of TF-binding sites in bacteriophage genomes.  
LPEG (Low Positional Entropy Genome) phages show a skewed distribution of predicted TF-binding sites.

## Pipeline used for our work on Alphaproteobacteria-infecting phages
### 1. Retrieve genomes
``retrieve_relevant_genomes.py``  

``retrieve_relevant_genomes_config.json``

### 2. PSSM-scan
``motif_matching_scoring.py``  

``motif_matching_scoring_config.json``

### 3. Statistical analysis
``analyze_sites_distribution.py``  

``analyze_sites_distribution_config.json``

### 4. Generate plots
``generate_plots.py``  

``generate_plots_config.json``


