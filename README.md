# LPEG phages
Code and data used to generate the results in Mascolo *et al*. (2022).

**Reference**

E. Mascolo *et al.*, “The transcriptional regulator CtrA controls gene expression in Alphaproteobacteria phages: Evidence for a lytic deferment pathway,” *Front Microbiol*, vol. 0, p. 2999, 2022, doi: 10.3389/FMICB.2022.918015.

[Link to paper](http://journal.frontiersin.org/article/10.3389/fmicb.2022.918015/full?&utm_source=Email_to_authors_&utm_medium=Email&utm_content=T1_11.5e1_author&utm_campaign=Email_publication&field=&journalName=Frontiers_in_Microbiology&id=918015)

## Short description
Code used for the analysis of TF-binding sites in bacteriophage genomes.  
Phages that show a skewed positional distribution of predicted TF-binding sites over their genome sequence are named LPEG (Low Positional Entropy Genome) phages.

## Pipeline

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


