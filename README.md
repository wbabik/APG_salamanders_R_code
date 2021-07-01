# APG_salamanders_R_code
R project with scripts performing most of the analyses reported in Palomar et al. *Coevolution between MHC class I and Antigen Processing Genes in salamanders*

The remaining, simple scripts plus all the data are on figshare [here](https://doi.org/10.6084/m9.figshare.c.5481927). Note that some files deposited on figshare are also in this R project for convenience.

## Usage
As some of the scripts use output from other scripts, they should be run in sequence:

1. `R_sequence_logos_classical_nonclassical.R`  - generates Fig. S1 (sequence logo) and assigns status to potentially functional MHC alleles, based on the number of conserved amino acids at key anchor residues; alleles with status **"class"** form the **"conserved anchor"** dataset.
2. `R_MHC_main.R` - performs most analyses on MHC (and *BRD2*), which was also studied via amplicon sequencing).
3. `R_microhaplotypes.R` - generates sequences of segment microhaplotypes for all species and individuals. Its inputs are included, they come from bash and R scripts that can be found [here](https://doi.org/10.6084/m9.figshare.c.5481927).
4. `R_microhaplotypes_diversity.R` takes output from `R_microhaplotypes.R` and calculates diversity for each segment.
5. `R_phylo_corr.R` takes output from `R_microhaplotypes_diversity.R` and `R_MHC_main.R`, calculates per gene diversities and does PGLS modeling.
6. `R_phylo_corr_summaries_diagnostics.R` summarises results of PGLS modeling and does some extra diagnostics.
7. `R_Tables_and_Figures.R` generates several Figures and Tables taking outputs from other scripts.
8. `R_intraspecific_corr.R` calculates intraspecific (individual based) correlations between MHC and APG diversity. These results were not reported in the paper.
9. `R_Fig_1.R` generates Figure 1.
10. `R_TAP_protein_polymorphism.R` summarises polymorphism at potentially functionally relevant residues of TAP proteins, producing part of Table S7. 

If you experience some unexpected errors when running these scripts, a possible reason are name conflicts between functions from different packages. In such a case try to run the scripts in a clean environment.

