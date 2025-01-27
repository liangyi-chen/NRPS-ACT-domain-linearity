# NRPS-ACT-domain-linearity
This script takes the input files downloaded from MIBiG database, which can be found in the subdirectory /MIBiG_data/, and output a list of the following (in .csv format):

1. MIBiG number;
2. Type, if the entry is a NRPS, or PKS, or both;
3. Compound name;
4. Molecular weight;
5. Numbers of A domain, C domain and T domain within.
   
The domain counts and the molecular weights will also be plotted in interactive html graphs, where the user can input their own compound to check if it might be a linear or nonlinear NRPS.

The output files are generated in the subdirectory /Output/. 

# Explanation of abbreviations

- **NRPS**: Nonribosomal peptide synthetase
- **A domain**: Adenylation domain
- **C domain**: Condensation domain
- **T domain**: Termination domain

# Publication

This data analysis method was applied in [Yang, J., Balutowski, A., Trivedi, M., & Wencewicz, T. A. (2025). Chemical Logic of Peptide Branching by Iterative Nonlinear Nonribosomal Peptide Synthetases. *Biochemistry*](https://pubs.acs.org/doi/full/10.1021/acs.biochem.4c00749).
