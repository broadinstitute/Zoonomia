# 200Mammals
Code used in Zoonomia Project analyses

# Diversity analysis (testing for correlations)

R script and input files for the diversity analysis in the Zoonomia paper

R script: diversity_analysis.correlation.R

Input files:
(1) diversity_analysis.input.txt
(2) diversity_analysis.pantheria.phenotypes.txt
(3) diversity_analysis.pantheria.phenotype_types.txt

Testing by regression with ordinal predictor for correlation of IUCN conservation status and SoH, heterozyogsity: IUCN.diversity.correlation.R

# Pipeline for calculating segments of heterozygosity from a mapped .bam file
# Scripts adapted to run on an SGE/bash environment
# paths are hardcoded in the scripts, matching directories need to be created 

1. Create a table of coverages per scaffold and discard anoumalous scaffolds:

qsub cov_and_disc.txt

2. Calculate average heterozygosity and windowed heterozygosity per sample

qsub hetcalc.txt

3. Calculate SoH per sample

for i in `ls -1 ./50kb_scaffolds_het/`; do 
  cat ${i}|python3 mammals.py > ${i}".soh"; 
done

