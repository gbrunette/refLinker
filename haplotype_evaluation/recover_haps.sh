#!/bin/bash

for chrom in {1..22} X;
do
python eagle2_recover.py ./cell_line_data/reflinker_haplotypes/pop_hap_solution_20231224_HCC1954_cancer_chr${chrom}.dat \
../haplotype_data/HCC1954BL_chr${chrom}_22-8-8.vcf.gz \
50000 \
chr${chrom} \
HCC1954_recovered_haplotype_chr${chrom}.dat

#python eagle2_recover.py ./cell_line_data/reflinker_haplotypes/pop_hap_solution_20240108_RPE1_SRS1045724_5_chr${chrom}.dat \
#../haplotype_data/RPE1_chr${chrom}_21-6-4.vcf.gz \
#50000 \
#chr${chrom} \
#RPE1_recovered_haplotype_chr${chrom}.dat
done
