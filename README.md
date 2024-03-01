# refLinker
Chromosome-length haplotype determination from Hi-C and external reference panels

refLinker resolves complete chromosomal haplotypes using Hi-C data and statistical phasing methods. For detailed installation instructions, see https://github.com/rwtourdot/mlinker.

refLinker performs haplotype refinement on an initial haplotype "guess," here, the statistical phasing haplotype generated by EAGLE2 (https://alkesgroup.broadinstitute.org/Eagle/). 

## Extracting Hi-C links

Starting from a statistical phasing haplotype, `{statistical haplotype}.vcf.gz` and aligned Hi-C reads, `{Hi-C reads}.bam`, extracts Hi-C links between variant genotypes for a single chromosome with the following command:

`reflinker extract -v {statistical haplotype}.vcf.gz -i {Hi-C reads}.bam -n {sample name} -e hic -c {chrom}`

This generates map from reads to variants, written to the file `graph_variant_{sample name}_{chrom}.dat`

## Haplotype refinement

Chromosome-length haplotype inference is then performed using the reads-to-variants map and statistical phasing haplotype with the following command:

`reflinker pop -v {statistical haplotype} -g graph_variant_{sample name}_{chrom}.dat -c {chrom} -e -10.0 -p 0.999`

Where the arguments `-e` and `-p` specify the block-flipping penalty and linkage pruning cut-offs respectively. `reflinker pop` writes whole-chromosome haplotype to the file `pop_hap_solution_{sample name}_{chrom}.dat`, where column 1 corresponds to variant position, column 4 corresponds to refLinker haplotype, and column 5 corresponds to the initial (e.g. EAGLE2) haplotype assignment.  

## Haplotype recovery

The python script `eagle2_recover.py` calculates the haplotype linkage between common, unlinked variants omitted from the refLinker Hi-C scaffold, based on their EAGLE2 linkage to the closest linked variant within 5 kb.
