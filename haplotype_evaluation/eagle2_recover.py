import numpy as np
import pandas as pd
import sys

hap_path = sys.argv[1]
vcf_path = sys.argv[2]
dist_cutoff = int(sys.argv[3])
chrom = sys.argv[4]
outfile = sys.argv[5]

def match_positions(unphased_positions, phased_positions):
    
    matched_positions = []
    unphased_counter = 0
    phased_counter = 0
    
    for i in range(len(unphased_positions) + len(phased_positions)):
        if ((unphased_counter < len(unphased_positions)) & (phased_counter < len(phased_positions))):
            if unphased_positions[unphased_counter] < phased_positions[phased_counter]:
                if phased_counter == 0:
                    matched_positions.append(phased_counter)
                elif abs(unphased_positions[unphased_counter] - phased_positions[phased_counter]) < abs(unphased_positions[unphased_counter] - phased_positions[phased_counter-1]):
                    matched_positions.append(phased_counter)
                else:
                    matched_positions.append(phased_counter-1)
                unphased_counter += 1
            else:
                phased_counter += 1
        elif phased_counter < len(phased_positions):
            return matched_positions
        else:
            matched_positions.append(len(phased_positions)-1)
            
    return matched_positions

def eagleRecover(hap_path, vcf_path, dist_cutoff, chrom):

    vcf = pd.read_csv(vcf_path,
                     comment='#', sep='\t', header=None)
    vcf = vcf[vcf[0]==chrom].reset_index(drop=True)

    parse = lambda x: x.split(':')[0].split('|')
    hap = pd.DataFrame((vcf[9].apply(parse)).to_list(), columns=[12,13]).astype('int')
    vcf['eagleHap'] = hap[12]-hap[13]
    vcf = vcf[vcf['eagleHap']!=0].reset_index(drop=True)
    vcf[1] = vcf[1]-1

    #Loading the refLinker haplotype

    linker = pd.read_csv(hap_path, header=None, sep='\t', index_col=0)
    phased = linker.copy()[linker[4] !=0].reset_index(drop=True)
    phased = phased.copy()[[1,4]]

    overlap = set(vcf[1]) & set(phased[1])
    phased = phased[phased[1].isin(overlap)].reset_index(drop=True)
    unphased = vcf[np.invert(vcf[1].isin(overlap))].reset_index(drop=True)

    unphased_positions = np.array(unphased[1])
    phased_positions = np.array(phased[1])

    nearest_ind = match_positions(unphased_positions, phased_positions)
    rescue_df = phased.copy().iloc[nearest_ind,:].reset_index(drop=True)
    distances = pd.Series(np.abs(unphased.copy()[1]-rescue_df.copy()[1]))
    a,b,c = np.intersect1d(vcf[1], rescue_df[1], return_indices=True)
    d, e = np.unique(rescue_df[1], return_counts=True)
    vcf_indices = np.repeat(b, e)
    recovered_linkage = vcf.copy().iloc[vcf_indices,list(vcf.columns).index('eagleHap')].reset_index(drop=True)*unphased.copy()['eagleHap']
    new_hap = rescue_df.copy()[4]*recovered_linkage
    unphased[4]  = new_hap.values
    unphased = unphased.copy()[[1,4]]
    unphased = unphased[distances <= dist_cutoff]
    phased = pd.concat([phased,unphased]).reset_index(drop=True).sort_values(by=1).reset_index(drop=True)
    phased = phased[phased[4]!=0]
    phased.columns = ["pos", "hap"]
    phased["chrom"] = chrom

    return phased[["chrom", "pos", "hap"]]


recovered_hap = eagleRecover(hap_path, vcf_path, dist_cutoff, chrom)

recovered_hap.to_csv(outfile, index=False, sep="\t")