import itertools
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os
from os import path
import pandas as pd
import pickle as pkl
import seaborn as sns


plt.rcParams["font.family"] = "Arial"
plt.rcParams['figure.figsize'] = [5, 1.8]
chrom_list = ['chr'+str(i) for i in range(1,23)]
chrom_list.append('chrX')


chr_bands = pd.read_csv('./genome_info/chromosome_banding_hg38.txt', sep='\t')
band_colors = []
for i in chr_bands['gieStain']:
    if i == 'acen':
        band_colors.append('Red')
    elif i == 'gneg':
        band_colors.append('White')
    elif i == 'gpos100':
        band_colors.append('Black')
    elif i == 'gpos25':
        band_colors.append('Silver')
    elif i == 'gpos50':
        band_colors.append('Gray')
    elif i == 'gpos75':
        band_colors.append('Gray')
    elif i == 'stalk':
        band_colors.append('None')
    elif i == 'gvar':
            band_colors.append('Blue')
chr_bands['color'] = band_colors

rpe1_truth = pd.read_csv('./cell_line_data/ground_truth/RPE1_Haplotype.txt.gz', sep='\t')
rpe1_truth = rpe1_truth[((rpe1_truth['linkage_filter']==1) & (rpe1_truth["haplotype_bulkSeq"]!=0))].reset_index(drop=True)

def getSize(chrom):
    chrom_sizes = pd.read_csv('./genome_info/hg38_chrom_sizes.txt', header=None, sep='\t')
    
    return int(chrom_sizes[chrom_sizes[0]==chrom][1].values[0])

def correctGC(depth_array, hg38_bins, bins_10kb):
    GC_bins = pd.read_csv('./genome_info/GCBins_10kb.GC.csv', header=None) #Average GC content for 50kb bins
    bins_10kb['GC_avg'] = GC_bins
    bins_10kb = bins_10kb[bins_10kb['NonNBase']>0]

    depth_mask = pd.concat([bins_10kb['NonNBase']<9000 for i in range(depth_array.shape[1])], axis=1)
    depth_mask.index = depth_array.index
    depth_mask.columns = depth_array.columns

    depth_array = depth_array.mask(depth_mask)
    depth_log2 = np.log2(depth_array)
    depth_log2 = depth_log2 - np.nanmedian(depth_log2, axis=0) #Normalizing each library by its genome-wide median cov

    GC_strata = pd.read_csv('./genome_info/GC_strata.csv', header=None)
    GC_bin_id = pd.read_csv('./genome_info/GC_bin_id.csv', header=None)
    good_bins = hg38_bins[(hg38_bins['Chr']!='chrY') & (hg38_bins['NonNBase']>0)].index
    GC_bin_id = GC_bin_id.iloc[good_bins, :]

    depth_log2.index = good_bins

    median_gc_strength = np.empty([GC_strata.shape[1]-1, depth_log2.shape[1]])
    for GC_id in range(GC_strata.shape[1]-1):
        binID = np.where(GC_bin_id[0]==GC_id)
        counts = depth_log2.iloc[binID[0].tolist(), :]
        no_gc_bins = counts.shape[0]-counts.isna().sum()
        median_gc_strength[GC_id, :] = np.nanmedian(counts, axis=0)
        median_gc_strength[GC_id, no_gc_bins<100] = np.nan

    GCcorr_logCov = depth_log2
    A = GCcorr_logCov.iloc[np.where([GC_bin_id[0]>0][0]==True)]
    B = pd.DataFrame(np.vstack(median_gc_strength[GC_bin_id.iloc[np.where([GC_bin_id[0]>0][0]==True)]]))
    B.index = A.index
    B.columns = A.columns

    GCcorr_logCov.iloc[np.where([GC_bin_id[0]>0][0]==True), :] = A-B
    C = pd.DataFrame(np.zeros((GCcorr_logCov.shape[0], GCcorr_logCov.shape[1]))+2)
    C.index = GCcorr_logCov.index
    C.columns = GCcorr_logCov.columns
    GCcorr_Cov = C.pow(GCcorr_logCov)
    
    return GCcorr_Cov, bins_10kb #return the GC-corrected read counts and filtered 10kb bins

def reformat_integrated_phasing(hapcut_path):
    hap = pd.read_csv(hapcut_path, sep='\t', comment='#', header=None)
    gt_entry = pd.Series([i.split(':')[0] for i in hap[9]])
    filter_hets = lambda x: x[gt_entry.isin(pd.Series(['0|1', '1|0']))].reset_index(drop=True)
    hap = filter_hets(hap)
    gt_entry = filter_hets(gt_entry)
    gt_hap = [int(i[0])-int(i[2]) for i in gt_entry]
    hap[10] = gt_hap
    hap[4] = hap[10]
    hap = hap[hap[4]!=0].reset_index(drop=True)
    
    return hap

def calculate_accuracy_hcc1954(test_hap, chrom):
    test_hap = test_hap[test_hap[4]!=0]
    truth_path = "./cell_line_data/ground_truth/hap_full_scaffold_nov10_BL1954_2000_"
    truth = pd.read_csv(truth_path+chrom+'.dat', sep='\t', header=None)
    truth = truth[(truth[5]!=0)]
    truth[2]=truth[2]+1
    overlap = list(set(test_hap[1]) & set(truth[2]))
    test_hap = test_hap[test_hap[1].isin(overlap)].reset_index(drop=True)
    truth = truth[truth[2].isin(overlap)].reset_index(drop=True)
    test_hap[11] = truth[5]
    test_hap[12] = [i*j for i,j in zip(test_hap[4].values, test_hap[11].values)]
    accuracy = max(list(test_hap[12]).count(1), list(test_hap[12]).count(-1))/len(list(test_hap[12]))
    print(chrom, accuracy)
    
    return accuracy, test_hap

def calculate_accuracy_recovered_hcc1954(test_hap, chrom):
    test_hap["pos"] = test_hap["pos"]+1
    test_hap = test_hap[test_hap["hap"]!=0]
    
    truth_path = "./cell_line_data/ground_truth/hap_full_scaffold_nov10_BL1954_2000_"
    truth = pd.read_csv(truth_path+chrom+'.dat', sep='\t', header=None)
    truth = truth[(truth[5]!=0)]
    truth[2]=truth[2]+1
    
    overlap = list(set(test_hap["pos"]) & set(truth[2]))
    test_hap = test_hap[test_hap["pos"].isin(overlap)]
    truth = truth[truth[2].isin(overlap)]
    
    test_hap["mlinker_hap"] = truth[5].values
    test_hap["agreement"] = [i*j for i,j in zip(test_hap["hap"].values, test_hap["mlinker_hap"].values)]
    accuracy = max(list(test_hap["agreement"]).count(1),list(test_hap["agreement"]).count(-1))/len(list(test_hap["agreement"]))
    print("refLinker recovered accuracy", chrom, accuracy)
    
    return accuracy, test_hap


def calculate_accuracy_rpe1(test_hap, chrom):
    test_hap[1] = test_hap[1]+1
    test_hap = test_hap[test_hap[4]!=0]

    chrom_truth = rpe1_truth[rpe1_truth["chr"]==chrom].copy()
    overlap = list(set(test_hap[1]) & set(chrom_truth["pos"]))
    test_hap = test_hap[test_hap[1].isin(overlap)]
    chrom_truth = chrom_truth[chrom_truth["pos"].isin(overlap)]
    test_hap[11] = chrom_truth["haplotype_bulkSeq"].values
    test_hap[12] = [i*j for i,j in zip(test_hap[4].values, test_hap[11].values)]
    accuracy = max(list(test_hap[12]).count(1), list(test_hap[12]).count(-1))/len(list(test_hap[12]))
    print("refLinker accuracy", chrom, accuracy)
    
    return accuracy, test_hap

def calculate_accuracy_recovered_rpe1(test_hap, chrom):
    test_hap["pos"] = test_hap["pos"]+1
    test_hap = test_hap[test_hap["hap"]!=0]
    
    chrom_truth = rpe1_truth[rpe1_truth["chr"]==chrom].copy()
    overlap = list(set(test_hap["pos"]) & set(chrom_truth["pos"]))
    test_hap = test_hap[test_hap["pos"].isin(overlap)]
    chrom_truth = chrom_truth[chrom_truth["pos"].isin(overlap)]
    test_hap["mlinker_hap"] = chrom_truth["haplotype_bulkSeq"].values
    test_hap["agreement"] = [i*j for i,j in zip(test_hap["hap"].values, test_hap["mlinker_hap"].values)]
    accuracy = max(list(test_hap["agreement"]).count(1), list(test_hap["agreement"]).count(-1))/len(list(test_hap["agreement"]))
    print("refLinker recovered accuracy", chrom, accuracy)
    
    return accuracy, test_hap
    
    

def calculate_eagle2_accuracy_rpe1(test_hap, chrom):
    test_hap[1] = test_hap[1]+1
    test_hap = test_hap[test_hap[5]!=0]

    chrom_truth = rpe1_truth[rpe1_truth["chr"]==chrom].copy()
    overlap = list(set(test_hap[1]) & set(chrom_truth["pos"]))
    test_hap = test_hap[test_hap[1].isin(overlap)]
    chrom_truth = chrom_truth[chrom_truth["pos"].isin(overlap)]
    test_hap[11] = chrom_truth["haplotype_bulkSeq"].values
    test_hap[12] = [i*j for i,j in zip(test_hap[5].values, test_hap[11].values)]
    accuracy = max(list(test_hap[12]).count(1), list(test_hap[12]).count(-1))/len(list(test_hap[12]))
    print("EAGLE2 accuracy", chrom, accuracy)
    
    return accuracy, test_hap


def calculate_acn(counts, hap, allelic_depths, chrom):
    
    counts = counts[counts["CONTIG"]==chrom].reset_index(drop=True)
    allelic_depths = allelic_depths[allelic_depths["contig"]==chrom].reset_index(drop=True)

    phased_overlap = list(set(hap[1]).intersection(set(allelic_depths["position"])))

    hap = hap[hap[1].isin(phased_overlap)]
    allelic_depths = allelic_depths[allelic_depths["position"].isin(phased_overlap)]

    allelic_depths["haplotype"] = hap[4].values

    allelic_depths["AR"] = ((1+allelic_depths.haplotype)*allelic_depths.refCount 
    + (1-allelic_depths.haplotype)*allelic_depths.altCount) / (2*(allelic_depths.refCount + allelic_depths.altCount))

    chrom_bins_10kb = range(0, getSize(chrom), int(1e4))
    count_ids = np.digitize(counts.START, chrom_bins_10kb)
    snp_ids = np.digitize(allelic_depths.position, chrom_bins_10kb)
    counts.index = count_ids
    allelic_depths['normCov'] = counts.loc[snp_ids].COUNT.values

    allelic_depths['CN_A'] = 4 * allelic_depths.AR.values * allelic_depths.normCov.values
    allelic_depths['CN_B'] = 4 * (1 - allelic_depths.AR.values) * allelic_depths.normCov.values
    
    return allelic_depths
    