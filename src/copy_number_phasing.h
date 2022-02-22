#ifndef COPY_NUMBER_PHASING_H
#define COPY_NUMBER_PHASING_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
using namespace std;

//////////////// linker include ////////////////
#include "read_linker_output.h"
#include "variant_site.h"
#include "coord_dict.h"
#include "bin_reference.h"
#include "mbin.h"

///////////// cn_cutoffs //////
#define min_switch 10  //10
#define min_switch_two 6  //10
#define minimum_snp 8
//#define first_pass_switch 4 // this should be small as it is used to phase between blocks
//#define first_pass_switch 0.3 // this should be small as it is used to phase between blocks
#define first_pass_switch 0.1 // this should be small as it is used to phase between blocks

/////////////// functions /////////////////////
void initialize_copy_num_map( cn_map& chromosome_map, coord_dictionary& pdict, variant_graph& vgraph, int binsize, std::vector<int>& bin_array );
void cn_phase_loop( cn_map& chromosome_map, std::vector<int>& good_bins );
double calc_switch( cn_map& chromosome_map, int index1, int index2 );
void cn_phasing( cn_map& chromosome_map, std::vector<int>& bin_array, std::vector<int>& good_bins, std::vector<int>& recovered_bins, std::vector<int>& merged_bins );
void cn_phasing_greg( coord_dictionary& pdict, variant_graph& vgraph, double delta, int range );
void flip_to_end( int next , cn_map& chromosome_map, std::vector<int>& good_bins );
void merge_bins( cn_map& chromosome_map, std::vector<int>& good_bins, std::vector<int>& merged_bins );
void cn_phase_between_blocks( cn_map& chromosome_map, std::vector<int>& good_bins );
double calc_switch_merge( cn_map& chromosome_map, std::vector<int> merged_index1, std::vector<int> merged_index2 );
void rescue_bad_bins( cn_map& chromosome_map, std::vector<int>& bin_array, std::vector<int>& good_bins, std::vector<int>& recovered_bins );

void cn_block_flip( coord_dictionary& pdict, variant_graph& vgraph );
void block_flip_recursive_cn( coord_dictionary& pdict, int range, int& prior, double cutoff, variant_graph& vgraph, double delta, std::vector<int> gbin );
void sliding_block_flip( coord_dictionary& pdict, int range, int& prior, double cutoff, variant_graph& vgraph, double delta, std::vector<int> gbin, bool max );

void bin_flip_recursive_cn( coord_dictionary& pdict, int range, int& prior, double cutoff, variant_graph& vgraph, double delta, std::vector<int> gbin, bool max );
void bin_flip_all( coord_dictionary& pdict, int range, int& prior, double cutoff, variant_graph& vgraph, double delta, std::vector<int> gbin );
void get_depth( coord_dictionary& pdict, int range, int& prior, double cutoff, variant_graph& vgraph, double delta, std::vector<int>& depth );

#endif  // COPY_NUMBER_PHASING_H
