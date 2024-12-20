#ifndef COORD_DICT_H
#define COORD_DICT_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
using namespace std;

//////////////// linker include ////////////////
#include "map_matrix.h"
#include "variant_site.h"

///////////////////// dictionary with general info
class coord_dictionary {
public:
    vector<double> deltaE,switchE;
    vector<double> deltaE_total,switchE_total;
    vector<int> block,block_total;
    vector<int> span_bound,span_bound_total;
    vector<int> haplotype;
    vector<int> pop_hap;
    vector<std::string> ref_handle,alt_handle;
    vector<bool> reload_bool;
    vector<bool> within_filter,within_filter_total;
    vector<int> sorted_paired_positions,double_positions,all_positions,sorted_all_positions;
    vector<int> up_bound_submatrix,low_bound_submatrix;
    vector<int> flip_up_bound,flip_low_bound;
    std::unordered_map<int,std::vector<std::string> > paired_dict;
    std::unordered_map<int,int> ref_index;
    std::unordered_map<int,int> pos_index;
    std::unordered_map<int,std::vector<int> > base_number;
    int num_paired,num_total;
    void initialize( std::unordered_map<std::string,variant_node>& , bool paired );
    void get_submatrix_bounds( map_matrix<int>& );
    void hap_random_initialization();
    void hap_pop_initialization();
    void hap_zero_initialization();
    void pop_hap_switch(int window, double prob );
};

//////////////// definitions //////////////////
//std::unordered_map<std::string,int> contig_dict;

#endif  // COORD_DICT_H
