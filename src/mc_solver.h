#ifndef MC_SOLVER_H
#define MC_SOLVER_H

//////////////// c++ include //////////////////
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <list> 
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <math.h>
#include <omp.h>
using namespace std;

//////////////// linker include ////////////////
#include "coord_dict.h"
#include "block_dict.h"
#include "map_matrix.h"
#include "sub_matrix.h"

///////////// solver_cutoffs //////
#define solver_loops 10                 // 10  // {val} the number of spin flip block flip loops
#define solver_loops_hic 10                 // 10  // {val} the number of spin flip block flip loops
#define pos_diff_cutoff 1000000          // {val} if the maximum delta genome distance - band width

/////////////// functions /////////////////////
void solver( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict,map_matrix<double> diff_matrix,map_matrix<int> num_matrix_second ); ///// fix this
static void set_switch_energy( coord_dictionary& pdict, vector<bool>& flip_maxima );
static void get_flip_positions( coord_dictionary& pdict, vector<bool>& flip_maxima );
static void energy_sum_min_max( coord_dictionary& pdict );
static void single_spin_flip_map( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix );
static void block_flip_map( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix, vector<bool>& flip_maxima );
static double energy_function_opt( submatrix_opt tempsub, vector<int> haplotype );
static double dot_product( vector<int> a, vector<double> b );
static vector<double> dot_product_matrix( vector<int> a, std::vector< std::vector<double> > b );
void call_blocks(coord_dictionary& pdict, int switch_cutoff);

/////////////////
void solver_recursive_pop( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<int> num_matrix, map_matrix<double> diff_matrix, int window_size, double cutoff, double prune );
void no_switch( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<int> num_matrix, map_matrix<double> diff_matrix, int window_size, double cutoff );

void greedy_cut( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<int> num_matrix, map_matrix<double> diff_matrix, int window_size, double cutoff, double prune );
double calculateOldCut( coord_dictionary pdict, set<int> S1, set<int> S2, vector<int> loop_hap, map_matrix<int> nmatrix, map_matrix<double> diff_matrix );
double calculateCut( coord_dictionary& pdict, set<int> vtcs, set<int> S1, set<int> S2, vector<int> loop_hap, map_matrix<int> nmatrix, map_matrix<double> diff_matrix );

void solver_recursive( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<double> diff_matrix, map_matrix<int> num_matrix_second );
void solver_recursive_hic( block_dictionary& bdict, std::vector<int> hic_limit_loop, map_matrix<int> block_matrix );


static void calculate_block_flip( coord_dictionary& pdict, int range, double& switch_min, bool do_flip );
static void block_flip_recursive_var_range( coord_dictionary& pdict, map_matrix<int>& nmatrix, int range, int& prior, int pop_weight, map_matrix<double> diff_matrix, double cutoff, int& bkp, bool switch_hap, double& switch_min );
static void block_flip_recursive_pop( coord_dictionary& pdict, map_matrix<int>& nmatrix, int dist, int pop_weight, int& prior, map_matrix<double> diff_matrix );

static void block_flip_scaffold( coord_dictionary& pdict, map_matrix<int>& nmatrix, int dist, int& prior, map_matrix<double> diff_matrix );

static void single_spin_flip_recursive_pop( coord_dictionary& pdict, map_matrix<int>& nmatrix );
static void flipE ( coord_dictionary& pdict, map_matrix<int>& nmatrix, int pop_weight, map_matrix<double> diff_matrix );
static void switchE( coord_dictionary& pdict, map_matrix<int>& nmatrix );

static void switchE_recursive( coord_dictionary& pdict, map_matrix<int>& nmatrix );
static void switchE_recursive_window( coord_dictionary& pdict, map_matrix<int>& nmatrix );
static void calculate_switchE( coord_dictionary& pdict, map_matrix<int>& nmatrix );
static void calculate_switchE_pop( coord_dictionary& pdict, map_matrix<int>& nmatrix, int pop_weight, map_matrix<double> diff_matrix );
static void switchE_recursive_pop( coord_dictionary& pdict, map_matrix<int>& nmatrix, int pop_weight, map_matrix<double> diff_matrix, double cutoff, int& bkp, bool flip, double& loop_min );
double calculate_flipE( coord_dictionary& pdict, map_matrix<int> nmatrix, int range, int pop_weight, map_matrix<double> diff_matrix, vector<int>& flip_range );

static void block_flip_fast( coord_dictionary& pdict, map_matrix<int>& nmatrix, int dist, int pop_weight, int& prior, map_matrix<double> diff_matrix );

static void single_spin_flip_recursive( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix );
static void block_flip_recursive( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix );
static void length_cutoff_nmatrix( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix, map_matrix<int>& nmatrix2 );
static void length_cutoff_nmatrix2( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix, map_matrix<int>& nmatrix2, double prune );

static void single_spin_flip_recursive_hic( block_dictionary& bdict );
static void block_flip_recursive_hic( block_dictionary& bdict );
static void energy_sum_min_max_hic( block_dictionary& bdict );
static void block_flip_brute_force_hic( block_dictionary& bdict );
//static void block_flip_brute_force2_hic( block_dictionary& bdict );

static void block_phasing_quick_flip( block_dictionary& bdict, vector<int>& graph_flip_list );
static void traverse_graph( std::unordered_map<int,std::unordered_map<int,int>> high_correlation_graph, vector<int>& traverse_list, int start_node, int subset_length );

#endif  // MC_SOLVER_H
