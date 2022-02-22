#ifndef RUN_POP_HIC_PHASE_H
#define RUN_POP_HIC_PHASE_H

//////////////// c++ include //////////////////
#include <getopt.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;
//using namespace opt;

//////////////// linker include ////////////////
#include "read_bam.h"
#include "read_vcf.h"
#include "read_linker_output.h"
#include "map_matrix.h"
#include "sub_matrix.h"
#include "mc_solver.h"
#include "coord_dict.h"
#include "variant_site.h"
#include "build_hap_matrix.h"
#include "write_linker_output.h"
#include "hic_links.h"

//////////////// definitions //////////////////

/////////////// functions /////////////////////
void run_pop_phaser( int argc, char** argv);
static void parse_pop_phaser_options( int argc, char** argv );
//void create_hic_links_pop( variant_graph& hic_vgraph, hic_link& hlink );
void init_hic_pop_matrix( map_matrix<double>& diff_matrix, map_matrix<int>& num_matrix_second, coord_dictionary& pdict, variant_graph& hic_vgraph, int centromere_pos );
void match_pdict( vcf_vector& vvec, coord_dictionary& pdict );

//static void parse_region_string( std::string samtools_region );

#endif  // RUN_POP_HIC_PHASE_H
