#include "run_cn_phaser.h"

///////////////////////////////////////////////
namespace opt {
        static std::string input_hap_file = "./output/hap_solution_default_tenx_chr20.dat";
        static std::string input_cov_file = "./output/het_coverage_default_tenx_chr20.dat";
        static std::string id_string = "default";
        static std::string input_vcf_file = "./sample_data/BL1954_PCRFree.hets.recalibrated.vcf";
        static std::string chr_choice = "chr18";
        static double delta = 4.0;
        static int range = 50;
	static int binsize = 10000;  // 10kb bins have 8-10 het sites on average
};

static const char* shortopts = "ho:l:m:n:v:b:c:d:r:";
static const struct option longopts[] = {
        { "hap-file",    no_argument, NULL, 'l' },
        { "cov-file",    no_argument, NULL, 'm' },
        { "id_string",   no_argument, NULL, 'n' },
        { "input_vcf_file", no_argument, NULL, 'v' },
        { "binsize",     no_argument, NULL, 'b' },
        { "chr_choice",     no_argument, NULL, 'c'},
        {"delta",   no_argument, NULL, 'd'}, 
        {"range", no_argument, NULL, 'r'}
};

///////////////////////////////////////////////
static const char *CN_PHASE_USAGE_MESSAGE =
"Usage: linker cn_phase [OPTION] -l /path/to/hap_solution.dat -m /path/to/het_coverage.dat \n\n"
"\n"
"  Options\n"
"  -l,      input haplotype solution path \n"
"  -m,      input het coverage path  \n"
"  -n,      id string for output files \n"
"  -b,      binsize - raw base number (default 10000 - 10kb) \n"
"\n";

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_cn_phaser_options( int argc, char** argv ) {
        bool die = false;
        if (string(argv[1]) == "help" || string(argv[1]) == "--help") {
                std::cerr << "\n" << CN_PHASE_USAGE_MESSAGE;
                exit(1);
        }
        for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
                std::istringstream arg(optarg != NULL ? optarg : "");
                switch (c) {
                case 'h': die = true; break;
                case 'l': arg >> opt::input_hap_file; break;
                case 'm': arg >> opt::input_cov_file; break;
                case 'n': arg >> opt::id_string; break;
                case 'v': arg >> opt::input_vcf_file; break;
                case 'b': arg >> opt::binsize; break;
                case 'c': arg >> opt::chr_choice; break;
                case 'd': arg >> opt::delta; break;
                case 'r': arg >> opt::range; break; 
                }
        }
        if (die) {
          std::cerr << "\n" << CN_PHASE_USAGE_MESSAGE;
          exit(1);
        }
        cout << endl;
        cout << "############### running linker copy number phase ### " << endl;
        cout << "== input hap  === " << opt::input_hap_file << endl;
        cout << "== input cov  === " << opt::input_cov_file << endl;
        cout << "== id string  === " << opt::id_string << endl;
        cout << "== binsize    === " << opt::binsize << endl;
        cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
        cout << endl;
};

////////////////////////////////////////////////////
// void run_cn_phaser(int argc, char** argv) {
//         parse_cn_phaser_options(argc, argv);
//         std::string outputFile = "./output/cn_phased_" + opt::id_string + ".dat";
//         variant_graph vgraph;
//         coord_dictionary pdict;
// 	//read_het_coverage( opt::input_cov_file, vgraph );
// 	//bool paired = true;
// 	//initialize_pdict( vgraph, pdict, paired );
// 	/////////////////////////////////////////
// 	//pdict.hap_zero_initialization();
// 	//read_hap_solution( opt::input_hap_file, pdict );
// 	/////////////////////////////////////////
// 	cn_map chromosome_map;
// 	std::vector<int> bin_array,good_bins,merged_bins,rec_bins;
// 	read_bin_haplotype_cn_data( opt::input_cov_file, chromosome_map, opt::binsize, bin_array );
// 	//initialize_copy_num_map( chromosome_map, pdict, vgraph, opt::binsize, bin_array );
// 	cout << " cn phasing section " << endl;
// 	cn_phasing( chromosome_map, bin_array, good_bins, rec_bins, merged_bins );
// 	write_cn_phased_bins( chromosome_map, outputFile, good_bins, merged_bins );
// 	//write_cn_phased( chromosome_map, outputFile, rec_bins, merged_bins );
//         //cout << " cn phasing coming soon " << endl;
//         return;
// };


void run_cn_phaser(int argc, char** argv) {
        parse_cn_phaser_options(argc, argv);
        std::string outputFile = "./output/cn_phased_" + opt::id_string + ".dat";
        variant_graph vgraph;
        coord_dictionary pdict;
        std::vector<int> pop_positions;
    vcf_vector vvec;
    vvec = load_vcf_file_pop_greg( opt::input_vcf_file, opt::chr_choice, pop_positions );

    read_het_coverage_greg( opt::input_cov_file, vgraph, vvec, pop_positions );
    
    bool paired = true;
    initialize_pdict( vgraph, pdict, paired );

    
    match_pdict_cn( vvec, pdict, pop_positions ); 
    pdict.hap_pop_initialization();
    
    /////////////////////////////////////////
    //pdict.hap_zero_initialization();
    //read_hap_solution( opt::input_hap_file, pdict );
    /////////////////////////////////////////
    cn_map chromosome_map;
    std::vector<int> bin_array,good_bins,merged_bins,rec_bins;

    //read_bin_haplotype_cn_data( opt::input_cov_file, chromosome_map, opt::binsize, bin_array );

    cn_phasing_greg( pdict, vgraph, opt::delta, opt::range);

    //initialize_copy_num_map( chromosome_map, pdict, vgraph, opt::binsize, bin_array );
    
    //cn_phasing( chromosome_map, bin_array, good_bins, rec_bins, merged_bins );
    write_cn_solution_greg( vgraph, outputFile, pdict, opt::chr_choice );
    //write_cn_phased( chromosome_map, outputFile, rec_bins, merged_bins );
        //cout << " cn phasing coming soon " << endl;
        return;
};

void match_pdict_cn( vcf_vector& vvec, coord_dictionary& pdict, std::vector<int> pop_positions ) {
        int count = 0;
        std::vector<int>::iterator pop_it;
        for (int i = 0; i < pdict.num_paired; i++) { pdict.ref_handle.push_back("ref"); pdict.alt_handle.push_back("alt"); }
        for (int i = 0; i < vvec.size(); i++) {
        pop_it = find(pop_positions.begin(),pop_positions.end(), vvec[i].pos);
        if ( pdict.pos_index.find(vvec[i].pos) != pdict.pos_index.end() && pop_it != pop_positions.end() ) {
            int j = pdict.pos_index[vvec[i].pos];
                        pdict.pop_hap[j] = vvec[i].pop_hap;
                        std::string ref_hash = std::to_string(vvec[i].pos) + "_" + vvec[i].ref_base + "_" + vvec[i].ref_base;
                        std::string alt_hash = std::to_string(vvec[i].pos) + "_" + vvec[i].ref_base + "_" + vvec[i].var_base;
                        pdict.ref_handle[j] = ref_hash;
                        pdict.alt_handle[j] = alt_hash;
                        count += 1;
                }
        }
        //#############################################################
        cout << "pdict number: " << pdict.num_paired << " matched count: " << count << endl;
};
















        //if ( opt::input_cov_file.substr(opt::input_cov_file.length() - 4) == ".vcf" ) {
        //   vcf_vector vvec;
        //   vvec = load_vcf_file_coverage(opt::input_cov_file,chromosome);
        //}


        //int chromosome = 0;  //2
        //if ( opt::input_cov_file.substr(opt::input_cov_file.length() - 4) == ".vcf" ) {
        //   vcf_vector vvec;
        //   vvec = load_vcf_file_coverage(opt::input_cov_file,chromosome);
        //}






