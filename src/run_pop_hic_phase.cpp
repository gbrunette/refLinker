#include "run_pop_hic_phase.h"

/////////////// namespaces ///////////////////
namespace opt {
        static std::string chr_choice = "chr20";
	static std::string input_graph_file = "./output/graph_variant_jun10_BL1954_tenx_chr20.dat";
	static std::string input_vcf_file = "./sample_data/BL1954_PCRFree.hets.recalibrated.vcf";
        static std::string id_string = "default";
        static double window_size = 0.5;
	static int start_bound = 0;
	static int end_bound = 300000000;
        static double cutoff = -10.0;
};

/////////////// structures ///////////////////
static const char* shortopts = "ho:g:c:v:n:b:w:e:";
static const struct option longopts[] = {
	{ "help",        no_argument, NULL, 'h' },
        { "graph-file",  no_argument, NULL, 'g' },
        { "vcf-file",    no_argument, NULL, 'v' },
        { "chr-choice",  no_argument, NULL, 'c' },
        { "id_string",   no_argument, NULL, 'n' },
        { "window_size", no_argument, NULL, 'w' },
        { "cutoff", no_argument, NULL, 'e'},
        
};

///////////////////////////////////////////////////////////////////////////////////////////////////////
static const char *POP_PHASER_USAGE_MESSAGE =
"Usage: linker pop [OPTION] -v ./vcf_data/sample.vcf -g ./output/graph_variant_hic_chr20.dat -c chr20 -n trial \n\n"
"\n"
"  Options\n"
"  -v,      population phased vcf file \n"
"  -g,      graph variant hic file \n"
"  -c,      chromosome name ( chr4 ) or contig name depending on bam \n"
"  -n,      id string for output files \n"
"\n";

///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
static void parse_pop_phaser_options( int argc, char** argv ) {
	bool die = false; //bool vcf_load = false; //bool cov_load = false;
	//if(argc <= 2) { die = true; }
        if (string(argv[1]) == "help" || string(argv[1]) == "--help") {
                std::cerr << "\n" << POP_PHASER_USAGE_MESSAGE;
                exit(1);
        }
        for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
                std::istringstream arg(optarg != NULL ? optarg : "");
                switch (c) {
		case 'h': die = true; break;
		case 'v': arg >> opt::input_vcf_file; break;
                case 'g': arg >> opt::input_graph_file; break;
                case 'c': arg >> opt::chr_choice; break;
                case 'n': arg >> opt::id_string; break;
                case 'w': arg >> opt::window_size; break;
                case 'e': arg >> opt::cutoff; break;

                }
        }
	if (die) {
	  std::cerr << "\n" << POP_PHASER_USAGE_MESSAGE;
	  exit(1);
	}
	//parse_region_string( opt::chr_choice );
        cout << endl;
        cout << "############### running linker pop ############### " << endl;
        cout << "== chromosome === " << opt::chr_choice << endl;
        cout << "== vcf file  === " << opt::input_vcf_file << endl;
        cout << "== input hic graph  === " << opt::input_graph_file << endl;
        cout << "== id string  === " << opt::id_string << endl;
        cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx " << endl;
        cout << endl;
};

/////////////////////////////////////////////////////////////////
void init_hic_pop_matrix( map_matrix<double>& diff_matrix, map_matrix<int>& num_matrix, coord_dictionary& pdict, variant_graph& hic_vgraph, int centromere_pos ) {
        map_matrix<int> temp_total_matrix(pdict.num_paired);
        map_matrix<int> count_matrix(pdict.num_paired);
        map_matrix<int> n_plus_matrix(pdict.num_paired);
        map_matrix<int> n_minus_matrix(pdict.num_paired);
        // ############################### fill block matrix #################################
        for (int i = 0; i < pdict.num_paired; i++) {

            int pos1 = pdict.double_positions[i];
            std::string arm1 = "p"; if ( pos1 > centromere_pos ) { arm1 = "q"; }
            int hap1 = pdict.haplotype[i]; //int pop_hap1 = pdict.pop_hap[i];

		if (hap1 != 0) {
                	std::string ref_hash = pdict.ref_handle[i];
                	std::string alt_hash = pdict.alt_handle[i];
			//cout << pos1 << "\t" << arm1 << "\t" <<  hap1 << "\t" << ref_hash << "\t" << alt_hash << endl;
			if ( hic_vgraph.find(ref_hash) != hic_vgraph.end() ) {
                                int alt1 = 1;
				for (auto& it : hic_vgraph[ref_hash].connections) {
					int pos2 = std::stoi(split_string_first(it.first,"_",0));
                        		std::string hic_ref_base = split_string_first(split_string_first(it.first,"_",1),"_",0);
                        		std::string hic_read_base = split_string_first(split_string_first(it.first,"_",1),"_",1);
					int alt2 = -1; if (hic_ref_base == hic_read_base) { alt2 = 1; }
					int j = pdict.pos_index[pos2];
					int hap2 = pdict.haplotype[j];
					if ( hap2 != 0 ) {
            					std::string arm2 = "p"; if ( pos2 > centromere_pos ) { arm2 = "q"; }
            					//int link_hap = alt1*hap1*alt2*hap2;
            					int link_hap = alt1*alt2;
                                                if ( link_hap == 1 ) { n_plus_matrix.add_to(i,j,1);  }//n_plus_matrix.add_to(j,i,1); }
                                                if ( link_hap == -1 ) { n_minus_matrix.add_to(i,j,1); }//n_minus_matrix.add_to(j,i,1); }
            					temp_total_matrix.add_to(i,j,link_hap);
            					temp_total_matrix.add_to(j,i,link_hap);
            					count_matrix.add_to(i,j,1);
            					count_matrix.add_to(j,i,1);
            					//cout << i << "\t" << pos1 << "\t" << alt1 << "\t" << j << "\t" << pos2 << "\t" << alt2 << "\t" << temp_total_matrix(i,j) << endl;
					}
				}
			}
			if ( hic_vgraph.find(alt_hash) != hic_vgraph.end() ) {
                                int alt1 = -1;
				for (auto& it : hic_vgraph[alt_hash].connections) {
					int pos2 = std::stoi(split_string_first(it.first,"_",0));
					std::string hic_ref_base = split_string_first(split_string_first(it.first,"_",1),"_",0);
                                        std::string hic_read_base = split_string_first(split_string_first(it.first,"_",1),"_",1);
					int alt2 = -1; if (hic_ref_base == hic_read_base) { alt2 = 1; }
					int j = pdict.pos_index[pos2];
					int hap2 = pdict.haplotype[j];	
					if ( hap2 != 0 ) {
    					std::string arm2 = "p"; if ( pos2 > centromere_pos ) { arm2 = "q"; }
    					//int link_hap = alt1*hap1*alt2*hap2;
    					int link_hap = alt1*alt2;
                                        if ( link_hap == 1 ) { n_plus_matrix.add_to(i,j,1); }//n_plus_matrix.add_to(j,i,1); }
                                        if ( link_hap == -1 ) { n_minus_matrix.add_to(i,j,1); }//n_minus_matrix.add_to(j,i,1); }
    					temp_total_matrix.add_to(i,j,link_hap);
    					temp_total_matrix.add_to(j,i,link_hap);
    					count_matrix.add_to(i,j,1);
    					count_matrix.add_to(j,i,1);
    					//cout << i << "\t" << pos1 << "\t" << alt1 << "\t" << j << "\t" << pos2 << "\t" << alt2 << "\t" << temp_total_matrix(i,j) << endl;
					}
				}
			}
		}
	}

    double eps0;
    double eps_sum = (double) 0;
    int eps_count = 0;
    for (int i=0; i < pdict.num_paired; i++) {
        for (auto const &ent1 : count_matrix.mat[i]){
            auto const &m = ent1.first;
            double eps = (double) min(n_plus_matrix(i,m), n_minus_matrix(i,m)) / count_matrix(i,m); //(n_plus_matrix(i,m) + n_minus_matrix(i,m));
            eps_sum = eps_sum+eps;
            eps_count++;
        }
    }
    eps0 = (double)eps_sum/eps_count;
    cout << "eps_sum" << "\t" << eps_sum << endl;
    cout << "eps_count" << "\t" << eps_count << endl;
    cout << "eps0" << "\t" << eps0  << endl;           


	for (int i=0; i < pdict.num_paired; i++) {
        for (auto const &ent1 : count_matrix.mat[i]){
            auto const &m = ent1.first;
            double prefactor = (double)temp_total_matrix(i,m)*abs(temp_total_matrix(i,m))/count_matrix(i,m);
            /*
            double eps_ij = (double) min(n_plus_matrix(i,m), n_minus_matrix(i,m)) / (n_plus_matrix(i,m) + n_minus_matrix(i,m));
            double eps = max(eps0, eps_ij);
            double prefactor = (double)(n_plus_matrix(i,m) - n_minus_matrix(i,m))*log((1-eps)/eps);
            */
            //cout << i << "\t" << m << "\t" << prefactor << endl;
            diff_matrix.set_val(i,m,prefactor);
            diff_matrix.set_val(m,i,prefactor);
        }
            }
	num_matrix = temp_total_matrix;
};

/////////////////////////////////////////////////////////////////
void match_pdict( vcf_vector& vvec, coord_dictionary& pdict ) {
        int count = 0;
        for (int i = 0; i < pdict.num_paired; i++) { pdict.ref_handle.push_back("ref"); pdict.alt_handle.push_back("alt"); }
        for (int i = 0; i < vvec.size(); i++) {
		if ( pdict.pos_index.find(vvec[i].pos) != pdict.pos_index.end() ) {
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

///////////////////////////////////////////////////////////////////////////////////////////////////////
void run_pop_phaser( int argc, char** argv ) {
        parse_pop_phaser_options(argc, argv);
        std::string pophapsolutionFile = "./output/pop_hap_solution_" + opt::id_string + "_" + opt::chr_choice + ".dat";
        cout << "- output hap file: " << pophapsolutionFile << endl;
        cout << endl;
        std::unordered_map<std::string,int> hg38_pq = {{"chr1",123400000},{"chr2",93900000},{"chr3",90900000},{"chr4",50000000},{"chr5",48800000},{"chr6",59800000},{"chr7",60100000},{"chr8",45200000},{"chr9",43000000},{"chr10",39800000},{"chr11",53400000},{"chr12",35500000},{"chr13",17700000},{"chr14",17200000},{"chr15",19000000},{"chr16",36800000},{"chr17",25100000},{"chr18",18500000},{"chr19",26200000},{"chr20",28100000},{"chr21",12000000},{"chr22",15000000},{"chrX",61000000},{"chrY",10400000}};
        //std::unordered_map<std::string,int> hg19_pq = {{"chr1",125000000},{"chr2",93300000},{"chr3",91000000},{"chr4",50400000},{"chr5",48400000},{"chr6",61000000},{"chr7",59900000},{"chr8",45600000},{"chr9",49000000},{"chr10",40200000},{"chr11",53700000},{"chr12",35800000},{"chr13",17900000},{"chr14",17600000},{"chr15",19000000},{"chr16",36600000},{"chr17",24000000},{"chr18",17200000},{"chr19",26500000},{"chr20",27500000},{"chr21",13200000},{"chr22",14700000},{"chrX",60600000},{"chrY",12500000}};
        int centromere_pos = hg38_pq[opt::chr_choice];
        coord_dictionary pdict;
        read_graph hic_rgraph;
        variant_graph hic_vgraph;
	vcf_vector vvec;
	hic_link hlink;
	bool paired;
	//#################### load_vcf_file ##########################
	vvec = load_vcf_file_pop( opt::input_vcf_file, opt::chr_choice );  // chromosome // opt::chr_choice
        //#################### start of code ##########################
	cout << "loading graph " << endl;
	read_variant_graph_file( opt::input_graph_file, opt::chr_choice, hic_vgraph, hic_rgraph );	
	cout << "linking graph " << endl;
	//#############################################################
        link_hashes( hic_vgraph, hic_rgraph );
	calc_coverage_unique_hash( hic_vgraph );
	//#############################################################
	for (auto& it : hic_vgraph) { 
		it.second.paired = true; 
	}
    initialize_pdict( hic_vgraph, pdict, paired );
	///////////// create global datastructures for haplotype solver
    map_matrix<int> num_matrix(pdict.num_paired);
    map_matrix<int> count_matrix(pdict.num_paired);
	map_matrix<double> diff_matrix(pdict.num_paired);
	////////////////////////////////////////////////////////////////
 	
	///////////////////////////////////////////////////////////////
	match_pdict( vvec, pdict );
    //pdict.pop_hap_switch(100, 0.0005);
	pdict.hap_pop_initialization();
	//pdict.hap_random_initialization();
	init_hic_pop_matrix( diff_matrix, num_matrix, pdict, hic_vgraph, centromere_pos );

        solver_recursive_pop( hic_vgraph, pdict, num_matrix, diff_matrix, opt::window_size, opt::cutoff );
	
	///////////////////////////////////////////////////////////////
	///////////// write haplotype output
	cout << " solver finished " << endl;
        write_pop_hap_solution( hic_vgraph, pophapsolutionFile, pdict, opt::chr_choice );
        return;
};





/////////////////////////////////////////////////////////////////
//void create_hic_links_pop( variant_graph& hic_vgraph, hic_link& hlink ) {
//        for (auto& it : hic_vgraph) { it.second.unique_hash(); }
//        int i=0;
//        for (auto& it : hic_vgraph) {
//                int pos1 = it.second.pos;
//                std::string var1 = it.first;   //cout << pos1 << "\t" << var1 << "\t" << it.second.unique_total << endl;
//                if (it.second.unique_total > 0) {
//                for (auto it2 : it.second.connections) {  //int pos2 = 100;  //var_dict[it2.first].pos;
//                        int pos2 = hic_vgraph[it2.first].pos;
//                        int depth = it2.second;
//                        std::string var2 = it2.first;
                        //cout << pos1 << "\t" << var1 << "\t" << pos2 << "\t" << var2 << "\t" << depth << endl;
//                        hlink.add_link(pos1,pos2,var1,var2,depth);  //anchor_pos1.push_back(pos1);
//                        i++;
//                } }
//        }
//        hlink.nlinks = i;
//}













                //}
                //ptrdiff_t j = get_index_var( pdict.sorted_paired_positions, vvec[i].pos );
                //if ( j < pdict.num_paired) {
                        //cout << i << " " << vvec[i].pos << " " << vvec[i].pop_hap << " " << vvec[i].ref_base << " " << vvec[i].var_base <<  " " << j << " " << pdict.sorted_paired_positions[j] << " " << pdict.pop_hap[j] << endl;



        //prune_graph( hic_vgraph );
        //#############################################################
        //write_link_network(vgraph,output_network_file,opt::chr_choice);
        //initialize_solver(vgraph,pdict,diff_matrix,num_matrix_second);


        //std::map<std::string,int> chr_str_map = {{"chr1", 1}, {"chr2", 2}, {"chr3", 3}, {"chr4", 4},{"chr5", 5}, {"chr6", 6}, {"chr7", 7}, {"chr8", 8},{"chr9", 9}, {"chr10", 10}, {"chr11", 11}, {"chr12", 12},{"chr13", 13}, {"chr14", 14}, {"chr15", 15}, {"chr16", 16},{"chr17", 17}, {"chr18", 18}, {"chr19", 19}, {"chr20", 20},{"chr21", 21}, {"chr22", 22}, {"chrX", 23}, {"chrY", 24}};
        //int chromosome = chr_str_map[opt::chr_choice];


