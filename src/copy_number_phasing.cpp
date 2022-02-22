#include "copy_number_phasing.h"

void initialize_copy_num_map( cn_map& chromosome_map, coord_dictionary& pdict, variant_graph& vgraph, int binsize, std::vector<int>& bin_array ) {
        int max_position = *std::max_element(pdict.sorted_paired_positions.begin(),pdict.sorted_paired_positions.end());
        int k = 0;
        for (int j = 0; j < max_position; j++) {
                if (j % binsize == 0) { //cout << j << endl;
                        cn_bin loop_bin;
                        loop_bin.initialize(j);
                        chromosome_map[k] = loop_bin;
			bin_array.push_back(k);
                        k += 1;
                }
        }
        ////////////////////////////////////////
        for (int j = 0; j < pdict.num_paired; j++) {
                if ( pdict.reload_bool[j] ) {
                        int bloc = floor(pdict.sorted_paired_positions[j]/binsize);
                        int loop_pos = pdict.sorted_paired_positions[j];
                        std::string het_string = pdict.paired_dict[loop_pos][0];
                        std::string rbase = vgraph[het_string].ref_base;
                        std::string vbase = vgraph[het_string].var_base;
                        char rbase_char = rbase[0];
                        char vbase_char = vbase[0];
                        int rnum = vgraph[het_string].base_dict[rbase_char];
                        int vnum = vgraph[het_string].base_dict[vbase_char];
                        //cout << pdict.sorted_paired_positions[j] << " " << pdict.reload_bool[j] << " " << bloc << " " << pdict.haplotype[j] << " " << pdict.deltaE[j] << " " << pdict.block[j] << "  " << rbase << "  " << vbase << " " << rnum << "  " <<  vnum << " " << endl;
                        chromosome_map[bloc].add_het(loop_pos,pdict.haplotype[j],pdict.block[j],rbase,vbase,rnum,vnum);
                }
        }
	for (auto i : bin_array) { chromosome_map[i].cov_calc(); }
	return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
double calc_switch( cn_map& chromosome_map, int index1, int index2 ) {
	double hapA_one,hapB_one,hapA_two,hapB_two;
	if ( chromosome_map[index1].switch_spin == 1 ) { hapA_one = chromosome_map[index1].hapA_avg; hapB_one = chromosome_map[index1].hapB_avg; }
	else {						 hapB_one = chromosome_map[index1].hapA_avg; hapA_one = chromosome_map[index1].hapB_avg; }
	if ( chromosome_map[index2].switch_spin == 1 ) { hapA_two = chromosome_map[index2].hapA_avg; hapB_two = chromosome_map[index2].hapB_avg; }
	else {						 hapB_two = chromosome_map[index2].hapA_avg; hapA_two = chromosome_map[index2].hapB_avg; }
	double diffA = abs(hapA_two - hapA_one);
	double diffB = abs(hapB_two - hapB_one);
	double diffB_switch = abs(hapA_two - hapB_one);
	double diffA_switch = abs(hapB_two - hapA_one);
	double min_func = std::min(diffA,diffB);
	double min_func_switch = std::min(diffA_switch,diffB_switch);
	double min_switch_prob = min_func - min_func_switch;
	return min_switch_prob;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
double calc_switch_merge( cn_map& chromosome_map, std::vector<int> merged_index1, std::vector<int> merged_index2 ) {
        double hapA_one_sum=0.0,hapB_one_sum=0.0;
	double hapA_two_sum=0.0,hapB_two_sum=0.0;
	int nbins_one = merged_index1.size();
	int nbins_two = merged_index2.size();
	for (auto index1 : merged_index1) {
        	if ( chromosome_map[index1].switch_spin == 1 ) { hapA_one_sum += chromosome_map[index1].hapA_avg; hapB_one_sum += chromosome_map[index1].hapB_avg; }
        	else {                                           hapB_one_sum += chromosome_map[index1].hapA_avg; hapA_one_sum += chromosome_map[index1].hapB_avg; }	
	}
	for (auto index2 : merged_index2) {
                if ( chromosome_map[index2].switch_spin == 1 ) { hapA_two_sum += chromosome_map[index2].hapA_avg; hapB_two_sum += chromosome_map[index2].hapB_avg; }
                else {                                           hapB_two_sum += chromosome_map[index2].hapA_avg; hapA_two_sum += chromosome_map[index2].hapB_avg; }
	}
	double hapA_one = hapA_one_sum/(double)nbins_one;
	double hapB_one = hapB_one_sum/(double)nbins_one;
	double hapA_two = hapA_two_sum/(double)nbins_two;
	double hapB_two = hapB_two_sum/(double)nbins_two;
        double diffA = abs(hapA_two - hapA_one);
        double diffB = abs(hapB_two - hapB_one);
        double diffB_switch = abs(hapA_two - hapB_one);
        double diffA_switch = abs(hapB_two - hapA_one);
        double min_func = std::min(diffA,diffB);
        double min_func_switch = std::min(diffA_switch,diffB_switch);
        double min_switch_prob = min_func - min_func_switch;
        return min_switch_prob;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void cn_phasing( cn_map& chromosome_map, std::vector<int>& bin_array, std::vector<int>& good_bins, std::vector<int>& rec_bins, std::vector<int>& merged_bins ) {
        //int k = 0;
        for (auto i : bin_array) {
		//cout << i << "\t" << chromosome_map[i].num_blocks << "\t" << chromosome_map[i].num_hets << endl; 
                if (chromosome_map[i].num_blocks == 1 && chromosome_map[i].num_hets > (minimum_snp-1)) {
			//good_bins.push_back(k); 
			good_bins.push_back(i); 
		}  //k++;
        }
	cout << "subset good bins: " << good_bins.size() << endl; 
	cn_phase_between_blocks( chromosome_map, good_bins ); 
	//for (int j = 0; j < 3; j++) { 
	//	cn_phase_loop( chromosome_map, good_bins ); 
	//}
	for (int j = 0; j < 1; j++) { merge_bins( chromosome_map, good_bins, merged_bins ); }
	////////////////////////////// updated bin printing
	//rescue_bad_bins( chromosome_map, bin_array, good_bins, rec_bins );
}
///////////////////////////////////////////////////////////////////////////////////////////////////
void cn_phasing_greg( coord_dictionary& pdict, variant_graph& vgraph, double delta, int range ) {
    int prior = 10;

	std::vector<int> depth;

    std::vector<int> gbin;
	int last_pos = 1;
   
	for (int i = 2; i < pdict.num_paired-2; i+=4) {
        
		if ( abs(pdict.sorted_paired_positions[i-2]-pdict.sorted_paired_positions[i+2] ) <= 10000) {
			for (int j = 0; j < 5; j++) {gbin.push_back(1);}
		}
		else { for (int j = 0; j < 5; j++) {gbin.push_back(0); } }
	}

    //for ( int i = 0; i < gbin.size(); i++ ) { if (gbin[i] == 0) { pdict.haplotype[i] = 0; } }
    
        get_depth(pdict, 100, prior, 1.0, vgraph, delta, depth); 

        //sliding_block_flip( pdict, range, prior, 1.0, vgraph, delta, gbin, false );
        //sliding_block_flip( pdict, range, prior, 1.0, vgraph, delta, gbin, true );
        //for ( int i = 0; i < 200; i++ ) { sliding_block_flip( pdict, range, prior, 1.0, vgraph, delta, gbin, true ); }
        //for ( int i = 0; i < 200; i++ ) { bin_flip_recursive_cn(pdict, range/2, prior, 1.0, vgraph, delta, depth, true); }
        //for ( int i = 0; i < 200; i++ ) { bin_flip_recursive_cn(pdict, range/5, prior, 1.0, vgraph, delta, depth, true); }
        //for ( int i = 0; i < 200; i++ ) { bin_flip_recursive_cn(pdict, range/10, prior, 1.0, vgraph, delta, depth, true); }
        //for ( int i = 0; i < 200; i++ ) { bin_flip_recursive_cn(pdict, 10, prior, 1.0, vgraph, delta, depth, true); }
        for ( int i = 0; i < 1800; i++ ) { bin_flip_recursive_cn(pdict, 200, prior, 1.0, vgraph, delta, depth, true); }
        sliding_block_flip( pdict, range, prior, 1.0, vgraph, delta, gbin, true );
        //for ( int i = 0; i < 20; i++ ) { bin_flip_recursive_cn(pdict, 300, prior, 1.0, vgraph, delta, depth, true); }
        //for ( int i = 0; i < 20; i++ ) { bin_flip_recursive_cn(pdict, 3000, prior, 1.0, vgraph, delta, depth, true); }
        //for ( int i = 0; i < 20; i++ ) { bin_flip_recursive_cn(pdict, 5000, prior, 1.0, vgraph, delta, depth, true); }
        


    
        //for ( int j = 0; j < 200; j++ ) { cout << j << endl; sliding_block_flip( pdict, range, prior, 1.0, vgraph, delta, gbin, true ); }
        //for ( int j = 0; j < 200; j++ ) { cout << j << endl; sliding_block_flip( pdict, 2*range, prior, 1.0, vgraph, delta, gbin, true ); }
        
        //for ( int j = 0; j < 10; j++ ) { cout << j << endl; sliding_block_flip( pdict, range/2, prior, 1.0, vgraph, delta, gbin, false ); }
        //for ( int j = 0; j < 10; j++ ) { cout << j << endl; sliding_block_flip( pdict, range/10, prior, 1.0, vgraph, delta, gbin, false ); }
        //for ( int j = 0; j < 500; j++ ) { bin_flip_recursive_cn( pdict, range/2, prior, 1.0, vgraph, 0.3, depth, true ); }
        //for ( int j = 0; j < 500; j++ ) { bin_flip_recursive_cn( pdict, range, prior, 1.0, vgraph, 0.3, depth, true ); }
        //for ( int j = 0; j < 500; j++ ) { bin_flip_recursive_cn( pdict, range/4, prior, 1.0, vgraph, 0.3, depth, true ); }
        //for ( int i = 0; i < 5; i++ ) { bin_flip_recursive_cn( pdict, range/10, prior, 1.0, vgraph, 0.3, depth, false ); }
        //for ( int i = 0; i < 5; i++ ) { bin_flip_recursive_cn( pdict, 10, prior, 1.0, vgraph, 0.3, depth, false ); }
        //for ( int i = 0; i < 5; i++ ) { bin_flip_recursive_cn( pdict, 5, prior, 1.0, vgraph, 0.3, depth, false ); }
        //for ( int j = 0; j < 5; j++ ) { bin_flip_recursive_cn( pdict, range, prior, 1.0, vgraph, delta, depth, false ); }
       

        

        //for ( int j = 0; j < 50; j++ ) { bin_flip_recursive_cn(pdict, 10, prior, 1.0, vgraph, delta, depth); }
        //for ( int j = 0; j < 50; j++ ) { bin_flip_recursive_cn(pdict, 5, prior, 1.0, vgraph, delta, depth); }

        //for (int i = 0; i < 2; i++ ) { block_flip_recursive_cn(pdict, range, prior, 1.0, vgraph, delta, gbin); }
        //for (int i = 0; i < 2; i++ ) { block_flip_recursive_cn(pdict, 10*range, prior, 1.0, vgraph, delta, gbin); }
        //for (int i = 0; i < 2; i++ ) { block_flip_recursive_cn(pdict, 100*range, prior, 1.0, vgraph, delta, gbin); }
        //for (int i = 0; i < 5; i++ ) { block_flip_recursive_cn(pdict, 5000, prior, 1.0, vgraph, delta, gbin); }

/*     
    //for ( int j = 0; j < 2; j++ ) {
        for ( int i = 0; i < 200; i++ ) { bin_flip_recursive_cn(pdict, 3*range/2, prior, 1.0, vgraph, delta, gbin); }
        //for ( int i = 0; i < 200; i++ ) { bin_flip_recursive_cn(pdict, range, prior, 1.0, vgraph, delta, gbin); }
        //for ( int i = 0; i < 500; i++ ) { bin_flip_recursive_cn(pdict, range/2, prior, 1.0, vgraph, delta, gbin); }
        //for ( int i = 0; i < 200; i++ ) { bin_flip_recursive_cn(pdict, 4*range, prior, 1.0, vgraph, delta, gbin); }
        //for ( int i = 0; i < 500; i++ ) { bin_flip_recursive_cn(pdict, range/5, prior, 1.0, vgraph, delta, gbin); }
        for (int i = 0; i < 5; i++ ) { block_flip_recursive_cn(pdict, 5000, prior, 1.0, vgraph, delta, gbin); }

    //}
*/
};
////////////////////////////////////////////////////////////////////////////////////////////////////
void block_flip_recursive_cn( coord_dictionary& pdict, int range, int& prior, double cutoff, variant_graph& vgraph, double delta, std::vector<int> gbin) {
       
        double min_switchE = 0.0;
        vector<int> min_haplotype = pdict.haplotype;
        vector<int> loop_haplotype = pdict.haplotype;

        double first_switchE = 0.0;
        double loop_switchE;
        double save_first = 0.0;
        int last_pos;

        int current_win;
        int min_pos;

        int loop_pos;
        std::string het_string;
        std::string rbase; 
        std::string vbase;
        char rbase_char;
        char vbase_char;
        int rnum;
        int vnum;
        int cnAx = 0; 
        int cnAy=0; 
        int cnBx=0; 
        int cnBy=0; 
        int cnAx_=0; 
        int cnAy_=0; 
        int cnBx_=0; 
        int cnBy_=0;

        //Calculating copy-number in the first block
        for ( int i = 0; i < range; i++ ) {
            loop_pos = pdict.sorted_paired_positions[i];
            het_string = pdict.paired_dict[loop_pos][0];
            rbase = vgraph[het_string].ref_base;
            vbase = vgraph[het_string].var_base;
            rbase_char = rbase[0];
            vbase_char = vbase[0];
            rnum = vgraph[het_string].base_dict[rbase_char];
            vnum = vgraph[het_string].base_dict[vbase_char];
          
            if (pdict.haplotype[i] == 1) {  cnAx += vnum;   cnBx += rnum;   cnAx_ += rnum;  cnBx_ += vnum;  }
            if (pdict.haplotype[i] == -1) {  cnAx += rnum;   cnBx += vnum;   cnAx_ += vnum;  cnBx_ += rnum;  }    
        }

        int last = range;
        //Calculating copy-number in the second block
        for ( int i = range; i < last+range; i++ ) {
            loop_pos = pdict.sorted_paired_positions[i];
            het_string = pdict.paired_dict[loop_pos][0];
            rbase = vgraph[het_string].ref_base;
            vbase = vgraph[het_string].var_base;
            rbase_char = rbase[0];
            vbase_char = vbase[0];
            rnum = vgraph[het_string].base_dict[rbase_char];
            vnum = vgraph[het_string].base_dict[vbase_char];

            if (pdict.haplotype[i] == 1) {  cnAy += vnum;   cnBy += rnum;   }
            if (pdict.haplotype[i] == -1) {  cnAy += rnum;   cnBy += vnum;   }
        }

        double Ax_avg = (double)cnAx/range;
        double Ay_avg = (double)cnAy/range;
        double diff1 = Ax_avg-Ay_avg;

        double Bx_avg = (double)cnBx/range;
        double By_avg = (double)cnBy/range;

        double diff2 = Bx_avg-By_avg;

        double x_high = max(Ax_avg, Bx_avg);
        double y_high = max(Ay_avg, By_avg);
        double x_low = min(Ax_avg, Bx_avg);
        double y_low = min(Ay_avg, By_avg);

        
        int d1;
        
        if ( min(abs(diff1), abs(diff2)) > delta ) { d1 = 1; }
        else { d1 = 0; }

        double diff3 = Ax_avg-By_avg;
        double diff4 = Bx_avg - Ay_avg;

        int d2;
        if ( min(abs(diff3), abs(diff4)) > delta ) { d2 = 1; }
        else { d2 = 0; }

        int exy = d1-d2;

        //cout << "x_high" << "\t" << x_high << "\t" << "x_low" << "\t" << x_low << "\t" << "y_high" << "\t" << y_high << 
        //"\t" << "y_low" << "\t" << y_low << "\t" << "exy" << "\t" << exy << endl;


        last_pos = range-1;
        
        pdict.switchE[0] = exy;

        //Calculating the subsequent block-flipping energies
        for ( int i = 1; i < (pdict.num_paired-2*range); i++ ) {
            last_pos = i+range-1;
            loop_haplotype[i-1] = pdict.haplotype[i-1];
            loop_haplotype[last_pos] = -1*pdict.haplotype[last_pos];          
            

            //Adding the new links from the last variant in the block
            loop_pos = pdict.sorted_paired_positions[last_pos];
            het_string = pdict.paired_dict[loop_pos][0];
            rbase = vgraph[het_string].ref_base;
            vbase = vgraph[het_string].var_base;
            rbase_char = rbase[0];
            vbase_char = vbase[0];
            rnum = vgraph[het_string].base_dict[rbase_char];
            vnum = vgraph[het_string].base_dict[vbase_char];

            if (pdict.haplotype[last_pos] == 1) {
                cnAx += vnum;
                cnBx += rnum;
                cnAx_ += rnum;
                cnBx_ += vnum;
                        }       
            if (pdict.haplotype[last_pos] == -1) {
                cnAx += rnum;
                cnBx += vnum;
                cnAx_ += vnum;
                cnBx_ += rnum;
            }

            loop_pos = pdict.sorted_paired_positions[last_pos+range];
            het_string = pdict.paired_dict[loop_pos][0];
            rbase = vgraph[het_string].ref_base;
            vbase = vgraph[het_string].var_base;
            rbase_char = rbase[0];
            vbase_char = vbase[0];
            rnum = vgraph[het_string].base_dict[rbase_char];
            vnum = vgraph[het_string].base_dict[vbase_char];

            if (pdict.haplotype[last_pos+range] == 1) {
                cnAy += vnum;
                cnBy += rnum;
                cnAy_ += rnum;
                cnBy_ += vnum;
                        }       
            if (pdict.haplotype[last_pos+range] == -1) {
                cnAy += rnum;
                cnBy += vnum;
                cnAy_ += vnum;
                cnBy_ += rnum;
            }

            //Removing terms that are now within the new block
            loop_pos = pdict.sorted_paired_positions[i-1];
            het_string = pdict.paired_dict[loop_pos][0];
            rbase = vgraph[het_string].ref_base;
            vbase = vgraph[het_string].var_base;
            rbase_char = rbase[0];
            vbase_char = vbase[0];
            rnum = vgraph[het_string].base_dict[rbase_char];
            vnum = vgraph[het_string].base_dict[vbase_char];
          
            if (pdict.haplotype[i-1] == 1) {
                cnAx += -vnum;
                cnBx += -rnum;

                cnAx_ += -rnum;
                cnBx_ += -vnum;
                        }       
            if (pdict.haplotype[i-1] == -1) {
                cnAx += -rnum;
                cnBx += -vnum;

                cnAx_ += -vnum;
                cnBx_ += -rnum;
            }

            loop_pos = pdict.sorted_paired_positions[i-1+range];
            het_string = pdict.paired_dict[loop_pos][0];
            rbase = vgraph[het_string].ref_base;
            vbase = vgraph[het_string].var_base;
            rbase_char = rbase[0];
            vbase_char = vbase[0];
            rnum = vgraph[het_string].base_dict[rbase_char];
            vnum = vgraph[het_string].base_dict[vbase_char];

            if (pdict.haplotype[i-1+range] == 1) {
                cnAy += -vnum;
                cnBy += -rnum;

                cnAy_ += -rnum;
                cnBy_ += -vnum;
                        }       
            if (pdict.haplotype[i-1+range] == -1) {
                cnAy += -rnum;
                cnBy += -vnum;

                cnAy_ += -vnum;
                cnBy_ += -rnum;
            }

            Ax_avg = (double)cnAx/range;
            Ay_avg = (double)cnAy/range;
            diff1 = Ax_avg-Ay_avg;

            Bx_avg = (double)cnBx/range;
            By_avg = (double)cnBy/range;

            diff2 = Bx_avg-By_avg;

            if ( min(abs(diff1), abs(diff2)) > delta ) { d1 = 1; }
            else { d1 = 0; }

            diff3 = Ax_avg-By_avg;
            diff4 = Bx_avg - By_avg;

            if ( min(abs(diff3), abs(diff4)) > delta ) { d2 = 1; }
            else { d2 = 0; }

            exy = d1-d2;
    
            pdict.switchE[i] = exy;

        }

        for(int i = range; i < pdict.num_paired; i++) {
            if (pdict.switchE[i] == 1) { //&& pdict.sorted_paired_positions[i+range]-pdict.sorted_paired_positions[i]<5000000) {
                for ( int j = 0; j<i+range; j++) {int initial_hap = pdict.haplotype[j]; pdict.haplotype[j] = -1*initial_hap;}
                i+=(range-1);
            }
        }

       
};

/////////////////////////////////////////////////////////////////////////////////////////////////////
void bin_flip_recursive_cn( coord_dictionary& pdict, int range, int& prior, double cutoff, variant_graph& vgraph, double delta, std::vector<int> gbin, bool max ) {
       
        int loop_pos;
        std::string het_string;
        std::string rbase; 
        std::string vbase;
        char rbase_char;
        char vbase_char;
        int rnum;
        int vnum;
        int cnAx=0; 
        int cnAy=0; 
        int cnBx=0; 
        int cnBy=0; 
       
        std::vector<int> A_cn;
        std::vector<int> B_cn;
        std::vector<int> depth;

        double Ax_avg;
        double Bx_avg;
        double Ay_avg;
        double By_avg;
        double diff1;
        double diff2;
        int d1;
        int d2;
        double diff3;
        double diff4;
        int exy;

        for(int i = 0; i < pdict.num_paired; i++) { pdict.switchE[i] = 0; }

        //Calculating copy-number in each block
        for( int j = 0; j < pdict.num_paired-range; j+=range ) {

            cnAx = 0; 
            cnAy = 0; 
            cnBx = 0; 
            cnBy = 0; 

            for ( int i = j; i < j + range; i++ ) {
                loop_pos = pdict.sorted_paired_positions[i];
                het_string = pdict.paired_dict[loop_pos][0];
                rbase = vgraph[het_string].ref_base;
                vbase = vgraph[het_string].var_base;
                rbase_char = rbase[0];
                vbase_char = vbase[0];
                rnum = vgraph[het_string].base_dict[rbase_char];
                vnum = vgraph[het_string].base_dict[vbase_char];
              
                if (pdict.haplotype[i] == 1) { cnAx += vnum;   cnBx += rnum; }
                if (pdict.haplotype[i] == -1) { cnAx += rnum;   cnBx += vnum; }    
            }

            A_cn.push_back(cnAx);
            B_cn.push_back(cnBx);
            depth.push_back(cnAx+cnBx);
        }

        
        //Calculating the switching energies
        for (int i = 0; i < A_cn.size(); i++ ) {


            //Ax_avg = (double)100*A_cn[i]/(depth[i]*range);
            Ax_avg = (double)A_cn[i]/depth[i];
            //Bx_avg = (double)100*B_cn[i]/(depth[i]*range);
            Bx_avg = (double)B_cn[i]/depth[i];

            for (int j = i+1; j < A_cn.size(); j++ ) {
                //Ay_avg = (double)100*A_cn[j]/(depth[j]*range);
                Ay_avg = (double)A_cn[j]/depth[j];
                //By_avg = (double)100*B_cn[j]/(depth[j]*range);
                By_avg = (double)B_cn[j]/depth[j];
                
                diff1 = Ax_avg-Ay_avg;
                diff2 = Bx_avg-By_avg;

                //cout << "Ax_avg " << Ax_avg << "\t" << "Bx_avg " << Bx_avg << "\t" << "Ay_avg " << Ay_avg << "\t" << "By_avg " << By_avg << endl;
               
                if ( min(abs(diff1), abs(diff2)) > delta) { d1 = 1; }
                else { d1 = 0; }
                diff3 = Ax_avg-By_avg;
                diff4 = Bx_avg-Ay_avg;
                if ( min(abs(diff3), abs(diff4)) > delta ) { d2 = 1; }
                else { d2 = 0; }
                exy = d1-d2;
                //if (gbin[range*i] == 1 && gbin[range*j] == 1) { 
                    pdict.switchE[range*i]+=exy; 
                    pdict.switchE[range*j]+=exy;
                //}
            }
        }

        int initial_hap;

        if (max == true) {
            //Finding the max switching energy
            int max_switchE = 0;
            int max_pos = 0;
            for (int i = 0; i < pdict.num_paired; i++) {
                if ( pdict.switchE[i] > max_switchE ) { max_switchE = pdict.switchE[i]; max_pos = i; }
            }

            cout << max_pos << "\t" << max_switchE << endl;
            if (max_switchE > 0) {
                for(int i = max_pos; i < max_pos+range; i++) {
                //for(int i = max_pos; i < pdict.num_paired; i++) {
                    initial_hap = pdict.haplotype[i]; 
                    pdict.haplotype[i] = -1*initial_hap;
                }
            }
        } else {

            for ( int j = 0; j < pdict.num_paired; j++ ) {
                if ( pdict.switchE[j] > 0 ) {
                    for(int i = j; i < j+range; i++) {
                        initial_hap = pdict.haplotype[i]; 
                        pdict.haplotype[i] = -1*initial_hap;
                    }
                    j+=range;
                }
            }

        }
};
//////////////////////////////////////////////////////////////////////////////////////////////////////
void sliding_block_flip( coord_dictionary& pdict, int range, int& prior, double cutoff, variant_graph& vgraph, double delta, std::vector<int> gbin, bool max ) {
       
        double min_switchE = 0.0;
        vector<int> min_haplotype = pdict.haplotype;
        vector<int> loop_haplotype = pdict.haplotype;

        double first_switchE = 0.0;
        double loop_switchE;
        double save_first = 0.0;
        int last_pos;

        int current_win;
        int min_pos;

        int loop_pos;
        std::string het_string;
        std::string rbase; 
        std::string vbase;
        char rbase_char;
        char vbase_char;
        int rnum;
        int vnum;
        int cnAx = 0; 
        int cnBx = 0; 

        std::vector<int> A_cn;
        std::vector<int> B_cn;
        std::vector<int> depth;

        double Ax_avg;
        double Bx_avg;
        double Ay_avg;
        double By_avg;
        double diff1;
        double diff2;
        double diff3;
        double diff4;
        int d1;
        int d2;
        int exy;

        for (int i = 0; i < pdict.num_paired; i++) { pdict.switchE[i] = 0; }

        //Calculating copy-number in the first block
        for ( int i = 0; i < range; i++ ) {
            loop_pos = pdict.sorted_paired_positions[i];
            het_string = pdict.paired_dict[loop_pos][0];
            rbase = vgraph[het_string].ref_base;
            vbase = vgraph[het_string].var_base;
            rbase_char = rbase[0];
            vbase_char = vbase[0];
            rnum = vgraph[het_string].base_dict[rbase_char];
            vnum = vgraph[het_string].base_dict[vbase_char];
          
            if (pdict.haplotype[i] == 1) {  cnAx += vnum;   cnBx += rnum; }
            if (pdict.haplotype[i] == -1) {  cnAx += rnum;   cnBx += vnum; }    
        }

        A_cn.push_back(cnAx);
        B_cn.push_back(cnBx);
        depth.push_back(cnAx + cnBx);

        int last = range;
        last_pos = range-1;


        //Calculating the subsequent block-flipping energies
        for ( int i = 1; i < (pdict.num_paired-range); i++ ) {
            last_pos = i+range-1;
            loop_haplotype[i-1] = pdict.haplotype[i-1];
            loop_haplotype[last_pos] = -1*pdict.haplotype[last_pos];          
            
            //Adding the new links from the last variant in the block
            loop_pos = pdict.sorted_paired_positions[last_pos];
            het_string = pdict.paired_dict[loop_pos][0];
            rbase = vgraph[het_string].ref_base;
            vbase = vgraph[het_string].var_base;
            rbase_char = rbase[0];
            vbase_char = vbase[0];
            rnum = vgraph[het_string].base_dict[rbase_char];
            vnum = vgraph[het_string].base_dict[vbase_char];

            if (pdict.haplotype[last_pos] == 1) {
                cnAx += vnum;
                cnBx += rnum;
                        }       
            if (pdict.haplotype[last_pos] == -1) {
                cnAx += rnum;
                cnBx += vnum;
            }

            //Removing terms that are now within the new block
            loop_pos = pdict.sorted_paired_positions[i-1];
            het_string = pdict.paired_dict[loop_pos][0];
            rbase = vgraph[het_string].ref_base;
            vbase = vgraph[het_string].var_base;
            rbase_char = rbase[0];
            vbase_char = vbase[0];
            rnum = vgraph[het_string].base_dict[rbase_char];
            vnum = vgraph[het_string].base_dict[vbase_char];
          
            if (pdict.haplotype[i-1] == 1) {
                cnAx += -vnum;
                cnBx += -rnum;

            }       
            if (pdict.haplotype[i-1] == -1) {
                cnAx += -rnum;
                cnBx += -vnum;
            }

            A_cn.push_back(cnAx);
            B_cn.push_back(cnBx);
            depth.push_back(cnAx+cnBx);
        }

        //Calculating the switching energies
        for (int i = 0; i < A_cn.size(); i++ ) {

            Ax_avg = (double)A_cn[i]/depth[i];
            //Bx_avg = (double)100*B_cn[i]/(depth[i]*range);
            Bx_avg = (double)B_cn[i]/depth[i];

            for (int j = i+1; j < A_cn.size(); j++ ) {
                //Ay_avg = (double)100*A_cn[j]/(depth[j]*range);
                Ay_avg = (double)A_cn[j]/depth[j];
                //By_avg = (double)100*B_cn[j]/(depth[j]*range);
                By_avg = (double)B_cn[j]/depth[j];
                
                diff1 = Ax_avg-Ay_avg;
                diff2 = Bx_avg-By_avg;

                if ( min(abs(diff1), abs(diff2)) > delta) { d1 = 1; }
                else { d1 = 0; }
                diff3 = Ax_avg-By_avg;
                diff4 = Bx_avg-Ay_avg;
                if ( min(abs(diff3), abs(diff4)) > delta ) { d2 = 1; }
                else { d2 = 0; }
                exy = d1-d2;
                //if (gbin[range*i] == 1 && gbin[range*j] == 1) { 
                    pdict.switchE[i]+=exy; 
                    pdict.switchE[j]+=exy;
                //}
            }
        }

        int initial_hap;

        if (max == true) {
            //Finding the max switching energy
            int max_switchE = 0;
            int max_pos = 0;
            for (int i = 0; i < pdict.num_paired; i++) {
                if ( pdict.switchE[i] > max_switchE ) { max_switchE = pdict.switchE[i]; max_pos = i; }
            }

            cout << max_pos << "\t" << max_switchE << endl;
            if (max_switchE > 0) {
                for(int i = max_pos; i < max_pos+range; i++) {
                //for(int i = max_pos; i < pdict.num_paired; i++) {
                    initial_hap = pdict.haplotype[i]; 
                    pdict.haplotype[i] = -1*initial_hap;
                }
            }
        } else {

            for ( int j = 0; j < pdict.num_paired; j++ ) {
                if ( pdict.switchE[j] > 0 ) {
                    for(int i = j; i < j+range; i++) {
                        initial_hap = pdict.haplotype[i]; 
                        pdict.haplotype[i] = -1*initial_hap;
                    }
                    j+=(range-1);
                }
            }

        }  
};
/////////////////////////////////////////////////////////////////////////////////////////////////////
void get_depth( coord_dictionary& pdict, int range, int& prior, double cutoff, variant_graph& vgraph, double delta, std::vector<int>& depth ) {
       
        int loop_pos;
        std::string het_string;
        std::string rbase; 
        std::string vbase;
        char rbase_char;
        char vbase_char;
        int rnum;
        int vnum;
        int cnAx=0; 
        int cnAy=0; 
        int cnBx=0; 
        int cnBy=0; 
       
        std::vector<int> A_cn;
        std::vector<int> B_cn;
        //std::vector<int> depth;

        double Ax_avg;
        double Bx_avg;
        double Ay_avg;
        double By_avg;
        double diff1;
        double diff2;
        int d1;
        int d2;
        double diff3;
        double diff4;
        int exy;

        //Calculating copy-number in each block
        for( int j = 0; j < pdict.num_paired-100; j+=100 ) {

            cnAx = 0; 
            cnAy = 0; 
            cnBx = 0; 
            cnBy = 0; 

            for ( int i = j; i < j + 100; i++ ) {
                loop_pos = pdict.sorted_paired_positions[i];
                het_string = pdict.paired_dict[loop_pos][0];
                rbase = vgraph[het_string].ref_base;
                vbase = vgraph[het_string].var_base;
                rbase_char = rbase[0];
                vbase_char = vbase[0];
                rnum = vgraph[het_string].base_dict[rbase_char];
                vnum = vgraph[het_string].base_dict[vbase_char];
              
                if (pdict.haplotype[i] == 1) { cnAx += vnum;   cnBx += rnum; }
                if (pdict.haplotype[i] == -1) { cnAx += rnum;   cnBx += vnum; }    
            }

            A_cn.push_back(cnAx);
            B_cn.push_back(cnBx);
            depth.push_back(cnAx+cnBx);
        }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void bin_flip_all( coord_dictionary& pdict, int range, int& prior, double cutoff, variant_graph& vgraph, double delta, std::vector<int> gbin ) {
       
        int loop_pos;
        std::string het_string;
        std::string rbase; 
        std::string vbase;
        char rbase_char;
        char vbase_char;
        int rnum;
        int vnum;
        int cnAx = 0; 
        int cnAy=0; 
        int cnBx=0; 
        int cnBy=0; 
       
        std::vector<int> A_cn;
        std::vector<int> B_cn;

        double Ax_avg;
        double Bx_avg;
        double Ay_avg;
        double By_avg;
        double diff1;
        double diff2;
        int d1;
        int d2;
        double diff3;
        double diff4;
        int exy;

        //Calculating copy-number in each block
        for( int j = 0; j < pdict.num_paired-range; j+=range ) {

            cnAx = 0; 
            cnAy = 0; 
            cnBx = 0; 
            cnBy = 0; 

            for ( int i = j; i < j + range; i++ ) {
                loop_pos = pdict.sorted_paired_positions[i];
                het_string = pdict.paired_dict[loop_pos][0];
                rbase = vgraph[het_string].ref_base;
                vbase = vgraph[het_string].var_base;
                rbase_char = rbase[0];
                vbase_char = vbase[0];
                rnum = vgraph[het_string].base_dict[rbase_char];
                vnum = vgraph[het_string].base_dict[vbase_char];
              
                if (pdict.haplotype[i] == 1) { cnAx += vnum;   cnBx += rnum; }
                if (pdict.haplotype[i] == -1) { cnAx += rnum;   cnBx += vnum; }    
            }

            A_cn.push_back(cnAx);
            B_cn.push_back(cnBx);
        }

        //Calculating the switching energies
        for (int i = 0; i < A_cn.size(); i++ ) {
            Ax_avg = (double)A_cn[i]/range;
            Bx_avg = (double)B_cn[i]/range;
            for (int j = i+1; j < A_cn.size(); j++ ) {
                Ay_avg = (double)A_cn[j]/range;
                By_avg = (double)B_cn[j]/range;
                diff1 = Ax_avg-Ay_avg;
                diff2 = Bx_avg-By_avg;
               
                if ( min(abs(diff1), abs(diff2)) > delta) { d1 = 1; }
                else { d1 = 0; }
                diff3 = Ax_avg-By_avg;
                diff4 = Bx_avg-Ay_avg;
                if ( min(abs(diff3), abs(diff4)) > delta ) { d2 = 1; }
                else { d2 = 0; }
                exy = d1-d2;
                //if (gbin[range*i] == 1 & gbin[range*j] == 1) { 
                    pdict.switchE[range*i]+=exy; 
                    pdict.switchE[range*j]+=exy;
                //}
            }
        }

        int initial_hap;

        for( int i = 0; i < pdict.num_paired; i+=range ) {
            if ( pdict.switchE[i] > 0 ) {
                cout << pdict.sorted_paired_positions[i] << endl;
                for ( int j = i; j < i+range; j++ ) {
                    initial_hap = pdict.haplotype[j]; 
                    pdict.haplotype[j] = -1*initial_hap;
                }
            }
        }
       
};
////////////////////////////////////////////////////////////////////////////////////////////////////
void cn_block_flip( coord_dictionary& pdict, variant_graph& vgraph ) {
	for (int j = 0; j < pdict.num_paired; j++) {
		int loop_pos = pdict.sorted_paired_positions[j];
        std::string het_string = pdict.paired_dict[loop_pos][0];
        std::string rbase = vgraph[het_string].ref_base;
        std::string vbase = vgraph[het_string].var_base;
        char rbase_char = rbase[0];
        char vbase_char = vbase[0];
        int rnum = vgraph[het_string].base_dict[rbase_char];
        int vnum = vgraph[het_string].base_dict[vbase_char];
        cout << pdict.sorted_paired_positions[j] << " " << rbase << "  " << vbase << " " << rnum << "  " <<  vnum << " " << endl;

	}
};
///////////////////////////////////////////////////////////////////////////////////////////////////
void rescue_bad_bins( cn_map& chromosome_map, std::vector<int>& bin_array, std::vector<int>& good_bins, std::vector<int>& rec_bins ) {
	std::vector<int> block_torescue;
	for (auto ind : good_bins) { int block_id = chromosome_map[ind].het_block[0]; block_torescue.push_back(block_id); }
	std::unordered_map<int,std::vector<int>> block_map_inclusive;
	std::unordered_map<int,int> block_map_flipped;
	std::unordered_map<int,bool> within_good_bin;
        for (auto i : bin_array) {
		if (chromosome_map[i].num_blocks == 1) { int block_id = chromosome_map[i].het_block[0];  block_map_inclusive[block_id].push_back(i); } 
                bool in_good_bins = std::find(std::begin(good_bins), std::end(good_bins), i) != std::end(good_bins);
		within_good_bin[i] = in_good_bins; if (in_good_bins) { block_map_flipped[i] = chromosome_map[i].switch_spin; }
        }
	std::vector<int> recovered_bins; int num_recovered_hets = 0;
	for ( auto& ind : block_map_inclusive ) {
		int num_good_bins = 0; int sum_flips = 0;
		for ( auto& ind2 : ind.second ) { if (within_good_bin[ind2]) { num_good_bins += 1; sum_flips += block_map_flipped[ind2]; } }
		if ( num_good_bins == abs(sum_flips) && num_good_bins != 0 ) {  //cout << "# in loop " << num_good_bins << "  " << sum_flips << endl;	
                	for ( auto& ind2 : ind.second ) {
				recovered_bins.push_back(ind2);
                        	if ( !within_good_bin[ind2] ) {
					num_recovered_hets += chromosome_map[ind2].num_hets; 
					if ( sum_flips < 0 ) { chromosome_map[ind2].flip_hap(); }
				}
			}
		}
	}
	std::set<int> previous(good_bins.begin(), good_bins.end());	
	std::set<int> recovered(recovered_bins.begin(), recovered_bins.end());
	std::vector<int> new_bins; set_union(previous.begin(),previous.end(),recovered_bins.begin(),recovered_bins.end(),back_inserter(new_bins));
	std::set<int> nbins(new_bins.begin(),new_bins.end());
	std::vector<int> rbins(nbins.begin(),nbins.end());
	std::sort(rbins.begin(),rbins.end(),std::less<int>());    // std::greater
	int num_total_hets = 0; for (auto i : bin_array) { num_total_hets += chromosome_map[i].num_hets; }
	int num_good_hets = 0; for (auto i : good_bins) { num_good_hets += chromosome_map[i].num_hets; }
	int num_new_hets = 0; for (auto i : rbins) { num_new_hets += chromosome_map[i].num_hets; rec_bins.push_back(i); }
	cout << "# total hets " << num_total_hets << endl;
	cout << "# good bin hets " << num_good_hets << endl;
	cout << "# recovered hets " << num_recovered_hets << endl;
	cout << "# combined hets " << num_new_hets << endl;	
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void flip_to_end( int next , cn_map& chromosome_map, std::vector<int>& good_bins ) {
	for (auto ind : good_bins) {
		if (ind >= next) {
			chromosome_map[ind].flip_hap();	
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void cn_phase_loop( cn_map& chromosome_map, std::vector<int>& good_bins ) {  //, std::vector<int>& bin_array
	int nbins = good_bins.size(); int m=0;
	for (auto ind : good_bins) {   //if ( m < 10) {
		if (nbins > (m+1)) {
			double min_switch_prob_plusone = 0.0;
			double min_switch_prob_plustwo = 0.0;
			double min_switch_prob_plusthree = 0.0;
			if (m >= 1 && m < nbins - 2) {
				int plusone = good_bins[m+2]; int minusone = good_bins[m-1];  //cout << "plus one" << " " << plusone << "  " << minusone << endl;
				min_switch_prob_plusone = calc_switch(chromosome_map,minusone,plusone);
			}
			if (m >= 2 && m < nbins - 3) {
				int plustwo = good_bins[m+3]; int minustwo = good_bins[m-2]; //cout << "plus two"<< " " << plustwo << "  " << minustwo << endl;
				min_switch_prob_plustwo = calc_switch(chromosome_map,minustwo,plustwo);
			}
                        if (m >= 3 && m < nbins - 4) {
                                int plusthree = good_bins[m+4]; int minusthree = good_bins[m-3]; //cout << "plus two"<< " " << plustwo << "  " << minustwo << endl;
                                min_switch_prob_plusthree = calc_switch(chromosome_map,minusthree,plusthree);
                        }
			//cout << m << " " << nbins << endl;
			int next = good_bins[m+1];
            		double min_switch_prob = calc_switch(chromosome_map,ind,next);
        		double hapA_one,hapB_one,hapA_two,hapB_two;
        		if ( chromosome_map[ind].switch_spin == 1 ) {    hapA_one = chromosome_map[ind].hapA_avg; hapB_one = chromosome_map[ind].hapB_avg; }
        		else {                                           hapB_one = chromosome_map[ind].hapA_avg; hapA_one = chromosome_map[ind].hapB_avg; }
        		if ( chromosome_map[next].switch_spin == 1 ) {   hapA_two = chromosome_map[next].hapA_avg; hapB_two = chromosome_map[next].hapB_avg; }
        		else {                                           hapB_two = chromosome_map[next].hapA_avg; hapA_two = chromosome_map[next].hapB_avg; } 
			cout << ind << "\t" << next << "\t" << min_switch_prob << "\t" << min_switch_prob_plusone << "\t" << min_switch_prob_plustwo << "\t" << min_switch_prob_plusthree << "\t" << hapA_one << "\t" << hapB_one << "\t" << hapA_two << "\t" << hapB_two << endl;
                        if ( min_switch_prob > min_switch ) { 
				cout << " switch: " << ind << endl;  
				flip_to_end(next,chromosome_map,good_bins);
			}
			//else {
			//	if (min_switch_prob_plusone > min_switch_two && min_switch_prob_plustwo > min_switch_two && min_switch_prob_plusthree > min_switch_two ) {
			//		cout << " switch plussthree: " << ind << endl;  
			//		flip_to_end(next,chromosome_map,good_bins);
			//	}
			//}
		}
		m+=1;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void cn_phase_between_blocks( cn_map& chromosome_map, std::vector<int>& good_bins ) {  //, std::vector<int>& bin_array
        int nbins = good_bins.size(); int m=0;
        for (auto ind : good_bins) {   //if ( m < 10) {
                if (nbins > (m+1)) { 
                        int next = good_bins[m+1];
                        double min_switch_prob = calc_switch(chromosome_map,ind,next);
                        double hapA_one,hapB_one,hapA_two,hapB_two;
			int block_ind = chromosome_map[ind].het_block[0];
			int block_next = chromosome_map[next].het_block[0];
			//cout << ind << "\t" << next << endl;
			//int num_block_ind = chromosome_map[ind].num_blocks;
			//int num_block_next = chromosome_map[next].num_blocks;
                        if ( chromosome_map[ind].switch_spin == 1 ) {    hapA_one = chromosome_map[ind].hapA_avg; hapB_one = chromosome_map[ind].hapB_avg; }
                        else {                                           hapB_one = chromosome_map[ind].hapA_avg; hapA_one = chromosome_map[ind].hapB_avg; }
                        if ( chromosome_map[next].switch_spin == 1 ) {   hapA_two = chromosome_map[next].hapA_avg; hapB_two = chromosome_map[next].hapB_avg; }
                        else {                                           hapB_two = chromosome_map[next].hapA_avg; hapA_two = chromosome_map[next].hapB_avg; }
                        cout << ind << "\t" << next << "\t" << "\t" << block_ind << "\t" << block_next << "\t" << min_switch_prob << "\t" << hapA_one << "\t" << hapB_one << "\t" << hapA_two << "\t" << hapB_two << endl;
                        if ( min_switch_prob > first_pass_switch && block_ind != block_next) {
                                cout << " switch: " << ind << endl;
                                flip_to_end(next,chromosome_map,good_bins);
                        }
                }
                m+=1;
        }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
void merge_bins( cn_map& chromosome_map, std::vector<int>& good_bins, std::vector<int>& merged_bins ) {
	int nbins = good_bins.size();
	int m = 0;
	double merge_energy_cutoff = -0.4;
	int array_index = 0;
	int merged_index = good_bins[0];
	std::map<int,mbin> merged_vec;	
	merged_vec[array_index].initialize_mbin(merged_index);
        for (auto ind : good_bins) {   //if ( m < 10) {
                if (nbins > (m+1)) {
			int next = good_bins[m+1];
			double min_switch_prob = calc_switch(chromosome_map,ind,next);
			cout << merged_index << " " << ind << " " << next << "  " << min_switch_prob << endl;
			if ( min_switch_prob < merge_energy_cutoff ) {
				merged_vec[array_index].set_energies(ind,min_switch_prob);	
			}
			else {
				array_index += 1;
				merged_index = next;
				merged_vec[array_index].initialize_mbin(next);
			}
		}
		m += 1;
		merged_bins.push_back(array_index);
	}
	for (auto const& x : merged_vec) {
		//cout << x.first << " " << x.second.bin_index << " " << x.second.num_mbins << endl;
		for (auto const& bin : x.second.merged_bins) {
			cout << bin << endl;
		}
	}
	std::vector<int> good_merged;
	for (auto const& x : merged_vec) {
		if (x.second.num_mbins > 2) {
			good_merged.push_back(x.first);
		}
	}
	int nbins_merged = good_merged.size();
	int l = 0;
	cout << "merged block phasing" << endl;
        for (auto ind : good_merged) {   //if ( m < 10) {
		//cout << ind << endl;
                if (nbins_merged > (l+1)) {
			int next = good_merged[l+1];
			//merged_vec[ind].merged_bins;  //merged_vec[next].merged_bins;
			int min_merge_bin = *std::min_element( merged_vec[next].merged_bins.begin(), merged_vec[next].merged_bins.end() );
			double min_switch_prob = calc_switch_merge(chromosome_map,merged_vec[ind].merged_bins,merged_vec[next].merged_bins);
			cout << merged_vec[ind].bin_index << " " << merged_vec[ind].num_mbins << " " << min_switch_prob << endl;
                        if ( min_switch_prob > first_pass_switch ) {
                                cout << " switch: " << ind << endl;
                                flip_to_end(min_merge_bin,chromosome_map,good_bins);
                        }
		}
		l += 1;
	}	
}







                //for ( auto& ind2 : ind.second ) {
                //   cout << ind.first << "\t" << ind2 << "\t" << within_good_bin[ind2] << "\t" << chromosome_map[ind2].start_pos << "\t" << chromosome_map[ind2].num_hets << endl;
                //}



//int num_goodbin_hets = 0;
                                //else { num_goodbin_hets += chromosome_map[ind2].num_hets; }


                //cout << "position " << chromosome_map[i].start_pos << " good bins: " << in_good_bins << " num blocks: " << chromosome_map[i].num_blocks << " switched: " << chromosome_map[i].switch_spin << " number hets: " << chromosome_map[i].num_hets << endl;













                //for (int k = 0; k < chromosome_map[i].num_hets; k++) {
                //        cout << "het pos " << chromosome_map[i].het_pos[k] << " " << chromosome_map[i].het_block[k] << endl;
                //      int hblock = chromosome_map[i].het_block[k]; int hpos = chromosome_map[i].het_pos[k];
                //      block_map[hblock].push_back(hpos);
                //}




		//cout << "######## " << num_good_bins << "  " << sum_flips << endl; 




                //for ( auto& ind2 : ind.second ) {
                //      if ( chromosome_map[ind2].num_hets > (minimum_snp-1) ) {
                //      }
                //      if ( num_good_bins == abs(sum_flips) && num_good_bins != 0 ) {
                //      }
                //}
                //cout << "######## " << num_good_bins << "  " << sum_flips << endl;





//min_switch_prob_plusone > min_switch_two


                        //if ( min_switch_prob_plusone > min_switch && min_switch_prob_plustwo > min_switch ) {
                        //if ( min_switch_prob > -min_switch ) {

                        //double min_switch_prob_back = -min_switch_prob;
                        //cout << ind << "\t" << next << "\t" << min_switch_prob << "\t" << min_switch_prob_plusone << "\t" << min_switch_prob_plustwo << endl;


                        //if ( min_switch_prob > min_switch ) { //cout << "switch block " << endl;
                        //        chromosome_map[next].flip_hap();
                        //}  //cout << endl;



  //cout << ind << "\t" << next << "\t" << min_switch_prob << "\t" << min_switch_prob_plusone << "\t" << min_switch_prob_plustwo << hapA_one << "\t" << hapB_one << hapA_two << "\t" << hapB_two << endl;

                        //cout << hapA_one << "\t" << hapB_one << "\t num blocks " << chromosome_map[ind].num_blocks << "\t num hets " << chromosome_map[ind].num_hets << endl;
                        //cout << hapA_two << "\t" << hapB_two << "\t num blocks " << chromosome_map[next].num_blocks << "\t num hets " << chromosome_map[next].num_hets << endl;

//cout << i << " " << chromosome_map[i].num_hets << endl;

