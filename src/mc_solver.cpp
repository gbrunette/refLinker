#include "mc_solver.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////
void solver( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<double> diff_matrix, map_matrix<int> num_matrix_second ) {
	int t = clock();
        static const std::size_t length = pdict.num_paired;
        map_matrix<int> num_matrix_third(length);
        cout << "========" << endl;
        cout << "monte carlo solver: " << endl;
        length_cutoff_nmatrix( pdict, diff_matrix, num_matrix_second, num_matrix_third );
        vector<bool> flip_maxima(length); 
	for (int j = 0; j < pdict.num_paired; j++) { flip_maxima[j] = true; }
        for (int i = 0; i < solver_loops; i++) {
                cout << "******************* solver loop     " << i << endl;
                t = clock();
                single_spin_flip_map( pdict, diff_matrix, num_matrix_third );  // num_matrix_second
                if (i > 0) { get_flip_positions( pdict, flip_maxima ); }
                t = clock() - t;
                cout << "spin flip --------- time: " << t << endl;
                t = clock();
                cout << "-- finished flip spins     " << endl;
                block_flip_map( pdict, diff_matrix, num_matrix_third, flip_maxima );  // num_matrix_second
                t = clock() - t;
                cout << "block flip -------- time: " << t << endl;
                t = clock();
                cout << "-- finished block flip     " << endl;
                energy_sum_min_max( pdict );
        }
        set_switch_energy( pdict, flip_maxima );   /// create a new haplotype block class and find the segments
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void solver_recursive( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<double> diff_matrix, map_matrix<int> num_matrix_second ) {
        int t = clock();
	static const std::size_t length = pdict.num_paired;
	map_matrix<int> num_matrix_third(length);
        cout << "========" << endl;
        cout << "monte carlo solver: " << endl;
	length_cutoff_nmatrix( pdict, diff_matrix, num_matrix_second, num_matrix_third );
        //vector<bool> flip_maxima(pdict.num_paired);
        //for (int j = 0; j < pdict.num_paired; j++) { flip_maxima[j] = true; }
        for (int i = 0; i < solver_loops; i++) {
                cout << "******************* solver loop     " << i << endl;
                t = clock();
		single_spin_flip_recursive( pdict, diff_matrix, num_matrix_third );
                t = clock() - t;
                cout << "spin flip --------- time: " << t << endl;
                t = clock();
		block_flip_recursive( pdict, diff_matrix, num_matrix_third );
                t = clock() - t;
                cout << "block flip -------- time: " << t << endl;
                energy_sum_min_max( pdict );
        }
        //set_switch_energy(pdict,flip_maxima);   /// create a new haplotype block class and find the segments
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void solver_recursive_pop( std::unordered_map<std::string,variant_node>& var_dict, coord_dictionary& pdict, map_matrix<int> num_matrix, map_matrix<double> diff_matrix, int window_size, double cutoff ) {
    
    static const std::size_t length = pdict.num_paired;

    int bkp = 0;
    
    map_matrix<int> num_matrix_cutoff(length);
    length_cutoff_nmatrix( pdict, diff_matrix, num_matrix, num_matrix_cutoff );

    int max_gap = 0;
    for (int i = 0; i < pdict.num_paired; i++) { 
            int loop_gap = pdict.sorted_paired_positions[i+1]-pdict.sorted_paired_positions[i];
            //if ( loop_gap > max_gap ) { max_gap = loop_gap; }
    }
    int max_loop = 2000;
    int count, prior_count;
    int window = window_size;
    if ( max_gap > 1e7 && pdict.num_paired < 1e5) { window = 17000; }
    int pop_weight = 4e5;
    int max_window = 0;
    double min_switch = 0;
    double loop_min;
    bool global;
    
    
  
    for (int i = 0; i < 50; i++) {
        cout << "global loop:" << "\t" << i << endl;
        switchE_recursive_pop( pdict, num_matrix, pop_weight, diff_matrix, cutoff, bkp, true, loop_min );      
    }

    count = 1e6;
    prior_count = count;
    //double cutoff = -10.0;


    for (int k = 0; k < 1000; k++) {
        max_window = 0;
        min_switch = 0;
        loop_min = 0;
        window = window_size;
        global = false;


        //switchE_recursive_pop( pdict, num_matrix, pop_weight, diff_matrix, cutoff, bkp, false, loop_min );
        //if (loop_min < min_switch) { min_switch = loop_min; max_window = window; global = true; } 

        while ( window > 0) {
            //cout << window << "\t" << "loop:" << "\t" << k << endl;
            block_flip_recursive_var_range( pdict, num_matrix, window, count, pop_weight, diff_matrix, 0, bkp, false, loop_min );
            if (loop_min < min_switch) { min_switch = loop_min; max_window = window; global = false; } 
            if (window > 200 ) { window = window - 20; }
            else { window = window - 1; }
            count = 1e6;
            prior_count = count;
        }

        window = (int)0.5*pdict.num_paired;

        while ( window > window_size ) {
            block_flip_recursive_var_range( pdict, num_matrix, window, count, pop_weight, diff_matrix, 0, bkp, false, loop_min );
            if (loop_min < min_switch) { min_switch = loop_min; max_window = window; global = false; } 
            window = (int)window/2;
        }
        if ( min_switch > cutoff ) { k = 10000; }
        cout << "loop" << "\t" << k << "\t" << "max_win" << "\t" << max_window << "\t" << min_switch << "\t" << global << endl; 
        if ( global == false ) { block_flip_recursive_var_range( pdict, num_matrix, max_window, count, pop_weight, diff_matrix, 0, bkp, true, min_switch ); }
        else { switchE_recursive_pop( pdict, num_matrix, pop_weight, diff_matrix, -10.0, bkp, true, min_switch ); }
    }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
void solver_recursive_hic( block_dictionary& bdict, std::vector<int> hic_limit_loop, map_matrix<int> block_matrix ) {
        int t = clock();
        //static const std::size_t length = pdict.num_paired;
        //map_matrix<int> num_matrix_third(length);
        cout << "========" << endl;
        cout << "monte carlo solver:" << endl;
        //length_cutoff_nmatrix( pdict, diff_matrix, num_matrix_second, num_matrix_third );
        //vector<bool> flip_maxima(pdict.num_paired);
        //vector<int> block_deltaE(length_subset);
        //for (int j = 0; j < pdict.num_paired; j++) { flip_maxima[j] = true; }
        //block_phasing_quick_flip( bdict );
	//
	vector<int> graph_flip_list;
	//block_phasing_quick_flip( bdict, graph_flip_list );
	for( int k = 0; k < hic_limit_loop.size(); k++ ) {
		int hlimit = hic_limit_loop[k];
		bdict.scale_matrix( hlimit, block_matrix );
        	for (int i = 0; i < solver_loops_hic; i++) {
                	cout << "******************* solver_loop" << "\t" << i << "\tscale_loop\t" << hlimit << endl;
                	t = clock();
                	single_spin_flip_recursive_hic( bdict );
        		//block_phasing_quick_flip( bdict, graph_flip_list );
                	t = clock() - t;
                	cout << "spin flip --------- time: " << t << endl;
                	t = clock();
			block_flip_brute_force_hic( bdict );  //block_flip_recursive_hic( bdict );
                	t = clock() - t;
                	cout << "block flip -------- time: " << t << endl;
                	energy_sum_min_max_hic( bdict );
        	} 
	}
        //set_switch_energy(pdict,flip_maxima);   /// create a new haplotype block class and find the segments
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_flip_recursive_var_range( coord_dictionary& pdict, map_matrix<int>& nmatrix, int range, int& prior, int pop_weight, map_matrix<double> diff_matrix, double cutoff, int& bkp, bool switch_hap, double& switch_min ) {
       
        double factor;
        int dist = 0;
        double min_switchE = 0.0;
        vector<int> min_haplotype = pdict.haplotype;
        vector<int> loop_haplotype = pdict.haplotype;


        double first_switchE = 0.0;
        double loop_switchE;
        double save_first = 0.0;
        double save_pop;
        int last_pos;

        int current_win;
        int pair_dist;
        int min_pos;

        int bin_distance;
        double bin_density;

        //Creating the first block flip
        for ( int i = 0; i < range; i++ )
        {
            loop_haplotype[i] = -1*pdict.haplotype[i];
        }

        last_pos = range-1;
        //Calculating the switching energy of the first block flip
        for ( int i = 0; i < range; i++ )
        {
            for ( auto const &ent1 : nmatrix.mat[i] ) 
            {
                auto const &m = ent1.first;
                if (m > last_pos)
                {
                    first_switchE += 2.0*pdict.haplotype[i]*pdict.haplotype[m]*diff_matrix(i,m);
                }
            }   
        }

//Adding the population term
        dist = abs(pdict.sorted_paired_positions[last_pos]-pdict.sorted_paired_positions[last_pos+1]);
        factor = (double)dist/pop_weight;
        save_pop = 2.0*exp(-1.0*factor)*pdict.haplotype[last_pos]*pdict.haplotype[last_pos+1]*pdict.pop_hap[last_pos]*pdict.pop_hap[last_pos+1];

        first_switchE += save_pop;

        pdict.switchE[0] = first_switchE;

        if (first_switchE < min_switchE) {
            min_switchE = first_switchE;
            min_haplotype = loop_haplotype;
        }
        loop_switchE = first_switchE;

        //Calculating the subsequent block-flipping energies
        for ( int i = 1; i < (pdict.num_paired-range); i++ ) {

            last_pos = i+range-1;

            loop_haplotype[i-1] = pdict.haplotype[i-1];
            loop_haplotype[last_pos] = -1*pdict.haplotype[last_pos];          
            

            //Subtracting the EAGLE2 term from the previous block
            loop_switchE += -1.0*save_pop;

            //Adding the new links from the last variant in the block
            for ( auto const &ent1 : nmatrix.mat[last_pos]) {
                auto const &m = ent1.first;

                 

                if ( m > last_pos ) {
                    loop_switchE += 2.0*pdict.haplotype[last_pos]*pdict.haplotype[m]*diff_matrix(last_pos,m);
                }
                if ( m < i ) {
                    loop_switchE += 2.0*pdict.haplotype[last_pos]*pdict.haplotype[m]*diff_matrix(last_pos,m);    
                }
                //Removing terms that are now within the new block
                if (m >= (i-1) ) {
                    if (m < last_pos) {
                        loop_switchE += -2.0*pdict.haplotype[last_pos]*pdict.haplotype[m]*diff_matrix(last_pos,m);
                    }
                }
            }

            //Taking care of the previous variant
            for (auto const &ent1 : nmatrix.mat[i - 1]){
                auto const &m = ent1.first;
                if ( m > last_pos) {
                    loop_switchE += -2.0*pdict.haplotype[i-1]*pdict.haplotype[m]*diff_matrix(i-1, m);
                }
                if ( m < (i-1)) {
                    loop_switchE += -2.0*pdict.haplotype[i-1]*pdict.haplotype[m]*diff_matrix(i-1,m);
                }
                if (m >= i) {
                    if (m < last_pos) {
                        loop_switchE += 2.0*pdict.haplotype[i-1]*pdict.haplotype[m]*diff_matrix(i-1,m);
                    }                
                }
            }

            //Adding population term
            dist = abs(pdict.sorted_paired_positions[last_pos]-pdict.sorted_paired_positions[last_pos+1]);
            factor = (double)dist/pop_weight;
            save_pop = 2.0*exp(-1.0*factor)*pdict.haplotype[last_pos]*pdict.haplotype[last_pos+1]*pdict.pop_hap[last_pos]*pdict.pop_hap[last_pos+1];
            loop_switchE += save_pop;
            pdict.switchE[i] = loop_switchE;

            bin_distance = pdict.sorted_paired_positions[last_pos] - pdict.sorted_paired_positions[i];
            bin_density = (double)range/bin_distance;
            
            if ( loop_switchE < min_switchE ) { 
            //if ( loop_switchE < min_switchE && i!= bkp ) { //&& max_win < 5e5) {  
                min_switchE = loop_switchE;
                min_haplotype = loop_haplotype;
                min_pos = i;
                switch_min = min_switchE;
                //cout << pdict.sorted_paired_positions[i+range-1]-pdict.sorted_paired_positions[i] << endl;
                //cout << min_switchE << endl;
            }
        }
        //cout << "min switching energy" << "\t" << min_switchE << "\t" << min_pos << endl;
        if ( min_switchE < -10.0 && switch_hap == true ) { pdict.haplotype = min_haplotype; prior = prior-1; bkp = min_pos; cout << "switch pos" << "\t" << min_pos << endl; }
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////

static void flipE( coord_dictionary& pdict, map_matrix<int>& nmatrix, int pop_weight, map_matrix<double> diff_matrix ) {
    int nsf = 0;
    double minE = 0.0;
    int min_pos = 0;
    int site_diff;
    double factor;
    for (int i = 0; i < pdict.num_paired; i++) {
            int init_i = pdict.haplotype[i]; //int final_i = -1*pdict.haplotype[i];
            double initial_energy = 0.0;     //double final_energy = 0.0;
            
            for (auto const &ent1 : nmatrix.mat[i]) {
                auto const &m = ent1.first;       
                //initial_energy += -1.0*init_i*pdict.haplotype[m]*diff_matrix(i,m);
                initial_energy += -1.0*init_i*pdict.haplotype[m]*nmatrix(i,m);
                
                //if (pdict.sorted_paired_positions[i] > 39000000 && pdict.sorted_paired_positions[i] < 40500000) {
                    //cout << "HiC_link" << "\t" << pdict.sorted_paired_positions[i] << "\t"<< pdict.sorted_paired_positions[m] << "\t" << abs(pdict.sorted_paired_positions[i]-pdict.sorted_paired_positions[m]) << "\t" << nmatrix(i,m) << "\t" << endl;
                //}
            }
            if (i < pdict.num_paired-1) {
                site_diff = abs(pdict.sorted_paired_positions[i]-pdict.sorted_paired_positions[i+1]);
                factor = (double)site_diff/pop_weight;
                //COMMENTED IN TESTING 12/01/2020
                initial_energy += -1.0*exp(-1.0*factor)*init_i*pdict.haplotype[i+1]*pdict.pop_hap[i]*pdict.pop_hap[i+1];
            }
            
            double final_energy = -1*initial_energy;
            double diff_energy = final_energy - initial_energy;
            //if (diff_energy < -1.0) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); nsf++; }
        
            if ( diff_energy < minE ) {minE = diff_energy; min_pos = i;} 
            pdict.deltaE[i] = (-1)*diff_energy;
    }
    //cout << "spin flips" << "\t" << nsf << endl;
    //if (minE < 0.0 ) { pdict.haplotype[min_pos] = pdict.haplotype[min_pos]*(-1); }
    for (int i = 0; i < pdict.num_paired; i++) {if (pdict.deltaE[i]<-100.0) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); nsf++; } }
    //cout << "min pos" << "\t" << min_pos << endl;
    cout << "number flips" << "\t" << nsf << endl;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_flip_scaffold( coord_dictionary& pdict, map_matrix<int>& nmatrix, int dist, int& prior, map_matrix<double> diff_matrix ) {
        int nbf = 0;
        double initial_energy = 0.0;
        int site_diff;
        double factor;
        int last_pos;
        
        for (int i = 0; i < pdict.num_paired; i++) {
            for (auto const &ent1 : nmatrix.mat[i]) {
                auto const &m = ent1.first;
                //initial_energy += -1.0*pdict.haplotype[i]*pdict.haplotype[m]*nmatrix(i,m); 
                initial_energy += -1.0*pdict.haplotype[i]*pdict.haplotype[m]*diff_matrix(i,m);    
            }
        }

        double energyA = initial_energy;
        double energyB = initial_energy;
        vector<int> loop_haplotype = pdict.haplotype;

        for (int i = 0; i < pdict.num_paired; i++) {
            vector<int> switched;
            double save_energyB = energyB;
            vector<int> switch_haplotype = loop_haplotype;
            int initial_pos = pdict.sorted_paired_positions[i];
            int j = i+1;
            bool window = true;
            while(window == true && j < pdict.num_paired) { 
                switch_haplotype[j] = pdict.haplotype[j]*(-1);
                switched.push_back(j);
                j+=1;
                if(pdict.sorted_paired_positions[j] > initial_pos+dist) {
                    window = false;
                }
            }
            
            int hap_spin  = loop_haplotype[i];    
            int flip_spin = -1*loop_haplotype[i];
            double spin_i_energy = 0.0;
            double spin_f_energy = 0.0; 
            for (int k = 0; k < switched.size(); k++) {
                for (auto const &ent1 : nmatrix.mat[switched[k]]) {
                    auto const &m = ent1.first;
                    site_diff = min(abs(pdict.sorted_paired_positions[i]-pdict.sorted_paired_positions[m]), 
                        abs(pdict.sorted_paired_positions[j]-pdict.sorted_paired_positions[m]));
                    if ( site_diff > 1e5 ) {
                        //spin_f_energy += -1.0*switch_haplotype[switched[k]]*switch_haplotype[m]*nmatrix(switched[k],m); 
                        spin_f_energy += -1.0*switch_haplotype[switched[k]]*switch_haplotype[m]*diff_matrix(switched[k],m); 
                        //spin_i_energy += -1.0*pdict.haplotype[switched[k]]*pdict.haplotype[m]*nmatrix(switched[k],m); 
                        spin_i_energy += -1.0*pdict.haplotype[switched[k]]*pdict.haplotype[m]*diff_matrix(switched[k],m); 
                    }
                }
            }
            
            double switch_energy = spin_f_energy - spin_i_energy;
            
            if ( switch_energy < 0.0) {
                loop_haplotype = switch_haplotype;
                nbf++;
                pdict.switchE[i] = (-1.0)*switch_energy;
            }
            i = j;
        }
        cout << "number block flips: " << "\t" << nbf << endl;
        if ( nbf < prior ) { pdict.haplotype = loop_haplotype; prior = nbf; }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////
static void switchE_recursive_pop( coord_dictionary& pdict, map_matrix<int>& nmatrix, int pop_weight, map_matrix<double> diff_matrix, double cutoff, int& bkp, bool flip, double& loop_min ) {//, map_matrix<double> diff_matrix, int pop_weight ) {
        int nbf = 0;
        int site_diff;
        int min_pos = 0;
        vector<int> h_plus;
        vector<int> h_minus;
        double initial_energy = 0.0;
        double factor;
        for (int i = 0; i < pdict.num_paired; i++) {
            int h_left = 0;
            int h_right = 0;
            for (auto const &ent1 : nmatrix.mat[i]) { 
                auto const &m = ent1.first;
                if ( pdict.sorted_paired_positions[m] < pdict.sorted_paired_positions[i] ) { h_left = h_left + diff_matrix(i,m)*pdict.haplotype[m]; }
                if ( pdict.sorted_paired_positions[m] > pdict.sorted_paired_positions[i] ) { h_right = h_right + diff_matrix(i,m)*pdict.haplotype[m]; }
            }
            h_minus.push_back(h_left);
            h_plus.push_back(h_right);
        }
        
        for (int i = 0; i < pdict.num_paired; i++) {

            for (auto const &ent1 : nmatrix.mat[i]) {
                auto const &m = ent1.first;
                initial_energy += -1.0*pdict.haplotype[i]*pdict.haplotype[m]*diff_matrix(i,m);
            }
            if (i < pdict.num_paired-1) {
                site_diff = abs(pdict.sorted_paired_positions[i]-pdict.sorted_paired_positions[i+1]);
                factor = (double)site_diff/pop_weight;
                //FOLLOWING COMMENTED IN TESTING 12/01/2020
                initial_energy += -1.0*exp(-1.0*factor)*pdict.haplotype[i]*pdict.haplotype[i+1]*pdict.pop_hap[i]*pdict.pop_hap[i+1];
            }
        }


        
        double final_energy = 0.0;
        vector<int> switch_haplotype = pdict.haplotype;
        for (int j = 1; j < pdict.num_paired; j++) { switch_haplotype[j] = pdict.haplotype[j]*(-1); }
        for (int i = 0; i < pdict.num_paired; i++) {
            for (auto const &ent1 : nmatrix.mat[i]){
                auto const &m = ent1.first;
                final_energy += -1.0*switch_haplotype[i]*switch_haplotype[m]*diff_matrix(i,m);
            }
            if (i < pdict.num_paired-1) {
                site_diff = abs(pdict.sorted_paired_positions[i]-pdict.sorted_paired_positions[i+1]);
                factor = (double)site_diff/pop_weight;
                //FOLLOWING COMMENTED IN TESTING 12/01/2020
                final_energy += -1.0*exp(-1.0*factor)*switch_haplotype[i]*switch_haplotype[i+1]*pdict.pop_hap[i]*pdict.pop_hap[i+1];
            }
        }
        double diffE_last = final_energy - initial_energy;
        double min_switch = diffE_last;
        pdict.switchE[0] = (-1.0)*diffE_last;
        for (int j = 1; j < pdict.num_paired-1; j++) {
            int site_diff1 = abs(pdict.sorted_paired_positions[j] - pdict.sorted_paired_positions[j+1]);
            int site_diff2 = abs(pdict.sorted_paired_positions[j-1] - pdict.sorted_paired_positions[j]);
            double factor1 = (double)site_diff1/pop_weight;
            double factor2 = (double)site_diff2/pop_weight;
            //FOLLOWING LINE COMMENTED IN TESTING 12/01/2020
            diffE_last = diffE_last + pdict.haplotype[j]*(h_plus[j]-h_minus[j]) + pdict.haplotype[j]*pdict.pop_hap[j]*(exp(-1.0*factor1)*pdict.haplotype[j+1]*pdict.pop_hap[j+1] - exp(-1.0*factor2)*pdict.haplotype[j-1]*pdict.pop_hap[j-1]);

            if ( diffE_last < min_switch && j != bkp ) { min_switch = diffE_last; min_pos = j; }
            //if ( diffE_last < min_switch ) { min_switch = diffE_last; min_pos = j; }
            pdict.switchE[j] = (-1.0)*diffE_last;
            
        }
        /*
        for ( int j = 0; j < pdict.num_paired; j++ ) {
            if ((-1.0)*pdict.switchE[j] < cutoff ){ for (int i = j; i < pdict.num_paired; i++) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); } }
        }
        */
        
        loop_min = min_switch; 

        if (min_switch < cutoff && flip == true ) { for (int i = min_pos+1; i < pdict.num_paired; i++) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); } bkp = min_pos; }

    cout << min_pos << endl;
};
////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_phasing_quick_flip( block_dictionary& bdict, vector<int>& graph_flip_list ) {
	//
	vector<int> loop_haplotype = bdict.subset_haplotype;
	std::unordered_map<int,std::unordered_map<int,int>> correlation_graph;
	std::unordered_map<int,std::unordered_map<int,int>> high_correlation_graph;
	for (int l = 0; l < bdict.length_subset; l++) {
		for ( auto const &ent1 : bdict.subset_matrix.mat[l] ) {
			auto const &m = ent1.first;
			if ( bdict.subset_matrix(l,m) != 0 & l != m ) {
			int lm_energy = loop_haplotype[l]*loop_haplotype[m]*bdict.subset_matrix(l,m);	
			for ( auto const &ent2 : bdict.subset_matrix.mat[m] ) {
				auto const &p = ent2.first;
				if ( bdict.subset_matrix(l,p) != 0 & l != p & m != p ) {
				int lp_energy = loop_haplotype[l]*loop_haplotype[p]*bdict.subset_matrix(l,p);
				int correlation = lm_energy*lp_energy;
				//cout << correlation << endl;	
				if ( correlation > 0 ) {
					//cout << l << "\t" << m << "\t" << p << "\t" << correlation << endl;
					correlation_graph[m][p] += 1;
					correlation_graph[p][m] += 1; //.push_back(m);
					if (correlation_graph[m][p] > 30) {  //2-
						high_correlation_graph[m][p] += 1;	
						high_correlation_graph[p][m] += 1;	
					}
				} }
			} }
                }
	}

	vector<int> covered_list = {};
	vector<int> traverse_list;
	for (int l = 0; l < bdict.length_subset; l++) {
		if (std::find(covered_list.begin(), covered_list.end(), l) != covered_list.end()) {
			cout << "covered" << endl;
		} else {
        		traverse_graph( high_correlation_graph, traverse_list, l, bdict.length_subset );
			cout << "start " << l << endl;
			for (int l = 0; l < traverse_list.size(); l++) {
				//cout << traverse_list[l] << " ";			
				covered_list.push_back(traverse_list[l]);
			}
			cout << endl;
		}	
	}

	for (int l = 0; l < traverse_list.size(); l++) {
		//cout << traverse_list[l] << " "; 
		bdict.subset_haplotype[l] = loop_haplotype[l]*-1;	
	}
	cout << endl;
	graph_flip_list = traverse_list;
}

static void traverse_graph( std::unordered_map<int,std::unordered_map<int,int>> high_correlation_graph, vector<int>& traverse_list, int start_node, int subset_length ) {
	traverse_list.clear(); 
        bool *visited = new bool[subset_length];
        for(int i = 0; i < subset_length; i++) { visited[i] = false; }
        // Create a queue for BFS
        int s = start_node; 
	visited[s] = true;
        std::list<int> queue;
        queue.push_back(s);
	traverse_list.push_back(s);

        /////////////////////////////////////////////////////////
        while(!queue.empty()) {
                s = queue.front();  //std::cout << s << " ";
                queue.pop_front();
                for (auto& it : high_correlation_graph[s]) {
                        int j = it.first;
                        if( !visited[j] ) {
				//cout << high_correlation_graph[s][j] << endl;
                                visited[j] = true;
                                queue.push_back(j);
                                traverse_list.push_back(j);
                        }
                }
        }  //cout << endl;
        std::sort( traverse_list.begin(),traverse_list.end() );
};
////////////////////////////////////////////////////////////////////////////////////////////////////////
static void length_cutoff_nmatrix( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix, map_matrix<int>& nmatrix2 ) {
        
        for (int i = 0; i < pdict.num_paired; i++) {
                for (auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;      
                        int diff_pos = abs(pdict.sorted_paired_positions[m] - pdict.sorted_paired_positions[i]);
                        //if (diff_pos < pos_diff_cutoff) { nmatrix2.add_to(i,m,nmatrix(i,m)); }
                        cout << i << " " << m << " " << diff_pos << " " << nmatrix(i,m) << endl;
                        if (diff_pos < 30000000) { nmatrix2.add_to(i,m,nmatrix(i,m)); nmatrix2.add_to(m,i,nmatrix(i,m)); }
                        if (diff_pos > 30000000) { nmatrix.add_to(i,m,-nmatrix(i,m)); diff_matrix.add_to(m,i,-nmatrix(m,i));}
                        if (diff_pos > 30000000) { diff_matrix.add_to(i,m,-diff_matrix(i,m)); diff_matrix.add_to(m,i,-diff_matrix(m,i));}
                }
        }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////
void call_blocks(coord_dictionary& pdict,int switch_cutoff ) {
	int block_num = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
		if (pdict.switchE[i] > switch_cutoff) { block_num += 1; }
		pdict.block[i] = block_num;
	}
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void set_switch_energy( coord_dictionary& pdict, vector<bool>& flip_maxima ) {
        for (int i = 0; i < pdict.num_paired; i++) {
                if (!flip_maxima[i]) { pdict.switchE[i] = -10000.0; }
        }
};
////////////////////////////////////////////////////////////////////////////////////////////////////////
static void get_flip_positions( coord_dictionary& pdict, vector<bool>& flip_maxima ) {
        for (int i = 0; i < pdict.num_paired; i++) {
                vector<int> neighbors{i-1,i+1};  //cout << "  " << i << "  ";
                double energy = pdict.deltaE[i];  bool check = true;
		if ( energy < -200 ) {
                for (int j=0; j < neighbors.size(); j++) {   //cout << neighbors[j] << "  ";
                        if ( neighbors[j] >= 0 && neighbors[j] < pdict.num_paired ) {
				if ( pdict.deltaE[neighbors[j]] > energy ) { check = false; } 
			}
                } }
                flip_maxima[i] = check;
        }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void energy_sum_min_max( coord_dictionary& pdict ) {
        double sum = 0.0; double max = 0.0; double min = 0.0; int az_num = 0;
        double sum_blockE = 0.0; double max_blockE = 0.0; double min_blockE = 0.0; int az_num_blockE = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
                sum += pdict.deltaE[i];
                sum_blockE += pdict.switchE[i];
                if (pdict.deltaE[i] > max) { max = pdict.deltaE[i]; }
                if (pdict.deltaE[i] < min) { min = pdict.deltaE[i]; }
                if (pdict.deltaE[i] > 0) { az_num += 1; }
                if (pdict.switchE[i] > max_blockE) { max_blockE = pdict.switchE[i]; }
                if (pdict.switchE[i] < min_blockE) { min_blockE = pdict.switchE[i]; }
                if (pdict.switchE[i] > 0) { az_num_blockE += 1; }
        }
        cout << "-- spin flip energy sum: " << sum << " minimum: " << min << " maximum: " << max << " #_above_zero: " << az_num << endl;
        cout << "-- block flip energy sum: " << sum_blockE << " minimum: " << min_blockE << " maximum: " << max_blockE << " #_above_zero: " << az_num_blockE << endl;
};

// added negative signs
////////////////////////////////////////////////////////////////////////////////////////////////////////
static void single_spin_flip_recursive( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix ) {
        int nsf = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
                int init_i = pdict.haplotype[i]; //int final_i = -1*pdict.haplotype[i];
                double initial_energy = 0.0;     //double final_energy = 0.0;
/*
                for (auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;
                        int diff_pos = abs(pdict.sorted_paired_positions[m] - pdict.sorted_paired_positions[i]);
                        initial_energy += -1.0*init_i*pdict.haplotype[m]*diff_matrix(i,m);  //final_energy += final_i*pdict.haplotype[m]*diff_matrix(i,m);
                        //cout << i << " " << m << " " << pdict.haplotype[i] << " " << diff_pos << " " << nmatrix(i,m) << " " << diff_matrix(i,m) << endl;
                }
                double final_energy = -1*initial_energy;
                double diff_energy = final_energy - initial_energy;
                //cout << "init e " << initial_energy << " final e " << final_energy << " e diff " << diff_energy << endl;
                if (diff_energy < 0.0) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); nsf++; }
                pdict.deltaE[i] = (-1.0)*diff_energy;
*/
		if (i > 0 && i < pdict.num_paired) {
                    int spin_prod = pdict.haplotype[i-1]*pdict.haplotype[i+1];
                    int pos_prod = pdict.haplotype[i]*pdict.haplotype[i+1];
                    if (spin_prod == 1 && pos_prod == -1){ pdict.haplotype[i] = pdict.haplotype[i]*(-1); nsf++; }
                }
        }
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////
static void single_spin_flip_recursive_pop( coord_dictionary& pdict, map_matrix<int>& nmatrix ) {
        int nsf = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
                int init_i = pdict.haplotype[i]; //int final_i = -1*pdict.haplotype[i];
                double initial_energy = 0.0;     //double final_energy = 0.0;
                for (auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;			
			initial_energy += -1.0*init_i*pdict.haplotype[m]*nmatrix(i,m);
                }
                double final_energy = -1*initial_energy;
                double diff_energy = final_energy - initial_energy;
                if (diff_energy < 0.0) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); nsf++; }
                pdict.deltaE[i] = (-1.0)*diff_energy;
        }
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////
static void switchE( coord_dictionary& pdict, map_matrix<int>& nmatrix ) {
        int nbf = 0;
        double initial_energy = 0.0;
        int site_diff;
        int min_pos = 0;
        double min_switch = 0;

        for (int i = 0; i < pdict.num_paired; i++) {

            for (auto const &ent1 : nmatrix.mat[i]) {
                auto const &m = ent1.first;
                initial_energy += -1.0*pdict.haplotype[i]*pdict.haplotype[m]*nmatrix(i,m); 
            }
        }

        for ( int i = 0; i < pdict.num_paired; i++ ) {
		
            int initial_pos  = pdict.sorted_paired_positions[i];
            double final_energy = 0.0;
            vector<int> switch_haplotype = pdict.haplotype;
            int j = i+1;
            bool window = true;
            while(window == true && j < pdict.num_paired) { 
                switch_haplotype[j] = pdict.haplotype[j]*(-1);
                j+=1;
                if(pdict.sorted_paired_positions[j] > initial_pos+1e5) {
                    window = false;
                }
            }
            for (int k = 0; k < pdict.num_paired; k++) {
                for (auto const &ent1 : nmatrix.mat[k]){
                    auto const &m = ent1.first;
                    final_energy += -1.0*switch_haplotype[k]*switch_haplotype[m]*nmatrix(k,m);
                }
            }
            double diffE = final_energy - initial_energy;
            pdict.switchE[i] = (-1.0)*diffE;
            //if ( diffE < min_switch ) { min_switch = diffE; min_pos = i; }
	if ( diffE < 0 ) {pdict.haplotype = switch_haplotype; initial_energy = final_energy; nbf++;}
        }
	cout << "nbf:" << "\t" << nbf << endl;

        //for (int i = min_pos+1; i < pdict.num_paired; i++) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); }
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
static void switchE_recursive( coord_dictionary& pdict, map_matrix<int>& nmatrix ) {
        int nbf = 0;
        int site_diff;
        int min_pos = 0;
        vector<int> h_plus;
        vector<int> h_minus;
        double initial_energy = 0.0;
        for (int i = 0; i < pdict.num_paired; i++) {
            int h_left = 0;
            int h_right = 0;
            for (auto const &ent1 : nmatrix.mat[i]) { 
                auto const &m = ent1.first;
                if (m < i) { h_left = h_left + nmatrix(i,m)*pdict.haplotype[m]; }
                else { h_right = h_right + nmatrix(i,m)*pdict.haplotype[m]; }
            }
            h_minus.push_back(h_left);
            h_plus.push_back(h_right);
        }
        
        for (int i = 0; i < pdict.num_paired; i++) {

            for (auto const &ent1 : nmatrix.mat[i]) {
                auto const &m = ent1.first;
                initial_energy += -1.0*pdict.haplotype[i]*pdict.haplotype[m]*nmatrix(i,m); 
            }
        }


        
        double final_energy = 0.0;
        vector<int> switch_haplotype = pdict.haplotype;
        for (int j = 1; j < pdict.num_paired; j++) { switch_haplotype[j] = pdict.haplotype[j]*(-1); }
        for (int i = 0; i < pdict.num_paired; i++) {
            for (auto const &ent1 : nmatrix.mat[i]){
                auto const &m = ent1.first;
                final_energy += -1.0*switch_haplotype[i]*switch_haplotype[m]*nmatrix(i,m);
            }
        }
        double diffE_last = final_energy - initial_energy;
        double min_switch = diffE_last;
        pdict.switchE[0] = (-1.0)*diffE_last;
        for (int j = 1; j < pdict.num_paired; j++) {
            diffE_last = diffE_last + pdict.haplotype[j]*(h_plus[j]-h_minus[j]);
            if (diffE_last < min_switch) { min_switch = diffE_last; min_pos = j; }
            pdict.switchE[j] = (-1.0)*diffE_last;
        }
        for (int i = min_pos+1; i < pdict.num_paired; i++) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); }
	cout << min_pos << endl;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////
static void calculate_switchE_pop( coord_dictionary& pdict, map_matrix<int>& nmatrix, int pop_weight, map_matrix<double> diff_matrix ) {//, map_matrix<double> diff_matrix, int pop_weight ) {
        int nbf = 0;
        int site_diff;
        int min_pos = 0;
        vector<int> h_plus;
        vector<int> h_minus;
        double initial_energy = 0.0;
        double factor;
        for (int i = 0; i < pdict.num_paired; i++) {
            int h_left = 0;
            int h_right = 0;
            for (auto const &ent1 : nmatrix.mat[i]) { 
                auto const &m = ent1.first;
                //if (pdict.sorted_paired_positions[m] < pdict.sorted_paired_positions[i]) { h_left = h_left + nmatrix(i,m)*pdict.haplotype[m]; }
                if ( pdict.sorted_paired_positions[m] < pdict.sorted_paired_positions[i] ) { h_left = h_left + diff_matrix(i,m)*pdict.haplotype[m]; }
                //if (pdict.sorted_paired_positions[m] > pdict.sorted_paired_positions[i]) { h_right = h_right + nmatrix(i,m)*pdict.haplotype[m]; }
                if ( pdict.sorted_paired_positions[m] > pdict.sorted_paired_positions[i] ) { h_right = h_right + diff_matrix(i,m)*pdict.haplotype[m]; }
                
            }
            h_minus.push_back(h_left);
            h_plus.push_back(h_right);
        }
        
        for (int i = 0; i < pdict.num_paired; i++) {

            for (auto const &ent1 : nmatrix.mat[i]) {
                auto const &m = ent1.first;
                //initial_energy += -1.0*pdict.haplotype[i]*pdict.haplotype[m]*nmatrix(i,m); 
                initial_energy += -1.0*pdict.haplotype[i]*pdict.haplotype[m]*diff_matrix(i,m);
            }
            if (i < pdict.num_paired-1) {
                site_diff = abs(pdict.sorted_paired_positions[i]-pdict.sorted_paired_positions[i+1]);
                factor = (double)site_diff/pop_weight;
                initial_energy += -1.0*exp(-1.0*factor)*pdict.haplotype[i]*pdict.haplotype[i+1]*pdict.pop_hap[i]*pdict.pop_hap[i+1];
            }
        }


        
        double final_energy = 0.0;
        vector<int> switch_haplotype = pdict.haplotype;
        for (int j = 1; j < pdict.num_paired; j++) { switch_haplotype[j] = pdict.haplotype[j]*(-1); }
        for (int i = 0; i < pdict.num_paired; i++) {
            for (auto const &ent1 : nmatrix.mat[i]){
                auto const &m = ent1.first;
                //final_energy += -1.0*switch_haplotype[i]*switch_haplotype[m]*nmatrix(i,m);
                final_energy += -1.0*switch_haplotype[i]*switch_haplotype[m]*diff_matrix(i,m);
            }
            if (i < pdict.num_paired-1) {
                site_diff = abs(pdict.sorted_paired_positions[i]-pdict.sorted_paired_positions[i+1]);
                factor = (double)site_diff/pop_weight;
                final_energy += -1.0*exp(-1.0*factor)*switch_haplotype[i]*switch_haplotype[i+1]*pdict.pop_hap[i]*pdict.pop_hap[i+1];
            }
        }
        double diffE_last = final_energy - initial_energy;
        double min_switch = diffE_last;
        pdict.switchE[0] = (-1.0)*diffE_last;
        for (int j = 1; j < pdict.num_paired-1; j++) {
            int site_diff1 = abs(pdict.sorted_paired_positions[j] - pdict.sorted_paired_positions[j+1]);
            int site_diff2 = abs(pdict.sorted_paired_positions[j-1] - pdict.sorted_paired_positions[j]);
            double factor1 = (double)site_diff1/pop_weight;
            double factor2 = (double)site_diff2/pop_weight;
            diffE_last = diffE_last + pdict.haplotype[j]*(h_plus[j]-h_minus[j]) + pdict.haplotype[j]*pdict.pop_hap[j]*(exp(-1.0*factor1)*pdict.haplotype[j+1]*pdict.pop_hap[j+1] - exp(-1.0*factor2)*pdict.haplotype[j-1]*pdict.pop_hap[j-1]);
            //diffE_last = diffE_last + pdict.haplotype[j]*(h_plus[j]-h_minus[j]);
            //if (j < pdict.num_paired-1) { diffE_last = diffE_last - exp(-1.0*factor)*pdict.haplotype[j+1]*pdict.haplotype[j]*pdict.pop_hap[j+1]*pdict.pop_hap[j]; }
            if ( diffE_last < min_switch ) { min_switch = diffE_last; min_pos = j; }
            pdict.switchE[j] = (-1.0)*diffE_last;
        }
        //for (int i = min_pos+1; i < pdict.num_paired; i++) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); }
    //scout << min_pos << endl;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////
static void calculate_switchE( coord_dictionary& pdict, map_matrix<int>& nmatrix ) {
        int nbf = 0;
        int site_diff;
        int min_pos = 0;
        vector<int> h_plus;
        vector<int> h_minus;
        double initial_energy = 0.0;
        
        for (int i = 0; i < pdict.num_paired; i++) {
            int h_left = 0;
            int h_right = 0;
            for (auto const &ent1 : nmatrix.mat[i]) { 
                auto const &m = ent1.first;
                if (m < i) { h_left = h_left + nmatrix(i,m)*pdict.haplotype[m]; }
                else { h_right = h_right + nmatrix(i,m)*pdict.haplotype[m]; }
            }
            h_minus.push_back(h_left);
            h_plus.push_back(h_right);
        }
        
        for (int i = 0; i < pdict.num_paired; i++) {

            for (auto const &ent1 : nmatrix.mat[i]) {
                auto const &m = ent1.first;
                initial_energy += -1.0*pdict.haplotype[i]*pdict.haplotype[m]*nmatrix(i,m); 
            }
        }


        
        double final_energy = 0.0;
        vector<int> switch_haplotype = pdict.haplotype;
        for (int j = 1; j < pdict.num_paired; j++) { switch_haplotype[j] = pdict.haplotype[j]*(-1); }
        for (int i = 0; i < pdict.num_paired; i++) {
            for (auto const &ent1 : nmatrix.mat[i]){
                auto const &m = ent1.first;
                final_energy += -1.0*switch_haplotype[i]*switch_haplotype[m]*nmatrix(i,m);
            }
        }
        double diffE_last = final_energy - initial_energy;
        double min_switch = diffE_last;
        pdict.switchE[0] = (-1.0)*diffE_last;
        
        for (int j = 1; j < pdict.num_paired; j++) {
            diffE_last = diffE_last + pdict.haplotype[j]*(h_plus[j]-h_minus[j]);
            if (diffE_last < min_switch) { min_switch = diffE_last; min_pos = j; }
            pdict.switchE[j] = (-1.0)*diffE_last;
        }
        //for (int i = min_pos+1; i < pdict.num_paired; i++) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); }
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
static void switchE_recursive_window( coord_dictionary& pdict, map_matrix<int>& nmatrix ) {
    int nbf = 0;
    calculate_switchE(pdict, nmatrix);
    for ( int i = 0; i < pdict.num_paired; i++ ) {
            int initial_pos  = pdict.sorted_paired_positions[i];
            vector<int> switch_haplotype = pdict.haplotype;
            int j = i+1;
            bool window = true;
            while(window == true && j < pdict.num_paired) { 
                switch_haplotype[j] = pdict.haplotype[j]*(-1);
                j+=1;
                if(pdict.sorted_paired_positions[j] > initial_pos+1e5) {
                    window = false;
                }
            }

            double windowE = (-1.0)*(pdict.switchE[i] - pdict.switchE[j]);
            if( windowE < 0.0 ) {pdict.haplotype = switch_haplotype; calculate_switchE(pdict, nmatrix); nbf++; cout << nbf << endl;}
    }
    cout << "window switches: " << nbf << endl;
};       
////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_flip_recursive( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix ) {
        int nbf = 0;
	int site_diff;
        double initial_energy = 0.0;
        for (int i = 0; i < pdict.num_paired; i++) {
                for (auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;
			site_diff = abs(pdict.sorted_paired_positions[i]-pdict.sorted_paired_positions[m]);
			
            		if (site_diff < 1e6) { initial_energy += -1.0*pdict.haplotype[i]*pdict.haplotype[m]*nmatrix(i,m); }
                        //initial_energy += -1.0*pdict.haplotype[i]*pdict.haplotype[m]*diff_matrix(i,m);
                }
        }
        double energyA = initial_energy;
        double energyB = initial_energy;
        vector<int> loop_haplotype = pdict.haplotype;
        for (int i = 0; i < pdict.num_paired; i++) {
                double save_energyB = energyB;
                vector<int> switch_haplotype = loop_haplotype;
                for (int j = 0; j < i+1; j++) { switch_haplotype[j] = loop_haplotype[j]*(-1); }
                int hap_spin  = loop_haplotype[i];
                int flip_spin = -1*loop_haplotype[i];
                double spin_i_energy = 0.0;
                double spin_f_energy = 0.0;
                for (auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;
			site_diff = abs(pdict.sorted_paired_positions[i]-pdict.sorted_paired_positions[m]);
			
			if (site_diff < 1e6) { 
				cout << "position: " << pdict.sorted_paired_positions[i] << endl;
				cout << "site diff: " << site_diff << endl;
                        	spin_i_energy += -1.0*hap_spin*switch_haplotype[m]*diff_matrix(i,m);  /// may need to add negative sign
                        	spin_f_energy += -1.0*flip_spin*switch_haplotype[m]*diff_matrix(i,m);
			}

                }
                double switch_energy = spin_f_energy - spin_i_energy;
                energyA = energyA + switch_energy;
                double diff_energy = energyA - energyB;
				
				int v1 = rand() % 100 + 1;
				int flipped = 0;
				cout << "Rand is: " << "\t" << v1 << endl;
                if ( diff_energy < 0.0 && v1 > 30) {
						cout << "Take the switch" << endl;
                        loop_haplotype = switch_haplotype;
                        energyB = energyA;
                        energyA = save_energyB;
                        nbf++;
						flipped++;
                }
				
				if ( diff_energy > 0.0 && v1 < 10 ) {
						loop_haplotype = switch_haplotype;
                        energyB = energyA;
                        energyA = save_energyB;
                        nbf++; 
						flipped++;
				}
				
                if ( flipped == 0 ) {	
                        energyB = save_energyB;
                        energyA = energyA;
                }
                pdict.switchE[i] = (-1.0)*diff_energy;
        }
        pdict.haplotype = loop_haplotype;
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// static void block_flip_recursive_var_range( coord_dictionary& pdict, map_matrix<int>& nmatrix, int range, int& prior, int pop_weight, map_matrix<double> diff_matrix ) {
//         int nbf = 0;
//         double initial_energy = 0.0;
//         int site_diff;
//         double factor;
//         int dist = 0;
//         double min_energy = 0.0;
//         vector<int> min_haplotype = pdict.haplotype;
        
//         for (int i = 0; i < pdict.num_paired; i++) {

//             for (auto const &ent1 : nmatrix.mat[i]) {
//             	auto const &m = ent1.first;
//                 initial_energy += -1.0*pdict.haplotype[i]*pdict.haplotype[m]*diff_matrix(i,m);     
//             }
//             if (i < pdict.num_paired-1) {
//             	site_diff = abs(pdict.sorted_paired_positions[i]-pdict.sorted_paired_positions[i+1]);
//             	factor = (double)site_diff/pop_weight;
//                 //COMMENTED IN TESTING 12/01/2020
//             	initial_energy += -1.0*exp(-1.0*factor)*pdict.haplotype[i]*pdict.haplotype[i+1]*pdict.pop_hap[i]*pdict.pop_hap[i+1];
//             }

//         }

//         double energyA = initial_energy;
//         double energyB = initial_energy;
//         vector<int> loop_haplotype = pdict.haplotype;

//         for (int i = 0; i < pdict.num_paired; i++) {
//             vector<int> switched;
//             double save_energyB = energyB;
//             vector<int> switch_haplotype = loop_haplotype;
//             int var_count =0;
//             int j = i+1;
//             while( var_count < range && j < pdict.num_paired )  {
//                 switch_haplotype[j] = pdict.haplotype[j]*(-1);
//                 switched.push_back(j);
//                 if ( pdict.pop_hap[j] != 0 ) { var_count+=1; }
//                 j+=1;
//             }

            
//             int hap_spin  = loop_haplotype[i];    
//             int flip_spin = -1*loop_haplotype[i];
//             double spin_i_energy = 0.0;
//             double spin_f_energy = 0.0; 
//             for (int k = 0; k < switched.size(); k++) {
//                 for (auto const &ent1 : nmatrix.mat[switched[k]]) {
//                     auto const &m = ent1.first; 
//                     spin_f_energy += -1.0*switch_haplotype[switched[k]]*switch_haplotype[m]*diff_matrix(switched[k],m);                      
//                     spin_i_energy += -1.0*pdict.haplotype[switched[k]]*pdict.haplotype[m]*diff_matrix(switched[k],m); 
// 				}
// 				if (switched[k] < pdict.num_paired-1) {
// 					site_diff = abs(pdict.sorted_paired_positions[switched[k]]-pdict.sorted_paired_positions[switched[k]+1]);
//                     factor = (double)site_diff/pop_weight;
//                     //FOLLOWING 2 LINES COMMENTED IN TESTING 12/01/2020
// 					spin_f_energy += -1.0*exp(-1.0*factor)*switch_haplotype[switched[k]]*switch_haplotype[switched[k]+1]*pdict.pop_hap[switched[k]]*pdict.pop_hap[switched[k]+1];
// 					spin_i_energy += -1.0*exp(-1.0*factor)*pdict.haplotype[switched[k]]*pdict.haplotype[switched[k]+1]*pdict.pop_hap[switched[k]]*pdict.pop_hap[switched[k]+1];
// 				}
//             }
            
//             double switch_energy = spin_f_energy - spin_i_energy;

//             //Block commented on 12/21
//             /*
//             if ( switch_energy < 0.0) {
//                 loop_haplotype = switch_haplotype;
//                 //energyB = energyA;
//                 //energyA = save_energyB;
//                 nbf++;
//                 pdict.switchE[i] = (-1.0)*switch_energy;
//                 i = j;
//             }
//             */
//             //pdict.switchE[i] = (-1.0)*switch_energy;
//             if ( switch_energy < min_energy ) {min_energy = switch_energy; min_haplotype = switch_haplotype; }
//         }
//         //cout << "number block flips: " << "\t" << nbf << endl;
//         cout << "min switching energy" << "\t" << min_energy << endl;
//         //if ( nbf < prior ) { pdict.haplotype = loop_haplotype; prior = nbf; } //ComMENTED on 12/21
//         if ( min_energy < 0.0 ) { pdict.haplotype = min_haplotype; prior = prior-1; }
// };
/////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_flip_recursive_pop( coord_dictionary& pdict, map_matrix<int>& nmatrix, int dist, int pop_weight, int& prior, map_matrix<double> diff_matrix ) {
        int nbf = 0;
        double initial_energy = 0.0;
        double min_energy = 0.0;
        int site_diff;
        double factor;
        int last_pos;
        vector<int> min_haplotype = pdict.haplotype;

        
        for (int i = 0; i < pdict.num_paired; i++) {

            for (auto const &ent1 : nmatrix.mat[i]) {
            	auto const &m = ent1.first;
                initial_energy += -1.0*pdict.haplotype[i]*pdict.haplotype[m]*diff_matrix(i,m);     
            }
            if (i < pdict.num_paired-1) {
            	site_diff = abs(pdict.sorted_paired_positions[i]-pdict.sorted_paired_positions[i+1]);
            	factor = (double)site_diff/pop_weight;
                //FOLLOWING LINE COMMENTED 12/01/2020
            	initial_energy += -1.0*exp(-1.0*factor)*pdict.haplotype[i]*pdict.haplotype[i+1]*pdict.pop_hap[i]*pdict.pop_hap[i+1];
            }

        }

        double energyA = initial_energy;
        double energyB = initial_energy;
        vector<int> loop_haplotype = pdict.haplotype;

        for (int i = 0; i < pdict.num_paired; i++) {
            vector<int> switched;
            double save_energyB = energyB;
            vector<int> switch_haplotype = loop_haplotype;
            int initial_pos = pdict.sorted_paired_positions[i];
            int j = i+1;
            bool window = true;
            while(window == true && j < pdict.num_paired) { 
                switch_haplotype[j] = pdict.haplotype[j]*(-1);
                switched.push_back(j);
                j+=1;
                if(pdict.sorted_paired_positions[j] > initial_pos+dist) {
                    window = false;
                }
            }

            
            int hap_spin  = loop_haplotype[i];    
            int flip_spin = -1*loop_haplotype[i];
            double spin_i_energy = 0.0;
            double spin_f_energy = 0.0;

            for (int k = 0; k < switched.size(); k++) {
                for (auto const &ent1 : nmatrix.mat[switched[k]]) {
                    auto const &m = ent1.first; 
                    spin_f_energy += -1.0*switch_haplotype[switched[k]]*switch_haplotype[m]*diff_matrix(switched[k],m); 
                    spin_i_energy += -1.0*pdict.haplotype[switched[k]]*pdict.haplotype[m]*diff_matrix(switched[k],m); 
				}

				if (switched[k] < pdict.num_paired-1) {
					site_diff = abs(pdict.sorted_paired_positions[switched[k]]-pdict.sorted_paired_positions[switched[k]+1]);
					factor = (double)site_diff/pop_weight;
                    //FOLLOWING 2 LINES COMMENTED IN TESTING 12/01/2020
					spin_f_energy += -1.0*exp(-1.0*factor)*switch_haplotype[switched[k]]*switch_haplotype[switched[k]+1]*pdict.pop_hap[switched[k]]*pdict.pop_hap[switched[k]+1];
					spin_i_energy += -1.0*exp(-1.0*factor)*pdict.haplotype[switched[k]]*pdict.haplotype[switched[k]+1]*pdict.pop_hap[switched[k]]*pdict.pop_hap[switched[k]+1];
				}
            }
            
            double switch_energy = spin_f_energy - spin_i_energy;

            pdict.switchE[i] = (-1.0)*switch_energy; //THIS WAS ADDED 12/21/2020

            if ( switch_energy < min_energy) {min_energy = switch_energy; min_haplotype = switch_haplotype; } //THIS WAS ADDED 12/21/2020
//FOLLOWING BLOCK COMMENTED ON 12/21/2020
           /*
            
            if ( switch_energy < 0.0) {
                loop_haplotype = switch_haplotype;
                //energyB = energyA;
                //energyA = save_energyB;
                nbf++;
                pdict.switchE[i] = (-1.0)*switch_energy;
                i = j;
            }
            */

            /*
            else {
                energyB = save_energyB;
                energyA = energyA;
            }
            */
            //pdict.switchE[i] = (-1.0)*switch_energy;
        }

        cout << "min energy: " << "\t" << min_energy << endl;
        //Following line commented on DEC 21
        //if ( nbf < prior ) { pdict.haplotype = loop_haplotype; prior = nbf; }
        if (min_energy < 0) { pdict.haplotype = min_haplotype; prior = prior-1; }
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_flip_fast( coord_dictionary& pdict, map_matrix<int>& nmatrix, int dist, int pop_weight, int& prior, map_matrix<double> diff_matrix ) {
        int nbf = 0;
        double initial_energy = 0.0;
        int site_diff;
        double factor;
        int last_pos;
        
        calculate_switchE_pop( pdict, nmatrix, pop_weight, diff_matrix );
        vector<int> loop_haplotype = pdict.haplotype;

        for (int i = 0; i < pdict.num_paired; i++) {
        	
            vector<int> switched;
            vector<int> switch_haplotype = pdict.haplotype;
            int initial_pos = pdict.sorted_paired_positions[i];
            int j = i+1;
            bool window = true;
            while(window == true && j < pdict.num_paired) { 
                switch_haplotype[j] = pdict.haplotype[j]*(-1);
                switched.push_back(j);
                j+=1;
                if(pdict.sorted_paired_positions[j] > initial_pos+dist) {
                    window = false;
                }
            }
            
            double switch_energy = (-1.0)*pdict.switchE[i] - (-1.0)*pdict.switchE[j];
            //double switch_energy = spin_f_energy - spin_i_energy;
            //double switch_energy = spin_i_energy - spin_f_energy;
            //energyA += switch_energy;
            //double diff_energy = energyA - energyB;
            
            if ( switch_energy < 0.0) {
                pdict.haplotype = switch_haplotype;
                calculate_switchE_pop( pdict, nmatrix, pop_weight, diff_matrix );
                //energyB = energyA;
                //energyA = save_energyB;
                nbf++;
                //pdict.switchE[i] = (-1.0)*switch_energy;
                cout << i << endl;
                i = j;
            }
        }
        cout << "number block flips: " << "\t" << nbf << endl;
        if ( nbf < prior ) { pdict.haplotype = loop_haplotype; prior = nbf; }
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////
static void energy_sum_min_max_hic( block_dictionary& bdict ) {
        double sum = 0.0; double max = 0.0; double min = 0.0;
        for (int i = 0; i < bdict.length_subset; i++) {
                sum += bdict.subset_deltaE[i];
                if (bdict.subset_deltaE[i] > max) { max = bdict.subset_deltaE[i]; }
                if (bdict.subset_deltaE[i] < min) { min = bdict.subset_deltaE[i]; }
        }
        cout << "-- energy sum: " << sum << " minimum: " << min << " maximum: " << max << endl;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////
static void single_spin_flip_recursive_hic( block_dictionary& bdict ) {
        int nsf = 0;
        for (int i = 0; i < bdict.length_subset; i++) {
                int init_i = bdict.subset_haplotype[i]; //int final_i = -1*pdict.haplotype[i];
                double initial_energy = 0.0;     //double final_energy = 0.0;
                for (auto const &ent1 : bdict.distance_matrix.mat[i]) {
                        auto const &m = ent1.first;
                        //int diff_pos = abs(pdict.sorted_paired_positions[m] - pdict.sorted_paired_positions[i]);
                        initial_energy += -1.0*init_i*bdict.subset_haplotype[m]*bdict.distance_matrix(i,m);  //final_energy += final_i*pdict.haplotype[m]*diff_matrix(i,m);
                }
                double final_energy = -1*initial_energy;
                double diff_energy = final_energy - initial_energy;
                //cout << "init e " << initial_energy << " final e " << final_energy << " e diff " << diff_energy << endl;
                if (diff_energy < 0.0) { bdict.subset_haplotype[i] = bdict.subset_haplotype[i]*(-1); nsf++; }
                bdict.subset_deltaE[i] = (-1.0)*diff_energy;
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_flip_recursive_hic( block_dictionary& bdict ) {
        int nbf = 0;
        double initial_energy = 0.0;
        for (int i = 0; i < bdict.length_subset; i++) {
                for (auto const &ent1 : bdict.distance_matrix.mat[i]) {
                        auto const &m = ent1.first;
                        initial_energy += -1.0*bdict.subset_haplotype[i]*bdict.subset_haplotype[m]*bdict.distance_matrix(i,m);
                }
        }
        double energyA = initial_energy;
        double energyB = initial_energy;
        //vector<int> loop_haplotype = pdict.haplotype;
        vector<int> loop_haplotype = bdict.subset_haplotype;
        for (int i = 0; i < bdict.length_subset; i++) {
                double save_energyB = energyB;
                vector<int> switch_haplotype = loop_haplotype;
                for (int j = 0; j < i+1; j++) { switch_haplotype[j] = loop_haplotype[j]*(-1); }
                int hap_spin  = loop_haplotype[i];
                int flip_spin = -1*loop_haplotype[i];
                double spin_i_energy = 0.0;
                double spin_f_energy = 0.0;
                for (auto const &ent1 : bdict.distance_matrix.mat[i]) {
                        auto const &m = ent1.first;
                        spin_i_energy += -1.0*hap_spin*switch_haplotype[m]*bdict.distance_matrix(i,m);  /// may need to add negative sign
                        spin_f_energy += -1.0*flip_spin*switch_haplotype[m]*bdict.distance_matrix(i,m);
                }
                double switch_energy = spin_f_energy - spin_i_energy;
                energyA = energyA + switch_energy;
                double diff_energy = energyA - energyB;
                if ( diff_energy < 0.0 ) {
                        loop_haplotype = switch_haplotype;
                        energyB = energyA;
                        energyA = save_energyB;
                        nbf++;
                }
                else {
                        energyB = save_energyB;
                        energyA = energyA;
                }
                bdict.subset_switchE[i] = (-1.0)*diff_energy;
        }
        bdict.subset_haplotype = loop_haplotype;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////
static void single_spin_flip_map( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix ) {
        int nsf = 0; clock_t t;
        for (int i = 0; i < pdict.num_paired; i++) {
                vector<int> sparse_map,sparse_hap,init_sub_haplotype,final_sub_haplotype;
                int l = 0; bool cross_i = true; int spin_pos;  double init_energy,final_energy,diff_energy;
                for(auto const &ent1 : nmatrix.mat[i]) {
                        auto const &m = ent1.first;
                        if( m>i && cross_i) { sparse_map.push_back(i); cross_i = false; spin_pos = l; l++; }
                        int diff_pos = abs(pdict.sorted_paired_positions[m] - pdict.sorted_paired_positions[i]);   ////// CHANGED THIS
                        //int diff_pos2 = abs(pdict.double_positions[m] - pdict.double_positions[i]);
                        //cout << "i " << i << "\tm " << m << "\t" << pdict.sorted_paired_positions[i] << "\t" << pdict.sorted_paired_positions[m] << "\t  sorted_paired " << diff_pos << "\t" << pdict.double_positions[i] << "\t" << pdict.double_positions[m] << "\t double_pos " << diff_pos2 << endl;      //////////////////////////
                        if (diff_pos < pos_diff_cutoff) { sparse_map.push_back(m); l++;}
                }
                if(cross_i) { sparse_map.push_back(i); cross_i = false; spin_pos = l; }
                for (int j = 0; j < sparse_map.size(); j++) { sparse_hap.push_back(pdict.haplotype[sparse_map[j]]); }
                submatrix_opt tempsub;  tempsub.initialize(sparse_map,diff_matrix);
                init_sub_haplotype = sparse_hap;
                final_sub_haplotype = init_sub_haplotype;
                //cout << "flip  position: " << i << " sparse_map_size " << sparse_map.size() << " spin_pos " << spin_pos << "  total " << l << endl;
                final_sub_haplotype[spin_pos] = final_sub_haplotype[spin_pos]*(-1);
                init_energy  = energy_function_opt(tempsub,init_sub_haplotype);
                final_energy = energy_function_opt(tempsub,final_sub_haplotype);
                diff_energy  = final_energy - init_energy;        //cout << "init e " << init_energy << " final e " << final_energy << " e diff " << diff_energy << endl;
                if (diff_energy < 0.0) { pdict.haplotype[i] = pdict.haplotype[i]*(-1); nsf++; }
                pdict.deltaE[i] = (-1.0)*diff_energy;
        }
        cout << "-- number of spin flips: " << nsf << endl;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_flip_map( coord_dictionary& pdict, map_matrix<double>& diff_matrix, map_matrix<int>& nmatrix, vector<bool>& flip_maxima ) {
        int nbf = 0;
        for (int i = 0; i < pdict.num_paired; i++) {
                if (flip_maxima[i]) {                   //if (pdict.switchE[i] > -switch_cut) {
                vector<int> sparse_map,sparse_hap,init_sub_haplotype,final_sub_haplotype;
                int flip_pos,l=0;   double init_energy,final_energy,diff_energy;
                for ( int j = pdict.flip_low_bound[i]; j <= pdict.flip_up_bound[i]; j++ ) {
                        if (j == i) { flip_pos = l; }
                        int diff_pos = abs(pdict.sorted_paired_positions[j] - pdict.sorted_paired_positions[i]);   ////// CHANGED THIS
                        //int diff_pos2 = abs(pdict.double_positions[j] - pdict.double_positions[i]);
                        //cout << "i " << i << "\tm " << j << "\t" << pdict.sorted_paired_positions[i] << "\t" << pdict.sorted_paired_positions[j] << "\t  sorted_paired " << diff_pos << "\t" << pdict.double_positions[i] << "\t" << pdict.double_positions[j] << "\t double_pos " << diff_pos2 << endl;
                        if (diff_pos < pos_diff_cutoff) { sparse_map.push_back(j); l++; }
                }
                for (int j = 0; j < sparse_map.size(); j++) { sparse_hap.push_back(pdict.haplotype[sparse_map[j]]); }
                //cout << "block position: " << i << " sparse_map_size " << sparse_map.size() << " flip_pos " << flip_pos << "  total " << l << endl;
                submatrix_opt tempsub; tempsub.initialize(sparse_map,diff_matrix);
                init_sub_haplotype = sparse_hap;
                final_sub_haplotype = init_sub_haplotype;
                for (int j = 0; j < flip_pos; j++) { final_sub_haplotype[j] = final_sub_haplotype[j]*(-1); }
                init_energy  = energy_function_opt(tempsub,init_sub_haplotype);
                final_energy = energy_function_opt(tempsub,final_sub_haplotype);
                diff_energy  = final_energy - init_energy;        //cout << "init e " << init_energy << " final e " << final_energy << " e diff " << diff_energy << endl;
                if (diff_energy < 0.0) { for (int l = 0; l < i; l++) { pdict.haplotype[l] = pdict.haplotype[l]*(-1); } nbf++; }
                pdict.switchE[i] = (-1.0)*diff_energy;
                }    //}
        }
        cout << "-- number of block flips: " << nbf << endl;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static void block_flip_brute_force_hic( block_dictionary& bdict ) {
        int nbf = 0;
        vector<int> loop_haplotype = bdict.subset_haplotype;
        for (int i = 0; i < bdict.length_subset; i++) {
                vector<int> switch_haplotype = loop_haplotype;
                for (int j = 0; j < i+1; j++) { switch_haplotype[j] = loop_haplotype[j]*(-1); }
                double spin_i_energy = 0.0;  
		double spin_f_energy = 0.0;
		for (int l = 0; l < bdict.length_subset; l++) {
                	for (auto const &ent1 : bdict.distance_matrix.mat[l]) {
                        	auto const &m = ent1.first;
                        	spin_i_energy += -1.0*loop_haplotype[l]*loop_haplotype[m]*bdict.distance_matrix(l,m);  /// may need to add negative sign
                        	spin_f_energy += -1.0*switch_haplotype[l]*switch_haplotype[m]*bdict.distance_matrix(l,m);
				//cout << l << "\t" << m << "\t"<<  bdict.distance_matrix(l,m) << "\t" << bdict.distance_matrix(m,l) << endl;
                	}
		}
                double switch_energy = spin_f_energy - spin_i_energy;
		//cout << "switch energy:\t" << i << "\t" << bdict.block_map_inverted[i] << "\t" << bdict.num_hets[bdict.block_map_inverted[i]] << "\t" << spin_i_energy << "\t" << spin_f_energy << "\t" << switch_energy << endl;
		bool flipped = false;
                if ( switch_energy < 0.0 ) {
                        loop_haplotype = switch_haplotype;
			flipped = true;
                        nbf++;
                }
		if ( flipped == true ) { bdict.subset_switchE[i] = switch_energy; }
		else { bdict.subset_switchE[i] = (-1.0)*switch_energy; } 
        }
        bdict.subset_haplotype = loop_haplotype;
        cout << "-- number of block flips: " << nbf << endl;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static double energy_function_opt( submatrix_opt tempsub, vector<int> haplotype ) {
        vector<double> temp_vector = dot_product_matrix(haplotype,tempsub.smat);
        double energy_return = (-0.5)*dot_product(haplotype,temp_vector);
        return energy_return;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static double dot_product( vector<int> a, vector<double> b ) {
        double dot = 0.0; for (int i = 0; i < a.size(); i++) { dot += (double)(a[i])*(b[i]); };
        return dot;
};

////////////////////////////////////////////////////////////////////////////////////////////////////////
static vector<double> dot_product_matrix( vector<int> a, std::vector< std::vector<double> > b ) {
        vector<double> dot_vector;
        for (int i = 0; i < a.size(); i++) {
                double i_value = 0.0; for (int j = 0; j < b.size(); j++) { i_value += (double)a[j]*b[i][j]; }
                dot_vector.push_back(i_value);
        }
        return dot_vector;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////
//static void energy_sum_min_max( coord_dictionary& pdict ) {
//        double sum = 0.0; double max = 0.0; double min = 0.0;
//        for (int i = 0; i < pdict.num_paired; i++) {
//                sum += pdict.deltaE[i];
//                if (pdict.deltaE[i] > max) { max = pdict.deltaE[i]; }
//                if (pdict.deltaE[i] < min) { min = pdict.deltaE[i]; }
//        }
//        cout << "-- energy sum: " << sum << " minimum: " << min << " maximum: " << max << endl;
//};

