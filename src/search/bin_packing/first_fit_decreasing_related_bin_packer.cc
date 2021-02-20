#include "first_fit_decreasing_related_bin_packer.h"

#include <iostream>
#include <deque>
#include <unordered_map>
#include "bin_packer.h"
#include "../task_proxy.h"
#include "../causal_graph.h"
#include "../globals.h"
#include "../utils/timer.h"
#include "../utils/math.h"
#include "../utils/rng.h"

using namespace pdbs;
using namespace std;

FirstFitDecreasingRelatedBinPacker::FirstFitDecreasingRelatedBinPacker(double pdb_size, int num_collections) :
        pdb_max_size (pdb_size),
        num_pbd_collections (num_collections) {
}

vector<Pattern> FirstFitDecreasingRelatedBinPacker::bin_packing(std::shared_ptr<TaskProxy> task_proxy) {

    vector<Pattern> pattern_collection;
    double current_size;

    //int temp = rand()%(max_target_size-min_target_size+1);
    //temp += min_target_size;

    //pdb_max_size=9*pow(10,temp);
    //pdb_max_size=min(pdb_max_size,pow(10,initial_max_target_size));
    //pdb_max_size=max(pdb_max_size,pow(10,min_target_size));

    cout << "Starting bin packing First Fit Decresing Related, pdb_max_size:" << pdb_max_size << endl;

    VariablesProxy variables = task_proxy->get_variables();

    const CausalGraph &cg = task_proxy->get_causal_graph();

    deque<pair<int,double>> sorted_vars;
	unordered_map<int,double> remaining_vars;
	remaining_vars.reserve(variables.size());
    for (size_t i = 0; i < variables.size(); ++i) {
        if (variables[i].get_domain_size() <= pdb_max_size) {
			pair<int,double> p=make_pair(i, variables[i].get_domain_size());
            sorted_vars.push_back(p);
			remaining_vars.insert(p);
        }
    }

    sort(sorted_vars.begin(),sorted_vars.end(),compare_domain_size_variable);

    //Init pattern
    vector<bool> pattern(variables.size(), false);
    current_size = 1;

    size_t current_var = 0;
	unordered_map<int,double>::iterator itr;
    while (remaining_vars.size() > 0) {
		auto&pos = sorted_vars.front();
		
		if((itr=remaining_vars.find(pos.first))!=remaining_vars.end()) {
			if(utils::is_product_within_limit(current_size, pos.second, pdb_max_size)) {

				current_var = pos.first;
				current_size *= pos.second;
				pattern[current_var] = true;

				sorted_vars.pop_front();
				remaining_vars.erase(itr);

				const vector<int> &rel_vars2 = cg.get_eff_to_pre(current_var);
				vector<int> rel_vars=rel_vars2;
			    g_rng()->shuffle(rel_vars);
				unordered_map<int,double>::iterator itr2;
				for (int r:rel_vars) {
					if ((itr2=remaining_vars.find(r))!=remaining_vars.end() &&
						utils::is_product_within_limit(current_size, itr2->second, pdb_max_size)) {
						current_size *= itr2->second;
						pattern[r] = true;
						remaining_vars.erase(itr2);
					}
				}
			}
			else {

				//cout << pdb_max_size << endl;
				//cout << current_size << endl;

				//Add pattern (bin)
				vector<int> trans_pattern=transform_to_pattern_normal_form(pattern);
				pattern_collection.push_back(trans_pattern);


				//Init pattern
				pattern.clear();
				pattern.resize(variables.size(), false);
				current_size = 1;
			}
		}
		else	sorted_vars.pop_front();
    }

    if (current_size > 1) {
        //Add pattern (bin)
        vector<int> trans_pattern=transform_to_pattern_normal_form(pattern);
        pattern_collection.push_back(trans_pattern);
    }

    //sort(pattern_collection.begin(), pattern_collection.end(), compare_pattern_length);

    //cout << " binpacking time: " << utils::g_timer << " with " << pattern_collection.back().size() << endl;
    return pattern_collection;
}
