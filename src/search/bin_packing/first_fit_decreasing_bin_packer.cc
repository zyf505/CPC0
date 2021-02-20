#include "first_fit_decreasing_bin_packer.h"

#include <iostream>

#include "../task_proxy.h"
#include "../utils/timer.h"
#include "../utils/math.h"

using namespace pdbs;
using namespace std;

FirstFitDecreasingBinPacker::FirstFitDecreasingBinPacker(double pdb_size, int num_collections) :
        pdb_max_size (pdb_size),
        num_pbd_collections (num_collections) {
}

vector<Pattern> FirstFitDecreasingBinPacker::bin_packing(std::shared_ptr<TaskProxy> task_proxy) {

    vector<Pattern> pattern_collection;
    double current_size;

    //int temp = rand()%(max_target_size-min_target_size+1);
    //temp += min_target_size;

    //pdb_max_size=9*pow(10,temp);
    //pdb_max_size=min(pdb_max_size,pow(10,initial_max_target_size));
    //pdb_max_size=max(pdb_max_size,pow(10,min_target_size));

    cout << "Starting bin packing First Fit Decresing, pdb_max_size:" << pdb_max_size << endl;

    VariablesProxy variables = task_proxy->get_variables();

    vector<pair<int,double>> remaining_vars;

    for (VariableProxy var : variables) {
        if (var.get_domain_size() <= pdb_max_size) {
            remaining_vars.push_back(make_pair(var.get_id(), var.get_domain_size()));
        }
    }

    sort(remaining_vars.begin(), remaining_vars.end(), compare_domain_size_variable);

    //Init pattern
    vector<bool> pattern(variables.size(), false);
    current_size = 1;

    for (size_t j =0; j < remaining_vars.size(); j++) {

        if(utils::is_product_within_limit(current_size, remaining_vars[j].second, pdb_max_size)) {
            current_size *= remaining_vars[j].second;
            pattern[remaining_vars[j].first] = true;
        }
        else {
            //Add pattern (bin)
            pattern_collection.push_back(transform_to_pattern_normal_form(pattern));

            vector<int> trans_pattern = transform_to_pattern_normal_form(pattern);

            //Init pattern
            pattern.clear();
            pattern.resize(variables.size(), false);
            current_size = 1;
        }
    }

    if (current_size > 1) {
        //Add pattern (bin)
        pattern_collection.push_back(transform_to_pattern_normal_form(pattern));
        vector<int> trans_pattern=transform_to_pattern_normal_form(pattern);
    }

    sort(pattern_collection.begin(), pattern_collection.end(), compare_pattern_length);

    //cout << " binpacking time: " << utils::g_timer << " with " << pattern_collection.back().size() << endl;
    return pattern_collection;
}
