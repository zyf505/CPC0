#include "pattern_database.h"

#include "match_tree.h"

#include "../priority_queue.h"
#include "../task_tools.h"

#include "../utils/collections.h"
#include "../utils/logging.h"
#include "../utils/math.h"
#include "../utils/timer.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <climits>
#include <string>
#include <vector>

using namespace std;

namespace pdbs3 {
AbstractOperator::AbstractOperator(const vector<pair<int, int>> &prev_pairs,
                                   const vector<pair<int, int>> &pre_pairs,
                                   const vector<pair<int, int>> &eff_pairs,
                                   int cost,
                                   const vector<size_t> &hash_multipliers)
    : cost(cost),
      regression_preconditions(prev_pairs) {
    regression_preconditions.insert(regression_preconditions.end(),
                                    eff_pairs.begin(),
                                    eff_pairs.end());
    // Sort preconditions for MatchTree construction.
    sort(regression_preconditions.begin(), regression_preconditions.end());
    for (size_t i = 1; i < regression_preconditions.size(); ++i) {
        assert(regression_preconditions[i].first !=
               regression_preconditions[i - 1].first);
    }
    hash_effect = 0;
    assert(pre_pairs.size() == eff_pairs.size());
    for (size_t i = 0; i < pre_pairs.size(); ++i) {
        int var = pre_pairs[i].first;
        assert(var == eff_pairs[i].first);
        int old_val = eff_pairs[i].second;
        int new_val = pre_pairs[i].second;
        assert(new_val != -1);
        size_t effect = (new_val - old_val) * hash_multipliers[var];
        hash_effect += effect;
    }
}

AbstractOperator::~AbstractOperator() {
}

void AbstractOperator::dump(const Pattern &pattern,
                            const TaskProxy &task_proxy) const {
    cout << "AbstractOperator:" << endl;
    cout << "Regression preconditions:" << endl;
    for (size_t i = 0; i < regression_preconditions.size(); ++i) {
        int var_id = regression_preconditions[i].first;
        int val = regression_preconditions[i].second;
        cout << "Variable: " << var_id << " (True name: "
             << task_proxy.get_variables()[pattern[var_id]].get_name()
             << ", Index: " << i << ") Value: " << val << endl;
    }
    cout << "Hash effect:" << hash_effect << endl;
}

PatternDatabase::PatternDatabase(
    const TaskProxy &task_proxy,
    const Pattern &pattern,
    bool dump,
    int time_limit,
    const vector<int> &operator_costs)
    : PatternDatabaseInterface (task_proxy, pattern, operator_costs) {
    verify_no_axioms(task_proxy);
    verify_no_conditional_effects(task_proxy);
    assert(operator_costs.empty() ||
           operator_costs.size() == task_proxy.get_operators().size());
    assert(utils::is_sorted_unique(pattern));
    assert(!pattern.empty());

    utils::Timer timer;
    hash_multipliers.reserve(pattern.size());
    num_states = 1;
    for (int pattern_var_id : pattern) {
        hash_multipliers.push_back(num_states);
        VariableProxy var = task_proxy.get_variables()[pattern_var_id];
        if (utils::is_product_within_limit(num_states, var.get_domain_size(),
                                           numeric_limits<int>::max())) {
            num_states *= var.get_domain_size();
        } else {
            cerr << "Given pattern is too large! (Overflow occured): " << endl;
            cerr << pattern << endl;
            utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
        }
    }
    //cout<<"calling create_pdb with num_states="<<num_states<<endl;fflush(stdout);dump=true;
    create_pdb(time_limit);
    if (dump)
        cout << "PDB construction time: " << timer << endl;
}

void PatternDatabase::multiply_out(
    int pos, int cost, vector<pair<int, int>> &prev_pairs,
    vector<pair<int, int>> &pre_pairs,
    vector<pair<int, int>> &eff_pairs,
    const vector<pair<int, int>> &effects_without_pre,
    vector<AbstractOperator> &operators) {
    if (pos == static_cast<int>(effects_without_pre.size())) {
        // All effects without precondition have been checked: insert op.
        if (!eff_pairs.empty()) {
            operators.push_back(
                AbstractOperator(prev_pairs, pre_pairs, eff_pairs, cost,
                                 hash_multipliers));
        }
    } else {
        // For each possible value for the current variable, build an
        // abstract operator.
        int var_id = effects_without_pre[pos].first;
        int eff = effects_without_pre[pos].second;
        VariableProxy var = task_proxy.get_variables()[pattern[var_id]];
        for (int i = 0; i < var.get_domain_size(); ++i) {
            if (i != eff) {
                pre_pairs.push_back(make_pair(var_id, i));
                eff_pairs.push_back(make_pair(var_id, eff));
            } else {
                prev_pairs.push_back(make_pair(var_id, i));
            }
            multiply_out(pos + 1, cost, prev_pairs, pre_pairs, eff_pairs,
                         effects_without_pre, operators);
            if (i != eff) {
                pre_pairs.pop_back();
                eff_pairs.pop_back();
            } else {
                prev_pairs.pop_back();
            }
        }
    }
}

void PatternDatabase::build_abstract_operators(
    const OperatorProxy &op, int cost,
    const std::vector<int> &variable_to_index,
    vector<AbstractOperator> &operators) {
    // All variable value pairs that are a prevail condition
    vector<pair<int, int>> prev_pairs;
    // All variable value pairs that are a precondition (value != -1)
    vector<pair<int, int>> pre_pairs;
    // All variable value pairs that are an effect
    vector<pair<int, int>> eff_pairs;
    // All variable value pairs that are a precondition (value = -1)
    vector<pair<int, int>> effects_without_pre;

    size_t num_vars = task_proxy.get_variables().size();
    vector<bool> has_precond_and_effect_on_var(num_vars, false);
    vector<bool> has_precondition_on_var(num_vars, false);

    for (FactProxy pre : op.get_preconditions())
        has_precondition_on_var[pre.get_variable().get_id()] = true;

    for (EffectProxy eff : op.get_effects()) {
        int var_id = eff.get_fact().get_variable().get_id();
        int pattern_var_id = variable_to_index[var_id];
        int val = eff.get_fact().get_value();
        if (pattern_var_id != -1) {
            if (has_precondition_on_var[var_id]) {
                has_precond_and_effect_on_var[var_id] = true;
                eff_pairs.push_back(make_pair(pattern_var_id, val));
            } else {
                effects_without_pre.push_back(make_pair(pattern_var_id, val));
            }
        }
    }
    for (FactProxy pre : op.get_preconditions()) {
        int var_id = pre.get_variable().get_id();
        int pattern_var_id = variable_to_index[var_id];
        int val = pre.get_value();
        if (pattern_var_id != -1) { // variable occurs in pattern
            if (has_precond_and_effect_on_var[var_id]) {
                pre_pairs.push_back(make_pair(pattern_var_id, val));
            } else {
                prev_pairs.push_back(make_pair(pattern_var_id, val));
            }
        }
    }
    multiply_out(0, cost, prev_pairs, pre_pairs, eff_pairs, effects_without_pre,
                 operators);
}

void PatternDatabase::create_pdb(int time_limit) {
  //cout<<"calling create_pdb_explicit"<<endl;
    VariablesProxy vars = task_proxy.get_variables();
    vector<int> variable_to_index(vars.size(), -1);
    finished=true;
    for (size_t i = 0; i < pattern.size(); ++i) {
        variable_to_index[pattern[i]] = i;
    }
    //cout<<"time_limit(ms):"<<time_limit<<endl;

    // compute all abstract operators
    vector<AbstractOperator> operators;
    for (OperatorProxy op : task_proxy.get_operators()) {
        int op_cost;
        if (operator_costs.empty()) {
            op_cost = op.get_cost();
        } else {
            op_cost = operator_costs[op.get_id()];
        }
        build_abstract_operators(op, op_cost, variable_to_index, operators);
    }

    // build the match tree
    MatchTree match_tree(task_proxy, pattern, hash_multipliers);
    for (const AbstractOperator &op : operators) {
        match_tree.insert(op);
    }

    // compute abstract goal var-val pairs
    vector<pair<int, int>> abstract_goals;
    for (FactProxy goal : task_proxy.get_goals()) {
        int var_id = goal.get_variable().get_id();
        int val = goal.get_value();
        if (variable_to_index[var_id] != -1) {
            abstract_goals.push_back(make_pair(variable_to_index[var_id], val));
        }
    }

    distances.reserve(num_states);
    // first implicit entry: priority, second entry: index for an abstract state
    AdaptiveQueue<size_t> pq;

    // initialize queue
    for (size_t state_index = 0; state_index < num_states; ++state_index) {
        if (is_goal_state(state_index, abstract_goals)) {
            pq.push(0, state_index);
            distances.push_back(0);
        } else {
            distances.push_back(numeric_limits<int>::max());
        }
    }

    int unseen_value=1;
    // Dijkstra loop
    int counter=0;
    //cout<<"pq empty:"<<pq.empty()<<",time_limit:"<<time_limit<<endl;
  
    //Best to start counting from here, otherwise we never get to do anything with PDBs expensive to set up
    //when using small time limits
    double start_create_time=utils::g_timer();
    while (!pq.empty()) {
        pair<int, size_t> node = pq.pop();
        int distance = node.first;
	//cout<<"\t\t distance:"<<distance<<endl;
        size_t state_index = node.second;
        if (distance > distances[state_index]) {
            continue;
        }
	if(time_limit){
	  //for those cases where there are so many goal states we never went past distance==0
	  if(distance==0){
	    if(counter%10000==0){
	      if((utils::g_timer()-start_create_time)*1000>4*time_limit){
		//cout<<"\t\tinterrupting explicit pdb generation,time taken:"<<utils::g_timer()-start_create_time<<",time_limit(ms):"<<time_limit<<",distance:"<<distance<<",counter:"<<counter<<endl;
		unseen_value=0;
		finished=false;
		break;
	      }
	    }
	  }
	  else if(distance>unseen_value){
	    unseen_value=distance;
	    if((utils::g_timer()-start_create_time)*1000>4*time_limit){
	      //cout<<"\t\tinterrupting explicit pdb generation,time taken:"<<utils::g_timer()-start_create_time<<",time_limit(ms):"<<time_limit<<",distance:"<<distance<<endl;
	      finished=false;
	      break;
	    }
	    //else{
	    //  cout<<"\t\tNot interrupting explicit pdb generation,time taken:"<<utils::g_timer()-start_create_time<<",time_limit(ms):"<<time_limit<<",distance:"<<distance<<endl;
	    //}
	  }
	}


        // regress abstract_state
        vector<const AbstractOperator *> applicable_operators;
        match_tree.get_applicable_operators(state_index, applicable_operators);
        for (const AbstractOperator *op : applicable_operators) {
            size_t predecessor = state_index + op->get_hash_effect();
            int alternative_cost = distances[state_index] + op->get_cost();
            if (alternative_cost < distances[predecessor]) {
                distances[predecessor] = alternative_cost;
                pq.push(alternative_cost, predecessor);
            }
        }
	counter++;
    }

    if(!finished){
      //cout<<"\t\tnow marking all unseen states to unseen value of:"<<unseen_value<<endl;
      for (auto& dist : distances){
	if(dist>unseen_value){
	  dist=unseen_value;
	}
      }
    }
    //else{
    //  cout<<"pdb took to fully generate:"<<(utils::g_timer()-start_create_time)*1000<<"ms"<<endl;
    //}
    //cout<<"pdb finished:"<<finished<<endl;
}

bool PatternDatabase::is_goal_state(
    const size_t state_index,
    const vector<pair<int, int>> &abstract_goals) const {
    for (pair<int, int> abstract_goal : abstract_goals) {
        int pattern_var_id = abstract_goal.first;
        int var_id = pattern[pattern_var_id];
        VariableProxy var = task_proxy.get_variables()[var_id];
        int temp = state_index / hash_multipliers[pattern_var_id];
        int val = temp % var.get_domain_size();
        if (val != abstract_goal.second) {
            return false;
        }
    }
    return true;
}

size_t PatternDatabase::hash_index(const State &state) const {
    size_t index = 0;
    for (size_t i = 0; i < pattern.size(); ++i) {
        index += hash_multipliers[i] * state[pattern[i]].get_value();
    }
    return index;
}

int PatternDatabase::get_value(const State &state) const {
  //if(distances[hash_index(state)]>1){
    //cout<<"Pattern:"<<pattern<<",initial_state:,"<<hash_index(state)<<",h:"<<distances[hash_index(state)]<<endl;
  //}
    return distances[hash_index(state)];
}

size_t PatternDatabase::hash_index(const vector<int> &state) const {
    size_t index = 0;
    for (size_t i = 0; i < pattern.size(); ++i) {
        index += hash_multipliers[i] * state[pattern[i]];
    }
    return index;
}

int PatternDatabase::get_value(const vector<int> &state) const {
    return distances[hash_index(state)];
}

double PatternDatabase::compute_mean_finite_h() const {
    double sum = 0;
    int size = 0;
    for (size_t i = 0; i < distances.size(); ++i) {
        if (distances[i] != numeric_limits<int>::max()) {
            sum += distances[i];
            ++size;
        }
    }
    if (size == 0) { // All states are dead ends.
        return numeric_limits<double>::infinity();
    } else {
        return sum / size;
    }
}

int PatternDatabase::compute_heuristic_id(size_t state_id) {
  //cout<<"calling offline_compute_heuristic_id"<<endl;fflush(stdout);
  //cout<<"state_id="<<state_id<<",entries:"<<num_states<<endl;fflush(stdout);
    int h = distances[state_id];
    //cout<<"h_offline:"<<h<<endl;fflush(stdout);
    if (h == numeric_limits<int>::max())
        return INT_MAX/2;//Better when doing maxes than DEAD_END=-1!!!!
    return h;
}

}
