//#include "pattern_collection_generator_complementary.h"
#include "pattern_collection_evaluator_RandWalk_Avg_h.h"
#include "../sampling.h"
#include "../utils/timer.h"
#include "../task_tools.h"
#include "../utils/countdown_timer.h"
#include "../options/option_parser.h"
#include "../options/plugin.h"

//#include "../causal_graph.h"
//#include "../globals.h"
//#include "../task_proxy.h"

//#include "../utils/markup.h"
//#include "../utils/math.h"
//#include "../utils/rng.h"
//#include "../utils/timer.h"
//#include "../heuristic.h"

//#include <algorithm>
//#include <cassert>
//#include <iostream>
//#include <unordered_set>
//#include <vector>
//#include <math.h>
//Hack to use SS get_type, it needs heuristic object in constructor
//#include "../heuristics/blind_search_heuristic.h"
//#include "../heuristics/lm_cut_heuristic.h"
#include "pdb_factory.h"
//#include "pattern_database_interface.h"
#include "../utils/debug_macros.h"
//#include <random>
//#include "../sampling.h"
#include "../state_registry.h"
#include <climits>
#include <boost/range/adaptor/reversed.hpp>
    
using namespace std;
//PatternCollectionEvaluator is to be the driver for 8 PDB-based options
//RandomCollectionGeneration: CBP, RBP, CGamer, the one used by *Pommerening et al. 
//Local Search: iPDB, gaPDB, CGamer, VPN, *none, *changing order of patterns in gaPDB mutation. (this one matters to 0-1 greedy cost partitioning).
//GenPDB: Symbolic, Explicit, Online, expressed on pdb_factory class
//PDBEval: AvgH, Random sampling, Stratified sampling, *original iPDB method
//CombPDBs->Canonical, hPO, Max
//CostPartition->*None, Saturated, 0-1 greedy 
//Learning: UCB1 to choose bin packing, pdb size. 
//Re-evaluate: None, RemovedPDBsDominated, Run GHS (note: I think if we run GHS here, it may help to generate PDBs that are complementary to LM-cut)
// And the pseudo-code logic for it:
//While (time < 900 seconds)
//     PC<-RandomCollectionGeneration(MaxSize,CostPartition)
//     setInterestingPCs<-LocalSearch(P,PDBEval,MaxSize,CostPartition)
//     selectedPCs<-SubsetSelection(setInterestingPCs,PDBEval,ComPDBs)
//     generatedPDBs <-GenPDB(selectedPCs)
//     H<-Re-evaluate (H¿generatedPCs)
//     Learning(BinPackingRewards)
//ComPDBs (Hl)
namespace pdbs3 {
PatterCollectionEvaluatorRandWalk_Avg_H::PatterCollectionEvaluatorRandWalk_Avg_H(const options::Options & opts) :
	time_limit (opts.get<int>("time_limit")){
    cout<<"hello EvaluatorRandWalk_Avg_H"<<flush<<endl;
}
  void PatterCollectionEvaluatorRandWalk_Avg_H::initialize(std::shared_ptr<AbstractTask> task) {
    /*int num_vars= task->get_num_variables();
    cout<<"num_vars:"<<num_vars<<flush<<endl;*/
    TaskProxy task_proxy_temp(*task);
    task_proxy=make_shared<TaskProxy>(task_proxy_temp);
    successor_generator=utils::make_unique_ptr<SuccessorGenerator>(task);
    //result=make_shared<PatternCollectionInformation>(task, make_shared<PatternCollection>());
//    cout << "Manual pattern collection: " << *patterns << endl;
//    return PatternCollectionInformation(task, patterns);
  }
   bool PatterCollectionEvaluatorRandWalk_Avg_H::evaluate(std::shared_ptr<ModularZeroOnePDBs> candidate_PC){
    //cout<<"candidate_pc.size:"<<candidate_PC->get_size()<<endl;
	increased_states=0;
	size_t additions=0;
	unsigned count_dead_ends=0;
	int new_val;
	std::map<size_t, std::pair<State,int> >::iterator itr = unique_samples.begin();
    while (itr != unique_samples.end()) {
		new_val=candidate_PC->get_value(itr->second.first);
		if(new_val!=numeric_limits<int>::max()) {
			if(new_val>itr->second.second){
				increased_states++;
				DEBUG_MSG(cout<<"\th improved from "<<itr->second.second<<" to "<<new_val<<endl;);
			}
			itr++;
			additions++;
		}
		else	{
			itr=unique_samples.erase(itr);//Found new dead_end	
			count_dead_ends++;
		}
	}
	if(additions==0){
        cerr<<"No additions for calculating ratio, debug me!!!"<<endl;
		exit(1);
      }
    //UPDATING UNIQUE_SAMPLES h VALUES ONLY IF COLLECTION WILL BE ADDED
	double ratio=increased_states/double(additions);
    if(ratio>=0.25){
      cout<<"time:"<<utils::g_timer()<<",Selecting PC,increased_states:"<<increased_states<<", out of "<<additions<<", ratio:"<<ratio<<endl;
      /*for(std::map<size_t,std::pair<State,int> >::iterator it=unique_samples.begin(); it!=unique_samples.end(); ++it){
        it->second.second=max(it->second.second,candidate_PC->get_value(it->second.first));
      }*/

      itr = unique_samples.begin();
      int new_h;
      while (itr != unique_samples.end()) {
        new_h=candidate_PC->get_value(itr->second.first);
        if (itr->second.second==numeric_limits<int>::max()){
          itr = unique_samples.erase(itr);
          count_dead_ends++;
        }
        else if(new_h==numeric_limits<int>::max()){
          itr = unique_samples.erase(itr);
          count_dead_ends++;
        }
        else {
          itr->second.second=max(itr->second.second,new_h);
          ++itr;
        }
      }
      cout<<"deleted "<<count_dead_ends<<" from list of sampled_states, remaining states:"<<unique_samples.size()<<endl; 
      return true;//Add collection
    }
    DEBUG_COMP(cout<<"time:"<<utils::g_timer()<<",Not_selecting PC,increased_states:"<<increased_states<<",threshold:"<<get_threshold()<<" out of "<<get_num_samples()<<endl;);
    //cout<<"time:"<<utils::g_timer()<<",Not_selecting PC,increased_states:"<<increased_states<<",threshold:"<<get_threshold()<<" out of "<<get_num_samples()<<endl;
	//cout<<"deleted "<<count_dead_ends<<" from list of sampled_states, remaining states:"<<unique_samples.size()<<endl; 
    return false;//Not adding collection
  }
  double PatterCollectionEvaluatorRandWalk_Avg_H::calculate_eval_score(std::shared_ptr<ModularZeroOnePDBs> candidate_PC){
    //cout<<"candidate_pc.size:"<<candidate_PC->get_size()<<endl;
    //increased_states=0;
    long sum_h=0;
    size_t additions=0;
	unsigned count_dead_ends=0;
	int new_val;
	std::map<size_t, State>::iterator itr = unique_samples_Gamer.begin();
    while(itr!=unique_samples_Gamer.end()) {
		new_val=candidate_PC->get_value(itr->second);
        //increased_states++;
        if (new_val!=numeric_limits<int>::max())	{
			sum_h+=new_val;
			additions++;
			itr++;
		}
		else	{
			itr=unique_samples_Gamer.erase(itr);//Found new dead_end
			count_dead_ends++;
		}
		//if(additions>=states_to_eval) break;
      }
	 if(additions==0){
        cerr<<"No additions for calculating avg_h_val, debug me!!!"<<endl;
		exit(1);
      }
	double result=sum_h/double(additions);
    //set_eval_score(result);
	//cout<<"deleted "<<count_dead_ends<<" from list of sampled_states for Gamer, remaining states:"<<unique_samples_Gamer.size()<<endl; 
	return result; 
  }

  void PatterCollectionEvaluatorRandWalk_Avg_H::sample_states(std::shared_ptr<PatternCollectionInformation> current_result){
    long sum_h=0;
    size_t additions=0;
    //samples.clear();//We only use as samples the states, ordered in the open they would be fetched, in the current open list
    //Need to keep pointer to result or sample_states... function will complain current_result is not captured
	unique_samples.clear();
    result=current_result;
    float start_time=utils::g_timer();
    //DEBUG_MSG(cout<<"adding to samples, unique_size prior:"<<unique_samples.size()<<",sampled_states:"<<samples.size()<<endl;);


// IF STATES_LOADED_FROM_OPEN_LIST IS EMPTY IT MEANS THIS IS FIRST CALL, REAL SEARCH HAS NOT STARTED
// FOR INITIALIZATION PURPOSES WE DO SAME RANDOM WALK AS IN RandWalk class

      const State &initial_state = task_proxy->get_initial_state();
      int init_h=current_result->get_value(initial_state);
      double average_operator_cost=get_average_operator_cost(*task_proxy);
      evaluator_timer = new utils::CountdownTimer(time_limit);
      cout<<"calling sample_states_with_random_walks"<<flush<<endl;
	  auto const &samples_just_states = sample_states_with_random_walks(
	      *task_proxy, *successor_generator, get_num_samples(), init_h,
	      average_operator_cost,
	      [this](const State &state) {
		  return result->is_dead_end(state);
	      },
	      evaluator_timer);
	  cout<<"We are finished,random_walk_time:"<<utils::g_timer()-start_time<<endl;
	  int h_val;
      for (auto state : samples_just_states){
		h_val=current_result->get_value(state);
        if (h_val==numeric_limits<int>::max())//Skipping dead_ends when doing avg_h
			continue;
		additions++;
		sum_h+=h_val;
		auto const &ret=unique_samples.insert(make_pair(state.hash(),make_pair(state,h_val)));
		if(ret.second==false)//Node already exist so update h_value
			ret.first->second.second=max(h_val,ret.first->second.second);
		//samples.push_back(make_pair(state,h_val));
      }
	  h_val=current_result->get_value(initial_state);  //insert initial state into samples
	  additions++;
	  sum_h+=h_val;
	  auto const &ret=unique_samples.insert(make_pair(initial_state.hash(),make_pair(initial_state,h_val)));
	  if(ret.second==false)//Node already exist so update h_value
		    ret.first->second.second=max(h_val,ret.first->second.second);
	  if(additions==0){
        cerr<<"Unique sample is empty, debug me!!!"<<endl;
		exit(1);
      }
      //set_num_samples(samples.size());
      set_sample_score(sum_h/double(additions));
      //cout<<"Sampled_score:"<<get_sample_score()<<endl;
      cout<<"adding to samples using RAND_WALK, unique_size:"<<unique_samples.size()<<flush<<",num_samples:,"<<get_num_samples()<<",sample_score:,"<<get_sample_score()<<endl;
      result.reset();
	  return;
  }
  
  void PatterCollectionEvaluatorRandWalk_Avg_H::sample_states_Gamer(std::shared_ptr<PatternCollectionInformation> current_result,int init_h){
	unique_samples_Gamer.clear();
	result=current_result;
    float start_time=utils::g_timer();
    //DEBUG_MSG(cout<<"adding to samples, unique_size prior:"<<unique_samples.size()<<",sampled_states:"<<samples.size()<<endl;);

      double average_operator_cost=get_average_operator_cost(*task_proxy);
      evaluator_timer = new utils::CountdownTimer(time_limit);
      cout<<"calling sample_states_with_random_walks_Gamer, init_h: "<<max(init_h,1)<<flush<<endl;
	  auto const &samples_just_states = sample_states_with_random_walks(
	      *task_proxy, *successor_generator, get_num_samples(), max(init_h,1),
	      average_operator_cost,
	      [this](const State &state) {
		  return result->is_dead_end(state);
	      },
	      evaluator_timer);
	  cout<<"We are finished,random_walk_time:"<<utils::g_timer()-start_time<<endl;
      for (auto state : samples_just_states)    unique_samples_Gamer.insert(make_pair(state.hash(),state)); 
	  const State &initial_state = task_proxy->get_initial_state();
	  unique_samples_Gamer.insert(make_pair(initial_state.hash(),initial_state));
      cout<<"adding to samples using RAND_WALK, unique_size:"<<unique_samples_Gamer.size()<<flush<<",num_samples:,"<<get_num_samples()<<endl;
      result.reset();
	  return;
  }
  
  void PatterCollectionEvaluatorRandWalk_Avg_H::clear_dominated_heuristics(const MaxAdditivePDBSubsets &current_max_additive_subsets,const std::shared_ptr<PatternCollectionInformation> &new_result,const PDBCollection &candidate_pdb) {
	  double start_time=utils::g_timer();
	  vector<int> current_best_h_values;
	  current_best_h_values.reserve(unique_samples.size());


	  cout<<"time:"<<utils::g_timer()<<",calling clear_dominated_heuristics with "<<current_max_additive_subsets.size()+1<<" best heuristics and unique_samples:"<<unique_samples.size()<<endl;
	  //shared_ptr<MaxAdditivePDBSubsets> current_max_additive_subsets=current_result->get_max_additive_subsets();
	  
	  //First we get all the values for sampled states with the latest PC
	  int i=0;
	  //const State &initial_state = task_proxy->get_initial_state();
	  //int init_h=calculate_max_additive_subset(candidate_pdb,initial_state);
	  for(map<size_t,pair<State,int> >::iterator it=unique_samples.begin(); it!=unique_samples.end();it++){
		current_best_h_values[i++]=calculate_max_additive_subset(candidate_pdb,it->second.first);
	  }
	  //cout<<"\ttime to populate current_best_h_values:"<<utils::g_timer()-start_time<<",unique_samples:"<<i<<endl;
		new_result->include_additive_pdbs(make_shared<PDBCollection>(candidate_pdb));

	  //Now go through each additive subset and check if they are dominated for the sampled set of states
	  //i=1;
	  int h;
	  //int new_init_h;
	  for(auto const &additive_subset : boost::adaptors::reverse(current_max_additive_subsets)){
		bool dominated_heur=true;
		/*new_init_h=calculate_max_additive_subset(additive_subset,initial_state);
		if(new_init_h>init_h) {
			init_h=new_init_h;
			dominated_heur=false;
		}*/
		int j=0;
	   for(map<size_t,pair<State,int> >::iterator it=unique_samples.begin(); it!=unique_samples.end();it++){
	  //If dead_end skip
		if(current_best_h_values[j]==INT_MAX){
			j++;
			continue;
		}
		h=calculate_max_additive_subset(additive_subset,it->second.first);
		  //NO BREAKS BECAUSE WE WANT TO CALCULATE ALL THE NEW HIGHER H VALUES
		  //IF HEUR IS NOT DOMINATED
		if (h == numeric_limits<int>::max()){
			dominated_heur=false;
			//cout<<"\tcolleciton ["<<i<<" is undominated because of dead_end, prev_val:"<<current_best_h_values[j]<<",h:"<<h<<endl;
			current_best_h_values[j]=INT_MAX;
		}
		else if(h>current_best_h_values[j]){
			dominated_heur=false;
			//cout<<"\tcolleciton ["<<i<<" is undominated because of higher_h, prev_val:"<<current_best_h_values[j]<<",h:"<<h<<endl;
			current_best_h_values[j]=h;
		}
		j++;
		}

		if(!dominated_heur){
		  //cout<<"adding heur["<<i<<"] to list of heurs"<<endl;	 
		  new_result->include_additive_pdbs(make_shared<PDBCollection>(additive_subset));
		}
		/*else{
		  cout<<"collection["<<i<<"] is dominated,eliminating "<<endl;
		}
		i++;*/
		}
		auto &new_max_additive_subsets=*(new_result->get_max_additive_subsets());
		std::reverse(new_max_additive_subsets.begin(),new_max_additive_subsets.end());
		cout<<"clear_dominated_heuristics,time_spent,"<<utils::g_timer()-start_time<<",Reduced number of collections from:,"<<current_max_additive_subsets.size()+1<<",to:,"<<new_result->get_max_additive_subsets()->size()<<endl;
  }

  static shared_ptr<PatternCollectionEvaluator>_parse(options::OptionParser &parser) {
    parser.add_option<int> ("time_limit", "If populated,stop construction on first node past boundary and time limit", "100");
    options::Options options = parser.parse();
    parser.document_synopsis(
        "Pattern Generator RBP",
        "RBP-stype selection of variables to generate Pattern Collection");
    options::Options opts = parser.parse();
    if (parser.dry_run())
        return 0;

    return make_shared<PatterCollectionEvaluatorRandWalk_Avg_H>(opts);
  }

  static options::PluginShared<PatternCollectionEvaluator> _plugin("rand_walk_evaluator_avg_h", _parse);
}
