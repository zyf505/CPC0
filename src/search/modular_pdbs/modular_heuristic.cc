#include "modular_heuristic.h"
#include "types.h"

//#include "pattern_generator.h"
#include "pattern_collection_generator_RBP.h"
#include "pattern_collection_generator_bin_packing_v1.h"
#include "pattern_collection_evaluator_RandWalk.h"
#include "pdb_factory_symbolic.h"
//#include "pattern_database_interface.h"

#include "../option_parser.h"
#include "../plugin.h"
#include "../task_proxy.h"
#include "../globals.h"

#include <limits>
#include <memory>
#include <climits>
#include "zero_one_pdbs.h"
#include "../utils/countdown_timer.h"
#include "learning.h"
#include "pattern_database_symbolic.h"
#include "pattern_collection_generator_GamerStyle.h"
#include <boost/range/adaptor/reversed.hpp>
//#include "pdb_factory_symbolic.h"
//#include "pdb_factory.h"

using namespace std;
std::ostream & operator<<(std::ostream &os,const set<int> &pattern){
    for (int v : pattern) os << "," << v;  
    return os;
}
std::ostream & operator<<(std::ostream &os,const vector<int> &pattern){
    for (int v : pattern) os << "," << v;  
    return os;
}
/*template<class T>
std::ostream & operator<<(ostream& stream, const std::vector<T>& values)
{
    stream << "[ ";
    copy( begin(values), end(values), ostream_iterator<T>(stream, " ") );
    stream << "] ";
    return stream;
}*/

//bool recompute_additive_sets=false;

//PatternCollectionGeneratorComplementary is to be the driver for 8 PDB-based options
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
	
//Pattern get_pattern_from_options(const shared_ptr<AbstractTask> task,
//				 const Options &opts) {
//    shared_ptr<PatternGenerator> pattern_generator =
//        opts.get<shared_ptr<PatternGenerator>>("pattern");
//    return pattern_generator->generate(task);
//}

ModularHeuristic::ModularHeuristic(const Options &opts)
    : Heuristic(opts),
    pattern_generator(opts.get<shared_ptr<PatternCollectionGeneratorComplementary>>("patterns")),
    pattern_evaluator(opts.get<shared_ptr<PatternCollectionEvaluator>>("evaluator")),
	  modular_time_limit(opts.get<int>("modular_time_limit")),
	  terminate_creation(opts.get<bool>("terminate_creation")),
	  create_perimeter(opts.get<bool>("create_perimeter")), 
	  only_gamer(opts.get<bool>("only_gamer")), 
	  only_CBP(opts.get<bool>("only_cbp")), 
      pdb_factory (opts.get<shared_ptr<PDBFactory>>("pdb_factory")) {
      cout<<"Hi modular_v1"<<endl;
      //bool unterminated_pdbs=false;
      cout<<"modular_time_limit:"<<modular_time_limit<<endl;
      cout<<"terminate_creation:"<<terminate_creation<<endl;
      cout<<"pdb_type:"<<pdb_factory->name()<<endl;
      cout<<"only_gamer:"<<only_gamer<<endl;
      cout<<"only_CBP:"<<only_CBP<<endl;
      int generator_choice=0;
	  double pdb_max_size=0;
      bool disjunctive_choice=false;
      unsigned num_goals_to_group=0;
      modular_heuristic_timer = new utils::CountdownTimer(modular_time_limit);
      TaskProxy task_proxy(*task);
	  //int time_to_prun=500;
      int initial_h=0;
      int new_initial_h;
	  int sampling_h=0;
      int num_episodes=0;
      int PC_counter=0;
      //vector<shared_ptr<ModularZeroOnePDBs> > pdb_ptr_collection;
      Options opts2=opts;
      //Options opts3=opts;
      Options opts4=opts;
      Options opts5=opts;
      opts2.set<int>("time_limit",40);
      //opts3.set<int>("packer_selection",0);
      opts4.set<int>("packer_selection",1);
      opts5.set<int>("packer_selection",2);
      PatternCollectionGeneratorGamer alternative_pattern_generator(opts2);
      


      vector<pair<double,Pattern>> improving_patterns;
      double best_value=0;//For Gamer-Style selection
	  int best_h;
	  int new_best_h;
      //0 means Random split into two patterns, 1 means CBP
      Learning UCB_generator; UCB_generator.insert_choice(1); UCB_generator.insert_choice(0);
      //Initialize reward to 2 so it does not go too fast for initial selection
      //UCB_generator.increase_reward(0);UCB_generator.increase_reward(1);

      Learning UCB_sizes;//PDB size selector
      //map<double,Learning> UCB_Disjunctive_patterns;//Disjunctive (or not) pattern selector,one per size;
      //UCB_sizes.insert_choice(pow(10,4));
      //UCB_sizes.insert_choice(8);
      //Learning binary_choice;binary_choice.insert_choice(1);binary_choice.insert_choice(0);
      //Initialize reward to 2 so it does not go too fast for initial selection
      //binary_choice.increase_reward(1);binary_choice.increase_reward(0);
      size_t overall_num_goals=task_proxy.get_goals().size();
      cout<<"overall_num_goals:"<<overall_num_goals<<endl;
      bool partial_gamer_run=false;//In case we did not process all available variabes due to taking too long
      //we need to know so that if Gamer come back with a partial list and none was used then we do not block
      //more Gamer turns unless none of those partial runs resulted on adding a variable
      //bool partial_gamer_improvement_found = true;
	  //bool backtrack=false;
      //bool finished=false;  //to decide whether terminate Gamer
	  //bool failure=false;
	  //bool terminate=false;  //to decide whether terminate CBP
	  //int failures_in_a_row=0;
	  //int CBP_extra=0;
      Learning goals_choice;
      for (size_t i=0;i<overall_num_goals;i++){
        goals_choice.insert_choice(i+1);
        //Initialize reward to 2 so it does not go too fast for initial selection
        //goals_choice.increase_reward(i+1);
      }
      //map<double,Learning> UCB_goals_to_group;//Disjunctive (or not) pattern selector,one per size;
      
      //Learning terminate_choice;terminate_choice.insert_choice(1);terminate_choice.insert_choice(0);
      //terminate_choice.increase_reward(0,10);//biasing towards not terminating in the begining

      UCB_sizes.insert_choice(9);
      /*UCB_sizes.insert_choice(10);
      UCB_sizes.insert_choice(11);
      UCB_sizes.insert_choice(12);
      UCB_sizes.insert_choice(13);
      UCB_sizes.insert_choice(14);
      UCB_sizes.insert_choice(15);
      UCB_sizes.insert_choice(16);
      UCB_sizes.insert_choice(17);
      UCB_sizes.insert_choice(18);
      UCB_sizes.insert_choice(19);
      UCB_sizes.insert_choice(20);
      UCB_sizes.insert_choice(25);
      UCB_sizes.insert_choice(30);
      UCB_sizes.insert_choice(35);*/

      /*UCB_Disjunctive_patterns[pow(10,8)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,9)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,10)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,11)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,12)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,13)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,14)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,15)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,16)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,17)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,18)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,19)]=binary_choice;
      //UCB_Disjunctive_patterns[pow(10,20)]=binary_choice;
      
      UCB_goals_to_group[pow(10,8)]=goals_choice;
      UCB_goals_to_group[pow(10,9)]=goals_choice;
      UCB_goals_to_group[pow(10,10)]=goals_choice;
      UCB_goals_to_group[pow(10,11)]=goals_choice;
      UCB_goals_to_group[pow(10,12)]=goals_choice;
      UCB_goals_to_group[pow(10,13)]=goals_choice;
      UCB_goals_to_group[pow(10,14)]=goals_choice;
      UCB_goals_to_group[pow(10,15)]=goals_choice;
      UCB_goals_to_group[pow(10,16)]=goals_choice;
      UCB_goals_to_group[pow(10,17)]=goals_choice;
      UCB_goals_to_group[pow(10,18)]=goals_choice;
      UCB_goals_to_group[pow(10,19)]=goals_choice;*/



      //need result here to store final PDB collection
      result=make_shared<PatternCollectionInformation>(task, make_shared<PatternCollection>());
      shared_ptr<ModularZeroOnePDBs> candidate_ptr;
      const State &initial_state = task_proxy.get_initial_state();
    
      //cout<<"initial pdb type:"<<pdb_factory->name()<<endl;

      //ModularZeroOnePDBs candidate(task_proxy, Initial_collection.get_PC(), *pdb_factory);
      //best_collection=Initial_collection;
      //generate sample states:
      pattern_evaluator->initialize(task);
      pattern_generator->initialize(task);
      PatternCollectionContainer perimeter_collection=pattern_generator->generate_perimeter();
      double overall_problem_size=perimeter_collection.get_overall_size();
      cout<<"overall_problem_size:"<<overall_problem_size<<endl;
      //Initializing alternative_pattern_generator as well
      alternative_pattern_generator.initialize();
      
      //result included in sample_states call because needed for dead_end detection
      //set_dead_ends add dead_ends for symbolic, NEED TO ASK ALVARO ABOUT THIS
      //OK, so this was a call to collapse de dead_ends so adding new symbolic PDBs
      //Would not take ages or something like that, needs to check on it and see what 
      //we do now.
      //result->set_dead_ends(pdb_factory->get_dead_ends());
      
      ///DISCUSS WITH ALVARO:adding pdbs to current set if evaluator says new collection is helpful
      //WHEN ADDING THE PDB, terminate_creation makes comparisons biased
      //because new candidate has less time to generate pdb, should we wait for 
      //terminate_pdb till the end of subset selection???
      //result->include_additive_pdbs(candidate_ptr->get_pattern_databases());
      float start_time=utils::g_timer();
      PatternCollectionContainer initial_Gamer_Collection=alternative_pattern_generator.get_PC();
      //cout<<"time:"<<utils::g_timer<<"Initial Gamer PDB:";initial_Gamer_Collection.print();
      candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, initial_Gamer_Collection.get_PC(), *pdb_factory);
	  best_h=candidate_ptr->get_value(initial_state);
	  cout<<"initial h value for Gamer-Style:"<<best_h<<endl;
	  pattern_evaluator->sample_states_Gamer(result,best_h/2+best_h%2);
      //best_value=pattern_evaluator->calculate_eval_score(candidate_ptr);
	  //cout<<"initial best_value for Gamer-Style:"<<best_value<<endl;
      float pdb_time=utils::g_timer()-start_time;
	  //UCB_generator.increase_cost(0,pdb_time);
	  //UCB_generator.increase_reward(0,pdb_time);

      //double initial_pdb_size=pow(10.0,UCB_sizes.make_choice());
	  double initial_pdb_size=pow(10.0,9);//start safely!
	  //double initial_disjunctive=UCB_Disjunctive_patterns[initial_pdb_size].make_choice();
	  bool initial_disjunctive=false;
	  double initial_goals_to_group=goals_choice.make_choice();
	  //Always adding first collection
	  start_time=utils::g_timer();
      if(create_perimeter){
        create_perimeter=false;
        cout<<"seeding with creating_perimeter,time:"<<utils::g_timer<<endl;
        cout<<"perimeter PC:";perimeter_collection.print();
        candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, perimeter_collection.get_PC(), *pdb_factory);
        result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases(), 250000, 50000, 10000000));
        //result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases(), 85000, 50000, 10000000));
        if(pdb_factory->is_solved()){
          cout<<"Perimeter solved the problem!"<<endl;
        }
        cout<<"seeding with creating_perimeter finished,time:"<<utils::g_timer<<",h_val:"<<candidate_ptr->get_value(initial_state)<<endl;
        //cout<<"second terminate finished,time:"<<utils::g_timer<<",h_val:"<<candidate_ptr->get_value(initial_state)<<endl;
        //pattern_evaluator->set_threshold(1);//If using perimeter, then we want all heuristics which can see further
      }
      else{
	    pattern_generator->set_pdb_max_size(initial_pdb_size);
        //pattern_generator->set_disjunctive_patterns(initial_disjunctive);
        pattern_generator->set_goals_to_add(initial_goals_to_group);
		cout<<"time:"<<utils::g_timer<<"Initial CBP,with pdb_max_size:"<<initial_pdb_size<<" disjunctive_choice:"<<initial_disjunctive<<" num_goals_to_group:"<<initial_goals_to_group<<endl;
        PatternCollectionContainer Initial_collection=pattern_generator->generate();
        cout<<"Initial PC:";
        Initial_collection.print();
        PatternCollection temp_pc=Initial_collection.get_PC();
        cout<<"temp_pc.size:"<<temp_pc.size()<<endl;
        candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, Initial_collection.get_PC(), *pdb_factory);
		pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases(),50000, 50000, 5000000);
        //cout<<"Initial h value :"<<candidate_ptr->get_value(initial_state)<<endl;
		//cout<<"Initial avg_h:"<<candidate_ptr->compute_approx_mean_finite_h();
        //pdb_factory->no_terminate_creation(candidate_ptr->get_pattern_databases());//To update pdb_factory in case problem is solved
        result->include_additive_pdbs(make_shared<PDBCollection>(candidate_ptr->get_pattern_databases())); 
      }
	  pdb_time=utils::g_timer()-start_time;
        
      if(pdb_factory->is_solved()){
          cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
          return;
      }
      else{
        cout<<"Sadly,solution not found yet after first PDB ;-) so we keep going"<<endl;
      }
      
      initial_h=result->get_value(initial_state);
	  sampling_h=initial_h;
      cout<<"Initial collection zero-one h value:"<<initial_h<<endl;
      
      PC_counter++;
      result->set_dead_ends(pdb_factory->get_dead_ends());//NOT SURE WHAT IT DOES, CHECK WITH ALVARO
       //NOTE:Skipped on purpose this for perimeter, because not want to give option to terminate for UBC
       //I think perimeter PDB would be too big for this to be effective
      pattern_evaluator->sample_states(result);
	  if(!create_perimeter)	{
		  UCB_generator.increase_cost(1,pdb_time);
		  UCB_sizes.increase_cost(9,pdb_time);
		  goals_choice.increase_cost(initial_goals_to_group,pdb_time);
		  UCB_generator.increase_reward(1,pdb_time);
		  UCB_sizes.increase_reward(9,pdb_time);
		  goals_choice.increase_reward(initial_goals_to_group,pdb_time);
	  }
      //cout<<"first pdb_max_size:"<<pattern_generator->get_pdb_max_size()<<endl;
      
      //Testing 1 Bin packing algorithms
      //PatternCollectionGeneratorBinPackingV1 Packer1(opts3);
      PatternCollectionGeneratorBinPackingV1 Packer2(opts4); 
      PatternCollectionGeneratorBinPackingV1 Packer3(opts5);

      //Packer1.initialize(task);
      Packer2.initialize(task);
      Packer3.initialize(task);
      
      //Testing 1st Packer  
      
      /*float temp_start_time=utils::g_timer();
      double temp_pdb_max_size=pow(10.0,7);
      while(utils::g_timer()-temp_start_time<50){
	temp_pdb_max_size=10*temp_pdb_max_size;
	Packer1.set_pdb_max_size(temp_pdb_max_size);
	PatternCollectionContainer Packer1Coll=Packer1.generate();
	shared_ptr<ModularZeroOnePDBs> BinPackerV1_ptr1=make_shared<ModularZeroOnePDBs>(task_proxy, Packer1Coll.get_PC(), *pdb_factory);
	cout<<"time:"<<utils::g_timer()<<"pdb_max_sizer:"<<temp_pdb_max_size<<",testing PC:";Packer1Coll.print();
	if(pattern_evaluator->evaluate(BinPackerV1_ptr1)){
	  cout<<"time:"<<utils::g_timer()<<",Packer1 improves initial PDB, adding to selection"<<endl;
	  result->include_additive_pdbs(pdb_factory->terminate_creation(BinPackerV1_ptr1->get_pattern_databases()));
	  result->set_dead_ends(pdb_factory->get_dead_ends());
	  pattern_evaluator->sample_states(result);
	  cout<<"After Packing1 initial h value:"<<result->get_value(initial_state)<<",max_pdb_size:"<<temp_pdb_max_size<<endl;
	}
      }
      cout<<"Packing1 initial h value:"<<result->get_value(initial_state)<<endl;*/
	  int best_size=11;
      //Testing 2nd Packer  
	  float eval_time=0,eval_start_time;
      float temp_start_time=utils::g_timer();
	  int cur_size=8;
      double temp_pdb_max_size=pow(10.0,8);
      while(utils::g_timer()-temp_start_time-eval_time<80){//time hardcoded to 50 secs because authors expect this method to do worse in average than next bin packing method, future verisons of this should be learnt!
		if(modular_heuristic_timer->is_expired()){
			cout<<"timer_is_expired,break-2"<<endl;
            break;
		}	
		if(overall_problem_size<=temp_pdb_max_size)
			break;//not going past problem size!
	    if(double(utils::get_current_memory_in_kb())/1024.0>memory_limit){
			cout<<"break-4,memory limit breached,current_memory(MB):"<<utils::get_current_memory_in_kb()/1024.0<<",memory_limit:"<<memory_limit<<endl;
			return;//Gone past memory limit, best to return with what we have
	    }
		Packer2.set_pdb_max_size(temp_pdb_max_size);
		PatternCollectionContainer Packer2Coll=Packer2.generate();
		shared_ptr<ModularZeroOnePDBs> BinPackerV1_ptr2=make_shared<ModularZeroOnePDBs>(task_proxy, Packer2Coll.get_PC(), *pdb_factory);
		auto const &candidate_pdb=BinPackerV1_ptr2->get_pattern_databases();
		pdb_factory->terminate_creation(candidate_pdb,50000, 50000, 5000000);
		PC_counter++;
		if(pdb_factory->is_solved()){
		  cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
		  result->include_additive_pdbs(make_shared<PDBCollection>(candidate_pdb));
		  return;
		}
		if(pattern_evaluator->evaluate(BinPackerV1_ptr2)){
		  cout<<"time:"<<utils::g_timer()<<",Packer2 improves initial PDB, adding to selection"<<endl;
		  result->include_additive_pdbs(make_shared<PDBCollection>(candidate_pdb));
		  result->set_dead_ends(pdb_factory->get_dead_ends());
		  new_initial_h=BinPackerV1_ptr2->get_value(initial_state);
		  if(new_initial_h>initial_h) { 
			  initial_h=new_initial_h;
			  cout<<"After Packing2 initial h value:"<<initial_h<<",max_pdb_size:"<<temp_pdb_max_size<<endl;
			   if(initial_h>sampling_h*1.1) {
				  sampling_h=initial_h;
				  eval_start_time=utils::g_timer();
				  pattern_evaluator->sample_states(result);
				  eval_time+=utils::g_timer()-eval_start_time;
			   }
		  }
		  auto &current_max_additive_subsets=*(result->get_max_additive_subsets());
		  if(current_max_additive_subsets.size()>=2){
			  auto const &new_result=make_shared<PatternCollectionInformation>(task, make_shared<PatternCollection>());
			  auto const candidate_pdb=current_max_additive_subsets.back();
			  current_max_additive_subsets.pop_back();
		      pattern_evaluator->clear_dominated_heuristics(current_max_additive_subsets,new_result,candidate_pdb);
			  result=new_result;
		  }
		  best_size=max(best_size,cur_size);
		}
		cur_size++;
		temp_pdb_max_size*=10;
      }
      cout<<"After Packing2 initial h value:"<<initial_h<<" total time:"<<utils::g_timer()-temp_start_time<<endl;
      pdb_factory->clean_up();
	  //Testing 3rd Packer  
	  eval_time=0;
      temp_start_time=utils::g_timer();
      cur_size=8;
      temp_pdb_max_size=pow(10.0,8);
      while(utils::g_timer()-temp_start_time-eval_time<80){//time hardcoded to  75 secs because authors expect this method to do better, future verisons of this should be learnt!
		 if(modular_heuristic_timer->is_expired()){
			cout<<"timer_is_expired,break-2"<<endl;
            break;
		 }	
		 if(overall_problem_size<=temp_pdb_max_size)
			break;//not going past problem size!
		 if(double(utils::get_current_memory_in_kb())/1024.0>memory_limit){
			cout<<"break-4,memory limit breached,current_memory(MB):"<<utils::get_current_memory_in_kb()/1024.0<<",memory_limit:"<<memory_limit<<endl;
			return;//Gone past memory limit, best to return with what we have
		 }
		Packer3.set_pdb_max_size(temp_pdb_max_size);
		PatternCollectionContainer Packer3Coll=Packer3.generate();
		shared_ptr<ModularZeroOnePDBs> BinPackerV1_ptr3=make_shared<ModularZeroOnePDBs>(task_proxy, Packer3Coll.get_PC(), *pdb_factory);
		auto const &candidate_pdb=BinPackerV1_ptr3->get_pattern_databases();
		pdb_factory->terminate_creation(candidate_pdb,50000, 50000, 5000000);
		PC_counter++;
		if(pdb_factory->is_solved()){
		  cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
		  result->include_additive_pdbs(make_shared<PDBCollection>(candidate_pdb));
		  return;
		}
		if(pattern_evaluator->evaluate(BinPackerV1_ptr3)){
		  cout<<"time:"<<utils::g_timer()<<",Packer3 improves initial PDB, adding to selection"<<endl;
		  result->include_additive_pdbs(make_shared<PDBCollection>(candidate_pdb));
		  result->set_dead_ends(pdb_factory->get_dead_ends());
		  new_initial_h=BinPackerV1_ptr3->get_value(initial_state);
		  if(new_initial_h>initial_h) { 
			  initial_h=new_initial_h;
			  cout<<"After Packing3 initial h value:"<<initial_h<<",max_pdb_size:"<<temp_pdb_max_size<<endl;
			  if(initial_h>sampling_h*1.1) {
				  sampling_h=initial_h;
				  eval_start_time=utils::g_timer();
				  pattern_evaluator->sample_states(result);
				  eval_time+=utils::g_timer()-eval_start_time;
				}
		  }
		  auto &current_max_additive_subsets=*(result->get_max_additive_subsets());
		  if(current_max_additive_subsets.size()>=2){
			  auto const &new_result=make_shared<PatternCollectionInformation>(task, make_shared<PatternCollection>());
			  auto const candidate_pdb=current_max_additive_subsets.back();
			  current_max_additive_subsets.pop_back();
		      pattern_evaluator->clear_dominated_heuristics(current_max_additive_subsets,new_result,candidate_pdb);
			  result=new_result;
		  }
		  best_size=max(best_size,cur_size);
		}
		cur_size++;
		temp_pdb_max_size*=10;
      }
	  cout << "After Packing3 initial h value:" <<initial_h<<" total time:"<<utils::g_timer()-temp_start_time<< endl;
      pdb_factory->clean_up();
	  ////////////////////////////////////
	  /*if(utils::g_timer()<400) {
		  PatternCollectionContainer initial_Gamer_Collection=alternative_pattern_generator.get_PC();
		  candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, initial_Gamer_Collection.get_PC(), *pdb_factory);
		  pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases(),50000, 50000, 5000000);
		  result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
		  new_initial_h=candidate_ptr->get_value(initial_state);
		   if(new_initial_h>initial_h) { 
			  initial_h=new_initial_h;
			  pattern_evaluator->sample_states(result);
		   }
	  }*/
	  
	  set<int> sizes={10,11,12,13,14,15,16,17,18,19,20,25,30,35};
	  //for(int i=-2;i<=2;i++) sizes.insert(best_size+i);
	  for (int s:sizes)	                                        //filter out sizes that is too large 
		  if(s<best_size+5)	UCB_sizes.insert_choice(s);	  
	  //UCB_sizes.increase_reward(11);
	  
 
	  //bool check_to_terminate=false;
      //bool terminate_or_not=true;
      PatternCollectionContainer candidate_collection;
      //separate so clearing it for other bin packing algorithms does not clear partial gamer collections which we might want to try
      PatternCollectionContainer candidate_collection_Gamer;
      PatternCollectionContainer selected_collection_Gamer;
	  /*deque<vector<Pattern>> candidate_collection_buffer; //for backtrcking
	  deque<map<size_t,State>> sample_states_buffer;   
	  deque<double> best_value_buffer;*/
      while(!modular_heuristic_timer->is_expired()){
		  if(double(utils::get_current_memory_in_kb())/1024.0>memory_limit){
			 cout<<"time:"<<utils::g_timer()<<",modular_heuristic constructor finished,time:"<<utils::g_timer()<<",episodes:"<<num_episodes<<",PC created:"<<PC_counter<<",final_pdbs:"<<result->get_patterns()->size()<<endl;
			 cout<<"break-4,memory limit breached,current_memory(MB):"<<utils::get_current_memory_in_kb()/1024.0<<",memory_limit:"<<memory_limit<<endl;
			 //Check if any collection is useless because it is dominated  //&&pattern_evaluator->get_num_unique_samples()>=100
			 auto &current_max_additive_subsets=*(result->get_max_additive_subsets());
			 if(current_max_additive_subsets.size()>=2){
			 /*cout<<"max_additive_subsets before pruning:"<<current_max_additive_subsets.size()<<endl;
			 for (size_t num=0; num<current_max_additive_subsets.size();num++){
				cout<<"\tcollection["<<num<<"]:"<<endl;
				const auto &pdbs=current_max_additive_subsets[num];
				for (size_t collection=0; collection<pdbs.size();collection++){
					cout<<"\tpdb["<<collection<<"]"<<*pdbs[collection]<<endl;
				} 
			 }*/
			 const auto& new_result=make_shared<PatternCollectionInformation>(task, make_shared<PatternCollection>());
			 auto const candidate_pdb=current_max_additive_subsets.back();
			 current_max_additive_subsets.pop_back();
			 pattern_evaluator->clear_dominated_heuristics(current_max_additive_subsets,new_result,candidate_pdb);
			 result=new_result;
			 }
		   pdb_factory->clean_up();
		   return;
		  }
        //First we decide whether to try improving existing pdbs or keep generating new ones
        //0 means no terminate, 1 means to try terminating the pdbs
          
       /* if(unterminated_pdbs){
          int initial_terminate_time=utils::g_timer();
          bool terminate_or_not=terminate_choice.make_choice();
          terminate_or_not=false;
          DEBUG_COMP(cout<<"TERMINATE WHILE SEARCHING FOR PATTERNS DISABLED,DEBUGGING GAMER-STYLE ON ITS OWN"<<endl;);
          DEBUG_COMP(cout<<"terminate_or_not:"<<terminate_or_not<<endl;);
          if(terminate_or_not){
            bool success=false;
            int pdb_counter=0;
            int finished_counter=0;
              auto pdb_collection=result->get_pdbs(); 
              for (auto pdb : *pdb_collection){
                pdb_counter++;
                if(pdb->is_finished()){//looking for first unfinished pdb
                  finished_counter++;
                  continue;
                }
                else{
                  cout<<"pdb is unfinished, trying to terminate it now"<<endl;
                }
                double initial_mean_finite_h=pdb->compute_mean_finite_h();
                cout<<"time:"<<utils::g_timer()<<",average h value before Terminate:"<<initial_mean_finite_h<<endl;
                pdb->terminate_creation(20000,10000,2000000,4000);
                cout<<"time:"<<utils::g_timer()<<",after Terminate average h value:"<<pdb->compute_mean_finite_h()<<endl;
                if(pdb->compute_mean_finite_h()>initial_mean_finite_h){//has terminate improved the pdb
                  if(pdb->is_finished()){
                    result->set_dead_ends(pdb_factory->get_dead_ends());//NOT SURE WHAT IT DOES, CHECK WITH ALVARO
                    cout<<"pdb_finished after terminate,time:"<<utils::g_timer()<<",after setting dead_ends "<<endl;
                    success=true;
                  }
                  //terminate_or_not=false;//already done with the terminate action
                  check_to_terminate=true;//need to check if it was the last pdb to terminate
                  success=true;
                  break;
                }
                //Check if we need to update initial_h_value after terminate
                int new_h_value=result->get_value(initial_state);
                if(result->get_value(initial_state)>initial_h){
                  cout<<"new inital_h_value raised to:"<<new_h_value<<endl;
                  initial_h=new_h_value;
                }
              }
              if(success){
                break;
              }
            if(pdb_counter==finished_counter){
              cout<<"Debug me, trying to terminate but pdb_counter:"<<pdb_counter<<"==finished_counter!!!"<<endl;
            }
            else{
              cout<<"pdb_counter:"<<pdb_counter<<",finished:"<<finished_counter<<endl;
            }
            //Always add cost for the choice,no matter the result!
            terminate_choice.increase_cost(terminate_or_not,utils::g_timer()-initial_terminate_time);
            if(success){//terminate helped
              terminate_choice.increase_reward(terminate_or_not,utils::g_timer()-initial_terminate_time);
              cout<<"terminate helped!"<<endl;
              //Now check if there are more unterminated pdbs
              unterminated_pdbs=false;
            for(auto i : pdb_ptr_collection){
              const PDBCollection & pdb_collection=i->get_pattern_databases(); 
              for (auto pdb : pdb_collection){
                if(pdb->is_finished()){
                  unterminated_pdbs=true;
                  check_to_terminate=false;
                  cout<<"pdb is not terminated,UCB to choose whether to terminate"<<endl;
                  break;
                }
                if(unterminated_pdbs)
                  break;
              }
            }
            continue;//Going back to decide whether to do more terminates
          }
        }
        }*/

        num_episodes++;

        /*if(pdb_factory->is_solved()){
          cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
          //best_pdb_collections.push_back(pdb_factory->terminate_creation(candidate.get_pattern_databases()));
          result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
          return;
        }*/
        //Now choosing between RandomSplit and CBP pattern generation
        if(only_CBP){
          generator_choice=1; cout<<"Forcing generator_choice to non-CGamer-style always, for debugging!"<<endl;
        }
        else if(only_gamer){
          generator_choice=0; cout<<"Forcing generator_choice to Gamer-style always, for debugging!"<<endl;
        }
		else {
          generator_choice=UCB_generator.make_choice();
		}
        DEBUG_COMP(cout<<"time:"<<utils::g_timer()<<",generator_choice:"<<generator_choice<<flush<<endl;);
        float start_time=utils::g_timer();
        float pdb_time=0;
        if(generator_choice==0){//Random split 
		   if(candidate_collection_Gamer.empty()){
			partial_gamer_run=false;
			candidate_collection_Gamer=alternative_pattern_generator.generate();
			cout<<"time:"<<utils::g_timer<<",Gamer,starting new_gamer_run with num_patterns:"<<candidate_collection_Gamer.get_size()<<endl;
			DEBUG_COMP(cout<<"modular_heuristic, new set of candidate_patterns from Gamer:";candidate_collection_Gamer.print();cout<<endl;);
		  }
		  else{
			partial_gamer_run=true;
			cout<<"time:"<<utils::g_timer<<",Gamer,starting partial_gamer_run with num_patterns:"<<candidate_collection_Gamer.get_size()<<endl;
		  }
          DEBUG_COMP(cout<<"modular_heuristic, continuing with remaining set of candidate_patterns from Gamer:";candidate_collection_Gamer.print();cout<<endl;);
        }
        else{
		  candidate_collection.clear();
          //Now getting next pdb size
          //cout<<"previos pdb_max_size:"<<pdb_max_size;
          pdb_max_size=pow(10.0,UCB_sizes.make_choice());
          pattern_generator->set_pdb_max_size(pdb_max_size);
          DEBUG_COMP(cout<<",new pdb_max_size:,"<<pattern_generator->get_pdb_max_size(););
          //bool disjunctive_choice=UCB_Disjunctive_patterns[pdb_max_size].make_choice();
          //disjunctive_choice=binary_choice.make_choice();
          DEBUG_COMP(cout<<",new disjunctive_choice:"<<disjunctive_choice;);
          //pattern_generator->set_disjunctive_patterns(disjunctive_choice);
          /*if(disjunctive_choice){//1 goal per pattern
            //Note that we may end up with more than one goal per disjunctive patern,
            //All this does is to disable greedy goal grouping, which makes sense to try when
            //variables are reusable among patterns but not if grouping all the goal variables 
            //in one patern means discarding most of variables not in first pattern due to inability
            //to reuse.  However, if goals happen to be dividide between a few patterns, that is not a 
            //problem either, can happen if goal variables get chosen randomly into patterns with existing goals
            num_goals_to_group=1;
          }
          else{
            num_goals_to_group=goals_choice.make_choice();
          }*/
		  num_goals_to_group=goals_choice.make_choice();
          DEBUG_COMP(cout<<",new goals_to_group:"<<num_goals_to_group;);
          pattern_generator->set_goals_to_add(num_goals_to_group);

          //Now generating next set of patterns
		  cout<<"time:"<<utils::g_timer<<",CBP,with pdb_max_size:"<<pdb_max_size<<" disjunctive_choice:"<<disjunctive_choice<<" num_goals_to_group:"<<num_goals_to_group<<endl;
          candidate_collection=pattern_generator->generate();
          candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, candidate_collection.get_PC(), *pdb_factory);
		  pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases(), 50000, 50000, 5000000);
          pdb_time=utils::g_timer()-start_time;
          DEBUG_COMP(cout<<"pdb_max_size:"<<pdb_max_size<<",pdb_time:"<<pdb_time<<endl;);
        }
          
		if(generator_choice==0){//check if avg_h_value was improved
            double new_best_value=-1;
            //Raising best_value by a thousandth so that it knows when 
            //DEBUG_COMP(cout<<"best_value raised by a thousandth to:"<<best_value<<",to ensure next set of variables added actually raises the avg_heuristic_value"<<endl;);
            double candidate_value=0;
			if(partial_gamer_run)	candidate_collection_Gamer.shuffle();
			else candidate_collection_Gamer.reverse();//So we can look at the same order as Gamer for comparison purposes
			vector<Pattern> patterns=candidate_collection_Gamer.get_PC();
			/*for(auto const &imprv : improving_patterns){
			  best_value=max(best_value,imprv.first);
			}*/
			improving_patterns.clear();
			cout<<"time:"<<utils::g_timer()<<",Gamer_iteration,partial_gamer_run="<<partial_gamer_run<<",starting best_value:"<<best_value<<endl;
			float start_gamer_iteration_time=utils::g_timer();
            while(patterns.size()>0){
			  cout<<"time:"<<utils::g_timer()<<",Gamer, patterns size:"<<patterns.size()<<endl;
			  Pattern pattern=patterns.back();patterns.pop_back();
			  cout<<"Gamer,working on pattern:"<<pattern<<endl;
			  candidate_collection_Gamer.pop_back();
              PatternCollectionContainer GamerCandidateContainer;GamerCandidateContainer.add_pc(pattern);
                
              //float indiv_start_time=utils::g_timer();
              shared_ptr<ModularZeroOnePDBs> Gamer_candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, GamerCandidateContainer.get_PC(),*pdb_factory);
              //Need to terminate creation to properly evaluate avg_h_value, limit is 30 secs and one million nodes, if bigger it is "deemed" too expensive to consider
              //Just for selection purposes for me
              pdb_factory->terminate_creation(Gamer_candidate_ptr->get_pattern_databases(),60000,20000,10000000);
              //float indiv_pdb_time=utils::g_timer()-indiv_start_time;
              //auto const &pdb=Gamer_candidate_ptr->get_pattern_databases()[0];
              //DEBUG_COMP(cout<<"\tGamer-style, indiv pdb gen time:"<<indiv_pdb_time<<",finished after terminate:"<<pdb->is_finished()<<",for pattern:"<<pattern<<endl;);
              candidate_value=pattern_evaluator->calculate_eval_score(Gamer_candidate_ptr);
              //DEBUG_COMP(cout<<"time:"<<utils::g_timer()<<",evaluating pattern:"<<pattern<<endl;);
              if(candidate_value>best_value){
                cout<<"\ttime:"<<utils::g_timer()<<",Gamer-Style,pattern:";for(auto const & i :pattern) cout<<i<<",";cout<<"raised value to:"<<candidate_value<<" vs prevous iteration best:"<<best_value<<endl;
                new_best_value=max(candidate_value,new_best_value);
                improving_patterns.push_back(make_pair(candidate_value,pattern));
                DEBUG_COMP(cout<<"modular_heuristic,Gamer candidate avg_h_val:"<<best_value<<",initial h:"<<Gamer_candidate_ptr->get_value(initial_state)<<endl;);
              }
              else{
                DEBUG_COMP(cout<<"candidate_pattern:";for(auto const & i :pattern) cout<<i<<",";cout<<",avg_h:"<<candidate_value<<",is smaller than current_avg_h:"<<best_value<<endl;);
              }

              /*if(indiv_pdb_time>150){//so if the time to generate all options is large, bettern to stop now with Gamer candidates
                cout<<"Stopping gamer-style selection, indiv pdb gen time:,"<<indiv_pdb_time<<",limit:150"<<",pattern:"<<pattern<<endl;
                UCB_generator.increase_cost(generator_choice,INT_MAX);//And make sure we do not try it again!
                break;
              }*/
	          if(utils::g_timer()-start_gamer_iteration_time>120){
				cout<<"break-3,interrupting gamer candidate generation, taking too long, invested time:"<<utils::g_timer()-start_gamer_iteration_time<<",limit:120"<<endl;
				break;
			  }
              if(float(utils::g_timer())/float(modular_time_limit)>0.80&&new_best_value>best_value){
		      //PDB generation in Gamer selection can take too long
                //We want to leave some time to add any improving patterns if found
                cout<<"modular_time_limit:"<<modular_time_limit<<",time:"<<utils::g_timer()<<",break-1,no more Gamer-style candidate PDB generation (or any other) because we are too close to time limit and we want remaining time to add any improving patterns from last iteration"<<endl;
                break;
              }
              if(modular_heuristic_timer->is_expired()){//PDB generation in Gamer selection can take too long
				cout<<"timer_is_expired,break-2"<<endl;
                break;
              }
            }
            cout<<"Gamer,removing all added patterns smaller than:"<<0.999*new_best_value<<",improving_patterns_size:"<<improving_patterns.size()<<endl;
              
	    /*for(auto const & imprv : improving_patterns){
	      new_best_value=max(new_best_value,imprv.first);
	    }*/
	
            improving_patterns.erase(std::remove_if(
                  improving_patterns.begin(), 
                  improving_patterns.end(),
                  [new_best_value](pair<double, Pattern> & x){
                  return x.first < 0.999*new_best_value;
                  }), improving_patterns.end());
            set<int> Gamer_Combined;
            for (auto const & i : improving_patterns){
              cout<<"\t\tpattern:";
              for(int j : i.second){ 
                cout<<j<<",";
                Gamer_Combined.insert(j);
              }
              cout<<"is improving avg_h:"<<i.first<<endl;
            }
            vector<int> Gamer_Combined_vect;for (auto const & i : Gamer_Combined) Gamer_Combined_vect.push_back(i);
            if(Gamer_Combined_vect.empty()) {
				if(candidate_collection_Gamer.empty()){
					/*if(!candidate_collection_buffer.empty()) {   //backtrack remaining candidates
						candidate_collection_Gamer.set_PC(candidate_collection_buffer.back());
						pattern_evaluator->set_unique_samples_Gamer(sample_states_buffer.back());
						best_value=best_value_buffer.back();
						candidate_collection_buffer.pop_back();
						sample_states_buffer.pop_back();
						best_value_buffer.pop_back();
						backtrack=true;
					}*/
				  	if(only_gamer) {   //||terminate
						cout<<"Finished, we stop Gamer-style as pattern generator"<<endl;
						break;
					}
					else {
						cout<<"Gamer_Combined_vect:"<<Gamer_Combined_vect<<",is empty,so we are finished with Gamer-style for this problem"<<endl;
						//finished=true;
						UCB_generator.increase_cost(generator_choice,INT_MAX);
					}
				}
			  //else   backtrack=false;
			  pdb_time=utils::g_timer()-start_time;
			  UCB_generator.increase_cost(generator_choice,pdb_time);  //make new choice
              continue;
            }
		  //backtrack=false;
		  cout<<"time:"<<utils::g_timer();
	      //if(partial_gamer_run)    cout<<"Gamer,confirming pattern:"<<Gamer_Combined<<"as improvements were found in a partial run"<<endl; //partial_gamer_improvement_found=true;	                   
	      if(candidate_collection_Gamer.empty())	cout<<",Gamer,confirming pattern:"<<Gamer_Combined<<"as there are no more patterns to test for this iteration"<<endl;
		  else {  
			  cout<<",Gamer,confirming pattern:"<<Gamer_Combined<<"discarding remaining candidates"<<endl;
			  /*if(candidate_collection_buffer.size()>=3) { //store remaining candidates for backtracking, max size 3
				  candidate_collection_buffer.pop_front();
				  sample_states_buffer.pop_front();
				  best_value_buffer.pop_front();
				}
			  candidate_collection_buffer.push_back(candidate_collection_Gamer.get_PC());	
			  sample_states_buffer.push_back(pattern_evaluator->get_unique_samples_Gamer());
			  best_value_buffer.push_back(best_value);*/
			  candidate_collection_Gamer.clear();
			}
	      alternative_pattern_generator.confirm(Gamer_Combined);
          selected_collection_Gamer.restart_pc(Gamer_Combined_vect);
          candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, selected_collection_Gamer.get_PC(), *pdb_factory);
	    //To get decent h values when combining all the variables need to terminate the combined vector as well
	      //cout<<"time:,"<<utils::g_timer()<<",finished:"<<pdb_factory->is_finished()<<endl;
	      auto const &pdb=candidate_ptr->get_pattern_databases()[0];
		  new_best_h=candidate_ptr->get_value(initial_state);
		  if (new_best_h>best_h) {
			  cout<<"best_h has been raised from "<<best_h<<" to "<<new_best_h<<endl;
			  if (new_best_h>best_h*1.1) {
				  best_h=new_best_h;
				  pattern_evaluator->sample_states_Gamer(result,best_h/2+best_h%2);
			  }
		  }
		  best_value=pattern_evaluator->calculate_eval_score(candidate_ptr);
	      pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases(), 50000, 50000, 5000000);
          pdb_time=utils::g_timer()-start_time;//So includes all the pdb generations for this Gamer step
          cout<<"time:"<<utils::g_timer()<<",Gamer_Combined_vect:"<<Gamer_Combined_vect<<",best_h:"<<best_h<<",best_value:"<<best_value<<",finished after terminate:"<<pdb->is_finished()<<endl;
          }
          else{
            UCB_sizes.increase_cost(log10(pdb_max_size),pdb_time);
            //binary_choice.increase_cost(double(disjunctive_choice),pdb_time);
            goals_choice.increase_cost(double(num_goals_to_group),pdb_time);
          }
          UCB_generator.increase_cost(generator_choice,pdb_time);
          //UCB_Disjunctive_patterns[pdb_max_size].increase_cost(double(pattern_generator->get_disjunctive_patterns()),pdb_time);
          //UCB_goals_to_group[pdb_max_size].increase_cost(double(pattern_generator->get_goals_to_add()),pdb_time);
         
 		  PC_counter++;
		  if(pdb_factory->is_solved()){
			  cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
			  result->include_additive_pdbs(make_shared<PDBCollection>(candidate_ptr->get_pattern_databases()));
			  return;
		  }	
		   if(modular_heuristic_timer->is_expired()){
			cout<<"timer_is_expired,break-2"<<endl;
            break;
		  }	
          new_initial_h=candidate_ptr->get_value(initial_state);
		  //initial_h=result->get_value(initial_state);
	  //cout<<"Initial value before evaluations:,"<<initial_h<<",new_initial_h:"<<new_initial_h<<endl;
          
         /* if(generator_choice==0){
            if(new_initial_h>initial_h){
              cout<<"Gamer,increased initial_h from:,"<<initial_h<<",to,"<<new_initial_h<<endl;
            }
          }*/
          

          if(new_initial_h>initial_h){//we always add the collection and re-sample if initial_h increased
            //int temp_initial_h=candidate_ptr->get_value(initial_state);
            //cout<<"time:"<<utils::g_timer()<<",Initial h value before terminate:"<<temp_initial_h<<endl;
            //if(generator_choice==0) cout<<"Gamer,increased initial_h from:,"<<initial_h<<",to,"<<new_initial_h<<endl;
			//else cout<<"CBP,increased initial_h from:,"<<initial_h<<",to,"<<new_initial_h<<endl;
            cout<<"time:"<<utils::g_timer()<<",pdb_max_size:"<<pdb_max_size<<",generator_choice:"<<generator_choice<<",goals_choice:"<<num_goals_to_group<<",disjoint:"<<disjunctive_choice<<",selecting PC because initial_h has been raised from "<<initial_h<<"to "<<new_initial_h<<endl;
			//Clean dominated PDBs after addition of improving patterns
			/*if(utils::g_timer()>=time_to_prun) {
				auto const &candidate_pdb=candidate_ptr->get_pattern_databases();
				pdb_factory->terminate_creation(candidate_pdb);
				auto const &current_max_additive_subsets=*(result->get_max_additive_subsets());
				const auto& new_result=make_shared<PatternCollectionInformation>(task, make_shared<PatternCollection>());
				pattern_evaluator->clear_dominated_heuristics(current_max_additive_subsets,new_result,candidate_pdb);
				result=new_result;
				time_to_prun+=200;
			}*/

		  result->include_additive_pdbs(make_shared<PDBCollection>(candidate_ptr->get_pattern_databases()));	
		  if(modular_heuristic_timer->is_expired()){
			cout<<"timer_is_expired,break-2"<<endl;
            break;
		  }		
	    //Clean dominated PDBs after addition of improving patterns
	    /*if(recompute_additive_sets){
	      cout<<"time:"<<utils::g_timer<<",calling recompute_max_additive_subsets"<<endl;
	      result->recompute_max_additive_subsets();
	      cout<<"time:"<<utils::g_timer<<",after recompute_max_additive_subsets"<<endl;
	    }*/
            
	    //check_to_terminate=true;
            //pdb_ptr_collection.push_back(candidate_ptr);

            //UCB_Disjunctive_patterns[pdb_max_size].increase_reward(double(pattern_generator->get_disjunctive_patterns()),pdb_time);
            //UCB_goals_to_group[pdb_max_size].increase_reward(double(pattern_generator->get_goals_to_add()),pdb_time);
            UCB_generator.increase_reward(generator_choice,pdb_time);
            if(generator_choice!=0){//Not Gamer options
              UCB_sizes.increase_reward(log10(pdb_max_size),pdb_time);
              //binary_choice.increase_reward(double(disjunctive_choice),pdb_time);              
			  goals_choice.increase_reward(double(num_goals_to_group),pdb_time);
			  //if(finished) CBP_extra=0;
			  /*if(failure) {
				  failure=false;
				  failures_in_a_row=0;
			  }*/
            }
            result->set_dead_ends(pdb_factory->get_dead_ends());//NOT SURE WHAT IT DOES, CHECK WITH ALVARO
            initial_h=candidate_ptr->get_value(initial_state);//Might be higher than simply new cadidate collection value due to max additive pattern combinations beyond current collection, unlikely but possible!
			if (initial_h == numeric_limits<int>::max()) {
				cout<<"initial state is dead_end according to PDB, problem unsolvable!!!"<<endl;
				exit(1);
             }
			if(initial_h>sampling_h*1.1) {
				  sampling_h=initial_h;
				  pattern_evaluator->sample_states(result);
			}
          }
				//OK,so lets check if candidate_PC is good enough to add to current collection
          else if(pattern_evaluator->evaluate(candidate_ptr)){
              //NEED TO CHECK WITHOUT TERMINATING CREATION UNTIL ALL PDBs ARE SELECTED
              cout<<"time:"<<utils::g_timer()<<",pdb_max_size:"<<pdb_max_size<<",generator_choice:"<<generator_choice<<",disjoint:"<<disjunctive_choice<<",goals_choice:"<<num_goals_to_group<<",modular_heuristic_selecting PC"<<endl;
			//Clean dominated PDBs after addition of improving patterns
              /*if(utils::g_timer()>=time_to_prun) {
				auto const &candidate_pdb=candidate_ptr->get_pattern_databases();
				pdb_factory->terminate_creation(candidate_pdb);
				auto const &current_max_additive_subsets=*(result->get_max_additive_subsets());
				const auto& new_result=make_shared<PatternCollectionInformation>(task, make_shared<PatternCollection>());
				pattern_evaluator->clear_dominated_heuristics(current_max_additive_subsets,new_result,candidate_pdb);
				result=new_result;
				time_to_prun+=200;
			}*/

			  result->include_additive_pdbs(make_shared<PDBCollection>(candidate_ptr->get_pattern_databases())); 
			  if(modular_heuristic_timer->is_expired()){
				cout<<"timer_is_expired,break-2"<<endl;
				break;
			  }	
	    /*if(recompute_additive_sets){
	      cout<<"time:"<<utils::g_timer<<",calling recompute_max_additive_subsets"<<endl;
	      result->recompute_max_additive_subsets();
	      cout<<"time:"<<utils::g_timer<<",after recompute_max_additive_subsets"<<endl;
	    }*/
              //check_to_terminate=true;
              UCB_generator.increase_reward(generator_choice,pdb_time);
              if(generator_choice!=0){
                UCB_sizes.increase_reward(log10(pdb_max_size),pdb_time);
                //UCB_Disjunctive_patterns[pdb_max_size].increase_reward(double(pattern_generator->get_disjunctive_patterns()),pdb_time);
                //UCB_goals_to_group[pdb_max_size].increase_reward(double(pattern_generator->get_goals_to_add()),pdb_time);
                //binary_choice.increase_reward(double(disjunctive_choice),pdb_time);
                goals_choice.increase_reward(double(num_goals_to_group),pdb_time);
				//if(finished) CBP_extra=0;
				/*if(failure) {
				  failure=false;
				  failures_in_a_row=0;
			  }*/
              }
              result->set_dead_ends(pdb_factory->get_dead_ends());//NOT SURE WHAT IT DOES, CHECK WITH ALVARO
			  new_initial_h=candidate_ptr->get_value(initial_state);
			  if(new_initial_h>initial_h) {
				  initial_h=new_initial_h;
				  if (initial_h == numeric_limits<int>::max()) {
					cout<<"initial state is dead_end according to PDB, problem unsolvable!!!"<<endl;
					exit(1);
				  }
				  //should we always re-sample or only if number of improved states is large enough?
				  if(initial_h>sampling_h*1.1) {
					  sampling_h=initial_h;
					  pattern_evaluator->sample_states(result);
				  }
			  }
            }
            else{
			  /*if(generator_choice!=0) {
				//  if(failure) {
				//	  failures_in_a_row++;
				//	  if(failures_in_a_row>=12) terminate=true;
				//  }
				//  else failure=true;
				  if(finished) {
					  CBP_extra++; 
					  if(CBP_extra>=5) terminate=true;
				  }
				  if(terminate) {
				  if(finished||only_CBP){
					cout<<"Finished, we stop CBP as pattern generator"<<endl;
					break;
					}
				  else {
					cout<<"We are finished with CBP for this problem"<<endl;
					UCB_generator.increase_cost(generator_choice,INT_MAX);
					}
				}
				}*/
              DEBUG_COMP(cout<<"time:"<<utils::g_timer()<<"pdb_max_size:"<<pdb_max_size<<",modular_heuristic_not_selecting_PC"<<endl;);
            }
       
          //Checking if any of the pdbs were not terminated
          //If they are, allow UCB to choose whether to terminate them instead of generating new pdbs
          /*if(check_to_terminate&&unterminated_pdbs==false){
            check_to_terminate=false;
            unterminated_pdbs=false;
              auto pdb_collection=result->get_pdbs(); 
              for (auto pdb : *pdb_collection){
                if(!pdb->is_finished()){
                  unterminated_pdbs=true;
                  check_to_terminate=false;
                  cout<<"pdb is not terminated,UCB to choose whether to terminate"<<endl;
                  break;
                }
              }
          }*/
      }
	  pdb_factory->clean_up();
      //Now terminate creation of unfinished selected PDBs
      start_time=utils::g_timer();
	  auto &current_max_additive_subsets=*(result->get_max_additive_subsets());
      for(auto const & pdbs : boost::adaptors::reverse(current_max_additive_subsets)){//better to start terminating the best PDBs
        //usually the last added, in case we run out of time
          if (utils::g_timer()>1320) {
		    cout<<"time:"<<utils::g_timer()<<",terminate_time:"<<utils::g_timer()-start_time<<",interrupting terminate because time is up"<<endl;
            break;
		  }
		  if(double(utils::get_current_memory_in_kb())/1024.0>memory_limit){
			cout<<"break-4,memory limit breached,current_memory(MB):"<<utils::get_current_memory_in_kb()/1024.0<<",memory_limit:"<<memory_limit<<endl;
			break;
		  }
          int temp_h=pattern_evaluator->calculate_max_additive_subset(pdbs,initial_state);
          cout<<"time:"<<utils::g_timer()<<",initial pdb collection h value before terminate:"<<temp_h<<endl;
		  for (size_t collection=0; collection<pdbs.size();collection++){
			cout<<"\tpdb["<<collection<<"]"<<*pdbs[collection]<<endl;
			} 
		  pdb_factory->terminate_creation(pdbs, 100000, 100000, 20000000);
          temp_h=pattern_evaluator->calculate_max_additive_subset(pdbs,initial_state);
          cout<<"time:"<<utils::g_timer()<<",initial pdb collection h value after terminate:"<<temp_h<<endl;
		  if(pdb_factory->is_solved()){
		    cout<<"Solution found while creating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
		    return;
			}
        /*if((utils::g_timer()-start_time)>200.0){
          cout<<"time:"<<utils::g_timer()<<",terminate_time:"<<utils::g_timer()-start_time<<",Interrupting terminate because it is taking too long"<<endl;
          break;
        }*/
      }
      //	    
      cout<<"time:"<<utils::g_timer()<<",modular_heuristic constructor finished,time:"<<utils::g_timer()<<",episodes:"<<num_episodes<<",PC created:"<<PC_counter<<",final_pdbs:"<<result->get_patterns()->size()<<",terminate_time:"<<utils::g_timer()-start_time<<endl;
      //result->recompute_max_additive_subsets();
      //cout<<"time:"<<utils::g_timer()<<",after recompute_max_additive_subset,Testing modular_heuristic constructor finished,time:"<<utils::g_timer()<<",episodes:"<<num_episodes<<",PC created:"<<PC_counter<<",final_pdbs:"<<result->get_patterns()->size()<<",terminate_time:"<<terminate_time<<endl;
	  
	  //Check if any collection is useless because it is dominated //&&pattern_evaluator->get_num_unique_samples()>=100
	  if(current_max_additive_subsets.size()>=2){
		auto const& new_result=make_shared<PatternCollectionInformation>(task, make_shared<PatternCollection>());
		auto const candidate_pdb=current_max_additive_subsets.back();
		current_max_additive_subsets.pop_back();
		pattern_evaluator->clear_dominated_heuristics(current_max_additive_subsets,new_result,candidate_pdb);
		result=new_result;
	  }
	  pdb_factory->clean_up();     //clean up stored pdbs to release memory
    }

int ModularHeuristic::compute_heuristic(const GlobalState &global_state) {
    State state = convert_global_state(global_state);
    int h = result->get_value(state);
    if (h == numeric_limits<int>::max()) {
        return DEAD_END;
    } else {
        return h;
    }
}

//int ModularHeuristic::compute_heuristic(const State &state) const {
//    int h = pdb.get_value(state);
//    if (h == numeric_limits<int>::max())
//        return DEAD_END;
//    return h;
//}

//int ModularHeuristic::compute_heuristic_id(size_t state_id) {
//  //cout<<"calling offline_compute_heuristic_id"<<endl;fflush(stdout);
//  //cout<<"state_id="<<state_id<<",entries:"<<num_states<<endl;fflush(stdout);
//    int h = pdb.distances[state_id];
//    //cout<<"h_offline:"<<h<<endl;fflush(stdout);
//    if (h == numeric_limits<int>::max())
//        return INT_MAX/2;//Better when doing maxes
//    return h;
//}

static Heuristic *_parse(OptionParser &parser) {
    parser.document_synopsis("Pattern database heuristic", "TODO");
    parser.document_language_support("action costs", "supported");
    parser.document_language_support("conditional effects", "not supported");
    parser.document_language_support("axioms", "not supported");
    parser.document_property("admissible", "yes");
    parser.document_property("consistent", "yes");
    parser.document_property("safe", "yes");
    parser.document_property("preferred operators", "no");


    parser.add_option<shared_ptr<PatternCollectionGeneratorComplementary>>(
        "patterns",
        "pattern Collection generation method",
        "modular_rbp");
    parser.add_option<shared_ptr<PatternCollectionEvaluator>>(
        "evaluator",
        "pattern Collection evaluation method",
        "rand_walk_evaluator_avg_h");
    parser.add_option<int>(
        "modular_time_limit",
        "time limit in seconds for modular_pdb_heuristic initialization cut off",
        "900");
    parser.add_option<shared_ptr<PDBFactory>>(
        "pdb_factory",
        "See detailed documentation for pdb factories. ",
	      "modular_symbolic");
    parser.add_option<bool>(
        "terminate_creation",
        "give extra generation time to selected PDBs but not to candidate PDBs",
	      "false");
    parser.add_option<bool>(
        "create_perimeter",
        "do initial perimeter",
	      "false");
    parser.add_option<bool>(
        "only_gamer",
        "only gamer-style w/wo Perimeter",
	      "false");
    parser.add_option<bool>(
        "only_cbp",
        "only CBP-style w/wo Perimeter",
	      "false");
    
    Heuristic::add_options_to_parser(parser);

    Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;

    return new ModularHeuristic(opts);
}

static Plugin<Heuristic> _plugin("modular_pdb", _parse);
}
