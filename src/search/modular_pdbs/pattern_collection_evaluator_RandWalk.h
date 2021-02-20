#ifndef MODULAR_PDBS_PATTERN_COLLECTION_EVALUATOR_RANDWALK_H
#define MODULAR_PDBS_PATTERN_COLLECTION_EVALUATOR_RANDWALK_H

#include "pattern_collection_evaluator.h"
#include "types.h"

#include <memory>
#include <vector>
#include <map>
#include <random>
#include "../global_state.h"
#include "../task_proxy.h"
#include "../successor_generator.h"
#include "pattern_collection_information.h"


class AbstractTask;

namespace options {
class OptionParser;
class Options;
}
namespace utils {
class CountdownTimer;
}

namespace pdbs3 {
//class PDBFactory;
//class PatternCollectionContainer;
class PatterCollectionEvaluatorRandWalk : public PatternCollectionEvaluator {
	int time_limit=20;
  unsigned increased_states=0;//We keep this value so we can pass overall progress 
  std::vector<State> samples;
  //std::shared_ptr<AbstractTask> task;
  std::shared_ptr<TaskProxy> task_proxy;
  std::unique_ptr<SuccessorGenerator> successor_generator;
  std::shared_ptr<PatternCollectionInformation> result;
  utils::CountdownTimer *evaluator_timer;
  std::map<size_t,std::pair<State,int> > unique_samples;
  
  public:
  
  virtual void initialize(std::shared_ptr<AbstractTask> task) override;
	explicit PatterCollectionEvaluatorRandWalk(const options::Options &options);
  virtual bool evaluate(std::shared_ptr<ModularZeroOnePDBs> candidate_PC) override;
  virtual void sample_states(std::shared_ptr<PatternCollectionInformation> current_result) override;
  virtual std::map<size_t,State> get_unique_samples_Gamer() override {return std::map<size_t,State>();}
  virtual void set_unique_samples_Gamer(std::map<size_t,State>& /*input_samples*/) override {}
  virtual void sample_states_Gamer(std::shared_ptr<PatternCollectionInformation> /*current_result*/,int /*init_h*/) override {}
  virtual double calculate_eval_score(std::shared_ptr<ModularZeroOnePDBs>/*candidate_PC*/) override {return 0;}
  virtual int get_reward() override {return increased_states;}//So metric is improved_states
  //virtual void set_reward(const PatternCollectionContainer & pc, double reward) = 0;
  virtual int get_time_limit() override {return time_limit;} 
  virtual int get_num_unique_samples() override {return unique_samples.size();} 
  virtual void clear_dominated_heuristics(const MaxAdditivePDBSubsets &/*current_max_additive_subsets*/,const std::shared_ptr<PatternCollectionInformation> &/*new_result*/,const PDBCollection &/*candidate_pdb*/) override {}
};

}

#endif
