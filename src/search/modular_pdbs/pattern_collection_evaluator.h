#ifndef MODULAR_PDBS_PATTERN_COLLECTION_EVALUATOR_H
#define MODULAR_PDBS_PATTERN_COLLECTION_EVALUATOR_H

//#include "types.h"
#include "pattern_collection_generator_complementary.h"
#include "pattern_collection_information.h"
#include "zero_one_pdbs.h"
#include <memory>
#include <vector>
#include <random>


class AbstractTask;

namespace options {
class Options;
}
namespace utils {
class CountdownTimer;
}

namespace pdbs3 {
//class PDBFactory;
//class PatternCollectionContainer;
class PatternCollectionEvaluator {
    //std::shared_ptr<PatternCollection> patterns;
  unsigned threshold=100;
  unsigned num_samples=10000;
  protected:
  double sample_score=0;
  double eval_score=0;
  public:
    virtual void initialize(std::shared_ptr<AbstractTask> task) = 0;
    virtual bool evaluate(std::shared_ptr<ModularZeroOnePDBs> candidate_PC)=0;
	virtual double calculate_eval_score(std::shared_ptr<ModularZeroOnePDBs> candidate_PC)=0;
    virtual void sample_states(std::shared_ptr<PatternCollectionInformation> current_result)=0;
	virtual void sample_states_Gamer(std::shared_ptr<PatternCollectionInformation> current_result,int init_h)=0;
	virtual std::map<size_t,State> get_unique_samples_Gamer()=0;
    virtual void set_unique_samples_Gamer(std::map<size_t,State>& input_samples)=0;
    void set_threshold(const unsigned thres){
      threshold=thres;
    }
    unsigned get_threshold(){
      return threshold;
    }
    void set_num_samples(const unsigned samples){
      num_samples=samples;
    }
    unsigned get_num_samples(){
      return num_samples;
    }
	virtual int get_num_unique_samples()=0;
    virtual int get_reward() = 0;//How much better is than current heuristic in whichever metric we are using
    virtual int get_time_limit()=0;
    //virtual void set_reward(const PatternCollectionContainer & pc, double reward) = 0;
    void set_sample_score(double _score){
      sample_score=_score;
    }
    double get_sample_score(){
      return sample_score;
    }
    void set_eval_score(double _score){
      eval_score=_score;
    }
    double get_eval_score(){
      return eval_score;
    }
	int calculate_max_additive_subset(const PDBCollection &max_subset,const State &current_state);
	virtual void clear_dominated_heuristics(const MaxAdditivePDBSubsets &current_max_additive_subsets,const std::shared_ptr<PatternCollectionInformation> &new_result,const PDBCollection &candidate_pdb)=0; 
};

}

#endif
