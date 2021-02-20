#include "pattern_collection_evaluator.h"
#include "pattern_database_interface.h"

#include "../option_parser.h"
#include "../plugin.h"
//#include "../task_proxy.h"

#include "../task_tools.h"
#include "../utils/timer.h"

    
using namespace std;
//PatterCollectionEvaluator is to be the driver for 8 PDB-based options
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
//    virtual void initialize(std::shared_ptr<AbstractTask> task) {
//    cout << "Manual pattern collection: " << *patterns << endl;
//    return PatternCollectionInformation(task, patterns);
//
//   } 
    
int PatternCollectionEvaluator::calculate_max_additive_subset(const PDBCollection &max_subset,const State &current_state){
	int h=0;
	int h_temp=0;
	for (auto const &pdb : max_subset){
	h_temp=pdb->get_value(current_state);
	if (h_temp == numeric_limits<int>::max()){
		h=numeric_limits<int>::max();
		break;
	}
	else{
		h+=h_temp;
	}
	}
	return h;
}

static options::PluginTypePlugin<PatternCollectionEvaluator> _type_plugin(
"PatterCollectionEvaluator",
"The various pattern evaluation algorithms usable by the complementary_modular_pdbs heurisitic.");
}
