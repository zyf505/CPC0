#include "pattern_collection_generator_bin_packing_v1.h"

//#include "../causal_graph.h"
//#include "../globals.h"
#include "../option_parser.h"
#include "../plugin.h"
#include "../bin_packing/first_fit_decreasing_bin_packer.h"
#include "../bin_packing/first_fit_decreasing_related_bin_packer.h"
#include "../bin_packing/first_fit_increasing_related_bin_packer.h"

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
//#include "../successor_generator.h"
//#include "../utils/countdown_timer.h"
//#include "pdb_factory.h"
//#include "pattern_database_interface.h"
//#include "../utils/debug_macros.h"
//#include <random>
//#include "../sampling.h"
#include "../task_tools.h"
#include "../options/plugin.h"


using namespace std;
using namespace pdbs;
//PatternCollectionGeneratorBinPackingV1 is to be the driver for 8 PDB-based options
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
//     H<-Re-evaluate (H�generatedPCs)
//     Learning(BinPackingRewards)
//ComPDBs (Hl)
namespace pdbs3 {

    void PatternCollectionGeneratorBinPackingV1::initialize(std::shared_ptr<AbstractTask> task) {
        task_proxy = make_shared<TaskProxy>(*task);
    }

    PatternCollectionGeneratorBinPackingV1::PatternCollectionGeneratorBinPackingV1(const options::Options & opts) :
    packer_selection(opts.get<int>("packer_selection")) {

        cout << "Packer_selection " << size_t(packer_selection) << endl;
        cout << "Packers " << packers.size() << endl;

        packers.resize(3);
        packers[0] = make_shared<FirstFitDecreasingBinPacker>(max_single_PDB_size,1 );
        packers[1] = make_shared<FirstFitDecreasingRelatedBinPacker>(max_single_PDB_size, 1);
        packers[2] = make_shared<FirstFitIncreasingRelatedBinPacker>(max_single_PDB_size, 1);

        if (size_t(packer_selection)>packers.size()) {
	        cout<<"selection can not be bigger than number of available packers!, pls DEBUG ME!"<<endl;
	        exit(1);
        }
    }

    PatternCollectionContainer PatternCollectionGeneratorBinPackingV1::generate(){

        PatternCollectionContainer PC;
        packers[packer_selection]->set_pdb_max_size(get_pdb_max_size());
        vector<Pattern> paterns = packers[packer_selection]->bin_packing(task_proxy);

        for (size_t i=0; i< paterns.size(); i++)
            PC.add_pc(paterns[i]);
        return PC;
    }


static options::PluginTypePlugin<PatternCollectionGeneratorBinPackingV1> _type_plugin(
"PatternCollectionGeneratorBinPackingV1",
"The various pattern selection algorithms usable by the complementary_modular_pdbs heurisitic.");

//static PluginTypePlugin<PatternGenerator> _type_plugin_single(
//    "PatternGenerator",
//    "Factory for single patterns");
//}
//    static PluginShared<PatternCollectionGenerator> _plugin("modular_complementary", _parse);
  static shared_ptr<PatternCollectionGeneratorComplementary>_parse(options::OptionParser &parser) {
    parser.add_option<int> ("packer_selection", "Which of the available packers are we calling.", "0");
    parser.document_synopsis(
        "Pattern Generator BinPackingV1",
        "Bin Packing selection of variables to generate Pattern Collections based on the work of ...");
    options::Options opts = parser.parse();
    if (parser.dry_run())
      return 0;
    return make_shared<PatternCollectionGeneratorBinPackingV1>(opts);
  }
  static options::PluginShared<PatternCollectionGeneratorComplementary> _plugin("BinPackingV1", _parse);
}
