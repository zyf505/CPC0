#ifndef MODULAR_PDBS_PATTERN_COLLECTION_GENERATOR_BIN_PACKING_V1
#define MODULAR_PDBS_PATTERN_COLLECTION_GENERATOR_BIN_PACKING_V1

#include "pattern_collection_generator_complementary.h"
#include "types.h"
#include "../globals.h"
#include "../task_tools.h"
#include "../task_proxy.h"
#include "../bin_packing/bin_packer.h"

#include <memory>
#include <vector>
#include <random>
#include <iostream>
#include <iterator>

using namespace pdbs;

class AbstractTask;

namespace options {
class Options;
}
namespace utils {
class CountdownTimer;
}

namespace pdbs3 {

    class PatternCollectionGeneratorBinPackingV1 : public PatternCollectionGeneratorComplementary {
    private:
        int packer_selection;
        std::vector<std::shared_ptr<BinPacker>> packers;
        std::shared_ptr<TaskProxy> task_proxy;
    public:
        PatternCollectionGeneratorBinPackingV1();
        explicit PatternCollectionGeneratorBinPackingV1(const options::Options &options);
        virtual void initialize(std::shared_ptr<AbstractTask> task) override;
        virtual PatternCollectionContainer generate() override;
    
  /*void transform_to_pattern_bitvector_form(std::vector<bool> &bitvector,const std::vector<int> &pattern) const {
      bitvector.assign(g_variable_name.size(), false);
      for (size_t i = 0; i < pattern.size(); ++i) {
        bitvector[pattern[i]]=true;
      }
    }

    Pattern transform_to_pattern_normal_form(const std::vector<bool> &bitvector) const {
      Pattern pattern;
      for (size_t i = 0; i < bitvector.size(); ++i) {
          if (bitvector[i]){
        pattern.push_back(i);
          }
      }
      return pattern;
    }*/

    };

}
#endif
