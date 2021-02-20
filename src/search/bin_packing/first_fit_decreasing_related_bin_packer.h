#ifndef FIRST_FIT_DECREASING_RELATED_BIN_PACKER_H
#define FIRST_FIT_DECREASING_RELATED_BIN_PACKER_H

#include "bin_packer.h"
#include "../abstract_task.h"

#include <vector>

using namespace std;

namespace pdbs {

    class FirstFitDecreasingRelatedBinPacker : public BinPacker {

    private:
        double pdb_max_size;
        int num_pbd_collections;

    public:
        FirstFitDecreasingRelatedBinPacker(double pdb_size, int num_collections);
        ~FirstFitDecreasingRelatedBinPacker() {}
        virtual vector<Pattern> bin_packing(std::shared_ptr<TaskProxy> task_proxy);
        virtual void set_pdb_max_size(double size) {pdb_max_size = size;}

        static bool compare_pattern_length(Pattern one, Pattern two) {
            return (count(two.begin(), two.end(), true)<count(one.begin(), one.end(), true)); }

        static bool compare_domain_size_variable(const pair<int,double> one, const pair<int,double> two) {
            return (one.second > two.second);
        }

    };
}

#endif //FIRST_FIT_DECREASING_BIN_PACKER_H
