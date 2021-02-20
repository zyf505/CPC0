#ifndef MODULAR_PDBS_TYPES_H
#define MODULAR_PDBS_TYPES_H

#include <memory>
#include <vector>

namespace pdbs3 {
class PatternDatabaseInterface;
using Pattern = std::vector<int>;
using PatternCollection = std::vector<Pattern>;
using PDBCollection = std::vector<std::shared_ptr<PatternDatabaseInterface>>;
using MaxAdditivePDBSubsets = std::vector<PDBCollection>;

/*std::ostream & operator<<(std::ostream &os,const PDBCollection &pdbs){
    for (size_t collection=0; collection<pdbs.size();collection++){
		os<<"\tpdb["<<collection<<"]"<<*pdbs[collection]<<endl;
	} 
    return os;
}*/
}

#endif
