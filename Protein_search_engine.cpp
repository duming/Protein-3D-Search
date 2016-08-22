#include "Protein_search_engine.hpp"


void Index_query::desQuery(std::vector<std::pair<float, unsigned> > & result, INDEX_DATA_TYPE *query)
{
    scanner->reset(query);
    lshIndex.query(query, *scanner);
    lshbox::Topk & tpk = scanner->topk();
    //tpk.genTopk();
    std::vector<std::pair<float, unsigned> > & tops = tpk.getTopk();
    result = tops;
}
