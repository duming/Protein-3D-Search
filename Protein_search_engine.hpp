#pragma once 
#include "lshbox/lshbox.h"
#include "Cathdata.hpp"
#include "utility.hpp"
#include "Protein_3D_Index.hpp"

class Index_query
{
    public:
        Index_query(int k, lshbox::Matrix<INDEX_DATA_TYPE> & desMat, std::string indexFile):
            desMatrix(desMat)
        {
            K = k;
            lshIndex.load(indexFile);
            lshbox::Matrix<INDEX_DATA_TYPE>::Accessor accessor(desMat);
            lshbox::Metric<INDEX_DATA_TYPE> metric(desMat.getDim(), L2_DIST);
            scanner = new lshbox::Scanner<lshbox::Matrix<INDEX_DATA_TYPE>::Accessor>(
                    accessor,
                    metric,
                    K
            );

        }


        Index_query(int k, lshbox::Matrix<INDEX_DATA_TYPE>& desMat, lshbox::kdbqLsh<INDEX_DATA_TYPE> & index):
            desMatrix(desMat)
        {
            K = k;
            lshIndex = index;
            lshbox::Matrix<INDEX_DATA_TYPE>::Accessor accessor(desMat);
            lshbox::Metric<INDEX_DATA_TYPE> metric(desMat.getDim(), L2_DIST);
            scanner = new lshbox::Scanner<lshbox::Matrix<INDEX_DATA_TYPE>::Accessor>(
                    accessor,
                    metric,
                    K
            );

        }

        ~Index_query()
        {
            delete scanner;
        }

        //query with descriptor
        void desQuery(std::vector<std::pair<float, unsigned> > & result,INDEX_DATA_TYPE * query);

        
        lshbox::kdbqLsh<INDEX_DATA_TYPE> & getIndex()
        {
            return lshIndex;
        }

    private:
        //number of nearest points
        int K;
        //descriptor matrix
        lshbox::Matrix<INDEX_DATA_TYPE>& desMatrix;
        //the actual index
        lshbox::kdbqLsh<INDEX_DATA_TYPE> lshIndex;

        //the scanner that required by lshbox
        lshbox::Scanner<lshbox::Matrix<INDEX_DATA_TYPE>::Accessor> *scanner;
};
