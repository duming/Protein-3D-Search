#ifndef PROTEIN_3D_INDEX_H
#define PROTEIN_3D_INDEX_H

#include "Utility.hpp"
#include "CathData.hpp"
#include "lshbox/lshbox.h"
//number of producer
#define P_NUM 1
//number of consumer
#define C_NUM 8

#define PROTEIN_DEFAULT_LENGTH 1024
#define BUFFER_DEFAULT_LENGTH 100

typedef double INDEX_DATA_TYPE;
using namespace std;



class Buf_Unit_Points
{
public:
    bool isEnd;
    int domain_idx;
    vector<POINT> points;
};


class Protein_3D_Index
{
public:

    Protein_3D_Index():cdata("null")
    {
    }

    Protein_3D_Index(string fileName):cdata(fileName),lshIndex()
    {
    }
    ~Protein_3D_Index()
    {
    }
    //one function for Protein descriptor calculation and lsh index making
    int BuildIndex();



    // Only calculate the descriptors for all pdb files under dataPath directory
    // Save the descriptors into desFileName.
    // Generate a list file that indicate the order of the pdb files
    void desCalculate(string dataPath, string desFileName, string listName)
    {
        cdata.fakeList(dataPath);
        MultiCalGI();
        cdata.saveDescriptor(desFileName);
        cdata.saveList(listName);
    }

    static void* PDBReader(void* ptr);

    static void* Consumer(void*ptr);


    lshbox::kdbqLsh<INDEX_DATA_TYPE>& getIndex()
    {
        return lshIndex;
    }

private:
    // create threads to read pdb files into buffer
    //  and threads that read the buffer and calculate the GaussIntegral
    //  the reader threads also responsible for writing the c-alpha atoms'
    //  coordinates into intermediate files
    int MultiCalGI();

    //build the Locality sensitive hashing index using the lshbox tool
    //input: 1. the data matrix that needs to be hashed
    //       2. hash table size
    //output: save the index file to indexfile
    void LSH(lshbox::Matrix<INDEX_DATA_TYPE> &data, string indexFile, int htSize);
private:
    CathData cdata;    
    ReadWriteBuffer<Buf_Unit_Points> *RWBuf;
    // descriptor data matrix
    lshbox::kdbqLsh<INDEX_DATA_TYPE> lshIndex;
};


struct threadPara
{
    Protein_3D_Index * P3Dptr;
    int startPos;
    int endPos;
    pthread_mutex_t * mt;
    //the number of the producer that read from file
    int p_num;
    //the number of consumer that read buffer and calculate the gauss integral
    int c_num;
};

#endif
