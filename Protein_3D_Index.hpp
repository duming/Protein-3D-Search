#ifndef PROTEIN_3D_INDEX_H
#define PROTEIN_3D_INDEX_H

#include "Utility.hpp"
#include "CathData.hpp"

//number of producer
#define P_NUM 1
//number of consumer
#define C_NUM 2

using namespace std;


struct Buf_Unit_Points
{
    bool isEnd;
    vector<double> points;
};

class Protein_3D_Index
{
public:
    Protein_3D_Index(string fileName):cdata(fileName)
    {
    }
    ~Protein_3D_Index();

    int BuildIndex();

    static void* PDBReader(void* ptr);

    static void* Consumer(void*ptr);
private:
    CathData cdata;    
    ReadWriteBuffer<Buf_Unit_Points> *RWBuf;
};

#endif
