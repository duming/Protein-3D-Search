#ifndef PROTEIN_3D_INDEX_H
#define PROTEIN_3D_INDEX_H

#include "Utility.hpp"
#include "CathData.hpp"

//number of producer
#define P_NUM 1
//number of consumer
#define C_NUM 8

#define PROTEIN_DEFAULT_LENGTH 1024
#define BUFFER_DEFAULT_LENGTH 100
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
    Protein_3D_Index(string fileName):cdata(fileName)
    {
    }
    ~Protein_3D_Index()
    {
    }

    int BuildIndex();

    static void* PDBReader(void* ptr);

    static void* Consumer(void*ptr);

private:
    // create threads to read pdb files into buffer
    //  and threads that read the buffer and calculate the GaussIntegral
    //  the reader threads also responsible for writing the c-alpha atoms'
    //  coordinates into intermediate files
    int MultiCalGI();
private:
    CathData cdata;    
    ReadWriteBuffer<Buf_Unit_Points> *RWBuf;
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
