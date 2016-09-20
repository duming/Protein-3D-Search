#ifndef UTILITY
#define UTILITY
#include <vector>
#include <string>
#include <iostream>
#include <sstream> 
#include <dirent.h>
#include <semaphore.h>
#include <pthread.h>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <sstream>

typedef
struct point
{
    double x;
    double y;
    double z;
}POINT;

std::ostream& operator <<(std::ostream& os, point & pt);


///////////////////////////////////////
//      semaphore
//      //////////////////////////////

//due to mac OS doen't support <semaphore.h> well. I write one by myself.
//TODO use boost to handle all the multithreading things.
class Semaphore
{
public:
    Semaphore(int val)
    {
        value = val;
        pthread_mutex_init(&m, NULL);
        pthread_cond_init(&c, NULL);
    }

    ~Semaphore()
    {
        pthread_mutex_destroy(&m);
        pthread_cond_destroy(&c);
    }

    void P()
    {
        pthread_mutex_lock(&m);
        while(value <= 0)
            pthread_cond_wait(&c,&m);
        value--;
        pthread_mutex_unlock(&m);
    }

    void V()
    {
        pthread_mutex_lock(&m);
        value++;
        if(value > 0)
            pthread_cond_signal(&c);
        pthread_mutex_unlock(&m);
    }

private:
    int value;
    pthread_mutex_t m;
    pthread_cond_t c;
};



////////////////////////////////////////
//      read write buffer
//      //////////////////////////////////
template<class T>
class ReadWriteBuffer
{
public:
    ReadWriteBuffer(int bufSize):full(0), empty(bufSize)
    {
        size = bufSize;
        //allocate memory
        isAllocateMemory = true;
        base = new T[bufSize];
        //set input and output pointer
        inptr = 0;
        outptr = 0;
        
    }

    ReadWriteBuffer(int bufSize, T* bufbase):full(0),empty(bufSize)
    {
        size = bufSize;
        //don't allocate memory 
        isAllocateMemory = false;
        base = bufbase;
        inptr = 0;
        outptr = 0;
    }

    ~ReadWriteBuffer()
    {   if(isAllocateMemory)
            delete [] base;
    }

    //return the next position the inptr and outptr should be
    inline int nextPos(int ptr)
    {
        return (ptr+1)%size;
    }

    //start to read the buffer
    //wait for full semaphre
    //use the mutex to create a critical area 
    //after update the pointer, release the mutex
    const T& startRead()
    {
        full.P();
        pthread_mutex_lock(&rw_mutex);
        int ret = outptr;
        outptr = nextPos(outptr);
        pthread_mutex_unlock(&rw_mutex);
        return base[ret];
    }

    void endRead()
    {
        empty.V();
    }


    //start to write 
    T& startWrite()
    {   
        empty.P();
        pthread_mutex_lock(&rw_mutex);
        int ret = inptr;
        inptr = nextPos(inptr);
        pthread_mutex_unlock(&rw_mutex);
        return base[ret];
    
    }

    void endWrite()
    {
        full.V();
    }


private:
    T *base;
    int inptr;
    int outptr;

    int size;
    bool isAllocateMemory;
    Semaphore full;
    Semaphore empty;
    pthread_mutex_t rw_mutex;
};


///////////////////////////////////////////
//      utility class
//      ///////////////////////////////////////


class Utility
{
    public:
    Utility();
    ~Utility();
    
    //standard version with memory allocation
    static std::vector<std::string> split(std::string &s, char delim =' ');

    //fast version of split store result in tokens
    static void split(std::vector<std::string> &tokens, std:: string &s, char delim = ' ') ;

    //even faster version process on char*
    static void split(std::vector<std::string> &tokens,const char* s, char delim = ' ');

    //if isFile == True return all the files in current directory through file
    //otherwise return all sub-directories
    static bool listFiles(std::string directory, std::vector<std::string> & files, bool isFile = true);



    // read a pdb file 
    // read all ATOM line in pdb file
    // save each line to a string vector
    static bool readPDB(std::string fileName , std::vector<std::string> & data, bool ca_only = true);

    // save the pdb data to disk
    // by default it will save the whole file
    // otherwise it will save the part from residue start to residue end
    static void writePDB(std::string fileName, std::vector<std::string> & data, int start = 1, int end = -1);

    // read a txt file save each line in line vector data
    static bool readFile(std::string fileName, std::vector<std::string> & data);

    // write a string vector to file
    // each unit in the vector is a line 
    static bool writeFile(std::string fileName, std::vector<std::string> & data);

    // delete all characters after last '.'
    // if there is no '.' in string, string remain unchanged
    static std::string removeSuff(std::string str);


    static bool readTable(std::string fileName
                , std::vector<std::vector<std::string> >&data
                , std::vector<std::string> colum_names);



    // return the first residue in a pdb line vector
    static int firstRes(std::vector<std::string> & dataPDB);

    // return the last residue in the pdb line vector
    static int lastRes(std::vector<std::string> & dataPDB);


    // split pdb data into segments each segment contains segLen residues
    static void splitPDB(std::vector<std::string> & dataPDB
            , std::vector<std::vector<std::string> >& segs, int segLen);

    template<class T>
    static double l1_norm(T* op1, T* op2, int length)
    {
        double sum = 0;
        for(int i=0; i < length; i++)
            sum += fabs(op1[i] - op2[i]);
        return sum/length;
    }

    template<class T>
    static double l2_norm(T* op1, T* op2, int length)
    {
        double sum = 0;
        for(int i=0; i < length; i++)
            sum += pow((op1[i] - op2[i]),2);
        return sqrt(sum)/length;
    }


    //vector length
    static double vectLen(POINT const &op);

    //vector dot product
    static double vectDot(POINT const &op1, POINT const &op2);

    //vector unit, return v/|v|
    static void vectUnit(POINT &result, POINT const &op);

    //vector subtraction generate the vector that from op1 to op2
    static void vectSub(POINT& result, POINT const &op1, POINT const &op2);

    //vector cross product
    static void vectCross(POINT& result, POINT const &op1, POINT const &op2);
    
    

    // significant digits
    static double round_to_digits(double value, int digits);


    //read gauss integral reasult from the original GI.c program
    static void readGIs(std::string fileName, std::vector<std::vector<double> > &result, bool isTrans = true);




};



#endif
