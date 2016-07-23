#ifndef UTILITY
#define UTILITY
#include <vector>
#include <string>
#include <iostream>
#include <sstream> 
#include <dirent.h>
#include <semaphore.h>
#include <pthread.h>

typedef
struct point
{
    double x;
    double y;
    double z;
}POINT;



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


template<class T>
class ReadWriteBuffer
{
public:
    ReadWriteBuffer(int bufSize):full(0), empty(bufSize)
    {
        size = bufSize;
        //allocate memory
        base = new T[bufSize];
        //set input and output pointer
        inptr = 0;
        outptr = 0;
        
    }
    ~ReadWriteBuffer()
    {
        delete [] base;
    }

    //return the next position the inptr and outptr should be
    inline int nextPos(int ptr)
    {
        return (ptr+1)%size;
    }

    //read buffer operation
    const T& startRead();
    void endRead();

    //write buffer operation
    T& startWrite();
    void endWrite();

private:
    T *base;
    int inptr;
    int outptr;

    int size;
    Semaphore full;
    Semaphore empty;
    pthread_mutex_t rw_mutex;
};





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
};



#endif
