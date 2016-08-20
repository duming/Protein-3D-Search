#include "Utility.hpp"
#include <iostream>
#include <fstream>
#include "CathData.hpp"
#include "Protein_3D_Index.hpp"
#include "GaussIntegral.hpp"
#include <time.h>
#include <unistd.h>
#include "lshbox/lshbox.h"

using namespace std;
///////////////////////////////////
//        semaphore test
//        ///////////////////
pthread_mutex_t main_mutex;

void *print_message_function( void *ptr )
{
    char *message;
    message = (char *) ptr;
    printf("%s \n", message);
    return NULL;
}

void * wait(void *ptr)
{
    Semaphore *sem;
    sem = (Semaphore*) ptr;
    pthread_t   tid;
    tid = pthread_self();
    
    pthread_mutex_lock(&main_mutex);
    cout<<tid<<":start to wait"<<endl;
    pthread_mutex_unlock(&main_mutex);
    
    sem->P();
   
    pthread_mutex_lock(&main_mutex);
    cout<<tid<<":finish waiting"<<endl;
    pthread_mutex_unlock(&main_mutex);
    
    return NULL;
}


void * sig(void *ptr)
{
    sleep(1);
    Semaphore *sem;

    sem = (Semaphore*) ptr; 
    pthread_t   tid;
    tid = pthread_self();
    
    pthread_mutex_lock(&main_mutex);
    cout<<tid<<":start to sig"<<endl;
    pthread_mutex_unlock(&main_mutex);
    
    sem->V();
    
    pthread_mutex_lock(&main_mutex);
    cout<<tid<<":finish sig"<<endl;
    pthread_mutex_unlock(&main_mutex);
    
    return NULL;
}

int semaphoreTest()
{
    pthread_mutex_init(&main_mutex,NULL);
    int p_num=3,v_num=2;
    Semaphore sem(3);
    pthread_t thread1[p_num], thread2[v_num];
	const char *message1 = "Thread 1";
	const char *message2 = "Thread 2";
	int  iret1, iret2;
    
    for(int i = 0; i<p_num;i++)
    {
	    iret1 = pthread_create( &thread1[i], NULL, wait, (void*) &sem);
	    if(iret1)
	    {
	        fprintf(stderr,"Error - pthread_create() return code: %d\n",iret1);
	        exit(EXIT_FAILURE);
	    }
    }

    for(int i=0; i<v_num;i++)
    {
	    iret2 = pthread_create( &thread2[i], NULL, sig, (void*) &sem);
	    if(iret2)
	    {
	        fprintf(stderr,"Error - pthread_create() return code: %d\n",iret2);
	        exit(EXIT_FAILURE);
	    }
    }
	
    for(int i=0;i<p_num;i++)	
	    pthread_join( thread1[i], NULL);

    for(int i=0;i<v_num;i++)
	    pthread_join( thread2[i], NULL);


    return 0;
}



/////////////////////////////////////////////
//    read write buffer test
//    ///////////////////////////////////
void* producer(void* ptr)
{
    //sleep(1);
    ReadWriteBuffer<int> * bptr = (ReadWriteBuffer<int>*) ptr;
    pthread_t   tid;
    tid = pthread_self();
 
    for(int i=0; i<100;i++)
    {
        //start writing
        pthread_mutex_lock(&main_mutex);
        cout<<tid<<":attempt to write "<<endl;
        pthread_mutex_unlock(&main_mutex);
        int &tmp = bptr->startWrite();
        tmp = i;
        bptr->endWrite();
        
        pthread_mutex_lock(&main_mutex);
        cout<<tid<<":write "<<i<<endl;
        pthread_mutex_unlock(&main_mutex);
 
    }
    return NULL;
}


void* consumer(void* ptr)
{
    pthread_t   tid;
    tid = pthread_self();
    

    ReadWriteBuffer<int>* bptr = (ReadWriteBuffer<int>*) ptr;
    for(int i =0;i<100;i++)
    {   
        pthread_mutex_lock(&main_mutex);
        cout<<tid<<":attempt to read "<<endl;
        pthread_mutex_unlock(&main_mutex);

        int tmp = bptr->startRead();
        bptr->endRead();
   
        pthread_mutex_lock(&main_mutex);
        cout<<tid<<":read "<<tmp<<endl;
        pthread_mutex_unlock(&main_mutex);
    }
    return NULL;
}


void producer_consumer_test()
{
    pthread_mutex_init(&main_mutex,NULL);
    ReadWriteBuffer<int> buffer(3);
    int c_num = 2;
    int p_num = 2;
    pthread_t threadc[c_num];
    pthread_t threadp[p_num];
    int iret1,iret2;

    for(int i = 0; i<c_num;i++)
    {
	    iret1 = pthread_create( &threadc[i], NULL, consumer, (void*) &buffer);
	    if(iret1)
	    {
	        fprintf(stderr,"Error - pthread_create() return code: %d\n",iret1);
	        exit(EXIT_FAILURE);
	    }
    }

    for(int i=0; i<p_num;i++)
    {
	    iret2 = pthread_create( &threadp[i], NULL, producer, (void*) &buffer);
	    if(iret2)
	    {
	        fprintf(stderr,"Error - pthread_create() return code: %d\n",iret2);
	        exit(EXIT_FAILURE);
	    }
    }

    for(int i=0;i<c_num;i++)
        pthread_join(threadc[i],NULL);
    for(int i=0;i<p_num;i++)
        pthread_join(threadp[i],NULL);

}


//////////////////////////////////
//test vector operations
///////////////////////////

void vectortest()
{
    POINT vct1 = {0,0,0};
    POINT vct2 = {-1,2,3};
    POINT result;

    cout<<Utility::vectLen(vct1)<<endl;
    cout<<Utility::vectDot(vct1,vct2)<<endl;
    Utility::vectUnit(result,vct1);
    cout<<result<<endl;

    cout<<Utility::vectLen(result)<<endl;

    Utility::vectSub(result,vct1,vct2);
    cout<<result<<endl;
    Utility::vectCross(result,vct1,vct2);
    cout<<result<<endl;

}



void DescriptorSavingTest()
{
    //Protein_3D_Index P3DIdx("CathDomainList.v4.0.0.txt");
    //P3DIdx.BuildIndex();

    
}


void lshtest()
{
    typedef double DATATYPE;
    std::cout << "LOADING DATA ..." << std::endl;
    lshbox::timer timer;
    lshbox::Matrix<DATATYPE> data("sampleData");
    std::cout << "LOAD TIME: " << timer.elapsed() << "s." << std::endl;

    std::cout << "CONSTRUCTING INDEX ..." << std::endl;
	timer.restart();
	std::string file("index/testindex");
	bool use_index = false;
	lshbox::psdLsh<DATATYPE> mylsh;
	if (use_index)
	{
	    mylsh.load(file);
	}
	else
	{
	    lshbox::psdLsh<DATATYPE>::Parameter param;
	    param.M = 1000;
	    param.L = 100;
	    param.D = 30;
	    param.T = 2;
	    param.W = 10;
	    mylsh.reset(param);
	    //mylsh.train(data);
	    mylsh.hash(data);
	}
	mylsh.save(file);
	std::cout << "CONSTRUCTING TIME: " << timer.elapsed() << "s." << std::endl;
	std::cout << "LOADING BENCHMARK ..." << std::endl;
}


int main()
{
    string ss("123  333 44 55");
/*
    const char *s1 = "123456";
    char s4[] = "123 : 333 44 5";

    clock_t t1,t2;
    vector<string> files;

   
    cc.test();
*/
   //cout<<t2-t1<<endl;
/*    
    vector<string> tks;
    Utility::split(tks, s4, ':');
    for(int i=0; i<tks.size(); i++)
        cout<<tks[i]<<endl;
*/
   // ReadWriteBuffer<POINT> pbf(10);
	
    //semaphoreTest();     


   // producer_consumer_test();
   

    Protein_3D_Index P3DIdx("CathDomainList.v4.0.0.txt");
    P3DIdx.BuildIndex();
    

//    GaussIntegral gi(1000,true);
//    gi.test();
    



//    vectortest();
    


//    CathData cd("CathDomainList.v4.0.0.txt");
//    cd.test();

//    lshtest();
    
   // DescriptorSavingTest();
   

    /*
    unsigned arry[3] = {8,281,30};
    unsigned arry2[3];
    cout<<arry[0]<<'\t'<<arry[1]<<'\t'<<arry[2]<<endl;

    ofstream of("test.data",ofstream::binary);
    of.write((char*)arry, sizeof arry);
    of.close();

    ifstream infile("test.data",std::ios::binary );
    infile.read((char*)arry2, sizeof(arry));
    infile.close(); 
    cout<<arry2[0]<<'\t'<<arry2[1]<<'\t'<<arry2[2]<<endl;

    typedef double DATATYPE;
    lshbox::timer timer;
    lshbox::Matrix<DATATYPE> data("test.data");
*/

    return 0;
}


