#include "Utility.hpp"
#include <iostream>
#include <fstream>
#include "CathData.hpp"
#include <time.h>
#include <unistd.h>
using namespace std;

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

int main()
{
    string ss("123  333 44 55");
/*
    const char *s1 = "123456";
    char s4[] = "123 : 333 44 5";

    clock_t t1,t2;
    vector<string> files;

   CathData cc("CathDomainList.v4.0.0.txt");

   t1 = clock();
   cc.readList();
   t2 = clock();

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
	
    semaphoreTest();     
    return 0;
}


