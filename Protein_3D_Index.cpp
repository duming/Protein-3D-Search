#include "Protein_3D_Index.hpp"


    

void* Protein_3D_Index::
PDBReader(void* ptr)
{
    threadPara *TPptr = (threadPara*)ptr;
    Protein_3D_Index* Pptr = TPptr->P3Dptr;
    CathData * Cptr = &Pptr->cdata;
    vector<POINT> tempPoints;
    tempPoints.reserve(PROTEIN_DEFAULT_LENGTH);

    int start = TPptr->startPos;
    int end = TPptr->endPos;
    //increase the producer thread number
    pthread_mutex_lock(TPptr->mt);
    TPptr->p_num++;
    pthread_mutex_unlock(TPptr->mt);


    for(int i = start; i<=end ; i++)
    {
        //read from pdb file
        if(!(*Cptr)[i].readPDB(tempPoints))
            continue;
        // write them into buffer 
        Buf_Unit_Points& Unit = Pptr->RWBuf->startWrite();
        Unit.domain_idx = i;
        // Unit.isEnd = false;
        Unit.points = tempPoints;
        Pptr->RWBuf->endWrite();
    }


    pthread_mutex_lock(TPptr->mt);
    //check if is the last PDBreader thread
    TPptr->p_num--;
    if(TPptr->p_num == 0)
        for(int i=0;i<TPptr->c_num;i++)
        {//put stop tokens for consumer threads;
            Buf_Unit_Points& Unit = Pptr->RWBuf->startWrite();
            Unit.isEnd = true;
            Pptr->RWBuf->endWrite();
        }
    pthread_mutex_unlock(TPptr->mt);
    return NULL;
}


void* Protein_3D_Index::
Consumer(void*ptr)
{
    threadPara *TPptr = (threadPara*)ptr;
    Protein_3D_Index* Pptr = TPptr->P3Dptr;
    CathData * Cptr = &Pptr->cdata;
    vector<POINT> tempPoints;
    tempPoints.reserve(PROTEIN_DEFAULT_LENGTH);
    int domain_idx;

    //increase the consumer thread number 
    pthread_mutex_lock(TPptr->mt);
    TPptr->c_num++;
    pthread_mutex_unlock(TPptr->mt);

    while(true)
    {
        const Buf_Unit_Points& unit = Pptr->RWBuf->startRead();
        if(true == unit.isEnd)
            break;
        tempPoints = unit.points; 
        domain_idx = unit.domain_idx;
        Pptr->RWBuf->endRead();

        //for test
        pthread_t tid;
        tid = pthread_self();
        pthread_mutex_lock(TPptr->mt);
        cout<<tid<<":read "<<endl;
        for(int i=0;i<tempPoints.size();i++)
            cout<<tempPoints[i]<<endl;
        cout<<endl;
        pthread_mutex_unlock(TPptr->mt);

    }
    //decrease the consumer thread number 
    pthread_mutex_lock(TPptr->mt);
    TPptr->c_num--;
    pthread_mutex_unlock(TPptr->mt);


    return NULL;
}


int Protein_3D_Index:: MultiCalGI()
{
    pthread_t threadc[C_NUM];
    pthread_t threadp[P_NUM];
    int iret1,iret2;
    
    //prepare the thread parameter
    threadPara TP;
    pthread_mutex_t mt;
    pthread_mutex_init(&mt,NULL);
    TP.P3Dptr = this;
    TP.mt = &mt;
    TP.p_num=0;
    TP.c_num=0;

    //init the buffer
    Buf_Unit_Points bufbase[BUFFER_DEFAULT_LENGTH];
    for(int i = 0; i< BUFFER_DEFAULT_LENGTH; i++)
    {
        bufbase[i].points.reserve(PROTEIN_DEFAULT_LENGTH);
        bufbase[i].isEnd = false;
        bufbase[i].domain_idx = -1;
    }
    RWBuf = new ReadWriteBuffer<Buf_Unit_Points>(BUFFER_DEFAULT_LENGTH, bufbase);


    for(int i = 0; i<C_NUM;i++)
    {
        //for test
        TP.startPos = 0;
        TP.endPos = 100;

        //for test end
	    iret1 = pthread_create( &threadc[i], NULL, Consumer, (void*) &TP);
	    if(iret1)
	    {
	        fprintf(stderr,"Error - pthread_create() return code: %d\n",iret1);
	        exit(EXIT_FAILURE);
	    }
    }

    for(int i=0; i<P_NUM;i++)
    {
	    iret2 = pthread_create( &threadp[i], NULL, PDBReader, (void*) &TP);
	    if(iret2)
	    {
	        fprintf(stderr,"Error - pthread_create() return code: %d\n",iret2);
	        exit(EXIT_FAILURE);
	    }
    }

    for(int i=0;i<C_NUM;i++)
        pthread_join(threadc[i],NULL);
    for(int i=0;i<P_NUM;i++)
        pthread_join(threadp[i],NULL);

    //free memory
    delete RWBuf;
return 1;
}

int Protein_3D_Index:: BuildIndex()
{
    // 1. load cath list file
    cdata.readList();
    cdata.printDomains(0,99);

    // 2. calculate GI
    MultiCalGI();
    
    // 3. construct the index according to the Gauss Integral

    return 0;
}
