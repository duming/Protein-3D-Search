#include "Protein_3D_Index.hpp"


int Protein_3D_Index:: BuildIndex()
{
    // 1. load cath list file
    cdata.readList();

    // 2. create threads to read pdb files into buffer
    //  and threads that read the buffer and calculate the GaussIntegral
    //  the reader threads also responsible for writing the c-alpha atoms'
    //  coordinates into intermediate files
    pthread_t threadc[C_NUM];
    pthread_t threadp[P_NUM];
    int iret1,iret2;

    for(int i = 0; i<C_NUM;i++)
    {
	    iret1 = pthread_create( &threadc[i], NULL, Consumer, (void*) RWBuf);
	    if(iret1)
	    {
	        fprintf(stderr,"Error - pthread_create() return code: %d\n",iret1);
	        exit(EXIT_FAILURE);
	    }
    }

    for(int i=0; i<P_NUM;i++)
    {
	    iret2 = pthread_create( &threadp[i], NULL, PDBReader, (void*) RWBuf);
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


    // 3. construct the index according to the Gauss Integral

    return 0;
}
