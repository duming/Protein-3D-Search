#include "GaussIntegral.hpp"

#include <time.h>
#include <iomanip>
using namespace std;

double GaussIntegral::Epsilon = 1e-16;



void GaussIntegral::test()
{
    /*
    //array test
    initAll();
    for(int i = 0; i < maxLen; i++)
        for(int j = 0; j < maxLen; j++)
            omega[i][j] = i*maxLen + j;

    printArray<double>(omega, maxLen);
    freeAll();
*/

    /*
    int testSize = 1000;
    int rtimes = 400;
    SquareArray<double> arraym(testSize);
    double arrayc[testSize*testSize];
    vector<vector<double> >arrayv;

    arrayv.reserve(testSize);
    for(int i=0;i<testSize;i++)
        arrayv[i].reserve(testSize);

    clock_t t1,t2;

    t1 = clock();
    for(int i=0; i < rtimes; i++)
    {
        for(int j = 0; j < testSize; j++)
            for(int k = 0; k < testSize; k++)
                arraym(j,k) = j*testSize +k;
    }
    t2 = clock();
    cout<<1000*(t2 -t1)/CLOCKS_PER_SEC<<endl;
 
    t1 = clock();
    for(int i=0; i < rtimes; i++)
    {
        for(int j = 0; j < testSize; j++)
            for(int k = 0; k < testSize; k++)
                arrayc[j*testSize+k] = j*testSize +k;
    }
    t2 = clock();
    cout<<1000*(t2 -t1)/CLOCKS_PER_SEC<<endl;
 
    t1 = clock();
    for(int i=0; i < rtimes; i++)
    {
        for(int j = 0; j < testSize; j++)
            for(int k = 0; k < testSize; k++)
                arrayv[j][k] = j*testSize +k;
    }
    t2 = clock();
    cout<<1000*(t2 -t1)/CLOCKS_PER_SEC<<endl;
    */
    

    //array speed test

    //correctness test
    string path = "data/testpdb/";
    vector<string> files;
    Utility::listFiles(path,files);

    vector<vector<double> > GIresult;

    vector<POINT> *crdPtr = new vector<POINT>[files.size()];
    

    for(int i=0;i<files.size();i++)
    {
        crdPtr[i].resize(1000);
    }

    string GIfile;
     int ds_len;
    if(CylinderTransform)
    {
        ds_len = 30;
        GIfile = "data/GI30out.txt";
    }
    else
    {
        GIfile = "data/GI29out.txt";
        ds_len = 29;
    }


    
    Utility::readGIs(GIfile, GIresult,CylinderTransform);
    cout<<"GIresults:"<<GIresult.size()<<endl;

    double dscrpt[files.size()][ds_len]; 
    
    clock_t t1,t2;

    t1 = clock();
    for(int i=0; i < files.size(); i++)
    {
        cout<<"reading:"<<path+files[i]<<endl;
        Cathdomain::readPDB(path+files[i], crdPtr[i]);
        //cout<<"calculating GI"<<endl;
        setProtein(&crdPtr[i]);
        GaussAll(dscrpt[i]);
    }
    t2 = clock();
    cout<<1000*(t2 -t1)/CLOCKS_PER_SEC<<endl;

    for(int i=0; i < GIresult.size(); i++)
    {
        for(int j = 0; j < ds_len ;j++)
        {
            double error,erat;
            error = Utility::round_to_digits(dscrpt[i][j],4) - GIresult[i][j];
            erat = abs(error)/(min(abs(dscrpt[i][j]),abs(GIresult[i][j]))); 
            if(error > 0.005 || erat > 0.001)
            {
                cout<<files[i]<<":";

                cout<<j<<"\t:\t"<<setw(20)<<GIresult[i][j]<<"\t"<<setw(10)<<dscrpt[i][j]<<"\t\t";

                cout<<error<<'\t'<<erat<<endl;
            }
        }

        cout<<endl;
    }
}

template<class T>
void GaussIntegral::printArray(T ** array, int size)
{
    for(int i=0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
            cout<<array[i][j]<<'\t';
        cout<<endl;
    }
}

template<class T>
void GaussIntegral::initArray(T *&base, T **&array, int size)
{
    base = new T[size*size];
    array = new T*[size];
    for(int i=0; i < size; i++)
        array[i] = base + size*i;
}

template<class T>
void GaussIntegral::freeArray(T *&base, T **&array)
{
    delete [] array;
    delete [] base;
}



void GaussIntegral::initAll()
{
    initArray<POINT>(unitvectorBase, unitvector, maxLen);
    initArray<double>(omegaBase, omega, maxLen);
    initArray<double>(absomegaBase, absomega, maxLen);
    initArray<double>(partsumBase, partsum, maxLen);
    initArray<double>(abspartsumBase, abspartsum, maxLen);
}


void GaussIntegral::freeAll()
{
    freeArray<POINT>(unitvectorBase, unitvector);
    freeArray<double>(omegaBase, omega);
    freeArray<double>(absomegaBase, absomega);
    freeArray<double>(partsumBase, partsum);
    freeArray<double>(abspartsumBase, abspartsum);
}


void GaussIntegral:: createUnitVector()
{
    POINT tmp;
    for(int i = 0; i < proteinLen-1 ; i++)
        for(int j = i+1 ; j < proteinLen ; j++)
        {
            Utility::vectSub(tmp, (*currentProtein)[i], (*currentProtein)[j]);
            Utility::vectUnit(unitvector[i][j],tmp);
        }
}




double GaussIntegral:: spangle(POINT const &op1, POINT const &op2, POINT const &op3)
{
    POINT vv;
    double ab,bc;
    Utility::vectCross(vv,op1,op2);
    ab = Utility::vectDot(op1, op2);
    bc = Utility::vectDot(op2, op3);
    return atan2(Utility::vectDot(vv, op3), ab*bc - Utility::vectDot(op1,op3));
}


void GaussIntegral:: createomega(void)
{
    double anglesum;
    POINT cross1, cross2, cross3;

    //the adjacent vectors has the value zero
    for(int i =0; i< proteinLen-1 ;i++)
    {
        omega[i][i+1] = 0;
        absomega[i][i+1] = 0;
    }

    for(int i=0;i<proteinLen-2;i++)
    {
        for(int j = i+2; j<proteinLen-1 ; j++)
        {
            Utility::vectCross(cross1, unitvector[i][j], unitvector[i][j+1]);
            Utility::vectCross(cross2, unitvector[i+1][j], unitvector[i+1][j+1]);
            Utility::vectCross(cross3, cross1, cross2);
            if(Utility::vectDot(cross3,cross3)<Epsilon)
            {
                omega[i][j] = 0;
                absomega[i][j] = 0;
            }
            else
            {
                anglesum =  -spangle(unitvector[i][j], unitvector[i][j+1], unitvector[i+1][j+1])
                            -spangle(unitvector[i][j+1], unitvector[i+1][j+1], unitvector[i+1][j])
                            -spangle(unitvector[i+1][j+1], unitvector[i+1][j], unitvector[i][j])
                            -spangle(unitvector[i+1][j], unitvector[i][j], unitvector[i][j+1]);

                if(anglesum >=0)
                    omega[i][j] = 2*M_PI - anglesum;
                else
                    omega[i][j] = -2*M_PI - anglesum;

                absomega[i][j] = fabs(omega[i][j]);
            }//end if
        }//end for j
    }//end for i
}



void GaussIntegral::createpartsum(void)
{
    //i and k in the comment
    int i, k;
    ////////////
    //partsum
    ////////
    //adjacent equals to zero
    for(k = 0; k < 2; k++)
        for(i = 0 ;i< proteinLen-k-1 ;i++)
            partsum[i][i+k] = 0;

    for(i = 0; i<proteinLen -2-1;i++)
        partsum[i][i+2] = omega[i][i+2];

    for(k = 3; k < proteinLen -1; k++)
        for(i = 0; i<proteinLen - k -1; i++)
            partsum[i][i+k] = partsum[i][i+k-1] + partsum[i+1][i+k] - partsum[i+1][i+k-1]
                            + omega[i][i+k];
    
    //////////
    //absolute partsum
    ////////
    for(k = 0; k < 2; k++)
        for(i = 0 ;i< proteinLen-k ;i++)
            abspartsum[i][i+k] = 0;

    for(i = 0; i<proteinLen -2-1;i++)
        abspartsum[i][i+2] = absomega[i][i+2];

    for(k = 3; k < proteinLen -1; k++)
        for(i = 0; i<proteinLen - k-1; i++)
            abspartsum[i][i+k] = abspartsum[i][i+k-1] + abspartsum[i+1][i+k] - abspartsum[i+1][i+k-1]
                            + absomega[i][i+k];

}


double GaussIntegral:: mixedsum(int a, int b, int c, int d, int isabs)
{
    double tempm, tempp;
    if(isabs != 0)
    {
        tempm = abspartsum[a][c-1] + abspartsum[b+1][d];
        tempp = abspartsum[a][d] + abspartsum[b+1][c-1];
    
    }
    else
    {
        tempm = partsum[a][c-1] + partsum[b+1][d];
        tempp = partsum[a][d] + partsum[b+1][c-1];
    }
    return tempp-tempm;
}


double GaussIntegral::int12(int isabs)
{
    if(isabs)
        return abspartsum[0][proteinLen-2]/(M_PI*2);
    else
        return partsum[0][proteinLen-2]/(M_PI*2);
}



double GaussIntegral::int12_34(int isabs1, int isabs2)
{
    int b;
    double temp;
    double** *arry;
    if(isabs1 == 1)
        arry = &abspartsum;
    else
        arry = &partsum;
    temp=0;
    for(b=3;b<proteinLen-3; b++)
        temp=temp+(*arry)[0][b-1]*mixedsum(b,b,b+2,proteinLen-2, isabs2);
    return temp/(M_PI*M_PI*4);
}


double GaussIntegral::int13_24(int isabs1, int isabs2)
{
    int a, b;
    double temp;
    double** * arry;
    if(isabs2 == 1)
        arry = & absomega;
    else
        arry = &omega;
    temp=0;
    for(a=0; a<proteinLen-3; a++)
        for(b=a+2; b<proteinLen -1; b++)
            temp = temp + mixedsum(0,a-1,a+1,b-1, isabs1)*(*arry)[a][b];
    return temp/(M_PI*M_PI*4);
}




double GaussIntegral::int14_23(int isabs1, int isabs2)
{
    int a, b;
    double temp;
    double ** *s_arry,***o_arry;
    if(isabs1 == 1)
        o_arry = & absomega;
    else
        o_arry = & omega;

    if(isabs2 == 1)
        s_arry = & abspartsum;
    else
        s_arry = & partsum;
    temp=0;
    for(a = 0;a<proteinLen-5; a++)
        for(b=a+4; b < proteinLen-1; b++)
            temp=temp+(*s_arry)[a+1][b-1]*(*o_arry)[a][b];
    return temp/(M_PI*M_PI*4);
}






double GaussIntegral::int12_34_56(void)
{
    int a, b;
    double temp1, temp2;
    temp1=0;
    for(a=3; a < proteinLen - 6; a++)
    {
        temp2=0;
        for(b = a+3; b < proteinLen-3; b++)
            temp2 = temp2 + mixedsum(a,a,a+2,b-1)*mixedsum(b,b,b+2,proteinLen - 2);
        temp1=temp1+temp2*partsum[0][a-1];
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}



double GaussIntegral::int12_35_46(void)
{
    int a, b;
    double temp1, temp2;
    temp1=0;
    for(a=3; a < proteinLen-4; a++)
    {
        temp2=0;
        for(b = a+2; b < proteinLen-2; b++)
            temp2 = temp2 + omega[a][b]*mixedsum(a+1,b-1,b+1,proteinLen-2);
        temp1 = temp1 + temp2*partsum[0][a-1];
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}



double GaussIntegral::int12_36_45(void)
{
    int a, b;
    double temp1, temp2;
    temp1=0;
    for(a=3; a < proteinLen-5; a++)
    {
        temp2=0;
        for(b=a+4; b < proteinLen-1 ; b++)
            temp2 = temp2 + omega[a][b]*partsum[a+1][b-1];
        temp1=temp1+temp2*partsum[0][a-1];
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}


double GaussIntegral::int13_24_56(void)
{
    int a, b;
    double temp1, temp2;
    temp1=0;
    for(a=3; a < proteinLen-4;a++)
    {
        temp2=0;
        for(b=1; b < a-1; b++)
            temp2 = temp2 + omega[b][a]*mixedsum(0,b-1,b+1,a-1);
        temp1=temp1+temp2*partsum[a+1][proteinLen-2];
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}


double GaussIntegral::int13_25_46(void)
{
    int a, b, c;
    double temp1, temp2;
    temp1=0;
    for(a=1; a < proteinLen-5; a++)
    {
        for(b = a+2; b < proteinLen-3; b++)
        {
            temp2=0;
            for(c=b+1; c < proteinLen-2; c++)
                temp2 = temp2 + omega[a][c]*mixedsum(b,b,c+1,proteinLen-2);
            temp1 = temp1 + temp2*mixedsum(0,a-1,a+1,b-1);
        }
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}


double GaussIntegral::int13_26_45(void)
{
    int a, b, c;
    double temp1, temp2;
    temp1=0;
    for(a = 1; a < proteinLen-6; a++)
    {
        for(b=a+5; b < proteinLen -1; b++)
        {
            temp2=0;
            for(c=a+1; c < b-3; c++)
                temp2 = temp2 + partsum[c+1][b-1]*mixedsum(0,a-1,c,c);
            temp1 = temp1 + temp2*omega[a][b];
        }
    }
        return temp1/(M_PI*M_PI*M_PI*8);
}



double GaussIntegral::int14_23_56(void)
{
    int a, b;
    double temp1, temp2;
    temp1=0;
    for(a=4; a < proteinLen-4; a++)
    {
        temp2=0;
        for(b=0; b<a-3; b++)
            temp2 = temp2 + omega[b][a]*partsum[b+1][a-1];
        temp1 = temp1 + temp2*partsum[a+1][proteinLen-2];
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}



double GaussIntegral::int14_25_36(void)
{
    int a, b, c;
    double temp1, temp2;
    temp1=0;
    for(a=1; a < proteinLen-5; a++)
    {
        for(b=a+3; b < proteinLen-2; b++)
        {
            temp2=0;
            for(c=a+1; c < b-1; c++)
                temp2 = temp2 + mixedsum(0,a-1,c+1,b-1)*mixedsum(c,c,b+1,proteinLen-2);
            temp1 = temp1 + temp2*omega[a][b];
        }
    }
          return temp1/(M_PI*M_PI*M_PI*8);
}

double GaussIntegral::int14_26_35(void)
{
    int a, b, c;
    double temp1, temp2;
    temp1=0;
    for(a=4; a < proteinLen-2; a++)
    {
        for(b=2; b<a-1; b++)
        {
            temp2 = 0;
            for(c=0; c<b-1; c++)
                temp2 = temp2 + mixedsum(c+1,b-1,a+1,proteinLen-2)*mixedsum(c,c,b+1,a-1);
            temp1=temp1+temp2*omega[b][a];
        }
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}


double GaussIntegral::int15_23_46(void)
{
    int a, b, c;
    double temp1, temp2;
    temp1=0;
    for(a=0; a < proteinLen-7; a++)
    {
        for(b=a+5; b < proteinLen-2; b++)
        {
            temp2 = 0;
            for(c=a+4; c < b; c++)
                temp2 = temp2 + partsum[a+1][c-1]*mixedsum(c,c,b+1,proteinLen-2);
            temp1 = temp1 + temp2*omega[a][b];
        }
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}



double GaussIntegral::int15_24_36(void)
{
    int a, b, c;
    double temp1, temp2;
    temp1=0;
    for(a = 1; a < proteinLen-5; a++)
    {
        for(b = a+2; b < proteinLen-3; b++)
        {
            temp2 = 0;
            for(c = b+2; c < proteinLen-1; c++)
                temp2 = temp2 + mixedsum(0,a-1,b+1,c-1)*mixedsum(a+1,b-1,c,c);
            temp1 = temp1 + temp2*omega[a][b];
        }
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}




double GaussIntegral::int15_26_34(void)
{
    int a, b, c;
    double temp1, temp2;
    temp1=0;
    for(a = 1; a < proteinLen-6; a++)
    {
        for(b = a+5; b < proteinLen-1; b++)
        {
            temp2 = 0;
            for(c = a+4; c < b; c++)
                temp2 = temp2 + partsum[a+1][c-1]*mixedsum(0,a-1,c,c);
            temp1 = temp1 + temp2*omega[a][b];
        }
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}



double GaussIntegral::int16_23_45(void)
{
    int a, b, c;
    double temp1, temp2;
    temp1=0;
    for(a=0; a < proteinLen-8; a++)
    {
        for(b=a+7; b < proteinLen -1; b++)
        {
            temp2 = 0;
            for(c = a+4; c < b-2; c++)
                temp2 = temp2 + partsum[a+1][c-1]*mixedsum(c,c,c+2,b-1);
            temp1 = temp1 + temp2*omega[a][b];
        }
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}


double GaussIntegral::int16_24_35(void)
{
    int a, b, c;
    double temp1, temp2;
    temp1 = 0;
    for(a=2; a<proteinLen-4; a++)
    {
        for(b=a+2; b < proteinLen-2; b++)
        {
            temp2=0;
            for(c = 0; c < a-1; c++)
                temp2 += mixedsum(c,c,b+1,proteinLen-2)*mixedsum(c+1,a-1,a+1,b-1);
            temp1+=temp2*omega[a][b];
        }
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}



double GaussIntegral::int16_25_34(void)
{
    int a, b;
    double temp1;
    temp1=0;
    for(a=1; a < proteinLen-6; a++)
        for(b = a+4; b < proteinLen-2; b++)
            temp1 = temp1 + mixedsum(0,a-1,b+1,proteinLen-2)*omega[a][b]*partsum[a+1][b-1];
    return temp1/(M_PI*M_PI*M_PI*8);
}


void GaussIntegral::GaussAll(double* gptr)
{

//    cout<<"protein length:"<<proteinLen<<endl;
//    cout<<"###################"<<endl;

//    for(int i=0; i< proteinLen; i++)
//        cout<<(*currentProtein)[i]<<endl;

//    cout<<"##########################"<<endl;
    createUnitVector();
   

//    unitvector.printArry();
    
    createomega();
//    omega.printArry();
//      absomega.printArry();
    createpartsum();

//    partsum.printArry();    

//    abspartsum.printArry();

    int offset;
    if(CylinderTransform)
    {
        gptr[0] = proteinLen;
        offset = 1;
    }
    else
        offset = 0;

    //"raw Guass Integral    
    gptr[offset+0] = int12(0);
    gptr[offset+1] = int12(1);
    //second order
    gptr[offset+2] = int12_34(0, 0);
    gptr[offset+3] = int12_34(1, 0);
    gptr[offset+4] = int12_34(0, 1);
    gptr[offset+5] = int12_34(1, 1);
    gptr[offset+6] = int13_24(0, 0);
    gptr[offset+7] = int13_24(1, 0);
    gptr[offset+8] = int13_24(0, 1);
    gptr[offset+9] = int13_24(1, 1);
    gptr[offset+10] = int14_23(0, 0);
    gptr[offset+11] = int14_23(1, 0);
    gptr[offset+12] = int14_23(0, 1);
    gptr[offset+13] = int14_23(1, 1);
    //third order
    gptr[offset+14] = int12_34_56();
    gptr[offset+15] = int12_35_46();
    gptr[offset+16] = int12_36_45();
   
    gptr[offset+17] = int13_24_56();
    gptr[offset+18] =  int13_25_46();
    gptr[offset+19] = int13_26_45();
    
    gptr[offset+20] = int14_23_56();
    gptr[offset+21] = int14_25_36();
    gptr[offset+22] = int14_26_35();
    
    gptr[offset+23] =  int15_23_46();
    gptr[offset+24] = int15_24_36();
    gptr[offset+25] = int15_26_34();
    
    gptr[offset+26] = int16_23_45();
    gptr[offset+27] = int16_24_35();
    gptr[offset+28] = int16_25_34();


    double pl = proteinLen;
    double pl2 = pl*pl;
    double pl3 = pl*pl2;

    // coefficient
    double coefficient[30]=
    {
        70.2135*1.339875
        ,pl*0.0401388*1.030110, pl*0.0730344*1.222925
        ,pl2*0.0014755*1.200300, pl2*0.00533629*1.050380, pl2*0.00513504*1.058926
        ,pl2*0.0157794*1.174822, pl2*5.7905e-05*1.462006, pl2*0.000460471*1.044327
        ,pl2*0.000485036*1.095845, pl2*0.00418937*1.282707, pl2*0.000188964*1.097088
        ,pl2*0.00269016*0.9816835, pl2*0.000892068*1.153978, pl2*0.00758828*1.241605
        ,pl3*3.7779e-05*1.394955, pl3*2.17081e-06*1.615541, pl3*2.98437e-06*1.098159
        ,pl3*2.21604e-06*1.592069, pl3*5.8425e-08*1.692098, pl3* 3.64209e-07*1.070658
        ,pl3*3.60442e-06*1.049276, pl3*6.07681e-08*1.313312, pl3*6.60069e-08*1.148898
        ,pl3*2.95267e-07*1.402297, pl3*5.84497e-08*1.369921, pl3*2.98807e-07*1.291011
        ,pl3*3.66585e-06*1.254290, pl3*2.25853e-07*1.359014, pl3*4.20434e-07*1.271251
    };
    if(CylinderTransform)
       for(int i = 0; i < 30; i++)
            gptr[i] /= coefficient[i];

}



