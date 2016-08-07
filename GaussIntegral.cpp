#include "GaussIntegral.hpp"

#include <time.h>
#include <iomanip>
using namespace std;

double GaussIntegral::Epsilon = 1e-16;

/*
template<class T>
void GaussIntegral::SquareArray<T>::printArry()
{
    for(int i=0; i < size; i++)
    {
        for(int j = 0; j < size; j++)
            cout<<base[i*size + j]<<end;
        cout<<endl;
    }
}
*/

void GaussIntegral::test()
{
    //array test
    /*
    int size = 10;
    SquareArray<int> sa(20);

    sa.resize(size);
    for(int i=0;i<size;i++)
        for(int j=0;j<size;j++)
            sa(i,j) = i*size +j;


    sa.printArry();

    sa.resize(100);
    sa(99,99) = 100;
    cout<<sa(99,99)<<endl;
    */
    string path = "data/UnitTestpdb/";
    vector<string> files;
    Utility::listFiles(path,files);

    vector<vector<double> > GIresult;

    vector<POINT> *crdPtr = new vector<POINT>[files.size()];
    

    for(int i=0;i<files.size();i++)
    {
        crdPtr[i].resize(1000);
    }


    Utility::readGIs("data/GIout.txt", GIresult);

    
    double dscrpt[files.size()][29]; 
    
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
    cout<<t2 -t1<<endl;

    for(int i; i < GIresult.size(); i++)
    {
        for(int j = 0; j < 29 ;j++)
        {
            cout<<j<<"\t:\t"<<setw(20)<<GIresult[i][j]<<"\t"<<setw(10)<<dscrpt[i][j]<<"\t\t";

            cout<<dscrpt[i][j] - GIresult[i][j]<<endl;
        }

        cout<<endl;
    }
    

}



void GaussIntegral:: createUnitVector()
{
    POINT tmp;
    for(int i = 0; i < proteinLen-1 ; i++)
        for(int j = i+1 ; j < proteinLen ; j++)
        {
            Utility::vectSub(tmp, (*currentProtein)[i], (*currentProtein)[j]);
            Utility::vectUnit(unitvector(i,j),tmp);
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
        omega(i,i+1) = 0;
        absomega(i,i+1) = 0;
    }

    for(int i=0;i<proteinLen-2;i++)
    {
        for(int j = i+2; j<proteinLen-1 ; j++)
        {
            Utility::vectCross(cross1, unitvector(i,j), unitvector(i,j+1));
            Utility::vectCross(cross2, unitvector(i+1,j), unitvector(i+1,j+1));
            Utility::vectCross(cross3, cross1, cross2);
            if(Utility::vectDot(cross3,cross3)<Epsilon)
            {
                omega(i,j) = 0;
                absomega(i,j) = 0;
            }
            else
            {
                anglesum =  -spangle(unitvector(i,j), unitvector(i,j+1), unitvector(i+1,j+1))
                            -spangle(unitvector(i,j+1), unitvector(i+1,j+1), unitvector(i+1,j))
                            -spangle(unitvector(i+1,j+1), unitvector(i+1,j), unitvector(i,j))
                            -spangle(unitvector(i+1,j), unitvector(i,j), unitvector(i,j+1));

                if(anglesum >=0)
                    omega(i,j) = 2*M_PI - anglesum;
                else
                    omega(i,j) = -2*M_PI - anglesum;

                absomega(i,j) = fabs(omega(i,j));
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
            partsum(i,i+k) = 0;

    for(i = 0; i<proteinLen -2-1;i++)
        partsum(i,i+2) = omega(i,i+2);

    for(k = 3; k < proteinLen -1; k++)
        for(i = 0; i<proteinLen - k -1; i++)
            partsum(i,i+k) = partsum(i,i+k-1) + partsum(i+1,i+k) - partsum(i+1,i+k-1)
                            + omega(i,i+k);
    
    //////////
    //absolute partsum
    ////////
    for(k = 0; k < 2; k++)
        for(i = 0 ;i< proteinLen-k ;i++)
            abspartsum(i,i+k) = 0;

    for(i = 0; i<proteinLen -2-1;i++)
        abspartsum(i,i+2) = absomega(i,i+2);

    for(k = 3; k < proteinLen -1; k++)
        for(i = 0; i<proteinLen - k-1; i++)
            abspartsum(i,i+k) = abspartsum(i,i+k-1) + abspartsum(i+1,i+k) - abspartsum(i+1,i+k-1)
                            + absomega(i,i+k);

}


double GaussIntegral:: mixedsum(int a, int b, int c, int d, int isabs)
{
    double tempm, tempp;
    if(isabs != 0)
    {
        tempm = abspartsum(a,c-1) + abspartsum(b+1,d);
        tempp = abspartsum(a,d) + abspartsum(b+1,c-1);
    
    }
    else
    {
        tempm = partsum(a,c-1) + partsum(b+1,d);
        tempp = partsum(a,d) + partsum(b+1,c-1);
    }
    return tempp-tempm;
}


double GaussIntegral::int12(int isabs)
{
    if(isabs)
        return abspartsum(0,proteinLen-2)/(M_PI*2);
    else
        return partsum(0,proteinLen-2)/(M_PI*2);
}



double GaussIntegral::int12_34(int isabs1, int isabs2)
{
    int b;
    double temp;
    SquareArray<double>* arry;
    if(isabs1 == 1)
        arry = &abspartsum;
    else
        arry = &partsum;
    temp=0;
    for(b=3;b<proteinLen-3; b++)
        temp=temp+(*arry)(0,b-1)*mixedsum(b,b,b+2,proteinLen-2, isabs2);
    return temp/(M_PI*M_PI*4);
}


double GaussIntegral::int13_24(int isabs1, int isabs2)
{
    int a, b;
    double temp;
    SquareArray<double>* arry;
    if(isabs2 == 1)
        arry = & absomega;
    else
        arry = &omega;
    temp=0;
    for(a=0; a<proteinLen-3; a++)
        for(b=a+2; b<proteinLen -1; b++)
            temp = temp + mixedsum(0,a-1,a+1,b-1, isabs1)*(*arry)(a,b);
    return temp/(M_PI*M_PI*4);
}




double GaussIntegral::int14_23(int isabs1, int isabs2)
{
    int a, b;
    double temp;
    SquareArray<double> *s_arry,*o_arry;
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
            temp=temp+(*s_arry)(a+1,b-1)*(*o_arry)(a,b);
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
        temp1=temp1+temp2*partsum(0,a-1);
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
            temp2 = temp2 + omega(a,b)*mixedsum(a+1,b-1,b+1,proteinLen-2);
        temp1 = temp1 + temp2*partsum(0,a-1);
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
            temp2 = temp2 + omega(a,b)*partsum(a+1,b-1);
        temp1=temp1+temp2*partsum(0,a-1);
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}


double GaussIntegral::int13_24_56(void)
{
    int a, b;
    double temp1, temp2;
    temp1=0;
    for(a=4; a < proteinLen-3;a++)
    {
        temp2=0;
        for(b=1; b < a-1; b++)
            temp2 = temp2 + omega(b,a)*mixedsum(0,b-1,b+1,a-1);
        temp1=temp1+temp2*partsum(a+1,proteinLen-2);
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
            for(c=b+1; c < proteinLen-1; c++)
                temp2 = temp2 + omega(a,c)*mixedsum(b,b,c+1,proteinLen-2);
            temp1 = temp1 + temp2*mixedsum(1,a-1,a+1,b-1);
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
                temp2 = temp2 + partsum(c+1,b-1)*mixedsum(0,a-1,c,c);
            temp1 = temp1 + temp2*omega(a,b);
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
            temp2 = temp2 + omega(b,a)*partsum(b+1,a-1);
        temp1 = temp1 + temp2*partsum(a+1,proteinLen-2);
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}



double GaussIntegral::int14_25_36(void)
{
    int a, b, c;
    double temp1, temp2;
    temp1=0;
    for(a=2; a < proteinLen-5; a++)
    {
        for(b=a+3; b < proteinLen-2; b++)
        {
            temp2=0;
            for(c=a+1; c < b-1; c++)
                temp2 = temp2 + mixedsum(1,a-1,c+1,b-1)*mixedsum(c,c,b+1,proteinLen-2);
            temp1 = temp1 + temp2*omega(a,b);
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
            temp1=temp1+temp2*omega(b,a);
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
                temp2 = temp2 + partsum(a+1,c-1)*mixedsum(c,c,b+1,proteinLen-2);
            temp1 = temp1 + temp2*omega(a,b);
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
            temp1 = temp1 + temp2*omega(a,b);
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
                temp2 = temp2 + partsum(a+1,c-1)*mixedsum(1,a-1,c,c);
            temp1 = temp1 + temp2*omega(a,b);
        }
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}



double GaussIntegral::int16_23_45(void)
{
    int a, b, c;
    double temp1, temp2;
    temp1=0;
    for(a=1; a < proteinLen-8; a++)
    {
        for(b=a+7; b < proteinLen -1; b++)
        {
            temp2 = 0;
            for(c = a+4; c < b-2; c++)
                temp2 = temp2 + partsum(a+1,c-1)*mixedsum(c,c,c+2,b-1);
            temp1 = temp1 + temp2*omega(a,b);
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
            temp1+=temp2*omega(a,b);
        }
    }
    return temp1/(M_PI*M_PI*M_PI*8);
}



double GaussIntegral::int16_25_34(void)
{
    int a, b;
    double temp1;
    temp1=0;
    for(a=2; a < proteinLen-6; a++)
        for(b = a+4; b < proteinLen-2; b++)
            temp1 = temp1 + mixedsum(0,a-1,b+1,proteinLen-2)*omega(a,b)*partsum(a+1,b-1);
    return temp1/(M_PI*M_PI*M_PI*8);
}




void GaussIntegral::GaussAll(double* gptr)
{
/*
    cout<<"protein length:"<<proteinLen<<endl;
    cout<<"###################"<<endl;


    for(int i=0; i< proteinLen; i++)
        cout<<(*currentProtein)[i]<<endl;

    cout<<"##########################"<<endl;
*/
    createUnitVector();
   

  //  unitvector.printArry();
    
    createomega();
//    omega.printArry();
//      absomega.printArry();
    createpartsum();

    //partsum.printArry();    

//    abspartsum.printArry();

    
    gptr[0] = int12(0);
    gptr[1] = int12(1);
    //second order
    gptr[2] = int12_34(0, 0);
    gptr[3] = int12_34(1, 0);
    gptr[4] = int12_34(0, 1);
    gptr[5] = int12_34(1, 1);
    gptr[6] = int13_24(0, 0);
    gptr[7] = int13_24(1, 0);
    gptr[8] = int13_24(0, 1);
    gptr[9] = int13_24(1, 1);
    gptr[10] = int14_23(0, 0);
    gptr[11] = int14_23(1, 0);
    gptr[12] = int14_23(0, 1);
    gptr[13] = int14_23(1, 1);
    //third order
    gptr[14] = int12_34_56();
    gptr[15] = int12_35_46();
    gptr[16] = int12_36_45();
   
    gptr[17] = int13_24_56();
    gptr[18] =  int13_25_46();
    gptr[19] = int13_26_45();
    
    gptr[20] = int14_23_56();
    gptr[21] = int14_25_36();
    gptr[22] = int14_26_35();
    
    gptr[23] =  int15_23_46();
    gptr[24] = int15_24_36();
    gptr[25] = int15_26_34();
    
    gptr[26] = int16_23_45();
    gptr[27] = int16_24_35();
    gptr[29] = int16_25_34();

}




