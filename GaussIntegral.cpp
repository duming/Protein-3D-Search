#include "GaussIntegral.hpp"

using namespace std;

void GaussIntegral::test()
{
    int size = 10;
    SquareArray<int> sa(10);
    for(int i=0;i<size;i++)
        for(int j=0;j<size;j++)
            sa(i,j) = i*size +j;

    for(int i=0;i<size;i++)
    {
        for(int j=-1;j<size;j++)
            cout<<sa(i,j)<<"\t";
        cout<<endl;
    }

    sa.resize(100);
    sa(99,99) = 100;
    cout<<sa(99,99)<<endl;
}
