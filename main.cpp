#include "Utility.hpp"
#include <iostream>
#include <fstream>
#include "CathData.hpp"
#include <time.h>
using namespace std;


int main()
{
    string ss("123  333 44 55");

    const char *s1 = "123456";
    char s4[] = "123 : 333 44 5";

    clock_t t1,t2;
    vector<string> files;

   CathData cc("CathDomainList.v4.0.0.txt");

   t1 = clock();
   cc.readList();
   t2 = clock();

    cc.test();

   //cout<<t2-t1<<endl;
/*    
    vector<string> tks;
    Utility::split(tks, s4, ':');
    for(int i=0; i<tks.size(); i++)
        cout<<tks[i]<<endl;
*/


    return 0;
}
