#include "Utility.hpp"
#include <iostream>
#include <fstream>
using namespace std;

int main()
{
    string ss("123  333 44 55");

    vector<string> files;

    Utility:: listFiles(string("./") , files, false);

    for(int i=0; i<files.size(); i++)
        cout<<files[i]<<endl;
       

    return 0;
}
