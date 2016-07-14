#ifndef UTILITY
#define UTILITY
#include <vector>
#include <string>
#include <iostream>
#include <sstream> 

class Utility
{
    Utility();
    ~Utility();

    static std::vector<std::string> split(std::string &s, char delim);
};



#endif
