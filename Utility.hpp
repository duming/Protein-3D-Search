#ifndef UTILITY
#define UTILITY
#include <vector>
#include <string>
#include <iostream>
#include <sstream> 
#include <dirent.h>

class Utility
{
    public:
    Utility();
    ~Utility();
    
    //standard version with memory allocation
    static std::vector<std::string> split(std::string &s, char delim =' ');

    //fast version of split store result in tokens
    static void split(std::vector<std::string> &tokens, std:: string &s, char delim = ' ') ;

    //if isFile == True return all the files in current directory through file
    //otherwise return all sub-directories
    static bool listFiles(std::string directory, std::vector<std::string> & files, bool isFile = true);
};



#endif
