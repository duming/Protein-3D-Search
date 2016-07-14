#include "Utility.hpp"


std::vector<std::string>
Utility:: split(std::string &s, char delim)
{
    std::stringstream ss(s); 
    std::string item;
    std::vector<std::string> tokens;
    while (getline(ss, item, delim)) 
    {
        tokens.push_back(item);
    }
    return tokens;
}
