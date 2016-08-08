#include "Utility.hpp"

std::ostream& operator <<(std::ostream& os, point & pt)
{
    os<<pt.x<<" "<<pt.y<<" "<<pt.z;
    return os;
}


std::vector<std::string>
Utility:: split(std::string &s, char delim)
{
    std::vector<std::string> tokens;
    split(tokens, s, delim);
    return tokens;
}    



void 
Utility:: split(std::vector<std::string> &tokens, std::string &s, char delim )
{
    tokens.resize(0);
    std::stringstream ss(s); 
    std::string item;
    while (getline(ss, item, delim)) 
    {
        if(delim == ' ' && item.empty())
            continue;
        tokens.push_back(item);
    }
}




void
Utility:: split(std::vector<std::string> &tokens, const char* s, char delim)
{
    tokens.resize(0);
    do
    {
        if(' ' ==delim)
            while( ' ' == *s)
                s++;
        const char*begin = s;
        while(*s != delim && *s)
            s++;
        tokens.push_back(std::string(begin, s));
    }while(0 != *s++);

}



bool 
Utility:: listFiles(std::string dir, std::vector<std::string> &files, bool isFile)
{
	DIR *dp;
	struct dirent *dirp;
	if((dp  = opendir(dir.c_str())) == NULL)
    {
        std::cout << "Can't open " << dir << std::endl;
	    return false;
	}
	
	while ((dirp = readdir(dp)) != NULL) 
    {
        if((isFile && dirp->d_type == DT_REG)
            || (!isFile && dirp->d_type == DT_DIR)    
           )
	        files.push_back(std::string(dirp->d_name));
	}
	closedir(dp);
	return true;
}



////////////////////////////////////
//    vector operation
///////////////////////////////////


double Utility::vectLen(POINT const &op)
{
   return sqrt(vectDot(op,op)); 
}


double Utility::vectDot(POINT const &op1, POINT const &op2)
{
    return  op1.x*op2.x + op1.y*op2.y + op1.z*op2.z;
}

void Utility::vectUnit(POINT &result, point const & op)
{   
    double len = vectLen(op);
    result.x = op.x/len;
    result.y = op.y/len;
    result.z = op.z/len;
}


void Utility::vectSub(POINT &result, POINT const &op1, POINT const &op2)
{
    result.x = op2.x - op1.x;
    result.y = op2.y - op1.y;
    result.z = op2.z - op1.z;
}


void Utility::vectCross(POINT &result, POINT const &op1, POINT const &op2)
{
    result.x = op1.y * op2.z - op2.y * op1.z;
    result.y = op1.z * op2.x - op2.z * op1.x;
    result.z = op1.x * op2.y - op2.x * op1.y;
}



double Utility::round_to_digits(double value, int digits)
{
    if (value == 0.0) // otherwise it will return 'nan' due to the log10() of zero
        return 0.0;

    double factor = pow(10.0, digits - ceil(log10(fabs(value))));
    return round(value * factor) / factor;   
}


void Utility:: readGIs(std::string fileName, std::vector<std::vector<double> > &result, bool isTrans)
{
    std::ifstream infile(fileName);
    if(!infile.is_open())
    {
        std::cout<<"error opening:"<<fileName<<std::endl;
        return;
    }
    
    int ds_len;
    if(isTrans)
        ds_len = 30;
    else
        ds_len = 29;

    std::vector<double> tmpvct;
    tmpvct.resize(ds_len);
    std::string tmpstr;
   
    double foo;

    while(true)
    {
        if(!(infile>>tmpstr))
            break;
        if(!(infile>>tmpstr))
            break;
        if(!(infile>>tmpstr))
            break;
        if(!(infile>>tmpstr))
            break;
 

   
        int i;
        for(i=0;i<ds_len;i++)
            if(!(infile>>tmpvct[i]))
                break;

        if(i==ds_len)
            result.push_back(tmpvct);
        else
            break;
    }

}
