#include "Utility.hpp"

using namespace std;


bool POINT::isTooFar(POINT & p1, POINT & p2)
{   

    if(dist(p1,p2) > 6.0)
        return true;
    return false;
}


 double POINT::dist(POINT &p1, POINT & p2)
{
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}


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



bool  Utility:: readPDB(std::string fileName, std::vector<std::string> & data, bool ca_only)
{
    std::ifstream infile(fileName);
    if(!infile.is_open())
    {
        std::cout<<"error opening file: "<<fileName<<std::endl;
        return false;
    }
    //setbuffer
    int LINELENGTH = 500;

    //read the file
    char line[LINELENGTH];
    char name[5];
    data.resize(0);

    while(infile.getline(line, LINELENGTH))
    {
        //Utility::split(tokens, line);
        sscanf(line,"%4c",name);
        
        if( strcmp(name, "ATOM") )
            continue;

        if(ca_only && !(line[13] =='C' && line[14] =='A'))
            continue;
      
        data.push_back(std::string(line));   
    }
    infile.close();
    return true;
}



int Utility::next_Gap(vector<POINT> & cords)
{
    int i;
    for(i = 0 ; i < cords.size() -1; i++)
        if(POINT::isTooFar(cords[i], cords[i+1]))
                break;
    if(i < cords.size() -1)
        return i;
    
    return -1;
}




bool  Utility:: readFile(std::string fileName, std::vector<std::string> & data)
{
    std::ifstream infile(fileName);
    if(!infile.is_open())
    {
        std::cout<<"error opening file: "<<fileName<<std::endl;
        return false;
    }
    //setbuffer
    int LINELENGTH = 500;

    //read the file
    char line[LINELENGTH];
    char name[5];
    data.resize(0);

    while(infile.getline(line, LINELENGTH))
    {
        data.push_back(std::string(line));   
    }
    infile.close();
    return true;
}



bool Utility::writeFile(string fileName, vector<string> & data)
{
    std::ofstream outfile(fileName);
    if(!outfile.is_open())
    {
        std::cout<<"error writing file: "<<fileName<<std::endl;
        return false;
    }

    for(int i=0; i < data.size(); i++)
        outfile<<data[i]<<'\n';
    outfile.close();
    return true;
}




int Utility::firstRes(vector<string>& dataPDB)
{
    if(dataPDB.size() == 0)
        return 0;
    int resnum;
    sscanf(dataPDB.front().c_str()+22, "%4i", &resnum);
    return resnum;
}

int Utility::lastRes(vector<string>& dataPDB)
{
    if(dataPDB.size() == 0)
        return 0;

    int resnum;
    sscanf(dataPDB.back().c_str()+22, "%4i", &resnum);
    return resnum;
}





void Utility::splitPDB(vector<string> & dataPDB, vector<vector<string> > & segs, int segLen)
{
    int segNum = 0;
    int curRes, lastRes;
    lastRes = Utility::firstRes(dataPDB);
    segs.resize(0);
    segs.resize(1);
    int resCount = 0;
    for(int i=0; i < dataPDB.size(); i++)
    {
        sscanf(dataPDB[i].c_str()+22,"%4i",&curRes); 
        if(curRes != lastRes)
        {
            resCount++;
            if(resCount >= segLen)
            {
                resCount = 0;
                segs.push_back(vector<string>());

            }
            lastRes = curRes;
        }
        segs.back().push_back(dataPDB[i]);
        
    }
    if(resCount < segLen )
        segs.pop_back();

}



void Utility:: writePDB(std::string fileName, std::vector<std::string> & data, int start, int end)
{
    std::ofstream outfile(fileName);
    int resSeq;
    if(end == -1)
        end = data.size();
    for(int i = 0; i < data.size(); i++)
    {
        sscanf(data[i].c_str()+22,"%4i",&resSeq);
        if(resSeq >= start && resSeq <=end)
            outfile<<data[i]<<'\n';
    }
    outfile.close();
}


std::string Utility:: removeSuff(const std::string &str)
{
    return str.substr(0, str.find_last_of("."));
}


bool Utility:: readTable(string fileName, vector<vector<string> >&dists, vector<string> colum_names)
{
    vector<string> data_str;
    if(!Utility::readFile(fileName, data_str))
        return false;
    int Pos[colum_names.size()];

    // search file header
    vector<string> header = Utility::split(data_str[0],' ');  
    for(int i=0; i < colum_names.size(); i++)
    {
        Pos[i] = find(header.begin(), header.end(), colum_names[i]) - header.begin();
        if(Pos[i] >= header.size())
        {
            cout<<"readDist: can't find: "<<colum_names[i]<<"in header"<<endl;
            return false;
        }
        //cout<<Pos[i]<<' ';
    }
    //cout<<endl;

    //fetch colum we want
    vector<string> row,row_data;
    dists.resize(0);
    for(int i=1; i < data_str.size(); i++)
    {
        if(data_str[i].length() == 0)
            continue;
        row = Utility::split(data_str[i], ' ');
        row_data.resize(0);
        for(int j = 0; j < colum_names.size(); j++)
            row_data.push_back(row[Pos[j]]);
        dists.push_back(row_data); 
    }
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
    if(len == 0)
    {//if the input is a zero vector return itself
        result = op;
        return;
    }
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
