#include "Cathdata.hpp"
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace std;
string CathData::dataPath = DATAPATH;
string CathData::domainPath = DOMAINPATH;
string Cathdomain::dataPath = DATAPATH;
string Cathdomain::domainPath = DOMAINPATH;

char* Cathdomain:: inputBuff = new char[BUFFLENGTH];
char* CathData:: inputBuff = new char[BUFFLENGTH];


void Cathdomain:: printDomain()
{
    cout<<domainName<<'\t';
    for(int i=0; i< CATE_LENGTH_DEFAULT; i++)
        cout<<category[i]<<'\t';
    cout<<domainLength<<'\t';
    cout<<resolution<<endl;
}




bool Cathdomain:: readPDB(string fileName, vector<POINT> &coords)
{
    coords.resize(0);
    ifstream infile(fileName);
    if(!infile.is_open())
    {
        cout<<"error opening file: "<<fileName<<endl;
        return false;
    }
    //setbuffer
    infile.rdbuf()->pubsetbuf(inputBuff, BUFFLENGTH);

    //read the file
    char line[LINELENGTH];
    ///vector<string> tokens;
    char lastChain[5] = "none", newChain[5], name[5], atom[3]; 
    int resSeq, lastResSeq = -1;
    POINT tempP;
    while(infile.getline(line, LINELENGTH))
    {
        //Utility::split(tokens, line);
        sscanf(line,"%4c",name);
        sscanf(line+22,"%4i",&resSeq);
        
        if( strcmp(name, "ATOM") || line[13]!= 'C' || line[14] != 'A'||resSeq == lastResSeq)
            continue;
      
       
        lastResSeq = resSeq;

        sscanf(line+30,"%lf %lf %lf",&tempP.x,&tempP.y,&tempP.z);
        coords.push_back(tempP);
    }

    return true;
}

bool CathData:: readList()
{
    ifstream infile(dataPath + listFileName);
    if(!infile.is_open())
    {
        cout<<"error opening file: "<<listFileName<<endl;
        return false;
    }

    //set buffer
   // infile.rdbuf()->pubsetbuf(inputBuff, BUFFLENGTH);

    //read the file
    char line[LINELENGTH];
    vector<string> tokens;
    Cathdomain tempd;
    while(infile.getline(line,LINELENGTH))
    {
        if(line[0] == '#')
            continue;
        Utility::split(tokens, line);
        if(tokens.size() != 12)
        {
            cout<<"CLF format error"<<endl;
            return false;
        }
        
        //read domain name
        tempd.domainName = tokens[0];
        //read classes
        for(int i=0; i<CATE_LENGTH_DEFAULT; i++)
            tempd.category[i] = atoi(tokens[i+1].c_str());

        //read domain length
        tempd.domainLength = atoi(tokens[CATE_LENGTH_DEFAULT +1].c_str());
        //read resolution
        tempd.resolution = atof(tokens[CATE_LENGTH_DEFAULT+2].c_str());

        domains.push_back(tempd);
    }
    domain_num = domains.size();

    return true;
}


void CathData::printData()
{
    for(int i=0; i< domain_num ; i++)
    {
        cout<<i<<"\t:\t";
        domains[i].printDomain();
    }
}




void CathData:: test()
{
    readList();
    string fileName = dataPath + domainPath + domains[0].domainName;
    cout<<fileName<<endl;
    vector<POINT> temp;
    domains[0].readPDB(fileName, temp);
    

    string filename = "data/UnitTestpdb/1a0hA02";
    Cathdomain::readPDB(filename, temp);

    for(int i=0; i< temp.size(); i++)
        cout<<i<<':'<<temp[i].x<<"\t\t"<<temp[i].y<<"\t\t"<<temp[i].z<<endl;

}
