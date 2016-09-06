#include "Cathdata.hpp"
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iomanip>
using namespace std;
//string CathData::dataPath = DATAPATH;
//string CathData::domainPath = DOMAINPATH;
//string Cathdomain::dataPath = DATAPATH;
//string Cathdomain::domainPath = DOMAINPATH;

char* Cathdomain:: inputBuff = new char[BUFFLENGTH];
char* CathData:: inputBuff = new char[BUFFLENGTH];


void Cathdomain:: printDomain(int option)
{
    if(option&1)
    {
        //print name
        cout<<domainName<<'\t';
    }

    if(option&2)
    {
        //print category
        for(int i=0; i< CATE_LENGTH_DEFAULT; i++)
            cout<<category[i]<<'\t';
        cout<<domainLength<<'\t';
        cout<<resolution<<endl;
    }

    if(option&4)
    {
        //print descriptor
        int num_per_line = 6;
        for(int i=0; i < DESCRIPTOR_LENGTH; i++)
        {
            cout<<setw(10)<<descriptor[i]<<'\t';
            if(i%num_per_line == num_per_line-1)
                cout<<endl;
        }
        cout<<endl;
    }
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
            
        name[4] = '\0';
        if( strcmp(name, "ATOM") || line[13]!= 'C' || line[14] != 'A'||resSeq == lastResSeq)
            continue;
      
       
        lastResSeq = resSeq;

        sscanf(line+30,"%lf %lf %lf",&tempP.x,&tempP.y,&tempP.z);
        coords.push_back(tempP);
    }

    infile.close();
    return true;
}

int CathData:: readList()
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
            infile.close();
            return 0;
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

    infile.close();
    return domain_num;
}


int CathData::saveList(string fileName)
{
    ofstream outfile(fileName);
    if(!outfile.is_open())
    {
        cout<<"error opening file: "<<listFileName<<endl;
        return false;
    }

    string whiteSpaces(4,' ');
    for(int i = 0; i < domain_num; i++)
    {
        //the name
        outfile<<domains[i].domainName<<whiteSpaces;
        //category numbers
        for(int j = 0; j < CATE_LENGTH_DEFAULT; j++)
            outfile<<domains[i].category[j]<<whiteSpaces;
        //domain length
        outfile<<domains[i].domainLength<<whiteSpaces;
        //resolution
        outfile<<domains[i].resolution;
        outfile<<endl;
    }

    outfile.close();

    return 1;
}




int CathData::fakeList(string directory)
{
    vector<string> fileList;
    Utility::listFiles(directory, fileList);
    Cathdomain tempd;
    for(int i = 0; i < fileList.size(); i++)
    {
        tempd.domainName = fileList[i];
        domains.push_back(tempd);
    }
    domain_num = fileList.size();
    return fileList.size();
}


void CathData::printData(int option)
{
    for(int i=0; i< domain_num ; i++)
    {
        cout<<i<<"\t:\t";
        domains[i].printDomain(option);
    }
}



void CathData::saveDescriptor(string fileName)
{
    ofstream outfile(fileName, ofstream::binary);
    unsigned dataheader[3];
    //output the size of data type
    dataheader[0] = (unsigned)sizeof(double);
    //output number of number and dimension of descriptor
    dataheader[1] = (unsigned)domains.size();
    dataheader[2] = (unsigned)DESCRIPTOR_LENGTH;
    outfile.write((char*)dataheader,sizeof dataheader);

    //output the actual data
    double *buffer = new double[domains.size() * DESCRIPTOR_LENGTH];
    for( int i=0; i < domains.size(); i++)
        for( int j=0; j < DESCRIPTOR_LENGTH; j++)
            buffer[i*DESCRIPTOR_LENGTH + j] = domains[i].descriptor[j];

    outfile.write((char*)buffer, sizeof(double)* domains.size() * DESCRIPTOR_LENGTH);
    //free memory
    delete [] buffer;
    //close file
    outfile.close();

}



void CathData:: test()
{
    /*
    readList();
    string fileName = dataPath + domainPath + domains[0].domainName;
    cout<<fileName<<endl;
    vector<POINT> temp;
    domains[0].readPDB(fileName, temp);
    

    string filename = "data/UnitTestpdb/1a0hA02";
    Cathdomain::readPDB(filename, temp);

    for(int i=0; i< temp.size(); i++)
        cout<<i<<':'<<temp[i].x<<"\t\t"<<temp[i].y<<"\t\t"<<temp[i].z<<endl;
    */
    string fileName =  "test.cathlist";
    readList();
    saveList(dataPath + fileName);
    CathData c1(fileName);
    c1.readList();
    cout<<(c1 == *this)<<endl;

}
