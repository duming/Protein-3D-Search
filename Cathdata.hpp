#ifndef CATHDATA_HPP
#define CATHDATA_HPP
#include "Utility.hpp"

#include "GaussIntegral.hpp"

//define the default length of how many hierarchical category will be used
#define CATE_LENGTH_DEFAULT 9
#define DESCRIPTOR_LENGTH 30
#define BUFFLENGTH 4096
#define LINELENGTH 128

#define DATAPATH "data/"
#define DOMAINPATH "dompdb/"

class Cathdomain
{
    public:
        Cathdomain()
        {
            coords = NULL;
        }
        ~Cathdomain()
        {}

        //read information from PDB file
        //store all C-alpha coordinates to variable coords
        static bool readPDB( std:: string fileName, std::vector<POINT> &coords);

        inline bool readPDB(std::vector<POINT> &coords, std::string path)
        {
            return readPDB(path + domainName, coords);
        }


        std::string getName()
        {
            return domainName;
        }
    

        void printDomain(int option);
    public:
        //Calculated Descriptors such as Gauss Integral
        double descriptor[DESCRIPTOR_LENGTH];


    private:

        static char* inputBuff;

                ////////////////////////
        //cath PDB file
        /////////////////////////
        //the coordinates of c-alpha atoms
        std::vector<POINT>* coords;


        //////////////////
        //cath list file
        ////////////////////
        //the name of the domain
        std::string domainName;

        //cath hierarchical category
        int category[CATE_LENGTH_DEFAULT];

        //the length of the domain
        int domainLength;

        //the structure resolution (999.000 for NMR and 1000 for obsolete PDB entries)
        double resolution;

        //static std::string dataPath;
        //static std::string domainPath;


        friend class CathData;
};


class CathData
{
    public:
        //need the Cath list file to initialize
        CathData(std::string fileName)
        {
            listFileName = fileName;
        }
        ~CathData()
        {}

        void setPath(std::string data_Path,std::string domain_Path)
        {
            dataPath = data_Path;
            domainPath = domain_Path;
        }
      
        // read cathList file and return the length of the list
        int readList();


        // pure for test fake a list of pdb file names by reading all files
        // in a directory
        int fakeList(std::string);


        // save the CathData to a cathList File
        int saveList(std::string fileName);

        void printData(int option);

        void test();

        Cathdomain& operator [](int idx)
        {
            return domains[idx];
        }


        bool readPDB(int i, std::vector<POINT> &coords)
        {
            return domains[i].readPDB(dataPath + domainPath + domains[i].domainName, coords);
        }

        void printDomains(int start = 0, int end = -1, int option= 1|2)
        {
            if(end == -1)
                end = domain_num -1;
            for(int i=start ; i<=end;i++)
                domains[i].printDomain(option);
        }

        int size()
        {
            return domain_num;
        }


        bool operator ==(const CathData & op2) const
        {
            //check length
            if(domain_num != op2.domain_num)
            {
                std::cout<<"unequal length\n";
            }
            for(int i = 0; i < domain_num; i++)
            {
                for(int j = 0; j < CATE_LENGTH_DEFAULT; j++)
                   if(domains[i].category[j] != op2.domains[i].category[j])
                    {
                        std::cout<<"number :"<<i<<" line category not equal\n";
                        return false;
                    }
                if(domains[i].domainLength != op2.domains[i].domainLength)
                {
                    std::cout<<"number :"<<i<<" line domain length not equal\n";
                    return false;
                }
                if(domains[i].resolution != op2.domains[i].resolution)
                {
                    std::cout<<"number :"<<i<<" line resolution not equal\n";
                    return false;
                }
            }
            return true;
        }

        // save all descriptors to binary file in order to 
        // build SLH index
        // format:
        // usigned int : sizeof(datatype)
        // unsigned int : "row number" one row for one descriptor
        // unsigned int : "descriptor length"
        // (row number) * (descriptor length) datatype
        void saveDescriptor(std::string fileName);


        void saveDescriptorText(std::string fileName);

    private:
        static char*inputBuff;

        std::vector<Cathdomain> domains;
        int domain_num;
        std::string listFileName;

        std::string dataPath;
        std::string domainPath;

};



#endif
