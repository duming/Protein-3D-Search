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

        inline bool readPDB(std::vector<POINT> &coords)
        {
            return readPDB(dataPath + domainPath + domainName, coords);
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

        static std::string dataPath;
        static std::string domainPath;


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

      
        // read cathList file and return the length of the list
        int readList();


        // pure for test fake a list of pdb file names by reading all files
        // in a directory
        int fakeList(std::string);

        void printData(int option);

        void test();

        Cathdomain& operator [](int idx)
        {
            return domains[idx];
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

    private:
        static char*inputBuff;

        std::vector<Cathdomain> domains;
        int domain_num;
        std::string listFileName;

        static std::string dataPath;
        static std::string domainPath;

};



#endif
