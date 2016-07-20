#ifndef CATHDATA_HPP
#define CATHDATA_HPP
#include "Utility.hpp"

//define the default length of how many hierarchical category will be used
#define CATE_LENGTH_DEFAULT 9
#define DESCRIPTOR_LENGTH 30
#define BUFFLENGTH 4096
#define LINELENGTH 128
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
        bool readPDB(std:: string fileName, std::vector<POINT> &coords);

        void printDomain();

    private:

        static char* inputBuff;

        //Calculated Descriptors such as Gauss Integral
        double descriptor[DESCRIPTOR_LENGTH];

        ////////////////////////
        //cath PDB file
        /////////////////////////
        //the coordinates of c-alpha atoms
        std::vector<double[3]>* coords;


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

       
        bool readList();

        void printData();

        void test();
    private:
        static char*inputBuff;

        std::vector<Cathdomain> domains;
        int domain_num;
        std::string listFileName;

        static std::string dataPath;
        static std::string domainPath;

};



#endif
