#include "Utility.hpp"
#include "Protein_3D_Index.hpp"


using namespace std;

class benchMark
{

    struct predict
    {
        string pname;
        double Descriptor[DESCRIPTOR_LENGTH];
        double l1_dist;
        double l2_dist;
        double RMSD;
        double GDT_TS;
        double TM_score;
    };

    struct target
    {
        string target_name;
        double Descriptor[DESCRIPTOR_LENGTH];
        int length;

        vector<predict> predictions;
    };

    public:
        benchMark():Target_Path(), Prediction_Path(), Result_Path()
        {}
        void DescriptorTest();


        void setParas(string tpath, string ppath, string rpath, string ddfile, string opath)
        {
            Target_Path = tpath;
            Prediction_Path = ppath;
            Result_Path = rpath;
            //domain_definition = ddfile;
            outputPath = opath;
        }

        // according to the domain_definition file 
        // trim all files in Target_Path to domain pdb files
        // for example T0759 to T0759-D1 T0759-D2 ...
        void split_to_domain(std::string domain_definition, std::string targetPath, std::string domainPath);

    private:
        string Target_Path;
        string Prediction_Path;
        string Result_Path;
        //string domain_definition;
        string outputPath;


};
