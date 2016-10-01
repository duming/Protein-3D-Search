#include "Utility.hpp"
#include "Protein_3D_Index.hpp"


using namespace std;

class benchMark
{

    struct predict
    {
        //string pname;
        int length;
        double Descriptor[DESCRIPTOR_LENGTH];
        double l1_dist;
        double l2_dist;
        double RMSD;
        double GDT_TS;
        double LGA_S3;
        double LGA_Q;
    };

    struct target
    {
        string target_name;
        double Descriptor[DESCRIPTOR_LENGTH];
        int length;

        map<string, predict > predictions;
    };

    public:
        benchMark():Target_Path(), Prediction_Path(), Result_Path()
        {}

        void DescriptorTest();


        // calculate RMSD between one PDB file and all PDB files under the path
        static void getRMSD(std::string fileName, std::string path, 
                std::vector<std::pair<std::string,double> > & result);


        // calculate rmsd between all <target,prediction> pairs
        // and save them into .SUMMARY.lga_sda.txt file
        void my_LGA_result();


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



        // split pdb to equal length segments
        // and save result in same directory

        void split_to_equal(std::string targetPath, std::string splitPath, int segLen);

        // prepare directories for pdb split
        void prepare_directory(std::string targetPath, std::string dirPath, int segLen);

        // split pdb to equal length segments 
        // and save them to different directories
        void split_to_equal_diff(std::string targetPath, std::string splitPath, int segLen);


        void split_all_prediction(std::string predictPath, std::string splitPath, int segLen);


        void pdb_range_copy(std::vector<std::string>& from
                        , std::vector<std::string>& to
                        , std::vector<std::pair<int,int> > ranges);

       
        // read all files includes:
        // 1. target.list target.descriptor
        // 2. outputPath/*.list outputPath/*.descriptor
        // 3. result_Path/*.txt files that contains all kinds of distances
        void readResult();



        // calculate l1_dist and l2_dist for every prediction
        void calDist();

        // check all the data that read from file is validate
        bool DataValidation();



        bool saveResult(std::string path);
    private:
        
    private:
        string Target_Path;
        string Prediction_Path;
        string Result_Path; //casp result file
        //string domain_definition;
        string outputPath;
        std::vector<target> result;

};
