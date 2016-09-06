#include "benchMark.hpp"
#include <fstream>
#include <iostream>
#include "Cathdata.hpp"

using namespace std;


void benchMark:: DescriptorTest()
{
    string TargetList = outputPath + "Target.list";
    string TargetDesc = outputPath + "Target.descriptor";
    //calculate descriptors for targets
    {
        Protein_3D_Index p3d;
        p3d.desCalculate(Target_Path, outputPath + "Target.descriptor", outputPath + "Target.list");
    }
    CathData targets(TargetList);
    targets.readList();
    //targets.printData(1|2);
    
    //calculate descriptor for each set of prediction
    for(int i = 0; i < targets.size(); i++)
    {
        //domainName in cathdata may contain a suffix
        string TargetName = targets[i].getName();
        TargetName = Utility::split( TargetName, '.')[0] +"\/";
        
        Protein_3D_Index p3d;
        p3d.desCalculate(Prediction_Path + TargetName
                        , outputPath + TargetName + ".descriptor"
                        , outputPath + TargetName + ".list"
                );


    }
}



void benchMark::
split_to_domain(std::string domain_definition, std::string targetPath, std::string domainPath)
{
    ifstream infile(domain_definition);
    if(!infile.is_open())
    {
        cout<<"error opening file: "<<domain_definition<<endl;
        return;
    }

    int line_len = 1000;
    int domain_col = 4;
    vector<string> domain_info;
    vector<string> pdbData;
    string curr_target, domain_name;
    string target;
    int start, end;
    char line[line_len];
    vector<string> tokens;
    while(infile.getline(line,line_len))
    {
        if(line[0] == '#')
            continue;
        Utility::split(tokens, line, '\t');
        //get target file name and domains    
        Utility::split(domain_info, tokens[domain_col].c_str(), ':');
        target = Utility::split(domain_info[0],'-')[0];
        if(curr_target != target)
        {
            curr_target = target;
            cout<<curr_target<<'\n';
            //read pdbs
            Utility::readPDB(targetPath + curr_target + ".pdb", pdbData);
        }
        domain_name = domain_info[0]+ ".pdb";
        auto range = Utility::split(domain_info[1], '-');
        start = stoi(range[0]);
        end = stoi(range[1]);
        //cout<<'\t'<<domain_name<<'\t'<<start<<'\t'<<end<<endl;
        int start_idx,end_idx;

        Utility::writePDB(domainPath + domain_name, pdbData, start, end);
        
        cout<<endl;

    }

}
