#include "benchMark.hpp"
#include <fstream>
#include <iostream>
#include "Cathdata.hpp"
#include "lshbox/matrix.h"
#include <stdlib.h>
#include "Procruste.hpp"

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
        TargetName = Utility::split( TargetName, '.')[0] ;
        
        Protein_3D_Index p3d;
        p3d.desCalculate(Prediction_Path + TargetName + "\/"
                        , outputPath + TargetName + ".descriptor"
                        , outputPath + TargetName + ".list"
                );
        
        cout<<TargetName<<endl;

    }
}




void benchMark::getRMSD(string fileName, string path, vector<pair<string ,double> > & result)
{
    result.resize(0);
    vector<POINT> target;
    vector<POINT> prediction;
    //read target
    Cathdomain::readPDB(fileName, target);

    //get file list
    vector<string> files;
    Utility::listFiles(path, files);

    result.resize(files.size());
    //iterate through files 
    for(int i = 0; i < files.size(); i++)
    {
        cout<<files[i]<<'\t';
        result[i].first = files[i];
        Cathdomain::readPDB(path + files[i], prediction);
        result[i].second = Procruste::procruste(target, prediction);
        cout<<endl;
    }
}




void  benchMark::my_LGA_result()
{
    vector<string> files;
    vector<string> data;
    vector<pair<string,double> > result;
    Utility::listFiles(Target_Path, files);

    //set line length
    int lineLen = 200;
    char line[lineLen];

    //iterate through the direcroty
    for(int i=0; i < files.size(); i++)
    {
        string path = Prediction_Path + files[i] + "/";
        string fname = Utility::removeSuff(files[i]);
        getRMSD(Target_Path + files[i], path, result);
        if(result.size() > 0)
        {//save result
            data.resize(1);
            //print the header
            data[0] = "NAME   N2   RMSD   GDT_TS   LGA_S3   LGA_Q";
            //print the data
            for(int j = 0; j < result.size(); j++)
            {
                sprintf(line,"%s     %d     %4f     %5f    %5f     %5f"
                                        , result[j].first.c_str(),0,  result[j].second
                                        , 0.0, 0.0, 0.0, 0.0);

                data.push_back(string(line));
            }
        
            Utility::writeFile(Result_Path + fname+".SUMMARY.lga_sda.txt", data);

        }
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
    vector<string> pdbData_trimmed;
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
            cout<<curr_target<<'\n';
            //read pdbs
            if(!Utility::readPDB(targetPath + target + ".pdb", pdbData))
                continue;
            curr_target = target;
        }
        domain_name = domain_info[0]+ ".pdb";
        //get all ranges i
        vector<pair<int,int> > ranges;
        auto ranges_str = Utility::split(domain_info[1], ',');
        for(int i = 0; i < ranges_str.size(); i++)
        {
            auto range_str = Utility::split(ranges_str[i], '-');
            start = stoi(range_str[0]);
            end = stoi(range_str[1]);
            ranges.push_back(pair<int, int>(start, end));
        }
        //cout<<'\t'<<domain_name<<'\t'<<start<<'\t'<<end<<endl;

        pdb_range_copy(pdbData, pdbData_trimmed, ranges);
        Utility::writePDB(domainPath + domain_name, pdbData_trimmed);
        
        //cout<<endl;dd

    }

}




void benchMark:: split_to_equal(string dataPath, string splitPath, int segLen)
{
    vector<string> files;
    Utility::listFiles(dataPath, files);
    vector<string> pdbData; 
    vector<vector<string> > segs;
    
    for(int i=0; i < files.size(); i++)
    {
        Utility::readPDB(dataPath + files[i], pdbData); 
        Utility::splitPDB(pdbData, segs, segLen);
        for(int j=0; j < segs.size(); j++)
        {
            if(segs[j].size() < segLen)
                continue;
            string newName = Utility::removeSuff(files[i]);
            newName = newName + "_" +  to_string(j);
            Utility::writeFile(splitPath + newName, segs[j]);
        }
    }
}



void benchMark::prepare_directory(string targetPath, string dir_path, int segLen)
{
    vector<string> files;
    Utility::listFiles(targetPath, files);
    vector<string> pdbData; 
    vector<vector<string> > segs;

    for(int i=0; i < files.size(); i++)
    {
        Utility::readPDB(targetPath + files[i], pdbData); 

        Utility::splitPDB(pdbData, segs, segLen);
        for(int j=0; j < segs.size(); j++)
        {
            if(segs[j].size() < segLen)
                continue;
            string newName = Utility::removeSuff(files[i]);
            newName = newName + "_" +  to_string(j);
            string cmd = "mkdir " + dir_path + newName;
            cout<<cmd<<endl;
            system(cmd.c_str());
        }
    }

}




void benchMark:: split_to_equal_diff(string dataPath, string splitPath, int segLen)
{
    vector<string> files;
    Utility::listFiles(dataPath, files);
    vector<string> pdbData; 
    vector<vector<string> > segs;

    string targetName = Utility::split(dataPath, '/').back();
    for(int i=0; i < files.size(); i++)
    {
        Utility::readPDB(dataPath + files[i], pdbData); 
        
        Utility::splitPDB(pdbData, segs, segLen);
        for(int j=0; j < segs.size(); j++)
        {
            if(segs[j].size() < segLen)
                continue;
            string tmp = Utility::removeSuff(files[i]);
            string newName = tmp.substr(0,5)+ "_" + to_string(j) + tmp.substr(5);
            cout<<newName<<endl;
            Utility::writeFile(splitPath + targetName + "_" + to_string(j)  + '/' + newName,segs[j]);  
        }
    }
}



void benchMark::split_all_prediction(string predPath, string splitPath, int segLen)
{
    vector<string> dirs;
    Utility::listFiles(predPath, dirs, 0);
    
    for(int i=0; i < dirs.size(); i++)
    {
        if(dirs[i][0] == '.')
            continue;
        split_to_equal_diff(predPath + dirs[i] + '/', splitPath, segLen);
    }
}




void benchMark::
pdb_range_copy(vector<string>& from, vector<string>& to, vector<pair<int,int> > ranges)
{
    to.resize(0);
    int resSeq;
    int j;
    for(int i=0; i < from.size(); i++)
    {
        sscanf(from[i].c_str()+22,"%4i",&resSeq);
        for(j = 0; j < ranges.size(); j++)
            if(resSeq >= ranges[j].first && resSeq <= ranges[j].second)
                break;
        if(j >= ranges.size())
            continue;

        to.push_back(from[i]);
    }
}



void benchMark:: readResult()
{
    //read All target files
    CathData cdata(outputPath + "Target.list");
    cdata.readList();
    typedef double DATATYPE;
    lshbox::Matrix<DATATYPE> target_des(outputPath + "Target.descriptor");

    //cout<<cdata.size()<<'\n'<<target_des.getSize()<<'\n'<<target_des.getDim()<<endl;

    
    for(int i=0; i < cdata.size(); i++)
    {
        target tt;
        string str = cdata[i].getName();
        tt.target_name = Utility::removeSuff(str);
        //cout<<cdata[i].getName()<<endl;
        for(int j=0; j < target_des.getDim(); j++)
        {
            tt.Descriptor[j] = target_des[i][j];
            // cout<<target_des[i][j]<<' ';
        }
        result.push_back(tt);
       //cout<<endl;
    }


    // read all prediction files and correspond distance file
    
    //choose the distance measures we want
    vector<string> cols;
    cols.push_back("NAME");
    cols.push_back("N2");
    cols.push_back("GDT_TS");
    cols.push_back("RMSD");
    cols.push_back("LGA_S3");
    cols.push_back("LGA_Q");

    for(int i=0; i < result.size(); i++)
    {
        string listName,desName;
        listName = outputPath + result[i].target_name + ".list";
        desName = outputPath + result[i].target_name + ".descriptor";
        CathData sublist(listName);
        if(!sublist.readList())
            continue;
        lshbox::Matrix<DATATYPE> sub_des(desName);
        
        vector<vector<string> > dist_data;
       
        string fileName = Result_Path + result[i].target_name + ".SUMMARY.lga_sda.txt";
        if(!Utility::readTable(fileName, dist_data, cols))
            continue;
        //cout<<listName<<'\n'<<desName<<endl;
        for(int j = 0; j < sublist.size(); j++)
        {
            string str = sublist[j].getName();
            str = Utility::removeSuff(str); 
            //cout<<str<<endl;
            predict tp;
            
            //copy descriptor
            for(int k = 0 ; k < DESCRIPTOR_LENGTH ; k++)
                tp.Descriptor[k] = sub_des[j][k];
            
            tp.length = numeric_limits<double>::quiet_NaN();
            tp.GDT_TS = numeric_limits<double>::quiet_NaN(); 
            tp.RMSD = numeric_limits<double>::quiet_NaN();
            tp.LGA_S3 = numeric_limits<double>::quiet_NaN();
            tp.LGA_Q = numeric_limits<double>::quiet_NaN();
            


            result[i].predictions.insert(pair<string, predict>(str, tp));
        }
        // since the distance file of cast may have different target order
        // a search is needed here
        for(int j = 0; j < dist_data.size(); j++)
        {
            // search with name
            //cout<<dist_data[j][0]<<endl;
            string key = dist_data[j][0];
            key = Utility::removeSuff(key);
            auto p_it = result[i].predictions.find(key);
            if(p_it == result[i].predictions.end())
                continue;
            predict & tp = (*p_it).second;
            tp.length = stoi(dist_data[j][1]);
            tp.GDT_TS = stod(dist_data[j][2]);
            tp.RMSD = stod(dist_data[j][3]);
            tp.LGA_S3 = stod(dist_data[j][4]);
            tp.LGA_Q = stod(dist_data[j][5]);
            
        }
    }//end of read predition descriptors and distances

}



void benchMark:: calDist()
{
    for(int i=0; i < result.size(); i++)
    {   
        for(auto it = result[i].predictions.begin(); it != result[i].predictions.end(); ++it)
        {
            predict& tp = (*it).second;        
            tp.l1_dist = Utility::l1_norm(result[i].Descriptor, tp.Descriptor, DESCRIPTOR_LENGTH);

            tp.l2_dist = Utility::l2_norm(result[i].Descriptor, tp.Descriptor, DESCRIPTOR_LENGTH);
        }
    }

}



bool benchMark :: DataValidation()
{
    bool ret = true;
    for(int i=0; i < result.size(); i++)
    {   
        cout<<result[i].target_name<<": "<<result[i].predictions.size()<<"predicitons"<<'\n';
        for(int j =0; j < DESCRIPTOR_LENGTH; j++)
            if(isnan(result[i].Descriptor[j]) || isinf(result[i].Descriptor[j]))
            {
                cout<<" target descriptor error"<<'\n';
                ret = false;
                break;
            }

        for(auto it = result[i].predictions.begin(); it != result[i].predictions.end(); ++it)
        {
            predict tp = (*it).second;
            if(isnan(tp.RMSD) || isinf(tp.RMSD)
                ||isnan(tp.LGA_S3) || isinf(tp.LGA_S3)
                ||isnan(tp.LGA_Q) || isinf(tp.LGA_Q)
                ||isnan(tp.GDT_TS) || isinf(tp.GDT_TS)
              )
            {
                cout<<(*it).first<<":mearuse error"<<'\n';
                ret =  false;
            }

            for(int j =0; j < DESCRIPTOR_LENGTH; j++)
            if(isnan(tp.Descriptor[j]) || isinf(tp.Descriptor[j]))
            {
                cout<<(*it).first<<" descriptor error"<<'\n';
                ret =  false;
                break;
            }

            if(fabs(tp.l1_dist) > 1000 || fabs(tp.l2_dist)> 1000)
            {
                cout<<(*it).first<<"  warning: distance too large "<<tp.l1_dist<<'\t'<<tp.l2_dist<<endl;
            }

            if(tp.l1_dist < 0 || tp.l2_dist <0 )
                cout<<(*it).first<<"  warning: error negtive value "<<tp.l1_dist<<'\t'<<tp.l2_dist<<endl;
   
        }
    }
    return ret;
}


// format : 1 prediction name
//          2 lenght
//          3 l1_dist
//          4 l2_dist
//          5 RMSD
//          6 GDT_TS
//          7 LGA_S3
//          8 LGA_Q
bool benchMark:: saveResult(string path)
{
    bool ret =true;
    vector<string> data;
    for(int i=0; i < result.size(); i++)
    {   
        data.resize(0);
        if(result[i].predictions.size() == 0)
            continue;
        for(auto it = result[i].predictions.begin(); it != result[i].predictions.end() ;++it)
        {
            ostringstream strs;
            strs<<(*it).first<<"  ";          
            predict tp = (*it).second;
            strs<<tp.length<<"  "<<tp.l1_dist<<"  "<<tp.l2_dist<<" ";
            strs<<tp.RMSD<<"  "<<tp.GDT_TS<<"  "<<tp.LGA_S3<<"  "<<tp.LGA_Q;
            data.push_back( strs.str());
        }

        if(Utility::writeFile(path + result[i].target_name + ".result", data))
            ret = false;
    }
    return ret;
}
