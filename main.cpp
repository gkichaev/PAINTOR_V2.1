//
//  main.cpp
//  TestTarget
//
//  Created by Gleb Kichaev on 4/16/15.
//  Copyright (c) 2015 Gleb Kichaev. All rights reserved.
//


#include "Header.h"
#include "Functions.cpp"

using namespace std;
using namespace Eigen;



void Welcome_Message(){
    cout << "Welcome to PAINTOR v2.0. Copyright Gleb Kichaev 2015." << endl;
    cout << "For questions or bug reports please contact: gkichaev@ucla.edu \n \n "<< endl;
    cout << "For required files and format specifications see User Manual \n \n" << endl;
    cout << "Usage: PAINTOR -input.files [input_filename] -in [input_directory] -out [output_directory] -Zhead [Zscore_header(s)] -LDname [LD_suffix(es)]  -annotations [annot_name1,annot_name2,...]  <other options> \n"<< endl;
    cout << "OPTIONS: -flag \t Description [default setting]  \n" << endl;
    cout << "-input \t (required) Filename of the input file containing the list of the fine-mapping loci [default: input.files]" << endl;
    cout << "-c \t The number of causal variants to consider per locus [default: 2]" << endl;
    cout << "-Zhead \t (required) The name(s) of the Zscores in the header of the locus file (comma separated) [default: N/A]" << endl;
    cout << "-LDname \t (required) Suffix(es) for LD files. Must match the order of Z-scores in locus file (comma separated) [Default:N/A]" << endl;
    cout << "-annotations \t The names of the annotations to include in model (comma separted) [default: N/A]" << endl;
    cout << "-in \t Input directory with all run files [default: ./ ]" << endl;
    cout << "-out \t Output directory where output will be written [default: ./ ]" << endl;
    cout << "-Gname \t Output Filename for enrichment estimates [default: Enrichment.Estimate]" << endl;
    cout << "-RESname \t Suffix for ouput files of results [Default: results] " << endl;
    cout << "-ANname \t Suffix for annotation files [Default: annotations]" << endl;
    cout << "-MI \t Maximum iterations for algorithm to run [Default: 10]" << endl;
    cout << "-post1CV \t fast conversion of Z-scores to posterior probabilites assuming a single casual variant and no annotations [Default: False]" << endl;
    cout << "-GAMinital \t inititalize the enrichment parameters to a pre-specified value (comma separated) [Default: 0,...,0]" << endl;

    cout << endl << endl ;
}

void Check_Flags(){

}


int main(int argc, const char * argv[])
{

    string in_dir = "./";
    string input_files = "input.files";
    string out_dir = "./";
    vector<string> annot_names;
    int max_causal = 2;
    string gammaName = "Enrichment.Values";
    string likeli_name;
    string single_post_flag;
    int maxIter= 10;
    string LD_suffix = "ld";
    string annot_suffix = "annotations";
    string results_suffix = "results";
    vector <string> LD_all_names;
    vector<string> z_headers;

    single_post_flag = "False";

    vector<VectorXd> z_score_loc;
    vector<string> snp_info;
    MatrixXd chol_factor;
    MatrixXd chol_factor_inv;
    MatrixXd out_annotations;
    vector<string> annot_header;
    string file_list_name = "input.files";
    string gamma_initial;

    if(argc < 2){
        Welcome_Message();
        return 0;
    }
       for(int i = 1; i < argc; i++){
           string argComp = argv[i];


           if(argComp.compare("-input")== 0){
               input_files = argv[i+1];
           }

           else if(argComp.compare("-in") == 0){
               in_dir = argv[i+1];
               string last_char = string(&in_dir[in_dir.size()-1]);
               if(last_char.compare("/")!=0){
                   in_dir = in_dir + "/";
               }
           }

           else if(argComp.compare("-out") == 0){
               out_dir = argv[i+1];
               string last_char = string(&out_dir[out_dir.size()-1]);
               if(last_char.compare("/")!=0){
                   out_dir = out_dir + "/";
               }
           }

           else if(argComp.compare("-c") == 0){
               max_causal = stoi(argv[i+1]);
           }


           else if(argComp.compare("-LDname") == 0){
               string temp_ldnames = argv[i+1];
               size_t n = count(temp_ldnames.begin(), temp_ldnames.end(), ',');
               if(n >0) {
                   LD_all_names = split(temp_ldnames, ',');
               }
               else{
                   LD_all_names = {temp_ldnames};
               }
           }

           else if(argComp.compare("-Zhead") == 0){
               string header_temp = argv[i+1];
               size_t n = count(header_temp.begin(), header_temp.end(), ',');
               if(n>0){
                   z_headers = split(header_temp, ',');
               }
               else{
                   z_headers = {header_temp};
               }
           }

           else if(argComp.compare("-annotations") == 0){
               string input_annotations = argv[i+1];
               size_t n = count(input_annotations.begin(), input_annotations.end(), ',');
               if(n >0){
                   annot_names = split(input_annotations, ',');
               }
               else{
                   annot_names= {input_annotations};
               }

           }

           else if(argComp.compare("-Gname") == 0){
               gammaName = argv[i+1];
           }

           else if(argComp.compare("-Lname") == 0){
               likeli_name = argv[i+1];
           }

           else if(argComp.compare("-post1CV") == 0){
               single_post_flag = argv[i+1];
           }

           else if(argComp.compare("-MI") == 0){
               maxIter = stoi(argv[i+1]);
           }

           else if(argComp.compare("-ANname") == 0){
               annot_suffix = argv[i+1];
           }

           else if(argComp.compare("-RESname") == 0){
               results_suffix = argv[i+1];
               results_suffix = "." + results_suffix;
           }

           else if(argComp.compare("-GAMinitial") == 0){
               gamma_initial = argv[i+1];
           }

       }

    /* initialize PAINTOR model parameters */

    vector<vector<VectorXd>> all_transformed_statistics;
    vector<vector<VectorXd>> all_lambdas;
    vector<vector<MatrixXd>> all_upper_cholesky;
    vector<vector<string>> all_snp_info;
    vector<MatrixXd> all_annotations;
    vector<string> all_headers;
    Get_All_Input(input_files, in_dir, z_headers, annot_names, all_transformed_statistics, all_lambdas, all_upper_cholesky, all_annotations, LD_all_names, annot_suffix,all_snp_info,all_headers);
    CausalProbs runProbs;
    VectorXd gamma_estimates(all_annotations[0].cols());
    gamma_estimates.setZero();

    if(gamma_initial.size() > 0){
        vector<string> gamma_initial_split = split(gamma_initial, ',');
        if(gamma_initial_split.size() != annot_names.size()+1){
            cout << "Warning: Incorrect number of Enrichment parameters specified. Pre-setting all paramters to zero" << endl;
        }
        else{

        }
    }

    /* run PAINTOR */

     if (single_post_flag.compare("True")==0) {
        if (z_headers.size() > 1) {
            cout << "Single Posterior Conversion not valid for more than 1 set of Z-scores per locus" << endl;
            return 0;
        }
        else{
            vector<VectorXd> all_single_post;
            VectorXd single_post;
            for(unsigned int i =0; i < all_transformed_statistics.size(); i++){
                single_post = Zscores2Post(all_transformed_statistics[i][0]);
                all_single_post.push_back(single_post);
            }
            Write_All_Output(input_files, out_dir, results_suffix, runProbs, all_snp_info, gamma_estimates,gammaName, 0, likeli_name, all_headers, annot_names);


        }
    }

    else{
        double Final_loglikeli = EM_Run_chol(runProbs, maxIter , all_transformed_statistics,all_lambdas, gamma_estimates, all_annotations, all_upper_cholesky, max_causal);

        Write_All_Output(input_files, out_dir, results_suffix, runProbs, all_snp_info, gamma_estimates,gammaName, Final_loglikeli, likeli_name, all_headers, annot_names);
    }

    return 0;
}

/*id Write_All_Output(string& input_files, string& out_dir,)

                      */