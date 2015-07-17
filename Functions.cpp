
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <iterator>
#include <string>
#include <iomanip>
#include <sstream>
#include <Eigen>
#include "Header.h"
#include <algorithm>

using namespace std;
using namespace Eigen;


void Write_Posterior(string& out_dir, string & out_name, VectorXd& locus_results, vector<string>& locus_info, string& header ){

    vector<double> probs_as_strings;
    ofstream myfile;
    string fname = out_dir+out_name;
    myfile.open(fname);
    for(int i =0; i<locus_results.size(); i++){
        probs_as_strings.push_back((locus_results[i]));
    }
    myfile << header + " Posterior_Prob" << endl;
    for(int i=0 ;i<locus_results.size(); i++){
        myfile << locus_info[i] + " ";
        myfile << probs_as_strings[i];
        myfile << "\n";
    }

    myfile.close();


}

void Write_All_Output(string& input_files, string& out_dir, string& out_suffix, CausalProbs& results, vector<vector<string>> & all_locus_info,
                      VectorXd & gamma_estimates, string& Gname, double log_likeli, string & Lname, vector<string>& all_headers, vector<string>& annot_names){
    ifstream input_list;
    input_list.open(input_files);
    string locus_name;
    unsigned int j=0;
    while(getline(input_list,locus_name)){
        string out_fname = locus_name + "." + out_suffix;
        Write_Posterior(out_dir, out_fname, results.probs_locs[j], all_locus_info[j], all_headers[j]);
        j++;
        if(j >= all_locus_info.size()) {
            break;
        }
    }
    input_list.close();

    ofstream gamma_out;
    string G_open = out_dir+Gname;
    gamma_out.open(G_open);

    gamma_out << "Baseline";
    for(unsigned int i =0; i < annot_names.size(); i++){
        gamma_out << " ";
        gamma_out << annot_names[i];
    }
    gamma_out << endl;

    gamma_out << gamma_estimates[0];
    for(unsigned int i=1; i < gamma_estimates.size(); i ++){
        gamma_out << " ";
        gamma_out << gamma_estimates[i];
    }
    gamma_out << endl;
    gamma_out.close();

    string L_open = out_dir+Lname;
    ofstream likeli_out;
    likeli_out.open(L_open);
    likeli_out<< std::setprecision(10) << log_likeli;
    likeli_out << "\n";
    likeli_out.close();

}


double EM_Run_chol(CausalProbs &probabilites, int iter_max, vector<vector<VectorXd>> &Zscores, vector<vector<VectorXd>> &Lambdas, VectorXd &beta_int, vector<MatrixXd> &Aijs, vector<vector<MatrixXd>> &upper_chol, int numCausal , vector<vector<MatrixXd>>& sigmas, string & version){
    vector<double> beta_run = eigen2vec(beta_int);
    //vector<MatrixXd> Invert_LD;
    //SigmaInvert(Sigmas, Invert_LD);
    ObjectiveData Opt_in;
    Opt_in.Aijs = Stack_Matrices(Aijs);

    int iterations = 0;
    VectorXd beta_update;
    double likeliOld = 0;
    double likeli =1;
    while(iterations < iter_max){
        likeli = Estep_chol(Zscores, Lambdas, beta_int, Aijs, upper_chol, probabilites, numCausal, sigmas,version);
        Opt_in.probs = probabilites.probs_stacked;

        void *opt_ptr = &Opt_in;
        Optimize_Nlopt(beta_run, -20, 20, beta_run[0], opt_ptr);
        beta_update = vector2eigen(beta_run);
        cout << std::setprecision(9) << "Log likelihood at iteration " << iterations << ": " <<likeli << endl;
        cout << "Parameter Estimates:" << endl << beta_update << endl << endl;
        if(abs(likeli - likeliOld) < 0.01){
            beta_int = beta_update;
            break;
        }
        else{
            beta_int = beta_update;
            likeliOld = likeli;
            iterations++;
        }
    }
    return(likeli);
}


double Estep_chol(vector<vector<VectorXd>> &Zscores, vector<vector<VectorXd>> &Lambdas, VectorXd &betas, vector<MatrixXd> &Aijs, vector<vector<MatrixXd>>& upper_chol, CausalProbs &E_out, int numberCausal, vector<vector<MatrixXd>>& sigmas, string& version){

    vector<VectorXd> marginal_i;
    VectorXd temp;
    VectorXd exp_temp;
    vector<double> stacker;
    vector<double> stack_temp;
    double fullLikeli = 0;

    for(unsigned i = 0; i < Zscores.size(); i ++){
        VectorXd temp(Zscores[i][0].size());

        if(version.compare("old")==0) {
            NewPost_chol(temp, Zscores[i], Lambdas[i], betas, Aijs[i], upper_chol[i], numberCausal, fullLikeli);
        }
        else  if(version.compare("default")==0){
            NewPost_chol_ncp(temp, Zscores[i], Lambdas[i], betas, Aijs[i], upper_chol[i], numberCausal, fullLikeli, sigmas[i]);
        }
        exp_temp = temp.array().exp();
        marginal_i.push_back(exp_temp);
        stack_temp =  eigen2vec(exp_temp);
        stacker.insert(stacker.end(), stack_temp.begin(), stack_temp.end());
    }

    E_out.probs_locs = marginal_i;
    E_out.probs_stacked = stacker;

    return(fullLikeli);
}




void NewPost_chol(VectorXd& Marginal, vector<VectorXd>& Zs, vector<VectorXd>& Lams, VectorXd& beta, MatrixXd& Aj,  vector<MatrixXd>& upper_chol, int NC, double& fullLikeli){
    int numsnps = Zs[0].size();
    int num_pops = Zs.size();
    double runsum = 0;
    for(int i =0 ; i < Marginal.size(); i++){
        Marginal(i) = -1e150;
    }
    double c_prob = 0;
    double z_prob = 0;
    double sum = 0;

    VectorXd causConfig(numsnps);
    causConfig.setZero();
    VectorXd causIndex(NC);
    causIndex.setZero();
    VectorXd evalPrior(numsnps);
    EvalAijs(Aj, beta, evalPrior);

    c_prob = Prob_CNew(evalPrior, beta, causConfig);
    double z_pop_prob = 0;
    for (int i = 0;  i<num_pops ; i++) {
        z_pop_prob=EvaluateLogMvn_Cholesky(Zs[i], causConfig, Lams[i], upper_chol[i]);
        z_prob += z_pop_prob;
    }

    sum = c_prob+z_prob;
    runsum = sum;
    int counter = 1;
    while(NextCombo(causIndex, NC, numsnps+1) == 1){
        BuildCausalVector(causConfig, causIndex);
        c_prob = Prob_CNew(evalPrior, beta, causConfig);

        z_prob = 0;
        for (int i = 0;  i<num_pops ; i++) {
            z_pop_prob=EvaluateLogMvn_Cholesky(Zs[i], causConfig, Lams[i], upper_chol[i]);
            z_prob += z_pop_prob;
        }
        sum = c_prob+z_prob;
        runsum = LogSum(runsum, sum);

        for(int j = 0; j < causIndex.size(); j++){
            if(causIndex[j] >0){
                Marginal[causIndex[j]-1] = LogSum(Marginal[causIndex[j]-1], sum);
            }

        }
        counter ++;
    }
    for(int f = 0 ; f < Marginal.size(); f++){
        Marginal[f] = Marginal[f]- runsum;
    }
    fullLikeli = fullLikeli+runsum;

}



void NewPost_chol_ncp(VectorXd& Marginal, vector<VectorXd>& Zs, vector<VectorXd>& Lams, VectorXd& beta, MatrixXd& Aj,  vector<MatrixXd>& upper_chol, int NC, double& fullLikeli, vector<MatrixXd> sigma_loc){
    int numsnps = Zs[0].size();
    int num_pops = Zs.size();
    double runsum = 0;
    for(int i =0 ; i < Marginal.size(); i++){
        Marginal(i) = -1e150;
    }
    double c_prob = 0;
    double z_prob = 0;
    double sum = 0;

    VectorXd causConfig(numsnps);
    causConfig.setZero();
    VectorXd causIndex(NC);
    causIndex.setZero();
    VectorXd evalPrior(numsnps);
    EvalAijs(Aj, beta, evalPrior);

    c_prob = Prob_CNew(evalPrior, beta, causConfig);
    double z_pop_prob = 0;
    for (int i = 0;  i<num_pops ; i++) {
        z_pop_prob=EvaluateLogMvn_NCP(Zs[i], causConfig, Lams[i], upper_chol[i], sigma_loc[i]);
        z_prob += z_pop_prob;
    }

    sum = c_prob+z_prob;
    runsum = sum;
    int counter = 1;
    while(NextCombo(causIndex, NC, numsnps+1) == 1){
        BuildCausalVector(causConfig, causIndex);
        c_prob = Prob_CNew(evalPrior, beta, causConfig);

        z_prob = 0;
        for (int i = 0;  i<num_pops ; i++) {
            z_pop_prob=EvaluateLogMvn_NCP(Zs[i], causConfig, Lams[i], upper_chol[i], sigma_loc[i]);
            z_prob += z_pop_prob;
        }
        sum = c_prob+z_prob;
        runsum = LogSum(runsum, sum);

        for(int j = 0; j < causIndex.size(); j++){
            if(causIndex[j] >0){
                Marginal[causIndex[j]-1] = LogSum(Marginal[causIndex[j]-1], sum);
            }

        }
        counter ++;
    }
    for(int f = 0 ; f < Marginal.size(); f++){
        Marginal[f] = Marginal[f]- runsum;
    }
    fullLikeli = fullLikeli+runsum;

}



// PAINTOR 2.1 Functions

void SigmaInvert(vector<vector<MatrixXd>> &Sigma_locs, vector<vector<MatrixXd>> &Inv_locs){

    for(unsigned i = 0; i < Sigma_locs.size(); i++){
        vector<MatrixXd> inv_ld;
        for(unsigned j = 0; j < Sigma_locs[i].size(); j++) {
            FullPivLU<MatrixXd> ldt(Sigma_locs[i][j]);
            inv_ld.push_back(ldt.inverse());
        }
        Inv_locs.push_back(inv_ld);
    }
}



void CorrectNCP(VectorXd & Z_vec, VectorXd & C_vec, VectorXd& NCP,  MatrixXd& Sig){
    vector<int> caus_ind;
    for(unsigned i = 0; i < C_vec.size(); i++){
        if(C_vec[i] ==1){
            caus_ind.push_back(i);
        }
    }
    if(caus_ind.size() > 1){
        VectorXd temp_NCP(caus_ind.size());
        MatrixXd temp_LD(caus_ind.size(),caus_ind.size());
        for(int i = 0; i < caus_ind.size(); i++){
            temp_NCP[i] = Z_vec[caus_ind[i]];
            for(int j = 0; j < caus_ind.size(); j++){
                if(i==j){
                    temp_LD(i,j) = 1.0001;
                }
                else{
                    temp_LD(i,j) = Sig(caus_ind[i],caus_ind[j]);
                }
            }
        }

        VectorXd transNCP;
        transNCP = temp_LD.inverse()*temp_NCP;
        for(unsigned i =0 ; i < caus_ind.size(); i++){
            NCP[caus_ind[i]] = transNCP[i];
        }
    }
    else if(caus_ind.size() == 1){
        NCP[caus_ind[0]] = Z_vec[caus_ind[0]];
    }

}

inline  double EvaluateLogMvn_NCP(VectorXd& Z_vec,  VectorXd C_vec,  VectorXd& Lam_vec,  MatrixXd& chol_factor,  MatrixXd& Sig){
    VectorXd NCP(Z_vec.size());
    NCP.setZero();
    CorrectNCP(Lam_vec, C_vec, NCP, Sig);
    VectorXd mu = chol_factor*NCP;
    VectorXd Z_mu = Z_vec-mu;
    double prob = (-.5)*(Z_mu.dot(Z_mu));
    return(prob);
}


//



inline  double EvaluateLogMvn_Cholesky(VectorXd& Z_vec,  VectorXd& C_vec,  VectorXd& Lam_vec, MatrixXd& upper_chol){

    VectorXd elm_wise = C_vec.cwiseProduct(Lam_vec);
    VectorXd mu = upper_chol*elm_wise;
    VectorXd Z_mu = Z_vec-mu;
    double prob = (-.5)*(Z_mu.dot(Z_mu));
    return(prob);
}



vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

vector<string> split(const string &s, char delim) {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

void Read_Locus(string &input_directory, string& fname, vector<string> & zname,  vector<VectorXd>& all_association_stat, vector<string>& snp_info, string & header){

    ifstream locus;
    locus.open(input_directory+fname);
    getline(locus, header);
    char delimiter = ' ';
    vector<string> split_header = split(header, delimiter);

    //vector<string> fields = {"SNP_ID", "CHR", "POS", "EFFECT_ALLELE"};
    vector<int> z_index;
    for(unsigned int i=0; i < split_header.size(); i++){
        for(unsigned int j=0; j < zname.size(); j++){
            if(zname[j].compare(split_header[i])==0){
                z_index.push_back(i);
            }
        }
    }
    string input_line;
    vector<string> split_line;
    vector<vector<string>> snp_info_split;
    string stat_holder;

    while(getline(locus,input_line)){
        split_line = split(input_line, ' ');
        snp_info.push_back(input_line);
        snp_info_split.push_back(split_line);
    }
    locus.close();
    unsigned long num_snps = snp_info.size();
    VectorXd association_stat;
    for(unsigned int i =0; i < z_index.size();i++){
        association_stat.resize(num_snps,1);
        for(unsigned int j=0; j < num_snps; j++){
            if(snp_info_split[j][z_index[i]].compare("NA") == 0){
                association_stat[j] = 0;
            }
            else{
                association_stat[j] = stod(snp_info_split[j][z_index[i]]);
            }
        }
        all_association_stat.push_back(association_stat);
    }

}

double Regularize_LD(MatrixXd & ld_mat){
    double reg_factor = 0.001;
    int reg_flag = 0;
    double det;
    VectorXd X(ld_mat.rows());
    MatrixXd Y = X.asDiagonal();
    while(reg_flag == 0){
        X.fill(reg_factor);
        Y = X.asDiagonal();
        det = (ld_mat+Y).determinant();
        if(det <= 1e-100){
            reg_factor = reg_factor*1.25;
        }
        else{
            reg_flag = 1;
        }

    }
    X.fill(reg_factor);
    Y = X.asDiagonal();
    ld_mat= ld_mat +Y;
    return reg_factor;
}



void Read_LD(string &input_directory, string& fname , MatrixXd& chol_factor, MatrixXd& chol_factor_inv){
    ifstream ld_file;
    ld_file.open(input_directory+fname);
    string ld_line;
    vector<string> ld_line_split;
    vector<vector<string>> all_ld_line_split;
    while(getline(ld_file, ld_line)){
        ld_line_split = split(ld_line, ' ');
        all_ld_line_split.push_back(ld_line_split);
    }
    ld_file.close();
    unsigned long num_snps = all_ld_line_split.size();
    MatrixXd numeric_ld(num_snps, num_snps);
    double val;
    for(unsigned int i = 0; i < num_snps; i++){
        for(unsigned int j = 0; j < num_snps; j++){
            val =  stod(all_ld_line_split[i][j]);
            if(isfinite(val)) {
                numeric_ld(i, j) = stod(all_ld_line_split[i][j]);
            }
            else{
               cout <<  "Error: Infinite value in LD matrix detected!!" << endl;
            }


        }
    }
    double reg_factor;
    reg_factor = Regularize_LD(numeric_ld);


    chol_factor = numeric_ld.llt().matrixL(); //cholesky decompostion
    chol_factor.transpose();
    chol_factor.triangularView<Eigen::Upper>();
    chol_factor_inv = chol_factor.inverse(); //cholesky decomp inverse
    chol_factor_inv.triangularView<Eigen::Lower>();

}


void Read_Annotations(string &input_directory, string& fname , vector<string>& model_annotations, MatrixXd& out_annotations, vector<string> & annotation_header_split){
    ifstream annotation_file;
    annotation_file.open(input_directory+fname);
    string annotation_header;
    getline(annotation_file,annotation_header);
    annotation_header_split = split(annotation_header, ' ');
    string annotation_line;
    vector<string> annotation_line_split;
    vector<vector<string>> all_annotation_line_split;
    while(getline(annotation_file, annotation_line)){
        annotation_line_split = split(annotation_line, ' ');
        all_annotation_line_split.push_back(annotation_line_split);
    }
    unsigned long num_annotations = annotation_header_split.size();
    unsigned long num_snps = all_annotation_line_split.size();
    MatrixXd annotation_matrix(num_snps, num_annotations);
    for(unsigned int i = 0; i < num_snps; i++){
        for(unsigned int j= 0; j < num_annotations; j++){
            annotation_matrix(i,j)= stod(all_annotation_line_split[i][j]);
        }
    }
    MatrixXd output;
    if(model_annotations.size()>0){
        vector<int> annotation_index;
        bool annotation_found;
        for(unsigned int i = 0; i < model_annotations.size();  i++){
            annotation_found=false;
            for(unsigned int j=0; j < annotation_header_split.size(); j++){
                if(annotation_header_split[j].compare(model_annotations[i])==0){
                    annotation_index.push_back(j);
                    annotation_found=true;
                }
            }
            if(!annotation_found){
                cout << "Error: Annotation " << model_annotations[i] << " not found in annotation file" << endl;
            }
        }
        output.resize(num_snps, model_annotations.size()+1);
        VectorXd A0(num_snps);
        A0.setOnes();
        output.col(0) = A0;
        for(unsigned int i = 0; i < model_annotations.size(); i++){
            output.col(i+1) = annotation_matrix.col(annotation_index[i]);
        }
    }
    else{
        output.resize(num_snps, 1);
        VectorXd A0(num_snps);
        A0.setOnes();
        output.col(0) = A0;
    }
    out_annotations=output;
}

void Generate_Lambda(VectorXd& Zscore, VectorXd& Lambda, double minNCP){
    Lambda.resize(Zscore.size(),1);
    for(int j=0; j <Zscore.size(); j++){
        if(Zscore[j]>0){
            if(Zscore[j] < minNCP){
                Lambda[j] = minNCP;
            }
            else{
                Lambda[j] = Zscore[j];
            }
        }
        else{
            if(Zscore[j] > -1*minNCP){
                Lambda[j] = -1*minNCP;
            }
            else{
                Lambda[j] = Zscore[j];
            }
        }
        if(Zscore[j]==0){
            Lambda[j] =0;
        }
    }
}

void Get_All_Input(string& file_list_name, string& directory, vector<string> znames, vector<string>& model_annotations, vector<vector<VectorXd>>& all_transformed_statistics, vector<vector<VectorXd>>& all_lambdas ,vector<vector<MatrixXd>> &all_upper_cholesky,  vector<MatrixXd>& all_annotations, vector<string>& LD_suffix, string& annotation_suffix,
                   vector<vector<string>>& all_SNP_info, vector<string>& all_headers) {
    ifstream input_list;
    input_list.open(file_list_name);
    string locus_name;
    vector<vector<MatrixXd>> all_lower_cholesky;
    string header;
    while(getline(input_list,locus_name)){
        cout << "Reading in files for: " << locus_name << endl;
        vector<VectorXd> locus_statistics;
        vector<string> snp_info;
        vector<MatrixXd> locus_ld;
        MatrixXd locus_annotations;
        Read_Locus(directory, locus_name, znames, locus_statistics, snp_info, header);
        all_headers.push_back(header);
        all_SNP_info.push_back(snp_info);
        vector<MatrixXd> locus_chol_factors;
        vector<VectorXd> transformed_locus_statistics;
        vector<VectorXd> locus_lambdas;

        for(unsigned int i = 0; i < LD_suffix.size(); i++){
            MatrixXd chol_factor;
            MatrixXd chol_factor_inv;
            string ld_name = locus_name+"."+LD_suffix[i];
            VectorXd transformed_z;
            VectorXd single_lambda;
            Read_LD(directory, ld_name, chol_factor, chol_factor_inv);
            locus_chol_factors.push_back(chol_factor);
            transformed_z = chol_factor_inv*locus_statistics[i];
            transformed_locus_statistics.push_back(transformed_z);
            Generate_Lambda(locus_statistics[i], single_lambda, 3.7);
            locus_lambdas.push_back(single_lambda);
        }
        all_transformed_statistics.push_back(transformed_locus_statistics);
        all_lambdas.push_back(locus_lambdas);
        all_lower_cholesky.push_back(locus_chol_factors);
        string annot_name = locus_name+"."+annotation_suffix;
        vector<string> annotation_header;
        Read_Annotations(directory, annot_name, model_annotations, locus_annotations, annotation_header);
        all_annotations.push_back(locus_annotations);
    }

    //Transpose the cholesky factors to make the upper triangular

    for(unsigned long i =0; i < all_lower_cholesky.size(); i++){
        vector<MatrixXd> trans_chol;
        for(unsigned long j =0; j < all_lower_cholesky[i].size(); j++){
            trans_chol.push_back(all_lower_cholesky[i][j].transpose());
        }
        all_upper_cholesky.push_back(trans_chol);
    }
}

bool Check_Locus_File(vector<VectorXd>& association_stats, vector<MatrixXd>& ld_matrices, MatrixXd& annotation_matrix, string& locus_name){
    bool pass = true;
    unsigned long numsnps_annots = annotation_matrix.rows();
    for(unsigned int i = 0; i< association_stats.size(); i++){
        unsigned long numsnps_assoc = association_stats[i].size();
        unsigned long numsnps_ld = ld_matrices[i].cols();
        if(numsnps_annots != numsnps_assoc || numsnps_assoc != numsnps_ld || numsnps_annots!=numsnps_ld) {
            pass = false;
            cout << "Error: Number of SNPs at Locus " + locus_name + " does not match:" << endl;
            cout << "Associations statistics in population "  <<  i+1 << numsnps_assoc << " SNPS";
            cout << "LD in population "  <<  i+1 << numsnps_ld << " SNPS";
            cout << "Annotations"  <<  numsnps_annots << " SNPS";
        }
    }
    return pass;
}



vector<int> ParseIndex(string& input){
    vector<int> output;
    string delimiter = ",";
    auto start = 0U;
    auto end = input.find(delimiter);
    while(end != std::string::npos){
        output.push_back(stoi(input.substr(start, end - start)));
        start = int(end + delimiter.length());
        end = input.find(delimiter, start);
    }
    output.push_back(stoi(input.substr(start,end)));
    return(output);
}



void MakeAnnotations(vector<MatrixXd>& annotsRun, vector<MatrixXd>& allAnnots,  vector<int>& indices){
    for(unsigned int i = 0; i < allAnnots.size(); i++){
        MatrixXd temp = allAnnots[i];
        MatrixXd output(temp.rows(), indices.size()+1);
        VectorXd A0(temp.rows());
        A0.setOnes();
        output.col(0) = A0;
        for(unsigned int j= 0; j < indices.size(); j++){
            output.col(j+1) = temp.col(indices[j]);
        }
        annotsRun.push_back(output);
    }
}

/*


 */


VectorXd GetVector(string filename){

    ifstream input_file(filename);
    istream_iterator<double> start(input_file), end;
    vector<double> vec(start, end);
    input_file.close();

    VectorXd outvec(vec.size());
    for(unsigned i = 0; i < vec.size(); i++){
        outvec(i) = vec[i];
    }

    return(outvec);

}



vector<double> eigen2vec(VectorXd &vec){
    vector<double> outvec;
    for(unsigned i = 0; i < vec.size(); i++){
        outvec.push_back(vec[i]);
    }
    return(outvec);
}

VectorXd vector2eigen(vector<double> &vec){
    VectorXd outVec(vec.size());
    for(unsigned i =0; i <vec.size(); i++){
        outVec[i] = vec[i];
    }
    return(outVec);
}

vector<vector<double> > eigen2Matrix(MatrixXd &mat){
    vector<vector<double> > outMat;
    for(int i = 0; i < mat.rows(); i++){
        vector<double> row;
        for(int j = 0; j < mat.cols(); j ++){
            row.push_back(mat(i,j));
        }
        outMat.push_back(row);
    }
    return(outMat);
}


inline double Prob_Cijs(VectorXd& beta , VectorXd& aij){

    double dotprod = beta.dot(aij);
    double prob = 1/(1+exp(dotprod));
    return(prob);
}



void EvalAijs(MatrixXd& Aj,VectorXd& beta, VectorXd& priorJ){
    for(int i= 0; i < Aj.rows(); i ++){
        VectorXd Aij_temp = Aj.row(i);
        priorJ[i] = Prob_Cijs(beta, Aij_temp);
    }
}


double Prob_CNew(VectorXd& priorJ,VectorXd& beta, VectorXd& C_vector){
    double probs = 0;
    for(int i = 0; i < C_vector.size(); i++){
        if(C_vector[i] == 1){
            probs += log(priorJ[i]);
        }
        else{
            probs += log(1-priorJ[i]);
        }
    }
    return(probs);
}



VectorXd Zscores2Post(VectorXd& Zs){
    VectorXd post(Zs.size());
    VectorXd Zsq  = Zs.array().square();
    for(int i = 0; i < Zsq.size(); i ++){
        VectorXd Ztemp = (Zsq.array() - Zsq[i])/2;
        VectorXd Zexp = Ztemp.array().exp();
        post[i] = 1/Zexp.sum();
    }
    return(post);
}



inline double LogSum(double val1, double val2){
    double logsum = 0;
    if(val1 > val2){
        logsum = log(1 + exp(val2-val1)) + val1;
    }
    else{
        logsum = log(1 + exp(val1-val2)) + val2;
    }

    return(logsum);
}



double dot_prod(vector<double>& vec1, vector<double>& vec2){
    double runsum = 0;
    for(unsigned i = 0; i < vec1.size(); i++){
        runsum += vec1[i]*vec2[i];
    }
    return(runsum);
}

vector<double> scalar_product(const double& scalar, vector<double>& vec){
    vector<double> scalvec(vec.size(),0);
    for(unsigned i = 0; i < vec.size(); i++){
        scalvec[i] = vec[i]*scalar;
    }
    return(scalvec);
}



vector<double> Vector_Sum(vector<double>& vec1, vector<double>& vec2){
    vector<double> outsum;
    for(unsigned i = 0; i< vec2.size(); i++){
        outsum.push_back(vec1[i] +vec2[i]);
    }
    return(outsum);
}



vector<double> GradientFxn(vector<double>& betas, ObjectiveData *in_data){
    ObjectiveData f_data;
    f_data.probs = in_data->probs;
    f_data.Aijs = in_data->Aijs;
    int numsnps = f_data.probs.size();
    double cij1 = 0;
    double cij0 = 0;
    double dp  = 0;
    vector<double> aij1(numsnps, 0);
    vector<double> aij0(numsnps, 0);
    vector<double> aij1_0(numsnps,0);
    vector<double> aij_out(numsnps,0);

    for(int i = 0; i < numsnps; i ++){
        dp = dot_prod(betas, f_data.Aijs[i]);
        cij1 = f_data.probs[i]*1/(1+exp(-1*dp));
        cij0 = (1-f_data.probs[i])*1/(1+exp(dp));
        aij1= scalar_product(-1*cij1, f_data.Aijs[i]);
        aij0 = scalar_product(cij0, f_data.Aijs[i]);
        aij1_0 = Vector_Sum(aij1, aij0);
        aij_out = Vector_Sum(aij_out, aij1_0);
    }

    return(aij_out);
}



double ObjectiveFxn(const vector<double> &x, vector<double> &grad, void *data){

    ObjectiveData *f_data = static_cast<ObjectiveData*>(data);
    int numsnps = f_data -> probs.size();
    double marginal_i;
    double runsum = 0;
    double cij1 = 0;
    double cij0 = 0;
    double dp  = 0;
    vector<double> temp;
    vector<double> betas = x;
    grad = GradientFxn(betas, f_data);

    for(int i = 0; i < numsnps; i++){
        temp = f_data -> Aijs[i];
        marginal_i = f_data -> probs[i];
        dp = dot_prod(betas, temp);
        cij1 = marginal_i*log(1+exp(dp));
        cij0 = (1-marginal_i)*log(1+exp(-dp));
        runsum = runsum - cij1 - cij0;

    }

    return(runsum);
}


void Optimize_Nlopt(vector<double>& x, double lower, double upper, double betaZero,  void* in_data){

    opt opt(LD_LBFGS, x.size());
    vector<double> lb;
    vector<double> ub;

    for(unsigned i = 0 ; i < x.size(); i++){
        lb.push_back(lower);
        ub.push_back(upper);
    }
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_max_objective(ObjectiveFxn, in_data);
    double minB;
    opt.optimize(x, minB);
}





vector<vector<double> > Stack_Matrices(vector<MatrixXd> &mats){
    vector<vector<double> > out_stack;

    for(unsigned i = 0; i < mats.size(); i++){
        MatrixXd tempMat = mats[i];
        for(unsigned j = 0; j < tempMat.rows(); j++){
            VectorXd temp = tempMat.row(j);
            out_stack.push_back(eigen2vec(temp));
        }
    }
    return(out_stack);
}


double CalcEuclidean(VectorXd &vec1 , VectorXd &vec2){
    VectorXd diff = vec1 - vec2;
    VectorXd sq = diff.array().square();
    return(sq.sum());
}





void BuildCausalVector(VectorXd& vec2Build , VectorXd& index){
    vec2Build.setZero();
    for(int i = 0; i <index.size(); i++){
        if(index[i] >0){
            vec2Build[index[i]-1] = 1;
        }
    }
}


int NextCombo(VectorXd& c, int k, int n) {
    for (int i= k; --i >= 0;) {
        if (++c[i] <= n-(k-i)) {
            while (++i < k)
                c[i]= c[i-1]+1;
            return 1;
        }
    }
    return 0;
}

void PrintVecXD(VectorXd& vec){
    for(int i = 0; i < vec.size(); i++){
        cout << vec[i] << " ";
    }
    cout << endl;
}


/* OLD FUNCTIONS NO LONGER USED


 void NewPost(VectorXd& Marginal, VectorXd& Zs, VectorXd& Lams, VectorXd& beta, MatrixXd& Aj, MatrixXd& InvLD, MatrixXd& LD, int NC, double& fullLikeli){
    int numsnps = Zs.size();
    double runsum = 0;
    for(int i =0 ; i < Marginal.size(); i++){
        Marginal(i) = -10000000;
    }
    double c_prob = 0;
    double z_prob = 0;
    double sum = 0;

    VectorXd causConfig(Zs.size());
    causConfig.setZero();
    VectorXd causIndex(NC);
    causIndex.setZero();
    VectorXd evalPrior(Zs.size());
    EvalAijs(Aj, beta, evalPrior);

    c_prob = Prob_CNew(evalPrior, beta, causConfig);
    z_prob = EvaluateLogMvn(Zs, causConfig, Lams, InvLD, LD);
    sum = c_prob+z_prob;
    runsum = sum;
    int counter = 1;
    while(NextCombo(causIndex, NC, numsnps+1) == 1){
        BuildCausalVector(causConfig, causIndex);
        c_prob = Prob_CNew(evalPrior, beta, causConfig);
        z_prob = EvaluateLogMvn(Zs, causConfig, Lams, InvLD, LD);
        sum = c_prob+z_prob;
        runsum = LogSum(runsum, sum);
        for(int j = 0; j < causIndex.size(); j++){
            if(causIndex[j] >0){
                Marginal[causIndex[j]-1] = LogSum(Marginal[causIndex[j]-1], sum);
            }

        }
        counter ++;
    }
    for(int f = 0 ; f < Marginal.size(); f++){
        Marginal[f] = Marginal[f]- runsum;
    }
    fullLikeli = fullLikeli+runsum;

}

 void EM_Run_No_Opt(CausalProbs &probabilities ,int iter_max, vector<VectorXd> &Zscores, vector<VectorXd> &Lambdas, VectorXd &beta_int, vector<MatrixXd> &Aijs, vector<MatrixXd> &Sigmas, int numCausal){
    vector<double> beta_run = eigen2vec(beta_int);
    vector<MatrixXd> Invert_LD;
    SigmaInvert(Sigmas, Invert_LD);
    ObjectiveData Opt_in;
    Opt_in.Aijs = Stack_Matrices(Aijs);

    int iterations = 0;
    VectorXd beta_update;

    while(iterations < iter_max){
        Estep(Zscores, Lambdas, beta_int, Aijs, Sigmas, Invert_LD, probabilities,numCausal);
        beta_update = vector2eigen(beta_run);
        cout << beta_update << endl << endl;
        if(CalcEuclidean(beta_update, beta_int) < .05){
            beta_int = beta_update;
            break;
        }
        else{
            beta_int = beta_update;

            iterations++;
        }
    }
}


double Estep(vector<VectorXd> &Zscores, vector<VectorXd> &Lambdas, VectorXd &betas, vector<MatrixXd> &Aijs, vector<MatrixXd> &Sigmas, vector<MatrixXd> &InvSigmas, CausalProbs &E_out, int numberCausal){

    vector<VectorXd> marginal_i;
    VectorXd temp;
    VectorXd exp_temp;
    vector<double> stacker;
    vector<double> stack_temp;
    double fullLikeli = 0;
    for(unsigned i = 0; i < Zscores.size(); i ++){
        VectorXd temp(Zscores[i].size());
        NewPost(temp, Zscores[i], Lambdas[i], betas, Aijs[i], InvSigmas[i], Sigmas[i], numberCausal, fullLikeli);
        exp_temp = temp.array().exp();
        marginal_i.push_back(exp_temp);
        stack_temp =  eigen2vec(exp_temp);
        stacker.insert(stacker.end(), stack_temp.begin(), stack_temp.end());
    }

    E_out.probs_locs = marginal_i;
    E_out.probs_stacked = stacker;


    return(fullLikeli);
}




double EM_Run(CausalProbs &probabilites, int iter_max, vector<VectorXd> &Zscores, vector<VectorXd> &Lambdas, VectorXd &beta_int, vector<MatrixXd> &Aijs, vector<MatrixXd> &Sigmas, int numCausal){
    vector<double> beta_run = eigen2vec(beta_int);
    vector<MatrixXd> Invert_LD;
    SigmaInvert(Sigmas, Invert_LD);
    ObjectiveData Opt_in;
    Opt_in.Aijs = Stack_Matrices(Aijs);

    int iterations = 0;
    VectorXd beta_update;
    double likeliOld = 0;
    double likeli =1;
    while(iterations < iter_max){
        likeli = Estep(Zscores, Lambdas, beta_int, Aijs, Sigmas, Invert_LD, probabilites, numCausal);
        Opt_in.probs = probabilites.probs_stacked;
        void *opt_ptr = &Opt_in;
        Optimize_Nlopt(beta_run, -20, 20, beta_run[0], opt_ptr);
        beta_update = vector2eigen(beta_run);
        cout << beta_update << endl << endl;
        if(abs(likeli - likeliOld) < 0.01){
            beta_int = beta_update;
            break;
        }
        else{
            beta_int = beta_update;
            likeliOld = likeli;
            iterations++;
        }
    }
    return(likeli);
}



inline  double EvaluateLogMvn(VectorXd& Z_vec,  VectorXd C_vec,  VectorXd& Lam_vec,  MatrixXd& Inv_Sig,  MatrixXd& Sig){

    VectorXd elm_wise = C_vec.cwiseProduct(Lam_vec);
    VectorXd mu = Sig*elm_wise;
    VectorXd Z_mu = Z_vec-mu;
    VectorXd RHSeval = Inv_Sig*Z_mu;
    double prob = (-.5)*(Z_mu.dot(RHSeval));

    return(prob);
}



double Prob_C(MatrixXd& Aj,VectorXd& beta, VectorXd& C_vector){
    double probs = 0;
    for(int i = 0; i < C_vector.size(); i++){
        if(C_vector[i] == 1){
            VectorXd Aij_temp = Aj.row(i);
            double temp_prob = Prob_Cijs(beta, Aij_temp);
            probs += log(temp_prob);
        }
        else{
            VectorXd Aij_temp = Aj.row(i);
            double temp_prob = (1- Prob_Cijs(beta, Aij_temp));
            probs += log(temp_prob);
        }
    }
    return(probs);
}


 void GetMultVectors(vector<VectorXd>& output, int locs, string root, string type){
    for(int i = 0 ; i < locs; i++){
        string fname = root +to_string((long long) (i+1) )+ "_" + type;
        output.push_back(GetVector(fname));
    }

}


MatrixXd GetMatrix(string filename){
    ifstream input_file(filename);
    int row, col;
    input_file >> row;
    input_file >> col;

    Eigen::MatrixXd outMat(row, col);
    for(int i = 0; i < row; i++){
        for(int j = 0; j < col; j++){
            input_file >> outMat(i,j);
        }
    }

    input_file.close();
    return(outMat);
}


void GenerateLambdas(vector<VectorXd>& allLams, vector<VectorXd>& allZscores, double minNCP){
    for(unsigned int i =0; i < allZscores.size(); i++){
        VectorXd zTemp =allZscores[i];
        VectorXd lamI(zTemp.size());
        for(int j =0; j <zTemp.size(); j++){
            if(zTemp[j]>0){
                if(zTemp[j] < minNCP){
                    lamI[j] = minNCP;
                }
                else{
                    lamI[j] = zTemp[j];
                }
            }
            else{
                if(zTemp[j] > -1*minNCP){
                    lamI[j] = -1*minNCP;
                }
                else{
                    lamI[j] = zTemp[j];
                }
            }
        }
        allLams.push_back(lamI);
    }
}



void GetCholeskys(vector<MatrixXd>& sigmas, vector<MatrixXd>& chol_factors, vector<MatrixXd>& chol_factors_inv){
    MatrixXd L;
    MatrixXd inv_L;
    for(int i=0; i < sigmas.size(); i++){
        VectorXd X(sigmas[i].rows());
        X.fill(.001);
        MatrixXd Y = X.asDiagonal();
        L = (sigmas[i]+Y).llt().matrixL();
        chol_factors.push_back(L.transpose().triangularView<Eigen::Upper>());
        inv_L = L.inverse();
        //cout << inv_L << endl;
        chol_factors_inv.push_back(inv_L.triangularView<Eigen::Lower>());
    }
}


void CholeskyTransform(vector<VectorXd>& input_vectors, vector<VectorXd>& output_vectors, vector<MatrixXd>& input_matrix){
    for(int i=0; i<input_vectors.size(); i++){
        output_vectors.push_back(input_matrix[i]*input_vectors[i]);
    }
}


 * */