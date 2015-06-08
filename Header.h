//
//  Header.h
//  PAINTOR.1.2
//
//  Created by Gleb Kichaev on 5/2/15.
//  Copyright (c) 2015 Gleb Kichaev. All rights reserved.
//

#ifndef PAINTOR_1_2_Header_h
#define PAINTOR_1_2_Header_h


#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen>
#include "nlopt.hpp"
#include <iterator>
#include <string>
#include <iomanip>
#include <algorithm>

using namespace Eigen;
using namespace std;
using namespace nlopt;

struct CausalProbs{
    vector<VectorXd> probs_locs;
    vector<double> probs_stacked;
};


struct ObjectiveData{
    vector<double> probs;
    vector<vector<double> > Aijs;
};


bool Check_Locus_File(vector<VectorXd>& association_stats, vector<MatrixXd>& ld_matrices, MatrixXd& annotation_matrix, string& locus_name);
vector<string> &split(const string &s, char delim, vector<string> &elems) ;
vector<string> split(const string &s, char delim);
void Read_Locus(string &input_directory, string& fname, vector<string> & zname,  vector<VectorXd>& all_association_stat, vector<string>& snp_info, string & header);
double Regularize_LD(MatrixXd & ld_mat);
void Read_LD(string &input_directory, string& fname , MatrixXd& chol_factor, MatrixXd& chol_factor_inv);
void Read_Annotations(string &input_directory, string& fname , vector<string>& model_annotations, MatrixXd& out_annotations, vector<string> & annotation_header_split);
void Generate_Lambda(VectorXd& Zscore, VectorXd& Lambda, double minNCP);
bool Check_Locus_File(vector<VectorXd>& association_stats, vector<MatrixXd>& ld_matrices, MatrixXd& annotation_matrix, string& locus_name);

void Get_All_Input(string& file_list_name, string& directory, vector<string> znames, vector<string>& model_annotations, vector<vector<VectorXd>>& all_transformed_statistics, vector<vector<VectorXd>>& all_lambdas ,vector<vector<MatrixXd>> &all_upper_cholesky,  vector<MatrixXd>& all_annotations, vector<string>& LD_suffix, string& annotation_suffix,
                   vector<vector<string>>& all_SNP_info, vector<string>& all_headers);


vector<double> eigen2vec(VectorXd &vec);
VectorXd vector2eigen(vector<double> &vec);
vector<vector<double> > eigen2Matrix(MatrixXd &mat);
inline double Prob_Cijs(VectorXd& beta , VectorXd& aij);
double Prob_C(MatrixXd& Aj,VectorXd& beta, VectorXd& C_vector);
double EM_Run_chol(CausalProbs &probabilites, int iter_max, vector<vector<VectorXd>> &Zscores, vector<vector<VectorXd>> &Lambdas, VectorXd &beta_int, vector<MatrixXd> &Aijs, vector<vector<MatrixXd>> &upper_chol, int numCausal);
double Estep_chol(vector<vector<VectorXd>> &Zscores, vector<vector<VectorXd>> &Lambdas, VectorXd &betas, vector<MatrixXd> &Aijs, vector<vector<MatrixXd>> &upper_chol, CausalProbs &E_out, int numberCausal);
vector<vector<double> > Stack_Matrices(vector<MatrixXd> &mats);
double CalcEuclidean(VectorXd &vec1 , VectorXd &vec2);
void EvalAijs(MatrixXd& Aj,VectorXd& beta, VectorXd& priorJ);
double Prob_CNew(VectorXd& priorJ,VectorXd& beta, VectorXd& C_vector);

void NewPost_chol(VectorXd& Marginal, vector<VectorXd>& Zs, vector<VectorXd>& Lams, VectorXd& beta, MatrixXd& Aj,  vector<MatrixXd>& upper_chol, int NC, double& fullLikeli);
inline  double EvaluateLogMvn_Cholesky(VectorXd& Z_vec,  VectorXd& C_vec,  VectorXd& Lam_vec, MatrixXd& upper_chol);
void BuildCausalVector(VectorXd& vec2Build , VectorXd& index);
int NextCombo(VectorXd& c, int k, int n) ;
inline double LogSum(double val1, double val2);
double dot_prod(vector<double>& vec1, vector<double>& vec2);
vector<double> scalar_product(const double& scalar, vector<double>& vec);
vector<double> Vector_Sum(vector<double>& vec1, vector<double>& vec2);
vector<double> GradientFxn(vector<double>& betas, ObjectiveData *in_data);
double ObjectiveFxn(const vector<double> &x, vector<double> &grad, void *data);
void Optimize_Nlopt(vector<double>& x, double lower, double upper, double betaZero,  void* in_data);
VectorXd Zscores2Post(VectorXd& Zs);

void NewPost_chol(VectorXd& Marginal, vector<VectorXd>& Zs, vector<VectorXd>& Lams, VectorXd& beta, MatrixXd& Aj,  vector<MatrixXd>& upper_chol, int NC, double& fullLikeli);


double Estep_chol(vector<vector<VectorXd>> &Zscores, vector<vector<VectorXd>> &Lambdas, VectorXd &betas, vector<MatrixXd> &Aijs, vector<vector<MatrixXd>>& upper_chol, CausalProbs &E_out, int numberCausal);
void NewPost(VectorXd& Marginal, VectorXd& Zs, VectorXd& Lams, VectorXd& beta, MatrixXd& Aj, MatrixXd& InvLD, MatrixXd& LD, int NC, double& fullLikeli);
void Write_Posterior(string& out_dir, string & out_name, VectorXd& locus_results, vector<string>& locus_info, string& header );
void Write_All_Output(string& input_files, string& out_dir, string& out_suffix, CausalProbs& results, vector<vector<string>> & all_locus_info, VectorXd & gamma_estimates, string& Gname, double log_likeli, string & Lname, vector<string>& all_headers, vector<string>& annot_names);
#endif
