#ifndef EM_CPP_Header_h
#define EM_CPP_Header_h
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen>
#include "nlopt.hpp"
#include <iterator>
#include <string>
#include <iomanip>
#include <omp.h>

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


VectorXd GetVector(string filename);
MatrixXd GetMatrix(string filename);
double Prob_Cijs(VectorXd& beta , VectorXd& aij);
double Prob_C(MatrixXd& Aj,VectorXd& beta, VectorXd& C_vector);
VectorXd Vec2EgVec(vector<double>& original);
double EvaluateLogMvn(const VectorXd& Z_vec, const VectorXd C_vec, const VectorXd& Lam_vec, const MatrixXd& Inv_Sig, const MatrixXd& Sig);

double LogSum(double val1, double val2);
VectorXd Calculate_Posterior4Caus(VectorXd& Zs, VectorXd& Lams, VectorXd& beta, MatrixXd& Aj, MatrixXd& InvLD, MatrixXd& LD);
VectorXd Calculate_Posterior1Caus(VectorXd& Zs, VectorXd& Lams, VectorXd& beta, MatrixXd& Aj, MatrixXd& InvLD);
vector<double> scalar_product(const double& scalar, vector<double>& vec);
vector<double> Vector_Sum(vector<double>& vec1, vector<double>& vec2);
MatrixXd GradientNR(VectorXd& betas,const VectorXd& Marginals, const MatrixXd& Aijs);
MatrixXd HessianNR(VectorXd& betas,const VectorXd& Marginals, const MatrixXd& Aijs);
VectorXd OptimizeNR(VectorXd& betas,const VectorXd& Marginals, const MatrixXd& Aijs);

vector<vector<double> > eigen2Matrix(MatrixXd &mat);
vector<double> eigen2vec(VectorXd &vec);

void ElementWiseMult(const vector<double>& el1, const vector<double>& el2, vector<double>& output);
vector<vector<double> > Stack_Matrices(vector<MatrixXd> &mats);
vector<VectorXd> GetMultVectors(int locs, string root, string type);
vector<MatrixXd> GetMultMatrices(int locs, string root, string type);
void MakeAnnotations(vector<MatrixXd>& annotsRun, vector<MatrixXd>& allAnnots, vector<int>& indices);

void GetMultVectors(vector<VectorXd>& output, int locs, string root, string type);
void GetMultMatrices(vector<MatrixXd>& output, int locs, string root, string type);
void GenerateLambdas(vector<VectorXd>& allLams, vector<VectorXd>& allZscores, double minNCP);
void SigmaInvert(vector<MatrixXd> &Sigma_locs, vector<MatrixXd> &Inv_locs);
void Optimize_Nlopt(vector<double>& x, vector<double> lower, vector<double> upper, double betaZero,  void* in_data);
VectorXd Zscores2Post(VectorXd& Zs);
int NextCombo(VectorXd& c, int k, int n);
void BuildCausalVector(VectorXd& vec2Build , VectorXd& index );
void PrintVecXD(VectorXd& vec);
void EvalAijs(MatrixXd& Aj,VectorXd& beta, VectorXd& priorJ);
double Prob_CNew(VectorXd& priorJ,VectorXd& beta, VectorXd& C_vector);
double GetLogLikeli(VectorXd& Zs, VectorXd& Lams, VectorXd& beta, MatrixXd& Aj, MatrixXd& InvLD, MatrixXd& LD, int NC);
double FullLogLikeli(vector<VectorXd> &Zscores, vector<VectorXd> &Lambdas, VectorXd &betas, vector<MatrixXd> &Aijs, vector<MatrixXd> &Sigmas, vector<MatrixXd> &InvSigmas, int numberCausal);

vector<int> ParseIndex(string& input);

void GetCholeskys(vector<MatrixXd>& sigmas, vector<MatrixXd>& chol_factors, vector<MatrixXd>& chol_factors_inv);
void CholeskyTransform(vector<VectorXd>& input_vectors, vector<VectorXd>& output_vectors, vector<MatrixXd>& input_matrix);

inline  double EvaluateLogMvn_Cholesky(VectorXd& Z_vec,  VectorXd C_vec,  VectorXd& Lam_vec, MatrixXd& upper_chol);
double Estep_chol(vector<VectorXd> &Zscores, vector<VectorXd> &Lambdas, VectorXd &betas, vector<MatrixXd> &Aijs, vector<MatrixXd> &upper_chol, CausalProbs &E_out, int numberCausal);
double NewPost_chol(VectorXd& Marginal, VectorXd& Zs, VectorXd& Lams, VectorXd& beta, MatrixXd& Aj,  MatrixXd& upper_chol, int NC);
double EM_Run_chol(CausalProbs &probabilites, int iter_max, vector<VectorXd> &Zscores, vector<VectorXd> &Lambdas, VectorXd &beta_int, vector<MatrixXd> &Aijs, vector<MatrixXd> &upper_chol, int numCausal);
vector<string> GetInputFiles(string input_directory, string fname, vector<VectorXd>& all_zscores, vector<MatrixXd>& all_LD, vector<MatrixXd>& all_annotations, string  LD_suffix, string annot_suffix);
#endif
