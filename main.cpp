#include "Header.h"

int main(int argc, const char * argv[])
{
    
    string root = "./";
    int number_loci = 2;
    int mallerFlag = 0;
    int OptBetaFlag = 0;
    int causalSNPs = 2;
    string outName = "_PaintorProb.txt";
    string gammaName = "EstimatedGamma.txt";
    int maxIter = 10;
    string outdir;
    string annotationIndex;
    string z_name = "Zscore.txt";
    string aij_name = "Aij.txt";
    string ld_name = "LD.txt";
    string likeli_name = "Likelihood.txt";
    vector<int> annot_Index;
    string beta_init = "NoInput";

    for(int i = 1; i < argc; i++){
        string argComp = argv[i];
        
        if(argComp.compare("-d") == 0){
            root = argv[i+1];
            string last_char = string(&root[root.size()-1]);
            if(last_char.compare("/")!=0){
                root = root + "/";
            }
        }
        else if(argComp.compare("-o") == 0){
            outdir = argv[i+1];
            string last_char = string(&outdir[outdir.size()-1]);
            if(last_char.compare("/")!=0){
                outdir = outdir + "/";
            }
        }
        else if(argComp.compare("-l") == 0){
            number_loci = stoi(argv[i+1]);
        }
        else if(argComp.compare("-i") == 0){
            annotationIndex = argv[i+1];
            annot_Index = ParseIndex(annotationIndex);
        }
        else if(argComp.compare("-c") == 0){
            causalSNPs = stoi(argv[i+1]);
        }
        else if(argComp.compare("-Gname") == 0){
            gammaName = argv[i+1];
        }
        else if(argComp.compare("-Lname") == 0){
            likeli_name = argv[i+1];
        }
        else if(argComp.compare("-m") == 0){
            mallerFlag = 1;
        }
        else if(argComp.compare("-MI") == 0){
            maxIter = stoi(argv[i+1]);
        }
        else if(argComp.compare("-Zname") == 0){
            z_name = argv[i+1];
        }
        else if(argComp.compare("-LDname") == 0){
            ld_name = argv[i+1];
        }
        else if(argComp.compare("-ANname") == 0){
            aij_name = argv[i+1];
        }
        else if(argComp.compare("-OUTname") == 0){
            outName = argv[i+1];
        }
	else if(argComp.compare("-GAMinitial") == 0){
            beta_init = argv[i+1];
        }
    }
    
    vector<VectorXd> Z_scores;
    vector<VectorXd> Lambdas;
    GetMultVectors(Z_scores, number_loci, root, z_name);
    GenerateLambdas(Lambdas, Z_scores, 3.7);
    
    vector<MatrixXd> Aij_Locs ;
    vector<MatrixXd> LD_locs;
    GetMultMatrices(Aij_Locs, number_loci, root, aij_name);
    GetMultMatrices(LD_locs,number_loci, root, ld_name);

    vector<MatrixXd> Aij_run;


    MakeAnnotations(Aij_run, Aij_Locs, annot_Index);
    VectorXd beta_est(Aij_run[0].cols());
    beta_est.setZero();
   
    if(beta_init.compare("NoInput")!= 0){
        beta_est = GetVector(root + beta_init);
    }
    
    if(mallerFlag == 1){
        for(int i = 0; i < number_loci; i++){
            VectorXd post1caus = Zscores2Post(Z_scores[i]);
            ofstream myfile;
            string fname = root + to_string((long long)(i+1)) + "_PosteriorMaller.txt";
            myfile.open(fname);
            myfile << post1caus;
            myfile.close();
        }
    }
    
    if(mallerFlag != 1){
        CausalProbs runProbs;
        if(OptBetaFlag == 0){
            double loglikelihod = EM_Run(runProbs, maxIter, Z_scores, Lambdas, beta_est, Aij_run, LD_locs, causalSNPs);
            for(int i = 0; i < (int) runProbs.probs_locs.size(); i++){
                ofstream myfile;
                string fname = outdir + to_string((long long)(i+1)) + outName;
                myfile.open(fname);
                myfile << runProbs.probs_locs[i];
                myfile << "\n";
                myfile.close();
            }
            
            ofstream myfile;
            myfile.open(outdir+likeli_name);
            myfile <<std::setprecision(9) << loglikelihod;
            myfile << "\n";
            myfile.close();
        }
        
        else{
            EM_Run_No_Opt(runProbs, maxIter, Z_scores, Lambdas, beta_est, Aij_Locs, LD_locs, causalSNPs);
            for(int i = 0; i < (int) runProbs.probs_locs.size(); i++){
                ofstream myfile;
                string fname = outdir + to_string((long long)(i+1)) + outName;
                myfile.open(fname);
                myfile << runProbs.probs_locs[i];
                myfile << "\n";
                myfile.close();
            }
        }
        
        ofstream myfile;
        myfile.open(outdir + gammaName);
        myfile << beta_est;
        myfile << "\n";
        myfile.close();
    }
    
return 0;

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


VectorXd GetVector(string filename){
    
    ifstream input_file(filename);
    istream_iterator<double> start(input_file), end;
    vector<double> vec(start, end);
    input_file .close();
    
    VectorXd outvec(vec.size());
    for(unsigned i = 0; i < vec.size(); i++){
        outvec(i) = vec[i];
    }
    
    return(outvec);
    
}

void GetMultVectors(vector<VectorXd>& output, int locs, string root, string type){
    for(int i = 0 ; i < locs; i++){
        string fname = root +to_string((long long) (i+1) )+ "_" + type;
        output.push_back(GetVector(fname));
    }

}

void GetMultMatrices(vector<MatrixXd>& output, int locs, string root, string type){
    for(int i = 0 ; i < locs; i++){
        string fname = root + to_string((long long)(i+1)) + "_" + type;
        output.push_back(GetMatrix(fname));
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



inline  double EvaluateLogMvn(VectorXd& Z_vec,  VectorXd C_vec,  VectorXd& Lam_vec,  MatrixXd& Inv_Sig,  MatrixXd& Sig){
    
    VectorXd elm_wise = C_vec.cwiseProduct(Lam_vec);
    VectorXd mu = Sig*elm_wise;
    VectorXd Z_mu = Z_vec-mu;
    VectorXd RHSeval = Inv_Sig*Z_mu;
    double prob = (-.5)*(Z_mu.dot(RHSeval));
    
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


void SigmaInvert(vector<MatrixXd> &Sigma_locs, vector<MatrixXd> &Inv_locs){
    for(unsigned i = 0; i < Sigma_locs.size(); i++){
        MatrixXd temp = Sigma_locs[i];
        for(unsigned j = 0; j < Sigma_locs[i].cols(); j++){
            
            temp(j,j) += 0.001;
        }
        Sigma_locs[i] = temp;
    }
	for(unsigned i = 0; i < Sigma_locs.size(); i++){
        FullPivLU<MatrixXd> ldt(Sigma_locs[i]);
        Inv_locs.push_back(ldt.inverse());
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
    
    cout << std::setprecision(9) << fullLikeli << endl;
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
