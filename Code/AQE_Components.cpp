#include <Rcpp.h>
#include </Users/dengdetian0603/Downloads/dlib-18.18/dlib/optimization.h>
#include <iostream>
#include <vector>


using namespace std;
using namespace dlib;
using namespace Rcpp;

typedef matrix<double,0,1> column_vector;

// [[Rcpp::export]]
double log_Prob_MSSBSigivenLj(const int& K, NumericVector& MSSi, NumericVector& MBSi, IntegerVector& Lj, 
      NumericVector& ss_tpr, NumericVector& bs_tpr, NumericVector& bs_fpr){
      double log_probs = 0.0;
      for(int k=0; k<K; ++k){
            double logP1 = 0.0;
            if (Lj[k] + MSSi[k] > 0.5) {
                  logP1 =  MSSi[k]*(log(Lj[k]) + Lj[k]*log(ss_tpr[k])) + Lj[k]*(1-MSSi[k])*log(1-ss_tpr[k]);
            }

            double logP2 = MBSi[k]*(Lj[k]*log(bs_tpr[k]) + (1-Lj[k])*log(bs_fpr[k])) + \
                        (1-MBSi[k])*(Lj[k]*log(1-bs_tpr[k]) + (1-Lj[k])*log(1-bs_fpr[k]));
        log_probs += logP1 + logP2;
      }
      return log_probs;
}

// [[Rcpp::export]]
NumericMatrix log_ProbMat_MSSBSgivenL(const int& K, IntegerMatrix& Lall_matrix,  NumericMatrix& MSS, 
                  NumericMatrix& MBS, NumericVector& ss_tpr, NumericVector& bs_tpr, NumericVector& bs_fpr){
      int J = Lall_matrix.nrow();
      int N = MSS.nrow();

      NumericMatrix condition_probs(N, J);
      for(int i=0; i<N; ++i){
            NumericVector mssi = MSS(i,_);
            NumericVector mbsi = MBS(i,_);
            for(int j=0; j<J; ++j){
                  IntegerVector lj = Lall_matrix(j,_);
                  condition_probs(i,j) = log_Prob_MSSBSigivenLj(K, mssi, mbsi, lj, \
                                                      ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr);
            }
      }
      return condition_probs;
}


/* ------------ log Likelihood for Control data ----------- */
// [[Rcpp::export]]
double log_Prob_ctrl(const int& K, const int& Nctrl, NumericMatrix& Mbs, NumericVector& bs_fpr){
      double log_prob_sum = 0.0;

      for(int i=0; i<Nctrl; ++i){
            for(int k=0; k<K; ++k){
                  log_prob_sum += Mbs(i,k)*log(bs_fpr[k]) + (1-Mbs(i,k))*log(1-bs_fpr[k]);
            }
      }
      return log_prob_sum;
}


/* ------------ Likelihood for Case data without GS ----------- */
// double log_Prob_SS_BS(int K, NumericMatrix MSS, NumericMatrix MBS, NumericVector I_GS, NumericVector Pr_Lj, 
//    NumericVector ss_tpr, NumericVector bs_tpr, NumericVector bs_fpr, IntegerMatrix Lall_matrix){
//    int J = Lall_matrix.nrow();
//    int N = I_GS.size();

//    double log_sum_probi = 0.0;
//    for(int i=0; i<N; ++i){
//          NumericVector mssi = MSS(i,_);
//          NumericVector mbsi = MBS(i,_);
            
//          if (I_GS[i]>0.5) {
//                log_sum_probi += 0.0;
//          }
//          else {
//                double prob_i = 0.0;
//                for(int j=0; j<J; ++j){
//                      IntegerVector lj = Lall_matrix(j,_);
//                      double prob_ij = Pr_Lj[j] * exp(log_Prob_MSSBSigivenLj(K=K, mssi, mbsi, lj, \
//                                                 ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr));
//                      prob_i += prob_ij;
//                }
//                log_sum_probi += log(prob_i);
//          }
//    }
//    return log_sum_probi;
// }



// [[Rcpp::export]]
NumericVector ThetaToMu(NumericVector& Theta1, NumericVector& Theta2, IntegerMatrix& LUmat, int& K, int& J1){
      NumericVector potentials(J1);
      NumericVector Mu(K, 0.0);
      double denominatorA = 1.0;
      int K2 = Theta2.size();
      for (int j=0; j<J1; ++j){
            double tmp = 0.0;
            for (int k=0; k<K; ++k){
                  tmp += LUmat(j,k)*Theta1(k);
            }
            for (int kk=0; kk<K2; ++kk){
                  tmp += LUmat(j,kk+K)*Theta2(kk);
            }     
            potentials(j) = exp(tmp);
            denominatorA = denominatorA + potentials(j);
            for (int k=0; k<K; ++k){
                  if (LUmat(j,k)>0) {
                        Mu(k) += potentials(j);
                  }
            }
      }
      for (int k=0; k<K; ++k){
            Mu(k) = Mu(k)/denominatorA;
      }
      return Mu;
}
// ThetaToMu(c(0,0,0), c(0,0,0), cbind(dmat$Lmat, dmat$Umat), 3,7)

NumericVector ThetaToMu(const NumericVector& Theta1, const NumericVector& Theta2, const IntegerMatrix& LUmat, 
                        const int& K, const int& J1) {
      NumericVector potentials(J1);
      NumericVector Mu(K, 0.0);
      double denominatorA = 1.0;
      int K2 = Theta2.size();
      for (int j=0; j<J1; ++j){
            double tmp = 0.0;
            for (int k=0; k<K; ++k){
                  tmp += LUmat(j,k)*Theta1(k);
            }
            for (int kk=0; kk<K2; ++kk){
                  tmp += LUmat(j,kk+K)*Theta2(kk);
            }     
            potentials(j) = exp(tmp);
            denominatorA = denominatorA + potentials(j);
            for (int k=0; k<K; ++k){
                  if (LUmat(j,k)>0) {
                        Mu(k) += potentials(j);
                  }
            }
      }
      for (int k=0; k<K; ++k){
            Mu(k) = Mu(k)/denominatorA;
      }
      return Mu;
}


// [[Rcpp::export]]
NumericVector ThetaToCellProb(NumericVector& Theta1, NumericVector& Theta2, IntegerMatrix& LUmat, 
                              int& K, int& J1, int& Option){ 
      //Option values 1: normalized prob; 2: log normalized probs 3: log_A, LU*Theta[-0] 4: unnormalized potentials 
      NumericVector potentials(J1+1);
      NumericVector LUTHETA(J1+1);
      double denominatorA = 1.0;
      potentials(0) = 1.0;
      LUTHETA(0) = 0;

      int K2 = Theta2.size();
      for (int j=1; j<J1+1; ++j){
            double tmp = 0.0;
            for (int k=0; k<K; ++k){
                  tmp += LUmat(j-1,k)*Theta1(k);
            }
            for (int kk=0; kk<K2; ++kk){
                  tmp += LUmat(j-1,kk+K)*Theta2(kk);
            }
            LUTHETA(j) = tmp;
            potentials(j) = exp(tmp);
            denominatorA = denominatorA + potentials(j);
      }

      if (Option == 1) {
            for (int j=0; j<J1+1; ++j){
                  potentials(j) = potentials(j)/denominatorA;
            }
            return potentials;
      } else if (Option == 2) {
            double log_A = log(denominatorA);
            for (int j=0; j<J1+1; ++j){
                  LUTHETA(j) = LUTHETA(j) - log_A;
            }
            return LUTHETA;
      } else if (Option == 3) {
            double log_A = log(denominatorA);
            LUTHETA(0) = log_A;
            return LUTHETA;
      } else {
            return potentials;
      }
}


/* -------------------------------------------------------------------------------------------- */

class QuadLoss
{
public:
      typedef ::column_vector column_vector;
      QuadLoss(const NumericVector& Mu, const NumericVector& Theta2, 
                  const IntegerMatrix& LUmat, const int& K, const int& J1)
      {
              this->mu = Mu;
              this->theta2 = Theta2;
              this->lumat = LUmat;
              this->K = K;
              this->J1 = J1;
      }

      double operator() ( const column_vector& Theta1) const
      {
            NumericVector x(this->K);
            for (int k=0; k< this->K; ++k){
                  x(k) = Theta1(k);
            }
            double output = 0.0;
            NumericVector tmp;
            tmp = ThetaToMu( x, this->theta2, this->lumat, this->K, this->J1 );
            for (int k=0; k< this->K; ++k){
                  output += pow(tmp(k) - this->mu(k), 2);
            }                 
            return output;    
      }

private:
      NumericVector mu;
      NumericVector theta2;
      IntegerMatrix lumat;
      int K;
      int J1;
};


// double qloss(NumericVector& Theta1, NumericVector& Mu, NumericVector& Theta2, IntegerMatrix& LUmat, int& K, int& J1) {
//       QuadLoss q(Mu, Theta2, LUmat, K, J1);
//       column_vector x(K);
//       for (int k=0; k< K; ++k){
//             x(k) = Theta1(k);
//       }
//       return q(x);
// }


// [[Rcpp::export]]
column_vector SolveForTheta1(NumericVector& Mu, NumericVector& Theta2, IntegerMatrix& LUmat, int& K, int& J1){
      column_vector starting_point(K);
      for (int k=0; k<K; ++k){
            starting_point(k) = 0.0;
      }

      find_min_using_approximate_derivatives(bfgs_search_strategy(),
                                               objective_delta_stop_strategy(1e-8),
                                               QuadLoss(Mu, Theta2, LUmat, K, J1), starting_point, -1, 1e-9);
      return starting_point;
}
// SolveForTheta1(c(0.5,0.5,0.5), c(0,0,0), cbind(dmat$Lmat, dmat$Umat), 3,7)



// [[Rcpp::export]]
NumericMatrix XbetaToProbL_unique(int& K, int& J1, NumericMatrix& X_unique, NumericVector& Beta, NumericVector& Theta2, 
                                    IntegerMatrix& LUmat, int& Option){
      //Option values 1: normalized prob; 2: log normalized probs 3: log_A, LU*Theta[-0] 4: unnormalized potentials 
      int D = Beta.size()/K;
      int Nunique = X_unique.nrow();
      NumericMatrix ProbL_unique(Nunique, J1+1);
      
      for (int u=0; u<Nunique; ++u){
            NumericVector Mu_unique(K);
            for (int k=0; k<K; ++k){
                  double tmp = 0.0;
                  for (int d=0; d<D; ++d){
                        tmp += X_unique(u, d)*Beta(k*D+d);
                  }
                  Mu_unique(k) = 1/(1+exp(-tmp));
            }
            column_vector theta1_u = SolveForTheta1(Mu_unique, Theta2, LUmat, K, J1);
            
            NumericVector Theta1_unique(K);
            for (int k=0; k<K; ++k){
                  Theta1_unique(k) = theta1_u(k);
            }

            ProbL_unique(u,_) = ThetaToCellProb(Theta1_unique, Theta2, LUmat, K, J1, Option);
      }
      return ProbL_unique;
} 
// XbetaToProbL_unique(3, 7, X.unique, as.vector(Data$Beta), c(0,0,0), lumat, 1)
// XbetaToProbL_unique(5, 31, X.unique, as.vector(Data$Beta), theta2, LUmat, 1)        


double log_Prior_BetaTheta2(NumericVector& Beta, NumericVector& Theta2, 
                              const double& var_beta, const double& var_theta2, const double& mean_theta2)
{
      double result = 0.0;
      for (int k=0; k<Beta.size(); ++k){
            result -= (Beta(k)*Beta(k))/(2*var_beta);
      }
      for (int k=0; k<Theta2.size(); ++k){
            result -= pow(Theta2(k)-mean_theta2, 2)/(2*var_theta2);
      }
      return result;
}


/* --------------------------------------- EM components --------------------------------------------*/


// [[Rcpp::export]]
NumericMatrix EM_UpdateRates(int& K, int& J1, int& N_case, int& N_ctrl, IntegerMatrix& MSS, IntegerMatrix& MBS_case, 
                        IntegerMatrix & MBS_ctrl, NumericMatrix& weights, IntegerVector& X_index, IntegerMatrix& Lmat_withZero,
                        NumericVector& aa, NumericVector& bb, NumericVector& cc, NumericVector& dd, NumericVector& ee, 
                        NumericVector& ff){
      //X_index starts from 0
      NumericMatrix PosRates(3, K);
      for (int k=0; k<K; ++k){
            double A_k = aa(k)-1; double B_k = bb(k)-1;
            double C_k = cc(k)-1; double D_k = dd(k)-1;
            double E_k = ee(k)-1; double F_k = ff(k)-1;
            for (int i=0; i<N_case; ++i){
                  for (int j=0; j<J1+1; ++j){
                        double WijLjk = weights(i,j)*Lmat_withZero(j,k);
                        double Wij_notLjk = weights(i,j)*(1-Lmat_withZero(j,k));
                        A_k += WijLjk * MSS(i,k);
                        B_k += WijLjk * (1-MSS(i,k));
                        C_k += WijLjk * MBS_case(i,k);
                        D_k += WijLjk * (1-MBS_case(i,k));
                        E_k += Wij_notLjk * MBS_case(i,k);
                        F_k += Wij_notLjk * (1-MBS_case(i,k));
                  }
            }
            for (int i=0; i<N_ctrl; ++i){
                  E_k += MBS_ctrl(i,k);
                  F_k += 1-MBS_ctrl(i,k);
            }
            PosRates(0,k) = A_k/(A_k + B_k);
            PosRates(1,k) = C_k/(C_k + D_k);
            PosRates(2,k) = E_k/(E_k + F_k);
      }
      return PosRates;
}




// /* ------------------------------------------------------------------------------------------ */
column_vector SolveForTheta1(NumericVector& Mu, NumericVector& Theta2, const IntegerMatrix& LUmat, 
                              const int& K, const int& J1)
{
      column_vector starting_point(K);
      for (int k=0; k<K; ++k){
            starting_point(k) = 0.0;
      }
      find_min_using_approximate_derivatives(bfgs_search_strategy(),
                                               objective_delta_stop_strategy(1e-8),
                                               QuadLoss(Mu, Theta2, LUmat, K, J1), starting_point, -1);
      return starting_point;
}

NumericVector ThetaToCellProb(NumericVector& Theta1, NumericVector& Theta2, const IntegerMatrix& LUmat, 
                              const int& K, const int& J1, const int& Option)
{ 
      //Option values 1: normalized prob; 2: log normalized probs 3: log_A, LU*Theta[-0] 4: unnormalized potentials 
      NumericVector potentials(J1+1);
      NumericVector LUTHETA(J1+1);
      double denominatorA = 1.0;
      potentials(0) = 1.0;
      LUTHETA(0) = 0;

      int K2 = Theta2.size();
      for (int j=1; j<J1+1; ++j){
            double tmp = 0.0;
            for (int k=0; k<K; ++k){
                  tmp += LUmat(j-1,k)*Theta1(k);
            }
            for (int kk=0; kk<K2; ++kk){
                  tmp += LUmat(j-1,kk+K)*Theta2(kk);
            }
            LUTHETA(j) = tmp;
            potentials(j) = exp(tmp);
            denominatorA = denominatorA + potentials(j);
      }

      if (Option == 1) {
            for (int j=0; j<J1+1; ++j){
                  potentials(j) = potentials(j)/denominatorA;
            }
            return potentials;
      } else if (Option == 2) {
            double log_A = log(denominatorA);
            for (int j=0; j<J1+1; ++j){
                  LUTHETA(j) = LUTHETA(j) - log_A;
            }
            return LUTHETA;
      } else if (Option == 3) {
            double log_A = log(denominatorA);
            LUTHETA(0) = log_A;
            return LUTHETA;
      } else {
            return potentials;
      }
}


NumericMatrix XbetaToProbL_unique(const int& K, const int& J1, const NumericMatrix& X_unique, 
                                    NumericVector& Beta, NumericVector& Theta2, 
                                    const IntegerMatrix& LUmat, const int& Option)
{//Option values 1: normalized prob; 2: log normalized probs 3: log_A, LU*Theta[-0] 4: unnormalized potentials 
      int D = Beta.size()/K;
      int Nunique = X_unique.nrow();
      NumericMatrix ProbL_unique(Nunique, J1+1);
      
      for (int u=0; u<Nunique; ++u){
            NumericVector Mu_unique(K);
            for (int k=0; k<K; ++k){
                  double tmp = 0.0;
                  for (int d=0; d<D; ++d){
                        tmp += X_unique(u, d)*Beta(k*D+d);
                  }
                  Mu_unique(k) = 1/(1+exp(-tmp));
            }

            // the following step is unstable.
            column_vector theta1_u = SolveForTheta1(Mu_unique, Theta2, LUmat, K, J1);
            
            NumericVector Theta1_unique(K);
            for (int k=0; k<K; ++k){
                  Theta1_unique(k) = theta1_u(k);
            }

            ProbL_unique(u,_) = ThetaToCellProb(Theta1_unique, Theta2, LUmat, K, J1, Option);
      }
      return ProbL_unique;
}  




class Qfunction
{
public:
      typedef ::column_vector column_vector;
      Qfunction(const NumericMatrix& Weights, const IntegerVector& X_index, const NumericMatrix& X_unique,
                  const IntegerMatrix& LUmat, const int& K, const int& J1, const int& D,
                  const double& prior_var_beta, const double& prior_var_theta2, const double& prior_mean_theta2)
      {
              this->weights = Weights;
              this->x_index = X_index; // index starts from 0
              this->lumat = LUmat;
              this->x_unique = X_unique;
              this->K = K;
              this->J1 = J1;
              this->D = D; // number of covariates + 1(intercept)
              this->Beta_size = D*K;
              this->Theta2_size = LUmat.ncol()-K;
              this->vbeta = prior_var_beta;
              this->vtheta2 = prior_var_theta2;
              this->mtheta2 = prior_mean_theta2;
      }

      double operator() ( const column_vector& Beta_Theta2) const
      {                  
            NumericVector Beta(this->Beta_size);
            for (int k=0; k<this->Beta_size; ++k){
                  Beta(k) = Beta_Theta2(k);
            }
            NumericVector Theta2(this->Theta2_size);
            for (int k=0; k<this->Theta2_size; ++k){
                  Theta2(k) = Beta_Theta2(this->Beta_size + k);
            }
            NumericMatrix qUnique(this->weights.nrow(), this->J1+1);

            double Qfunc_1 = -1e+26;
            int iter = 0;
            while( (abs(Qfunc_1) > 1e+25) && iter <10){ 
                  // unstable step to get qUnique   
                  try {      
                        qUnique = XbetaToProbL_unique(this->K, this->J1, this->x_unique, Beta, Theta2, this->lumat, 3);
                  }
                  catch(...) {
                        //cout << "Error in solving Theta1." << endl; 
                        ++iter; 
                        continue;
                  }

                  Qfunc_1 = 0.0;
                  for (int i=0; i<this->x_index.size(); ++i){
                        Qfunc_1 -= qUnique(this->x_index(i),0);
                        for (int j=1; j<this->J1+1; ++j){
                              Qfunc_1 += qUnique(this->x_index(i),j) * this->weights(i, j);
                        }
                  }
                  //cout << Qfunc_1 <<endl; 
                  ++iter;
            }
            return Qfunc_1 + log_Prior_BetaTheta2(Beta, Theta2, vbeta, vtheta2, mtheta2);    
      }

private:
      NumericMatrix weights;
      IntegerVector x_index;
      IntegerMatrix lumat;
      NumericMatrix x_unique;
      int K;
      int J1;
      int D;
      int Beta_size;
      int Theta2_size;
      double vbeta;
      double vtheta2;
      double mtheta2;
};


// [[Rcpp::export]]
double qfunc(NumericMatrix& Weights, IntegerVector& X_index, NumericMatrix& X_unique,
                  IntegerMatrix& LUmat, int& K, int& J1, int& D,
                  double& prior_var_beta, double& prior_var_theta2, double& prior_mean_theta2,
                  NumericVector& initial_value) {
      Qfunction q(Weights, X_index, X_unique, LUmat, K, J1, D,
                      prior_var_beta, prior_var_theta2, prior_mean_theta2);

      column_vector x(initial_value.size());
      for (int k=0; k< initial_value.size(); ++k){
            x(k) = initial_value(k);
      }
      return q(x);
}


// [[Rcpp::export]]
NumericMatrix EM_GetWeights(int& K, IntegerMatrix& Lmat_withZero,  NumericMatrix& MSS, NumericMatrix& MBS,
                              NumericVector& ss_tpr, NumericVector& bs_tpr, NumericVector& bs_fpr,
                              IntegerVector& X_index, NumericMatrix& X_unique, IntegerMatrix& LUmat,
                              NumericVector& Beta, NumericVector& Theta2)
{
      int J = Lmat_withZero.nrow();
      int N = MSS.nrow();

      NumericMatrix log_ProbMat_a(N,J);
      NumericMatrix unnormed_ProbMat_b_unique(X_unique.nrow(),J);

      log_ProbMat_a = log_ProbMat_MSSBSgivenL(K, Lmat_withZero, MSS, MBS, ss_tpr, bs_tpr, bs_fpr);

      int iter = 0;
      while (iter < 10){
            try {
                  unnormed_ProbMat_b_unique = XbetaToProbL_unique(K, J-1, X_unique, Beta, Theta2, LUmat, 4);
            }
            catch(...) {cout << "Error getting ProbMat_b." << endl; ++iter; continue;}
            break;
      }

      NumericMatrix weights(N, J);
      for (int i=0; i<N; i++){
            double denominatorA = 0.0;
            for (int j=0; j<J; ++j){
                  weights(i,j) = exp(log_ProbMat_a(i,j)) * unnormed_ProbMat_b_unique(X_index(i), j);
                  denominatorA += weights(i,j);
            }
            for (int j=0; j<J; ++j){
                  weights(i,j) = weights(i,j)/denominatorA;
            }
      }
      return weights;
}


// [[Rcpp::export]]
column_vector EM_UpdateBetaTheta2(NumericMatrix& Weights, IntegerVector& X_index, NumericMatrix& X_unique,
                  IntegerMatrix& LUmat, int& K, int& J1, int& D,
                  double& prior_var_beta, double& prior_var_theta2, double& prior_mean_theta2,
                  NumericVector& initial_value)
{
      int Par_size = initial_value.size();
      column_vector starting_point(Par_size);
      for (int k=0; k<Par_size; ++k){
            starting_point(k) = initial_value(k);
      }

      find_max_using_approximate_derivatives(bfgs_search_strategy(), objective_delta_stop_strategy(1e-8),
            Qfunction(Weights, X_index, X_unique, LUmat, K, J1, D,
                      prior_var_beta, prior_var_theta2, prior_mean_theta2), 
            starting_point, 0.0, 1e-9);
      return starting_point;
}


/* ------------------------------------------------------------------------------------------------------- */

class log_full_distn
{
public:
      typedef ::column_vector column_vector;
      log_full_distn( const IntegerVector& X_index, const NumericMatrix& X_unique,
                  const IntegerMatrix& LUmat, const int& K, const int& J1, const int& D, IntegerMatrix& L_withZero_mat,
                  NumericMatrix& MSS, NumericMatrix& MBS, NumericMatrix& MBS_ctrl,
                  const double& prior_var_beta, const double& prior_var_theta2, const double& prior_mean_theta2, 
                  const NumericVector& aa, const NumericVector& bb, const NumericVector& cc, const NumericVector& dd, 
                  const NumericVector& ee, const NumericVector& ff)
      {
              this->x_index = X_index; // index starts from 0
              this->lumat = LUmat;
              this->x_unique = X_unique;
              this->K = K;
              this->J1 = J1;
              this->D = D; // number of covariates + 1(intercept)
              this->Lmat_withZero = L_withZero_mat;
              this->mss = MSS;
              this->mbs = MBS;
              this->mbs_ctrl = MBS_ctrl;

              this->Beta_size = D*K;
              this->Theta2_size = LUmat.ncol()-K;
              this->vbeta = prior_var_beta;
              this->vtheta2 = prior_var_theta2;
              this->mtheta2 = prior_mean_theta2;

              this->aa = aa;
              this->bb = bb;
              this->cc = cc;
              this->dd = dd;
              this->ee = ee;
              this->ff = ff;
      }

      double operator() ( const column_vector& All_Par)
      {                  
            NumericVector ss_tpr(this->K);
            for (int k=0; k<this->K; ++k){
                  ss_tpr(k) = All_Par(k);
            }
            NumericVector bs_tpr(this->K);
            for (int k=0; k<this->K; ++k){
                  bs_tpr(k) = All_Par(k+this->K);
            }
            NumericVector bs_fpr(this->K);
            for (int k=0; k<this->K; ++k){
                  bs_fpr(k) = All_Par(k+this->K*2);
            }

            NumericVector Beta(this->Beta_size);
            for (int k=0; k<this->Beta_size; ++k){
                  Beta(k) = All_Par(k+this->K*3);
            }
            NumericVector Theta2(this->Theta2_size);
            for (int k=0; k<this->Theta2_size; ++k){
                  Theta2(k) = All_Par(this->Beta_size + this->K*3 + k);
            }

            NumericMatrix log_ProbMat_MgivenL(this->mss.nrow(), this->J1+1);
            log_ProbMat_MgivenL = log_ProbMat_MSSBSgivenL(this->K, this->Lmat_withZero,  this->mss, this->mbs,
                                    ss_tpr, bs_tpr, bs_fpr);


            NumericMatrix Prob_L_unique(this->mss.nrow(), this->J1+1);

            double log_case = -1e+26;
            int iter = 0;
            while( (abs(log_case) > 1e+25) && iter <10){ 
                  // unstable step to get qUnique   
                  try {      
                        Prob_L_unique = XbetaToProbL_unique(this->K, this->J1, this->x_unique, Beta, Theta2, this->lumat, 1);
                  }
                  catch(...) {
                        //cout << "Error in solving Theta1." << endl; 
                        ++iter; 
                        continue;
                  }

                  log_case = 0.0;
                  for (int i=0; i<this->x_index.size(); ++i){
                        double lik_i = 0;
                        for (int j=0; j<this->J1+1; ++j){
                              lik_i += Prob_L_unique(this->x_index(i),j) * exp(log_ProbMat_MgivenL(i, j));
                        }
                        log_case += log(lik_i);
                  }
                  ++iter;
            }

            double rates_logprior = 0.0;
            for (int k=0; k< this->K; ++k){
                  rates_logprior += (this->aa(k)-1)*log(ss_tpr(k)) + (this->bb(k)-1)*log(1-ss_tpr(k)) + 
                                    (this->cc(k)-1)*log(bs_tpr(k)) + (this->dd(k)-1)*log(1-bs_tpr(k)) + 
                                    (this->ee(k)-1)*log(bs_fpr(k)) + (this->ff(k)-1)*log(1-bs_fpr(k));
            }
            return log_case + log_Prob_ctrl(this->K, this->mbs_ctrl.nrow(), this->mbs_ctrl, bs_fpr)+ 
                    rates_logprior + log_Prior_BetaTheta2(Beta, Theta2, vbeta, vtheta2, mtheta2);    
      }

private:
      IntegerVector x_index;
      IntegerMatrix lumat;
      NumericMatrix x_unique;
      IntegerMatrix Lmat_withZero;
      NumericMatrix mss;
      NumericMatrix mbs;
      NumericMatrix mbs_ctrl;
      NumericMatrix ss_tpr;
      NumericMatrix bs_tpr;
      NumericMatrix bs_fpr;
      int K;
      int J1;
      int D;
      int Beta_size;
      int Theta2_size;
      double vbeta;
      double vtheta2;
      double mtheta2;
      NumericVector aa;
      NumericVector bb;
      NumericVector cc;
      NumericVector dd;
      NumericVector ee;
      NumericVector ff;
};


// [[Rcpp::export]]
double Evaluate_log_distn(NumericVector& all_par, const IntegerVector& X_index, const NumericMatrix& X_unique,
                  const IntegerMatrix& LUmat, const int& K, const int& J1, const int& D, IntegerMatrix& L_withZero_mat,
                  NumericMatrix& MSS, NumericMatrix& MBS, NumericMatrix& MBS_ctrl,
                  const double& prior_var_beta, const double& prior_var_theta2, const double& prior_mean_theta2, 
                  const NumericVector& aa, const NumericVector& bb, const NumericVector& cc, const NumericVector& dd, 
                  const NumericVector& ee, const NumericVector& ff)
{     
      log_full_distn logD( X_index, X_unique, LUmat, K, J1, D, L_withZero_mat, 
                   MSS, MBS, MBS_ctrl, prior_var_beta, prior_var_theta2, prior_mean_theta2, 
                   aa, bb, cc, dd, ee, ff);

      column_vector x(all_par.size());
      for (int k=0; k< all_par.size(); ++k){
            x(k) = all_par(k);
      }
      return logD(x);
}


// column_vector Maximize_Posterior(NumericVector initial_value, const IntegerVector& X_index, const NumericMatrix& X_unique,
//                   const IntegerMatrix& LUmat, const int& K, const int& J1, const int& D, IntegerMatrix& L_withZero_mat,
//                   NumericMatrix& MSS, NumericMatrix& MBS, NumericMatrix& MBS_ctrl,
//                   const double& prior_var_beta, const double& prior_var_theta2, const double& prior_mean_theta2, 
//                   const NumericVector& aa, const NumericVector& bb, const NumericVector& cc, const NumericVector& dd, 
//                   const NumericVector& ee, const NumericVector& ff)
// {
//       int Par_size = initial_value.size();
//       column_vector starting_point(Par_size);
//       for (int k=0; k<Par_size; ++k){
//             starting_point(k) = initial_value(k);
//       }

//       find_max_using_approximate_derivatives(bfgs_search_strategy(), objective_delta_stop_strategy(1e-8),
//             log_full_distn(X_index, X_unique, LUmat, K, J1, D, L_withZero_mat, 
//                    MSS, MBS, MBS_ctrl, prior_var_beta, prior_var_theta2, prior_mean_theta2, 
//                    aa, bb, cc, dd, ee, ff), 
//             starting_point, 0.0, 1e-9);
//       return starting_point;
// }




