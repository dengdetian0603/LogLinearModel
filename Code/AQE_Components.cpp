#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
double log_Prob_MSSBSigivenLj(int K, NumericVector MSSi, NumericVector MBSi, IntegerVector Lj, 
	NumericVector ss_tpr, NumericVector bs_tpr, NumericVector bs_fpr){
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
NumericMatrix log_ProbMat_MSSBSgivenL(int K, IntegerMatrix Lall_matrix,  NumericMatrix MSS, NumericMatrix MBS,
	NumericVector ss_tpr, NumericVector bs_tpr, NumericVector bs_fpr){
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





/* ------------ Likelihood for Case data without GS ----------- */
// double log_Prob_SS_BS(int K, NumericMatrix MSS, NumericMatrix MBS, NumericVector I_GS, NumericVector Pr_Lj, 
// 	NumericVector ss_tpr, NumericVector bs_tpr, NumericVector bs_fpr, IntegerMatrix Lall_matrix){
// 	int J = Lall_matrix.nrow();
// 	int N = I_GS.size();

// 	double log_sum_probi = 0.0;
// 	for(int i=0; i<N; ++i){
// 		NumericVector mssi = MSS(i,_);
// 		NumericVector mbsi = MBS(i,_);
		
// 		if (I_GS[i]>0.5) {
// 			log_sum_probi += 0.0;
// 		}
// 		else {
// 			double prob_i = 0.0;
// 			for(int j=0; j<J; ++j){
// 				IntegerVector lj = Lall_matrix(j,_);
// 				double prob_ij = Pr_Lj[j] * exp(log_Prob_MSSBSigivenLj(K=K, mssi, mbsi, lj, \
//                                                 ss_tpr=ss_tpr, bs_tpr=bs_tpr, bs_fpr=bs_fpr));
// 				prob_i += prob_ij;
// 			}
// 			log_sum_probi += log(prob_i);
// 		}
// 	}
// 	return log_sum_probi;
// }

/* ------------ log Likelihood for Control data ----------- */
// [[Rcpp::export]]
double log_Prob_ctrl(int K, int Nctrl, NumericMatrix Mbs, NumericVector bs_fpr){
	double log_prob_sum = 0.0;

	for(int i=0; i<Nctrl; ++i){
		for(int k=0; k<K; ++k){
			log_prob_sum += Mbs(i,k)*log(bs_fpr[k]) + (1-Mbs(i,k))*log(1-bs_fpr[k]);
		}
	}
	return log_prob_sum;
}