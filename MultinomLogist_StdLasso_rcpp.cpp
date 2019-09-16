#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double softmax(double a, double delta){
  double val=std::max(0.0,fabs(a)-delta);
  double sign=0;
  if (a!=0.0){
    sign=a/fabs(a);
  }
  return(sign*val);
}

// [[Rcpp::export]]
List multinomial_logistic_lasso_vecteur_optimal(NumericVector beta,
                                                NumericMatrix X,
                                                NumericVector Y,
                                                int pp,
                                                int L,
                                                int n,
                                                NumericVector P,
                                                NumericVector beta_resamp,
                                                NumericVector XY,
                                                NumericVector Xbeta,
                                                NumericVector E,
                                                NumericVector S,
                                                NumericVector S0,
                                                NumericVector B,
                                                List NN,
                                                NumericVector lambda,
                                                double epsilon,
                                                int max_iter,
                                                bool reference,
                                                double tolerance,
                                                double decay_rate,
                                                double decay_min
) {
  //  << std::numeric_limits<int>::max()
  double incr=std::numeric_limits<double>::max();
  double grad=0.0;
  
  for (int iter=0; iter<max_iter; ++iter){
    double norme2_beta_diff=0;
    double norme2_beta_old=0;
    
    for (int l=0; l<pp; ++l){
      double beta_old=beta[l];
      
      if ((beta_old!=0) || (R::runif(0,1) <= beta_resamp[l])){
        double XdotP=0.0;
        for (int i=0; i<L; ++i){
          XdotP += X(i,l)*P[i];
        }
        
        grad=XY[l]-XdotP;
        
        double beta_new=softmax(beta_old+grad/(B[l]),(lambda[l])/(B[l]));
        
        double beta_diff=beta_new-beta_old;
        if (beta_diff!=0){
          IntegerVector classe=NN[l];
          
          for (int cl=0; cl<classe.size(); ++cl){
            int cla=classe[cl]-1;
            Xbeta[cla]=Xbeta[cla]+X(cla,l)*beta_diff;
            E[cla]=exp(Xbeta[cla]);
          }
          
          
          for (int li=0; li<n; ++li){
            S[li]=epsilon;
            for (int cl=li; cl<L;(cl=cl+n) ){
              S[li] += E[cl];
            }
          }
          
          for (int li=0; li<n; ++li){
            for (int cl=li; cl<L;(cl=cl+n) ){
              S0[cl]=S[li];
            }
          }
          P=E/S0;
          beta[l]=beta_new;
          
          norme2_beta_diff += beta_diff*beta_diff;
        }
      
        if (beta_new!=0 && beta_old==0){
          beta_resamp[l]=1;
        }
        
        else if (beta_new==0 && beta_old==0){
          beta_resamp[l]=(beta_resamp[l]-decay_min)*decay_rate+decay_min;
        }
        
        
        norme2_beta_old += beta_old*beta_old;
      }
    }
    
    incr=sqrt(norme2_beta_diff)/(sqrt(norme2_beta_old)+std::numeric_limits<double>::epsilon());
    if (incr<tolerance){
      break;
    }
  }
  

  return(List::create(
      Rcpp::Named("beta")=beta,
      Rcpp::Named("Xbeta")=Xbeta,
      Rcpp::Named("S")=S)
  );
}
