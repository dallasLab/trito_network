#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

List Erlang_Malaria_Cut(
                double t,  
                NumericVector y, 
                int patch_num,
                NumericVector params,
                NumericMatrix disp_mat) {
        
        //Demographic parameters
        NumericVector b_H = rep(params["b_H"],patch_num); //Human birth rate
        NumericVector b_P = rep(params["b_P"],patch_num); //P.vector birth rate
        NumericVector b_S =  rep(params["b_S"],patch_num);  //S. vector birth rate
        NumericVector mu_H = rep(params["mu_H"],patch_num);  //Human death rate
        NumericVector mu_P =  rep(params["mu_P"],patch_num); //P. vector death rate
        NumericVector mu_S = rep(params["mu_S"],patch_num); // S. vector death rate
        
        //Force of infection parameters
        NumericVector a_P = rep(params["a_P"],patch_num); //biting rate of the p. vector
        NumericVector a_S =  rep(params["a_S"],patch_num);  //biting rate of the s.vector
        
        NumericVector phi_P =   rep(params["phi_P"],patch_num);  //transmission probability of p. vector
        NumericVector phi_S =    rep(params["phi_S"],patch_num);   //transmission probability of s. vector
        NumericVector phi_H =   rep(params["phi_H"],patch_num);   //transmission probability of human
        
        // Recovery rate
        NumericVector gamma <- rep(params["gamma"],patch_num);   //recovery rate of infected human
        
        //#competition coefficient
        NumericVector c_PS <-  rep(params["c_PS"],patch_num);   //competitition effect of p.vector on s.vector
        NumericVector c_SP <-  rep(params["c_SP"],patch_num);  //competition effect of s.vector on p.vector
        
        // FOI
        NumericVector FOI_P <- a_P * phi_P; // #FOI for a primary vector
        NumericVector FOI_S <- a_S * phi_S  //#FOI for a secondary vector
        NumericVector FOI_H_P <-  a_P * phi_H //#FOI for a human to a primary vector
        NumericVector FOI_H_S <- a_S * phi_S //#FOI for a human to a secondary vector
        
        
        
        
        Rcpp::NumericVector dy(patch_num * 7);
       
       //H_S
       for (int i=0; i < patch_num; ++i) {
               dy[i] = (b_H[i] * y[i])-  ;  
       }

        ////////////////////////////
        //Infected Red Blood Cells//
        ////////////////////////////
        
        dy[1] = ((1.0 - c) * 1.0 * p * on_off* y[0] * y[n1+1]) - (n1 * alpha1 * y[1])- (muI * y[1]);
        
        for (int i=1; i < n1; ++i) {
                dy[1+i] = (n1 * alpha1 * y[i]) - (n1 * alpha1 * y[i+1])- (muI  * y[1+i]);  
        }
        
        //////////////
        //Merozoites//
        //////////////
        
        dy[n1+1]  =  B*(n1*alpha1*y[n1]) - (p * y[n1+1]* y[0]) - muM * y[n1+1]; 
        
        ////////////////////////
        //immature gametocytes//
        ////////////////////////
        
        dy[n1+2] =  (c * p * y[0] * y[n1+1]) - muR * y[n1+2] - alpha2*n2*y[n1+2];
        
        for (int i=1; i < n2; ++i) {
                dy[n1+2+i] = alpha2*n2*y[n1+1+i] - n2*alpha2*y[n1+2+i] - muR*y[n1+2+i];
        }
        
        //gametocytes
        dy[n1+n2+2] = (n2 * alpha2 * y[n1+n2+1]) - muG*y[n1+n2+2] ;
        
        
        return List ::create(dy);
        
}