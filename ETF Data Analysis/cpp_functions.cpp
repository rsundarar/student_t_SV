// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>
#include <math.h>

#define _USE_MATH_DEFINES

#include <cmath>


using namespace Numer;
using namespace Rcpp;
using namespace std;

double hypergeometric( double a, double b, double c, double x )
{
  const double TOLERANCE = 1.0e-10;
  double term = a * b * x / c;
  double value = 1.0 + term;
  int n = 1;
  
  while ( abs( term ) > TOLERANCE )
  {
    a++, b++, c++, n++;
    term *= a * b * x / c / n;
    value += term;
  }
  
  return value;
}


// Obtain environment containing function
Rcpp::Environment package_env("package:hypergeo"); 

// Make function callable from C++
Rcpp::Function hypg = package_env["hypergeo"];  


//............................................................................//



// [[Rcpp::export]]
Rcpp::List conditional_exp(NumericVector parms, NumericVector yt1, NumericMatrix xt1)
{
  
  int n1 = xt1.nrow(), p1 = xt1.ncol();
  
  NumericVector beta_clem (p1);
  for(int i=0; i<p1; i++) beta_clem[i] = parms[i];

  double sig1 = parms[p1];
  double nu1 = parms[p1+1];
  double rho1 = parms[p1+2];
  
  
  
  // For getting v_t, t=1,2,...,n
  NumericVector vt1 (n1);
  for(int i=0 ; i<n1 ; i++){
    
  NumericVector curr_xt = xt1(i,_);  
  double sq_term = pow((yt1[i] - inner_product(curr_xt.begin(), curr_xt.end() , beta_clem.begin(), 0.0 ) ),2);
    
    vt1[i] = 1 + (1-rho1)/nu1/sig1*sq_term;
  }
  // .................................................................................
  
  
  
  //For getting w_t, t=2,3,...,n-1
  double term1 =  pow( 1.0*(1-rho1) , 1.0*(nu1/2 + 1.0) )/M_PI/nu1/sig1;
    
  double term2 = pow( tgamma( 1.0*(nu1/2+0.5) ) , 2 )/pow(tgamma(1.0*nu1/2) , 2);
  
  NumericVector wterm (n1-1);
  for(int i=1; i<n1; i++){
    NumericVector curr_xt = xt1(i,_); 
    double inst1 = 1+pow((yt1[i] - inner_product(curr_xt.begin(), curr_xt.end() , beta_clem.begin(), 0.0 ) ),2)*(
      (1-rho1)/nu1/sig1);
    
    curr_xt = xt1(i-1,_); 
    double inst2 = 1+pow((yt1[i-1] - inner_product(curr_xt.begin(), curr_xt.end() , beta_clem.begin(), 0.0 ) ),2)*(
      (1-rho1)/nu1/sig1);
      
    wterm[i-1] = pow(inst1,-1)*pow(inst2,-1);  
  }
  //....................................................................................
  

  
  // For getting f_y_t, t=2,3,...,n
  
  NumericVector hgmt (n1-1);
  for(int i=1; i<n1; i++){
   // NumericVector curr_hgm = hypg( 1*(nu1/2 + 0.5) , 1*(nu1/2+0.5) , nu1/2 , rho1*wterm[i]);
    hgmt[i-1] = hypergeometric( 1.0*(nu1/2 + 0.5) , 1.0*(nu1/2+0.5) , 1.0*nu1/2 , rho1*wterm[i-1]);//curr_hgm[0];
 // cout<<"a,b="<<1.0*(nu1/2 + 0.5)<<",c="<<1.0*nu1/2<<",z="<<rho1*wterm[i]<<"--"<<"hgm="<<hgmt[i]<<"--end--";
  }
  
  NumericVector fyt (n1-1);
  for(int i=1; i<n1; i++){
    fyt[i-1] = term1*term2*pow(wterm[i-1] , 1.0*(nu1/2+0.5) )*hgmt[i-1];
  }
  
  
  //..........................................................................................
  
  // zeta_t, t=2,3,...,n
  NumericVector zeta1 (n1-1);
  
  double cons1 = 2.0*pow( (1-rho1)*1.0 , (nu1/2+2.0)*1.0 )*tgamma( (nu1/2+0.5)*1.0 )*tgamma( (nu1/2+1.5)*1.0 );
  double cons2 = 1.0/M_PI/nu1/nu1/sig1/pow( tgamma(nu1/2*1.0)   , 2  );
  
  for(int i=1; i<n1; i++){
    
  zeta1[i-1] = hypergeometric( 1.0*(nu1/2 + 0.5) , 1.0*(nu1/2+1.5) , 1.0*nu1/2 , rho1/vt1[i]/vt1[i-1]);  
  
  double dnm1 = pow( vt1[i] , (nu1/2+1.5)*1.0  )*pow( vt1[i-1] , (nu1/2+0.5)*1.0  )*fyt[i-1];
  
  zeta1[i-1] = cons1*cons2*zeta1[i-1]/dnm1;
    
  }
  
  // zeta_{t-1}, t=2,3,...,n
  NumericVector zeta2 (n1-1);
  
  for(int i=1; i<n1; i++){
    
    zeta2[i-1] = hypergeometric( 1.0*(nu1/2 + 0.5) , 1.0*(nu1/2+1.5) , 1.0*nu1/2 , rho1/vt1[i]/vt1[i-1]);  
    
    double dnm2 = pow( vt1[i-1] , (nu1/2+1.5)*1.0  )*pow( vt1[i] , (nu1/2+0.5)*1.0  )*fyt[i-1];
    
    zeta2[i-1] = cons1*cons2*zeta2[i-1]/dnm2;
    
  }
  
  //.....................................................................................
  
  
  // for tau_t, t=2,3,...,n
  
  NumericVector tau (n1-1);
  
  for(int i=1; i<n1; i++){
    
    double curr_y1 = yt1[i];
    NumericVector curr_xt1 = xt1(i,_); 
    double mu_curr1 = inner_product(curr_xt1.begin(), curr_xt1.end() , beta_clem.begin(), 0.0 );
    
    double v1 = 1+(1-rho1)*pow( (curr_y1 - mu_curr1)*1.0 , 2.0)/nu1/sig1;
    
    
    double curr_y2 = yt1[i-1];
    NumericVector curr_xt2 = xt1(i-1,_); 
    double mu_curr2 = inner_product(curr_xt2.begin(), curr_xt2.end() , beta_clem.begin(), 0.0 );
    
    double v2 = 1+(1-rho1)*pow( (curr_y2 - mu_curr2)*1.0 , 2.0)/nu1/sig1;
    
    
    tau[i-1] = 2.0*rho1*pow( (1-rho1)*1.0 , nu1/2+1.0 )/M_PI/nu1/nu1/sig1*(
      tgamma((nu1/2 + 1.5) )*tgamma((nu1/2 + 1.5) )/tgamma(nu1/2)/tgamma(nu1/2))*(
    hypergeometric( 1.0*(nu1/2 + 1.5) , 1.0*(nu1/2+1.5) , 1.0*(nu1/2+1.0) , rho1/v1/v2))/(
      pow( 1.0*v1*v2 , 1.0*(nu1/2+1.5))*fyt[i-1]);  
      
  }
  
  
  return Rcpp::List::create(
    Rcpp::Named("wt") = wterm,
    Rcpp::Named("zeta1")=zeta1,
    Rcpp::Named("zeta2")=zeta2,
    Rcpp::Named("tau") = tau
  );
  
}
