# Functions needed for the CL algorithm #


# Load required libraries #
library(hypergeo)
library(gsl)
library(ordinal)
library(BMS)
library(optimParallel)
library(fAsianOptions)
library(cubature)


#########################################################################################################

# Initial guess: function to output initial guesses of (\beta,\sigma^2,\nu)

# yt: Input vector of the response variable
# xt: A (n \times p) matrix of covariates
# Note that n = series length, p = dimension of covariates

initial.guess = function(yt,xt,nuval){
  
  n = length(yt)
  p = dim(xt)[2]
  
  
  obj.fn.initial.guess = function( param , p.in = p , n.in = n ){
    
    beta.in = param[1:p.in] #\beta
    sigma2.in = param[p.in+1]  #\sigma^2
    
    term1 = lgamma((nuval + 1)/2)  - lgamma(nuval/2) - 0.5*log(sigma2.in) - 0.5*log(nuval)
    
    sum1 = 0
    for(i in 1:n.in) { sum1 = sum1 + log( 1 + (1/nuval/sigma2.in)*(yt[i] -  
                                                                     (matrix(xt[i,],1,p.in)%*%matrix(beta.in,p.in,1)))^2 )
    }
    
    term2 = (nuval + 1)/2*sum1
    
    return(-term1*n.in + term2 ) # returns negative log-likelihood
  } # end internal function obj.fn.initial.guess
  
  
  # constrained optimization 
  optim.par = optim( par = c( rnorm(p) , 1 ) , 
                             fn = obj.fn.initial.guess , 
                             lower = c(rep(-5,p) , 0.01 ) ,
                             upper = c(rep(5,p) , 10  ) ,
                             method = "L-BFGS-B" )$par

  names(optim.par) = c(paste0("beta",1:p),"sigma^2"  ) 
  
  return(optim.par)
  
} # end function initial.guess()


#####################################################################################################


# Bivariate Composite likelihood function
obj.fn.mcl.alt.1 = function( param , p.in , n.in  , yt1 , xt1 , nuinp ){
  
  beta.in = param[1:p.in]  #\beta
  
  sigma2.in = param[p.in+1] #\sigma^2
  
  nu.in = nuinp #\nu, known input
  
  rho.in = param[p.in+2] # \rho
  
  term1 =  (1-rho.in)^(nu.in/2+1)/pi/nu.in/sigma2.in  #scalar
  
  term2 = (gamma(nu.in/2+1/2))^2/(gamma(nu.in/2))^2 #scalar
  
  # (n-1) x 1 vector
  wterm = ( ( yt1[2:n.in] - xt1[2:n.in,]%*%matrix(beta.in,p.in,1) )^2*(1-rho.in)/nu.in/sigma2.in + 
              1 )^(-1)*( 
                ( yt1[1:(n.in-1)] - xt1[1:(n.in-1),]%*%matrix(beta.in,p.in,1))^2*(1-rho.in)/nu.in/sigma2.in + 
                  1 )^(-1)
  
  # (n-1) x 1 vector
  hgmterm =  Re( hypergeo(  
    
    A = nu.in/2 + 1/2 , B = nu.in/2+1/2 , C = nu.in/2 , 
    
    z = rho.in*wterm
    
  ) # end hyper()
  ) # end Re()
  
  term3 = wterm^(nu.in/2+1/2)*hgmterm
  
  fyt = term1*term2*term3  # (n-1) x 1 vector
  # f(y_t,y_{t-1}) for t=2,3,...,n  
  
  return(-sum(log(fyt)))
  
} # end function obj.fn.mcl



# student-t regression: function to output student-t time series regression coefficients
# composite likelihood inference is implemented: MCL(m+1) #

# yt: Input vector of the response variable
# xt: A (n \times p) matrix of covariates
# Note that n = series length, p = dimension of covariates
# m: input for composite likelihood

student.t.mcl = function(yt,xt,nuval){
  
  n = length(yt)
  p = dim(xt)[2]
  
  # Obtain initial guess for the optimization of negative log likelihood of t distribution 
  initial.value = initial.guess(yt,xt,nuval)
 
  #print(initial.value)
  
  # starting value for \rho
  rho.seq = seq(0.05,0.95,0.05)
  #n.rho = length(rho.seq)
  
  cl1 <- makeCluster(detectCores()-1)
  setDefaultCluster(cl=cl1)
  clusterEvalQ(cl1,{library(BMS)
    library(hypergeo)
  })
  
  #optim.value.c = lapply( rho.seq , optimize.over.rho , initial1 = initial.value , cl1 )
  optim.value.c <- optimParallel( par = c( initial.value , rho.seq[1])  , 
                                  fn = obj.fn.mcl.alt.1 , 
                                  lower = c( rep(-5,p) , 0.01 ,  0.05 ) ,
                                  upper = c( rep(5,p) , 10 ,  0.95 ) ,
                                  method = "L-BFGS-B", 
                                  p.in = p , n.in = n ,
                                  yt1=yt,xt1=xt , nuinp = nuval  )
  
  optim.value.begin = optim.value.c$value
  optim.par = optim.value.c$par
  optim.value = optim.value.begin
  
  for(rr in 2:length(rho.seq)){
    optim.value.c <- optimParallel( par = c(initial.value , rho.seq[rr])  , 
                                    fn = obj.fn.mcl.alt.1 , 
                                    lower = c( rep(-5,p) , 0.01 ,  0.05 ) ,
                                    upper = c( rep(5,p) , 10  ,  0.95 ) ,
                                    method = "L-BFGS-B", 
                                    p.in = p , n.in = n ,
                                    yt1=yt,xt1=xt , nuinp = nuval  )
    
    if(optim.value.c$value < optim.value){
      optim.value = optim.value.c$value
      optim.par = optim.value.c$par
    }
    #print(rr) 
  }
  stopCluster(cl1)
  
  optim.out = c( optim.par[1:(p+1)] , nuval , optim.par[p+2] )
  
  names(optim.out) = c(paste0("beta",1:p),"sigma^2" , "nu" , "rho" )
  
  return( optim.out )
  
} # end function student.t.mcl



#########################################################################################################


# Generation of the GAR(1) process by Sim (1990)
# Obs. # rgamma: mgf (1 - t*s)^(-a);  s=scale and a=shape

rgar<-function(n,rho,phi){
  
  n1 = n
  x<-array(0,c(n1,1))
  
  epsilon<-rgamma(n1, shape=phi,  scale = (1-rho)/phi)
  
  x[1]<-rgamma(1, shape=phi, scale = 1/phi)
  
  for(i in 2:n1){
    x[i]<-rgamma(1,shape=rpois(1,lambda=rho*phi*x[i-1]/(1-rho)), scale = (1-rho)/phi ) +
      epsilon[i]
  }
  
  return(x)
}

####################################################################################

initial.guess.student.t = function(yt,xt){
  
  n = length(yt)
  p = dim(xt)[2]
  
  
  obj.fn.initial.guess = function( param , p.in = p , n.in = n ){
    
    beta.in = param[1:p.in] #\beta
    sigma2.in = param[p.in+1]  #\sigma^2
    nuval = param[p.in+2] # \nu
    
    term1 = lgamma((nuval + 1)/2)  - lgamma(nuval/2) - 0.5*log(sigma2.in) - 0.5*log(nuval)
    
    sum1 = 0
    for(i in 1:n.in) { sum1 = sum1 + log( 1 + (1/nuval/sigma2.in)*(yt[i] -  
                                                                     (matrix(xt[i,],1,p.in)%*%matrix(beta.in,p.in,1)))^2 )
    }
    
    term2 = (nuval + 1)/2*sum1
    
    return(-term1*n.in + term2 ) # returns negative log-likelihood
  } # end internal function obj.fn.initial.guess
  
  
  # constrained optimization 
  optim.par = optim( par = c( rnorm(p) , 1 , 1 ) , 
                     fn = obj.fn.initial.guess , 
                     lower = c(rep(-5,p) , 0.01 , 0.01 ) ,
                     upper = c(rep(5,p) , 10 , 10  ) ,
                     method = "L-BFGS-B" )$par
  
  names(optim.par) = c(paste0("beta",1:p),"sigma^2"  , "nu" ) 
  
  return(optim.par)
  
} # end function initial.guess()
