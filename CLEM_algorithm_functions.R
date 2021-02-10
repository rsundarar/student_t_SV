# Functions needed for the CLEM algorithm #


# Load required libraries #
library(Rcpp)
library(RcppNumerical)
library(RcppEigen)
library(BB)
library(lbfgs)


source("CL_algorithm_functions.R") 
# Functions for getting initial guess (Optimizing Eq. (10) of the paper)


sourceCpp("cpp_functions.cpp")
# C++ functions used in the CLEM algorithm

##########################################################################################################


# Finding conditional expectations #

clem_cexp = function( old.pars.value , yt , xt  ){
  
  n = length(yt)
  p = dim(xt)[2]
  
  # old paramater values
  beta.clem = old.pars.value[1:p]
  sig.clem = old.pars.value[p+1]
  nu.clem = old.pars.value[p+2]
  rho.clem = old.pars.value[p+3]
  
  
  # v_t term, n x 1 vector
  vt = 1 + (1-rho.clem)/nu.clem/sig.clem*(yt -  xt%*%matrix(beta.clem,p,1))^2
  
  
  
  # Computations of f(y_t , y_{t-1})######################################################
  
  term1 =  (1-rho.clem)^(nu.clem/2+1)/pi/nu.clem/sig.clem  #scalar
  
  term2 = (gamma(nu.clem/2+1/2))^2/(gamma(nu.clem/2))^2 #scalar
  
  # (n-1) x 1 vector
  wterm = ( ( yt[2:n] - xt[2:n,]%*%matrix(beta.clem,p,1) )^2*(1-rho.clem)/nu.clem/sig.clem + 
              1 )^(-1)*( 
                ( yt[1:(n-1)] - xt[1:(n-1),]%*%matrix(beta.clem,p,1))^2*(1-rho.clem)/nu.clem/sig.clem + 
                  1 )^(-1)
  
  
  # (n-1) x 1 vector
  hgmterm =  Re( hypergeo(  
    
    A = (nu.clem/2 + 1/2)*1 , B = (nu.clem/2+1/2)*1 , C = nu.clem/2 , 
    
    z = rho.clem*wterm
    
  ) # end hyper()
  ) # end Re()
  
  term3 = wterm^(nu.clem/2+1/2)*hgmterm
  
  fyt = term1*term2*term3  # (n-1) x 1 vector
  # f(y_t,y_{t-1}) for t=2,3,...,n
  
 
  #########################################################################################
  
  
  ### \zeta_t terms ########################################################################
  
  # scalar
  z11 = 2*(1-rho.clem)^(nu.clem/2+2)*gamma(nu.clem/2+1/2)*gamma(nu.clem/2+3/2)/pi/nu.clem^2/sig.clem/(gamma(nu.clem/2))^2
  
  # (n-1) x 1 vector
  z12 = Re( hypergeo(  
    A = nu.clem/2 + 1/2 , B = nu.clem/2+3/2 , C = nu.clem/2 , 
    z = rho.clem/vt[2:n]/vt[1:(n-1)]
  ) )
  
  # denominator term, (n-1) x 1 vectors
  dnm1 = (vt[2:n])^(nu.clem/2+3/2)*(vt[1:(n-1)])^(nu.clem/2+1/2)*fyt 
  dnm2 = (vt[2:n])^(nu.clem/2+1/2)*(vt[1:(n-1)])^(nu.clem/2+3/2)*fyt
  
  # (n-1) x 1 vectors
  zeta1 = z11*z12/dnm1 # \zeta_{t,i}, t=2,3,...,n
  zeta2 = z11*z12/dnm2 # \zeta_{t-1,i}, t=1,2,...,n-1
  
  
  # w code
  zeta.wrapper.1 = function(indxt){
    y1 = yt[indxt]
    mu1 = xt[indxt,]%*%matrix(beta.clem)
    y2 = yt[indxt-1]
    mu2 = xt[indxt-1,]%*%matrix(beta.clem)
    
    return(zeta.f.1(y1,y2,mu1,mu2,sig.clem,nu.clem,rho.clem))
  }
  
  zeta.wrapper.2 = function(indxt){
    y1 = yt[indxt]
    mu1 = xt[indxt,]%*%matrix(beta.clem)
    y2 = yt[indxt-1]
    mu2 = xt[indxt-1,]%*%matrix(beta.clem)
    
    return(zeta.f.2(y1,y2,mu1,mu2,sig.clem,nu.clem,rho.clem))
  }
  
  zeta1 = sapply(2:n , zeta.wrapper.1)   
  zeta2 = sapply(2:n , zeta.wrapper.2)   
  
  ### \tau_t terms #########################################################################
  
  tau.term = sapply(1:(n-1) , tau.wrapper , beta.c = beta.clem , sig.c =sig.clem , 
                    nu.c = nu.clem , rho.c = rho.clem , yt2 = yt , xt2 = xt)  
  
  out.lis = list( zeta1 , zeta2 , tau.term)
  names(out.lis) = c("zeta1","zeta2"  , "tau" )
  return( out.lis  )
  
} # end clem_cexp function


########################################################################################################
########################################################################################################
########################################################################################################


# Update function: inputs old parameters, returns the new updated parameters #
# input a value for nu
update_pars = function(old.pars.value , yt , xt   ){
  
  p = dim(xt)[2]
  n = dim(xt)[1]
  
  beta.clem = old.pars.value[1:p]
  sig.clem = old.pars.value[p+1]
  nu.clem = old.pars.value[p+2]
  rho.clem = old.pars.value[p+3]
  
  
  # read conditional expectations
  cexp = conditional_exp( old.pars.value , yt , xt  )
  

  # Update beta #####################################################################
  

  beta.function.solve = function(betas){
    p = length(betas)
    out.f = rep(0,p)
    for(j in 1:p){
      out.f[j] = sum(  cexp$zeta2*(yt[1:(n-1)] - xt[1:(n-1),]%*%matrix(betas))*xt[1:(n-1),j]  ) + 
        sum(  cexp$zeta1*(yt[2:(n)] - xt[2:(n),]%*%matrix(betas))*xt[2:(n),j]  )
    }
    return( (out.f) )
  } # end function betas.function.solve()
  
  p0 = rnorm(p)
  beta.new = BBsolve(p0 , beta.function.solve , quiet = TRUE)$par
  
  
  #####################################################################################
  
  
  # sigma^2 updated #
  sigma.new = 1/2/(n-1)*sum(  cexp$zeta1*( yt[2:(n)] - xt[2:(n),]%*%matrix(beta.new) )^2  )  +
    1/2/(n-1)*sum(  cexp$zeta2*( yt[1:(n-1)] - xt[1:(n-1),]%*%matrix(beta.new) )^2  )
  
  # Finding value of Q^1
  q1.value = -log(sigma.new) - 
    1/2/sigma.new*( sum(cexp$zeta1*(yt[2:(n)] - xt[2:(n),]%*%matrix(beta.new))^2)  ) - 
    1/2/sigma.new*( sum(cexp$zeta2*(yt[1:(n-1)] - xt[1:(n-1),]%*%matrix(beta.new))^2)  )
  
  
  # \rho update function #
  rho.function.optim = function(rp){
    
    
    t11 = -1*( 2*sum(cexp$tau) + (n-1)*nu.clem/2 )*log( 1 - rp )
    
    
    t12 = -nu.clem*( sum(cexp$zeta1 + cexp$zeta2) )/2/(1-rp)  
    
    
    t13 = sum(cexp$tau)*log(rp)
    
    
    return( -1*( t11 +  t12 + t13 ) )
    
  } # end function rho.function.optim()
  
  # \rho update function #
  gradient.rho.function.optim = function(rp){
    
    
    t11 = 1*( 2*sum(cexp$tau) + (n-1)*nu.clem/2 )*1/( 1 - rp )
    
    
    t12 = -nu.clem*( sum(cexp$zeta1 + cexp$zeta2) )/2/(1-rp)/(1-rp)  
    
    
    t13 = sum(cexp$tau)*1/(rp)
    
    
    return( -1*( t11 +  t12 + t13 ) )
    
  } # end function rho.function.optim()
  
  # optimize function for getting \rho
    #nr.curr = optim( 0.5 , rho.function.optim, 
    #             lower = c( 0.05 ) , upper=c( 0.95 ) , method="L-BFGS-B" )
  

  nr.curr = nlm( rho.function.optim ,  0.5) 

    
  nr.new = nr.curr
  
  start.rho.vec = c(0.1,0.3,0.7,0.9)
  for(rho.start in start.rho.vec) {
    
    #nr.curr = optim( rho.start , rho.function.optim, 
    #                 lower = c( 0.05 ) , upper=c( 0.95 ) , method="L-BFGS-B" )
    
    
    # Using nlm() function 
    nr.curr = nlm( rho.function.optim , rho.start )
                     
    
    if(nr.curr$minimum<nr.new$minimum) nr.new = nr.curr
    
  }
  

  
  out.lis = c(beta.new , sigma.new , nu.clem , nr.new$estimate , q1.value , nr.new$minimum  )
  #out.lis = c(beta.new , sigma.new , nu.clem , nr.new$xopt , q1.value , nr.new$fopt  )
  names(out.lis) = c(paste0("beta",1:p),"sigma^2" , "nu" , "rho" , "q1v" , "q2v")
  
  
  return(out.lis) # returns updated set of parameters 
} # end function update_pars()


################################################################################################


# final wrapper function, inputs a value for nu #
clem_final_pars = function(yt,xt,nuval){
  
  n = length(yt)
  p =  dim(xt)[2]
    
  estimates.initial <- initial.guess(yt,xt,nuval) # initial guess for beta and sigma^2
  estimates.initial <- c(estimates.initial , nuval , 0.5 )
  
  pars.old = estimates.initial[1:(p+3)] # Initialize
  pars.new = estimates.initial[1:(p+3)] # Initialize
  
  
  stop.value = 100 # initialize stopping criterion
  
  niter = 1
  
  while(stop.value>0.00000001 & niter<20000){ # iterate until stopping criterion
    
    pars.new = update_pars(pars.old , yt , xt  ) # update function
    
    pars.old = pars.new[1:(p+3)]
    
    # stopping criterion
    if(niter>1) stop.value = max( abs(q1.old - pars.new[p+4]) , abs(q2.old - pars.new[p+5])  )
    
    q1.old = pars.new[p+4] # old value of Q^1
    q2.old = pars.new[p+5] # old value of Q^2
    
    
    #print(niter)
    #print(pars.new[1:(p+3)])
    niter = niter+1
  }
  
  names(pars.new[1:(p+3)]) = c(paste0("beta",1:p),"sigma^2" , "nu" , "rho") 
  return(pars.new[1:(p+3)]) # returns final set of parameters 
  
} # end function clem_final_pars()  
