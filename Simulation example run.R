# Robust Student-t stochastic volatility model #

# Simulation Example #


source("CL_algorithm_functions.R")

source("CLEM_algorithm_functions.R")


p = 5 # Dimension of the covariates

nu.inp = 5 # \nu parameter

rho.inp = 0.3 # \rho parameter

sigma2 = 0.5 # \sigma^2 parameter

phi.inp =  nu.inp/2 # \phi parameter from GAR(1) model

# Set the seed
set.seed(1012)  

# Generate the regression coefficients \beta
beta.coef = runif(p,-1,1)


# Generate ARMA(1,1) coefficients for regressors x_t
phi1 = runif(1,0.3,0.7)
theta1 = runif(1,0.3,0.7)

# Sample size
n = 500
    
    # Initialize matrix of covariates 
    xt = matrix(0,n,p)
    
    # Initialize the mean vector 
    mut = matrix(0,n,1)
    
    # Generate the covariate matrix x_t 
    for(j in 1:n){
      xt[j,] =   arima.sim( model = list(ar = phi1 , ma = theta1  ) , n = p  )
      
      mut[j] = sum( xt[j,]%*%matrix(beta.coef) )
    }
    
    
    # Generate n observations from GAR(1) model by Sim (1990)
    # These are the z_t in Definition 2 of the paper
    eps = rgar(n , rho.inp , phi.inp  )
    
    
    # Generating the response process (Definition 2 of the paper)
    yt = rep(0,n) # Initialize vector of size n
    for(j in 1:n) yt[j]  = rnorm( 1 , xt[j,]%*%matrix(beta.coef) , sqrt(sigma2/eps[j]) )
    
    
    # Initial Guess of all the parameters
    # Optimizing Eq. (10) of the paper
    initial.estimates <- initial.guess.student.t( yt , xt )
    
    # \nu parameter estimate
    nu.initial.estimate = initial.estimates[p+2]
    
    # CL method estimates  
    estimates.cl <- student.t.mcl(yt , xt , nu.initial.estimate ) 
    

    # CLEM method estimates
    estimates.clem <- clem_final_pars(yt,xt , nu.initial.estimate ) 
    

  