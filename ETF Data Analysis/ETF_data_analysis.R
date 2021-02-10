# Application of the NSVt regression model #

# Regression of Stock price returns (vs) US Sector ETF returns #


# Load the y and x for the two regressions in Eq. (15) of the paper
# This data corresponds to the training period June 11-14, 2019
load("returns_data_june_11_2019.RData")

n = dim(x)[1] # sample size

p = dim(x)[2] # no. of US sector ETFs chosen

# ETF Tickers for the p=11 US sectors plus Chevron (CVX) and Exxon (XOM) #
etf_tkr <- c("XLF" , "XLE" , "OIH" , "XLK" , "XLP" ,
             "XLV" , "XLU" , "GDX"   , "XLI" , "IYE" , "XME" ,
             "CVX" , "XOM" ) 

####################################################################################

  #  Run the CL and CLEM methods to get parameter estimates # 

  source("CL_algorithm_functions.R") # source functions for the CL method
  
  source("CLEM_algorithm_functions.R") # source functions for the CLEM method
  
  
  # initial guess for \nu for Chevron CVX
  initial.guess.CVX = initial.guess.student.t(y1,x)
  # For Exxon XOM, use y2 above
  
  # \nu parameter for Chevron CVX
  nu.initial.estimate.CVX = tail(initial.guess.CVX,1)
  # For Exxon XOM, use y2 above
  
  
  # CL method estimates for Chevron CVX
  estimates.cl.CVX <- student.t.mcl(y1 , x , nu.initial.estimate.CVX ) 
  names(estimates.cl.CVX)[1:p] = etf_tkr[1:p]
  # For Exxon XOM, use y2
  
  
  # CLEM method estimates for Chevron CVX
  estimates.clem.CVX <- clem_final_pars(y1 , x , nu.initial.estimate.CVX ) 
  names(estimates.clem.CVX)[1:p] = etf_tkr[1:p]
  # For Exxon XOM, use y2 above
  
####################################################################################

