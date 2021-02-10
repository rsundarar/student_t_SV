#-------------------------------------------------------------
# Joint density function

f_joint<-function(y,x,muy,mux,sigma2,nu,rho){
  omega<-(1+(1-rho)*(y-muy)^2/(nu*sigma2))^(-1)*(1+(1-rho)*(x-mux)^2/(nu*sigma2))^(-1)
  (1-rho)^(nu/2+1)/(pi*nu*sigma2)*(gamma((nu+1)/2)/gamma(nu/2))^2*
    omega^((nu+1)/2)*Re(hypergeo((nu+1)/2,(nu+1)/2,nu/2,rho*omega))
}




#-------------------------------------------------------------
# Conditional Expectations


#-------------------------------------------------------------
# tau=E(U|Yobs)  

tau.f<-function(y,x,muy,mux,sigma2,nu,rho){
  
  vy<-1+(1-rho)*(y-muy)^2/(nu*sigma2)
  vx<-1+(1-rho)*(x-mux)^2/(nu*sigma2)
  
  return( 2*rho*(1-rho)^(nu/2+1)/(pi*nu^2*sigma2)*(gamma((nu+3)/2)/gamma(nu/2))^2*
            Re(hypergeo((nu+3)/2,(nu+3)/2,nu/2+1,rho/(vy*vx)))/
            ((vy*vx)^((nu+3)/2)*f_joint(y,x,muy,mux,sigma2,nu,rho))
  )
}

tau.wrapper = function(indxt,beta.c,sig.c,nu.c,rho.c,yt2,xt2){
  y1t = yt2[indxt+1]
  mu1t = xt2[indxt+1,]%*%matrix(beta.c)
  y2t = yt2[indxt]
  mu2t = xt2[indxt,]%*%matrix(beta.c)
  
  return(tau.f(y1t,y2t,mu1t,mu2t,sig.c,nu.c,rho.c))
}

#-------------------------------------------------------------
# zeta=E(Z|Yobs)  

# for zeta_t, t=2,3,...n
zeta.f.1 <- function(y,x,muy,mux,sigma2,nu,rho){
  vy<-1+(1-rho)*(y-muy)^2/(nu*sigma2)
  vx<-1+(1-rho)*(x-mux)^2/(nu*sigma2)
  2*(1-rho)^(nu/2+2)/(pi*nu^2*sigma2)*gamma((nu+1)/2)*gamma((nu+3)/2)/gamma(nu/2)^2*
    Re(hypergeo((nu+1)/2,(nu+3)/2,nu/2,rho/(vy*vx)))/
    (vy^((nu+3)/2)*vx^((nu+1)/2)*f_joint(y,x,muy,mux,sigma2,nu,rho))
}


# for zeta_t, t=1,2,...n-1
zeta.f.2 <- function(y,x,muy,mux,sigma2,nu,rho){
  vy<-1+(1-rho)*(y-muy)^2/(nu*sigma2)
  vx<-1+(1-rho)*(x-mux)^2/(nu*sigma2)
  2*(1-rho)^(nu/2+2)/(pi*nu^2*sigma2)*gamma((nu+1)/2)*gamma((nu+3)/2)/gamma(nu/2)^2*
    Re(hypergeo((nu+1)/2,(nu+3)/2,nu/2,rho/(vy*vx)))/
    (vx^((nu+3)/2)*vy^((nu+1)/2)*f_joint(y,x,muy,mux,sigma2,nu,rho))
}

