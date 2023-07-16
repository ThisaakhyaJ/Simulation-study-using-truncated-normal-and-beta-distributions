# This module has functions needed for simulations from truncated normal distributions

# Libraries
library(Metrics)

# Mean
Tr.mean<-function(a,b,mu,sigma){
  alpha<-(a-mu)/sigma
  beta<-(b-mu)/sigma
  dalpha<-dnorm(alpha,0,1)
  dbeta<-dnorm(beta,0,1)
  palpha<-pnorm(alpha,0,1)
  pbeta<-pnorm(beta,0,1)
  mubar<-mu-sigma*((dbeta-dalpha)/(pbeta-palpha))
  mubar
}

# Standard deviation
Tr.sd<-function(a,b,mu,sigma){
  alpha<-(a-mu)/sigma
  beta<-(b-mu)/sigma
  dalpha<-dnorm(alpha,0,1)
  dbeta<-dnorm(beta,0,1)
  palpha<-pnorm(alpha,0,1)
  pbeta<-pnorm(beta,0,1)
  v1<-(alpha*dalpha-beta*dbeta)/(pbeta-palpha)
  v2<-((dalpha-dbeta)/(pbeta-palpha))^2
  sigmabar<-(sigma^2)*(1+v1-v2)
  sigmabar
  sqrt(sigmabar)
}

# Median
Tr.median<-function(a,b,mu,sigma){
  alpha<-(a-mu)/sigma
  beta<-(b-mu)/sigma
  palpha<-pnorm(alpha,0,1)
  pbeta<-pnorm(beta,0,1)
  quant<-(palpha+pbeta)/2
  medbar<-mu+(qnorm(quant,0,1))*sigma
  medbar
}

# Mode
Tr.mode<-function(a,b,mu,sigma, na.rm=TRUE){
  if(mu<a){
    a
  }else if(mu>b){
    b
  }else{
    mu
  }
}

# Simulation function
Simulation<-function(Response,n, nsim=10000){
  Parameters<-read.csv("Parameters.csv")
  Samp.mean=c()
  Samp.median=c()
  Samp.mode=c()
  Samp.freq.adj.mode=c()
  tr.mean=Parameters$Mean[Parameters$Response==Response & Parameters$Fixed=="Mean"]
  tr.median=Parameters$Median[Parameters$Response==Response & Parameters$Fixed=="Median"]
  tr.mode=Parameters$Mode[Parameters$Response==Response & Parameters$Fixed=="Mode"]
  
  ## Simulations
  
  # Simulation 1 - Mean
  set.seed(1)
  for(i in 1:nsim){
    X1=rtnorm(n,Parameters$mu[Parameters$Response==Response & Parameters$Fixed=="Mean"],
              Parameters$sigma[Parameters$Response==Response  & Parameters$Fixed=="Mean"], 1, 5)
    D1=round(X1)
    Samp.mean[i]=mean(D1)
  }
  
  # Simulation 2 - Median
  set.seed(1)
  for(j in 1:nsim){
    X2=rtnorm(n,Parameters$mu[Parameters$Response==Response & Parameters$Fixed=="Median"],
              Parameters$sigma[Parameters$Response==Response & Parameters$Fixed=="Median"], 1, 5)
    D2=round(X2)
    Samp.median[j]=median(D2)
  }
  
  # Simulation 3 - Mode
  set.seed(1)
  for(k in 1:nsim){
    X3=rtnorm(n,Parameters$mu[Parameters$Response==Response & Parameters$Fixed=="Mode"],
              Parameters$sigma[Parameters$Response==Response & Parameters$Fixed=="Mode"], 1, 5)
    D3=round(X3)
    Samp.mode[k]=mode(D3)
    Samp.freq.adj.mode[k]=adjmode(D3)
  }
  
  
  # Metrics 
  bias.samp.mean=round(mean(Samp.mean-tr.mean),2)
  bias.samp.median=round(mean(Samp.median-tr.median),2)
  bias.samp.mode=round(mean(Samp.mode-tr.mode),2)
  bias.samp.freq.adj.mode=round(mean(Samp.freq.adj.mode-tr.mode),2)
  se.samp.mean=round( sqrt(mean(Samp.mean-tr.mean)^2),4)
  se.samp.median=round(sqrt(mean(Samp.median-tr.median)^2),4)
  se.samp.mode=round(sqrt(mean(Samp.mode-tr.mode)^2),4)
  se.samp.freq.adj.mode=round(sqrt(mean(Samp.freq.adj.mode-tr.mode)^2),4)
  list(Response=Response, tr.mean=tr.mean,tr.median=tr.median,tr.mode=tr.mode, n=n, bias.samp.mean=bias.samp.mean,
       bias.samp.median=bias.samp.median, bias.samp.mode=bias.samp.mode, bias.samp.freq.adj.mode=bias.samp.freq.adj.mode,
       se.samp.mean=se.samp.mean, se.samp.median=se.samp.median, se.samp.mode=se.samp.mode, se.samp.freq.adj.mode=se.samp.freq.adj.mode)
}



