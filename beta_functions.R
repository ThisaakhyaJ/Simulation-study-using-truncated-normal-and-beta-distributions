# This module has functions for simulations from the beta distribution

# Mean
Beta.mean<-function(alpha,beta){
  Beta.mean<-(5*alpha + beta)/(alpha + beta)
  Beta.mean
}

# Variance
Beta.var<-function(alpha,beta){
  Beta.variance<-(16*alpha*beta)/((alpha+beta)^2+(alpha + beta + 1))
  Beta.variance
}

# Median
Beta.med<-function(alpha,beta){
  Beta.median<-(15*alpha+3*beta-6)/(3*alpha+3*beta-2)
  Beta.median
}

# Mode
Beta.mod<-function(alpha,beta){
  Beta.mode<-(5*alpha+beta-6)/(alpha+beta-2)
  Beta.mode
}

# Finding parameter beta when alpha and mean is fixed
mean_set_beta<-function(alpha, mean){
  beta<-(alpha*(5-mean))/(mean-1)
  beta
}

# Finding parameter beta when alpha and median is fixed
median_set_beta<-function(alpha, median){
  beta<-(3*alpha*(5-median)+2*(median-3))/(3*(median-1))
  beta
}

# Finding parameter beta when alpha and mode is fixed
mode_set_beta<-function(alpha, mode){
  beta<-(alpha*(5-mode)+2*(mode-3))/(mode-1)
  beta
}

# Finding parameter alpha when beta and median is fixed 
median_set_alpha<-function(beta, median){
  alpha<-(3*beta*(1-median)+2*(median-3))/(3*(median-5))
  alpha
}

# Finding parameter alpha when beta and mode is fixed
mode_set_alpha<-function(beta, mode){
  alpha<-(beta*(1-mode)+2*(mode-3))/(mode-5)
  alpha
}

# Simulation function
Simulation<-function(Response,n, nsim=1000){
  Parameters<-read.csv("Parameters.csv")
  Samp.mean=c()
  Samp.median=c()
  Samp.mode=c()
  Samp.freq.adj.mode=c()
  beta.mean=Parameters$Mean[Parameters$Response==Response & Parameters$Fixed=="Mean"]
  beta.median=Parameters$Median[Parameters$Response==Response & Parameters$Fixed=="Median"]
  beta.mode=Parameters$Mode[Parameters$Response==Response & Parameters$Fixed=="Mode"]
  
  ## Simulations
  
  # Simulation 1 - Mean
  
  set.seed(1)
  for(i in 1:nsim){
    X1=rbeta(n,Parameters$Alpha[Parameters$Response==Response & Parameters$Fixed=="Mean"],
             Parameters$Beta[Parameters$Response==Response & Parameters$Fixed=="Mean"])
    X1<-4*X1+1
    D1=round(X1)
    Samp.mean[i]=mean(D1)
  }
  
  # Simulation 2 - Median
  
  set.seed(1)
  for(j in 1:nsim){
    X2=rbeta(n,Parameters$Alpha[Parameters$Response==Response & Parameters$Fixed=="Median"],
             Parameters$Beta[Parameters$Response==Response & Parameters$Fixed=="Median"])
    X2<-4*X2+1
    D2=round(X2)
    Samp.median[j]=median(D2)
  }
  
  # Simulation 3 - Mode
  
  set.seed(1)
  for(k in 1:nsim){
    X3=rbeta(n,Parameters$Alpha[Parameters$Response==Response & Parameters$Fixed=="Mode"],
             Parameters$Beta[Parameters$Response==Response & Parameters$Fixed=="Mode"])
    X3<-4*X3+1
    D3=round(X3)
    Samp.mode[k]=mode(D3)
    Samp.freq.adj.mode[k]=adjmode(D3)
  }
  
  
  
  # Metrics 
  bias.samp.mean=round(mean(Samp.mean-beta.mean),2)
  bias.samp.median=round(mean(Samp.median-beta.median),2)
  bias.samp.mode=round(mean(Samp.mode-beta.mode),2)
  bias.samp.freq.adj.mode=round(mean(Samp.freq.adj.mode-beta.mode),2)
  se.samp.mean=round( sqrt(mean(Samp.mean-beta.mean)^2),4)
  se.samp.median=round(sqrt(mean(Samp.median-beta.median)^2),4)
  se.samp.mode=round(sqrt(mean(Samp.mode-beta.mode)^2),4)
  se.samp.freq.adj.mode=round(sqrt(mean(Samp.freq.adj.mode-beta.mode)^2),4)
  list(Response=Response, beta.mean=beta.mean,beta.median=beta.median,beta.mode=beta.mode, n=n, bias.samp.mean=bias.samp.mean,
       bias.samp.median=bias.samp.median, bias.samp.mode=bias.samp.mode, bias.samp.freq.adj.mode=bias.samp.freq.adj.mode,
       se.samp.mean=se.samp.mean, se.samp.median=se.samp.median, se.samp.mode=se.samp.mode, se.samp.freq.adj.mode=se.samp.freq.adj.mode)
}
