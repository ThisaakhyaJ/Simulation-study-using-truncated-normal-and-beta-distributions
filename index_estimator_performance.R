# This module has all functions needed to explore the performance of the index estimator


# Libraries

library(msm)
library(tidyverse)

##  General functions

# Actual Index

Index=function(Indices,Nitems){
  index=sum(Nitems*Indices)/sum(Nitems)
  index  
}

# Estimated index

Estimated.index=function(Estimated.Indices,Nitems){
  estimated.index=sum(Nitems*Estimated.Indices)/sum(Nitems)
  estimated.index  
}

# Rounds decimal numbers 

My.round=function(Data){
  X=Data
  n=length(X)
  Rx=X
  for(i in 1:n){
    y=X[i]
    if(y>0){
      fy=floor(y)
      dy=y-fy
      if(dy<0.5){ry=fy}
      else {ry=fy+1}
      Rx[i]=ry}
    
    if(y<0){
      fy=ceiling(y)
      dy=fy-y
      if(dy<0.5){ry=fy}
      else {ry=fy-1}
      Rx[i]=ry}
  }
  
  Rx
}

# Mode based on density (Adjusted sample mode)

My.dmode <- function(X) {
  X=as.numeric(unlist(X))
  Freq=c()
  U=1:5
  for(i in 1:5){
    Freq[i]=sum((X==i)*1)
  }
  Freq[is.na(Freq)]=0
  Freq[1]=2*Freq[1]
  Freq[5]=2*Freq[5]
  modes=U[Freq == max(Freq)]
  dmode=mean(modes)
  adjusted.freq=max(Freq)
  list(dmode=dmode, adjusted.freq=adjusted.freq)
}

# Mode based on frequency

My.fmode <- function(X) {
  #This function estimates the mode based on frequencies
  Freq=c()
  U=1:5
  for(i in 1:5){
    Freq[i]=sum((X==i)*1)
  }
  Freq[is.na(Freq)]=0
  modes=U[Freq == max(Freq)]
  fmode=mean(modes)
  freq=max(Freq)
  list(fmode=fmode, freq=freq)
}


# Generating items for a construct

Generate.items=function(c, range=1, nitems){
  if((c==1)){Items=c+rbinom(nitems,2,0.1)}
  if((c==5)){Items=c-rbinom(nitems,2,0.1)}
  if((c>1) & (c<5) ){Items=round(c+rnorm(nitems,0,1/2))
  Items[Items<1]=1
  Items[Items>5]=5
  }
  Items
}

## Index simulations from truncated normal distribution

# Calculates the parameter sigma of a trnorm distribution when sd of data is specified

My.tnormsig=function(mod,std,  a=1, b=5){
  if(std>1){stop("std should be less than 1")}
  sig1=0.1
  sig2=6
  tol=10^(-4)
  
  za=(a-mod)/sig1
  zb=(b-mod)/sig1
  z=pnorm(zb)-pnorm(za)
  f1=(1-(zb*dnorm(zb)-za*dnorm(za))/z-((dnorm(zb)-dnorm(za))/z)^2)*sig1^2-std^2
  
  za=(a-mod)/sig2
  zb=(b-mod)/sig2
  z=pnorm(zb)-pnorm(za)
  f2=(1-(zb*dnorm(zb)-za*dnorm(za))/z-((dnorm(zb)-dnorm(za))/z)^2)*sig2^2-std^2
  if(f1*f2>0){stop("no solution")}
  repeat{
    sig3=(sig1+sig2)/2
    za=(a-mod)/sig3
    zb=(b-mod)/sig3
    z=pnorm(zb)-pnorm(za)
    f3=(1-(zb*dnorm(zb)-za*dnorm(za))/z-  ((dnorm(zb)-dnorm(za))/z)^2)*sig3^2-std^2
    #print(c(sig3, f3))
    if(abs(f3)<tol){break}
    if(f1*f3<0){
      sig2=sig3
      f2=f3}
    if(f3*f2<0){sig1=sig3
    f1=f3
    }}
  sig3
}

# Calculates the theoretical parameters of the truncated normal distribution for each item

My.tnormpara=function(mod, std, a=1, b=5){
  sig=My.tnormsig(mod,std)
  del=ptnorm(mod+0.5, mod, sig, a,b)-ptnorm(mod-0.5, mod, sig, a,b)
  list(std=std,  lam=mod, del=del)
}

# Actual construct value - Trnorm

Index.for.construct.trnorm=function(Items,std){
  nitems=length(Items)
  Lambda=c()
  W=c()
  for(i in 1:nitems){
    mod=Items[i]
    OUT=My.tnormpara(mod,std,a=1, b=5)
    Lambda[i]=OUT$lam
    W[i]=OUT$del
  }
  ic=sum(W*Lambda)/sum(W)
  ic
}

# Construct estimate for - Trnorm

Estimated.index.for.construct.trnorm=function(Items,std, ss){
  nitems=length(Items)
  
  Dmode=c()
  W=c()
  for(i in 1:nitems){
    mod=Items[i]
    sigma=My.tnormsig(mod,std,  a=1, b=5)
    X=rtnorm(ss,mod, sigma)
    R=My.round(X)
    OUT=My.dmode(R)
    Dmode[i]=OUT$dmode
    OUT=My.fmode(R)
    W[i]=OUT$freq
    
  }
  #print(cbind(Dmode, W))
  eic=sum(W*Dmode)/sum(W)
  eic
}

# Index estimate - Trnorm

Trnorm.simulation.index=function(C,Stds=c(0.5,0.5), Nitems=c(8,6), ss=50, nsim=10000){
  
  Items1=Generate.items(C[1],1, Nitems[1])
  Items2=Generate.items(C[2],1, Nitems[2])
  Ic1=Index.for.construct.trnorm(Items1,std=Stds[1])
  Ic2=Index.for.construct.trnorm(Items2,std=Stds[2])
  Indices=c(Ic1, Ic2)
  index=Index(Indices,Nitems) # Actual index
  
  Econstruct1=c()
  Econstruct2=c()
  Eindex=c()
  
  for(sim in 1:nsim){
    Econstruct1[sim]=Estimated.index.for.construct.trnorm(Items1,std=Stds[1], ss)
    Econstruct2[sim]=Estimated.index.for.construct.trnorm(Items2,std=Stds[2], ss)
    Eindex[sim]=Estimated.index(c(Econstruct1[sim],Econstruct2[sim] ), Nitems)
  }
  
  df<-data.frame(Actual=rep(index,nsim), index_estimate=Eindex)
  return(df)
  
}

## Index simulations from beta distribution

# Calculates the parameters alpha and beta for a beta distribution

My.alphabeta=function(mod, std, a=1, b=5){

  if(std>1){stop("std should be 1 or smaller")}
  f=function(alpha,mod, std, a=1, b=5){
    m=(mod-a)/(b-a)
    beta=(alpha-1+2*m-m*alpha)/m
    v=alpha*beta/((alpha+beta+1)*(alpha+beta)^2)
    if(v<0){stop("negative variance")}
    (b-a)*sqrt(v)-std
  }
  g=function(beta,std, a=1, b=5){
    #This function determines the beta to get mode=1
    alpha=1
    v=alpha*beta/((alpha+beta+1)*(alpha+beta)^2)
    if(v<0){stop("negative variance")}
    (b-a)*sqrt(v)-std
  }
  
  h=function(alpha,std, a=1, b=5){
    #This function determines the beta to get mode=5
    beta=1
    v=alpha*beta/((alpha+beta+1)*(alpha+beta)^2)
    if(v<0){stop("negative variance")}
    (b-a)*sqrt(v)-std
  }
  
  if((mod>1) & (mod< 5)){
    alpha1=1.0001 #should be greater than 1
    alpha2=75
    f1=f(alpha1, mod, std)
    f2=f(alpha2, mod, std)
    if(f1*f2>0){stop(paste("Bad initial interval for alpha: alpha1=", alpha1, "alpha2=", alpha2, "mod= ",mod,"std= ",std) )}
    repeat{
      alpha3=(alpha1+alpha2)/2
      f3=f(alpha3, mod, std)
      if(abs(f3)<10^(-4)){break}
      if((f1*f3)<0){alpha2=alpha3
      f2=f3}
      if((f3*f2)<0){alpha1=alpha3
      f1=f3}
    }
    alpha=alpha3
    m=(mod-a)/(b-a)
    beta=(alpha-1+2*m-m*alpha)/m
  }
  
  if(mod==1){alpha=1
  beta1=1.0001
  beta2=75
  
  f1=g(beta1,  std)
  f2=g(beta2,  std)
  if(f1*f2>0){stop(paste("Bad initial interval for beta: beta1=", beta1, "beta2=", beta2, "mod= ",mod,"std= ",std) )}
  repeat{
    beta3=(beta1+beta2)/2
    f3=g(beta3, std)
    if(abs(f3)<10^(-4)){break}
    if((f1*f3)<0){beta2=beta3
    f2=f3}
    if((f3*f2)<0){beta1=beta3
    f1=f3}
  }
  alpha=alpha
  beta=beta3
  
  }

  if(mod==5){
    alpha1=1.001 #should be less than 1
    alpha2=75
    beta=1
    
    f1=h(alpha1,  std)
    f2=h(alpha2,  std)
    if(f1*f2>0){stop(paste("Bad initial interval for alpha: alpha1=", alpha1, "alpha2=", alpha2, "mod= ",mod,"std= ",std) )}
    repeat{
      alpha3=(alpha1+alpha2)/2
      f3=h(alpha3,  std)
      if(abs(f3)<10^(-4)){break}
      if((f1*f3)<0){alpha2=alpha3
      f2=f3}
      if((f3*f2)<0){alpha1=alpha3
      f1=f3}
    }
    alpha=alpha3
    beta=beta
  }
  
  list(alpha=alpha, beta=beta)
}

# Calculates the theoretical parameters of a general beta distribution

My.betapara=function(mod, std, a=1, b=5){
  alphabeta= My.alphabeta(mod, std, a=1, b=5)
  alpha=alphabeta$alpha
  beta=alphabeta$beta
  m=(mod-a)/(b-a)
  del=pbeta(m+0.5, alpha, beta)-pbeta(m-0.5, alpha, beta)
  list(std=std, lam=mod, del=del)
}

# Actual construct value - Beta

Index.for.construct.beta=function(Items,std){
  nitems=length(Items)
  Lambda=c()
  W=c()
  for(i in 1:nitems){
    mod=Items[i]
    OUT=My.betapara(mod,std,a=1, b=5)
    Lambda[i]=OUT$lam
    W[i]=OUT$del
  }
  ic=sum(W*Lambda)/sum(W)
  ic
  
}

# Construct estimate for - Beta

Estimated.index.for.construct.beta=function(Items,std, ss){
  nitems=length(Items)
  
  Dmode=c()
  W=c()
  for(i in 1:nitems){
    mod=Items[i]
    alphabeta=My.alphabeta(mod, std,  a=1, b=5)
    alpha=alphabeta$alpha
    beta=alphabeta$beta
    Y=rbeta(ss, alpha, beta)
    a=1
    b=5
    X=(b-a)*Y+a
    R=My.round(X)
    OUT=My.dmode(R)
    Dmode[i]=OUT$dmode
    OUT=My.fmode(R)
    W[i]=OUT$freq
    
  }
  eic=sum(W*Dmode)/sum(W)
  eic
}

# Index estimate - Beta

Beta.simulation.index=function(C,Stds=c(0.75,0.75), Nitems=c(8,6), ss=50, nsim=10000){
  
  Items1=Generate.items(C[1],1, Nitems[1])
  Items2=Generate.items(C[2],1, Nitems[2])
  
  Ic1=Index.for.construct.beta(Items1,std=Stds[1])
  Ic2=Index.for.construct.beta(Items2,std=Stds[2])
  Indices=c(Ic1, Ic2)
  index=Index(Indices,Nitems) # Actual index
  
  Econstruct1=c()
  Econstruct2=c()
  Eindex=c()
  
  for(sim in 1:nsim){
    Econstruct1[sim]=Estimated.index.for.construct.beta(Items1,std=Stds[1], ss)
    Econstruct2[sim]=Estimated.index.for.construct.beta(Items2,std=Stds[2], ss)
    Eindex[sim]=Estimated.index(c(Econstruct1[sim],Econstruct2[sim] ), Nitems)
  }
  
  df<-data.frame(Actual=rep(index,nsim), index_estimate=Eindex)
  return(df)
  
}


