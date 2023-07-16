# This module has functions needed for parameter estimation for fitting truncated normal and beta distributions

# Libraries
library(R.utils)
library(msm)
library(fitdistrplus)
library(rootSolve)
library(Metrics)

## Support Functions

# Generating exact values from responses using a uniform proposal model
simulation_uniform_prop<-function(R){
  n<-length(R)
  captureOutput(for (i in 1:n){
    if(R[i]==1){
      print(runif(1,R[i],R[i]+0.5))
    }else if(R[i]==5){
      print(runif(1,R[i]-0.5,R[i]))
    }else{
      print(runif(1,R[i]-0.5, R[i]+0.5))
    }
    
  }) %>% gsub("\\[1] ","",.) %>% as.numeric()
}

## Parameter estimation functions

# Parameter estimation for Trnorm fit 1
sigma.given.tr.var<-function(x, a=1, b=5, Interval = c(0.1,5)){
  rounded_x<-round(x)
  mu<-adjmode(rounded_x)
  sample.var<-var(x)
  f<-function(z){sample.var-Tr.var(a,b,mu,z)}
  Out=uniroot(f, interval=Interval, extendInt = "yes")
  sigma=Out$root
  sigma
}

#  Parameter estimation for Trnorm fit 2
sigma.given.tr.mean<-function(x, a=1, b=5, Interval = c(0.1,5)){
  rounded_x<-round(x)
  mu<-adjmode(rounded_x)
  sample.mean<-mean(x)
  f<-function(z){sample.mean-Tr.mean(a,b,mu,z)}
  Out=uniroot(f, interval=Interval, extendInt = "yes")
  sigma=Out$root
  sigma
}

# Parameter estimation for Trnorm fit 3
sigma.by.min.ssd2<-function(x, a=1, b=5, Interval = c(0.1,5)){
  
  mu<-adjmode(round(x))
  
  # Finding sigma based on Newton Raphson method
  sigma1<-sigma.given.tr.var(x)
  
  # Getting a sequence of sigmas
  sigma_seq<-seq(sigma1-1, sigma1+1, length.out=1000)
  sigma_seq <-sigma_seq[sigma_seq > 0] 
  
  # Empirical densities 
  emp.dens<-ecdf(x)
  emp_dens_x<-emp.dens(x)
  
  # Calculating theoretical densities for each sigma value
  
  ssd<-c()
  for(i in sigma_seq){
    theorotical_dens_x<-ptnorm(x, mean=mu, sd=i, lower=1, upper=5)
    diff_squared<- (theorotical_dens_x-emp_dens_x)^2
    sum_of_squares<-sum(diff_squared)
    ssd<-c(ssd, sum_of_squares)
  }
  
  # Dataframe of sigmas and SSDs
  df<-data.frame(sigma=sigma_seq, ssd=ssd)
  
  # Minimum ssd
  min_ssd<-min(ssd)
  
  
  # Selecting the best sigma
  optimal_sigma<-df$sigma[df$ssd==min_ssd]
  
  optimal_sigma
  
}

# Parameter estimation for Beta fit 1
transformed_beta_parameters<-function(x){
  y<-(x-1)/4
  fit <- fitdist(y, "beta")
  shape1<-fit$estimate[1]
  shape2<-fit$estimate[2]
  l<-list(shape1, shape2)
  return(l)
}

# Parameter estimation for Beta fit 2
transformed_beta_parameters_mode.beta<-function(x){
  mode<-adjmode(round(x))
  y<-(x-1)/4
  fit <- fitdist(y, "beta")
  shape1<-fit$estimate[1]
  shape2<-mode_given_beta(alpha=shape1, mode=mode)
  l<-list(shape1, shape2)
  return(l)
}


# Parameter estimation for Beta fit 3
transformed_beta_parameters_mean.beta<-function(x){
  mean<-mean(round(x))
  y<-(x-1)/4
  fit <- fitdist(y, "beta")
  shape1<-fit$estimate[1]
  shape2<-mean_given_beta(alpha=shape1, mean=mean)
  l<-list(shape1, shape2)
  return(l)
}

## Plot functions

# Plot  for Trnorm fit 1
Plot_truncated_fit1<-function(x, Method,Fit){
  mu<-adjmode(round(x))
  sigma<-sigma.given.tr.mean(x)
  hist(x,prob = TRUE, breaks=c(1,1.5,2.5,3.5,4.5,5), main = paste("Trnorm fit ", Fit), xlab = "Data")
  curve(dtnorm(x, mean=mu, sd=sigma, lower=1, upper=5), col = "blue", add = TRUE)
  text(x = 2, y = 0.3, paste("mu = ", mu ))
  text(x = 2, y = 0.2, paste("sigma = ", round(sigma,2) ))
}

# Plot  for Trnorm fit 1
Plot_truncated_fit2<-function(x, Method,Fit){
  mu<-adjmode(round(x))
  sigma<-sigma.given.tr.var(x)
  hist(x,prob = TRUE, breaks=c(1,1.5,2.5,3.5,4.5,5), main = paste("Trnorm fit ", Fit), xlab = "Data")
  curve(dtnorm(x, mean=mu, sd=sigma, lower=1, upper=5), col = "blue", add = TRUE)
  text(x = 2, y = 0.3, paste("mu = ", mu ))
  text(x = 2, y = 0.2, paste("sigma = ", round(sigma,2) ))
}

# Plot for Trnorm fit 3
Plot_truncated_fit3<-function(x, Fit){
  mu<-adjmode(round(x))
  sigma<-sigma.by.min.ssd2(x)
  hist(x,prob = TRUE, breaks=c(1,1.5,2.5,3.5,4.5,5), main = paste("Trnorm fit ", Fit), xlab = "Data")
  curve(dtnorm(x, mean=mu, sd=sigma, lower=1, upper=5), col = "blue", add = TRUE)
  text(x = 2, y = 0.3, paste("mu = ", mu ))
  text(x = 2, y = 0.2, paste("sigma = ", round(sigma,2) ))
}


# Plot for Beta fit 1
Plot_beta_fit1<-function(x, Fit){
  hist(x,prob = TRUE, breaks=c(1,1.5,2.5,3.5,4.5,5), main = paste("Beta fit ", Fit), xlab = "Data")
  p<-transformed_beta_parameters(x)
  shape1<-p[1] %>% unlist()
  shape2<-p[2] %>% unlist()
  curve(dbeta((x-1)/4 , shape1, shape2)/4, col = "blue", add = TRUE)
  text(x = 2, y = 0.3, paste("alpha = ", round(shape1,2) ))
  text(x = 2, y = 0.2, paste("beta = ", round(shape2,2) ))
}


# Plot for Beta fit 1
Plot_beta_fit2<-function(x,Fit){
  hist(x,prob = TRUE, breaks=c(1,1.5,2.5,3.5,4.5,5), main = paste("Beta fit ", Fit), xlab = "Data")
  p<-transformed_beta_parameters_mode.beta(x)
  shape1<-p[1] %>% unlist()
  shape2<-p[2] %>% unlist()
  curve(dbeta((x-1)/4 , shape1, shape2)/4, col = "blue", add = TRUE)
  text(x = 2, y = 0.3, paste("alpha = ", round(shape1,2) ))
  text(x = 2, y = 0.2, paste("beta = ", round(shape2,2) ))
}


# Plot for Beta fit 1
Plot_beta_fit3<-function(x,Fit){
  hist(x,prob = TRUE, breaks=c(1,1.5,2.5,3.5,4.5,5), main = paste("Beta fit ", Fit), xlab = "Data")
  p<-transformed_beta_parameters_mean.beta(x)
  shape1<-p[1] %>% unlist()
  shape2<-p[2] %>% unlist()
  curve(dbeta((x-1)/4 , shape1, shape2)/4, col = "blue", add = TRUE)
  text(x = 2, y = 0.3, paste("alpha = ", round(shape1,2) ))
  text(x = 2, y = 0.2, paste("beta = ", round(shape2,2) ))
}


## Calculating SSD for each fit

# SSD for Trnorm fit 1
truncated_fit1<-function(x){
  # Empirical densities 
  emp.dens<-ecdf(x)
  emp_dens_x<-emp.dens(x)
  
  # mu and sigma
  mu<-adjmode(round(x))
  sigma<-sigma.given.tr.var(x)
  
  # Theoretical densities
  theorotical_dens_x<-ptnorm(x, mean=mu, sd=sigma, lower=1, upper=5)
  
  # Calculating SSD
  diff_squared<- (theorotical_dens_x-emp_dens_x)^2
  sum_of_squares<-sum(diff_squared)
  
  return(sum_of_squares)
}

# SSD for Trnorm fit 2
truncated_fit2<-function(x){
  # Empirical densities 
  emp.dens<-ecdf(x)
  emp_dens_x<-emp.dens(x)
  
  # mu and sigma
  mu<-adjmode(round(x))
  sigma<-sigma.given.tr.mean(x)
  
  # Theoretical densities
  theorotical_dens_x<-ptnorm(x, mean=mu, sd=sigma, lower=1, upper=5)
  
  # Calculating SSD
  diff_squared<- (theorotical_dens_x-emp_dens_x)^2
  sum_of_squares<-sum(diff_squared)
  
  return(sum_of_squares)
}

# SSD for Trnorm fit 3
truncated_fit3<-function(x){
  # Empirical densities 
  emp.dens<-ecdf(x)
  emp_dens_x<-emp.dens(x)
  
  # mu and sigma
  mu<-adjmode(round(x))
  sigma<-sigma.by.min.ssd2(x)
  
  # Theoretical densities
  theorotical_dens_x<-ptnorm(x, mean=mu, sd=sigma, lower=1, upper=5)
  
  # Calculating SSD
  diff_squared<- (theorotical_dens_x-emp_dens_x)^2
  sum_of_squares<-sum(diff_squared)
  
  return(sum_of_squares)
}

# SSD for Beta fit 1
beta_fit1<-function(x){
  # Empirical densities 
  emp.dens<-ecdf(x)
  emp_dens_x<-emp.dens(x)
  
  # alpha and beta
  p<-transformed_beta_parameters(x)
  shape1<-p[1] %>% unlist()
  shape2<-p[2] %>% unlist()
  
  # Theoretical densities
  theorotical_dens_x<-pbeta((x-1)/4 , shape1, shape2)
  
  # Calculating SSD
  diff_squared<- (theorotical_dens_x-emp_dens_x)^2
  sum_of_squares<-sum(diff_squared)
  
  return(sum_of_squares)
}

# SSD for Beta fit 2
beta_fit2<-function(x){
  # Empirical densities 
  emp.dens<-ecdf(x)
  emp_dens_x<-emp.dens(x)
  
  # alpha and beta
  p<-transformed_beta_parameters_mode.beta(x)
  shape1<-p[1] %>% unlist()
  shape2<-p[2] %>% unlist()
  
  # Theoretical densities
  theorotical_dens_x<-pbeta((x-1)/4 , shape1, shape2)
  
  # Calculating SSD
  diff_squared<- (theorotical_dens_x-emp_dens_x)^2
  sum_of_squares<-sum(diff_squared)
  
  return(sum_of_squares)
}

# SSD for Beta fit 3
beta_fit3<-function(x){
  # Empirical densities 
  emp.dens<-ecdf(x)
  emp_dens_x<-emp.dens(x)
  
  # alpha and beta
  p<-transformed_beta_parameters_mean.beta(x)
  shape1<-p[1] %>% unlist()
  shape2<-p[2] %>% unlist()
  
  # Theoretical densities
  theorotical_dens_x<-pbeta((x-1)/4 , shape1, shape2)
  
  # Calculating SSD
  diff_squared<- (theorotical_dens_x-emp_dens_x)^2
  sum_of_squares<-sum(diff_squared)
  
  return(sum_of_squares)
}


