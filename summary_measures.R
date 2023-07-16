# This module has all functions needed to calculate summary measures mode and adjusted sample mode

# Mode
mode <- function(X) {
  Freq=c()
  U=1:5
  for(i in 1:5){
    Freq[i]=sum((X==i)*1)
  }
  Freq[is.na(Freq)]=0
  modes=U[Freq == max(Freq)]
  mean(modes)
}

# Adjusted sample mode
adjmode <- function(X) {
  Freq=c()
  U=1:5
  for(i in 1:5){
    Freq[i]=sum((X==i)*1)
  }
  Freq[is.na(Freq)]=0
  Freq[1]=2*Freq[1]
  Freq[5]=2*Freq[5]
  modes=U[Freq == max(Freq)]
  mean(modes)
}