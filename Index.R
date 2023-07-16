# This module has functions needed to calculate the index for a given lecturer

# Calculating the value of a construct for a given course
non_binary_construct<-function(proportion, adj_sample_mode){
  num<-sum(proportion*adj_sample_mode)
  den<-sum(proportion)
  construct<-num/den
  construct
}

# Calculating the index value of a given course
Index<-function(m, construct){
  num<-sum(m*construct)
  den<-sum(m)
  index<-num/den
  index
}

