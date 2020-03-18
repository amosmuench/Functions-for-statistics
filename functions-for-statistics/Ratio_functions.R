#calculates the SEM of the ratio of 2 arrays
#source: http://www.stat.cmu.edu/~hseltman/files/ratio.pdf
ratio_sem<-function(A,B) {
  ratio_standarderrorofmean<-abs(mean(A)/mean(B))* sqrt( (sd(A) /sqrt(length(A))/mean(A))^2 + (sd(B) /sqrt(length(B))/mean(B))^2)
  return(ratio_standarderrorofmean)
}


#calculates the variance of the ratio of 2 arrays
#source: http://www.stat.cmu.edu/~hseltman/files/ratio.pdf
ratio_variance<-function(A,B) {
  ratio_var<-sqrt(mean(A))/sqrt(mean(B))*(var(B)/sqrt(mean(B))+var(B)/sqrt(mean(B)))
  return(ratio_var)
}