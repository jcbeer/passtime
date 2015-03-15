#' Simulate Data, One Group Time Trend  
#' 
#' Generates data for the one group time trend simple linear regression model.
#' @param n number of subjects at each time point
#' @param t.max number of time points
#' @param pi1 proportion differentially expressed genes
#' @param FC fold change at time t.max
#' @param variance gene expression variance at each time point 
#' @param m number of genes
#' @return one m*n*t.max array of gene expression data

SimulateData1 <- function(n, t.max, pi1, FC, variance, m){
  data <- array(,dim=c(m, n, t.max))
  for (t in 1:t.max){
    data[,,t] <- matrix(c(rnorm(pi1*m*n, mean=log2(FC)*((t-1)/(t.max-1)), sd=sqrt(variance)), 
                          rnorm((1-pi1)*m*n, mean=0, sd=sqrt(variance))), nrow=m, ncol=n, byrow=TRUE)
  }
  return(data)
}