#' Simulate Data, Two Groups
#' 
#' Generates data for two groups single time point (t.max = 1) or the two group time trend simple linear regression model.
#' @param n number of subjects per group at each time point
#' @param t.max number of time points
#' @param pi1 proportion differentially expressed genes
#' @param FC fold change at time t.max
#' @param variance gene expression variance at each time point 
#' @param m number of genes
#' @return one m*n*t.max array of gene expression data

SimulateData2 <- function(n, t.max, pi1, FC, variance, m){
  data <- array(,dim=c(m, 2*n, t.max))
  if(t.max==1){
    control.data <- matrix(rnorm(m*n, mean=0, 
                                 sd=sqrt(variance)), nrow=m, ncol=n)
    treatment.data <- matrix(c(
      ### expressed genes
      rnorm(pi1*m*n, mean=log2(FC), sd=sqrt(variance)), 
      ### unexpressed genes
      rnorm((1-pi1)*m*n, mean=0, sd=sqrt(variance))
    ), nrow=m, ncol=n, byrow=TRUE)
    data <- cbind(control.data, treatment.data)
    return(data)    
  } else {
    for (t in 1:t.max){
      control.data <- matrix(rnorm(m*n, mean=0,
                                   sd=sqrt(variance)), nrow=m, ncol=n)
      treatment.data <- matrix(c(
        ### expressed genes
        rnorm(pi1*m*n, mean=log2(FC)*((t-1)/(t.max-1)), sd=sqrt(variance)),
        ### unexpressed genes
        rnorm((1-pi1)*m*n, mean=0, sd=sqrt(variance))
      ), nrow=m, ncol=n, byrow=TRUE)
      data[,,t] <- cbind(control.data, treatment.data)
    }
  }
  return(data)
}