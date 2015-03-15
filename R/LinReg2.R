#' Linear Regression, Two Groups Time Trend  
#' 
#' For two group, time trend simple linear regression model of gene expression. Function fits OLS linear regression model to each gene (j = 1, ..., m): 
#' Y.j = beta.0.j + beta.1.j(t) + beta.2.j(I.group) + beta.3.j(t*I.group)
#' Returns t-statistic for test of different slopes. Assumes group sizes are equal, i.e. n1 = n2.
#' @param data a m*n*t.max array of gene expression data, where m is number of genes, n is total number of subjects at each time point (i.e. the sum of the number of subjects in each group, n = n1 + n2, assumes n1 = n2), and t.max is number of time points
#' @param I.group vector of 0/1 variables indicating group membership for each subject at each time point, length n*t.max. For example, I.group = rep(c(rep(0, n1), rep(1, n2)), t.max).
#' @return vector of t-statistics for interaction term coefficient

LinReg2 <- function(data, I.group){
  m <- dim(data)[1]
  n <- dim(data)[2]
  ni <- n/2
  t.max <- dim(data)[3]
  ones <- rep(1, t.max*n)
  t <- rep(1:t.max, each=n)
  t.I.group <- t*I.group
  y.jt <- t(apply(data, 1, as.vector))
  X <- matrix(c(ones, t, I.group, t.I.group), nrow=t.max*n, ncol=4)
  beta.hat <- t(apply(y.jt, 1, function(y){
    solve(t(X) %*% X) %*% t(X) %*% y}))
  y.jt.hat <- t(X %*% t(beta.hat))
  sigma.j.hat.2 <- rowSums((y.jt - y.jt.hat)^2)/(t.max*n - 4)
  ### Get variance of beta.3.j from the variance-covariance matrix
  ### of the least squares parameter estimates
  var.beta.3.j <- solve(t(X) %*% X)[4,4]*sigma.j.hat.2
  Tj <- beta.hat[,4]/sqrt(var.beta.3.j)
  return(Tj)
}