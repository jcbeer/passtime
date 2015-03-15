#' Linear Regression, One Group Time Trend  
#' 
#' For one group, time trend simple linear regression model of gene expression. Function fits OLS linear regression model to each gene. Returns t-statistic for test of zero slope.
#' @param data a m*n*t.max array of gene expression data, where m is number of genes, n is number of subjects at each time point, and t.max is number of time points
#' @return vector of t-statistics for slope coefficient

LinReg1 <- function(data){
  m <- dim(data)[1]
  n <- dim(data)[2]
  t.max <- dim(data)[3]
  t <- matrix(rep(1:t.max, each=n), nrow=m, ncol=(n*t.max), byrow=TRUE)
  y.jt <- t(apply(data, 1, as.vector))
  sd.t <- apply(t, 1, sd)
  sd.y.jt <- apply(y.jt, 1, sd)
  r <- sapply(seq.int(m), function(i) cor(t[i,], y.jt[i,]))
  beta.1.hat <- r*(sd.y.jt/sd.t)
  beta.0.hat <- sapply(seq.int(m), function(i) mean(y.jt[i,]) - beta.1.hat[i]*mean(t[i,]))
  fitted.values <- beta.0.hat + t*beta.1.hat 
  residuals.sqrd <- (y.jt - fitted.values)^2
  s.beta.1.hat <- sapply(seq.int(m), function(i){ 
    sqrt(sum(residuals.sqrd[i,])/((n*t.max-2)*(n*t.max-1)*sd.t[i]^2))})
  Tj <- beta.1.hat/s.beta.1.hat
  return(Tj)
}