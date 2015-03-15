#' T-statistics
#' 
#' Function to calculate t-statistics to test difference of two group means.
#' @param group.1 group 1 gene expression matrix, rows are genes, columns are subjects
#' @param group.2 group 2 gene expression matrix, rows are genes, columns are subjects
#' @return vector of t-statistics

Tstats <- function(group.1, group.2){
  mean.1 <- rowMeans(group.1)
  mean.2 <- rowMeans(group.2)
  sd.1 <- apply(group.1,1,sd)
  sd.2 <- apply(group.2,1,sd)
  n1 <- dim(group.1)[2]
  n2 <- dim(group.2)[2]
  pooled.variance <- ((n1-1)*sd.1^2+(n2-1)*sd.2^2)/(n1+n2-2)
  pooled.sd <- sqrt(pooled.variance)
  t <- as.matrix((mean.1-mean.2)/(pooled.sd*sqrt(1/n1 + 1/n2)))
  return(t)
}