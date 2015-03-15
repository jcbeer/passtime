#' Significance test and power calculation
#' 
#' For simulated microarray data with proportion pi1 of differentially expressed (DE) genes, this function takes a vector of q-values and a vector of the corresponding p-values, sorted such that all differentially expressed genes come first, and returns s (number of true positives), r (number of genes called significant), power (proportion of DE genes that are called significant), estimated false discovery rate (FDR) (maximum q-value <= desired FDR), and true FDR ((r - s)/r if r not equal to zero, otherwise zero).
#' @param q.values vector of q-values sorted such that all DE genes come first
#' @param p.values vector of corresponding p-values in the same order as the q.values
#' @param f desired false discovery rate (between 0 and 1)
#' @param pi1 the proportion of genes that are truly differentially expressed
#' @return vector of s, r, power, estimated FDR, and true FDR

TestPower <- function(q.values, p.values, f, pi1){
  m <- length(p.values)
  q.val.index <- which(q.values <= f)
  if (length(q.val.index)==0){
    s <- 0
    r <- 0
    power <- 0
    estimated.FDR <- NA
    true.FDR <- NA
  } else {
    ordered.p.values <- sort(p.values,decreasing=TRUE,index.return=TRUE) 
    gene.index <- ordered.p.values[[2]][q.val.index]
    s <- sum(gene.index <= pi1*m)
    r <- length(gene.index)
    power <- s/(pi1*m)
    estimated.FDR <- max(q.values[q.val.index])
    true.FDR <- (r - s)/r    
  }
  out <- c(s,r,power,estimated.FDR,true.FDR)
  return(out) 
}