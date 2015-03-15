#' Q-value function
#' 
#' Applies positive FDR algorithm to obtain q-value for each gene.
#' @param p.values vector of p-values
#' @param lambda pFDR tuning parameter (between 0 and 1)
#' @return vector of q-values

Qval <- function(p.values, lambda=0.5) {
  ### Estimate number of null genes
  m <- length(p.values)
  m0hat <- sum(p.values > lambda)/(lambda) 
  pi0hat <- m0hat/m
  ### Sort p-values from largest to smallest
  ordered.p.values <- sort(p.values,decreasing=TRUE,index.return=TRUE) 
  p_m <- ordered.p.values[[1]][1]                ### Value of largest p-value
  q.p_m <- ((m0hat*p_m)/sum(p.values <= p_m))    ### Obtain q-value estimate for the largest p-value
  q.values <- matrix(NA, nrow=m, ncol=1)
  q.values[1,] <- q.p_m
  for (i in 2:m) {
    q.values[i,] <- min(q.values[i-1,], ((m0hat*ordered.p.values[[1]][i])/
                                           sum(p.values<=ordered.p.values[[1]][i])))
  }
  return(q.values)
}