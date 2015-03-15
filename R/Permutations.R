#' Permutations
#' 
#' Function to generate random permutations. For two groups, permutes labels independently at each time point. For one group time trend, permutes time points. 
#' @param groups: number of groups (1 or 2)
#' @param n sample size per group
#' @param t.max number of time points
#' @param B number of permutations per significance test, maximum 2000. Recommend 200 for two group single time point, 100 for two group time course, 1000 for one group time course.
#' @return matrix of permutations 

Permutations <- function(groups, n, t.max, B=200){
  if(groups==2 && t.max==1){   
    ### two groups single time point
    ### generate 2000 permutations
    permutations <- t(replicate(2000, sample(2*n, n)))
    ### remove duplicates
    permutations <- permutations[!duplicated(permutations),]
    ### take the first B unique permutations
    if(dim(permutations)[1] > B) permutations <- permutations[1:B,]
  } else if(groups==2 && t.max > 1){ 
    ### two groups time course
    ### does B unique permutations of group assignment indicator variable
    ### generate 2000 permutations of indicator variable
    permutations <- matrix(
      replicate(2000, 
                as.vector(replicate(t.max, sample(c(rep(0,n), rep(1,n)))))),
      nrow=2000, ncol=t.max*2*n, byrow=TRUE)
    ### remove duplicates
    permutations <- permutations[!duplicated(permutations),]
    ### take the first B unique permutations
    permutations <- permutations[1:B,]
  } else if(groups==1){ 
    ### one group time course 
    ### generate 2000 permutations
    permutations <- t(replicate(2000, sample(t.max, t.max)))
    ### remove duplicates
    permutations <- permutations[!duplicated(permutations),]
    ### take the first B unique permutations
    if(dim(permutations)[1] > B) permutations <- permutations[1:B,]
  }
  return(permutations)
}