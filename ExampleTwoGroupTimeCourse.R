### Example simulation script
### Two groups time trend

### Set working directory
setwd('working directory')

### Load library snowfall (depends on snow) for parallel computing
library(snowfall)

### Create blocks for snowfall parallel implementation
### 20 blocks of 50 
### each node will do 50 simulations at a time
### 1000 total simulations for a given set of parameters
### the numbers in the blocks matrix serve as the seeds
### in set.seed() for each simulation
blocks <- matrix(1:1000, nrow=20, byrow=TRUE)

### Set values for global variables
groups <- 2
n <- 5
t.max <- 5
m <- 5000
pi1 <- 0.1
FC <- 3
variance <- 0.1
B <- dim(permutations)[2]

### Generate permutations
### Permutations(groups, n, t.max, B=200)
permutations <- Permutations(groups, n, t.max, B)

### Generate group indicator variable for data
I.group.data <- rep(c(rep(0, n), rep(1, n)), t.max)

### INITIALIZE SNOWFALL ###
### Enter number of cpus according to machine
sfInit(parallel=TRUE, cpus=12, type="SOCK")  
sfExportAll()   ### Export all objects in workspace to each node

### Record start time.
start.time <- proc.time() 

### BEGIN ITERATIONS ###
output <- sfApply(blocks, 1, function(k){      # snowfall: each row distributed to different node
  sapply(k, function(l){                       # entries of each row specify the seed for random number generation
    set.seed(l)
    ### Generate gene expression data
    data <- SimulateData2(n, t.max, pi1, FC, variance, m)    
    ### Generate empirical null distribution by permutation
    perms.index <- 1:B
    null.t.stats <- sapply(perms.index, function(i){
      LinReg2(data, permutations[i,])})
    ### Calculate the test statistics for the data
    Tj <- LinReg2(data, I.group.data)  
    ### Calculate raw p-values using pooled null test statistics
    p.values <- sapply(Tj,function(x){(sum(abs(null.t.stats) 
                                           >= abs(x)))/(B*m)}) 
    ### Positive FDR algorithm to obtain q-value for each gene 
    q.values <- Qval(p.values)       
    ### Do significance test and power calculation 
    ### for FDR levels f = 0.05 and f = 0.10
    result <- c(TestPower(q.values, p.values, 0.05, pi1),
                TestPower(q.values, p.values, 0.10, pi1))
    return(result)
  })
})      ### END ITERATIONS ###

sfStop()        ### Stop cluster


### Record end time
end.time <- proc.time()

### Restructure output matrix
output <- matrix(output, nrow=1000, ncol=10, byrow=TRUE)
dimnames(output) <- list(NULL, c("s.05","r.05","power.05",
                                 "estimated.fdr.05","true.fdr.05","s.10","r.10",
                                 "power.10","estimated.fdr.10","true.fdr.10"))

### Save output matrix in csv file
filename <- paste("t.max_",t.max,"_pi1_",pi1,"_FC_",FC,
                  "_n.it_",n,".csv",sep="")
write.csv(output, filename)

### Print out filename, time, means
print(filename)
print(structure(end.time - start.time, class = "proc_time"))
print(c("means",colMeans(output)))