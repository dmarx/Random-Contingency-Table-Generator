mcmc_contingency_tables = function(root_table, k=1e3, thin=1){
  
  flatten = function(mat){matrix(mat, 1, d1*d2)} # convert to row vector
  
  d1 = dim(root_table)[1]; d2 = dim(root_table)[2]
  d1_seq = seq(1,d1); d2_seq = seq(1,d2)
  
  # pre-allocate matrix to store MCMC samples
  hx = matrix(NA_real_, k, d1*d2)
  last_table = root_table
  
  hx_ix=1
  hx[hx_ix,] = flatten(root_table)
  
  U=runif(k)
  for(iter in 2:k){
    i = sample(d1, 1)
    j = sample(d1_seq[-i],1)
    p = sample(d2,1)
    q = sample(d2_seq[-p],1)
    #print(c(i,j,p,q))
    u=U[iter]
        
    v = min(u, last_table[i,p], last_table[j,q])
    
    last_table[i,p] = last_table[i,p] - v
    last_table[j,p] = last_table[j,p] + v
    last_table[i,q] = last_table[i,q] + v
    last_table[j,q] = last_table[j,q] - v
    
    if(iter%%thin==0){
      hx_ix=hx_ix+1
      hx[hx_ix,] = flatten(last_table)      
    }
  }
  hx[1:floor(k/thin),]
}

###############################################################################

set.seed(123)

# Let's' build a simple 2x10 contingency table
N = 1e3
n_values = 10
class = rbinom(N, 1, .7) 
values = sample(n_values,N,replace=TRUE, prob=seq(1,10))
tab = table(class, values)/N

# Thin=30 does a good job of knocking down the autocorrelation 
random_tables = mcmc_contingency_tables(tab, 3e4, 30)

acf(random_tables[,1])
acf(random_tables[,3])
acf(random_tables[,5])

rowSums(matrix(random_tables[1e3,], 2)) # .3 .7, original probabilities maintained after 3e4 iterations
colSums(matrix(random_tables[1e3,], 2)) # ditto
(1:10)/sum(1:10) # Sanity check column marginal probabilities

###############################################################################

# Now let's try this with 5 classes instead of 2

n_values = 10
class = rbinom(N, 4, (1:4)/sum(1:4)) 
values = sample(n_values,N,replace=TRUE, prob=seq(1,10))
tab2 = table(class, values)/N

random_tables2 = mcmc_contingency_tables(tab2, 3e4, 100)

# We could maybe thin it a bit more, but the autocorrelation here seems tolerable to me
acf(random_tables2[,1])
acf(random_tables2[,2])
acf(random_tables2[,3])
acf(random_tables2[,4])
acf(random_tables2[,5])
