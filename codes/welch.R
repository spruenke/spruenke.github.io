# Welch-Test
welch.sim = function(x.1, x.2, nboot = 1000, alpha = 0.05){ #x.1 and x.2 are data, nboot the number of simulations and alpha is significance level
  # length of groups
  n.1 = length(x.1)
  n.2 = length(x.2)
  
  x = c(x.1, x.2)
  n = length(x)
  ### permutation
  
  x.perm = matrix(replicate(nboot, sample(x)), byrow = T, nrow = nboot)
  x.1.perm = x.perm[,1:n.1]
  x.2.perm = x.perm[,(n.1+1):length(x)]
  x.1.bar = rowMeans(x.1.perm)
  x.2.bar = rowMeans(x.2.perm)
  x.1.var = (rowSums(x.1.perm^2)-n.1*x.1.bar^2)/(n.1-1)
  x.2.var = (rowSums(x.2.perm^2) - n.2*x.2.bar^2)/(n.2-1)
  T.star.perm = (x.1.bar - x.2.bar)/(sqrt((x.1.var/n.1) + (x.2.var/n.2)))
  
  ### non-parametric
  
  x.npb = matrix(sample(x, n*nboot, replace = T), byrow = T, ncol = n)
  x.1.npb = x.npb[,1:n.1]
  x.2.npb = x.npb[,(n.1+1):length(x)]
  x.1.bar = rowMeans(x.1.npb)
  x.2.bar = rowMeans(x.2.npb)
  x.1.var = (rowSums(x.1.npb^2)-n.1*x.1.bar^2)/(n.1-1)
  x.2.var = (rowSums(x.2.npb^2) - n.2*x.2.bar^2)/(n.2-1)
  T.star.npb = (x.1.bar - x.2.bar)/(sqrt((x.1.var/n.1) + (x.2.var/n.2)))
  
  ### Parametric Bootstrap
  
  x.pb = matrix(rnorm(nboot * length(x), 0, sd(x)), ncol = n, byrow = T)
  x.1.pb = x.pb[,1:n.1]
  x.2.pb = x.pb[,(n.1+1):length(x)]
  x.1.bar = rowMeans(x.1.pb)
  x.2.bar = rowMeans(x.2.pb)
  x.1.var = (rowSums(x.1.pb^2)-n.1*x.1.bar^2)/(n.1-1)
  x.2.var = (rowSums(x.2.pb^2) - n.2*x.2.bar^2)/(n.2-1)
  T.star.pb = (x.1.bar - x.2.bar)/(sqrt((x.1.var/n.1) + (x.2.var/n.2)))
  
  ### Wild Bootstrap
  
  z = c(x.1 - mean(x.1), x.2 - mean(x.2))
  w = matrix(sample(c(1,-1), size = n*nboot, replace=T), ncol = nboot)
  x.star = w*z
  x.star.sq = x.star^2
  
  x.1.wb = x.star[1:n.1,]
  x.2.wb = x.star[(n.1+1):n,]
  
  x.star.sq.1 = x.star.sq[1:n.1,]
  x.star.sq.2 = x.star.sq[(n.1+1):n,]
  
  x.1.bar = colMeans(x.1.wb)
  x.2.bar = colMeans(x.2.wb)
  
  x.1.var = (colSums(x.star.sq.1) - n.1*(x.1.bar)^2)/(n.1-1)
  x.2.var = (colSums(x.star.sq.2) - n.2*(x.2.bar)^2)/(n.2-1)
  
  T.star.wb = (x.1.bar - x.2.bar)/(sqrt((x.1.var/n.1) + (x.2.var/n.2)))
  
  #compute critical values from empirical/bootstrapped distribution
  c.perm  = quantile(T.star.perm, c(alpha/2, 1-alpha/2))
  c.npb   = quantile(T.star.npb, c(alpha/2, 1-alpha/2))
  c.pb    = quantile(T.star.pb, c(alpha/2, 1-alpha/2))
  c.wb    = quantile(T.star.wb, c(alpha/2, 1-alpha/2))
  
  #T-Value of original data
  T.true = (mean(x.1) - mean(x.2))/(sqrt((var(x.1)/n.1) + (var(x.2)/n.2)))
  
  #Check whether H_0 is rejected or not (evaluates to TRUE if H_0 is rejected!)
  perm  = T.true < c.perm[1] | T.true > c.perm[2]
  npb   = T.true < c.npb[1] | T.true > c.npb[2]
  pb    = T.true < c.pb[1] | T.true > c.pb[2]
  wb    = T.true < c.wild[1] | T.true > c.wild[2]
  
  #Results
  res           = data.frame(c(c.perm[1], c.npb[1], c.pb[1], c.wb[1]), c(c.perm[2], c.npb[2], c.pb[2], c.wb[2]), c(perm, npb, pb, wb))
  colnames(res) = c(paste0(alpha/2, "-Quantil"), paste0("1-", alpha/2, "-Quantil"), "Reject H_0:" ) #set column names
  rownames(res) = c("Permutation Bootstrap", "Non-Parametric Bootstrap", "Parametric Bootstrap", "Wild Bootstrap")
  abc           = t.test(x.1, x.2)
  orig          = abc$p.value < 0.05 # If p-value lower than 0.05 than returns TRUE, thus H_0 is rejected
  
  return((list("Resampling" = res, "Classic Two-Sample t-test" = list("P-Value" = abc$p.value, "H_0 rejected?" = orig))))
  
}
