# t-test (one sample)
t.sim = function(x, mu.t = 0, nboot = 4000, alpha = 0.05, direction = "equal"){ #test whether x.bar = mu.t, x.bar > mu.t (greater) or x.bar < mu.t (lower)
  n = length(x)
  mu = mean(x)
  
  # non-parametric bootstrap
  x.npb = matrix(sample(x, n*nboot, replace = T), byrow = T, ncol = n)
  x.bar = rowMeans(x.npb)
  x.var = (rowSums(x.npb^2)-n*x.bar^2)/(n-1)
  T.star.npb = (x.bar -mu)/(sqrt(x.var/n))
  
  # parametric bootstrap
  x.pb = matrix(rnorm(nboot * n, 0, sd(x)), ncol = n, byrow = T)
  x.bar = rowMeans(x.pb)
  x.var = (rowSums(x.pb^2)-n*x.bar^2)/(n-1)
  T.star.pb = (x.bar - mu)/(sqrt(x.var/n))
  
  # wild bootstrap
  z = x - mean(x)
  w = matrix(sample(c(1,-1), size = n*nboot, replace=T), ncol = nboot)
  x.star = w*z
  x.star.sq = x.star^2
  x.bar = colMeans(x.star)
  x.var = (colSums(x.star.sq) - n*x.bar^2)/(n-1)
  T.star.wb = x.bar * sqrt(n/x.var)
  
  
  c.wb    = quantile(T.star.wb, c(alpha/2, 1-alpha/2))
  c.npb   = quantile(T.star.npb, c(alpha/2, 1-alpha/2))
  c.pb    = quantile(T.star.pb, c(alpha/2, 1-alpha/2))
  
  c.npb.2 = quantile(T.star.npb, c(alpha, 1-alpha))
  c.pb.2  = quantile(T.star.pb, c(alpha, 1-alpha))
  c.wb.2  = quantile(T.star.wb, c(alpha, 1-alpha))
  
  #T-Value of original data
  T.true = (mu - mu.t)/(sqrt((var(x)/n)))
  
  #Check whether H_0 is rejected or not (evaluates to TRUE if H_0 is rejected!)
  switch(direction,
         equal = {
           npb   = T.true < c.npb[1] | T.true > c.npb[2]
           pb    = T.true < c.pb[1] | T.true > c.pb[2]
           wb    = T.true < c.wb[1] | T.true > c.wb[2]
           res   = data.frame(c(c.npb[1], c.pb[1], c.wb[1]), c(c.npb[2], c.pb[2], c.wb[2]), c(npb, pb, wb))
           abc   = t.test(x, mu = mu.t)},
         greater = {
           npb   = T.true < c.npb.2[1]
           pb    = T.true < c.pb.2[1]
           wb    = T.true < c.wb.2[1]
           res   = data.frame(c(c.npb.2[1], c.pb.2[1], c.wb.2[1]), c(c.npb.2[2], c.pb.2[2], c.wb.2[2]), c(npb, pb, wb))
           abc   = t.test(x, mu = mu.t, alternative = "greater")
         },
         lower    = {
           npb   = T.true > c.npb.2[2]
           pb    = T.true > c.pb.2[2]
           wb    = T.true > c.wb.2[2]
           res   = data.frame(c(c.npb.2[1], c.pb.2[1], c.wb[1]), c(c.npb.2[2], c.pb.2[2], c.wb.2[2]), c(npb, pb, wb))
           abc   = t.test(x, mu = mu.t, alternative = "lower")
         }
  )
  
  
  colnames(res) = c(paste0("Lower Quantile"), paste0("Upper Quantile"), "Reject H_0:" ) #set column names
  rownames(res) = c("Non-Parametric Bootstrap", "Parametric Bootstrap", "Wild Bootstrap")
  abc           = t.test(x, mu = mu.t)
  orig          = abc$p.value < 0.05 # If p-value lower than 0.05 than returns TRUE, thus H_0 is rejected
  
  return((list("Resampling" = res, "Classic t-test" = list("P-Value" = abc$p.value, "H_0 rejected?" = orig))))
}