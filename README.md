[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="888" alt="Visit QuantNet">](http://quantlet.de/)

## [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **simulations** [<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/)

```yaml

Name of Quantlet: simulations

Published in: submitted to N/A 

Description: ‘Performs the Penalized Adaptive Method (PAM), a combination of propagation-separation approach and a SCAD penalty, to fit a model to a simulated data with two pre-defined change points. Computes the percentage of correctly identifing the homogeneous intervals.’

Keywords: ‘linear model, regression, SCAD penalty, bic, change points, simulations, bootstrap’

Author: Lenka Zboňáková; Xinjue Li

Submitted:  25 Nov. 2018 

Input: 
- n.obs   : Number of observations
- n.par   : Number of parameters
- n.sim   : Number of simulated scenarios
- r       : Correlation parameter of the design matrix X
- sd.eps  : Standard deviation of the error term
- cp1.seq : Sequence of change points
- n.boot  : Number of bootstrap loops
- K       : Sequence of increments between adjacent subintervals
- mb.type : Distribution of multipliers (Bound, Exp or Pois)

```

### R Code
```r

# Clear all variables
rm(list = ls(all = TRUE))
graphics.off()

# Set directory
setwd("C:/Users/xinjueli/PAM/simulations")

# Install and load packages
libraries = c("MASS", "bbmle", "glmnet", "doParallel", "LaplacesDemon", "optimx",
              "lars", "scales", "tilting")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)} )
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

source("PAMfunc.r")

# Simulation setup
n.obs       = 500                   # No of observations
n.par       = 10                    # No  of parameters
n.sim       = 50                    # No of simulations
seed1       = 20171110              # Seed simulation X
seed2       = 20170602              # Seed simulation epsilon
r           = 0.5                   # Correlation parameter of X
sd.eps      = 1                     # Standard deviation of epsilon
cp1.seq     = matrix(c(50,150, 100,150), nrow = 2, ncol = 2, byrow = TRUE)              # Define change points
a           = 3.7                   # Recommended value of parameter a for SCAD
n.boot      = 10                    # Number of bootstrapped multipliers
K           = 50                    # Increment between adjacent sub-intervals
M           = n.obs/K               # Number of sub-intervals
mb.type     = "Exp"               # Bound, Exp or Pois
sel         = 0                     # Position of s

# Initiate cluster for parallel computing
n.cores = detectCores()           # Number of cores to be used
cl      = makeCluster(n.cores)
registerDoParallel(cl)
getDoParWorkers()

# Define number of scenarios in every core
n.simcores = rep((n.sim %/% n.cores), n.cores)
h.simcores = n.sim %% n.cores
if (h.simcores != 0){
  n.simcores[1:h.simcores] = n.simcores[1:h.simcores] + 1
}

# True beta coefficients with change in t = cp1
tmp1.1  = rep(1, 5)
tmp1.2  = rep(0, n.par - length(tmp1.1))
b1      = c(tmp1.1, tmp1.2)

tmp2.1  = rep(1, 3)
tmp2.2  = rep(0, n.par - length(tmp2.1))
b2      = c(tmp2.1, tmp2.2)

tmp3.1  = rep(1, 5)
tmp3.2  = rep(0, n.par - length(tmp1.1))
b3      = c(tmp3.1, tmp3.2)

# Simulation of the design matrix
mu    = rep(0, n.par)
Sigma = matrix(0, nrow = n.par, ncol = n.par)

for (i in 1:n.par) {
  for (j in 1:n.par) {
    if (i == j){
      Sigma[i, j] = 1
    }else {
      Sigma[i, j] = r^abs(i - j)
    }
  }
}

X = list()
for (i in 1:n.sim){
  set.seed(seed1)
  X[[i]] = mvrnorm(n = n.obs, mu, Sigma)
}  

# Simulation of the error term for t = 1, ..., n.obs
eps  = list()
set.seed(seed2)
for (i in 1:n.sim){
  eps[[i]] = rnorm(n.obs, mean = 0, sd = sd.eps)
} 

# Function for computing Y for t = 1, ..., n.obs
Y.sim = function(cp1, cp2){
  # Computation of Y for t = 1, ..., cp1
  Y1    = list()
  for (i in 1:n.sim){
    Y.tmp = numeric(0)
    for (j in 1:cp1){
      Y.tmp = c(Y.tmp, b1 %*% X[[i]][j, ] + eps[[i]][j])
    }
    Y1[[i]] = Y.tmp
  }
  
  # Computation of Y for t = cp1, ..., cp2
  Y2    = list()
  for (i in 1:n.sim){
    Y.tmp = numeric(0)
    for (j in (cp1 + 1):cp2){
      Y.tmp = c(Y.tmp, b2 %*% X[[i]][j, ] + eps[[i]][j])
    }
    Y2[[i]] = Y.tmp
  }
  
  # Computation of Y for t = cp2, ..., n.obs
  Y3    = list()
  for (i in 1:n.sim){
    Y.tmp = numeric(0)
    for (j in (cp2 + 1):n.obs){
      Y.tmp = c(Y.tmp, b3 %*% X[[i]][j, ] + eps[[i]][j])
    }
    Y3[[i]] = Y.tmp
  }
  
  Y   = list()
  for (i in 1:n.sim){
    Y[[i]] = c(Y1[[i]], Y2[[i]], Y3[[i]])
  }
  Y
}

# Simulating mutlipliers from a bounded distribution
runif.mod = function(n){
  unif1    = runif(n, min = 0, max = 12)
  unif.a   = unif1[which(unif1 <= 1)]
  index.a  = which(unif1 <= 1)
  unif.b   = unif1[which((1 < unif1) & (unif1 <= 4))]
  index.b  = which((1 < unif1) & (unif1 <= 4))
  unif.c   = (unif1[which(unif1 > 4)] - 4)/8
  index.c  = which(unif1 > 4)
  unif2    = cbind(c(index.a, index.b, index.c), c(unif.a, unif.b, unif.c))
  unif.fin = unif2[order(unif2[, 1], decreasing = FALSE), 2]
  unif.fin
}



sim.change.points = function(interval.in, s, mb.type = c("Bound", "Exp", "Pois")){
  
    b.pts          = numeric()
    ratio.sele     = numeric()
    Beta.r         = numeric()
    interval.index = numeric()
    k           = interval.in+1
    k1          = interval.in
    k2          = interval.in+1
    beta        = list()
    beta0.tmp   = list()
    a0          = list()
    lambda      = list()
    b.pts.tmp   = 0
    
    # 1.step: Fit the model for assumed homogeneous interval I_t^(1)
    X.tmp.L        = X[[s]][((k1 * K) + 1):(k2 * K), ]
    Y.tmp.L        = Y[[s]][((k1 * K) + 1):(k2 * K)]
    n.obs.tmp.L    = length(Y.tmp.L)
    object.tmp.L   = Onestep.SCAD(X.tmp.L, Y.tmp.L, a, n.obs.tmp.L)
    beta[[k]]      = object.tmp.L$beta
    a0[[k]]        = object.tmp.L$a
    lambda[[k]]    = object.tmp.L$lambda
    beta0.tmp[[k]] = object.tmp.L$beta.0
    
    while (k < M){
      # Fit the model over the next interval I_t^(k + 1) - I_t^(k)
      X.tmp.R      = X[[s]][((k2 * K) + 1):((k2 + 1) * K), ]
      Y.tmp.R      = Y[[s]][((k2 * K) + 1):((k2 + 1) * K)]
      n.obs.tmp.R  = length(Y.tmp.R)
      object.tmp.R = Onestep.SCAD(X.tmp.R, Y.tmp.R, a, n.obs.tmp.R)
      beta.tmp.R   = object.tmp.R$beta
      a.tmp.R      = object.tmp.R$a
      lambda.tmp.R = object.tmp.R$lambda
      beta0.tmp.R  = object.tmp.R$beta.0
      
      # Fit the model over the whole interval I_t^(k + 1)
      X.tmp.T      = X[[s]][((k1 * K) + 1):((k2 + 1) * K), ]
      Y.tmp.T      = Y[[s]][((k1 * K) + 1):((k2 + 1) * K)]
      n.obs.tmp.T  = length(Y.tmp.T)
      object.tmp.T = Onestep.SCAD(X.tmp.T, Y.tmp.T, a, n.obs.tmp.T)
      beta.tmp.T   = object.tmp.T$beta
      a.tmp.T      = object.tmp.T$a
      lambda.tmp.T = object.tmp.T$lambda
      beta0.tmp.T  = object.tmp.T$beta.0
      
      lik.T        = loglik.pen(a.tmp.T, beta0.tmp.T, beta.tmp.T, lambda.tmp.T, 
                                X.tmp.T, Y.tmp.T)
      
      # Simulation of multipliers u_i (i = 1, ..., n.obs.tmp.T)
      set.seed        = 20170424 * k * s * 10
      
      if (mb.type == "Bound"){
        multipliers   = matrix(runif.mod(n.boot * n.obs.tmp.T),      
                               ncol = n.obs.tmp.T, nrow = n.boot) # Bounded distribution
      }
      if (mb.type == "Exp"){
        multipliers   = matrix(rexp(n.boot * n.obs.tmp.T, rate = 1), 
                               ncol = n.obs.tmp.T, nrow = n.boot) # Exp(1) distribution
      }
      if (mb.type == "Pois"){
        multipliers   = matrix(rpois(n.boot * n.obs.tmp.T, lambda = 1), 
                               ncol = n.obs.tmp.T, nrow = n.boot) # Pois(1) distribution
      }
      
      MB.ratios      = numeric(0)
      lik.sel        = numeric(0)
      
      
        lik.ratio    = numeric(0)
        
        # 1.step: Fit the model for assumed homogeneous interval I_t^(k,s)
        X.sel.L      = X[[s]][((k1 * K) + 1):(k2 * K ), ]
        Y.sel.L      = Y[[s]][((k1 * K) + 1):(k2 * K )]
        n.obs.sel.L  = length(Y.sel.L)
        object.sel.L = Onestep.SCAD(X.sel.L, Y.sel.L, a, n.obs.sel.L)
        beta.sel.L   = object.sel.L$beta
        a.sel.L      = object.sel.L$a
        lambda.sel.L = object.sel.L$lambda
        beta0.sel.L  = object.sel.L$beta.0
        
        # 2. a) step: Fit the model over the next interval I_t^(k + 1,s) - I_t^(k,s)
        X.sel.R      = X[[s]][((k2 * K) + 1 ):((k2 + 1) * K), ]
        Y.sel.R      = Y[[s]][((k2 * K) + 1 ):((k2 + 1) * K)]
        n.obs.sel.R  = length(Y.sel.R)
        object.sel.R = Onestep.SCAD(X.sel.R, Y.sel.R, a, n.obs.sel.R)
        beta.sel.R   = object.sel.R$beta
        a.sel.R      = object.sel.R$a
        lambda.sel.R = object.sel.R$lambda
        beta0.sel.R  = object.sel.R$beta.0
        
        # 2. b) step: Fit the model over the whole interval I_t^(k + 1,s)
        X.sel.T      = X.tmp.T
        Y.sel.T      = Y.tmp.T
        n.obs.sel.T  = n.obs.tmp.T
        
        # 3.step: Evaluate the test statistic
        # Real likelihood ratio
        lik.sel.L    = loglik.pen(a.sel.L, beta0.sel.L, beta.sel.L, lambda.sel.L, 
                                  X.sel.L, Y.sel.L)
        lik.sel.R    = loglik.pen(a.sel.R, beta0.sel.R, beta.sel.R, lambda.sel.R, 
                                  X.sel.R, Y.sel.R)
        lik.sel.T    = lik.T
        
        lik.sel      = c(lik.sel, (((n.obs.sel.L/n.obs.sel.T) * lik.sel.L)
                                   + ((n.obs.sel.R/n.obs.sel.T) * lik.sel.R) 
                                   - lik.sel.T))
        
        # 4.step: Find quantiles of LR with use of MB
        for (l in 1:(n.boot)){
          # Multiplier bootstrap for left-hand interval I_t^(k)
          multipl.L      = multipliers[l, 1:n.obs.sel.L]
          
          boot.object.L  = Onestep.SCAD.MB(X.sel.L, Y.sel.L, a, n.obs.sel.L, 
                                           as.numeric(lambda.sel.L), multipl.L)
          
          boot.beta.L    = boot.object.L$beta
          boot.a.L       = boot.object.L$a
          boot.a0.L      = boot.object.L$a.0
          boot.beta0.L   = boot.object.L$beta.0
          
          lik.boot.L     = loglik.pen.MB(boot.a.L, boot.beta0.L, boot.beta.L, 
                                         as.numeric(lambda.sel.L), X.sel.L, Y.sel.L, 
                                         multipl.L)
          
          # Multiplier bootstrap for right-hand interval I_t^(k + 1) - I_t^(k)
          multipl.R      = multipliers[l, (n.obs.sel.L + 1):n.obs.sel.T]
          
          boot.object.R  = Onestep.SCAD.MB(X.sel.R, Y.sel.R, a, n.obs.sel.R, 
                                           as.numeric(lambda.sel.R), multipl.R)
          
          boot.beta.R    = boot.object.R$beta
          boot.a.R       = boot.object.R$a
          boot.a0.R      = boot.object.R$a.0
          boot.beta0.R   = boot.object.R$beta.0
          
          lik.boot.R     = loglik.pen.MB(boot.a.R, boot.beta0.R, boot.beta.R, 
                                         as.numeric(lambda.sel.R), X.sel.R, Y.sel.R, 
                                         multipl.R)
          
          # Multiplier bootstrap for the whole interval I_t^(k + 1) (with shift)
          multipl        = multipliers[l, ]
          
          boot.object.T  = Onestep.SCAD.MB.shift(X.sel.L, Y.sel.L, X.sel.R, Y.sel.R, a, 
                                                 lambda.sel.L, lambda.sel.R, multipl.L, 
                                                 multipl.R, beta.sel.L, beta.sel.R, 
                                                 a.sel.L, a.sel.R)
          
          boot.beta.T    = boot.object.T$beta
          boot.a.T       = boot.object.T$a
          boot.a0.T      = boot.object.T$a.0
          boot.lambda.T  = boot.object.T$lambda
          boot.beta0.T   = boot.object.T$beta.0
          
          Y.boot.T       = c(Y.sel.L, (Y.sel.R - rep((a.sel.R - a.sel.L), n.obs.sel.R) 
                                       - X.sel.R %*% (beta.sel.R - beta.sel.L)))
          
          lik.boot.T     = loglik.pen.MB(boot.a.T, boot.beta0.T, boot.beta.T, 
                                         boot.lambda.T, X.sel.T, Y.boot.T, multipl)
          
          lik.ratio      = c(lik.ratio, ((n.obs.sel.L/n.obs.sel.T) * lik.boot.L 
                                         + (n.obs.sel.R/n.obs.sel.T) * lik.boot.R 
                                         - lik.boot.T))
        }
        
        MB.ratios = cbind(MB.ratios, lik.ratio)
      
      
      # Find maximum over the real and bootstrapped likelihood ratios
      t.stat        = max(lik.sel)
      MB.ratios.max = apply(MB.ratios, 1, max)
      
      q90 = quantile(MB.ratios.max, probs = 1-(0.1*k/M))
      q95 = quantile(MB.ratios.max, probs = 1-(0.05*k/M))
      
      k2  = k2 + 1
      k   = k + 1
      
      # 5.step: Test for homogeneity and evaluate current beta estimator
      if (t.stat <= q95){
        for (k3 in (k1 + 1):k){
          beta[[k3]]   = beta.tmp.T
          a0[[k3]]     = a.tmp.T
          lambda[[k3]] = lambda.tmp.T
        }
      } else {
        b.pts.tmp      = (k - 1) * K + 1
        k1             = k2 - 1
        beta[[k]]      = beta.tmp.R
        a0[[k]]        = a.tmp.R
        lambda[[k]]    = lambda.tmp.R
        index          = which(beta[[(k)]]<= 0.2*max(beta[[(k)]]))
        BETA           = beta[[(k)]]
        BETA[index]    = 0
        if(k<=2){
          ratio        = 1-(length(which(BETA[6:10]!= 0))/length(b1))}
        if(k>2 && k<=3){
          ratio        = 1-(length(which(BETA[4:10]!= 0))/length(b1))}
        else{
          ratio          = 1-(length(which(BETA[6:10]!= 0))/length(b1))}
        break
      } 
    }
    
    if(b.pts.tmp == 0){
     sim.no                    = s
     b.pts[sim.no]           = 500
     index                     = which(beta.tmp.T <= 0.4*max(beta.tmp.T))
     Beta.r[sim.no]          = beta.tmp.T
     Beta.r[sim.no][index]   = 0
     ratio.sele[sim.no]      = 1-(length(which((Beta.r[[sim.no]][6:10])!= 0))/length(b1))
     interval.index[sim.no]  = k
     value                     = list(b.pts[sim.no], ratio.sele[sim.no], Beta.r[sim.no], interval.index[sim.no])
     return(value)
    }
    else{
    Sys.time()
    sim.no                   = s
    b.pts[sim.no]          = b.pts.tmp 
    ratio.sele[sim.no]     = ratio
    Beta.r[sim.no]         = BETA
    interval.index[sim.no] = k-1
    value                    = list(b.pts[sim.no], ratio.sele[sim.no], Beta.r[sim.no], interval.index[sim.no])
     return(value)
    }
}

results = function(h){
value.interval.1 = list()
value.ratio.1    = list()
a1               = matrix(0, 50, 2)
a2               = matrix(0, 50, 2)
s1 = ifelse(h == 1, 1, sum(n.simcores[1:(h - 1)],1))
s2 = sum(n.simcores[1:h])
for(i in s1:s2){
value.interval        = list()
value.ratio           = list()
values1               = sim.change.points(0, i, mb.type)
j                     = values1[[4]]
value.interval        = j
value.ratio           = values1[[2]]
while(j < 9){
values1               = sim.change.points(j, i, mb.type)
value.interval        = c(value.interval, values1[[4]])
value.ratio           = c(value.ratio, values1[[2]])
j                     = values1[[4]]}
value.interval.1[[i]] = value.interval
if(c(2, 3)%in%value.interval == c(T,T)){a1[i, ]=c(1,1)}
else{a1[i, ]=c(0,0)}
value.ratio.1[[i]]    = value.ratio
if(c(2, 3)%in%value.interval == c(T,T)){
a2[i, ] = value.ratio[c(1,2)]}
else{a2[i, ]=c(0,0)}
}
Value.final           = list(a1[s1:s2, ], a2[s1:s2, ])
names(Value.final)    = list("interval.index", "selection.v.ratio")
return(Value.final)
}

out.sum    = list()
final.list = list() 
Sys.time()
Y   = Y.sim(100, 150)
  final.chpts       = foreach(ipar = 1:n.cores, 
                              .packages = c("MASS", "bbmle", "glmnet", "doParallel", 
                                            "LaplacesDemon", "optimx", "lars", "scales", 
                                            "tilting"), 
                              .combine = 'c') %dopar% results(ipar)   
  final.list = final.chpts
Sys.time()
stopCluster(cl)

subinterval     = final.list[[1]]
for(k in 2:n.cores){
  subinterval   = rbind(subinterval, final.list[[(2*k-1)]])
}

Ratio.select    = final.list[[2]]
for(k in 2:n.cores){
  Ratio.select  = rbind(Ratio.select, final.list[[(2*k)]])
}

corr.ratio.structural.break.1 = sum(subinterval[, 1])/n.sim
corr.ratio.structural.break.2 = sum(subinterval[, 2])/n.sim

corr.ratio.v.selection.1      = mean(Ratio.select[, 1])
corr.ratio.v.selection.2      = mean(Ratio.select[, 2])

corr.ratio.structural.break.1
corr.ratio.structural.break.2
corr.ratio.v.selection.1 
corr.ratio.v.selection.2


```

automatically created on 2019-11-26