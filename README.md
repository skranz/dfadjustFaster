Test of a modified version of Michal Kolesar's `dfadjust` package

https://github.com/kolesarm/Robust-Small-Sample-Standard-Errors

This version achieves quite substantial speed-ups simply by using Sebastian Krantz's (not me, there is an extra t!) `collapse` package (https://github.com/SebKrantz/collapse) and replacing
the `tapply(x, clustervariable, sum)` and `apply(some_matrix,2,tapply(x,clustervaiable,sum))` calls with `collapse::fsum`.

## Example and Speed comparision

Speed gains seem largest if there are many x variables. That is because the call to `df0` inside the code loops through all `x` variables and `collapse::fsum` substantially speeds up `df0`.

Since it looks as if currently only `lm` is supported and not yet direct fixed effects models like `feols` from `fixest` the number of x variables may indeed be large in applications that use dummy encoding for fixed effects. 

In the following example with many x, the modification yields roughly a 10x speed-up.

```r
N = 1000
S = 20 # num_cluster
K = 200 # num_x

X = matrix(rnorm(N*K),N,K)
colnames(X) = paste0("x",1:K)
clustervar = as.factor(rep(1:S, each = N/S))
y = rnorm(N)
dat = cbind(data.frame(y=y, as.data.frame(X)))
reg = lm(y~., data=dat)

# Original dfadjust
start.time <- Sys.time()
org = dfadjust::dfadjustSE(reg, clustervar)
Sys.time()-start.time
```

```
## Time difference of 8.170841 secs
```

```r
# Faster version using collapse::fsum
start.time <- Sys.time()
fast = dfadjustFaster::dfadjustSE(reg, clustervar)
```

```r
Sys.time()-start.time
```

```
## Time difference of 0.7872219 secs
```

```r
# Roughly a 10x speed improvement in this example

# Results are essentially equal
all.equal(org$coefficients,fast$coefficients)
```

```
## [1] TRUE
```
