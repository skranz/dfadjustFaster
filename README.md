Test of a modified version of Michal Kolesar's `dfadjust` package

https://github.com/kolesarm/Robust-Small-Sample-Standard-Errors

This version achieves speed ups using two modifications:

1. Use Sebastian Krantz's (not me, there is an extra t!) `collapse` package (https://github.com/SebKrantz/collapse) and replace the `tapply(x, clustervariable, sum)` and `apply(some_matrix,2,tapply(x,clustervaiable,sum))` calls with `collapse::fsum`.

2. Add new argument `inv_mode`. Only used for cluster robust se. If "K" then behaviour is as in `dfadjust`:  for each cluster invert (i.e. perform eigenvalue decomposition) of a K x K matrix where K is the number of columns in X.  If "n_s" then for each cluster s invert a n_s x n_s matrix, where n_s is the number of rows in cluster s. If "min" check for each cluster  whether K or n_s is smaller and invert the smaller matrix. For backward compatibility the default is "K", but using "min" optimizes for speed.

## 1. Speed comparision change to `collapse::fsum`

Speed gains seem largest if there are many x variables. That is because the call to `df0` inside the code loops through all `x` variables and `collapse::fsum` substantially speeds up `df0`.

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

## 2. Speed comparison different `inv_mode`


```r
# Example with fe
# comparing different inv_mode options

N = 10000
S = 500 # num_cluster
K = 200 # num_x
num_fe = 200
fe_vals = rnorm(num_fe)
fe_int = sample(1:num_fe,N, replace=TRUE)
fe_val = fe_vals[fe_int]

library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
dat <- tibble(
  x1=rnorm(N,0,1)+fe_val,
  x2=rnorm(N,0,1)+x1,
  cl = as.factor(sample(1:S, N, replace=TRUE, prob = (1:S)^20)),
  #cl = as.factor(rep(1:S, each = N/S)),
  fe = as.factor(fe_int),
  y = 0*x1 + 0*x2 + fe_val + rnorm(N, 0,1)
)
dat = arrange(dat, cl)
sort(table(dat$cl))
```

```
## 
## 335 346 352 360 364 366 369 371 375 376 377 380 393 402 345 356 359 374 383 384 
##   1   1   1   1   1   1   1   1   1   1   1   1   1   1   2   2   2   2   2   2 
## 390 391 400 357 372 382 385 386 387 395 396 398 399 401 378 392 394 397 407 408 
##   2   2   2   3   3   3   3   3   3   3   3   3   3   3   4   4   4   4   4   4 
## 388 403 405 409 414 406 418 404 411 416 423 424 412 413 417 419 422 429 420 421 
##   5   6   6   7   7   8   8   9   9   9   9  10  11  12  12  12  12  12  13  13 
## 425 427 410 434 415 432 433 426 428 431 430 435 438 440 437 445 442 447 436 444 
##  13  13  14  14  15  15  15  16  17  18  22  25  30  30  31  32  33  33  37  37 
## 441 446 439 450 443 452 448 449 453 459 451 457 456 454 455 458 460 462 461 465 
##  38  38  42  43  44  48  51  55  55  62  65  65  66  68  68  70  80  88  89  91 
## 463 471 464 466 468 467 470 469 472 475 474 476 478 477 473 479 480 481 486 482 
##  92  95 100 104 111 114 114 115 126 128 131 144 153 154 159 174 174 185 201 209 
## 483 484 485 489 487 488 490 491 493 494 496 492 495 498 500 497 499 
## 211 215 219 252 262 278 278 298 309 326 326 332 364 379 388 404 411
```

```r
library(dfadjustFaster)

# Spec 1: just add fe which is not nested in cl
reg = lm(y~x1+x2+fe, data=dat)
dfadjustSE(reg,for_coefs = 2:3, dat$cl,inv_mode="n_s", timing=TRUE)
```

```
## Time difference of 5.550497 secs
```

```
## 
## Coefficients:
##       Estimate      HC1 se      HC2 se     Adj. se       df    p-value
## x1 -0.02579377 0.010743285 0.010751923 0.011050633 44.74868 0.02066906
## x2  0.01063598 0.008973349 0.008994776 0.009251056 43.66271 0.24342177
```

```r
dfadjustSE(reg,for_coefs = 2:3, dat$cl,inv_mode="K", timing=TRUE)
```

```
## Time difference of 4.986795 secs
```

```
## 
## Coefficients:
##       Estimate      HC1 se      HC2 se     Adj. se       df    p-value
## x1 -0.02579377 0.010743285 0.010751923 0.011050633 44.74868 0.02066906
## x2  0.01063598 0.008973349 0.008994776 0.009251056 43.66271 0.24342177
```

```r
dfadjustSE(reg,for_coefs = 2:3, dat$cl,inv_mode="min",timing=TRUE)
```

```
## Time difference of 2.796929 secs
```

```
## 
## Coefficients:
##       Estimate      HC1 se      HC2 se     Adj. se       df    p-value
## x1 -0.02579377 0.010743285 0.010751923 0.011050633 44.74868 0.02066906
## x2  0.01063598 0.008973349 0.008994776 0.009251056 43.66271 0.24342177
```

```r
# Spec 2: cluster variable cl is also a fixed effect
reg = lm(y~x1+x2+fe+cl, data=dat)
dfadjustSE(reg,for_coefs = 2:3, dat$cl,inv_mode="n_s", timing=TRUE)
```

```
## Time difference of 6.984944 secs
```

```
## 
## Coefficients:
##        Estimate      HC1 se      HC2 se     Adj. se       df    p-value
## x1 -0.024680185 0.010855348 0.010787876 0.011093674 43.88105 0.02702734
## x2  0.009218437 0.009015457 0.008974077 0.009234686 42.86081 0.31007496
```

```r
dfadjustSE(reg,for_coefs = 2:3, dat$cl,inv_mode="K", timing=TRUE)
```

```
## Time difference of 17.71424 secs
```

```
## 
## Coefficients:
##        Estimate      HC1 se      HC2 se     Adj. se       df    p-value
## x1 -0.024680185 0.010855348 0.010787876 0.011093674 43.88105 0.02702734
## x2  0.009218437 0.009015457 0.008974077 0.009234686 42.86081 0.31007496
```

```r
dfadjustSE(reg,for_coefs = 2:3, dat$cl,inv_mode="min",timing=TRUE)
```

```
## Time difference of 6.393841 secs
```

```
## 
## Coefficients:
##        Estimate      HC1 se      HC2 se     Adj. se       df    p-value
## x1 -0.024680185 0.010855348 0.010787876 0.011093674 43.88105 0.02702734
## x2  0.009218437 0.009015457 0.008974077 0.009234686 42.86081 0.31007496
```
