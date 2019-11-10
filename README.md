# ebpm
R package to fit tthe Empirical Bayes Poisson Means model.
See model details and derivation in https://github.com/stephenslab/ebpm/blob/master/derivations/ebpm.pdf


## Quick Start

Install the ebpm package using `devtools`:

```R
library(devtools)
install_github("stephenslab/ebpm")
```

```R
beta = c(rep(0,50),rexp(50))
x = rpois(100,beta) # simulate Poisson observations
s = replicate(100,1)
out = ebpm_exponential_mixture(x,s,m=2)
plot(x, out$posterior$mean)
```
vignettes: 
https://zihao12.github.io/ebpmf_demo/vignette_ebpm.html
