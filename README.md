# ebpm
R package to fit tthe Empirical Bayes Poisson Means model.

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
