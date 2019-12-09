# ebpm
R package to fit the Empirical Bayes Poisson Means model:
\begin{equation}
 x_i \sim Pois(\lambda_i)
\end{equation}

See model details and derivation in https://zihao12.github.io/ebpmf_demo/ebpm.pdf


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
## there are the following options: `ebpm_exponential_mixture`, `ebpm_gamma_mixture_single_scale`,`ebpm_point_gamma`, `ebpm_two_gamma`
out = ebpm_exponential_mixture(x,s,m=2)
plot(x, out$posterior$mean)
```
vignettes: 
https://zihao12.github.io/ebpmf_demo/vignette_ebpm.html

a comparison with the method proposed in "Bayesian inference on quasi-sparse count data": https://zihao12.github.io/ebpmf_demo/compare_GH.html

