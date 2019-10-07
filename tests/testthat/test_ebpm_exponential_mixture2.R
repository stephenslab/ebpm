# ## test with  data from ebpmf failed cases
# 

context("Previous Issues in ebpmf::ebpmf_exponential_mixture_experiment")
data = readRDS("../../data/ebpmf_issue.Rds")
x = data$x
s = data$s
out = ebpm::ebpm_exponential_mixture(x, s)
expect_false(is.nan(out$log_likelihood))