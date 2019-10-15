# ##
# data = readRDS("../../data/ebpmf_point_issue3.Rds")
# x = data$x
# s =  data$s
# 
# #browser()
# out = ebpm_point_gamma(x, s)
# 
# 
# ## if use default init, get Inf in one of teh gradient, and result in error
# ## if use c(0.5, 100, 100) as init, get message 
# 
# # Last global step failed to locate a point lower than x.
# # Either x is an approximate local minimum of the function,
# # the function is too non-linear for this algorithm,
# # or steptol is too large.
# 
# ##  affter setting steptol to 1e-20, get Infinite value again. 
# 
# ## Fixed
# 


data = readRDS("../../data/ebpmf_point_issue4.Rds")
x = data$x
s =  data$s

#browser()
out = ebpm_point_gamma(x, s)