# # ## test with  data from ebpmf failed cases
# #
# 
# context("Previous Issues in ebpmf::ebpmf_point_gamma")
# data = readRDS("../../data/ebpmf_gamma_issue.Rds")
# x = data$x
# s = data$s
# init_par = c(0.1,10,10)
# #browser()
# out = ebpm::ebpm_point_gamma(x, s, init_par)
# print(out$log_likelihood)
# #expect_false(is.nan(out$log_likelihood))