context("Previous Issues in ebpmf::ebpmf_point_gamma")
data = readRDS("../../data/ebpmf_gamma_issue.Rds")
x = data$x
s = data$s
init_par = point_gamma(pi0 = 0.1,shape = 10,scale= 1/10)
out = ebpm::ebpm_point_gamma(x, s, g_init = init_par)
expect_false(is.nan(out$log_likelihood))



data = readRDS("../../data/ebpmf_point_issue3.Rds")
x = data$x
s =  data$s
out = ebpm_point_gamma(x, s)
expect_false(is.nan(out$log_likelihood))


data = readRDS("../../data/ebpmf_point_issue4.Rds")
x = data$x
s =  data$s
out = ebpm_point_gamma(x, s)
expect_false(is.nan(out$log_likelihood))


data = readRDS("../../data/ebpmf_point_issue5.Rds")
x = data$x
s =  data$s
out1 = ebpm_point_gamma(x, s, point_gamma(0.5,1,1))
out2 = ebpm_point_gamma(x, s, point_gamma(0.5,1.64,1))
expect_equal(out1$log_likelihood, out2$log_likelihood)







