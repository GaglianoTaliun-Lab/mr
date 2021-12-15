library('TwoSampleMR')

ao <- available_outcomes()
#exposure ebi-a-GCST90002412: LDL cholesterol levels
exposure_dat <- extract_instruments(c('ebi-a-GCST000759'))
exposure_dat <- clump_data(exposure_dat)
#outcome ieu-a-798: myocardial infarction
outcome_dat <- extract_outcome_data(exposure_dat$SNP, c('ieu-a-798'), proxies = 1, rsq = 0.8, align_alleles = 1, palindromes = 1, maf_threshold = 0.3) 
dat <- harmonise_data(exposure_dat, outcome_dat, action = 2)

#Summary MR results
res<-mr(dat)
res

#Analysis by single SNPs
#Default methods: single_method = "mr_wald_ratio", all_method = c("mr_ivw", "mr_egger_regression")
res_single<-mr_singlesnp(dat)
res_single

#Leave-one-out Analysis
#Default method: method = mr_ivw
res_loo<-mr_leaveoneout(dat)
res_loo

#Inverse variance weighted
mr_ivw(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)

#Weighted median
mr_weighted_median(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)

#Weighted mode
mr_weighted_mode(dat$beta.exposure, dat$beta.outcome, dat$se.exposure, dat$se.outcome)

#Plots
#Scatter plot (lines with methods used in res)
pdf("scatter-male.pdf")
mr_scatter_plot(res, dat)
dev.off()

#Forest
pdf("forest.pdf")
mr_forest_plot(res_single)
dev.off()

#Funnel plot
pdf("funnel.pdf")
mr_funnel_plot(res_single)
dev.off()

#Leave-one-out
pdf("leaveoneoute.pdf")
mr_leaveoneout_plot(res_loo)
dev.off()

#Testing for horizontal pleiotropy using MR-Egger intercept test
mr_pleiotropy_test(dat)


