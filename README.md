# OLP-vs-OC
Mendelian randomization between Oral lichen planus and Oral Cancer
#install.packages("ieugwasr")
library(ieugwasr)
library(TwoSampleMR)

tok <- "TOKEN"
Sys.setenv(OPEN_GWAS_JWT = tok)
Sys.setenv(OPENGWAS_JWT   = tok)
options(opengwas_jwt       = tok)
options(ieugwasr_api_token = tok)

substring(Sys.getenv("OPEN_GWAS_JWT"), 1, 25)
substring(getOption("opengwas_jwt"), 1, 25)

ieugwasr::gwasinfo(id = "finn-b-L12_LICHENPLANUS") #OLP
ieugwasr::gwasinfo(id = "ebi-a-GCST90025967") #niveles de vit D
ieugwasr::gwasinfo(id = "ieu-b-5132") #oral cancer
ieugwasr::gwasinfo(id = "ieu-b-4834") #ALCOHOL

##Contrastar
exp_IVs <- extract_instruments(
  outcomes = "ieu-b-5132",
  p1 = 5e-8,  
  clump = TRUE, r2 = 0.001, kb = 10000
)
outcome_overlap <- extract_outcome_data(
  snps = exp_IVs$SNP,
  outcomes = "finn-b-L12_LICHENPLANUS",
  proxies = TRUE, rsq = 0.8
)

nrow(outcome_overlap)
head(outcome_overlap)

##Armonizard
dat <- TwoSampleMR::harmonise_data(exp_IVs,outcome_overlap)
head(dat)
#mr analisis 
mr(dat, method_list = c("mr_egger_regression", "mr_ivw", "mr_weighted_median"))

#sensitivity analyses

mr_heterogeneity(dat)
het_results <- mr_heterogeneity(dat)  #para visualizar en tabla y no en consola
View(het_results)

mr_pleiotropy_test(dat)
plei_results <- mr_pleiotropy_test(dat)
View(plei_results)

#single SNP analysis 
##para obtener MR estimado para cada SNP individual
res_single <- mr_singlesnp(dat)
#se obtiene un data frame con resultados del analisis separado para cada par exp-outcome, usando un SNP diferente a la vez
res_single <- mr_singlesnp(dat, single_method = "mr_meta_fixed")

#plot
res <- mr(dat)
scatterplot <- mr_scatter_plot(res, dat)
scatterplot[[1]]
