suppressMessages({
  library("sleuth")
})

sample_id <- dir(file.path("myResults"))
sample_id
##[1] "CRI1" "CRI2" "CRI3" "CRI4" "KO1"  "KO2"  "KO3"  "KO4" 

kal_dirs <- file.path("myResults", sample_id, "kalipso")
kal_dirs
 

s2c <- read.table(file.path("my_meta.tsv"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = folder, cond)
s2c
##[1] "myResults/CRI1/kalipso" "myResults/CRI2/kalipso" "myResults/CRI3/kalipso"
##[4] "myResults/CRI4/kalipso" "myResults/KO1/kalipso"  "myResults/KO2/kalipso" 
##[7] "myResults/KO3/kalipso"  "myResults/KO4/kalipso" 

s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)
 sample cond                   path
##1   CRI1  CRI myResults/CRI1/kalipso
##2   CRI2  CRI myResults/CRI2/kalipso
##3   CRI3  CRI myResults/CRI3/kalipso
##4   CRI4  CRI myResults/CRI4/kalipso
##5    KO1   KO  myResults/KO1/kalipso
##6    KO2   KO  myResults/KO2/kalipso
##7    KO3   KO  myResults/KO3/kalipso
##8    KO4   KO  myResults/KO4/kalipso

#here
so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

# fit the full model
so <- sleuth_fit(so, ~cond, 'full')

#What this has accomplished is to “smooth” the raw kallisto abundance estimates for each sample using a linear model with a parameter that represents the experimental condition (in this case scramble vs. HOXA1KD). To test for transcripts that are differential expressed between the conditions, sleuth performs a second fit to a “reduced” model that presumes abundances are equal in the two conditions. To identify differential expressed transcripts sleuth will then identify transcripts with a significantly better fit with the “full” model.

#In the typical Sleuth's workflow, the likelihood ratio test (LRT) is applied. Briefly, the LRT models the likelihood of the data given 2 models:

#full: transcript abundance affected on one or more dependent variables (here just being treated or not)
#reduced: transcript abundance unaffected by the treatment (null hypothesis)
#The LRT then estimates for each transcript the ratio of the 2 likelihoods and produces a q-value (i.e. p-value adjusted for multiple testing by means of false discovery rate, FDR) which can be used as measure of significance. Note that the LRT does not produce any metric equivalent to the fold change, which indicates the maginitude of the change in expression between the 2 conditions and is commonly reported in DEA. That's why the script below also applies the Wald test (WT), which is somewhat related to the LRT and is also used to test for differential expression. WT is used becase it generates the beta statistic, which approximates to the log2 fold change in expression between the 2 condition tested. However, LRT is considerd a better test than the WT (see here) and thus significance filtering is based on LRT's q-values.

#The reduced model is fit with

so <- sleuth_fit(so, ~1, 'reduced')

# then the test
so <- sleuth_lrt(so, 'reduced', 'full')

# can view model with 
models(so)
##[  full  ]
##formula:  ~cond 
##coefficients:
##	(Intercept)
## 	condKO
##[  reduced  ]
##formula:  ~1 
##coefficients:
##	(Intercept)
#
#HERE
# finally view the results
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
head(sleuth_significant, 20)

write.csv(sleuth_table, 'sleuth.csv')
write.csv(sleuth_significant, 'sleuth_significant.csv')

#
plot_bootstrap(so, "ENST00000301480.4", units = "est_counts", color_by = "cond")

# add gene names from ensembl
#source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
library('biomaRt')

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# add these to the tabl
so <- sleuth_prep(s2c, target_mapping = t2g)

so <- sleuth_fit(so, ~cond, 'full')
## fitting measurement error models
## shrinkage estimation
## computing variance of betas
so <- sleuth_fit(so, ~1, 'reduced')
## fitting measurement error models
## shrinkage estimation
## computing variance of betas
so <- sleuth_lrt(so, 'reduced', 'full')
results_table_lrt <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
#head(sleuth_significant, 20)

#write.csv(sleuth_significant, 'sleuth_significant_withnames.csv')

#plot_pca(so, color_by = 'cond')

# add log fold changes by running wald test
# the treatment labe is 'cond' and the non-control value is KO
# combine these to be condKO, doh (f'n R)
so <- sleuth_wt(so, 'condKO')
results_table_wt <- sleuth_results(so, 'condKO')

d5.lrt.sig_ids <- results_table_lrt$target_id[which(results_table_lrt$qval < 0.05)]
d5.wt.sig_ids <- results_table_wt$target_id[which(results_table_wt$qval < 0.05)]
shared_ids <- d5.wt.sig_ids[d5.wt.sig_ids %in% d5.lrt.sig_ids]

#Grab just the rows of the Wald test that correspond to the shared transcripts, then write them to a file for posterity, downstream analysis in another program, etc..

shared_results <- results_table_wt[results_table_wt$target_id %in% shared_ids,]
write.csv(shared_results, file="kallisto_wald_test_lrt_passed_d5.csv")

#To make sure we aren't throwing the baby out with the bathwater, let's write out the Wald-pass-but-LRT-failed genes to another file in case there's something interesting there.

wald_only_ids <- d5.wt.sig_ids[!(d5.wt.sig_ids %in% d5.lrt.sig_ids)]
wald_only_results <- results_table_wt[results_table_wt$target_id %in% wald_only_ids,]
write.csv(wald_only_results, file="kallisto_wald_test_only_passed_d5.csv")

# but we still need to add the ids

library('biomaRt')
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "hsapiens_gene_ensembl",
  host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
    "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
  ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# add these to the tabl
so <- sleuth_prep(s2c, target_mapping = t2g)

so <- sleuth_fit(so, ~cond, 'full')
## fitting measurement error models
## shrinkage estimation
## computing variance of betas
so <- sleuth_fit(so, ~1, 'reduced')
## fitting measurement error models
## shrinkage estimation
## computing variance of betas
so <- sleuth_lrt(so, 'reduced', 'full')
results_table_lrt <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

so <- sleuth_wt(so, 'condKO')
results_table_wt <- sleuth_results(so, 'condKO')

d5.lrt.sig_ids <- results_table_lrt$target_id[which(results_table_lrt$qval < 0.05)]
d5.wt.sig_ids <- results_table_wt$target_id[which(results_table_wt$qval < 0.05)]
shared_ids <- d5.wt.sig_ids[d5.wt.sig_ids %in% d5.lrt.sig_ids]


