library("sleuth")
setwd("/Users/vikrants/Desktop/OIR_New/mir143_mir150_mir126/kallisto/")
base_dir <- "/Users/vikrants/Desktop/OIR_New/mir143_mir150_mir126/kallisto"
sample_id <- dir(file.path(base_dir,"results_1/"))
sample_id


kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "results_1", id, "kallisto"))
kal_dirs



s2c <- read.table(file.path(base_dir, "hiseq_info_1.txt"), header = TRUE, stringsAsFactors=FALSE)
s2c <- dplyr::select(s2c, sample = run_accession, condition)

s2c$condition <- as.factor(s2c$condition)
s2c$condition <- relevel(s2c$condition, "Control")
s2c

s2c <- dplyr::mutate(s2c, path = kal_dirs)
print(s2c)
## design matrix to select reference from the sample
#cond <- factor(s2c$condition)
#cond <- relevel(cond, ref = 'HRPE')
#md <- model.matrix(~cond, s2c)


#so <- sleuth_prep(s2c, ~condition,extra_bootstrap_summary = TRUE,read_bootstrap_tpm=TRUE)
#so <- sleuth_prep(s2c, ~condition, target_mapping = t2g, aggregation_column = 'ens_gene',
#                 extra_bootstrap_summary = TRUE,
#                transformation_function_reads=function(x) log2(x + 0.5),
#               transformation_function_tpm=function(x) log2(x + 0.5),
#              read_bootstrap_tpm = TRUE,  ## required if sleuth_fit uses which_var = "obs_tpm"
#             gene_mode = TRUE) 

#so <- sleuth_fit(so)

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'www.ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id","external_gene_name"), mart = mart)

t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

#so <- sleuth_prep(s2c, ~condition, target_mapping = t2g,extra_bootstrap_summary = TRUE,read_bootstrap_tpm=TRUE)
so <- sleuth_prep(s2c, ~condition, target_mapping = t2g,
                  extra_bootstrap_summary = TRUE,
                  transform_fun_counts=function(x) log2(x + 0.5),
                  transform_fun_tpm=function(x) log2(x + 0.5),
                  read_bootstrap_tpm = TRUE  ## required if sleuth_fit uses which_var = "obs_tpm"
) 

so <- sleuth_fit(so)
#so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_wt(so, 'conditionmiR150')

#so <- sleuth_wt(so)

results_table <- sleuth_results(so, 'conditionmiR150')

results_ordered <- results_table[order(results_table$qval),]

table(results_ordered$qval <= 0.01)
write.table(results_ordered, file = 'Sleuth.DE_Control_vs_miR150.txt', sep="\t",row.names=F, quote=F)
write.table( subset(results_ordered, qval <= 0.01), file='sleuth.DE_transcripts.qval_0.01_control_vs_150_new.txt', sep="\t",row.names=F, quote=F)

#sleuth_live(so)
