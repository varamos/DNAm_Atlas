### Generate betas using SeSAMe
library(sesame)
library(readxl)
library(dplyr)
library(data.table)
library(tidyr)
library(parallel)

IDs <- read_excel("~/Downloads/scripts/MouseArrayMaster.xlsx", sheet = "Data") %>% dplyr::select(Prep, IDAT) 

pfxes = searchIDATprefixes("~/zhou_lab/projects/20230125_nonGEO_IDATs/MM285/") ## IDAT files location
sub_pfxes <- pfxes[grepl(paste(IDs$IDAT, collapse = "|"), pfxes)]

all(IDs$IDAT %in% names(sub_pfxes))
all(names(sub_pfxes) %in% IDs$IDAT)

betas <- openSesame(sub_pfxes, BPPARAM = BiocParallel::MulticoreParam(20)) ## 265 BS | 265 bACE
saveRDS(betas, "~/path/betas.rds")


IDs_5hmC <- IDs %>% filter(grepl("_A", Prep))
IDs_5modC <- IDs %>% filter(grepl("_B", Prep))

betas_5hmC <- betas[, IDs_5hmC$IDAT]
saveRDS(betas_5hmC, "~/path/betas_5hmC.rds")
betas_5modC <- betas[, IDs_5modC$IDAT]
saveRDS(betas_5modC, "~/path/betas_5modC.rds")

## Methylation standards (GSE184410)
IDs <- read.csv("~/Downloads/data/DNAm_standards.csv")
# download IDAT files from GSE184410
pfxes = searchIDATprefixes("~/zhou_lab/projects/20230125_nonGEO_IDATs/MM285/") ## location IDAT files downloaded
sub_pfxes <- pfxes[IDs$IDAT]
all(IDs$IDAT %in% names(sub_pfxes))
setdiff(names(sub_pfxes), IDs$IDAT)
betas = openSesame(sub_pfxes,BPPARAM = BiocParallel::MulticoreParam(6))
betas_df <- as.data.frame(betas, col.names= names(x)) %>% filter(grepl("cg", rownames(betas))) 
betas_df <- betas_df %>% filter(!grepl("ctl-", rownames(betas_df))) %>% na.omit()
betas_df$probeID <- rownames(betas_df)
m_ref <- reshape2::melt(betas_df)
setDT(m_ref)
m_ref_1 <- merge(m_ref,IDs %>% select(-c(Sample_Type)), by.x="variable", by.y="IDAT")
m_ref_1 <- m_ref_1 %>% select(-variable)
m_ref_1$info <- m_ref_1$info*0.01

## 5hmC, interpolation
result <- m_ref_1 %>%
  pivot_wider(names_from = probeID, values_from = value, values_fn = mean) %>%
  arrange(info)

exp_5hmC <- readRDS("~/path/betas_5hmC.rds") %>% as.data.frame()
exp_5hmC$probeID <- rownames(exp_5hmC)
exp_5hmC <- exp_5hmC %>% filter(probeID %in% unique(m_ref_1$probeID)) %>% na.omit() 
exp_5hmC <- exp_5hmC %>% select(probeID, everything())

result <- result %>% select(info, rownames(exp_5hmC))
setdiff(rownames(exp_5hmC)[1:nrow(exp_5hmC)], colnames(result)[2:ncol(result)]) # check probes order

rownames(exp_5hmC) <- NULL

calculate_corrected_methylation <- function(signal, CpG) {
  approx(result[[CpG]], result$info, xout = signal, rule=2, method = "linear")$y
}

corrected_methylation_values <- mclapply(2:ncol(exp_5hmC), function(sample_col) {
  sapply(1:nrow(exp_5hmC), function(i) {
    CpG <- exp_5hmC$probeID[i]
    signal <- exp_5hmC[i, sample_col]
    calculate_corrected_methylation(signal, CpG)
  })
}, mc.cores = 48)

corrected_methylation_values <- do.call(cbind, corrected_methylation_values)
corrected_methylation_df <- data.frame(corrected_methylation_values)
colnames(corrected_methylation_df) <- paste("corrected_", colnames(exp_5hmC)[2:ncol(exp_5hmC)], sep = "")

final_results <- cbind(exp_5hmC, corrected_methylation_df)
saveRDS(final_results, "~/path/final_results_5hmC.rds")

## 5modC, interpolation
result <- m_ref_1 %>%
  pivot_wider(names_from = probeID, values_from = value, values_fn = mean) %>%
  arrange(info)

exp_5modC <- readRDS("~~/path/betas_5modC.rds") %>% as.data.frame()
exp_5modC$probeID <- rownames(exp_5modC)
exp_5modC <- exp_5modC %>% filter(probeID %in% unique(m_ref_1$probeID)) %>% na.omit() 
exp_5modC <- exp_5modC %>% select(probeID, everything())

result <- result %>% select(info, rownames(exp_5modC))
setdiff(rownames(exp_5modC)[1:nrow(exp_5modC)], colnames(result)[2:ncol(result)]) # check probes order

rownames(exp_5modC) <- NULL

calculate_corrected_methylation <- function(signal, CpG) {
  approx(result[[CpG]], result$info, xout = signal, rule=2, method = "linear")$y
}

corrected_methylation_values <- mclapply(2:ncol(exp_5modC), function(sample_col) {
  sapply(1:nrow(exp_5modC), function(i) {
    CpG <- exp_5modC$probeID[i]
    signal <- exp_5modC[i, sample_col]
    calculate_corrected_methylation(signal, CpG)
  })
}, mc.cores = 48)

corrected_methylation_values <- do.call(cbind, corrected_methylation_values)
corrected_methylation_df <- data.frame(corrected_methylation_values)
colnames(corrected_methylation_df) <- paste("corrected_", colnames(exp_5modC)[2:ncol(exp_5modC)], sep = "")

final_results <- cbind(exp_5modC, corrected_methylation_df)
saveRDS(final_results, "~/path/final_results_5modC.rds")

final_results <- readRDS("~/path/final_results_5modC.rds")
rownames(final_results) <- final_results$probeID
final_results1 <- final_results %>% select(-probeID)
final_before <- final_results1[,1:(ncol(final_results1)/2)]
final_after <- final_results1[,(1+(ncol(final_results1)/2)):ncol(final_results1)]
colnames(final_after) <- gsub("^corrected_", "", colnames(final_after))
saveRDS(final_after, "~/path/betas_5modC_interpolated.rds") # 5modC betas

final_results <- readRDS("~/path/final_results_5hmC.rds")
rownames(final_results) <- final_results$probeID
final_results1 <- final_results %>% select(-probeID)
final_before <- final_results1[,1:(ncol(final_results1)/2)]
final_after <- final_results1[,(1+(ncol(final_results1)/2)):ncol(final_results1)]
colnames(final_after) <- gsub("^corrected_", "", colnames(final_after))
saveRDS(final_after, "~/path/betas_5hmC_interpolated.rds") # 5hmC betas

## Generate 5mC profiles (5modC - 5hmC)
IDs <- read_excel("~/Downloads/data/MouseArrayMaster.xlsx", sheet = "Data"
)  %>% dplyr::select(Prep, IDAT, IDAT, Tissue, Sex) 

modC <- readRDS("~/path/betas_5modC_interpolated.rds")
hmC <- readRDS("~/path/betas_5hmC_interpolated.rds")

dictionary_df <- setNames(IDs$Prep, IDs$IDAT)
names(modC) <- dictionary_df[names(modC)]
names(hmC) <- dictionary_df[names(hmC)]

saveRDS(modC, "~/path/betas_5modC_interpolated_prep.rds")
saveRDS(hmC, "~/path/betas_5hmC_interpolated_prep.rds")

modC_only <- setdiff(rownames(modC), rownames(hmC))
hmC_only <- setdiff(rownames(hmC), rownames(modC))
common_rows <- intersect(rownames(modC), rownames(hmC))
modC_1 <- modC[common_rows, , drop = FALSE]
hmC_1 <- hmC[common_rows, , drop = FALSE]

mC_1 <- modC_1 - hmC_1
mC_1 <- apply(mC_1, c(1, 2), function(x) ifelse(is.na(x), NA, ifelse(x < 0, 0, x)))
mC_1 <- mC_1 %>% as.data.frame()
saveRDS(mC_1,"~/path/betas_5mC_substracted_prep.rds")

modC_1$probeID <- rownames(modC_1)
modC1 <- melt(modC_1)
hmC_1$probeID <- rownames(hmC_1)
hmC1 <- melt(hmC_1)
setDT(modC1)
setDT(hmC1)
all <- merge(modC1, hmC1, by=c("probeID", "variable"))

mC_1$probeID <- rownames(mC_1)
mC1 <- melt(mC_1)
setDT(mC1)
rm(modC_1, hmC_1, hmC1, modC1)
final <- merge(all, mC1, by=c("probeID", "variable"))

## add tissue data
IDs_1 <- IDs%>% select(c(Prep, Tissue, Sex)) %>% distinct() %>% na.omit()
colnames(IDs_1)[1] <-"variable"
setDT(final)
setDT(IDs_1)
final_1 <- merge(final, IDs_1, by="variable") %>% na.omit()
colnames(final_1) <- c("sample", "probeID", "modC", "hmC", "mC","Tissue", "Sex")
final_1a <- final_1 %>% filter(grepl("cg", probeID))

### Calculate the mean cytosine modification level for each sample
tissue_type <- final_1a %>% select(c(sample,Tissue)) %>% distinct()
mean_C_modi <- final_1a %>%
  group_by(sample,Tissue) %>%
  dplyr::summarise_at(vars("modC", "hmC", "mC"), ~mean(., na.rm = TRUE))

## Calculate the Spearman correlation coefficient between 5hmC and tissue turnover days
Turnover <- read_excel("~/Downloads/data/turnover.xlsx")
IDs_s <- read_excel("~/Downloads/data/MouseArrayMaster.xlsx", sheet = "Data") %>% dplyr::select(Prep, Tissue, Sex, AgeInWeeks)  
IDs_s <- IDs_s%>% distinct()

mean_hmC_sample <- final_1a %>%
  group_by(sample) %>%
  dplyr::summarise_at(vars("hmC"), ~mean(., na.rm = TRUE))

mean_hmC_sample_1 <- merge(mean_hmC_sample,IDs_s, by.x="sample", by.y="Prep")
mean_hmC_sample_2 <- mean_hmC_sample_1 %>% filter(Tissue %in% Turnover$Tissue)

mean_hmC_sample_3 <- merge(Turnover, mean_hmC_sample_2, by="Tissue")
cor(mean_hmC_sample_3$hmC, mean_hmC_sample_3$`Estimated_turnover (days)`, method="spearman")
# 0.8170366

## Calculate mean DNAm per sample, grouped by chromHMM
IDs_s <- read_excel("~/Downloads/data/MouseArrayMaster.xlsx", sheet = "Data") %>%  
  dplyr::filter(bACE=="1")%>% dplyr::select(Prep, AgeInWeeks)
chromHMM <- fread("~/Downloads/data/ChromHMM.20220318.gz")

setDT(final_1a)
final_chrom <- merge(final_1a, chromHMM, by.x="probeID", by.y="Probe_ID")
final_chrom$Knowledgebase <- gsub("ChromHMM;", "", final_chrom$Knowledgebase)

