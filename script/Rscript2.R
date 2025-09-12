library(Rtsne)
library(dbscan)
library(ggplot2)
library(DescTools)

cleanMatrixForClusterW <- function(mtx, f_row = 0.5, f_col = 0.5) {
  cat(sprintf("Filter rows with >%1.2f missingness and columns with >%1.2f missingness.\n",
              f_row, f_col))
  cat("Before: ", nrow(mtx), "rows and ", ncol(mtx),"columns.\n")
  namtx = is.na(mtx)
  good_row = rowSums(namtx) <= ncol(mtx) * (1-f_row)
  good_col = colSums(namtx) <= nrow(mtx) * (1-f_col)
  cat("After: ", sum(good_row), "rows and ", sum(good_col),"columns.\n")
  mtx[good_row, good_col]
}

imputeRowMean = function(mtx) {
  k <- which(is.na(mtx), arr.ind=TRUE)
  mtx[k] <- rowMeans(mtx, na.rm=TRUE)[k[,1]]
  mtx
}

## 5modC, unsupervised clustering
betas <- readRDS("/path/betas_5modC_interpolated_prep.rds")
tsne = Rtsne(t(betas), dims=2, perplexity=20)
df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")
df$barcode = colnames(betas)

df <- readRDS("~/Downloads/data/5modC_tSNE.rds") # the file is in the data directory
db <- dbscan(df[,1:2], eps=1.1, minPts=5)
df$cluster <- db$cluster

ggplot(df) + geom_point(aes(tSNE1, tSNE2, color=as.factor(cluster)))
ggplot(df) + geom_point(aes(tSNE1, tSNE2, color=as.factor(Tissue)))

meta = read_excel("~/Downloads/data/MouseArrayMaster.xlsx", sheet = "Data") %>%
  dplyr::filter(bACE=="0")%>%  dplyr::select(Prep, IDAT, Tissue, AgeInWeeks, Sex)

df$Tissue <- meta$Tissue[match(df$barcode, meta$Prep)]
df$Sex <- meta$Sex[match(df$barcode, meta$Prep)]
df$AgeInWeeks <- meta$AgeInWeeks[match(df$barcode, meta$Prep)]
df$AgeInWeeks <- as.numeric(df$AgeInWeeks)
uncertcoefs <- sapply(c("Tissue", "Sex", "AgeInWeeks"), function(x) UncertCoef(df$cluster, df[[x]], direction="row"))

dfuc <- data.frame(uc = uncertcoefs, variable = names(uncertcoefs))
dfuc$variable = factor(dfuc$variable, levels=dfuc$variable[order(dfuc$uc)])

## 5hmC, unsupervised clustering
betas <- readRDS("/path/betas_5hmC_interpolated_prep.rds")
tsne = Rtsne(t(betas), dims=2, perplexity=20)
df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")
df$barcode = colnames(betas)

df <- readRDS("~/Downloads/data/5hmC_tSNE.rds") # the file is in the data directory
db <- dbscan(df[,1:2], eps=1.1, minPts=5)
df$cluster <- db$cluster

ggplot(df) + geom_point(aes(tSNE1, tSNE2, color=as.factor(cluster)))
ggplot(df) + geom_point(aes(tSNE1, tSNE2, color=as.factor(Tissue)))

meta = read_excel("~/Downloads/data/MouseArrayMaster.xlsx", sheet = "Data") %>%
  dplyr::filter(bACE=="0")%>%  dplyr::select(Prep, IDAT, Tissue, AgeInWeeks, Sex)

df$Tissue <- meta$Tissue[match(df$barcode, meta$Prep)]
df$Sex <- meta$Sex[match(df$barcode, meta$Prep)]
df$AgeInWeeks <- meta$AgeInWeeks[match(df$barcode, meta$Prep)]
df$AgeInWeeks <- as.numeric(df$AgeInWeeks)
uncertcoefs <- sapply(c("Tissue", "Sex", "AgeInWeeks"), function(x) UncertCoef(df$cluster, df[[x]], direction="row"))

dfuc <- data.frame(uc = uncertcoefs, variable = names(uncertcoefs))
dfuc$variable = factor(dfuc$variable, levels=dfuc$variable[order(dfuc$uc)])

## 5mC, unsupervised clustering
betas <- readRDS("/path/betas_5mC_substracted_prep.rds")
tsne = Rtsne(t(betas), dims=2, perplexity=20)
df = as.data.frame(tsne$Y)
colnames(df) = c("tSNE1", "tSNE2")
df$barcode = colnames(betas)

df <- readRDS("~/Downloads/data/5mC_tSNE.rds") # the file is in the data directory
db <- dbscan(df[,1:2], eps=1.1, minPts=5)
df$cluster <- db$cluster

ggplot(df) + geom_point(aes(tSNE1, tSNE2, color=as.factor(cluster)))
ggplot(df) + geom_point(aes(tSNE1, tSNE2, color=as.factor(Tissue)))

meta = read_excel("~/Downloads/data/MouseArrayMaster.xlsx", sheet = "Data") %>%
  dplyr::filter(bACE=="0")%>%  dplyr::select(Prep, IDAT, Tissue, AgeInWeeks, Sex)

df$Tissue <- meta$Tissue[match(df$barcode, meta$Prep)]
df$Sex <- meta$Sex[match(df$barcode, meta$Prep)]
df$AgeInWeeks <- meta$AgeInWeeks[match(df$barcode, meta$Prep)]
df$AgeInWeeks <- as.numeric(df$AgeInWeeks)
uncertcoefs <- sapply(c("Tissue", "Sex", "AgeInWeeks"), function(x) UncertCoef(df$cluster, df[[x]], direction="row"))

dfuc <- data.frame(uc = uncertcoefs, variable = names(uncertcoefs))
dfuc$variable = factor(dfuc$variable, levels=dfuc$variable[order(dfuc$uc)])
