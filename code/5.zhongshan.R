library(maftools)
library(tidyverse)

# snv
site <-
  list.files("D:/task/zhongshanOS/funcotator/", full.names = T)
site <- site[grep("*T.tsv", site)]
site
snvlist <- list()
for (i in site) {
  snvlist <- append(snvlist, list(read.csv(i, sep = "\t")))
}
names(snvlist) <-
  sub(".*/(.*?)\\..*", "\\1", site[seq(1, length(site))])

# tumor sample barcode
for (i in seq_along(snvlist)) {
  snvlist[[i]]$Tumor_Sample_Barcode <-  names(snvlist)[i]
}

# filter
# for (i in seq_along(snvlist)) {
#   snvlist[[i]] <-   snvlist[[i]][snvlist[[i]]$t_alt_count >= 5,]
# }

for (i in seq_along(snvlist)) {
  snvlist[[i]] <- snvlist[[i]] |> filter(
    Variant_Classification %in% c(
      "Missense_Mutation",
      "Nonsense_Mutation",
      "Frame_Shift_Del",
      "Frame_Shift_Ins",
      "In_Frame_Del",
      "In_Frame_Ins",
      "Splice_Site"
    )
  )
}


# cnv
site <-
  list.files("D:/task/zhongshanOS/facets/", full.names = T)
site <- site[seq(2, 52, 2)]
site
cnvlist <- list()
for (i in site) {
  cnvlist <- append(cnvlist, list(read.csv(i, sep = "\t")))
}
names(cnvlist) <-
  sub(".*/(.*?)\\..*", "\\1", site[seq(1, length(site))])

cnvlist <- lapply(cnvlist, function(x) {
  x <-
    separate(
      data = x,
      col = INFO,
      into = paste0("INFO", 1:14),
      sep = ";"
    )
})

cnvlist <- lapply(cnvlist, function(x) {
  x <- x[, c(1, 2, 10, 19, 20)]
})

cnvlist <- lapply(cnvlist, function(x) {
  x %>% mutate_all(~ str_replace_all(., "END=", ""))
})

cnvlist <- lapply(cnvlist, function(x) {
  x %>% mutate_all(~ str_replace_all(., "TCN_EM=", ""))
})

cnvlist <- lapply(cnvlist, function(x) {
  x %>% mutate_all(~ str_replace_all(., "LCN_EM=", ""))
})

for (i in seq_along(cnvlist)) {
  cnvlist[[i]]$INFO13[cnvlist[[i]]$INFO13 == "."] <- "0"
}

for (i in seq_along(cnvlist)) {
  cnvlist[[i]]$POS <- as.numeric(cnvlist[[i]]$POS)
  cnvlist[[i]]$INFO3 <- as.numeric(cnvlist[[i]]$INFO3)
  cnvlist[[i]]$INFO12 <- as.numeric(cnvlist[[i]]$INFO12)
  cnvlist[[i]]$INFO13 <- as.numeric(cnvlist[[i]]$INFO13)
}

for (i in seq_along(cnvlist)) {
  cnvlist[[i]]$INFO12 <- cnvlist[[i]]$INFO12 - cnvlist[[i]]$INFO13
}

for (i in seq_along(cnvlist)) {
  colnames(cnvlist[[i]]) <-
    c("chr", "start", "end", "major", "minor")
}


# RETCHER
library(fuzzyjoin)

datalist <- snvlist

for (i in seq_along(datalist)) {
  datalist[[i]] <- fuzzy_left_join(
    datalist[[i]],
    cnvlist[[i]],
    by = c(
      "Chromosome" = "chr",
      "Start_Position" = "start",
      "Start_Position" = "end"
    ),
    match_fun = list(`==`, `>=`, `<=`)
  ) %>%
    select(Chromosome,
           Start_Position,
           Hugo_Symbol,
           t_ref_count,
           t_alt_count,
           major,
           minor)
}

for (i in seq_along(datalist)) {
  datalist[[i]]$purity <- 1
}

for (i in seq_along(datalist)) {
  colnames(datalist[[i]]) <-
    c("chr", "pos", "gene", "ref", "alt", "major", "minor", "purity")
}

for (i in seq_along(datalist)) {
  datalist[[i]][is.na(datalist[[i]])] <- 1
}


# purity
p <-
  c(
    0.32,
    0.32,
    0.62,
    0.23,
    0.7,
    0.39,
    0.42,
    0.41,
    0.35,
    0.19,
    0.47,
    0.39,
    0.74,
    0.62,
    0.37,
    0.64,
    0.3,
    0.54,
    0.35,
    0.31
  )

median(p)

datalist$OS2T$purity <- 0.32
datalist$OS3T$purity <- 0.39
datalist$OS4T$purity <- 0.32
datalist$OS7T$purity <- 0.62
datalist$OS8T$purity <- 0.23
datalist$OS9T$purity <- 0.39
datalist$OS10T$purity <- 0.7
datalist$OS11T$purity <- 0.39
datalist$OS12T$purity <- 0.39
datalist$OS13T$purity <- 0.42
datalist$OS14T$purity <- 0.41
datalist$OS15T$purity <- 0.35
datalist$OS16T$purity <- 0.19

# for (i in seq_along(datalist)) {
#   write.csv(datalist[[i]],
#             paste("data/zhongshan/RETCHER/",
#                   names(datalist)[i],
#                   ".csv",
#                   sep = ""),
#             row.names = F)
# }

inputList <- datalist
save(inputList, file = "inputList_no_filter.Rdata")


# RETCHER
library(RETCHER)

load("inputList_no_filter.Rdata")

for (i in seq_along(inputList)) {
  inputList[[i]] <- inputList[[i]][inputList[[i]]$major != 0, ]
}

set.seed(1207)

for (i in c(1:13)) {
  result <-
    RunPipeline(
      inputList = inputList[i],
      saveDir = paste0("result/zhongshan/RETCHER/filter/", names(inputList[i])),
      clusterMinSnvNum = 3,
      minimumDepth = 0,
      removeOutliers = F,
      maxCCFValue = 1.2,
      maximumClusters = 10,
      useSexChrs = T,
      calcCCFMethod = "threshold",
      usePercentage = T,
      LimitCCF = T
    )
  
  save(result,
       file = paste0(
         "result/zhongshan/RETCHER/filter/",
         names(inputList[i]),
         "/",
         names(inputList[i]),
         ".Rdata"
       ))
}

x <- lapply(inputList, function(x) {
  x <- unique(x$purity)
})

write.table(unlist(x), file = "result/zhongshan/pyclone/input/purity.txt")

# maftools
maf <- do.call(rbind, snvlist)

laml = read.maf(maf = maf)

plotmafSummary(
  maf = laml,
  rmOutlier = T,
  dashboard = T,
  titvRaw = F,
  addStat = 'median',
  log_scale = F,
  showBarcodes = T,
  color = vc_cols,
  fs = 1,
  textSize = 1,
  titleSize = c(1.2, 1),
  top = 10
)

oncoplot(
  maf = laml,
  top = 30,
  showTumorSampleBarcodes = F,
  bgCol = "grey90",
  borderCol = "white",
  gene_mar = 6.3,
  legendFontSize = 1.4
)

titv(maf = laml, plot = T, useSyn = TRUE)

# cnv segment
#BiocManager::install("GenVisR")
library(GenVisR)

site <-
  list.files("D:/task/zhongshanOS/cnvkit/seg/", full.names = T)
site
cnvlist <- list()
for (i in site) {
  cnvlist <- append(cnvlist, list(read.csv(i, sep = "\t")))
}
names(cnvlist) <-
  sub(".*/(.*?)\\..*", "\\1", site[seq(1, length(site))])

for (i in seq_along(cnvlist)) {
  cnvlist[[i]] <- cnvlist[[i]][, -5]
}

cnv <- do.call(rbind, cnvlist)

colnames(cnv) <-
  c("sample", "chromosome", "start", "end", "segmean")
cnv <- cnv[, c(2:5, 1)]

cnSpec(cnv, genome = "hg38", CNscale = "relative")