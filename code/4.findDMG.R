library(maftools)
library(tidyverse)

# Read in data
maf <- read.csv("data/wes/maf.csv")

cnv <- read.csv("data/wes/cnv.csv")

clinical <- read.csv("data/wes/clinical.csv")
clinical$metastasis_at_diagnosis[clinical$metastasis_at_diagnosis == "Metastasis, NOS"] <-
  "Metastasis"

# copy number variation laml
cnTable <-
  merge(maf[, c(1, 5, 6, 16)],
        cnv,
        by = c("Tumor_Sample_Barcode", "Chromosome"),
        all.x = T)

cnTable <- cnTable %>%
  filter(Start_Position >= Start & Start_Position <= End)

cnTable$CN <-
  ifelse(cnTable$Copy_Number == 2,
         "None",
         ifelse(cnTable$Copy_Number > 2,
                "Amp",
                "Del"))

cnTable <- cnTable[, c(1, 3, 10)]
cnTable <- cnTable |> filter(CN != "None") |> distinct()
cnTable <- cnTable[, c(2, 1, 3)]

# Extract non-synonymous mutations
maf <- maf |> filter(
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

# Grouping based on whether transfer occurs at diagnosis
no_met_barcode <-
  clinical[clinical$metastasis_at_diagnosis == "No Metastasis", ]$Tumor_Sample_Barcode
no_met_maf <- maf[maf$Tumor_Sample_Barcode %in% no_met_barcode, ]
met_maf <- maf[!maf$Tumor_Sample_Barcode %in% no_met_barcode, ]

# met_gene <-
#   setdiff(unique(met_maf$Hugo_Symbol), unique(no_met_maf$Hugo_Symbol))
# met_maf <- met_maf[met_maf$Hugo_Symbol %in% met_gene, ]

maf <- rbind(no_met_maf, met_maf)

laml = read.maf(maf = maf,
                cnTable = cnTable,
                clinicalData = clinical)

# waterfall
vc_cols = c(
  "#2372a9",
           "#E84C42",
           "#62A97E",
           "#e0c092",
           "#f8d7e7",
           "#ee7c6f",
           "#91baC2",
           "#6f6baa",
           "grey30",
           "#435C68"
           
)

names(vc_cols) = c(
  'Missense_Mutation',
  'Nonsense_Mutation',
  'Frame_Shift_Ins',
  'Frame_Shift_Del',
  'In_Frame_Del',
  'Splice_Site',
  "Amp",
  "Del",
  'Multi_Hit',
  "Complex_Event"
)

# left bar VAF
genes_vaf <-
  data.frame(
    Hugo_Symbol = maf$Hugo_Symbol,
    VAF = maf$t_alt_count / (maf$t_alt_count + maf$t_ref_count)
  )
genes_vaf <- genes_vaf %>%
  group_by(Hugo_Symbol) %>%
  summarise(VAF = mean(VAF, na.rm = TRUE))

# metastasis annotations
met_colors = c("#8AB2D1", "#C2A7EB")
names(met_colors) = c("Metastasis", "No Metastasis")
status_colors = c("#C2DF8A", "#FB9A99", "#FDBF6F")
names(status_colors) = c("Alive", "Dead", "Not Reported")

anno_colors = list(metastasis_at_diagnosis = met_colors, vital_status =
                     status_colors)

# plot
oncoplot(
  maf = laml,
  top = 30,
  colors = vc_cols,
  # minMut = 2,
  drawRowBar = T,
  drawColBar = T,
  draw_titv = F,
  leftBarData = genes_vaf,
  leftBarLims = c(0, 1),
  clinicalFeatures = c("metastasis_at_diagnosis", "vital_status"),
  sortByAnnotation = T,
  annotationColor = anno_colors,
  # pathways = "smgbp", # sigpw,smgbp
  # topPathways = 5
  showTumorSampleBarcodes = F,
  gene_mar = 6.3,
  anno_height = 1,
  legend_height = 4,
  bgCol = "grey90",
  borderCol = "white",
  annoBorderCol = NULL,
  drawBox = F,
  fontSize = 0.8,
  legendFontSize = 1.4,
  annotationFontSize = 1.4,
  writeMatrix = F
)

# TCGA TMB
tcgaCompare(
  laml,
  cohortName = "OS",
  logscale = T,
  capture_size = 42.1
)

# TMB
tmb(laml, captureSize = 42.1, logScale = F)

# mutational signature
library(BSgenome.Hsapiens.UCSC.hg38)
library(sigminer)

mats <- sig_tally(laml,
                  ref_genome = "BSgenome.Hsapiens.UCSC.hg38",
                  use_syn = F,
                  mode = "SBS")

show_catalogue(mats$all_matrices$SBS_96 %>% t(),
               mode = "SBS",
               style = "cosmic")

est <- sig_estimate(mats$all_matrices$SBS_96,
                    range = 2:5,
                    nrun = 30)

sigs <- sig_auto_extract(mats$all_matrices$SBS_96)
p <- show_sig_profile(sigs, mode = "SBS", style = "cosmic")
sim <- get_sig_similarity(sigs, sig_db = "SBS")
add_labels(
  p,
  x = 0.75,
  y = 0.78,
  y_end = 0.43,
  n_label = 2,
  labels = sim
)

# Mutation gene correlation
library(stringr)

rel_g <-
  somaticInteractions(
    maf = laml,
    top = 50,
    showSum = F,
    fontSize = 0.8,
    nShiftSymbols = 2.5,
    sigSymbolsSize = 2,
    sigSymbolsFontSize = 0.8
  )

# Differential mutation analysis
met <- read.maf(maf = met_maf)
no_met <- read.maf(maf = no_met_maf)

met_compare <- mafCompare(
  m1 = no_met,
  m2 = met,
  m1Name = "no_Met",
  m2Name = "Met",
  pseudoCount = T,
  minMut = 2
)

met_compare

forestPlot(
  mafCompareRes = met_compare,
  pVal = 1,
  color = c('maroon', 'royalblue'),
  geneFontSize = 0.8,
  titleSize = 1.1,
  lineWidth = 1.1
)

# Extract differential mutations
maf1 <-
  maf[maf$Hugo_Symbol %in% met_compare$results[met_compare$results$pval < 1, ]$Hugo_Symbol,]

# load CCF
site <-
  list.files(
    path = "result/RETCHER",
    pattern = "\\.Rdata$",
    recursive = TRUE,
    full.names = TRUE
  )
site <- site[-1]
data_list <- lapply(site, function(x) {
  env <- new.env()
  load(x, envir = env)
  as.list(env)
})
names(data_list) <-
  sub(".*/(.*?)\\..*", "\\1", site[seq(1, length(site))])

# match CCF
maf1$ccf <- NA

for (i in 1:nrow(maf1)) {
  sample_barcode <- maf1$Tumor_Sample_Barcode[i]
  gene_name <- maf1$Hugo_Symbol[i]
  
  sample_data <- data_list[[sample_barcode]]
  
  if (!is.null(sample_data)) {
    mat <- sample_data$result$cluster_res_post$mat
    matching_row <- mat[mat$gene == gene_name, ]
    if (nrow(matching_row) > 0) {
      maf1$ccf[i] <-
        matching_row[matching_row[, 4] == maf1$t_alt_count[i],][, 6]
      
    }
  }
}

# Extract key columns
ccf <-
  maf1[, c("Hugo_Symbol",
           "Tumor_Sample_Barcode",
           "ccf")]

# Convert to a data frame in which sample is the row,
# gene is the column, and ccf is the value.
df <- ccf %>%
  pivot_wider(names_from = Hugo_Symbol,
              values_from = ccf,
              values_fn = mean) %>%
  mutate(across(everything(), ~ replace_na(., 0)))

df <- as.data.frame(df)

Group <-
  data.frame(sample = c(
    unique(met_maf$Tumor_Sample_Barcode),
    unique(no_met_maf$Tumor_Sample_Barcode)
  ),
  Group = c(rep("met", length(
    unique(met_maf$Tumor_Sample_Barcode)
  )), rep("no_met", length(
    unique(no_met_maf$Tumor_Sample_Barcode)
  ))))

df <-
  merge(df, Group, by.x = "Tumor_Sample_Barcode", by.y = "sample")

df <- df[order(df$Group),]

rownames(df) <- df[, 1]
group <- df$Group
df <- df[,-c(1, ncol(df))]

library(pheatmap)

annotation_col = data.frame(group = group)
rownames(annotation_col) = rownames(df)
ann_colors = list(group = c(met = '#7AC5CD',
                            no_met = '#CDB5CD'))

pheatmap(
  t(df),
  color = colorRampPalette(c("grey99", "#CC2727"))(100),
  show_rownames = T,
  show_colnames = F,
  cluster_rows = T,
  cluster_cols = F,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  fontsize = 9,
  angle_col = "90",
  #legend_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)
)


# lasso
library(glmnet)

# metastasis_at_diagnosis
df$Tumor_Sample_Barcode <- rownames(df)
os <-
  clinical[, c("Tumor_Sample_Barcode", "metastasis_at_diagnosis")]
lasso_data <- merge(df, os, by = "Tumor_Sample_Barcode")
lasso_data <-
  lasso_data[, c(1, ncol(lasso_data), 2:(ncol(lasso_data) - 1))]
rownames(lasso_data) <- lasso_data[, 1]
lasso_data <- lasso_data[, -1]
colnames(lasso_data)[1] <- "met"
lasso_data$met <- ifelse(lasso_data$met == "Metastasis", 1, 0)

x <- as.matrix(lasso_data[, c(2:ncol(lasso_data))])
y <- as.matrix(lasso_data$met)

lasso_fit <- glmnet(x,
                    y,
                    alpha = 1,
                    family = "binomial",
                    #lambda = c(0.03408),
                    nlambda = 100)

plot(lasso_fit, xvar = "lambda", label = TRUE)

# feature <-
#   as.data.frame(as.matrix(coef(lasso_fit)))
# feature <-  feature %>% filter(abs(s0) > 0)
# feature

# cross validation
set.seed(1614)
lasso_fit_cv <-
  cv.glmnet(
    x,
    y,
    nfolds = 10,
    type.measure = "mse",
    alpha = 1,
    family = "binomial",
    nlambda = 100
  )
lasso_fit_cv
plot(lasso_fit_cv)

# Extract features filtered by lasso
feature <-
  as.data.frame(as.matrix(coef(lasso_fit_cv, s = lasso_fit_cv$lambda.min)))
feature <-  feature %>% filter(abs(s1) > 0)
feature

# Save key genes
feature$genes <- rownames(feature)
feature <- feature[, c(2, 1)]
write.csv(feature, file = "result/feature/keygenes.csv", row.names = F)