library(maftools)
library(tidyverse)

maf <- read.csv("data/wes/maf.csv")

cnv <- read.csv("data/wes/cnv.csv")

clinical <- read.csv("data/wes/clinical.csv")
clinical$metastasis_at_diagnosis[clinical$metastasis_at_diagnosis == "Metastasis, NOS"] <-
  "Metastasis"

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

no_met_barcode <-
  clinical[clinical$metastasis_at_diagnosis == "No Metastasis", ]$Tumor_Sample_Barcode
no_met_maf <- maf[maf$Tumor_Sample_Barcode %in% no_met_barcode, ]
met_maf <- maf[!maf$Tumor_Sample_Barcode %in% no_met_barcode, ]

library(ggvenn)

venn_list <-
  list(Metastasis = met_maf$Hugo_Symbol,
       Non_Metastasis = no_met_maf$Hugo_Symbol)

ggvenn(
  venn_list,
  c("Non_Metastasis", "Metastasis"),
  stroke_color = c('black'),
  stroke_size = 0.5,
  stroke_linetype = 0,
  stroke_alpha = 1,
  #,"#CCCFFF",'#FFCCCC'
  fill_color = c('#D1EBF2', '#F6D1CA'),
  fill_alpha = 0.8,
  set_name_color = "black",
  set_name_size = 8,
  show_percentage = F,
  auto_scale = T,
  text_size = 10,
  text_color = "black",
)

inner_genes <-
  intersect(unique(met_maf$Hugo_Symbol), unique(no_met_maf$Hugo_Symbol))
length(inner_genes)

met_genes <- setdiff(met_maf$Hugo_Symbol, inner_genes)
length(met_genes)

no_met_genes <- setdiff(no_met_maf$Hugo_Symbol, inner_genes)
length(no_met_genes)

inner_maf <- maf[maf$Hugo_Symbol %in% inner_genes, ]
met_maf <- maf[maf$Hugo_Symbol %in% met_genes, ]
no_met_maf <- maf[maf$Hugo_Symbol %in% no_met_genes, ]

# Add cnv to maf and extract key columns
getCNV <- function(maf) {
  cnTable <-
    merge(maf[, c(1, 5, 6, 16)],
          cnv,
          by = c("Tumor_Sample_Barcode", "Chromosome"),
          all.x = T)
  
  cnTable <- cnTable %>%
    filter(Start_Position >= Start & Start_Position <= End)
  
  cnTable <- cnTable[, c(1, 3, 7)]
  cnTable <- cnTable |> distinct()
  
  maf <-
    merge(
      maf,
      cnTable,
      by = c("Tumor_Sample_Barcode", "Hugo_Symbol"),
      all.x = T
    )
  
  maf <- maf[, c(6, 7, 2, 41, 42, 141, 1)]
  
  maf$Copy_Number[is.na(maf$Copy_Number)] <- 2
  
  return(maf)
}

inner_maf <- getCNV(inner_maf)
met_maf <- getCNV(met_maf)
no_met_maf <- getCNV(no_met_maf)

# add purity
load("result/RETCHER/RETCHER_input/inputList.Rdata")

purity_df <-
  data.frame(Tumor_Sample_Barcode = character(61), purity = 0)

for (i in seq_along(inputList)) {
  purity_df$Tumor_Sample_Barcode[i] <- names(inputList)[i]
  purity_df$purity[i] <- unique(inputList[[i]]$purity)
}

inner_maf <-
  merge(inner_maf, purity_df, by = "Tumor_Sample_Barcode", all.x = T)

met_maf <-
  merge(met_maf, purity_df, by = "Tumor_Sample_Barcode", all.x = T)

no_met_maf <-
  merge(no_met_maf, purity_df, by = "Tumor_Sample_Barcode", all.x = T)

# Function to calculate CCF in maf
calcCCF <- function(df, LimitCCF = TRUE) {
  m <- pmax(1, round(((
    df$t_alt_count / (df$t_ref_count + df$t_alt_count)
  ) *
    ((df$Copy_Number * df$purity) + 2 * (1 - df$purity)
    )) / df$purity))
  
  ccf <- ((df$t_alt_count / (df$t_ref_count + df$t_alt_count)) *
            ((df$Copy_Number * df$purity) + 2 * (1 - df$purity))) / (df$purity * m)
  
  df$ccf <- ccf
  df$mutMulti <- m
  df <- df[order(df$Tumor_Sample_Barcode), ]
  
  if (LimitCCF) {
    df$ccf <- ifelse(df$ccf > 1, 1, df$ccf)
  }
  
  return(df)
}

inner_maf <- calcCCF(inner_maf, LimitCCF = T)
met_maf <- calcCCF(met_maf, LimitCCF = T)
no_met_maf <- calcCCF(no_met_maf, LimitCCF = T)

# Calculate weighted average ccf
weight_average_ccf <- function(maf) {
  x <- as.data.frame(
    maf %>%
      group_by(Hugo_Symbol) %>%
      summarize(
        weighted_ccf = sum(ccf * (t_ref_count + t_alt_count)) / sum(t_ref_count + t_alt_count),
        .groups = 'drop'
      )
  )
  return(x)
}

inner_maf_weighted <- weight_average_ccf(inner_maf)
met_maf_weighted <- weight_average_ccf(met_maf)
no_met_maf_weighted <- weight_average_ccf(no_met_maf)

# Draw density plots and scatter plots
library(ggthemes)

plot_Density_Point <- function(weighted) {
  # density
  p1 <- ggplot(weighted, aes(weighted_ccf)) +
    geom_density(color =  "black", fill =  "gray") +
    theme_clean() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    ) +
    labs(x = "Weighted_CCF",
         y = "Density") +
    scale_x_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0.01, 0.01)),
      breaks = seq(0, 1, 0.2),
      labels = seq(0, 1, 0.2)
    )
  
  print(p1)
  
  
  # scatter
  # Randomly sample 250X depth for easy display
  weighted$depth <-
    round(rnorm(nrow(weighted), 250, 100))
  
  p2 <- ggplot(weighted,
               aes(weighted_ccf,
                   depth)) +
    geom_point(size = 2,
               shape = 15,
               alpha = 0.8) +
    theme_classic() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    ) +
    labs(x = "CCF",
         y = "Depth") +
    scale_x_continuous(
      limits = c(0, 1),
      expand = expansion(mult = c(0.01, 0.01)),
      breaks = seq(0, 1, 0.2),
      labels = seq(0, 1, 0.2)
    )
  
  print(p2)
}

plot_Density_Point(inner_maf_weighted)
plot_Density_Point(met_maf_weighted)
plot_Density_Point(no_met_maf_weighted)

# 12 gene evolutionary trees
genes <-
  c(
    "TP53",
    "ATRX",
    "MUC16",
    "RB1",
    "TTN",
    "CNTNAP5",
    "ALMS1",
    "CSMD2",
    "HELZ2",
    "PCDH15",
    "RYR2",
    "ZFHX3"
  )

genes_12_met <-
  met_maf_weighted[met_maf_weighted$Hugo_Symbol %in% genes, ]

genes_12_no_met <-
  no_met_maf_weighted[no_met_maf_weighted$Hugo_Symbol %in% genes, ]

plot_Density_Point(genes_12_met)
plot_Density_Point(genes_12_no_met)

# Queue-based analysis of evolutionary paths of high-frequency mutations------
maf <- read.csv("data/wes/maf.csv")

cnv <- read.csv("data/wes/cnv.csv")

clinical <- read.csv("data/wes/clinical.csv")
clinical$metastasis_at_diagnosis[clinical$metastasis_at_diagnosis == "Metastasis, NOS"] <-
  "Metastasis"

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

no_met_barcode <-
  clinical[clinical$metastasis_at_diagnosis == "No Metastasis", ]$Tumor_Sample_Barcode
no_met_maf <- maf[maf$Tumor_Sample_Barcode %in% no_met_barcode, ]
met_maf <- maf[!maf$Tumor_Sample_Barcode %in% no_met_barcode, ]

inner_genes <-
  intersect(unique(met_maf$Hugo_Symbol), unique(no_met_maf$Hugo_Symbol))
length(inner_genes)
maf <- maf[maf$Hugo_Symbol %in% inner_genes, ]

# Extract mutations that appear in more than or equal to three samples,
# a total of 16, of which 12 are previously shared mutations.
genes <-
  as.data.frame(maf |> group_by(Hugo_Symbol) |>
                  filter(n_distinct(Tumor_Sample_Barcode) >= 3) |>
                  summarise() |>  ungroup())[, 1]

maf <- maf[maf$Hugo_Symbol %in% genes, ]

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

laml = read.maf(maf = maf,
                cnTable = cnTable,
                clinicalData = clinical)

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

# oncoplot
oncoplot(
  maf = laml,
  top = 30,
  colors = vc_cols,
  drawRowBar = F,
  drawColBar = F,
  draw_titv = F,
  #leftBarData = genes_vaf,
  leftBarLims = c(0, 1),
  clinicalFeatures = c("metastasis_at_diagnosis"),
  sortByAnnotation = T,
  annotationColor = anno_colors,
  showTumorSampleBarcodes = F,
  gene_mar = 6.3,
  anno_height = 0.5,
  legend_height = 4,
  bgCol = "grey90",
  borderCol = "white",
  annoBorderCol = NULL,
  drawBox = F,
  fontSize = 0.8,
  legendFontSize = 1.4,
  annotationFontSize = 1.4,
  writeMatrix = T,
)

library(reshape2)

df <-
  readxl::read_xlsx("result/Figure/onco_matrix.xlsx")
df <- melt(df, id = "Gene")
df$Gene <- factor(df$Gene, levels = df$Gene[1:12])

df1 <-
  readxl::read_xlsx("result/Figure/onco_matrix_cnv.xlsx")
df1 <- melt(df1, id = "Gene")
df1$Gene <- factor(df1$Gene, levels = df1$Gene[1:12])

ggplot(df, aes(x = variable, y = Gene)) +
  geom_tile(aes(fill = value),
            color = "grey50",
            linewidth = 0.4) +
  geom_dotplot(
    data = df1,
    aes(variable, Gene, fill = value),
    na.rm = T,
    dotsize = 0.8,
    binwidth = 0.6,
    stackdir = "center"
  ) +
  scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "bottom") +
  theme_minimal() +
  theme(
    panel.border = element_rect(
      fill = NA,
      color = "white",
      size = 1,
      linetype = "solid"
    ),
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    # axis.text.x = element_text(
    #   angle = 45,
    #   colour = 'black',
    #   size = 12,
    #   hjust = 0,
    #   vjust = 0
    # ),
    axis.text.x = element_blank(),
    axis.text.y = element_text(colour = 'black', size = 12),
    #axis.text.y = element_blank(),
    plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), units = , "cm"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.position = "right",
    #legend.margin = margin(5, 5, 5, 5)
  ) +
  scale_fill_manual(
    name = " ",
    values = c(
      "#0072B5CC",
                 "#FFDC91CC",
                 "#E18727CC",
                 "#E4BBBA",
                 "#7876B1CC",
                 "#6F99ADCC",
                 "#DCDCDC",
                 "#B7BFA0",
                 "#8B6969"
    ),
    na.value = "white",
    labels = c(
      "Missense_Mutation",
      "Frame_Shift_Del",
      "In_Frame_Del",
      "Frame_Shift_Ins",
      "Nonsense_Mutation",
      "Splice_Site",
      "Multi_Hit",
      "Amp",
      "Del"
    ),
    breaks = c(
      "Missense_Mutation",
      "Frame_Shift_Del",
      "In_Frame_Del",
      "Frame_Shift_Ins",
      "Nonsense_Mutation",
      "Splice_Site",
      "Multi_Hit",
      "Amp",
      "Del"
    )
  )




# Divide into transfer and non-transfer, respectively infer the evolutionary
# relationship between mutations
no_met_barcode <-
  clinical[clinical$metastasis_at_diagnosis == "No Metastasis", ]$Tumor_Sample_Barcode
no_met_maf <- maf[maf$Tumor_Sample_Barcode %in% no_met_barcode, ]
met_maf <- maf[!maf$Tumor_Sample_Barcode %in% no_met_barcode, ]

# Non-Metastasis
oncoplot(
  read.maf(no_met_maf),
  colors = vc_cols,
  writeMatrix = T,
  drawRowBar = F,
  drawColBar = F
)

no_met_maf <- read.csv("onco_matrix.txt", sep = "\t")
no_met_maf[no_met_maf != ""] <- 1
no_met_maf[no_met_maf == ""] <- 0
rname <- rownames(no_met_maf)
no_met_maf <- as.data.frame(apply(no_met_maf, 2, as.numeric))
rownames(no_met_maf) <- rname

# Metastasis
oncoplot(
  read.maf(met_maf),
  colors = vc_cols,
  writeMatrix = T,
  drawRowBar = F,
  drawColBar = F
)

met_maf <- read.csv("onco_matrix.txt", sep = "\t")
met_maf[met_maf != ""] <- 1
met_maf[met_maf == ""] <- 0
rname <- rownames(met_maf)
met_maf <- as.data.frame(apply(met_maf, 2, as.numeric))
rownames(met_maf) <- rname


# Infer mutational evolutionary relationships
library(cate)
no_met_evo <- mutEvolution(no_met_maf)
met_evo <- mutEvolution(met_maf)

plotGraph(no_met_evo,
          layoutType = "sugiyama",
          main = "No_Met Mutation Evolutionary Trails")

plotGraph(met_evo,
          layoutType = "sugiyama",
          main = "Met Mutation Evolutionary Trails")
