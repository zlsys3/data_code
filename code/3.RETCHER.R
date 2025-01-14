# library(RETCHER)

#load("result/RETCHER/RETCHER_input/inputList.Rdata")
#
# for (i in seq_along(inputList)) {
#   inputList[[i]] <- inputList[[i]][inputList[[i]]$major != 0, ]
# }

# set.seed(1207)
#
# for (i in c(1:61)) {
#   result <-
#     RunPipeline(
#       inputList = inputList[i],
#       saveDir = paste0("result/RETCHER/", names(inputList[i])),
#       clusterMinSnvNum = 3,
#       minimumDepth = 0,
#       removeOutliers = F,
#       maximumClusters = 10,
#       calcCCFMethod = "threshold",
#       error_buf = 0.05,
#       usePercentage = T,
#       LimitCCF = T
#     )
#
#   save(result,
#        file = paste0(
#          "result/RETCHER/",
#          names(inputList[i]),
#          "/",
#          names(inputList[i]),
#          ".Rdata"
#        ))
#
# }


# Distribution box plot of ccf

# load clinical information
clinical <- read.csv("data/wes/clinical.csv")
clinical$metastasis_at_diagnosis[clinical$metastasis_at_diagnosis == "Metastasis, NOS"] <-
  "Metastasis"

# Load CCF data
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

# Get the data frame for drawing ccf box plot
box_df <- data.frame(sample = character(0), ccf = numeric(0))

for (i in seq_along(data_list)) {
  ccf <- data_list[[i]]$result$cluster_res_post$inferTree[, 4]
  sampleName <- rep(names(data_list[i]), length(ccf))
  box_df <-
    rbind(box_df, data.frame(sample = sampleName, ccf = ccf))
}

df <- read.csv("Figure2.csv")

box_df <- merge(box_df, df, by = "sample")[, 1:3]

box_df <-
  box_df[order(box_df$metastasis_at_diagnosis, decreasing = T),]

box_df$sample <-
  factor(box_df$sample, levels = unique(box_df$sample))

library(ggplot2)

ggplot(box_df, aes(x = sample, y = ccf)) +
  geom_boxplot(
    aes(fill = metastasis_at_diagnosis, color = metastasis_at_diagnosis),
    outliers = F,
    outlier.alpha = 1,
    outlier.fill = "black",
    outlier.color = "black",
    outlier.shape = 16,
    outlier.size = 0.6,
    size = 0.8,
    width = 0.7,
    alpha = 0.1
  ) +
  geom_jitter(
    position = position_jitter(0.3),
    alpha = 1,
    size = 0.6,
    shape = 16
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 14),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02)),
    breaks = seq(0, 1, 0.2),
    labels = seq(0, 1, 0.2)
  )


# Draw clone cluster histogram
df <- read.csv("Figure2.csv")

table(df[df$metastasis_at_diagnosis == "No Metastasis",]$cloneNum)

table(df[df$metastasis_at_diagnosis != "No Metastasis",]$cloneNum)

df <-
  data.frame(
    clusterNum = 1:7,
    no_met = c(12, 10, 5, 5, 4, 2, 2),
    met = c(7, 9, 0, 1, 1, 1, 0)
  )

library(tidyverse)

df <-
  df %>% pivot_longer(cols = c(no_met:met),
                      names_to = 'met',
                      values_to = 'sampleNum')

library(ggplot2)
library(ggpubr)

ggplot(df, aes(x = clusterNum,
               y = sampleNum,
               fill = met)) +
  geom_bar(
    stat = 'identity',
    width = 0.7,
    colour = 'white',
    position = 'dodge'
  ) +
  theme_pubclean() +
  scale_fill_manual(values = c("#246b93", "#F39B7F")) +
  labs(x = "clusterNum",
       y = "sampleNum",
       fill = "Metastasis") +
  scale_x_continuous(
    limits = c(0.5, 7.5),
    breaks = seq(1, 7, 1),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  scale_y_continuous(
    limits = c(0, 12.5),
    breaks = seq(0, 12, 2),
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )