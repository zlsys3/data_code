library(dbscan)
library(ggplot2)

# x is a data frame that must have an id ref alt major minor column
calcPurity <- function(df,
                       mergePvalue = 0.05) {
  # Extract copy number neutral genes
  x <- df[df$major == 1 & df$minor == 1,]
  if (nrow(x) == 0 | nrow(x) == 1) {
    cat("No copy number neutral mutations\n")
    return(0)
  }
  x$vaf <- x$alt / (x$alt + x$ref)
  x$depth <- x$alt + x$ref
  
  # Calculate k-nearest neighbor distance, k = 2*feature dimension-1 = 1
  k_distances <-
    data.frame(x = seq_along(x$vaf),
               y = sort(kNNdist(as.data.frame(x$vaf), 1)))
  
  # Visualize k distance
  plot(k_distances$y,
       type = "l",
       xlab = "Points sorted by distance",
       ylab = "Distance")
  
  # Find the point with the largest gradient
  max_index <-
    which(diff(k_distances$y) == max(diff(k_distances$y)))
  
  # Get the eps corresponding to the inflection point
  eps <- k_distances$y[max_index]
  
  # cat("The inflection point is", max_index, "\n")
  # cat("The eps is", eps, "\n")
  
  abline(h = eps,
         col = "blue",
         lty = 2)
  
  # dbscan
  # minPts = feature dimension + 1 = k + 1 = 2
  dbscan_result <-
    dbscan(as.data.frame(x$vaf),
           eps = eps,
           minPts = 2)
  
  # print(dbscan_result)
  
  x$cluster <- dbscan_result$cluster
  
  # Remove outliers and merge indifference clusters
  x_merged <- mergeCluster(x, pvalue = mergePvalue)
  
  # Twice the vaf of the cloned cluster is the purity
  purity <-
    mean(x_merged[x_merged$cluster == max(x_merged$cluster), ]$vaf) * 2
  
  # cat("\nTumor purity is", purity, "\n")
  
  df$purity <- purity
  
  # Visualization
  p <- ggplot(x, aes(
    x = vaf,
    y = depth,
    color = as.factor(cluster)
  )) +
    geom_point(size = 3) +
    labs(
      title = "DBSCAN Cluster Result",
      x = "VAF",
      y = "Depth",
      color = "Cluster"
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    ) +
    scale_x_continuous(
      limits = 0:1,
      expand = expansion(mult = c(0.01, 0.01)),
      breaks = seq(0, 1, 0.2),
      labels = seq(0, 1, 0.2)
    )
  
  # print(p)
  
  # Merged clustering results
  p <- ggplot(x_merged, aes(
    x = vaf,
    y = depth,
    color = as.factor(cluster)
  )) +
    geom_point(size = 3) +
    labs(
      title = "Processed Cluster Result",
      x = "VAF",
      y = "Depth",
      color = "Cluster"
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12)
    ) +
    scale_x_continuous(
      limits = 0:1,
      expand = expansion(mult = c(0.01, 0.01)),
      breaks = seq(0, 1, 0.2),
      labels = seq(0, 1, 0.2)
    )
  
  vaf_line <-
    mean(x_merged[x_merged$cluster == max(x_merged$cluster),]$vaf)
  
  p <- p + geom_vline(aes(xintercept = vaf_line),
                      linetype = 5,
                      col = "grey50")
  # print(p)
  
  return(df)
}

# Remove outliers and merge indifference clusters
mergeCluster <- function(df, pvalue = 0.05) {
  df <- df[df$cluster != 0, ]
  label <- 0
  for (i in sort(unique(df$cluster), decreasing = T)) {
    for (j in sort(unique(df$cluster), decreasing = T)) {
      if (j < i) {
        if (length(df[df$cluster == i, ]$vaf) != 0 &
            length(df[df$cluster == j, ]$vaf) != 0) {
          wilcox <-
            wilcox.test(df[df$cluster == i, ]$vaf,
                        df[df$cluster == j, ]$vaf,
                        exact = F)
          
          if (wilcox$p.value >= pvalue) {
            df[df$cluster == j, ]$cluster <- i
            # cat(
            #   "Cluster",
            #   i,
            #   "and Cluster",
            #   j,
            #   "have been merged, they are not statistical difference\n"
            # )
            label <- label + 1
          }
        }
      }
    }
  }
  
  if (label == 0) {
    # cat("All clusters are statistically different\n")
  }
  return(df)
}

# load data
load("result/purity/inputList.Rdata")

for (i in seq_along(inputList)) {
  df <- inputList[[i]]
  df <- calcPurity(df)
  if (length(df) == 1) {
    next
  }
  inputList[[i]] <- df
}

purity <- c()
for (i in seq_along(inputList)) {
  x <- unique(inputList[[i]]$purity)
  if (!is.null(x)) {
    purity <- c(purity, x)
  }
}

barplot(purity,
        ylim = c(0, 1),
        ylab = "Tumor Purity")

median(purity)

for (i in seq_along(inputList)) {
  if (ncol(inputList[[i]]) == 7) {
    inputList[[i]]$purity <- median(purity)
  }
}

all_purity <- numeric()
for (i in seq_along(inputList)) {
  all_purity <- c(all_purity, unique(inputList[[i]]$purity))
}

all_purity
barplot(all_purity,
        ylim = c(0, 1),
        ylab = "Tumor Purity")

save(inputList, file = "data/wes/RETCHER/inputList.Rdata")

# demo
df <- inputList[[17]]
df <- df[,-8]
res <- calcPurity(df)


library(ggplot2)

purity_df <-
  data.frame(sample = character(), purity = numeric())

for (i in seq_along(inputList)) {
  sample <- names(inputList)[i]
  purity <- unique(inputList[[i]]$purity)
  if (!is.null(purity)) {
    purity_df <- rbind(purity_df, data.frame(sample, purity))
  }
}

ggplot(purity_df, aes(x = sample, y = purity)) +
  geom_bar(stat = "identity",
           width = 0.8,
           fill = "grey55") +
  theme_classic() +
  theme(
    # legend.position = "right",
    panel.grid.major.x = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  labs(x = "Samples",
       y = "Purity") +
  scale_x_discrete(expand = expansion(mult = c(0.03, 0.03))) +
  scale_y_continuous(
    limits = c(0:1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.01))
  ) + geom_hline(aes(yintercept = 0.4926884),
                 linetype = 5,
                 col = "grey40")



# grouped histogram
library(tidyverse)
library(ggpubr)

df <- readxl::read_xlsx("result/purity/zhongshan.xlsx", sheet = 1)

df <-
  df %>% pivot_longer(cols = c(facets:our),
                      names_to = 'method',
                      values_to = 'Purity')

ggplot(df, aes(x = sample,
               y = Purity,
               fill = method)) +
  geom_bar(
    stat = 'identity',
    width = 0.7,
    colour = 'white',
    position = 'dodge'
  ) +
  theme_pubclean() +
  scale_fill_manual(values = c("#246b93", "#F39B7F")) +
  labs(x = "Sample",
       y = "Purity",
       fill = "Method") +
  scale_y_continuous(
    limits = c(0:1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  theme(
    legend.position = "top",
    #panel.grid.major.x = element_blank(),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
