# SNV ---------------------------------------------------------------------
all.maf <-
  list.files(
    path = "data/TARGET-OS/maf/",
    pattern = ".gz",
    full.names = T,
    recursive = T
  )

all.maf[1:3]

maf.list <- lapply(
  all.maf,
  data.table::fread,
  sep = "\t",
  header = T,
  skip = 7
)

maf.merge <- do.call(rbind, maf.list)

dim(maf.merge)
write.csv(maf.merge, "maf.csv", row.names = F, na = "")


# Clinical -----------------------------------------------------------------
library("rjson")
json <-
  jsonlite::fromJSON("data/TARGET-OS/metadata.cart.2024-08-09.json")

entity_submitter_id <-
  sapply(json$associated_entities, function(x) {
    x[, 1]
  })
case_id <- sapply(json$associated_entities, function(x) {
  x[, 3]
})
sample_case <- t(rbind(entity_submitter_id, case_id))

clinical <-
  read.delim('data/TARGET-OS/clinical.cart.2024-08-09/clinical.tsv',
             header = T)
clinical <- as.data.frame(clinical[!duplicated(clinical$case_id), ])

clinical_matrix <-
  merge(sample_case, clinical, by.x = "V3", by.y = "case_id")
clinical_matrix <- clinical_matrix[, -1]


# snv match clinical--------------------------------------------------------------------
rearrange_rows <- function(df) {
  for (i in 1:nrow(df)) {
    if (grepl("-01A-", df$V1[i])) {
      if (grepl("-10[AB]-", df$V2[i])) {
        next
      } else {
        temp <- df$V1[i]
        df$V1[i] <- df$V2[i]
        df$V2[i] <- temp
      }
    } else if (grepl("-10[AB]-", df$V1[i])) {
      if (grepl("-01A-", df$V2[i])) {
        temp <- df$V1[i]
        df$V1[i] <- df$V2[i]
        df$V2[i] <- temp
      }
    }
  }
  return(df)
}

clinical_matrix <- rearrange_rows(clinical_matrix)
clinical_matrix <- clinical_matrix[, -c(2, 3)]
clinical_matrix <-
  clinical_matrix[, sapply(clinical_matrix, function(col)
    any(col != "'--"))]
colnames(clinical_matrix)[1] <- "TumorSampleID"
clinical_matrix[clinical_matrix == "'--"] <- NA

write.csv(clinical_matrix, "clinical.csv", row.names = F)

a <-  read.csv("clinical.csv")
b <- read.csv("maf.csv")
a <- a[a$TumorSampleID %in% unique(b$Tumor_Sample_Barcode),]
write.csv(a, "clinical.csv", row.names = F, na = "")
write.csv(b, "maf.csv", row.names = F, na = "")



# CNV ---------------------------------------------------------------------
all.cnv <-
  list.files(
    path = "data/TARGET-OS/cnv/",
    pattern = ".txt",
    full.names = T,
    recursive = T
  )

all.cnv[1:3]

cnv.list <- lapply(all.cnv,
                   data.table::fread,
                   sep = "\t",
                   header = T)

cnv.list <- cnv.list[-c(10, 83)]
cnv.merge <- do.call(rbind, cnv.list)

dim(cnv.merge)
write.csv(cnv.merge, "cnv.csv", row.names = F, na = "")


# merge SNV,CNV,clinical ----------------------------------------------------------
snv <- read.csv("maf.csv")
cnv <- read.csv("cnv.csv")
clinal <- read.csv("clinical.csv")

snv <-
  snv[snv$Tumor_Sample_UUID %in% cnv$GDC_Aliquot,]

cnv <-
  cnv[cnv$GDC_Aliquot %in% snv$Tumor_Sample_UUID,]

clinal <-
  clinal[clinal$TumorSampleID %in% snv$Tumor_Sample_Barcode,]

write.csv(snv, "data/wes/maf.csv", row.names = F, na = "")
write.csv(cnv, "data/wes/cnv.csv", row.names = F, na = "")
write.csv(clinal,
          "data/wes/clinical.csv",
          row.names = F,
          na = "")

snv <- read.csv("data/wes/maf.csv")
cnv <- read.csv("data/wes/cnv.csv")

barcode <-
  data.frame(
    Tumor_Sample_Barcode = snv$Tumor_Sample_Barcode,
    Tumor_Sample_UUID = snv$Tumor_Sample_UUID
  )

barcode <- unique(barcode)

cnv <- merge(cnv,
             barcode,
             by.x = "GDC_Aliquot",
             by.y = "Tumor_Sample_UUID",
             all.x = T)

cnv <- cnv[, -1]
cnv <- cnv[, c(7, 1:6)]

write.csv(cnv, "data/cnv.csv", row.names = F, na = "")


# clonal format---------------------------------------------------------
snv <- read.csv("data/wes/maf.csv")

snv_list <- split(snv, snv$Tumor_Sample_Barcode)

snv_list <- lapply(snv_list, function(x) {
  x[, c(1, 5, 6, 41, 42)]
})

cnv <- read.csv("data/wes/cnv.csv")

cnv_list <- split(cnv, cnv$Tumor_Sample_Barcode)

library(dplyr)

# return copy number
get_copy_number <- function(chromosome, start_position, cnv_df) {
  matched_cnv <- cnv_df %>%
    filter(Chromosome == chromosome,
           Start <= start_position,
           End >= start_position)
  
  if (nrow(matched_cnv) > 0) {
    list(
      Major_Copy_Number = matched_cnv$Major_Copy_Number[1],
      Minor_Copy_Number = matched_cnv$Minor_Copy_Number[1]
    )
  } else {
    list(Major_Copy_Number = 1,
         Minor_Copy_Number = 1)
  }
}

# get copy number
process_snv_cnv <- function(snv_df, cnv_df) {
  snv_df %>%
    rowwise() %>%
    mutate(copy_info = list(get_copy_number(Chromosome, Start_Position, cnv_df))) %>%
    mutate(
      Major_Copy_Number = copy_info$Major_Copy_Number,
      Minor_Copy_Number = copy_info$Minor_Copy_Number
    ) %>%
    select(-copy_info) %>%
    ungroup()
}

if (length(snv_list) != length(cnv_list)) {
  stop("Error,length of snv_list is not match cnv_list")
}

result_list <-
  mapply(process_snv_cnv, snv_list, cnv_list, SIMPLIFY = FALSE)

for (i in seq_along(result_list)) {
  result_list[[i]] <-
    as.data.frame(result_list[[i]][, c(2, 3, 1, 4:7)])
  colnames(result_list[[i]]) <-
    c("chr", "pos", "gene", "ref", "alt", "major", "minor")
}

# remove major and minor == 0
for (i in seq_along(result_list)) {
  result_list[[i]] <- result_list[[i]][result_list[[i]]$major != 0, ]
}

# save
inputList <- result_list
save(inputList, file = "result/purity/inputList.Rdata")


# bar plot
barplot(
  height = c(19, 42),
  names.arg = c('Metastasis', 'Non-Metastasis'),
  ylim = c(0, 50)
)
