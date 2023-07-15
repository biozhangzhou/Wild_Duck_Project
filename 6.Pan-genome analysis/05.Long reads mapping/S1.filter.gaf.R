##### Parameter Settings #####
### Load packages
library(optparse)
### Describes how parameters are parsed
option_list <- list(
  make_option(c("-i", "--input"), type = "character", default = F,
              action = "store", help = "path/to/input.gaf"
  ),
  make_option(c("-o", "--output"), type = "character", default = FALSE,
              action = "store", help = "path/to/output_filtered.gaf"
  ),
  make_option(c("-q", "--minmapq"), type = "integer", default = FALSE,
              action = "store", help = "Minimum MAPQ to keep (option)"
  )
)
### Parsing parameters
opt = parse_args(OptionParser(option_list = option_list, usage = "Rscript script.R <options>"))

##### Parameters  #####
INPUT <- opt$input
OUTPUT <- opt$output
MINMAPQ <- opt$minmapq

### test
#INPUT <- "data/test.gaf"
#INPUT <- "data/A.6.graph_rm_s"
#OUTPUT <- "cache/A.6.graph_map_MA_n100_filter.gaf"
#MINMAPQ <- 50

##### Running #####
### Load packages
library(tidyverse)
library(data.table)

### try to read gaf file
df<- fread(INPUT, sep="\t", fill = TRUE)

### rename col 
#Because colnames from different read functions are fifferent
colnames(df) <- paste("V", 1:ncol(df), sep = "")

### choose highest alignment score
df_filter1 <- df %>% distinct(V1, .keep_all=T)

### discard reads which <80% read length mapped to graph (ALRL)
# Residue matches/Path length
df_filter1$V10 <- as.numeric(df_filter1$V10)
df_filter1$V7 <- as.numeric(df_filter1$V7)

df_filter1$ALRL <- df_filter1$V10/df_filter1$V7
df_filter2 <- df_filter1 %>% filter(ALRL >= 0.8)



if(is.null(MINMAPQ) == T){
  df_filter3 <- df_filter2 %>% select(-ALRL)
} else
  {
  # Filter MAPQ
  df_filter3 <- df_filter2 %>% filter(V12 >= MINMAPQ) %>% select(-ALRL)
  }


write.table(df_filter3,
            OUTPUT, 
            sep="\t", 
            col.names = F,
            row.names = F,
            quote = F)
