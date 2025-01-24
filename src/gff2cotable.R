library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(stringr)

args <- commandArgs(trailingOnly=T)

input <- args[1]
output <- args[2]

gff <- read_tsv(input, skip=3, col_names=F)

#========= example
#input <- read_tsv(file="/datasets/data_4/nison/colcen/liftoff_pipeline/results/02_liftoff/Selected_inhouse_WT_Col-CC.gff", skip=3, col_names=F)
#gff <- input
#===========================

cotable <- gff %>%
	tidyr::separate_wider_delim(X9, delim=";", names=c("ID", "coverage", "sequence_ID", "extra_copy_number", "copy_num_ID"), too_many="merge") %>%
	mutate(ID=str_replace(ID, "ID=", "")) %>%
	separate_wider_delim(ID, delim="_", names=c("lib", "counter")) %>%
	dplyr::select(c(X1, X4, X5, lib)) %>%
	relocate(lib, .before=X1) %>%
	mutate(cos=X4) %>%
	mutate(start = cos-1000, stop = cos + 1000) %>%
	mutate(width = stop - start) %>%
	dplyr::select(-c(X4, X5)) %>%
	relocate(cos, .after=stop)

colnames(cotable) <- c("lib", "chrs", "start", "stop", "cos", "width")

#write_csv(cotable, file="Selected_inhouse_WT_cotable_Col-CC.liftoff.txt", col_names=T)
write_csv(cotable, file=output, col_names=T)
