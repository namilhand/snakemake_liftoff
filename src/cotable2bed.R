#library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(tibble)
library(readr)
library(stringr)
#library(tidyverse)

args <- commandArgs(trailingOnly=T)
input <- args[1]
output <- args[2]

cotable <- read_csv(input, col_names=T) %>%
	mutate(lib = str_replace(lib, "results/06_tiger/lib", "")) %>%
  	mutate(lib = str_replace(lib, "_MappedOn_tair10", "")) %>%
    mutate(lib = as.numeric(lib)) %>%
	arrange(lib, chrs, start)


#===== test data =======
#cotable <- read_csv("/datasets/data_4/nison/colcen/liftoff/data/Selected_inhouse_WT_cotable.txt", col_names=T) %>%
#	mutate(lib = str_replace(lib, "results/06_tiger/lib", "")) %>%
#  	mutate(lib = str_replace(lib, "_MappedOn_tair10", "")) %>%
#    mutate(lib = as.numeric(lib)) %>%
#	arrange(lib, chrs, start)
#
#output <- "/datasets/data_4/nison/colcen/liftoff_pipeline/cotable.gff"
#========================



cotable.indexed <- group_by(cotable, lib) %>% mutate(counter = row_number(lib)) %>% ungroup()
colnames(cotable.indexed)[2] <- "seqnames"


cotable.indexed.bed <- dplyr::select(cotable.indexed, c(seqnames, cos, lib, counter)) %>%
	mutate(start = cos) %>%
	mutate(end = cos + 3000) %>%
	mutate(attributes = paste0("ID=", lib, "_", counter)) %>%
	#mutate(width = end - start)
	dplyr::select(-c(cos, lib, counter))

cotable.bed.gr <- makeGRangesFromDataFrame(df = cotable.indexed.bed, starts.in.df.are.0based=T, keep.extra.columns=T)

chr_lengths <- c(30427671, 19698289, 23459830, 18585056, 26975502)
names(chr_lengths) <- paste0("Chr", 1:5)

seqinfo(cotable.bed.gr) <- Seqinfo(seqnames = names(chr_lengths), seqlengths = chr_lengths)

cotable.bed.gr.trim <- trim(cotable.bed.gr)
cotable.bed.trim <- as_tibble(as.data.frame(cotable.bed.gr.trim))
cotable.gff3 <- cotable.bed.trim %>%
	add_column(source = ".",
			   type = "gene",
			   score = ".",
			   phase = "."
			   ) %>%
	dplyr::select(-width) %>%
	relocate(source, .after=seqnames) %>%
	relocate(type, .after=source) %>%
	relocate(strand, .after=score) %>%
	relocate(attributes, .after=phase)

sink(output)
cat("##gff-version 3\n")
sink()
write_tsv(cotable.gff3, file=output, col_names=F, append=T)
