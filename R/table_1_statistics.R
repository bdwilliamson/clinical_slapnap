# download CATNAP data, get statistics for table 1:
# (a) number of total viruses, (b) number since 2005, (c) number subtype C since 2005

# set up libraries and list of nabs for CATNAP queries
library("here")
library("tidyverse")
library("seqinr")

source(here("R", "compute_table_1_stats_nab.R"))
# these next two come from SLAPNAP
source(here("R", "00_utils.R"))
source(here("R", "02_multi_ab.Rlib"))

nabs <- c("VRC01", "VRC07-523-LS", "PGT121", "PGDM1400", "VRC26.25", 
          "VRC07-523-LS;PGT121", "VRC07-523-LS;PGDM1400", "VRC07-523-LS;VRC26.25",
          "VRC07-523-LS;10-1074", "VRC07-523-LS;PGDM1400;PGT121", "VRC01/PGDM1400-10E8v4")

# read in the assay, virus, mAb data
data.assay <- read_delim(here("dat", "catnap", "assay.txt"), delim = "\t", quote = "\"")
data.viruses <- read_delim(here("dat", "catnap", "viruses.txt"), delim = "\t", quote = "\"")
data.abs <- read_delim(here("dat", "catnap", "abs.txt"), delim = "\t", quote = "\"")

# remove outliers, same as in SLAPNAP
data.assay.reduced <- data.assay[data.assay$IC50 != ">1", ]

# load and process virus info and sequences
data.seq <- read.fasta(here("dat", "catnap", "virseqs_aa.fasta"), seqtype="AA")
seqname.full <- names(data.seq)
header.info <- strsplit(names(data.seq), split=".", fixed = TRUE)
subtype <- unlist(lapply(header.info, function(x) return(x[1])))
country <- unlist(lapply(header.info, function(x) return(x[2])))
year <- unlist(lapply(header.info, function(x) return(x[3])))
seqname.db <- unlist(lapply(header.info, function(x) return(x[4])))

table_1_statistics <- NULL
for (i in seq_len(length(nabs))) {
    antibody_string <- nabs[i]
    antibodies <- strsplit(antibody_string, split = ";")[[1]]
  
    these_summary_statistics <- suppressWarnings(
      compute_stats_nab(antibodies = antibodies,
                        data.assay.reduced = data.assay.reduced,
                        data.seq = data.seq, seqname.db = seqname.db,
                        subtype = subtype, country = country, 
                        year = year, measure = "ic80")
    )
    table_1_statistics <- dplyr::bind_rows(
      table_1_statistics, 
      tibble::tibble(nab = antibody_string, num_total = these_summary_statistics$tot,
                     num_since_2005 = these_summary_statistics$since_2005,
                     num_c_since_2005 = these_summary_statistics$c_since_2005)
    )
    
}
if (!dir.exists(here("R_output"))) {
  dir.create(here("R_output"))
}
readr::write_csv(table_1_statistics, here("R_output", "table_1_summary_statistics.csv"))
