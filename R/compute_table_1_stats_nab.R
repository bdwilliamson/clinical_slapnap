# for a given set of antibodies, compute the summary statistics

# @param antibodies the list of antibodie
# @param data.assay.reduced the set of assay data
# @param data.seq the sequence data
# @param seqname.db the sequence name database
# @param subtype the subtype info
# @param country the country info
# @param year the year info
# @param measure "ic80" or "ic50" (which measurement to count)
compute_stats_nab <- function(antibodies, data.assay.reduced,
                              data.seq, seqname.db, subtype, country, year, 
                              measure = "ic80") {
  # get initial sequence, country, year info
  ab.seqs <- list()
  for(ab.index in 1:length(antibodies)) {
    ab.seqs[[ab.index]] <- unique(as.character(data.assay.reduced[data.assay.reduced[, 1] == antibodies[ab.index], ] %>% pull(Virus)))
  }
  seqs.include <- Reduce(intersect, ab.seqs)
  
  seqs.selected <- data.seq[seqname.full[seqname.db %in% seqs.include]]
  seqname.selected.full <- names(seqs.selected)
  seqname.selected.db <- seqname.db[seqname.full %in% seqname.selected.full]
  subtype.selected <- subtype[seqname.full %in% seqname.selected.full]
  country.selected <- country[seqname.full %in% seqname.selected.full]
  year.selected <- year[seqname.full %in% seqname.selected.full]
  
  # process IC50 and IC80 readouts
  # initialize our readout vectors
  imputed.ic50 <- list()
  imputed.ic80 <- list()
  censored.ic50 <- list()
  censored.ic80 <- list()
  for(ab.tmp in antibodies) {
    ab_tmp_str <- paste0("nab_", ab.tmp)
    imputed.ic50[[ab_tmp_str]] <- rep(NA, length(seqs.selected))
    imputed.ic80[[ab_tmp_str]] <- rep(NA, length(seqs.selected))
    censored.ic50[[ab_tmp_str]] <- rep(0, length(seqs.selected))
    censored.ic80[[ab_tmp_str]] <- rep(0, length(seqs.selected))
  }
  
  # collect and process our imputed/censored information
  for(ab.tmp in antibodies) {
    ab_tmp_str <- paste0("nab_", ab.tmp)
    for(seqname.index in 1:length(seqname.selected.db)) {
      
      # isolate our readouts of interest
      seqname.tmp <- seqname.selected.db[seqname.index]
      data.assay.reduced.tmp <- data.assay.reduced[data.assay.reduced$Antibody == ab.tmp & data.assay.reduced$Virus == seqname.tmp, 5:6]
      
      # IC50:  let's confirm that we actually have readout data for this sequence
      if(length(data.assay.reduced.tmp[data.assay.reduced.tmp[, 1] != "", 1]) > 0) {
        
        # make the binary notation that we have a right-censored variable
        if(">" %in% substr(data.assay.reduced.tmp[, 1], 1, 1)) {
          censored.ic50[[ab_tmp_str]][seqname.index] <- 1
        }
        
        # merge all of our readouts, both censored and numeric
        imputed.ic50[[ab_tmp_str]][seqname.index] <- merge.readouts(data.assay.reduced.tmp[, 1])
        
        # if no assay readouts exist, "Mark it zero"
      } else {
        censored.ic50[[ab_tmp_str]][seqname.index] <- NA
        imputed.ic50[[ab_tmp_str]][seqname.index] <- NA
      }
      
      # IC80:  let's confirm that we actually have readout data for this sequence
      if(length(data.assay.reduced.tmp[data.assay.reduced.tmp[, 2] != "", 2])) {
        
        # make the binary notation that we have a right-censored variable
        if(">" %in% substr(data.assay.reduced.tmp[, 2], 1, 1)) {
          censored.ic80[[ab_tmp_str]][seqname.index] <- 1
        }
        
        # merge all of our readouts, both censored and numeric
        imputed.ic80[[ab_tmp_str]][seqname.index] <- merge.readouts(data.assay.reduced.tmp[, 2])
        
        # if no assay readouts exist, "Mark it zero"
      } else {
        censored.ic80[[ab_tmp_str]][seqname.index] <- NA
        imputed.ic80[[ab_tmp_str]][seqname.index] <- NA
      }
    }
  }
  
  # combine results into new dataframe
  readouts <- data.frame(seq.id.catnap=seqname.selected.db)
  for(ab.tmp in antibodies) {
    ab_tmp_str <- paste0("nab_", ab.tmp)
    readouts.tmp <- data.frame(censored.ic50[[ab_tmp_str]], censored.ic80[[ab_tmp_str]],
                               imputed.ic50[[ab_tmp_str]], imputed.ic80[[ab_tmp_str]])
    names(readouts.tmp) <- c(paste0(ab_tmp_str, ".ic50.censored"),
                             paste0(ab_tmp_str, ".ic80.censored"),
                             paste0(ab_tmp_str, ".ic50.imputed"),
                             paste0(ab_tmp_str, ".ic80.imputed"))
    readouts <- data.frame(readouts, readouts.tmp)
  }
  
  if(length(antibodies) > 1) {
    # run additive model to compute (Bliss-Hill will be non-NA if this is)
    readouts$pc.ic50 <- apply(readouts[, grep("ic50.imputed", names(readouts), fixed=T)], 1, wagh.additive.method)
    readouts$pc.ic80 <- apply(readouts[, grep("ic80.imputed", names(readouts), fixed=T)], 1, wagh.additive.method)
  } else { # if only one antibody, define them as the single-ab version
    readouts$pc.ic50 <- readouts[, grep("ic50.imputed", names(readouts), fixed = TRUE)]
    readouts$log10.pc.ic50 <- log10(readouts$pc.ic50)
    readouts$pc.ic80 <- readouts[, grep("ic80.imputed", names(readouts), fixed = TRUE)]
    readouts$log10.pc.ic80 <- log10(readouts$pc.ic80)
  }
  readouts$log10.pc.ic50 <- log10(readouts$pc.ic50)
  readouts$log10.pc.ic80 <- log10(readouts$pc.ic80)
  
  # assemble our data set
  binary_subtype <- data.frame(bin.subtype(subtype.selected))
  data.final <- data.frame(seq.id.lanl=seqname.selected.full,
                           seq.id.catnap=seqname.selected.db,
                           readouts[, -1],
                           year=year.selected,
                           subtype_c = binary_subtype$subtype.is.C)
  # count the total number, total number since 2005, and the number of 
  # subtype C since 2005
  num_year <- as.numeric(data.final$year)
  if (measure == "ic80") {
    num_tot <- sum(!is.na(data.final$pc.ic80))
    num_since_2005 <- sum(!is.na(data.final$pc.ic80) & (num_year >= 2005 & !is.na(num_year)))
    num_c_since_2005 <- sum(!is.na(data.final$pc.ic80) & (num_year >= 2005 & !is.na(num_year)) &
                              (data.final$subtype_c == 1))
  } else {
    num_tot <- sum(!is.na(data.final$pc.ic50))
    num_since_2005 <- sum(!is.na(data.final$pc.ic50) & (num_year >= 2005 & !is.na(num_year)))
    num_c_since_2005 <- sum(!is.na(data.final$pc.ic50) & (num_year >= 2005 & !is.na(num_year)) &
                              (data.final$subtype_c == 1))
  }
  
  return(list(tot = num_tot, since_2005 = num_since_2005, c_since_2005 = num_c_since_2005))
}