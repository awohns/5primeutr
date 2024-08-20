rm(list=ls())

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

setwd('/home/users/tami/5utr_extended_panel/data/')


#------------#
# FUNCTIONS  #
#------------#
process_ref_data = function(raw_df){
  
  df = raw_df %>%
    
    # Length refers to the length of 5' UTR
    rename(Length_UTR = 'Length',
           CHR = 'Chr',
           REF_sequence_250bp = 'UTR5_sequence_250') %>%
    
    # Sequence lengths (250bp)
    mutate(ref_seq_len_250bp = nchar(REF_sequence_250bp)) %>%
    
    # Shorten sequences longer than 180nt 
    mutate(REF_sequence_180bp = ifelse(ref_seq_len_250bp > 180,
                                       substr(REF_sequence_250bp,
                                              ref_seq_len_250bp-179, ref_seq_len_250bp),
                                       REF_sequence_250bp)) %>%
    
    # Sequence lengths (180bp)
    mutate(ref_seq_len_180bp = nchar(REF_sequence_180bp)) %>%
    
    # Sequence lengths (adjusted for stars) %>%
    mutate(ref_seq_len_180bp_adjusted = 
             ref_seq_len_180bp - str_count(REF_sequence_180bp, "\\*"),
           
           ref_seq_len_250bp_adjusted = 
             ref_seq_len_250bp - str_count(REF_sequence_250bp, "\\*"))
  
  return(df)
  
}


#------------------------#
# LOAD AND PROCESS DATA  #
#------------------------#
raw_ref_df = data.frame(
  fread('./processed/5utr_ref_sequences_250bp_gencode_v45.csv'))

ref_df = process_ref_data(raw_ref_df)

ref_180_df = ref_df %>%
  mutate(ALT_sequence = REF_sequence_180bp)

ref_250_df = ref_df %>%
  mutate(ALT_sequence = REF_sequence_250bp)


#------------#
# SAVE DATA  #
#------------#
write.csv(ref_180_df, 
          str_interp("./processed/references_180_annotated.csv"),
          row.names=FALSE)

write.csv(ref_250_df, 
          str_interp("./processed/references_250_annotated.csv"),
          row.names=FALSE)

