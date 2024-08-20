rm(list=ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)


#------------#
# VARIABLES  #
#------------#
bp_limits = c(250)

stop_codons = c('TAA', 'TAG', 'TGA')


#------------#
# LOAD DATA  #
#------------#
gencode_df = data.frame(fread('../../../data/processed/gencode_v45_protein_coding_MANE_clean.csv'))
gencode2_df = data.frame(fread('../../../data//processed/gencode2_v45_protein_coding_MANE_clean.csv'))

gencode_transcripts_df =
  data.table(fread('../../../data/processed/gencode_v45_MANE_transcripts_5utr_sequences.csv'))


#------------------#
# EXTRACT 5' UTRs  #
#------------------#
df_five_prime_UTR = gencode2_df %>%
  filter(Element == 'five_prime_UTR') %>%
  select(-Element)

# Merge gencode data and sequences
df_five_prime_UTR = merge(df_five_prime_UTR, gencode_transcripts_df,
                          by='transcript_name')

# Construct sequences with length limits
for (bp_limit in bp_limits){

  limit = bp_limit-6

  colname = str_interp("UTR5_sequence_${limit}")

  length_w_limit = ifelse(df_five_prime_UTR$Length > limit,
                          limit, df_five_prime_UTR$Length)

  df_five_prime_UTR = df_five_prime_UTR %>%
    mutate(length_w_limit = ifelse(df_five_prime_UTR$Length > limit,
                                   limit, df_five_prime_UTR$Length)) %>%
    mutate(start = Length - length_w_limit + 1)

  colname = str_interp("UTR5_sequence_${bp_limit}")
  df_five_prime_UTR[,colname] = substr(df_five_prime_UTR$UTR5_sequence_and_6bp,
                                       df_five_prime_UTR$start,
                                       df_five_prime_UTR$Length+6)

  df_five_prime_UTR = df_five_prime_UTR %>%
    select(-length_w_limit, -start)

}


#--------------------#
# NOTE START CODONS  #
#--------------------#
df_start_codons = gencode_transcripts_df %>%
  mutate(length = nchar(UTR5_sequence_and_6bp),
         start = length-5,
         end = length-3) %>%
  mutate(start_codon = substr(UTR5_sequence_and_6bp, start, end)) %>%
  select(transcript_name, start_codon)

write.csv(df_start_codons,
          '../../../data/processed/start_codons_gencode_v45.csv', row.names=FALSE)

# Merge start codons
df_five_prime_UTR = merge(df_five_prime_UTR, df_start_codons, by='transcript_name')


#------------#
# SAVE DATA  #
#------------#
write.csv(df_five_prime_UTR,
          '../../../data/processed/5utr_ref_sequences_250bp_gencode_v45.csv', row.names=FALSE)
