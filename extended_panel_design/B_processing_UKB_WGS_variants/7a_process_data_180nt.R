#--------------------------------#
# DATA PROCESSING FOR < 180 NT   #
#--------------------------------#

rm(list=ls())

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)

set.seed(123456)

#setwd('~/Desktop/5prime_utr/data/')
setwd('/home/users/tami/5utr_extended_panel/data/')


#------------#
# FUNCTIONS  #
#------------#
remove_reporters_within_10p_of_length = function(df, ref_df){

  tmp_df = merge(df, ref_df %>%
                   select(gene_name,
                          ref_seq_len_180bp,
                          ref_seq_len_180bp_adjusted),
                 by='gene_name', all.x=TRUE, all.y=FALSE)

  tmp_df = tmp_df %>%

    # Calculate boundaries
    mutate(ref_lower = ref_seq_len_180bp_adjusted*0.9,
           ref_upper = ref_seq_len_180bp_adjusted*1.1,

           alt_lower = alt_seq_len_180bp_adjusted*0.9,
           alt_upper = alt_seq_len_180bp_adjusted*1.1) %>%

    # Filtering
    filter((alt_seq_len_180bp_adjusted > ref_lower &
              alt_seq_len_180bp_adjusted < ref_upper) |
             (ref_seq_len_180bp_adjusted > alt_lower &
                ref_seq_len_180bp_adjusted < alt_upper)) %>%

    # Remove columns
    select(-ref_seq_len_180bp,
           -ref_seq_len_180bp_adjusted,

           -ref_lower, -ref_upper, -alt_lower, -alt_upper)

  return(tmp_df)

}

remove_duplicates = function(main_df, df, ref_df){

  duplicates = main_df$ALT_sequence[duplicated(main_df$ALT_sequence)]

  main_df = main_df %>%

    # Remove duplicates within
    filter(!(ALT_sequence %in% duplicates)) %>%

    # Remove duplicates with invariant reporters
    filter(!(ALT_sequence %in% df$ALT_sequence)) %>%

    # Remove duplicates with references
    filter(!(ALT_sequence %in% ref_df$ALT_sequence))

  return(main_df)

}


#------------------------#
# LOAD AND PROCESS DATA  #
#------------------------#
ukb_180_df = data.frame(fread(str_interp('./processed/UKB_variants_180_annotated.csv')))
invar_180_df = data.frame(fread(str_interp('./processed/invariant_180_annotated.csv')))
ref_180_df = data.frame(fread(str_interp('./processed/references_180_annotated.csv')))


# UKB sumstats
ukb_sumstats_df = data.frame(
    fread(str_interp('./processed/filtering_ukb_sumstats.csv')))

ukb_sumstats_df[nrow(ukb_sumstats_df)+1, ] =
    c('Keep variants less than 180nt',
      ukb_sumstats_df[nrow(ukb_sumstats_df),3],
      nrow(ukb_180_df), length(unique(ukb_180_df$gene_name)),
      ukb_sumstats_df[nrow(ukb_sumstats_df),3] - nrow(ukb_180_df))


# Invariant sumstats
invar_sumstats_df = data.frame(
    fread(str_interp('./processed/filtering_invariant_sumstats.csv')))

invar_sumstats_df[nrow(invar_sumstats_df)+1, ] =
    c('Keep variant less than 180nt',
      invar_sumstats_df[nrow(invar_sumstats_df),3],
      nrow(invar_180_df),
      invar_sumstats_df[nrow(invar_sumstats_df),3]-nrow(invar_180_df))


#------------------------------------------------------------#
# QUALITY CONTROL                                            #
#   1. Keep reporters within 10% of length (REF and ALT)     #
#   2. Remove duplicates (within and between)                #
#   3. Remove references that are 24nt                       #
#------------------------------------------------------------#

# 1.
ukb_180_df = remove_reporters_within_10p_of_length(ukb_180_df, ref_180_df)
invar_180_df = remove_reporters_within_10p_of_length(invar_180_df, ref_180_df)

ukb_sumstats_df[nrow(ukb_sumstats_df)+1, ] =
    c('Remove reporters >10% length difference REF and ALT',
      ukb_sumstats_df[nrow(ukb_sumstats_df),3],
      nrow(ukb_180_df), length(unique(ukb_180_df$gene_name)),
      nrow(ukb_180_df) - as.numeric(ukb_sumstats_df[nrow(ukb_sumstats_df),3]))

invar_sumstats_df[nrow(invar_sumstats_df)+1, ] =
    c('Remove reporters >10% length difference REF and ALT',
      invar_sumstats_df[nrow(invar_sumstats_df),3],
      nrow(invar_180_df),
      nrow(invar_180_df) - as.numeric(invar_sumstats_df[nrow(invar_sumstats_df),3]))


# 2.
ukb_180_df = remove_duplicates(ukb_180_df, invar_180_df, ref_180_df)
invar_180_df = remove_duplicates(invar_180_df, ukb_180_df, ref_180_df)

ukb_sumstats_df[nrow(ukb_sumstats_df)+1, ] =
    c('Remove duplicates',
      ukb_sumstats_df[nrow(ukb_sumstats_df),3],
      nrow(ukb_180_df), length(unique(ukb_180_df$gene_name)),
      nrow(ukb_180_df) - as.numeric(ukb_sumstats_df[nrow(ukb_sumstats_df),3]))

invar_sumstats_df[nrow(invar_sumstats_df)+1, ] =
    c('Remove duplicates',
      invar_sumstats_df[nrow(invar_sumstats_df),3],
      nrow(invar_180_df),
      nrow(invar_180_df) - as.numeric(invar_sumstats_df[nrow(invar_sumstats_df),3]))


# 3. 
ukb_180_df = ukb_180_df %>% 
  filter(ref_seq_len >= 25)

invar_180_df = invar_180_df %>%
  filter(ref_seq_len >= 25)

ukb_sumstats_df[nrow(ukb_sumstats_df)+1, ] =
  c('Remove variants with 24nt references',
    ukb_sumstats_df[nrow(ukb_sumstats_df),3],
    nrow(ukb_180_df), length(unique(ukb_180_df$gene_name)),
    nrow(ukb_180_df) - as.numeric(ukb_sumstats_df[nrow(ukb_sumstats_df),3]))

invar_sumstats_df[nrow(invar_sumstats_df)+1, ] =
  c('Remove variants with 24nt references',
    invar_sumstats_df[nrow(invar_sumstats_df),3],
    nrow(invar_180_df),
    nrow(invar_180_df) - as.numeric(invar_sumstats_df[nrow(invar_sumstats_df),3]))


#-----------------------------------------------#
# SANITY CHECKS                                 #
#   1. No duplicates                            #
#   2. Duplicates within references             #
#   3. Lengths of reporters should be 25-180    #
#-----------------------------------------------#

# 1.
reporters = c(ukb_180_df$ALT_sequence, invar_180_df$ALT_sequence)
n_duplicates = length(reporters[duplicated(reporters)])
print(str_interp("There are ${n_duplicates} duplicates."))

# 2.
ref_reporters = ref_180_df$ALT_sequence
n_duplicates = length(ref_reporters[duplicated(ref_reporters)])
print(str_interp("There are ${n_duplicates} duplicates."))

# 3.
range((ukb_180_df %>%
         mutate(reporter_len = nchar(ALT_sequence) -
                  str_count(ALT_sequence, fixed("*"))))$reporter_len)
range((invar_180_df %>%
         mutate(reporter_len = nchar(ALT_sequence) -
                  str_count(ALT_sequence, fixed("*"))))$reporter_len)


#------------#
# SAVE DATA  #
#------------#
write.csv(ukb_180_df, './processed/ukb_180_annotated_processed.csv', row.names=FALSE)
write.csv(invar_180_df, './processed/invar_180_annotated_processed.csv', row.names=FALSE)
write.csv(ref_180_df, './processed/ref_180_annotated_processed.csv', row.names=FALSE)

write.csv(ukb_sumstats_df, './processed/filtering_ukb_180_sumstats.csv', row.names=FALSE)
write.csv(invar_sumstats_df, './processed/filtering_invariant_180_sumstats.csv', row.names=FALSE)
