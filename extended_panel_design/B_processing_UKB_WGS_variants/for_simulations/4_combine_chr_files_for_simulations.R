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
process_filtering_sumstats = function(filtering_sumstats){

  df = data.frame(matrix(nrow=3, ncol=5))
  names(df) = c('Filter', 'Before', 'After', 'N_vars_lost', 'Perc_var_lost')

  df$Filter = c('Keep only single nucleotide variants',
                'Remove variants overlapping >1 gene',
                'Within the 250bp sequence')

  for (i in 1:nrow(df)){

    df[i, 2:3] = c(sum(filtering_sumstats[,i+1]),
                   sum(filtering_sumstats[,i+2]))

  }

  df$N_vars_lost = df$Before-df$After
  df$Perc_var_lost = round((df$N_vars_lost / df$Before)*100, 0)

  return(df)
}


#--------------------#
# LOAD GENCODE DATA  #
#--------------------#
gencode_exons_df = data.frame(fread('./processed/gencode_v45_protein_coding_MANE_clean.csv'))

gencode_6bp_df = data.frame(fread('./processed/gencode_v45_6bp_coordinates.csv'))
gencode_6bp_df = gencode_6bp_df %>%
  select(-exon_length) %>%
  filter(gene_name %in% gencode_exons_df$gene_name)

gencode_df = rbind(gencode_exons_df, gencode_6bp_df)

gencode_sum_df = gencode_df %>%
  filter(Element == 'five_prime_UTR') %>%
  mutate(exon_length = abs(End-Start)+1) %>%
  group_by(gene_name) %>%
  mutate(utr_length = sum(exon_length)) %>%
  select(gene_name, Chr, utr_length)

gencode_genes = unique(gencode_sum_df$gene_name)


#-------------------#
# INVARIANT SITES   #
#-------------------#
for (i in c(1:22, 'X')){

  file_path =  str_interp("./processed/UKB_invariant_variants_annotated_preliminary_panel_chr${i}.csv")

  if (file.exists(file_path)){
    df = fread(file_path)

    df = df %>%
      mutate(seq_len = nchar(UTR5_sequence_250)) %>%
      select(gene_name, CHR, POS, REF, ALT, seq_len,
             position_within_reporter_seq,
             start_codon, UTR5_sequence_250,
             orf_ann_250bp, orf_ann_250bp_ALT,
             kozak_strength, kozak_strength_ALT,
             variant_type, gained_main_ORF:kozak_stronger,
             predicted_category,
             posterior_hs:hs_percentile)

    if (i==1){
      invar_df = df
    } else {
      invar_df = rbind(invar_df, df)
    }

  } else {
    print(paste("File does not exist for i =", i))
  }

}

write.csv(invar_df, "./processed/UKB_invariant_variants_annotated_preliminary_panel.csv",
          row.names=FALSE)


#---------------------------------#
# FILTERING SUMSTATS - INVARIANT  #
#---------------------------------#
for (i in c(1:22, 'X')){

  file_path =  str_interp("./processed/filtering_invariant_sumstats_preliminary_panel_chr${i}.csv")

  if (file.exists(file_path)) {

    df = fread(file_path)
    df$Chr = i

    if (i==1){
      sumstats_df = df
    } else {
      sumstats_df = rbind(sumstats_df, df)
    }

  } else {
    print(paste("File does not exist for i =", i))
  }

}

sumstats = sumstats_df %>%
  group_by(Filter) %>%
  summarize(Before_total = sum(Before),
            After_total = sum(After)) %>%
  arrange(desc(Before_total)) %>%
  mutate(Variants_lost = After_total - Before_total)

write.csv(sumstats, "./processed/filtering_invariant_sumstats_preliminary_panel.csv",
          row.names=FALSE)
