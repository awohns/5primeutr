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


#-----------#
# UKB DATA  #
#-----------#
for (n in c(180, 250)){

  for (i in c(1:22, 'X')){

    df = fread(str_interp("./processed/UKB_variants_${n}_annotated_chr${i}.csv"))

    df = df %>%
      mutate(AF = as.numeric(AC)/as.numeric(AN)) %>%
      mutate(alt_seq_len = nchar(ALT_sequence),
             ref_seq_len = nchar(REF_sequence)) %>%

      select(gene_name:AN, AF,
             ALT_sequence, REF_sequence,
             position_within_reporter_seq,
             start_codon, recessive,
             ref_seq_len, alt_seq_len,
             alt_seq_len_180bp_adjusted,
             alt_seq_len_250bp_adjusted,
             orf_ann_ref, orf_ann_alt,
             variant_type, gained_main_ORF:changed_pos_ORFs_total,
             variant_type_adjusted,
             predicted_category,
             posterior_hs:hs_percentile)

    if (i==1){
      ukb_df = df
    } else {
      ukb_df = rbind(ukb_df, df)

    }
  }

  write.csv(ukb_df,
            str_interp("./processed/UKB_variants_${n}_annotated.csv"),
            row.names=FALSE)

}


#---------------------#
# FILTERING SUMSTATS  #
#---------------------#
for (i in c(1:22, 'X')){

  sumstats_df = fread(str_interp("./processed/filtering_sumstats_chr${i}.csv"))
  sumstats_df$Chr = i

  if (i==1){
    main_df = sumstats_df
  } else {
    main_df = rbind(main_df, sumstats_df)
  }
}

sumstats = main_df %>%
  group_by(Filter) %>%
  summarize(Before_total = sum(Before),
            After_total = sum(After),
            Genes_total = sum(N_Genes)) %>%
  arrange(desc(Before_total)) %>%
  mutate(Variants_lost = After_total - Before_total)

write.csv(sumstats, "./processed/filtering_ukb_sumstats.csv",
          row.names=FALSE)


#-------------------#
# INVARIANT SITES   #
#-------------------#
for (n in c(180, 250)){

  for (i in c(1:22, 'X')){

    file_path =  str_interp("./processed/invariant_${n}_annotated_chr${i}.csv")

    if (file.exists(file_path)){
      df = fread(file_path)

      df = df %>%
        mutate(alt_seq_len = nchar(ALT_sequence),
               ref_seq_len = nchar(REF_sequence)) %>%

        select(gene_name, CHR, POS, REF, ALT,
               ALT_sequence, REF_sequence,
               position_within_reporter_seq,
               start_codon, recessive,
               ref_seq_len, alt_seq_len,
               orf_ann_ref, orf_ann_alt,
               alt_seq_len_180bp_adjusted,
               alt_seq_len_250bp_adjusted,
               variant_type, gained_main_ORF:changed_pos_ORFs_total,
               variant_type_adjusted,
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

  write.csv(invar_df, str_interp("./processed/invariant_${n}_annotated.csv"),
            row.names=FALSE)

}


#---------------------------------#
# FILTERING SUMSTATS - INVARIANT  #
#---------------------------------#
for (i in c(1:22, 'X')){

  file_path =  str_interp("./processed/filtering_invariant_sumstats_chr${i}.csv")

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
  arrange(desc(After_total)) %>%
  mutate(Variants_lost = After_total - Before_total)

write.csv(sumstats, "./processed/filtering_invariant_sumstats.csv",
          row.names=FALSE)
