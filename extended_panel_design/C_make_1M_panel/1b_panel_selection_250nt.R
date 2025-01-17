#----------------------------------#
# PANEL SELECTION FOR 180-250 NT   #
#----------------------------------#

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
# VARIABLES  #
#------------#
n = 250000
n_variants = n*0.9
n_invar = n*0.1


#------------#
# FUNCTIONS  #
#------------#
keep_one_per_duplicate = function(ref_5utrs_df){

  duplicate_refs = ref_5utrs_df$ALT_sequence[duplicated(ref_5utrs_df$ALT_sequence)]

  refs_duplicates_df = ref_5utrs_df %>%
    filter(ALT_sequence %in% duplicate_refs) %>%
    group_by(ALT_sequence) %>%
    mutate(alternative_annotations = sapply(gene_name, function(g) {
      other_genes <- setdiff(gene_name, g)
      paste(other_genes, collapse = ', ')})) %>%
    ungroup()

  ref_not_duplicated_df = ref_5utrs_df %>%
    filter(!(ALT_sequence %in% duplicate_refs)) %>%
    mutate(alternative_annotations = NA)

  for (ref in duplicate_refs){

    tmp_df = refs_duplicates_df %>%
      filter(ALT_sequence == ref)

    keep_tmp_df = tmp_df[1,]

    if (ref == duplicate_refs[1]){
      refs_keep_one_duplicate_df = keep_tmp_df
    } else {
      refs_keep_one_duplicate_df = rbind(refs_keep_one_duplicate_df, keep_tmp_df)
    }
  }

  ref_df = rbind(ref_not_duplicated_df, refs_keep_one_duplicate_df)

  return(ref_df)

}

add_NA_missing_columns = function(ref_5utrs_panel_df, selected_vars_df){

  missing_columns = setdiff(names(selected_vars_df), names(ref_5utrs_panel_df))

  for (col_name in missing_columns) {
    ref_5utrs_panel_df[[col_name]] = NA
  }

  return(ref_5utrs_panel_df)
}


#-------------#
# LOAD DATA   #
#-------------#
ukb_df = data.frame(fread('./processed/ukb_250_annotated_processed.csv'))
invar_df = data.frame(fread('./processed/invar_250_annotated_processed.csv'))
ref_df = data.frame(fread('./processed/ref_250_annotated_processed.csv'))

# Process data
ukb_df = ukb_df %>%
  rename(seq_len = 'alt_seq_len_250bp_adjusted') %>%
  select(-alt_seq_len_180bp_adjusted)

invar_df = invar_df %>%
  rename(seq_len = 'alt_seq_len_250bp_adjusted') %>%
  select(-alt_seq_len_180bp_adjusted)

ref_df = ref_df %>%
  rename(seq_len = 'ref_seq_len_250bp_adjusted') %>%
  select(gene_name, CHR, start_codon, ALT_sequence, seq_len)


####################################################################################

#-------------------#
# PANEL SELECTION   #
#-------------------#
n_vars_to_select = n_variants

selected_vars_df = data.frame(matrix(nrow=0, ncol=ncol(ukb_df)))
names(selected_vars_df) = names(ukb_df)

summary_df = data.frame(matrix(nrow=0, ncol=3))
names(summary_df) = c('Category', 'N_variants', 'N_genes')


#-------------------------------------------------------------------#
# 1. KEEP ALL OBSERVED VARIANTS IN HIGH S-HET AND RECESSIVE GENES   #
#-------------------------------------------------------------------#

# Split by high and low shet genes
df_low_hs = ukb_df %>%
  filter(hs_percentile <= 65) %>%
  filter(recessive == 'no')

df_high_hs = ukb_df %>%
  filter(hs_percentile > 65 | recessive == 'yes')

mean_utr_len_low_hs = round(mean((df_low_hs %>%
                                    select(gene_name, seq_len) %>%
                                    distinct())$seq_len), 1)
mean_utr_len_high_hs = round(mean((df_high_hs %>%
                                     select(gene_name, seq_len) %>%
                                     distinct())$seq_len), 1)

print(str_interp("Average 5' UTR length in low shet genes is ${mean_utr_len_low_hs}."))
print(str_interp("Average 5' UTR length in high shet and recessive genes is ${mean_utr_len_high_hs}."))

n_low_hs = nrow(df_low_hs)
n_high_hs = nrow(df_high_hs)

n_genes_low_hs = length(unique(df_low_hs$gene_name))
n_genes_high_hs = length(unique(df_high_hs$gene_name))

print(str_interp("There are ${n_low_hs} variants in low shet genes, and ${n_genes_low_hs} genes."))
print(str_interp("There are ${n_high_hs} variants in high shet genes, and ${n_genes_high_hs} genes."))


##################################################
# KEEP all observed variants in high hs genes    #

n_vars_to_select = n_vars_to_select-n_high_hs
selected_vars_df = rbind(selected_vars_df, df_high_hs)

summary_df[nrow(summary_df)+1, ] =
  c('Observed variants in high s-het and recessive genes', n_high_hs, n_genes_high_hs)
                                                 #
##################################################


#--------------------------------------------#
# 2. KEEP ALL MAF > 0.1% IN LOW S-HET GENES  #
#--------------------------------------------#
df_low_hs_maf_1p = df_low_hs %>%
  filter(AF >= 0.001)

n_vars = nrow(df_low_hs_maf_1p)
print(str_interp("There are ${n_vars} variants in low s-het genes, MAF >= 0.1%."))


##################################################
# KEEP all MAF > 1% variants in low shet genes   #

n_vars_to_select = n_vars_to_select-n_vars
selected_vars_df = rbind(selected_vars_df, df_low_hs_maf_1p)

summary_df[nrow(summary_df)+1, ] =
  c('MAF > 0.1% in low s-het genes', n_vars, NA) #
##################################################


#-----------------------------------------------------------------#
# 3. RANDOM SAMPLING MAF < 0.1% IN LOW S-HET GENES UP TO 900,000  #
#-----------------------------------------------------------------#
set.seed(123456)

df_sampled = df_low_hs %>%
  filter(AF < 0.001) %>%
  sample_n(n_vars_to_select, replace=FALSE)

n_vars = nrow(df_sampled)
print(str_interp("We sampled ${n_vars} variants in low s-het genes, MAF < 0.1%."))


###################################################
# SAMPLE REST FROM MAF < 0.1%  in low shet genes  #

n_vars_to_select = n_vars_to_select-n_vars
selected_vars_df = rbind(selected_vars_df, df_sampled)

summary_df[nrow(summary_df)+1, ] =
  c('MAF < 0.1% in low s-het genes', n_vars, NA)  #
###################################################


#----------------------------#
# 4. ADD REFERENCE SEQUENCES #
#----------------------------#
panel_genes = unique(selected_vars_df$gene_name)
n_panel_genes = length(panel_genes)

ref_5utrs_panel_df = ref_df %>%
  filter(gene_name %in% panel_genes) %>%
  mutate(gene_name = paste(gene_name, "_ref", sep=''))

ref_5utrs_panel_df = keep_one_per_duplicate(ref_5utrs_panel_df)

n_reporters = nrow(ref_5utrs_panel_df)
print(str_interp("There are ${n_panel_genes} genes and ${n_reporters} unique reporters."))

# Identify missing columns in ref_5utrs_panel_df and add NA
ref_5utrs_panel_df = add_NA_missing_columns(ref_5utrs_panel_df, selected_vars_df)

# Add information to summary_df
summary_df[nrow(summary_df)+1, ] =
  c('Reference genes', n_reporters, n_panel_genes)


#--------------------------------#
# 5. ADDING UNOBSERVED VARIANTS  #
#--------------------------------#
n_invar_sites = n_invar - nrow(ref_5utrs_panel_df)

set.seed(123456)
selected_invar_null_df = invar_df %>%
  filter(predicted_category == "Null") %>%
  sample_n(0.1*n_invar_sites)

set.seed(123456)
selected_invar_pos_df = invar_df %>%
  filter(predicted_category == "Positive") %>%
  sample_n(0.5*n_invar_sites+2)

set.seed(123456)
selected_invar_neg_df = invar_df %>%
  filter(predicted_category %in% c("Negative", 'Total loss')) %>%
  sample_n(0.4*n_invar_sites)

selected_invar_df = rbind(selected_invar_null_df,
                          selected_invar_pos_df,
                          selected_invar_neg_df)

selected_invar_df = add_NA_missing_columns(selected_invar_df, selected_vars_df)

n_null = nrow(selected_invar_null_df)
n_pos = nrow(selected_invar_pos_df)
n_neg = nrow(selected_invar_neg_df)

summary_df[nrow(summary_df)+1, ] = c('Positive variants', n_pos, NA)
summary_df[nrow(summary_df)+1, ] = c('Null variants', n_null, NA)
summary_df[nrow(summary_df)+1, ] = c('Negative variants', n_neg, NA)

################################################################################################

#----------------------#
# RESULTING 250K PANEL #
#----------------------#
final_panel_df = rbind(selected_vars_df %>%
                         mutate(alternative_annotations = NA),
                       selected_invar_df %>%
                         mutate(alternative_annotations = NA),
                       ref_5utrs_panel_df)

final_panel_df = final_panel_df %>%
  rename(insert_seq = 'ALT_sequence')

n_variants = nrow(final_panel_df)
n_unique_reporters = length(unique(final_panel_df$insert_seq))

n_variants == n_unique_reporters

summary_df[nrow(summary_df)+1, ] = c('Final panel', n_variants, n_panel_genes)

print(str_interp("The total number of variants is ${n_variants} and of reporters is ${n_unique_reporters}."))


#------------------#
# SAVE 250K PANEL  #
#------------------#
write.csv(final_panel_df, './processed/panel_250k.csv', row.names=FALSE)
write.csv(summary_df, './processed/panel_250k_sumstats.csv', row.names=FALSE)
