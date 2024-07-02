rm(list=ls())

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)

set.seed(123456)

setwd('~/Desktop/5prime_utr/data/')


#------------#
# FUNCTIONS  #
#------------#
remove_variants_not_within_10p_of_ref_alt_lengths = function(df, sumstats_df, alternative_seq, data_name){
  
  n_before = nrow(df)
  
  df$alt_seq_len = nchar(df[,alternative_seq])
  
  # Remove variants where REF and ALT sequence are not within a 10% size from each other 
  df$min_len = df$seq_len*0.9
  df$max_len = df$seq_len*1.1
  
  df$remove = ifelse(df$alt_seq_len < df$min_len | df$alt_seq_len > df$max_len, 'yes', 'no')
  
  df = df %>%
    filter(remove == 'no') %>%
    select(-min_len, -max_len, -remove)
  
  n_after = nrow(df)
  n_genes = length(unique(df$gene_name))

  if (data_name == "invariant"){
    sumstats_df[nrow(sumstats_df)+1, ] = c('Remove variants where REF/ALT seq are not within a 10% size difference',
                                           n_before, n_after, n_after-n_before)
  } else {
    sumstats_df[nrow(sumstats_df)+1, ] = c('Remove variants where REF/ALT seq are not within a 10% size difference',
                                           n_before, n_after, n_genes, n_after-n_before)
  }
  
  return(list(df, sumstats_df))
  
}

annotate_vartype = function(final_panel_df){
  
  final_panel_df$vartype1 = ifelse((nchar(final_panel_df$REF) == 1 & nchar(final_panel_df$ALT) == 1),
                                   'SNP', 'Larger variant')
  
  final_panel_df$vartype2 = ifelse(final_panel_df$vartype1 == 'SNP' & 
                                     final_panel_df$ALT == "*", 'Deletion', '')
  
  final_panel_df$vartype = paste(final_panel_df$vartype2, ", ", 
                                 final_panel_df$vartype1, sep='')
  
  # Replace NA, NA with reference 
  final_panel_df$vartype = gsub("NA, NA", "Reference", final_panel_df$vartype)
  final_panel_df$vartype = gsub(", Larger variant", "Larger variant", final_panel_df$vartype)
  final_panel_df$vartype = gsub("Deletion, SNP", "SNP, indel", final_panel_df$vartype)
  final_panel_df$vartype = gsub(", SNP", "SNP", final_panel_df$vartype)
  
  final_panel_df = final_panel_df %>%
    select(-vartype1, -vartype2)
  
  return(final_panel_df)
}

add_tami_theme = function(p){
  
  p = p +
    theme_light() +
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=16),
          plot.title = element_text(hjust=0.5, face='bold', size=18),
          strip.text = element_text(size=14),
          strip.background = element_rect(fill = "black", color = "black"),
          legend.position = 'none')
  
  print(p)
  
}

make_freq_df = function(df, colname, order=NA){
  
  freq_df = data.frame(table(df[,colname]))
  names(freq_df) = c('Variable', 'Count')
  
  freq_df$Percentage = round(freq_df$Count/sum(freq_df$Count), 3)*100
  
  if (length(order) != 1){
    
    freq_df$Variable = factor(freq_df$Variable, levels=order)
    
    freq_df = freq_df[order(freq_df$Variable), ]
  }
  
  return(freq_df)
  
}

remove_duplicate_reporters = function(df, sumstats_df, data_txt){
  
  n_variants_before = nrow(df)
  n_genes_before = length(unique(df$gene_name))
  
  duplicates = df$ALT_sequence_250bp[duplicated(df$ALT_sequence_250bp)]
  
  df = df %>% 
    filter(!(ALT_sequence_250bp %in% duplicates))
  
  n_variants_after = nrow(df)
  n_genes_after = length(unique(df$gene_name))
  
  if (data_txt == 'ukb'){
    sumstats_df[nrow(sumstats_df)+1, ] = 
      c('Remove duplicate reporters', n_variants_before, n_variants_after,
        n_genes_after, n_variants_after-n_variants_before)
    
  } else if (data_txt == 'invariant'){
    
    sumstats_df[nrow(sumstats_df)+1, ] = 
      c('Remove duplicate reporters', n_variants_before, n_variants_after,
        n_variants_after-n_variants_before)
  }
  
  return(list(df, sumstats_df))
  
}

remove_multinucleotide_variants_outside_25_250nt = function(df, sumstats_df, data_txt){

  n_variants_before = nrow(df)
  n_genes_before = length(unique(df$gene_name))
  
  # Compute reporter length and filter out variants with reporters >250nt 
  df = df %>%
    mutate(reporter_len = nchar(ALT_sequence_250bp) - str_count(ALT_sequence_250bp, fixed("*"))) %>%
    filter(reporter_len <= 250,
           reporter_len >= 25)
  
  n_variants_after = nrow(df)
  n_genes_after = length(unique(df$gene_name))
  
  if (data_txt == 'ukb'){
    sumstats_df[nrow(sumstats_df)+1, ] = 
      c('Remove variants that produce reporters >250nt and <25nt', 
        n_variants_before, n_variants_after,
        n_genes_after, n_variants_after-n_variants_before)
    
  } else if (data_txt == 'invariant'){
    
    sumstats_df[nrow(sumstats_df)+1, ] = 
      c('Remove variants that produce reporters >250nt and <25nt', 
        n_variants_before, n_variants_after,
        n_variants_after-n_variants_before)
    
  }
  
  return(list(df, sumstats_df))
  
}

shorten_ALT_reporters_250nt_to_249nt = function(df, sumstats_df){
  
  n_before = nrow(df)
  
  # Subset affected reporters 
  tmp_reporters_df = df %>% 
    mutate(reporter_len = nchar(ALT_sequence_250bp)) %>%
    filter(reporter_len == 251)
  
  # Identify affected genes (since you need to shorten reference )
  genes_affected = unique(tmp_reporters_df$gene_name)
  
  tmp_reporters_df = df %>% 
    filter(gene_name %in% genes_affected)
  
  tmp_unaffected_reporters_df = df %>%
    filter(!(gene_name %in% genes_affected))
  
  # 1. Remove variants affecting position 250 
  tmp_reporters_df = tmp_reporters_df %>%
    filter(position_within_reporter_seq != 250)
  
  # 2. Shorten alternative reporters to 249nt 
  tmp_reporters_df = tmp_reporters_df %>%
    mutate(ALT_sequence_250bp = substr(ALT_sequence_250bp, 2, nchar(ALT_sequence_250bp)))
  
  # 3. Merge affected and unaffected reporters
  df = rbind(tmp_reporters_df, 
             tmp_unaffected_reporters_df)
  
  n_after = nrow(df)
  
  sumstats_df[nrow(sumstats_df)+1,] = c('Variants affecting position 250',
                                        n_before, n_after, n_after-n_before)
  
  return(list(df, sumstats_df))
  
}

keep_one_per_duplicate = function(ref_5utrs_df){
  
  duplicate_refs = ref_5utrs_df$ALT_sequence_250bp[duplicated(ref_5utrs_df$ALT_sequence_250bp)]
  
  refs_duplicates_df = ref_5utrs_df %>%
    filter(ALT_sequence_250bp %in% duplicate_refs) %>%
    group_by(ALT_sequence_250bp) %>%
    mutate(alternative_annotations = sapply(gene_name, function(g) {
      other_genes <- setdiff(gene_name, g)
      paste(other_genes, collapse = ', ')})) %>%
    ungroup()
  
  ref_not_duplicated_df = ref_5utrs_df %>%
    filter(!(ALT_sequence_250bp %in% duplicate_refs)) %>%
    mutate(alternative_annotations = NA)
  
  for (ref in duplicate_refs){
    
    tmp_df = refs_duplicates_df %>%
      filter(ALT_sequence_250bp == ref)
    
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


#------------#
# LOAD DATA  #
#------------#
ukb_df = data.frame(fread("./processed/UKB_observed_variants_annotated.csv"))
invar_df = data.frame(fread("./processed/UKB_invariant_variants_annotated_top40hs_AR.csv"))

ukb_sumstats_df = data.frame(fread("./processed/filtering_sumstats_all.csv"))
invar_sumstats_df = data.frame(fread("./processed/filtering_sumstats_all_top40hs_AR.csv"))

recessive_genes = data.frame(fread('./raw/all_recessive_genes.tsv', header=FALSE))
names(recessive_genes) = "gene_name"

ref_5utrs_df = data.frame(fread('./raw/5utr_ref_sequences_gencode_v45_annotated.csv'))
# The annotations in here are incorrect. In annotate_ukb_variants, I'd reannotated them! 

#---------------------#
# PROCESS REFERENCES  #
#---------------------#
ref_5utrs_df = ref_5utrs_df %>%
  rename(CHR = 'Chr',
         ALT_sequence_250bp = 'UTR5_sequence_250') %>%
  select(gene_name, CHR, start_codon, ALT_sequence_250bp, orf_ann_250bp, kozak_strength) %>%
  
  # Shorten references to 249nt if they are 250nt long
  mutate(reporter_len = nchar(ALT_sequence_250bp) - str_count(ALT_sequence_250bp, fixed("*"))) %>%
  mutate(ALT_sequence_250bp = ifelse(reporter_len == 250,
                                     substr(ALT_sequence_250bp, 2, nchar(ALT_sequence_250bp)),
                                     ALT_sequence_250bp))

# Note that there are duplicates!


#---------------#
# PROCESS DATA  #
#---------------#

# Annotate recessive genes 
ukb_df$recessive = ifelse(ukb_df$gene_name %in% recessive_genes$gene_name, "yes", 'no')

# Clean up invariant sites
invar_df$predicted_category = ifelse(invar_df$variant_type != 'No change' & 
                                       invar_df$predicted_category == 'Null',
                                     'Ambiguous', invar_df$predicted_category)

invar_df = invar_df %>% filter(hs_percentile >= 35)


#----------------------------------------------------#
# REMOVE VARIANTS NOT WITHIN 10% OF REF/ALT LENGTHS  #
#----------------------------------------------------#
ukb_res = remove_variants_not_within_10p_of_ref_alt_lengths(
  ukb_df, ukb_sumstats_df, 'ALT_sequence_250bp', 'ukb')

ukb_df = ukb_res[[1]]
ukb_sumstats_df = ukb_res[[2]]


invar_res = remove_variants_not_within_10p_of_ref_alt_lengths(
  invar_df, invar_sumstats_df, 'ALT_sequence_250bp', 'invariant')

invar_df = invar_res[[1]]
invar_sumstats_df = invar_res[[2]]


#-----------------------------#
# REMOVE DUPLICATE REPORTERS  #
#-----------------------------#
ukb_res = remove_duplicate_reporters(ukb_df, ukb_sumstats_df, 'ukb')

ukb_df = ukb_res[[1]]
ukb_sumstats_df = ukb_res[[2]]


invar_res = remove_duplicate_reporters(invar_df, invar_sumstats_df, 'invariant')

invar_df = invar_res[[1]]
invar_sumstats_df = invar_res[[2]]


#----------------------------------------------------------#
# UKB — Remove multinucleotide variants outside of 250nt   #
#----------------------------------------------------------#
ukb_res = remove_multinucleotide_variants_outside_25_250nt(ukb_df, ukb_sumstats_df, 'ukb')

ukb_df = ukb_res[[1]]
ukb_sumstats_df = ukb_res[[2]]


invar_res = remove_multinucleotide_variants_outside_25_250nt(invar_df, invar_sumstats_df, 'invariant')

invar_df = invar_res[[1]]
invar_sumstats_df = invar_res[[2]]


#---------------------------------------------------------#
# SHORTEN 250nt to 249nt if they have insertion variants  #
#---------------------------------------------------------#
invar_res = shorten_ALT_reporters_250nt_to_249nt(invar_df, invar_sumstats_df)

invar_df = invar_res[[1]]
invar_sumstats_df = invar_res[[2]]


#---------------------------------------------------------------------------------#
# REMOVE ANY INVARIANT VARIANTS PRODUCING A REFERENCE OR A UKB-OBSERVED SEQUENCE  #
#---------------------------------------------------------------------------------#
n_before = nrow(invar_df)

invar_df = invar_df %>%
  filter(!(ALT_sequence_250bp %in% c(ref_5utrs_df$ALT_sequence_250bp, ukb_df$ALT_sequence_250bp)))

invar_sumstats_df[nrow(invar_sumstats_df)+1, ] = 
  c('Reporters duplicated in UKB or references', n_before, nrow(invar_df), n_before-nrow(invar_df))


#--------------------------------------------#
# REMOVE UKB VARIANTS PRODUCING A REFERENCE  #
#--------------------------------------------#
n_before = nrow(ukb_df)

ukb_df = ukb_df %>%
  filter(!(ALT_sequence_250bp %in% as.character(ref_5utrs_df$ALT_sequence_250bp)))

ukb_sumstats_df[nrow(ukb_sumstats_df)+1, ] = 
  c('Reporters duplicated in references', n_before, nrow(ukb_df), length(unique(ukb_df$gene_name)),
    n_before - nrow(ukb_df))
  

#----------------#
# SANITY CHECKS  #
#----------------#

# 1. No duplicates anywhere [DONE]
n_unique_reporters_ukb = length(unique(ukb_df$ALT_sequence_250bp))
n_unique_rows_ukb = nrow(ukb_df)

n_unique_reporters_ukb == n_unique_rows_ukb

n_unique_reporters_invar = length(unique(invar_df$ALT_sequence_250bp))
n_unique_rows_invar = nrow(invar_df)

n_unique_reporters_invar == n_unique_rows_invar


# 2. Lengths of reporters (UKB and invariant, REF and ALT) between 25-250 [DONE]
range((ukb_df %>%
  mutate(reporter_len = nchar(ALT_sequence_250bp) - str_count(ALT_sequence_250bp, fixed("*"))))$reporter_len)

range((invar_df %>%
         mutate(reporter_len = nchar(ALT_sequence_250bp) - str_count(ALT_sequence_250bp, fixed("*"))))$reporter_len)

# 3. Combine ALL reporters 
all_reporters = c(ukb_df$ALT_sequence_250bp,
                  invar_df$ALT_sequence_250bp)

all_refs = ref_5utrs_df$ALT_sequence_250bp
# The only duplicates are within the reference sequences 

# 4. Check manually a few +/- stranded genes that:
#     - alternative sequences are legit             [DONE]
#     - position in reporter seq is correct         [DONE]
#     - length is accurate (especially if *)        [DONE]

# AIDA (-) CDC7 (+) FARSB (25, -) ABL2 (250, -) 

# tmp = ukb_df %>%
#   filter(gene_name %in% c("AIDA", "CDC7", "FARSB", "ABL2"))
# 
# tmp_ref = ref_5utrs_df %>%
#   filter(gene_name %in% c("AIDA", "CDC7", "FARSB", "ABL2"))

###################################################################################################

#-----------------------------------#
# VARIABLES ON # of VARS TO SELECT  #
#-----------------------------------#
n_vars_to_select = 900000

selected_vars_df = data.frame(matrix(nrow=0, ncol=ncol(ukb_df)))
names(selected_vars_df) = names(ukb_df)

summary_df = data.frame(matrix(nrow=0, ncol=3))
names(summary_df) = c('Category', 'N_variants', 'N_genes')


#------------------------------------------------------------------#
# 1. KEEP ALL OBSERVED VARIANTS IN HIGH S-HET AND RECESSIVE GENES  #
#------------------------------------------------------------------#

# Split by high and low shet genes
df_low_hs = ukb_df %>%
  filter(hs_percentile <= 65) %>%
  filter(recessive == 'no')

df_high_hs = ukb_df %>%
  filter(hs_percentile > 65 | recessive == 'yes')

mean_utr_len_low_hs = round(mean((df_low_hs %>% select(gene_name, seq_len) %>% distinct())$seq_len), 1)
mean_utr_len_high_hs = round(mean((df_high_hs %>% select(gene_name, seq_len) %>% distinct())$seq_len), 1)

print(str_interp("Average 5' UTR length in low shet genes is ${mean_utr_len_low_hs}."))
print(str_interp("Average 5' UTR length in high shet and recessive genes is ${mean_utr_len_high_hs}."))

n_low_hs = nrow(df_low_hs)
n_high_hs = nrow(df_high_hs)

n_genes_low_hs = length(unique(df_low_hs$gene_name))
n_genes_high_hs = length(unique(df_high_hs$gene_name))

print(str_interp("There are ${n_low_hs} variants in low shet genes, and ${n_genes_low_hs} genes."))
print(str_interp("There are ${n_high_hs} variants in high shet genes, and ${n_genes_high_hs} genes."))


#################################################
# KEEP all observed variants in high hs genes   #

n_vars_to_select = n_vars_to_select-n_high_hs
selected_vars_df = rbind(selected_vars_df, df_high_hs)

summary_df[nrow(summary_df)+1, ] = 
  c('Observed variants in high s-het and recessive genes', n_high_hs, n_genes_high_hs)

#################################################


#--------------------------------------------#
# 2. KEEP ALL MAF > 0.1% IN LOW S-HET GENES  #
#--------------------------------------------#
df_low_hs_maf_1p = df_low_hs %>%
  filter(AF >= 0.001)

n_vars = nrow(df_low_hs_maf_1p)
print(str_interp("There are ${n_vars} variants in low s-het genes, MAF >= 0.1%."))


#################################################
# KEEP all MAF > 0.1% variants in low shet genes  #

n_vars_to_select = n_vars_to_select-n_vars
selected_vars_df = rbind(selected_vars_df, df_low_hs_maf_1p)

summary_df[nrow(summary_df)+1, ] = 
  c('MAF > 0.1% in low s-het genes', n_vars, NA)
#################################################


#-----------------------------------------------------------------#
# 3. RANDOM SAMPLING MAF < 0.1% IN LOW S-HET GENES UP TO 900,000  #
#-----------------------------------------------------------------#
set.seed(123456)

df_sampled = df_low_hs %>%
  filter(AF < 0.001) %>%
  sample_n(n_vars_to_select, replace=FALSE)

n_vars = nrow(df_sampled)
print(str_interp("We sampled ${n_vars} variants in low s-het genes, MAF < 0.1%."))


#################################################
# SAMPLE REST FROM MAF < 1%  in low shet genes  #

n_vars_to_select = n_vars_to_select-n_vars
selected_vars_df = rbind(selected_vars_df, df_sampled)

summary_df[nrow(summary_df)+1, ] = 
  c('MAF < 0.1% in low s-het genes', n_vars, NA)
#################################################


#----------------------------#
# 4. ADD REFERENCE SEQUENCES #
#----------------------------#
panel_genes = unique(selected_vars_df$gene_name)
n_panel_genes = length(panel_genes)

ref_5utrs_panel_raw_df = ref_5utrs_df %>%
  filter(gene_name %in% panel_genes) %>%
  mutate(orf_ann_250bp = NA,
         gene_name = paste(gene_name, "_ref", sep=''))

ref_5utrs_panel_df = keep_one_per_duplicate(ref_5utrs_panel_raw_df)
n_reporters = nrow(ref_5utrs_panel_df)

print(str_interp("There are ${n_panel_genes} genes and ${n_reporters} unique reporters."))

# Identify missing columns in ref_5utrs_panel_df and add NA
missing_columns = setdiff(names(selected_vars_df), names(ref_5utrs_panel_df))

for (col_name in missing_columns) {
  ref_5utrs_panel_df[[col_name]] = NA
}

summary_df[nrow(summary_df)+1, ] = 
  c('Reference genes', n_reporters, n_panel_genes)


#--------------------------------#
# 5. ADDING UNOBSERVED VARIANTS  #
#--------------------------------#
n_invar_sites = 100000 - nrow(ref_5utrs_panel_df)

set.seed(123456)
selected_invar_null_df = invar_df %>%
  filter(predicted_category == "Null") %>%
  sample_n(0.1*n_invar_sites)

set.seed(123456)
selected_invar_pos_df = invar_df %>%
  filter(predicted_category == "Positive") %>%
  sample_n(0.5*n_invar_sites+2)

set.seed(123456)
selected_invar_nonnull_df = invar_df %>%
  filter(!(predicted_category %in% c("Null", 'Positive'))) %>%
  sample_n(0.4*n_invar_sites)

selected_invar_df = rbind(selected_invar_null_df, 
                          selected_invar_pos_df,
                          selected_invar_nonnull_df)

n_null = nrow(selected_invar_null_df)
n_pos = nrow(selected_invar_pos_df)
n_nonnull = nrow(selected_invar_nonnull_df)

summary_df[nrow(summary_df)+1, ] = c('Positive variants', n_pos, NA)
summary_df[nrow(summary_df)+1, ] = c('ORF/Kozak-disrupting variants', n_nonnull, NA)
summary_df[nrow(summary_df)+1, ] = c('Null variants', n_null, NA)


#-------------------------#
# ANALYZE RESULTING PANEL #
#-------------------------#
final_panel_df = rbind(selected_vars_df %>%
                         rename(insert_seq = 'ALT_sequence_250bp') %>%
                         mutate(alternative_annotations = NA), 
                       selected_invar_df %>%
                         rename(insert_seq = 'ALT_sequence_250bp') %>%
                         mutate(AC = NA, AN = NA, AF = NA, recessive=NA,
                                alternative_annotations = NA),
                       ref_5utrs_panel_df %>%
                         rename(insert_seq = 'ALT_sequence_250bp'))

n_variants = nrow(final_panel_df)
n_unique_reporters = length(unique(final_panel_df$insert_seq))

summary_df[nrow(summary_df)+1, ] = c('Final panel', n_variants, n_panel_genes)

print(str_interp("The total number of variants is ${n_variants} and of reporters is ${n_unique_reporters}."))


## 1. Number of variants per gene ------------------------------------
gene_df = final_panel_df %>%
  group_by(gene_name) %>%
  dplyr::mutate(n_variants = n(),
                hs_high_low = ifelse(hs_decile <=6, 'low shet', 'high shet')) %>%
  ungroup() %>%
  select(gene_name, contains("hs"), seq_len, n_variants, hs_high_low) %>%
  distinct()

p = ggplot(gene_df, aes(x=n_variants)) +
  geom_histogram(binwidth = 5) +
  xlab("Number of observed variants per gene") + 
  ylab("Count")

add_tami_theme(p)

range(gene_df$n_variants)

png("./panel_outputs/n_vars_per_gene.png", width=450, height=350)
add_tami_theme(p)
dev.off()


## 2. Fraction of positive/negative/null variants ----------------------
freq_df = make_freq_df(final_panel_df, 'predicted_category',
                       c('Positive', 'Negative', 'Ambiguous', 'Null'))

p = ggplot(freq_df, aes(x=Variable, y=Percentage, fill=Variable)) +
  geom_bar(stat='identity') +
  geom_text(aes(label=comma(Count)), 
            vjust=-0.5, size=6) +
  xlab("\nPredicted effect category") + 
  ylab("Percentage (%)") +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_manual(values=c('red', 'blue', 'dark grey', 'light grey')) +
  ylim(c(0,100))

add_tami_theme(p)

png("./panel_outputs/variant_types_all.png", width=450, height=350)
add_tami_theme(p)
dev.off()

## Only observed 
freq_df = make_freq_df(selected_vars_df, 'predicted_category',
                       c('Positive', 'Negative', 'Ambiguous', 'Null'))

p = ggplot(freq_df, aes(x=Variable, y=Percentage, fill=Variable)) +
  geom_bar(stat='identity') +
  geom_text(aes(label=comma(Count)), 
            vjust=-0.5, size=6) +
  xlab("\nPredicted effect category") + 
  ylab("Percentage (%)") +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_manual(values=c('red', 'blue', 'dark grey', 'light grey')) +
  ylim(c(0,100)) +
  ggtitle("Observed UKB variants (N=900,000)")

add_tami_theme(p)

png("./panel_outputs/variant_types_observed.png", width=450, height=350)
add_tami_theme(p)
dev.off()

## Only unobserved
freq_df = make_freq_df(selected_invar_df, 'predicted_category',
                       c('Positive', 'Negative', 'Ambiguous', 'Null'))

p = ggplot(freq_df, aes(x=Variable, y=Percentage, fill=Variable)) +
  geom_bar(stat='identity') +
  geom_text(aes(label=comma(Count)), 
            vjust=-0.5, size=6) +
  xlab("\nPredicted effect category") + 
  ylab("Percentage (%)") +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_manual(values=c('red', 'blue', 'dark grey', 'light grey')) +
  ylim(c(0,100)) +
  ggtitle("Non-observed variants (N=100,000)")

add_tami_theme(p)

png("./panel_outputs/variant_types_unobserved.png", width=450, height=350)
add_tami_theme(p)
dev.off()


## 3. Fractions of indels 
final_panel_df = annotate_vartype(final_panel_df)

freq_df = make_freq_df(final_panel_df %>% filter(is.na(vartype) == FALSE), 'vartype',
                       rev(c('SNP', 'SNP, indel', 'Larger variant', 'Reference')))

p = ggplot(freq_df, aes(x=Variable, y=Percentage)) +
  geom_bar(stat='identity') +
  geom_text(aes(label=comma(Count)), 
            vjust=0.5, hjust=-0.2, size=6) +
  xlab("") + 
  ylab("Percentage (%)") +
  theme_classic() +
  theme(legend.position = 'none') +
  ylim(c(0,100)) +
  coord_flip()

add_tami_theme(p)

png("./panel_outputs/snps_indels_large_variants.png", width=450, height=350)
add_tami_theme(p)
dev.off()


## 4. Frequencies
selected_vars_df$freq = NA
selected_vars_df$freq = ifelse(selected_vars_df$AF <= 0.05, 'MAF 1-5%', selected_vars_df$freq)
selected_vars_df$freq = ifelse(selected_vars_df$AF <= 0.01, 'MAF 0.1-1%', selected_vars_df$freq)
selected_vars_df$freq = ifelse(selected_vars_df$AF <= 0.001, 'MAF 0.01-0.1%', selected_vars_df$freq)
selected_vars_df$freq = ifelse(selected_vars_df$AF <= 0.0001, 'MAF 0.001-0.01%', selected_vars_df$freq)
selected_vars_df$freq = ifelse(selected_vars_df$AF <= 0.00001, 'MAF < 0.001%', selected_vars_df$freq)
selected_vars_df$freq = ifelse(selected_vars_df$AC == 1, 'Singleton', selected_vars_df$freq)

freq_df = make_freq_df(selected_vars_df, 'freq',
                       order=c('Singleton', 'MAF < 0.001%', 'MAF 0.001-0.01%', 'MAF 0.01-0.1%',
                               'MAF 0.1-1%', 'MAF 1-5%'))

p = ggplot(freq_df, aes(x=Variable, y=Percentage)) +
  geom_bar(stat='identity', fill='#ECB159') +
  geom_text(aes(label=comma(Count)), 
            vjust=0.5, hjust=-0.2, size=6) +
  xlab("") + 
  ylab("Percentage (%)") +
  theme_classic() +
  theme(legend.position = 'none') +
  ylim(c(0,105)) +
  coord_flip()

add_tami_theme(p)

png("./panel_outputs/frequencies.png", width=500, height=300)
add_tami_theme(p)
dev.off()


#-----------------------#
# SAVE RESULTING PANEL  #
#-----------------------#

# Make a clean file to send to Yale
final_panel_clean_df = final_panel_df %>%
  rename(ref_length = 'seq_len',
         alt_length = 'reporter_len',
         shet = 'posterior_hs') %>%
  select(gene_name:AF, posterior_hs, start_codon, insert_seq, 
         position_within_reporter_seq, ref_length, alt_length,
         orf_ann_250bp, orf_ann_250bp_ALT, 
         variant_type, predicted_category, 
         alternative_annotations)

write.table(final_panel_clean_df, './panel_outputs/naptrap_final_panel_1m_variants.txt', sep='\t',
            row.names=FALSE, quote=FALSE)

# Full NaP-TRAP data and summary 
write.table(final_panel_df, './panel_outputs/naptrap_final_panel_1m_variants_annotated.txt', sep='\t',
            row.names=FALSE, quote=FALSE)

write.table(summary_df, './panel_outputs/naptrap_final_panel_sumstats.txt', sep='\t',
            row.names=FALSE, quote=FALSE)

# Updated UKB and invariant sumstats
write.csv(ukb_sumstats_df, './panel_outputs/filtering_sumstats_all.csv',
          row.names=FALSE, quote=FALSE)

write.csv(invar_sumstats_df, './panel_outputs/filtering_sumstats_all_top40hs_AR.csv',
          row.names=FALSE, quote=FALSE)


