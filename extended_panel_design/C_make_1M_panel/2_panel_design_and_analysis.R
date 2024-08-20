#-----------------#
# PANEL ANALYSIS  #
#-----------------#
# NOTE that this code was run locally.

rm(list=ls())

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)

set.seed(123456)

setwd('/home/users/tami/5utr_extended_panel/data')


#------------#
# FUNCTIONS  #
#------------#
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

make_freq_df = function(df, colname, order=NA, add_refs=FALSE){
  
  if (add_refs){
    df = df %>% 
      mutate(predicted_category = ifelse(is.na(predicted_category),
                                         'Reference', predicted_category))
    order = c(order, 'Reference')
  }
  
  freq_df = data.frame(table(df[,colname]))
  names(freq_df) = c('Variable', 'Count')
  
  freq_df$Percentage = round(freq_df$Count/sum(freq_df$Count), 3)*100
  
  if (length(order) != 1){
    
    freq_df$Variable = factor(freq_df$Variable, levels=order)
    
    freq_df = freq_df[order(freq_df$Variable), ]
  }
  
  return(freq_df)
  
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


#------------#
# LOAD DATA  # 
#------------#
panel_750k = data.frame(fread('./processed/panel_750k.csv'))
panel_250k = data.frame(fread('./processed/panel_250k.csv'))

final_panel_df = rbind(panel_750k %>% mutate(panel = '750k'), 
                       panel_250k %>% mutate(panel = '250k'))

sumstats_750k = data.frame(fread('./processed/panel_750k_sumstats.csv'))
sumstats_250k = data.frame(fread('./processed/panel_250k_sumstats.csv'))


## 0a. Basic numbers
n_reporters = nrow(final_panel_df)
n_genes = length(unique((final_panel_df %>%
                          filter(grepl("ref", gene_name)))$gene_name))
n_refs = nrow(final_panel_df %>% filter(grepl("ref", gene_name)))

print(str_interp("${n_reporters} reporters, ${n_genes} genes, ${n_refs} references"))

## 0b. References 
refs_df = data.frame(table((final_panel_df %>%
                              filter(grepl("ref", gene_name)))$gene_name))

n_once = as.numeric(table(refs_df$Freq)[1])
n_twice = as.numeric(table(refs_df$Freq)[2])

print(str_interp("There are ${n_refs} reporters and ${n_genes} genes."))
print(str_interp("${n_once} genes have one, and ${n_twice} genes have two reporters."))


## 1. Number of variants per gene ------------------------------------
gene_df = final_panel_df %>%
  group_by(gene_name) %>%
  dplyr::mutate(n_variants = n(),
                hs_high_low = ifelse(hs_decile <=6, 'low shet', 'high shet')) %>%
  ungroup() %>%
  filter(!(grepl('ref', gene_name))) %>%
  select(gene_name, contains("hs"), n_variants, hs_high_low) %>%
  distinct()

p = ggplot(gene_df, aes(x=n_variants)) +
  geom_histogram(binwidth = 5) +
  xlab("Number of observed variants per gene") + 
  ylab("Count")

add_tami_theme(p)

nrow(gene_df)
range(gene_df$n_variants)

png("./naptrap_panel/n_vars_per_gene.png", width=450, height=350)
add_tami_theme(p)
dev.off()


## 2. Fraction of positive/negative/null variants ----------------------
freq_df = make_freq_df(final_panel_df, 'predicted_category',
                       c('Positive', 'Negative', 'Total loss', 'Null'), add_refs = TRUE)

p = ggplot(freq_df, aes(x=Variable, y=Percentage, fill=Variable)) +
  geom_bar(stat='identity') +
  geom_text(aes(label=comma(Count)), 
            vjust=-0.5, size=6) +
  xlab("\nPredicted effect category") + 
  ylab("Percentage (%)") +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_manual(values=c('red', 'blue', '#00008B', 'grey', 'black')) +
  ylim(c(0,100))

add_tami_theme(p)

png("./naptrap_panel/variant_types_all.png", width=450, height=350)
add_tami_theme(p)
dev.off()


## 2a. Only observed ---------------------------------------------------
freq_df = make_freq_df(final_panel_df %>%
                         filter(is.na(AF) == FALSE), 
                       'predicted_category',
                       c('Positive', 'Negative', 'Total loss', 'Null'))

p = ggplot(freq_df, aes(x=Variable, y=Percentage, fill=Variable)) +
  geom_bar(stat='identity') +
  geom_text(aes(label=comma(Count)), 
            vjust=-0.5, size=6) +
  xlab("\nPredicted effect category") + 
  ylab("Percentage (%)") +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_manual(values=c('red', 'blue', '#00008B', 'grey')) +
  ylim(c(0,100)) +
  ggtitle("Observed UKB variants (N=900,000)")

add_tami_theme(p)

png("./naptrap_panel/variant_types_observed.png", width=450, height=350)
add_tami_theme(p)
dev.off()


## 2b. Only unobserved -------------------------------------------------
freq_df = make_freq_df(final_panel_df %>%
                         filter(is.na(AF)), 
                       'predicted_category',
                       c('Positive', 'Negative', 'Total loss', 'Null'))

n_unobserved = sum(freq_df$Count)

p = ggplot(freq_df, aes(x=Variable, y=Percentage, fill=Variable)) +
  geom_bar(stat='identity') +
  geom_text(aes(label=comma(Count)), 
            vjust=-0.5, size=6) +
  xlab("\nPredicted effect category") + 
  ylab("Percentage (%)") +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_manual(values=c('red', 'blue', '#00008B', 'grey')) +
  ylim(c(0,100)) +
  ggtitle(str_interp("Non-observed variants (N=${n_unobserved})"))

add_tami_theme(p)

png("./naptrap_panel/variant_types_unobserved.png", width=450, height=350)
add_tami_theme(p)
dev.off()


## 3. Fractions of indels ----------------------------------------------
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

png("./naptrap_panel/snps_indels_large_variants.png", width=450, height=350)
add_tami_theme(p)
dev.off()


## 4. Frequencies ------------------------------------------------------
selected_vars_df = final_panel_df %>%
  filter(is.na(AF) == FALSE)

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

png("./naptrap_panel/frequencies.png", width=500, height=300)
add_tami_theme(p)
dev.off()


#-----------------------#
# SAVE RESULTING PANEL  #
#-----------------------#

# Make a clean file to send to Yale
final_panel_clean_df = final_panel_df %>%
  rename(reporter_length = 'seq_len',
         shet = 'posterior_hs') %>%
  select(gene_name:AF, shet, start_codon, insert_seq, 
         position_within_reporter_seq, reporter_length,
         orf_ann_ref, orf_ann_alt, 
         variant_type, predicted_category, 
         alternative_annotations, panel)


# NaP-TRAP panel for Yale 
write.table(final_panel_clean_df, './naptrap_panel/naptrap_final_panel_1m_variants.txt', sep='\t',
            row.names=FALSE, quote=FALSE)

# Full NaP-TRAP data and summary 
write.table(final_panel_df, './naptrap_panel/naptrap_final_panel_1m_variants_annotated.txt', sep='\t',
            row.names=FALSE, quote=FALSE)
