rm(list=ls())

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)

setwd('/home/users/tami/5utr_extended_panel/data/')
#setwd('~/Desktop/5prime_utr/data/')


#--------------------#
# PROCESS ARGUMENTS  #
#--------------------#
args = commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  stop("You must specify the fraction of individuals that should be in the reference list", call.=FALSE)
  
} else if (length(args) != 0) {
  chr = args[1]
}


#------------#
# FUNCTIONS  #
#------------#

# Annotate gene names
check_overlap <- function(chrom, pos, start, end) {
  overlap = pos >= start & pos <= end
  return(overlap)
}

annotate_gene_names = function(ukb_df, gencode_df){
  
  n_before = nrow(ukb_df)
  
  gencode_df <- as.data.frame(gencode_df)
  ukb_df <- as.data.frame(ukb_df)
  
  for (i in 1:nrow(ukb_df)){
    
    CHR = ukb_df[i,'CHR']
    POS = ukb_df[i,'POS']
    
    subset_genes = subset(gencode_df, Chr == CHR &
                            check_overlap(Chr, POS, Start, End))
    
    if (nrow(subset_genes) > 0) {
      ukb_df[i,'gene_name'] = paste(subset_genes$gene_name, collapse=", ")
    } else {
      ukb_df[i,'gene_name'] = NA
    }
  }
  
  ukb_df = ukb_df %>%
    filter(is.na(gene_name) == FALSE)
  
  n_after = nrow(ukb_df)
  
  # Variants
  print("Remove variants that couldnt be linked to genes")
  print(str_interp("Before filtering: ${n_before} variants"))
  print(str_interp("After filtering: ${n_after} variants"))
  print(str_interp("Dropped ${abs(n_after-n_before)} variants"))
  
  return(ukb_df)
}


# Variants that fall in 5' UTRs of two genes
remove_vars_in_two_genes = function(ukb_df){
  
  # How many variants
  n_before = nrow(ukb_df)
  overlapped_genes = unique(unlist(str_split(ukb_df[grepl(",", ukb_df$gene_name), 'gene_name'], ", ")))
  
  # Divide single and double genes
  df_single_genes = ukb_df %>%
    filter(!(grepl(",", gene_name)))
  
  n_after = nrow(df_single_genes)
  
  # Variants
  print("Remove variants overlapping multiple gene 5' UTRs")
  print(str_interp("Before filtering: ${n_before} variants"))
  print(str_interp("After filtering: ${n_after} variants"))
  print(str_interp("Dropped ${abs(n_after-n_before)} variants"))
  
  print("------------------")
  
  # Genes
  for (gene in overlapped_genes){
    n_vars = sum(grepl(gene, df_single_genes$gene_name))
    print(str_interp("For ${gene}, there are ${n_vars} remaining variants."))
  }
  
  return(df_single_genes)
}


# Position AUGs in contiguous sequence
add_adjusted_start_and_end_to_gencode = function(gencode_df){
  
  new_gencode_df = data.frame(matrix(nrow=0, ncol=ncol(gencode_df)+2))
  names(new_gencode_df) = c(names(gencode_df), 'Adjusted_start', 'Adjusted_end')
  
  for (gene in unique(gencode_df$gene_name)){
    
    # Subset GENCODE exons for gene
    gen_df = gencode_df %>%
      filter(gene_name == gene) %>%
      filter(Element == 'five_prime_UTR') %>%
      mutate(exon_length = abs(End-Start)+1) %>%
      arrange(Start)
    
    strand = unique(gen_df$Strand)
    
    # Positively stranded genes
    if (strand == '+'){
      
      previous_end = 0
      
      for (i in 1:nrow(gen_df)){
        gen_df[i, 'Adjusted_start'] = previous_end
        gen_df[i, 'Adjusted_end'] = previous_end + gen_df[i,'exon_length'] - 1
        
        previous_end = gen_df[i, 'Adjusted_end']+1
      }
      
      # Negatively stranded genes
    } else if (strand == '-'){
      
      previous_end = 0
      
      for (i in nrow(gen_df):1){
        gen_df[i, 'Adjusted_start'] = previous_end
        gen_df[i, 'Adjusted_end'] = previous_end + gen_df[i,'exon_length'] - 1
        
        previous_end = gen_df[i, 'Adjusted_end']+1
      }
    }
    
    new_gencode_df = rbind(new_gencode_df, gen_df)
  }
  
  return(new_gencode_df)
  
}

compute_position_within_reporter_seq = function(ukb_df, gencode_df){
  
  new_ukb_df = data.frame(matrix(nrow=0, ncol=ncol(ukb_df)+1))
  names(new_ukb_df) = c(names(ukb_df), 'position_within_reporter_seq')
  
  for (gene in unique(ukb_df$gene_name)){
    
    tmp_ukb_df = ukb_df %>%
      filter(gene_name == gene)
    
    tmp_gen_df = gencode_df %>%
      filter(gene_name == gene) %>%
      arrange(Start)
    
    # Identify position_from_AUG for UKB chromosomal position
    strand = unique(tmp_gen_df$Strand)
    
    exons_seq = c()
    
    for (i in 1:nrow(tmp_gen_df)){
      exons_seq = c(exons_seq, seq(as.numeric(tmp_gen_df[i,'Start']),
                                   as.numeric(tmp_gen_df[i,'End'])))
    }
    
    expanded_tmp_gen_df = data.frame(exons_seq)
    names(expanded_tmp_gen_df) = 'POS'
    
    if (strand == "+"){
      expanded_tmp_gen_df$adjusted_pos = seq(nrow(expanded_tmp_gen_df), 1)
    } else if (strand == "-"){
      expanded_tmp_gen_df$adjusted_pos = seq(1, nrow(expanded_tmp_gen_df))
    }
    
    for (j in 1:nrow(tmp_ukb_df)){
      tmp_ukb_df[j,'position_within_reporter_seq'] =
        expanded_tmp_gen_df[which(expanded_tmp_gen_df$POS == tmp_ukb_df[j,'POS']), 'adjusted_pos']
    }
    
    # Attach gene data to new UKB df
    new_ukb_df = rbind(new_ukb_df, tmp_ukb_df)
    
  }
  
  new_ukb_df = remove_multiposition_variants_outside_utr_len(new_ukb_df)
  
  return(new_ukb_df)
  
}

# For multi-position variants, remove them if they fall outside 5' UTR
remove_multiposition_variants_outside_utr_len = function(new_ukb_df){
  
  n_before = nrow(new_ukb_df)
  
  new_ukb_snps_df = new_ukb_df %>%
    filter(nchar(REF) == 1 & nchar(ALT) == 1)
  
  new_ukb_longer_df = new_ukb_df %>%
    filter(!(nchar(REF) == 1 & nchar(ALT) == 1))
  
  new_ukb_longer_df = merge(new_ukb_longer_df,
                            ref_5utrs_df %>%
                              mutate(seq_len = nchar(UTR5_sequence_250)) %>%
                              select(gene_name, Strand, seq_len), by='gene_name')
  
  for (i in 1:nrow(new_ukb_longer_df)){
    
    strand = new_ukb_longer_df[i,'Strand']
    REF = new_ukb_longer_df[i,'REF']
    ALT = new_ukb_longer_df[i,'ALT']
    position_within_reporter_seq = new_ukb_longer_df[i,'position_within_reporter_seq']
    seq_len = new_ukb_longer_df[i,'seq_len']
    
    len_REF = nchar(REF)
    len_ALT = nchar(ALT)
    
    end_REF = ifelse(strand == "+", position_within_reporter_seq + len_REF,
                     position_within_reporter_seq - len_ALT)
    end_ALT = ifelse(strand == "+", position_within_reporter_seq + len_ALT,
                     position_within_reporter_seq - len_ALT)
    
    if (end_REF %in% seq(1, seq_len) && end_ALT %in% seq(1, seq_len)){
      new_ukb_longer_df[i,'out_of_range'] = 'No'
    } else {
      new_ukb_longer_df[i,'out_of_range'] = 'Yes'
    }
    
  }
  
  new_ukb_longer_df = new_ukb_longer_df %>%
    filter(out_of_range == "No") %>%
    select(-out_of_range, -Strand, -seq_len)
  
  new_ukb_df = rbind(new_ukb_snps_df, new_ukb_longer_df)
  
  n_after = nrow(ukb_df)
  
  print(str_interp("Removing multiposition varinats outside of UTR length:"))
  print(str_interp("Before filtering: ${n_before} variants"))
  print(str_interp("After filtering: ${n_after} variants"))
  print(str_interp("Dropped ${abs(n_after-n_before)} variants"))
  
  return(new_ukb_df)
  
}


# Complementary REF/ALT for - stranded genes
replace_acgt <- function(input_sequence) {
  output_sequence <- chartr("ACTG", "TGAC", input_sequence)
  return(output_sequence)
}

complementary_ref_alt_neg_strand = function(ukb_df, ref_5utrs_df){
  
  tmp_df = merge(ukb_df, ref_5utrs_df %>%
                   select(gene_name, Strand), by='gene_name')
  
  ukb_pos_df = tmp_df %>%
    filter(Strand == "+") %>%
    select(-Strand)
  
  ukb_neg_df = tmp_df %>%
    filter(Strand == "-") %>%
    mutate(REF = replace_acgt(REF),
           ALT = replace_acgt(ALT)) %>%
    select(-Strand)
  
  tmp_df = rbind(ukb_pos_df, ukb_neg_df)
  
  return(tmp_df)
  
}


# Remove variants >250 bp from start codon
remove_vars_further_than_limit = function(ukb_df, limit){
  
  n_before = nrow(ukb_df)
  
  # Split to SNPs and otherwise
  ukb_snps_df = ukb_df %>% filter(nchar(REF) == 1 & nchar(ALT) == 1)
  ukb_larger_df = ukb_df %>% filter(!(nchar(REF) == 1 & nchar(ALT) == 1))
  
  # Filter among SNPs
  ukb_snps_df = ukb_snps_df %>%
    filter(position_within_reporter_seq <= limit)
  
  # Filter among longer variants
  ukb_larger_df = calculate_max_end(ukb_larger_df)
  ukb_larger_df = ukb_larger_df %>%
    filter(max_end <= limit) %>%
    select(-max_end)
  
  # Merge them back
  ukb_df = rbind(ukb_snps_df, ukb_larger_df)
  
  n_after = nrow(ukb_df)
  
  print(str_interp("Removing variants further than ${limit} bp from start codon:"))
  print(str_interp("Before filtering: ${n_before} variants"))
  print(str_interp("After filtering: ${n_after} variants"))
  print(str_interp("Dropped ${abs(n_after-n_before)} variants"))
  
  return(ukb_df)
}

calculate_max_end = function(ukb_larger_df){
  
  for (i in 1:nrow(ukb_larger_df)){
    position_within_reporter_seq = ukb_larger_df[i, 'position_within_reporter_seq']
    REF = ukb_larger_df[i,'REF']
    ALT = ukb_larger_df[i,'ALT']
    
    end_ref = position_within_reporter_seq + nchar(REF)-1
    end_alt = position_within_reporter_seq + nchar(ALT)-1
    
    max_end = max(end_ref, end_alt)
    
    ukb_larger_df[i,'max_end'] = max_end
    
  }
  
  return(ukb_larger_df)
  
}


# Construct alternative sequences
construst_alt_sequences = function(ukb_df, ref_5utrs_df){
  
  ukb_df = merge(ukb_df, ref_5utrs_df %>%
                   select(gene_name, UTR5_sequence_250, start_codon, orf_ann_250bp, kozak_strength),
                 by='gene_name')
  
  for (i in 1:nrow(ukb_df)){
    
    ref = ukb_df[i,'REF']
    alt = ukb_df[i,'ALT']
    seq = ukb_df[i,'UTR5_sequence_250']
    gene = ukb_df[i,'gene_name']
    
    utr_len = nchar(seq)
    pos = utr_len - ukb_df[i,'position_within_reporter_seq'] + 1
    
    # Gene strand
    strand = ref_5utrs_df[which(ref_5utrs_df$gene_name == gene), 'Strand']
    
    add = nchar(ref)-1
    base_at_pos = ifelse(strand == "+", substr(seq, pos, pos+add),
                         substr(seq, pos-add, pos))
    
    if (strand == "-" && (nchar(ref) != 1 | nchar(alt) != 1)){
      base_at_pos = intToUtf8(rev(utf8ToInt(base_at_pos)))
    }
    
    if (base_at_pos == ref) {
      
      len_ref = nchar(ref)
      
      alt_seq = ifelse(strand == "+",
                       paste0(substr(seq, 1, pos-1), alt, substr(seq, pos+len_ref, nchar(seq))),
                       paste0(substr(seq, 1, pos-len_ref), alt, substr(seq, pos+1, nchar(seq))))
      
    } else {
      
      ukb_df[i,'ALT_sequence_250bp'] = 'mismatch UKB REF to GENCODE REF'
    }
    
    ukb_df[i,'ALT_sequence_250bp'] = alt_seq
  }
  
  ukb_df = ukb_df %>%
    filter(ALT_sequence_250bp != "mismatch UKB REF to GENCODE REF")
  
  return(ukb_df)
  
}



#------------#
# LOAD DATA  #
#------------#

# UKB variants
ukb_df = read.csv('./processed/ukb_variants_all.csv')
ukb_df = ukb_df %>%
  distinct() %>%
  filter(CHR == chr)

# Reference 5' UTRs
ref_5utrs_df = data.frame(fread('./processed/5utr_ref_sequences_gencode_v45_annotated.csv'))

# GENCODE exons (5' UTRs)
gencode_exons_df = data.frame(fread('./processed/gencode_v45_protein_coding_MANE_clean.csv'))
gencode_exons_df = gencode_exons_df %>%
  filter(Element == 'five_prime_UTR') %>%
  filter(Chr %in% c(str_interp('chr${chr}')))

# GENCODE exons (6bp of CDS)
gencode_6bp_df = data.frame(fread('./processed/gencode_v45_6bp_coordinates.csv'))
gencode_6bp_df = gencode_6bp_df %>%
  filter(Chr %in% c(str_interp('chr${chr}'))) %>%
  select(-exon_length) %>%
  filter(gene_name %in% gencode_exons_df$gene_name)

gencode_df = rbind(gencode_exons_df, gencode_6bp_df)
gencode_df$Chr = str_replace_all(gencode_df$Chr, "chr", "")


#--------------------------#
# LOAD GENE FEATURES DATA  #
#--------------------------#
df_ann_raw = data.frame(fread('./processed/ann_gencode_v45_protein_coding_MANE_clean.csv'))

df_hs_raw = data.frame(fread('./raw/s_het_estimates.genebayes.tsv'))
df_gene_features_raw = read.delim2('./raw/pc_genes.txt')


# Process other files to merge
df_ann = df_ann_raw %>%
  filter(gene_type == 'protein_coding') %>%
  select(hgnc_id, gene_name) %>%
  dplyr::rename(hgnc = hgnc_id) %>%
  unique()

df_hs = df_hs_raw %>%
  dplyr::rename(posterior_hs = 'post_mean') %>%
  select(hgnc, posterior_hs) %>%
  mutate(hs_decile = ntile(posterior_hs, 10),
         hs_by_median = ntile(posterior_hs, 2),
         hs_percentile = ntile(posterior_hs, 100))

# TSSD = gene density (# of TSSs within 1Mb window around a gene's TSS)
df_gene_features = df_gene_features_raw %>%
  dplyr::rename(hgnc = hgnc_id) %>%
  select(hgnc, TSSD, length, CDS_length, pLI, LOEUF,
         promoter_count, TF)

df_hs_and_features = merge(df_hs, df_gene_features, by='hgnc', all=TRUE)
df_hs_and_features = df_hs_and_features %>% distinct()

df_gene_info = merge(df_ann, df_hs_and_features, by='hgnc', all.x=FALSE, all.y=TRUE)

df_gene_info = df_gene_info %>%
  filter(is.na(gene_name) == FALSE) %>%
  distinct()


#------------------#
# BASIC FILTERING  #
#------------------#
summary_df = data.frame(matrix(nrow=0, ncol=4))
names(summary_df) = c('Filter', 'Before', 'After', 'N_Genes')

n_before = nrow(ukb_df)

# 0a. Remove varants with AC == 0
ukb_df = ukb_df %>% filter(AC != 0) 
n_genes = length(unique(ukb_df$gene_name))

summary_df[nrow(summary_df)+1, ] = c("Variants with allele count = 0",
                                     n_before, nrow(ukb_df), n_genes)
n_before = nrow(ukb_df)

# 0b. Remove varants with AN covering 50% of individuals 
max_AN = max(ukb_df$AN)
ukb_df = ukb_df %>% filter(AN >= max_AN*0.9)

summary_df[nrow(summary_df)+1, ] = c("Variants with allele number < 90% of individuals",
                                     n_before, nrow(ukb_df), n_genes)
n_before = nrow(ukb_df)

# 0c. Calculate allele frequencies
ukb_df$AF = as.numeric(ukb_df$AC) / as.numeric(ukb_df$AN)


#--------------------#
# ANNOTATE VARIANTS  #
#--------------------#

# 1. Annotate gene names
ukb_df = annotate_gene_names(ukb_df, gencode_df)
n_genes = length(unique(ukb_df$gene_name))

write.csv(ukb_df, str_interp('./processed/ukb_variants_all_w_gene_annotations_chr${chr}.csv'))

summary_df[nrow(summary_df)+1, ] = c("Variants in genes without 5' UTRs",
                                     n_before, nrow(ukb_df), n_genes)
n_before = nrow(ukb_df)

# 2. Remove variants that fall in multiple genes
ukb_df = remove_vars_in_two_genes(ukb_df)
n_genes = length(unique(ukb_df$gene_name))

summary_df[nrow(summary_df)+1, ] = c("Variants overlapping two gene 5' UTRs",
                                     n_before, nrow(ukb_df), n_genes)
n_before = nrow(ukb_df)

# 3. Adjust ref/alt for negatively stranded variants (put them all on + strand)
ukb_df = complementary_ref_alt_neg_strand(ukb_df, ref_5utrs_df)

# 4. Determine position within reporter sequence
ukb_df = compute_position_within_reporter_seq(ukb_df, gencode_df)


#--------------------------------#
# ANNOTATE 5' UTR VARIANT TYPES  #
#--------------------------------#

# 5. Remove variants further than 250bp
ukb_df = remove_vars_further_than_limit(ukb_df, 250)
n_genes = length(unique(ukb_df$gene_name))

summary_df[nrow(summary_df)+1, ] = c("Variants further than 250bp sequence",
                                     n_before, nrow(ukb_df), n_genes)
n_before = nrow(ukb_df)


# 6. Construct alternative reporter sequences
ukb_df = construst_alt_sequences(ukb_df, ref_5utrs_df)
n_genes = length(unique(ukb_df$gene_name))

summary_df[nrow(summary_df)+1, ] = c("Mismatch between UKB REF and GENCODE REF",
                                     n_before, nrow(ukb_df), n_genes)
n_before = nrow(ukb_df)


# 7. Remove UTRs < 25nt
ukb_df = ukb_df %>%
  mutate(utr_length = nchar(UTR5_sequence_250)) %>%
  filter(utr_length >= 25)

n_genes = length(unique(ukb_df$gene_name))

summary_df[nrow(summary_df)+1, ] = c("Variants in UTRs shorter than 25bp",
                                     n_before, nrow(ukb_df), n_genes)


#----------------------------------#
# MERGE WITH OTHER GENIC FEATURES  #
#----------------------------------#

# 8. Merge with other genic features [hs, length, OMIM, etc.]
ukb_df = merge(ukb_df, df_gene_info, by='gene_name', all.x=TRUE, all.y=FALSE)

write.csv(ukb_df, str_interp('./processed/UKB_variants_chr${chr}.csv'),
          row.names=FALSE)

write.csv(summary_df, str_interp('./processed/filtering_sumstats_chr${chr}.csv'),
          row.names=FALSE)
