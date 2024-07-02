rm(list=ls())

#install.packages(c('data.table', 'dplyr'))
library(data.table)
library(dplyr)
library(stringr)

#------------#
# LOAD DATA  #
#------------#
gencode_df = fread('../../data/raw/gencode.v45.annotation.gff3.no.comments',
                   header=FALSE)

names(gencode_df) = c('Chr', 'Annotation', 'Element', 'Start', 'End',
                      'idk', 'Strand', 'idk2', 'Details')

#------------#
# FUNCTIONS  #
#------------#

# GENCODE file processing
split_details_col = function(gencode_df){
  
  details_column = gencode_df$Details
  
  # Split the data by semicolon
  split_data <- strsplit(details_column, ";")
  
  # Extract unique keys
  all_keys <- unique(unlist(lapply(split_data, function(x) sapply(strsplit(x, "="), `[`, 1))))
  
  # Create a matrix to store values
  values_matrix <- matrix(NA, ncol = length(all_keys), nrow = length(split_data))
  
  # Fill in the matrix with values
  for (i in seq_along(split_data)) {
    entry_pairs <- strsplit(split_data[[i]], "=")
    entry_keys <- sapply(entry_pairs, `[`, 1)
    entry_values <- sapply(entry_pairs, `[`, 2)
    values_matrix[i, match(entry_keys, all_keys)] <- entry_values
  }
  
  # Create a data frame
  details_df <- as.data.frame(values_matrix, stringsAsFactors = FALSE)
  
  # Assign column names
  colnames(details_df) <- all_keys
  
  gencode_df = cbind(gencode_df %>%
                       select(Chr:idk2), details_df)
  
  return(gencode_df)
}

# Data wrangling
put_elements_on_exonic_position = function(gencode_df){

  colnames = c(names(gencode_df), 'Adjusted_start', 'Adjusted_end')

  gencode2_df = data.frame(matrix(nrow=0, ncol=length(colnames)))
  names(gencode2_df) = colnames

  for (gene in unique(gencode_df$gene_name)){
    
    print(gene)

    df_gene = gencode_df %>%
      filter(gene_name == gene)

    transcripts = as.character(unique(df_gene$transcript_name))

    for (transcript in transcripts){

      df_gene_transcript = df_gene %>%
        filter(transcript_name == as.character(transcript))

      main_gene_transcript_df = compute_adjusted_start_and_end_exonic(df_gene_transcript)

      ## FOR EACH CATEGORY — CHECK IF CONTIGUOUS -------------------
      categories = unique(main_gene_transcript_df$Element)
      categories = categories[-(which(categories == 'exon'))]

      for (cat in categories){

        gene_cat_df = main_gene_transcript_df %>%
          filter(Element == cat) %>%
          filter(transcript_name == transcript)

        ## If there are multiple rows for the category
        if (nrow(gene_cat_df) != 1){

          contiguous = check_if_contiguous(gene_cat_df, cat)

          ## Record with adjusted start/end, and calculate length 
          if (contiguous){
            gencode2_df = record_contiguous_entires(gene_cat_df, gencode2_df)

          } else {
            print("These elements are discontiguous. This shouldn't happen!")
          }

        } else {

          ## Record the element as such, and calculate length
          gencode2_df = record_element_as_such(gene_cat_df, gencode2_df)

        }
      }
    }
  }

  return(gencode2_df)

}

check_if_contiguous = function(gene_cat_df, cat){

  gene_cat_df = gene_cat_df[order(gene_cat_df$Adjusted_start),]

  contiguous = TRUE

  for (j in 2:nrow(gene_cat_df)){
    previous_end = gene_cat_df[j-1, 'Adjusted_end']
    start = gene_cat_df[j,'Adjusted_start']

    if (start != previous_end+1){
      contiguous = FALSE
    }
  }

  return(contiguous)
}

compute_adjusted_start_and_end_exonic = function(df){

  df_gene_transcript = df %>%
    filter(Element != 'transcript')

  exons_df = df_gene_transcript %>%
    filter(Element == 'exon')

  # Order them by smallest start (if + strand) or largest end (if - strand)
  exons_df = exons_df[order(exons_df$Start),]
  offset = 0

  for (i in 1:nrow(exons_df)){

    start = as.numeric(exons_df[i,'Start'])
    end = as.numeric(exons_df[i,'End'])
    length = abs(start-end)

    elements_df = df_gene_transcript %>%
      filter(Start %in% seq(start, end))

    elements_df$Adjusted_start = elements_df$Start - start + offset
    elements_df$Adjusted_end = elements_df$End - start + offset

    if (i==1){
      main_gene_df = elements_df
    } else {
      main_gene_df = rbind(main_gene_df, elements_df)
    }

    if (i == 1){
      offset = length+1
    } else {
      offset = offset + length+1
    }

  }

  return(main_gene_df)

}

record_contiguous_entires = function(gene_cat_df, gencode2_df){

  # Note that this is independent of strand as it records it on the 5'>3' logic
  start = min(gene_cat_df$Start)
  end = max(gene_cat_df$End)

  # Record data by adjusting Start/End
  gene_cat_df$Length = abs(gene_cat_df$Adjusted_end-gene_cat_df$Adjusted_start)+1
  element_length = sum(gene_cat_df$Length)
  
  df = gene_cat_df %>%
    select(-Start, -End, -Length) %>%
    mutate(Start = start,
           End = end,
           Adjusted_start = min(Adjusted_start),
           Adjusted_end = max(Adjusted_end),
           Length = element_length) %>%
    distinct()

  gencode2_df = rbind(gencode2_df, df)

  return(gencode2_df)

}

record_element_as_such = function(gene_cat_df, gencode2_df){
  
  # Calculate length
  length = length(seq(unique(gene_cat_df$Adjusted_start), 
                      unique(gene_cat_df$Adjusted_end)))
  
  df = gene_cat_df %>%
    mutate(Length = length) %>%
    distinct()
  
  gencode2_df = rbind(gencode2_df, df)
  
  return(gencode2_df)
  
}


#-----------------------------#
# SPLITTING "DETAILS" column  #
#-----------------------------#
gencode_df = split_details_col(gencode_df)


#----------------------------#
# EXTRACT STRANDEDNESS INFO  #
#----------------------------#
df = gencode_df %>%
  filter(gene_type == 'protein_coding') %>%
  select(hgnc_id, gene_name, Strand) %>%
  distinct()

write.csv(df, '../../data/processed/gencode_v45_genes_and_strandedness.csv', row.names=FALSE)


#--------------#
# FILTER DATA  #
#--------------#
gencode_df_filtered = gencode_df %>%
  filter(gene_type == 'protein_coding') %>%
  filter(grepl('MANE_Select', tag))

gencode_final_df = gencode_df_filtered %>%
  select(gene_name, Element, Chr,
         Start, End, Strand, transcript_name)

ann_df = gencode_df %>%
  select(hgnc_id, gene_name, gene_id, gene_type,
         transcript_name, transcript_id, transcript_type) %>%
  distinct()


#------------#
# SAVE DATA  #
#------------#
write.csv(gencode_final_df,
          '../../data/processed/gencode_v45_protein_coding_MANE_clean.csv', row.names=FALSE)

write.csv(ann_df,
          '../../data/processed/ann_gencode_v45_protein_coding_MANE_clean.csv', row.names=FALSE)


#------------------------------#
# PROCESSING TO REMOVE INTRONS #
#------------------------------#
gencode_df = gencode_final_df

gencode2_df = put_elements_on_exonic_position(gencode_df)


#--------------------------------------#
# COUNT NUMBER OF UTRs AND ADD COLUMN  #
#--------------------------------------#

# Number of different elements per gene
n_utrs_df = gencode2_df %>%
  select(gene_name, Element) %>%
  group_by(gene_name, Element) %>%
  summarize(n_element_per_gene=n()) %>%
  distinct() %>%
  spread(Element, n_element_per_gene) %>%
  select(gene_name, five_prime_UTR) %>%
  rename(n_five_prime_UTR = 'five_prime_UTR')

n_utrs_df[is.na(n_utrs_df)] = 0

gencode2_df = merge(gencode2_df, n_utrs_df, by='gene_name')

write.csv(gencode2_df, "../../data/processed/gencode2_v45_protein_coding_MANE_clean.csv",
          row.names=FALSE)


#--------------------------------------------------------------------#
# LIST OF MANE SELECT TRANSCRIPTS FOR GENCODE TRANSCRIPT PROCESSING  #
#--------------------------------------------------------------------#
mane_select_transcripts = unique(gencode_df$transcript_name)

write.csv(mane_select_transcripts, 
          '../../data/processed/MANE_select_transcript_names_GENCODE_v45.csv',
          row.names=FALSE)


#-------------------#
# PRINT SOME STATS  #
#-------------------#
n_rows = nrow(gencode_df)
n_rows_post_filtering = nrow(gencode_df_filtered)
p_post_filtering = round(n_rows_post_filtering/n_rows, 3)*100

print(str_interp("There were ${n_rows} before filtering."))
print(str_interp("There were ${n_rows_post_filtering} after keeping protein-coding genes and MANE_select transcripts."))
print(str_interp("That's ${p_post_filtering}%."))
