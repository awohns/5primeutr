rm(list=ls())

library(dplyr)
library(tidyr)
library(stringr)


#------------#
# FUNCTIONS  #
#------------#
get_positions_for_ukb_extraction_6bp = function(gencode_df){
  
  df_coordinates = data.frame(matrix(nrow=0, ncol=8))
  names(df_coordinates) = c(names(gencode_df), 'exon_length')
  
  for (gene in unique(gencode_df$gene_name)){
    
    ## Add the first 6 bp of CDS [NOTE: CDS includes the start codon.]
    cds_df = gencode_df %>%
      filter(gene_name == gene) %>%
      filter(Element == 'CDS') %>%
      mutate(exon_length = abs(End-Start)+1)
    
    strand = unique(cds_df$Strand)
    
    cds_6bp_df = data.frame(matrix(nrow=0, ncol=ncol(cds_df)))
    names(cds_6bp_df) = names(cds_df)
    
    len_left = 6
    
    if (strand == "+"){
      cds_df = cds_df %>%
        arrange(Start)
      
      for (j in 1:nrow(cds_df)){
        if (len_left != 0){
          
          exon_len = cds_df[j,'exon_length']
          
          if (exon_len >= len_left){
            cds_df[j,'End'] = cds_df[j,'Start']+len_left-1
            cds_6bp_df = rbind(cds_6bp_df, cds_df[j,])
            len_left = 0
            
          } else {
            cds_6bp_df = rbind(cds_6bp_df, cds_df[j,])
            len_left = len_left-exon_len
          }
        }
      }
      
    } else if (strand == "-"){
      
      cds_df = cds_df %>%
        arrange(-End)
      
      for (j in 1:nrow(cds_df)){
        if (len_left != 0){
          
          exon_len = cds_df[j,'exon_length']
          
          if (exon_len >= len_left){
            cds_df[j,'Start'] = cds_df[j,'End']-len_left+1
            cds_6bp_df = rbind(cds_6bp_df, cds_df[j,])
            len_left = 0
            
          } else {
            cds_6bp_df = rbind(cds_6bp_df, cds_df[j,])
            len_left = len_left-exon_len
          }
        }
      }
    }
    
    cds_6bp_df$exon_length = abs(cds_6bp_df$End - cds_6bp_df$Start)+1
    df_coordinates = rbind(df_coordinates, cds_6bp_df)
    
  }
  
  return(df_coordinates)
  
}


#------------#
# LOAD DATA  #
#------------#
gencode_df = data.frame(fread('./raw/gencode_v45_protein_coding_MANE_clean.csv'))

df_coordinated_6bp = get_positions_for_ukb_extraction_6bp(gencode_df)


#------------#
# SAVE DATA  #
#------------#
write.csv(df_coordinated_6bp, 
          '../../data/processed/5utrs_hg38_coordinates_6bp.csv',
          row.names=FALSE)