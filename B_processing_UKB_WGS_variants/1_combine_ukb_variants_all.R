rm(list=ls())

library(dplyr)
library(stringr)
library(data.table)

#-----------#
# FUNCTIONS #
#-----------#
extract_chr = function(filename){
  
  parts = unlist(strsplit(filename, "_"))
  
  # Remove b, .txt and chr from parts
  parts = str_replace(parts, "chr", "")
  chr = parts[3]
  
  return(chr)
  
}

split_multiallelic_vars_into_rows_chr = function(df){
  
  multiallelic_df = df %>%
    filter(grepl(",", ALT))
  
  biallelic_df = df %>%
    filter(!grepl(",", ALT))
  
  final_df = data.frame(matrix(nrow=0, ncol=ncol(df)))
  names(final_df) = names(df)
  
  for (i in 1:nrow(multiallelic_df)){
    
    row = multiallelic_df[i,]
    
    alt <- unlist(strsplit(as.character(row$ALT), ","))
    ac <- unlist(strsplit(as.character(row$AC), ","))
    
    new_rows = data.frame(matrix(nrow=length(alt), ncol=ncol(final_df)))
    names(new_rows) = names(df)
    
    new_rows$POS = row$POS
    new_rows$REF = row$REF
    new_rows$AN = row$AN
    new_rows$CHR = row$CHR
    
    new_rows$ALT = alt
    new_rows$AC = ac
    
    final_df = rbind(final_df, new_rows)
    
  }
  
  df = rbind(biallelic_df, final_df)
  df = df[order(df$POS),]
  
  return(df)
  
}


#-----------#
# LOAD DATA #
#-----------#
setwd('/home/users/tami/5utr_extended_panel/data/ukb_variants')

files = list.files(".")
files = files[grepl(".txt", files)]


#--------------#
# PROCESS DATA #
#--------------#
df = data.frame(matrix(nrow=0, ncol=6))
names(df) = c('CHR', 'POS', 'REF', 'ALT', 'AC', 'AN')

for (file in files){
  
  # Load data
  tmp_df = fread(file)
  
  if (nrow(tmp_df) != 0){
    
    # Add chromosome
    tmp_df = cbind(extract_chr(file), tmp_df)
    names(tmp_df) = c('CHR', 'POS', 'REF', 'ALT', 'AC', 'AN')
    
    # Add to main df
    df = rbind(df, tmp_df)
  }
}

#-----------------------------#
# SEPARATE VARIANTS INTO ROWS  #
#------------------------------#
df = split_multiallelic_vars_into_rows_chr(df)

n_vars = nrow(df)
print(n_vars)

write.csv(df, "/home/users/tami/5utr_extended_panel/data/processed/ukb_variants_all.csv", row.names=FALSE)
