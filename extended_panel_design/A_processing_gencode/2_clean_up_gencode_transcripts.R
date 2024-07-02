rm(list=ls())

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)


#-----------#
# LOAD DATA #
#-----------#
df_raw = data.frame(fread('../../data/raw/gencode.v45.pc_transcripts.fa',
                          header=FALSE))
mane_transcripts = 
  read.csv('../../data/processed/MANE_select_transcript_names_GENCODE_v45.csv')


# Extract 5' UTR sequences from GENCODE file -------------------------------
start_lines = grep("ENST", df_raw$V1)

df = data.frame(matrix(nrow=0, ncol=2))
names(df) = c('transcript_name', 'UTR5_sequence_and_6bp')

for (i in 1:length(start_lines)){
  
  # Define start and end row for sequence 
  info_line = start_lines[i]
  seq_start = start_lines[i]+1
  seq_end = start_lines[i+1]-1
  
  # If this is the last entry in data, sequence ends at last row
  if (i == length(start_lines)){
    seq_end = nrow(df_raw)
  }
  
  # Extract important info (transcript name and UTR positions)
  info = df_raw[info_line,1]
  info_list = unlist(str_split(str_replace(info, ">", ""), "\\|"))
  
  transcript_name = info_list[5]
  
  # If transcript is a MANE_Select
  if (transcript_name %in% mane_transcripts$x){
    
    utr_positions = info_list[grep("UTR5", info_list)]
    
    # If transcript has a 5' UTR
    if (length(utr_positions) != 0){
      
      # Combine sequence across rows into a string 
      sequence = paste(df_raw[seq_start:seq_end,1], collapse = '')
      
      # Obtain 5' UTR relative positions in sequence   
      utr5_positions = as.numeric(unlist(str_split(
        unlist(str_split(utr_positions, ":"))[2], "-")))
      
      start_utr = utr5_positions[1]
      end_utr = utr5_positions[2]+6   # This covers first 6bp from CDS
      
      # Subset 5' UTR sequence from string 
      sequence_5utr = substr(sequence, start_utr, end_utr)
      
      # Add to reformatted df 
      df[nrow(df)+1, ] = c(transcript_name, sequence_5utr)
      
    } else {
      df[nrow(df)+1, ] = c(transcript_name, NA)
    }
  }
  
}

write.csv(df, '../../data/processed/gencode_v45_MANE_transcripts_5utr_sequences.csv', 
          row.names=FALSE)



