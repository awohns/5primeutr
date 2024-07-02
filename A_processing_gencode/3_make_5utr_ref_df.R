rm(list=ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)


#------------#
# VARIABLES  #
#------------#
bp_limits = c(180)

stop_codons = c('TAA', 'TAG', 'TGA')


#-----------#
# FUNCTIONS #
#-----------#

# uORF annotation (& composite functions)
classify_orfs = function(seq, start_codon){

  main_orf_position = determine_main_orf_position(seq, start_codon)+1

  # Determine if AUG exists in sequence
  aug_positions = unlist(gregexpr("ATG", seq))

  orf_types = c()
  skip = c()

  for (start_pos in aug_positions){

    if (!(start_pos %in% skip)){

      # If position is main start codon
      if (start_pos == main_orf_position){
        main_ORF = str_interp("${main_orf_position}: 'main ORF'")
        orf_types = c(orf_types, main_ORF)

        # Otherwise, annotate u/oORF
      } else {

        # In or out of frame?
        frame = determine_frame(seq, start_pos, main_orf_position)

        # Upstream or overlapping
        stream = determine_stream(seq, start_pos, main_orf_position)[1]
        stop = as.numeric(determine_stream(seq, start_pos, main_orf_position)[2])

        # Remove all the AUG positions between start and stop
        skip = aug_positions[which(aug_positions %in% seq(start_pos, stop))]

        orf_types = c(orf_types, str_interp("${start_pos}: '${frame} ${stream}ORF'"))
      }

    }

  }

  orf_types_txt = str_interp("{${paste(orf_types, collapse=', ')}}")

  return(orf_types_txt)

}

determine_main_orf_position = function(seq, start_codon){

  # Find the last occurrence of start codon
  last_atg_index = max(gregexpr(start_codon, seq)[[1]])-1
  main_orf_pos = last_atg_index

  return(main_orf_pos)
}

determine_frame = function(sequence, index, main_orf_pos){

  in_frame_indices = unique(c(seq(main_orf_pos, 1, by=-3),
                              seq(main_orf_pos, nchar(sequence), by=3)))

  if (index %in% in_frame_indices){
    frame = 'inframe'
  } else {
    frame = 'outframe'
  }

  return(frame)

}

determine_stream = function(sequence, pos, main_orf_pos){

  # List all the triplets from start codon till end
  triplets = c()

  for (i in 1:floor((nchar(sequence) - pos) / 3)) {
    start = pos + (i - 1) * 3 + 1
    end = start + 2
    triplets = c(triplets, substr(sequence, start, end))
  }

  # Check if there is a STOP codon among those triplets
  found_stop_codon = any(triplets %in% stop_codons)

  # If so, check if it's before or after main ORF
  if (found_stop_codon){
    n_triplet_stop = which(sapply(triplets, `%in%`, stop_codons))[1]
    pos_stop = as.numeric((pos+1)+3*(n_triplet_stop-1))

    # If before, ORF is upstream
    if (pos_stop-3 < main_orf_pos){
      stream = 'u'

      # If after, ORF is overlapping
    } else if (pos_stop-3 >= main_orf_pos){
      stream = 'o'

      # Error
    } else {
      print("Check what's up!")
    }

    # If there is no stop codon, it must come after, therefore ORF is overlapping
  } else {
    pos_stop = nchar(sequence)-2
    stream = 'o'
  }

  return(c(stream, pos_stop))

}


# Kozak annotation
determine_kozak_strength = function(seq, start_codon){

  main_orf_position = determine_main_orf_position(seq, start_codon)

  kozak_seq = substr(seq, main_orf_position-2, main_orf_position+4)

  # Kozak strength classifications (Whiffin et al 2021):
  #   - Strong Kozak - [A/G]NNAUGG
  #   - Moderate Kozak - NNNAUGG or [A/G]NNAUGN
  #   - Weak Kozak - NNNAUGN

  if (substr(kozak_seq, 7, 7) == 'G'){

    if (substr(kozak_seq, 1, 1) %in% c('A', 'G')){
      kozak_strength = 'Strong'

    } else {
      kozak_strength = 'Moderate'
    }
  } else {

    if (substr(kozak_seq, 1, 1) %in% c('A', 'G')){
      kozak_strength = 'Moderate'
    } else {
      kozak_strength = 'Weak'
    }
  }

  return(kozak_strength)
}


#------------#
# LOAD DATA  #
#------------#
gencode_df = data.frame(fread('../../data/processed/gencode_v45_protein_coding_MANE_clean.csv'))
gencode2_df = data.frame(fread('../../data//processed/gencode2_v45_protein_coding_MANE_clean.csv'))

gencode_transcripts_df =
  data.table(fread('../../data/processed/gencode_v45_MANE_transcripts_5utr_sequences.csv'))


#------------------#
# EXTRACT 5' UTRs  #
#------------------#
df_five_prime_UTR = gencode2_df %>%
  filter(Element == 'five_prime_UTR') %>%
  select(-Element)

# Merge gencode data and sequences
df_five_prime_UTR = merge(df_five_prime_UTR, gencode_transcripts_df,
                          by='transcript_name')

# Construct sequences with length limits
for (bp_limit in bp_limits){

  limit = bp_limit-6

  colname = str_interp("UTR5_sequence_${limit}")

  length_w_limit = ifelse(df_five_prime_UTR$Length > limit,
                          limit, df_five_prime_UTR$Length)

  df_five_prime_UTR = df_five_prime_UTR %>%
    mutate(length_w_limit = ifelse(df_five_prime_UTR$Length > limit,
                                   limit, df_five_prime_UTR$Length)) %>%
    mutate(start = Length - length_w_limit + 1)

  colname = str_interp("UTR5_sequence_${bp_limit}")
  df_five_prime_UTR[,colname] = substr(df_five_prime_UTR$UTR5_sequence_and_6bp,
                                       df_five_prime_UTR$start,
                                       df_five_prime_UTR$Length+6)

  df_five_prime_UTR = df_five_prime_UTR %>%
    select(-length_w_limit, -start)

}


#--------------------#
# NOTE START CODONS  #
#--------------------#
df_start_codons = gencode_transcripts_df %>%
  mutate(length = nchar(UTR5_sequence_and_6bp),
         start = length-5,
         end = length-3) %>%
  mutate(start_codon = substr(UTR5_sequence_and_6bp, start, end)) %>%
  select(transcript_name, start_codon)

write.csv(df_start_codons,
          '../../data/processed/start_codons_gencode_v45.csv', row.names=FALSE)


#----------------------------#
# ANNOTATE 5' UTR SEQUENCES  #
#----------------------------#
df = merge(df_five_prime_UTR, df_start_codons, by='transcript_name')

df_annotated = df %>%
  rowwise() %>%
  mutate(orf_ann_full = classify_orfs(UTR5_sequence_and_6bp, start_codon),
         orf_ann_180bp = classify_orfs(UTR5_sequence_180, start_codon),
         kozak_strength = determine_kozak_strength(UTR5_sequence_and_6bp, start_codon)) %>%
  ungroup()


#------------#
# SAVE DATA  #
#------------#
write.csv(df_five_prime_UTR,
          '../../data/processed/5utr_ref_sequences_180bp_gencode_v45.csv', row.names=FALSE)

write.csv(df_annotated,
          '../../data/processed/5utr_ref_sequences_180bp_gencode_v45_annotated.csv', row.names=FALSE)
