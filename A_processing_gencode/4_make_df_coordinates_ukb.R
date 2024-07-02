rm(list=ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)


#------------#
# VARIABLES  #
#------------#
bp_limits = c(250)


#------------#
# FUNCTIONS  #
#------------#
get_positions_for_ukb_extraction = function(gencode_df, gencode2_df, limit){

  # If the 5' UTR is encoded in one exon (limit=200 as illustrative):
  #     (a) if exon < 200 --> record it
  #     (b) if exon > 200 --> cut it
  #
  # If the 5' UTR is encoded in multiple exons:
  #     (a) LOOP, and cut the last.

  limit = limit-6

  df_coordinates = data.frame(matrix(nrow=0, ncol=8))
  names(df_coordinates) = c(names(gencode_df), 'exon_length')

  for (gene in unique(gencode_df$gene_name)){

    # Subset GENCODE exons for gene
    gen_df = gencode_df %>%
      filter(gene_name == gene) %>%
      filter(Element == 'five_prime_UTR') %>%
      mutate(exon_length = abs(End-Start)+1) %>%
      arrange(Start)

    # Subset GENCODE2 element
    gen2_df = gencode2_df %>%
      filter(gene_name == gene) %>%
      filter(Element == 'five_prime_UTR')

    ## RECORD EXONS AS COORDINATES -----------------------------------------------

    # If the gene has a 5' UTR:
    if (nrow(gen_df) != 0 & nrow(gen2_df) != 0){

      utr_length = unique(gen2_df$Length)

      len_limit = ifelse(utr_length > limit, limit, utr_length)

      # IF 5' UTR is encoded in one exon ----------------------
      if (nrow(gen_df) == 1){
        exon_len = gen_df$exon_length

        # If exon > LIMIT, cut it (adjust end)
        if (exon_len > len_limit){
          gen_df$End = gen_df$Start+(len_limit-1)
          gen_df$exon_length = abs(gen_df$End-gen_df$Start)+1
        }

        # Record exon in coordinates
        df_coordinates = rbind(df_coordinates, gen_df)

        # IF 5' UTR is encoded in several exons ---------------
      } else {

        # Cut the exon in loops
        len_left = len_limit

        # Loop through exons
        for (i in 1:nrow(gen_df)){

          # If there is still length of the UTR left
          if (len_left != 0){

            exon_length = gen_df[i,'exon_length']
            exon_start = gen_df[i,'Start']
            exon_end = gen_df[i,'End']

            # If exon if shorter than length left, record it.
            if (exon_length < len_left){

              # Subtract exon length from length left
              len_left = len_left - exon_length

              # If exon is longer than length left, CUT (update end)
            } else {

              # Record new end, and drop subsequent exons
              gen_df[i,'End'] = exon_start+(len_left-1)
              gen_df[i,'exon_length'] = abs(gen_df[i,'End'] - gen_df[i,'Start']) + 1

              gen_df = gen_df[1:i, ]

              # Update length left
              len_left = 0
            }

            # Add exon to coordinates
            df_coordinates = rbind(df_coordinates, gen_df[i,])

          }
        }
      }

      # Check that recorded exons sum up to expected length
      df_coordinates_gene = df_coordinates %>%
        filter(gene_name == gene)

      total_exonic_bp = sum(df_coordinates_gene$exon_length)

      if (len_limit != total_exonic_bp){
        print(gene)
        print("Ayy you messed it up!")

        if (total_exonic_bp > len_limit){
          print("More exonic bp's than allowed")
        } else {
          print("Less exonic bp's than allowed")
        }
      }

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
      df_coordinates = rbind(df_coordinates, cds_df = gencode_df %>%
                               filter(gene_name == gene) %>%
                               filter(Element == 'start_codon') %>%
                               mutate(exon_length = abs(End-Start)+1))

    } else {
      print(str_interp("${gene} has no 5' UTRs! "))
    }

  }

  return(df_coordinates)

}


#------------#
# LOAD DATA  #
#------------#
gencode_df = data.frame(fread('../../data/processed/gencode_v45_protein_coding_MANE_clean.csv'))
gencode2_df = data.frame(fread('../../data/processed/gencode2_v45_protein_coding_MANE_clean.csv'))


#-------------------------------------------#
# MAKE DF_COORDINATES FOR DIFFERENT LIMITS  #
#-------------------------------------------#
for (limit in bp_limits){

  df_coordinates = get_positions_for_ukb_extraction(gencode_df, gencode2_df, limit=limit)

  df_coordinates = df_coordinates %>%
    filter(Element == "five_prime_UTR")
  
  write.csv(df_coordinates,
            str_interp('../../data/processed/5utrs_hg38_coordinates_limit_${limit}bp.csv'),
            row.names=FALSE)

	print("Done! File saved!")

}
