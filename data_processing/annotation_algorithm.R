# ANNOTATE PREDICTED EFFECTS — REVISED

## Annotate data for predicted positive, predicted negative, or null effect

rm(list=ls())

library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(stringr)

setwd('~/Desktop/5prime_utr')

##################################################################################################

#-----------#
# VARIABLES #
#-----------#
stop_codons = c('TAA', 'TAG', 'TGA')

orf_types = c('main_ORF', 'inframe_oORF', 'inframe_uORF',
              'outframe_oORF', 'outframe_uORF', 'stop_first_codon')

# Combinations of variant types that act like replacements
replacements_df = data.frame(
  
  rbind(c('Loss of main ORF, Gain of inframe oORF',
          'Loss of main ORF',
          'lost_main_ORF'),
        
        c('Loss of main ORF, Gain of outframe oORF',
          'Loss of main ORF',
          'lost_main_ORF'),
        
        c('Loss of inframe oORF, Gain of main ORF',
          'Gain of main ORF',
          'gained_main_ORF'),
        
        c('Loss of outframe oORF, Gain of main ORF',
          'Gain of main ORF',
          'gained_main_ORF'),
        
        c('Change of position for outframe uORF, Gain of outframe uORF',
          'Gain of outframe uORF',
          'gained_outframe_uORF'),
        
        c('Change of position for inframe uORF, Gain of inframe uORF',
          'Gain of inframe uORF',
          'gained_inframe_uORF'),
        
        c('Change of position for outframe oORF, Gain of outframe oORF', 
          'Gain of outframe oORF',
          'gained_outframe_oORF'),
        
        c('Change of position for inframe oORF, Gain of inframe oORF',
          'Gain of inframe oORF',
          'gain_inframe_oORF'))
)

names(replacements_df) = c('pattern_txt', 'drop_txt', 'col_drop')


# Allocation of variant types to positive/negative/null
pos_var_types = c("lost_inframe_uORF", "lost_outframe_uORF", 
                  "lost_outframe_oORF", "gained_main_ORF")

neg_var_types = c("lost_main_ORF", "gained_inframe_uORF", 
                  "gained_outframe_oORF", "gained_outframe_uORF")

null_var_types = c("gained_inframe_oORF", "changed_pos_main_ORF",
                   "changed_pos_inframe_uORF", "changed_pos_outframe_uORF",
                   "changed_pos_outframe_oORF", "changed_pos_inframe_oORF", 
                   "lost_inframe_oORF")

total_loss_var_types = c('Complete loss of ORFs', 
                         'Loss of main ORF', 
                         'Gain of outframe oORF', 
                         'Gain of stop first codon')

##################################################################################################

#------------#
# FUNCTIONS  #
#------------#

# Data cleaning (may need adjustments based on raw data)
transform_to_df__gene_humvar_ref_insertseq = function(df_raw, df_sequences){
  
  # Remove spike-ins 
  df_raw = df_raw %>%
    rename(humvar='X') %>%
    filter(!(grepl("spk", humvar))) %>%
    filter(!(grepl("NFASC", humvar)))
  
  # Replace U's with T's 
  df_sequences$insert_seq = str_replace_all(df_sequences$insert_seq, "U", "T")
  
  # Merge df_raw and df_sequences 
  df = merge(df_raw %>%
               select(humvar),
             df_sequences %>%
               select(-X), by='humvar', all.x=TRUE, all.y=FALSE)
  
  # Clean up raw df 
  df = df %>%
    separate(humvar, into = c('gene_id', 'humvar'), sep='_') %>%
    select(gene_id, humvar, insert_seq)
  
  # Split references and variants 
  df_refs = df %>%
    filter(grepl("ref", humvar)) %>%
    rename(ref_seq = 'insert_seq') %>%
    select(gene_id, ref_seq)
  
  genes = df_refs$gene_id
  
  df_vars = df %>%
    filter(!grepl("ref", humvar))
  
  n_vars_before = nrow(df_vars)
  
  df = merge(df_vars, df_refs, by='gene_id')
  n_vars_after = nrow(df)
  
  print(str_interp("Before: ${n_vars_before} variants"))
  print(str_interp("After removing variants w/o ref: ${n_vars_after} variants"))
  
  df = df %>%
    select(gene_id, humvar, ref_seq, insert_seq)
  
  return(df)
  
}

merge_df_with_start_codons = function(df, df_gencode){
  
  match = intersect(df_gencode$gene_name, df$gene_id)
  no_match = df$gene_id[!(df$gene_id %in% match)]
  
  # Merge data with perfect gene matches 
  match_df = merge(df, df_gencode, by.x='gene_id', by.y='gene_name')
  
  # Merge data with imperfect gene matches 
  isoforms_match_df = df %>% 
    filter(gene_id %in% no_match) %>%
    mutate(gene_id_to_match = str_replace_all(gene_id, 'iso2', ''))
  
  isoforms_match_df = merge(isoforms_match_df, df_gencode, 
                            by.x='gene_id_to_match', by.y='gene_name',
                            all.x=TRUE, all.y=FALSE)
  
  # Genes without start codons 
  genes_no_start_codon = unique(isoforms_match_df[
    is.na(isoforms_match_df$start_codon), 'gene_id'])
  
  # Add start codon ATG 
  isoforms_match_df$start_codon[is.na(isoforms_match_df$start_codon)] = 'ATG'
  
  print(str_interp("For genes ${paste(genes_no_start_codon, collapse=', ')}"))
  print("start codon was set to ATG (no start codon detected in GENCODE file).")
  
  # Merge the two parts of the data
  df = rbind(match_df, isoforms_match_df %>%
               select(-gene_id_to_match))
  
  return(df)
  
}


# uORF annotation (& composite functions)
classify_orfs = function(seq, start_codon){
  
  # Determine main ORF position
  main_orf_position = determine_main_orf_position(seq, start_codon)+1
  
  if (is.na(main_orf_position)){
    
    orf_types_txt = "{}"
    
  } else {
    
    # Determine if AUG exists in sequence
    aug_positions_ATG = unlist(gregexpr("ATG", seq))
    
    # If there are AUGs in the sequence, process them:
    if (!(length(aug_positions_ATG) == 1 && aug_positions_ATG == -1)){
      
      orf_types = c()
      skip = c()
      
      for (start_pos in aug_positions_ATG){
        
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
            skip = aug_positions_ATG[which(aug_positions_ATG %in% seq(start_pos, stop))]
            
            orf_types = c(orf_types, str_interp("${start_pos}: '${frame} ${stream}ORF'"))
          }
          
        }
        
      }
      
    } else {
      orf_types = c()
    }
    
    # If start codon is NOT "AUG", add that
    if (start_codon != 'ATG'){
      orf_types = c(orf_types, str_interp("${main_orf_position}: 'main ORF'"))
    }
    
    # ANNOTATE premature stop in first codon 
    seq_len = nchar(seq)
    first_codon = substr(seq, seq_len-2, seq_len)
    
    if (first_codon %in% stop_codons){
      pos = seq_len-1
      orf_types = c(orf_types, str_interp("${pos}: 'stop first codon'"))
    }
    
    orf_types_txt = str_interp("{${paste(orf_types, collapse=', ')}}")
    
  }
  
  return(orf_types_txt)
  
}

classify_orfs_alt = function(seq, seq_ref, start_codon){
  
  main_orf_position_ref = determine_main_orf_position(seq_ref, start_codon)+1
  
  # Sometimes, indels shift the main ORF start position by 1.
  # In this case, use the position of the main ORF based on the insert_seq
  if (nchar(seq) != nchar(seq_ref)){
    main_orf_position_alt = determine_main_orf_position(seq, start_codon)+1
    
    if (main_orf_position_alt %in% seq(main_orf_position_ref-1, 
                                       main_orf_position_ref+1)){
      main_orf_position = main_orf_position_alt
      
      # Likely loss of ORF 
    } else {
      main_orf_position = main_orf_position_ref 
    }
    
  } else {
    main_orf_position = main_orf_position_ref 
  }
  
  if (is.na(main_orf_position)){
    
    orf_types_txt = "{}"
    
  } else {
    
    # Determine if AUG exists in sequence
    aug_positions_ATG = unlist(gregexpr("ATG", seq))
    
    # If there are AUGs in the sequence, process them:
    if (!(length(aug_positions_ATG) == 1 && aug_positions_ATG == -1)){
      
      orf_types = c()
      skip = c()
      
      for (start_pos in aug_positions_ATG){
        
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
            skip = aug_positions_ATG[which(aug_positions_ATG %in% seq(start_pos, stop))]
            
            orf_types = c(orf_types, str_interp("${start_pos}: '${frame} ${stream}ORF'"))
          }
          
        }
        
      }
      
    } else {
      orf_types = c()
    }
    
    # If start codon is NOT "AUG", add that
    if (start_codon != 'ATG'){
      orf_types = c(orf_types, str_interp("${main_orf_position}: 'main ORF'"))
    }
    
    ## ANNOTATE premature stop in first codon 
    seq_len = nchar(seq)
    first_codon = substr(seq, seq_len-2, seq_len)
    
    if (first_codon %in% stop_codons){
      
      # Reference first codon 
      ref_len = nchar(seq_ref)
      ref_first_codon = substr(seq_ref, ref_len-2, ref_len)
      
      # Compare REF and ALT first codon
      tmp_df = data.frame(
        cbind(unlist(strsplit(first_codon, "")),
              unlist(strsplit(ref_first_codon, ""))))
      
      # Find position of variant 
      pos = ref_len-3+which(tmp_df$X1 != tmp_df$X2)
      
      if (length(pos) != 1){
        pos = pos[1]
      }
      
      orf_types = c(orf_types, str_interp("${pos}: 'stop first codon'"))
    }
    
    orf_types_txt = str_interp("{${paste(orf_types, collapse=', ')}}")
    
  }
  
  return(orf_types_txt)
}

determine_main_orf_position = function(seq, start_codon){
  
  # Find the last occurrence of start codon
  matches = gregexpr(start_codon, seq)[[1]]
  last_atg_index = ifelse(all(matches == -1), NA, max(matches) - 1)
  
  main_orf_pos = last_atg_index
  
  return(main_orf_pos)
}

determine_frame = function(sequence, index, main_orf_pos){
  
  if (index < main_orf_pos){
    in_frame_indices = seq(main_orf_pos, 1, by=-3)
  } else if (index >= main_orf_pos){
    in_frame_indices = seq(main_orf_pos, nchar(sequence), by=3)
  }
  
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
    start = pos + (i - 1) * 3
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
    if (pos_stop-2 < main_orf_pos){
      stream = 'u'
      
      # If after, ORF is overlapping
    } else if (pos_stop-2 >= main_orf_pos){
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


# Annotate variant changes
annotate_variant_types = function(df_annotated, orf_types,
                                  ref_orf_colname, alt_orf_colname){
  
  # Annotate variants as "Some change" vs "No change"
  df_annotated$variant_type = 'Some change'
  df_annotated[which(df_annotated[,ref_orf_colname] == df_annotated[,alt_orf_colname]),
               'variant_type'] = 'No change'
  
  # Make empty columns for variant types 
  for (var_type in c('gained', 'lost', 'changed_pos')){
    for (orf_type in c(orf_types, "ORFs_total")){
      df_annotated[,str_interp("${var_type}_${orf_type}")] = 0
    }
  }
  
  for (i in 1:nrow(df_annotated)){
    print(i)
    
    if (df_annotated[i,'variant_type'] == 'Some change'){
      
      # Ref and alt
      ref = df_annotated[i,ref_orf_colname]
      alt = df_annotated[i, alt_orf_colname]
      
      ## Extract non-shared ORFs between ref and alt
      orfs_ref = strsplit(gsub("[{}]", "", as.character(ref)), ", ")[[1]]
      orfs_alt = strsplit(gsub("[{}]", "", as.character(alt)), ", ")[[1]]
      
      # Check if there is a shared variant in ref and alt, and if so, remove it.
      no_change_var = intersect(orfs_ref, orfs_alt)
      
      if (length(no_change_var) > 0){
        orfs_ref = orfs_ref[!(orfs_ref %in% no_change_var)]
        orfs_alt = orfs_alt[!(orfs_alt %in% no_change_var)]
      }
      
      # Make an empty variant type list
      variant_type_list = c()
      
      
      ## TYPE 1: TOTAL LOSS OF ORFs -----------------------------------------------------------------
      if (alt == "{}"){
        
        orf_change_type = 'Complete loss of ORFs'
        variant_type_list = orf_change_type
        
        if (grepl('main ORF', orfs_ref)){
          df_annotated[i,'lost_main_ORF'] = 1
        } else if (grepl('inframe oORF', orfs_ref)){
          df_annotated[i,'lost_inframe_oORF'] = 1
        } else {
          print("Reference doesn't have a main ORF or inframe oORF. Adjust annotation accordingly!")
          print(i)
        }
        
        ## TYPES 2-3: LOSS OR CHANGE OF POSITION ----------------------------------------------------
      } else {
        
        # Process reference ORFs (must be losses or change of position)
        for (orf_ref in orfs_ref){
          
          ref_type = extract_orf_types_from_txt(orf_ref)
          alt_types = extract_orf_types_from_txt(orfs_alt)
          
          ## TYPE 2: LOSS OF ORF --------------------------------------------------------------------
          if (!(ref_type %in% alt_types)){
            
            col_name = str_interp("lost_${ref_type}")
            df_annotated[i,col_name] = df_annotated[i,col_name] + 1
            
            orf_change_type = str_interp("Loss of ${ref_type}")
            
            ## TYPE 3: CHANGE OF ORF POSITION -------------------------------------------------------
          } else if (ref_type %in% alt_types){
            
            col_name = str_interp("changed_pos_${ref_type}")
            df_annotated[i,col_name] = df_annotated[i,col_name] + 1
            
            orf_change_type = str_interp("Change of position for ${ref_type}")
            
            # Remove variant from alt
            orfs_alt = remove_orf_from_orfs_alt(orfs_alt, orf_ref, ref_type)
          }
          
          # Add variant to list
          variant_type_list = c(variant_type_list, orf_change_type)
          
        }
        
        # 2. Process alt ORFs (which must be gains)
        if (length(orfs_alt) != 0){
          
          for (orf_alt in orfs_alt){
            
            alt_type = extract_orf_types_from_txt(orf_alt)
            
            ## TYPE 4: GAIN OF ORF POSITION -------------------------------------------------------
            col_name = str_interp("gained_${alt_type}")
            df_annotated[i,col_name] = df_annotated[i,col_name] + 1
            
            orf_change_type = str_interp("Gain of ${alt_type}")
            variant_type_list = c(variant_type_list, orf_change_type)
            
          }
        }
      }
      
      variant_type_txt = gsub("_", " ", paste(variant_type_list, collapse=', '))
      df_annotated[i,'variant_type'] = variant_type_txt
      
    }
  }
  
  return(df_annotated)
  
}

extract_orf_types_from_txt = function(txt){
  
  types = gsub(" ", "_",
               gsub("'", "", unlist(regmatches(txt, gregexpr("'(.*?)'", txt)))))
  
  return(types)
  
}

remove_orf_from_orfs_alt = function(orfs_alt, orf_ref, ref_type){
  
  ref_type_txt = str_replace(ref_type, "_", " ")
  
  matches = orfs_alt[grepl(ref_type_txt, orfs_alt)]
  n_matches = length(matches)
  
  # If ref_type exists in alternative types only once, match it and remove it 
  if (n_matches  == 1){
    orfs_alt = orfs_alt[!grepl(ref_type_txt, orfs_alt)]
    
  } else {
    ref_index = as.numeric(sub(":.*", "", orf_ref)) 
    
    for (match in matches){
      alt_index = as.numeric(sub(":.*", "", match))
      
      # If index is within +/- 1 of reference index, match and remove  
      if (alt_index %in% seq(ref_index-1, ref_index+1)){
        orfs_alt = orfs_alt[-(which(orfs_alt == match))]
      }
    }
  }
  
  return(orfs_alt)
}


# Replacements
adjust_replacements_in_variant_ann = function(df, replacements_df){
  
  df$variant_type_adjusted = df$variant_type 
  
  for (k in 1:nrow(replacements_df)){
    
    patterns_txt = unlist(strsplit(replacements_df[k,'pattern_txt'], ", "))
    drop_txt = replacements_df[k,'drop_txt']
    col_drop = replacements_df[k,'col_drop']
    
    for (i in 1:nrow(df)){
      variant_ann = unlist(strsplit(df[i,'variant_type'], ', '))
      
      # If the variant annotation contains a pattern for replacement, adjust it. 
      if (all(patterns_txt %in% variant_ann)){
        variant_ann_revised = variant_ann[!(variant_ann %in% drop_txt)]
        
        variant_ann_revised_txt = paste(variant_ann_revised, collapse=', ')
        
        df[i,'variant_type_adjusted'] = variant_ann_revised_txt
        df[i,col_drop] = 0
        
        # Otherwise, keep as is. 
      }
    }
    
  }
  
  return(df)
  
}


# Annotate predicted categories
annotate_pos_neg_null_loss = function(df_annotated){
  
  # Calculate # of positive variants
  df_annotated$n_pos_var_types = ifelse(df_annotated$variant_type == 'No change',
                                        NA, rowSums(df_annotated[,pos_var_types]))
  
  # Calculate # of negative variants
  df_annotated$n_neg_var_types = ifelse(df_annotated$variant_type == 'No change',
                                        NA, rowSums(df_annotated[,neg_var_types]))
  
  # Calculate # of null variants
  df_annotated$n_null_var_types = ifelse(df_annotated$variant_type == 'No change',
                                         NA, rowSums(df_annotated[,null_var_types]))
  
  for (i in 1:nrow(df_annotated)){
    
    if (df_annotated[i,'variant_type'] == 'No change'){
      df_annotated[i,'predicted_category'] = 'Null'
      
    } else {
      
      n_pos = df_annotated[i,'n_pos_var_types']
      n_neg = df_annotated[i,'n_neg_var_types']
      n_null = df_annotated[i,'n_null_var_types']
      
      if (n_pos != 0 & n_pos > n_neg){
        df_annotated[i,'predicted_category'] = 'Positive'
      } else if (n_neg !=0 & n_neg > n_pos){
        df_annotated[i,'predicted_category'] = 'Negative'
      } else if (n_pos == 0 & n_neg == 0 & n_null != 0){
        df_annotated[i,'predicted_category'] = 'Null'
      } else if (n_pos == n_neg){
        df_annotated[i,'predicted_category'] = 'Null'
      } else {
        df_annotated[i,'predicted_category'] = 'Unresolved'
      }
      
    }
  }
  
  # Annotate total losses 
  df_annotated = df_annotated %>%
    
    mutate(predicted_category = ifelse(
      variant_type_adjusted %in% 
        total_loss_var_types,
      'Total loss', predicted_category)) %>%
    
    select(-starts_with("n_"))
  
  return(df_annotated)
  
}



##################################################################################################

#----------------------#
# LOAD & PROCESS DATA  #
#----------------------#

# NaP-TRAP data  
df_raw = read.csv('./data/raw/humvar_5utr_ntrap_v6_pm_reads_counts_041824.csv')
df_sequences = read.csv('./data/raw/humvar_insert_seq_all.csv')

df = transform_to_df__gene_humvar_ref_insertseq(df_raw, df_sequences)

# Start codons
df_gencode = data.frame(fread('./data/processed/5utr_ref_sequences_gencode_v45_annotated.csv'))
df_gencode = df_gencode %>% 
  select(gene_name, start_codon) %>% distinct()

df = merge_df_with_start_codons(df, df_gencode)


#----------------------#
# SEQUENCE ANNOTATION  #
#----------------------#
df$ref_len = nchar(df$ref_seq)
df$alt_len = nchar(df$insert_seq)

for (i in 1:nrow(df)){
  
  ref_sequence = df[i, 'ref_seq']
  alt_sequence = df[i, 'insert_seq']
  start_codon = df[i, 'start_codon']
  
  df[i,'main_orf_position'] = determine_main_orf_position(ref_sequence, start_codon)
  df[i,'orf_ann_ref'] = classify_orfs(ref_sequence, start_codon)
  df[i,'orf_ann_alt'] = classify_orfs_alt(alt_sequence, ref_sequence, start_codon)
  
}


#-------------------------#
# ANNOTATE VARIANT TYPES  #
#-------------------------#
df = annotate_variant_types(df, orf_types,
                            'orf_ann_ref', 'orf_ann_alt')

df = adjust_replacements_in_variant_ann(df, replacements_df)


#--------------------------------#
# ANNOTATE PREDICTED CATEGORIES  #
#--------------------------------#
df = annotate_pos_neg_null_loss(df)


#------------#
# SAVE DATA  #
#------------#
write.csv(df, './data/processed/annotated_variants_90k_revised.csv')
