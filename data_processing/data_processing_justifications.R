rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(stringr)

library(plotrix)
library(devtools)
library(ashr)
library(viridis)
library(DescTools)

setwd('~/Desktop/5prime_utr/')


#-----------------#
# SPECIFICATIONS  #
#-----------------#
read_filter_spec = 0
epsilon_spec = 1
normalization_spec = 'rpm'
reference_spec = 'reference'
remove_input_0_spec = TRUE
remove_median_0_spec = TRUE


#------------#
# FUNCTIONS  #
#------------#

# Recording filtering effects
add_sumstat = function(df, sumstats_df, filter_txt){
  
  n_variants = nrow(df)
  n_genes = length(unique(df$gene))
  
  sumstats_df[nrow(sumstats_df)+1, ] = c(filter_txt, n_variants, n_genes)
  
  return(sumstats_df)
  
}

# Data cleaning
keep_relevant_cols_from_df_categories = function(df_categories_raw){
  
  df_categories = df_categories_raw %>%
    rowwise() %>%
    select(gene_id, humvar, 
           variant_type, variant_type_adjusted,
           predicted_category,
           orf_ann_ref, orf_ann_alt,
           ref_seq, insert_seq)
  
  return(df_categories)
}


# Data processing
read_filter = function(df, n_min_reads){
  
  df = df %>%
    rowwise() %>%
    mutate(min_reads_input = min(input_12h_B1, input_12h_B2,
                                 input_12h_B3, input_12h_B4),
           min_reads_pulldown = min(pulldown_12h_B1, pulldown_12h_B2,
                                    pulldown_12h_B3, pulldown_12h_B4))
  
  if (n_min_reads >= 0){
    
    df = df %>%
      filter(min_reads_input >= n_min_reads | min_reads_pulldown >= n_min_reads) %>%
      ungroup()
    
  } 
  
  return(df)
  
}

remove_reporters_with_any_zeros = function(df){
  
  df = df %>%
    mutate(to_remove = ifelse(input_12h_B1 == 0 | input_12h_B2 == 0 |
                                input_12h_B3 == 0 | input_12h_B4 == 0, 1, 0)) %>%
    filter(to_remove == 0) %>%
    select(-to_remove)
  
  return(df)
  
}

remove_genes_with_0_median = function(df){
  
  # Calculate translation
  df = df %>%
    mutate(translation1 = pulldown_12h_B1/input_12h_B1,
           translation2 = pulldown_12h_B2/input_12h_B2,
           translation3 = pulldown_12h_B3/input_12h_B3,
           translation4 = pulldown_12h_B4/input_12h_B4)
  
  # Calculate median 
  df = df %>%
    group_by(gene) %>%
    
    # Compute median translation 
    mutate(median1 = median(translation1),
           median2 = median(translation2),
           median3 = median(translation3),
           median4 = median(translation4)) %>%
    
    ungroup() %>%
    
    filter(median1 != 0,
           median2 != 0,
           median3 != 0,
           median4 != 0)
  
  return(df)
  
}

add_epsilon = function(df, epsilon){
  
  df = df %>%
    mutate(input_12h_B1 = input_12h_B1 + epsilon,
           input_12h_B2 = input_12h_B2 + epsilon,
           input_12h_B3 = input_12h_B3 + epsilon,
           input_12h_B4 = input_12h_B4 + epsilon,
           
           pulldown_12h_B1 = pulldown_12h_B1 + epsilon, 
           pulldown_12h_B2 = pulldown_12h_B2 + epsilon,
           pulldown_12h_B3 = pulldown_12h_B3 + epsilon,
           pulldown_12h_B4 = pulldown_12h_B4 + epsilon)
  
  return(df)
  
}

normalize_replicates = function(df_raw, normalization, df_spikeins=NA){
  
  if (normalization == 'rpm'){
    
    df = df_raw %>%
      
      # Calculate RPM as ((reporter read count/Total Reads)*1,000,000)
      mutate(input_12h_B1_norm = (input_12h_B1/sum(input_12h_B1))*1000000,
             input_12h_B2_norm = (input_12h_B2/sum(input_12h_B2))*1000000,
             input_12h_B3_norm = (input_12h_B3/sum(input_12h_B3))*1000000,
             input_12h_B4_norm = (input_12h_B4/sum(input_12h_B4))*1000000,
             
             pulldown_12h_B1_norm = (pulldown_12h_B1/sum(pulldown_12h_B1))*1000000,
             pulldown_12h_B2_norm = (pulldown_12h_B2/sum(pulldown_12h_B2))*1000000,
             pulldown_12h_B3_norm = (pulldown_12h_B3/sum(pulldown_12h_B3))*1000000,
             pulldown_12h_B4_norm = (pulldown_12h_B4/sum(pulldown_12h_B4))*1000000)
    
  } else if (normalization == 'spikein'){
    
    df = df_raw %>%
      
      # Calculate spike-in normalized as (reporter read count/sum(spike-ins))
      mutate(input_12h_B1_norm = (input_12h_B1/sum(df_spikeins$input_12h_B1))*1000000,
             input_12h_B2_norm = (input_12h_B2/sum(df_spikeins$input_12h_B2))*1000000,
             input_12h_B3_norm = (input_12h_B3/sum(df_spikeins$input_12h_B3))*1000000,
             input_12h_B4_norm = (input_12h_B4/sum(df_spikeins$input_12h_B4))*1000000,
             
             pulldown_12h_B1_norm = (pulldown_12h_B1/sum(df_spikeins$pulldown_12h_B1))*1000000,
             pulldown_12h_B2_norm = (pulldown_12h_B2/sum(df_spikeins$pulldown_12h_B2))*1000000,
             pulldown_12h_B3_norm = (pulldown_12h_B3/sum(df_spikeins$pulldown_12h_B3))*1000000,
             pulldown_12h_B4_norm = (pulldown_12h_B4/sum(df_spikeins$pulldown_12h_B4))*1000000)
    
  }
  
  return(df)
  
}


# Calculate constructs 
calculate_translation = function(df, epsilon){
  
  df = df %>%
    
    ungroup() %>%
    
    mutate(translation1 = pulldown_12h_B1_norm/input_12h_B1_norm,
           translation2 = pulldown_12h_B2_norm/input_12h_B2_norm,
           translation3 = pulldown_12h_B3_norm/input_12h_B3_norm,
           translation4 = pulldown_12h_B4_norm/input_12h_B4_norm)
  
  return(df)
  
}

calculate_delta_translation = function(df, reference_txt){
  
  # Calculate median for all
  df = df %>%
    group_by(gene) %>%
    
    # Compute median translation 
    mutate(median1 = median(translation1),
           median2 = median(translation2),
           median3 = median(translation3),
           median4 = median(translation4)) %>%
    
    ungroup()
  
  if (reference_txt == 'median'){

    # Calculate delta translation 
    df = df %>%
      mutate(delta_t1 = (translation1/median1),
             delta_t2 = (translation2/median2),
             delta_t3 = (translation3/median3),
             delta_t4 = (translation4/median4)) %>%
      filter(humvar != 'ref')
        
  } else if (reference_txt == 'reference'){
    
    df_refs = df %>%
      filter(humvar == 'ref') %>%
      select(gene, translation1:translation4) %>%
      rename(ref1 = 'translation1',
             ref2 = 'translation2',
             ref3 = 'translation3',
             ref4 = 'translation4')
    
    df = df %>%
      filter(humvar != 'ref')
    
    genes_without_refs = df$gene[!(df$gene %in% df_refs$gene)]
    n_vars_affected = length(genes_without_refs)
    
    print(str_interp("${n_vars_affected} variants removed due to missing reference."))
    
    df = merge(df, df_refs, by='gene', all.x=TRUE, all.y=FALSE)
    
    df = df %>%
      
      filter(!(gene %in% genes_without_refs)) %>%
      
      # Calculate delta translation 
      mutate(delta_t1 = (translation1/ref1),
             delta_t2 = (translation2/ref2),
             delta_t3 = (translation3/ref3),
             delta_t4 = (translation4/ref4)) %>%
      
      select(-ref1, -ref2, -ref3, -ref4)
      
  }
  
  # Calculate log2_delta_t
  df = df %>%
    mutate(log2_delta_t1 = log2(delta_t1),
           log2_delta_t2 = log2(delta_t2),
           log2_delta_t3 = log2(delta_t3),
           log2_delta_t4 = log2(delta_t4))
  
  return(df)
  
}

calculate_fold_changes = function(df){
  
  df = df %>%
    
    # Haploid FC 
    mutate(haploid_fc1 = delta_t1-1,
           haploid_fc2 = delta_t2-1,
           haploid_fc3 = delta_t3-1,
           haploid_fc4 = delta_t4-1,
           
           # Diploid FC 
           diploid_fc1 = 1+0.5*haploid_fc1, 
           diploid_fc2 = 1+0.5*haploid_fc2, 
           diploid_fc3 = 1+0.5*haploid_fc3, 
           diploid_fc4 = 1+0.5*haploid_fc4, 
           
           # log2(Diploid FC)
           log2_diploid_fc1 = log2(diploid_fc1),
           log2_diploid_fc2 = log2(diploid_fc2),
           log2_diploid_fc3 = log2(diploid_fc3),
           log2_diploid_fc4 = log2(diploid_fc4)) %>%
    
    select(-starts_with("haploid"), -starts_with("diploid"))
  
  return(df)
  
}


# Calculate SEs 
solve_for_se_values = function(var_diffs) {
  
  var_values = solve(coef_matrix, var_diffs)
  se_values = sqrt(var_values)
  
  return(se_values)
  
}

calculate_SE_diploid = function(df){
  
  # Bin variants based on mean expected pulldown  
  df = df %>%
    
    # Calculate expected pulldown
    mutate(expected_pulldown_1 = input_12h_B1_norm*median1,
           expected_pulldown_2 = input_12h_B2_norm*median2,
           expected_pulldown_3 = input_12h_B3_norm*median3,
           expected_pulldown_4 = input_12h_B4_norm*median4) %>%
    
    rowwise() %>%
    
    mutate(mean_expected_pulldown = sum(expected_pulldown_1,
                                        expected_pulldown_2, 
                                        expected_pulldown_3,
                                        expected_pulldown_4)/4) %>%
    
    ungroup() %>%
    
    # Bin reporters based on expected pulldown 
    mutate(bins = ntile(mean_expected_pulldown, n_bins)) %>%
    select(-starts_with("expected_pulldown"))
  
  
  # Calculate var(diff_[1-4]) within bins
  df = df %>% 
    
    # Calculate var(beta-mean_beta)
    rowwise() %>%
    mutate(mean = (log2_diploid_fc1 + log2_diploid_fc2 +
                     log2_diploid_fc3 + log2_diploid_fc4)/4,
           
           diff_1 = log2_diploid_fc1 - mean,
           diff_2 = log2_diploid_fc2 - mean,
           diff_3 = log2_diploid_fc3 - mean,
           diff_4 = log2_diploid_fc4 - mean) %>%
    
    ungroup() %>% 
    group_by(bins) %>%
    
    mutate(var_diff_1 = var(diff_1),
           var_diff_2 = var(diff_2),
           var_diff_3 = var(diff_3),
           var_diff_4 = var(diff_4)) %>%
    
    ungroup()
  
  
  # Calculate SE values by solving for the 4 equations:
  #   var_diff_1 = (3/4)^2*se_1 + 1/(4^2)*(se_2 + se_3 + se_4)
  #   var_diff_2 = (3/4)^2*se_2 + 1/(4^2)*(se_1 + se_3 + se_4)
  #   var_diff_3 = (3/4)^2*se_3 + 1/(4^2)*(se_1 + se_2 + se_4)
  #   var_diff_4 = (3/4)^2*se_4 + 1/(4^2)*(se_1 + se_2 + se_3)
  
  df = df %>%
    
    rowwise() %>%
    mutate(se_log2_diploid = list(solve_for_se_values(c(var_diff_1, var_diff_2, 
                                                        var_diff_3, var_diff_4)))) %>%
    unnest_wider(se_log2_diploid, names_sep = "_") %>%
    ungroup()
  
  # Remove unnecessary columns 
  df = df %>%
    select(-bins, -mean, -starts_with("diff"), -starts_with("var_diff"))
  
  return(df)
  
}

calculate_SE_delta = function(df){
  
  # Bin variants based on mean expected pulldown  
  df = df %>%
    
    # Calculate expected pulldown
    mutate(expected_pulldown_1 = input_12h_B1_norm*median1,
           expected_pulldown_2 = input_12h_B2_norm*median2,
           expected_pulldown_3 = input_12h_B3_norm*median3,
           expected_pulldown_4 = input_12h_B4_norm*median4) %>%
    
    rowwise() %>%
    
    mutate(mean_expected_pulldown = sum(expected_pulldown_1,
                                        expected_pulldown_2, 
                                        expected_pulldown_3,
                                        expected_pulldown_4)/4) %>%
    
    ungroup() %>%
    
    # Bin reporters based on expected pulldown 
    mutate(bins = ntile(mean_expected_pulldown, n_bins)) %>%
    select(-starts_with("expected_pulldown"))
  
  
  # Calculate var(diff_[1-4]) within bins
  df = df %>% 
    
    # Calculate var(beta-mean_beta)
    rowwise() %>%
    mutate(mean = (log2_delta_t1 + log2_delta_t2 +
                     log2_delta_t3 + log2_delta_t4)/4,
           
           diff_1 = log2_delta_t1 - mean,
           diff_2 = log2_delta_t2 - mean,
           diff_3 = log2_delta_t3 - mean,
           diff_4 = log2_delta_t4 - mean) %>%
    
    ungroup() %>% 
    group_by(bins) %>%
    
    mutate(var_diff_1 = var(diff_1),
           var_diff_2 = var(diff_2),
           var_diff_3 = var(diff_3),
           var_diff_4 = var(diff_4)) %>%
    
    ungroup()
  
  
  # Calculate SE values by solving for the 4 equations:
  #   var_diff_1 = (3/4)^2*se_1 + 1/(4^2)*(se_2 + se_3 + se_4)
  #   var_diff_2 = (3/4)^2*se_2 + 1/(4^2)*(se_1 + se_3 + se_4)
  #   var_diff_3 = (3/4)^2*se_3 + 1/(4^2)*(se_1 + se_2 + se_4)
  #   var_diff_4 = (3/4)^2*se_4 + 1/(4^2)*(se_1 + se_2 + se_3)
  
  df = df %>%
    
    rowwise() %>%
    mutate(se_log2_delta_t = list(solve_for_se_values(c(var_diff_1, var_diff_2,
                                                        var_diff_3, var_diff_4)))) %>%
    unnest_wider(se_log2_delta_t, names_sep = "") %>%
    ungroup()
  
  # Remove unnecessary columns 
  df = df %>%
    select(-bins, -mean, -starts_with("diff"), -starts_with("var_diff"))
  
  return(df)
  
}


# SE QC plots
plot_mean_expected_pulldown_vs_SE = function(df, beta_variable){
  
  if (beta_variable == 'diploid'){
    
    add_theme = function(p){
      p = p +
        xlab("Mean expected pulldown") +
        ylab(expression(paste("SE of log"[2]*"(diploid FC)"))) +
        theme(axis.text = element_text(size=15),
              axis.title = element_text(size=15))
      
      return(p)
    }
    
    pA = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_diploid_1)) + 
      geom_point() + theme_light()
    pB = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_diploid_2)) + 
      geom_point() + theme_light()
    pC = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_diploid_3)) +
      geom_point() + theme_light()
    pD = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_diploid_4)) +
      geom_point() + theme_light()
    
  } else if (beta_variable == 'delta_t'){
    
    add_theme = function(p){
      p = p +
        xlab("Mean expected pulldown") +
        ylab(expression(paste("SE of log"[2]*" "*Delta*" T"))) +
        theme(axis.text = element_text(size=15),
              axis.title = element_text(size=15))
      
      return(p)
    }
    
    pA = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_delta_t1)) +
      geom_point() + theme_light()
    pB = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_delta_t2)) +
      geom_point() + theme_light()
    pC = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_delta_t3)) +
      geom_point() + theme_light()
    pD = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_delta_t4)) +
      geom_point() + theme_light()
    
  }
  
  # Save plot 
  png(str_interp('./output_plots/${beta_variable}_SE_vs_expected_pulldown.png'), 
      width=600, height=450)
  gridExtra::grid.arrange(add_theme(pA), add_theme(pB), 
                          add_theme(pC), add_theme(pD), nrow=2)  
  dev.off()
  
  # Print plot
  gridExtra::grid.arrange(add_theme(pA), add_theme(pB), 
                          add_theme(pC), add_theme(pD), nrow=2) 
  
}


# Calculate precision-weighted mean and SE 
calculate_prec_mean_and_se_diploid = function(df){
  
  df = df %>%
    
    rowwise() %>%
    
    # Calculate precision-weighted mean and SE 
    mutate(num = sum((log2_diploid_fc1)/(se_log2_diploid_1^2), 
                     (log2_diploid_fc2)/(se_log2_diploid_2^2), 
                     (log2_diploid_fc3)/(se_log2_diploid_3^2), 
                     (log2_diploid_fc4)/(se_log2_diploid_4^2)),
           
           denom = sum(1/(se_log2_diploid_1^2), 
                       1/(se_log2_diploid_2^2), 
                       1/(se_log2_diploid_3^2), 
                       1/(se_log2_diploid_4^2))) %>%
    
    mutate(pw_mean_log2_diploid = (num/denom),
           pw_se_log2_diploid = sqrt(1/denom)) %>%
    
    select(-num, -denom) %>%
    ungroup()
  
  return(df)
  
}

calculate_prec_mean_and_se_delta_t = function(df){
  
  df = df %>%
    
    rowwise() %>%
    
    mutate(num = sum(log2_delta_t1/(se_log2_delta_t1^2), 
                     log2_delta_t2/(se_log2_delta_t2^2), 
                     log2_delta_t3/(se_log2_delta_t3^2), 
                     log2_delta_t4/(se_log2_delta_t4^2)),
           
           denom = sum(1/(se_log2_delta_t1^2), 
                       1/(se_log2_delta_t2^2), 
                       1/(se_log2_delta_t3^2), 
                       1/(se_log2_delta_t4^2))) %>%
    
    mutate(pw_mean_log2_delta_t = num/denom,
           pw_se_log2_delta_t = sqrt(1/denom)) %>%
    
    select(-num, -denom) %>%
    ungroup()
  
  return(df)
  
}


# Ash
run_ash = function(df, n_categories, mean_col, se_col, 
                   posterior_mean_col, posterior_sd_col){
  
  df = data.frame(df)
  
  # Run ash on all the data 
  if (n_categories == 1){
    
    ash_results = 
      ash(df[,mean_col], df[,se_col],
          mixcompdist = 'halfnormal', pointmass=TRUE, pi_thresh=0)
    
    ash_results_df = ash_results$result
    ash_results_df = ash_results_df[,c('PosteriorMean', 'PosteriorSD')]
    names(ash_results_df) = c(posterior_mean_col, posterior_sd_col)
    
    df_post_ash = cbind(df, ash_results_df)
    
    # Run ash on predicted categories 
  } else if (n_categories == 4){
    
    # Categories 
    categories = unique(na.omit(df$predicted_category))
    
    # Run ash within each category 
    for (category in categories){
      
      # Subset variants in category 
      tmp_df = df %>%
        filter(predicted_category == category)
      
      if (category %in% c('Total loss')){
        
        # Run ash 
        ash_results = 
          ash(tmp_df[,mean_col], tmp_df[,se_col],
              mixcompdist = 'normal', pi_thresh=0, 
              pointmass=TRUE, mode='estimate')
        
      } else {
        
        # Run ash 
        ash_results = 
          ash(tmp_df[,mean_col], tmp_df[,se_col],
              mixcompdist = 'halfnormal', pointmass=TRUE, pi_thresh=0)
      }
      
      ash_results_df = ash_results$result
      ash_results_df = ash_results_df[,c('PosteriorMean', 'PosteriorSD')]
      names(ash_results_df) = c(posterior_mean_col, posterior_sd_col)
      
      df_merged_results = cbind(tmp_df, ash_results_df)
      
      if (category == categories[1]){
        df_post_ash = df_merged_results
      } else {
        df_post_ash = rbind(df_post_ash, df_merged_results)
      }
    }
    
  }
  
  return(df_post_ash)
  
}


# Make directory name from specifications 
make_subdirectory_name = function(){
  
  subdir = paste(read_filter_spec,
                 epsilon_spec,
                 normalization_spec,
                 reference_spec,
                 remove_input_0_spec,
                 remove_median_0_spec, sep='_')
  
  return(subdir)
  
}


# Ash plots
plot_pre_vs_post_ash_scatter = function(df, beta){
  
  if (beta == 'diploid'){
    
    pA = ggplot(df %>%
                  filter(is.na(predicted_category) == FALSE), 
                aes(x=pw_mean_log2_diploid, y=ash_posterior_mean_log2_diploid, 
                    color=pw_se_log2_diploid)) +
      geom_abline(intercept = 0, slope = 1, linetype='dashed', color='grey') +
      geom_point(alpha=0.5) + 
      facet_wrap(~predicted_category, nrow=1) +
      theme_light() +
      theme(axis.text = element_text(size=17),
            axis.title = element_text(size=17),
            legend.text = element_text(size=15),
            legend.title = element_text(size=17),
            strip.text = element_text(size=17, color='white'),
            strip.background = element_rect(fill = "black", color = "black")) +
      scale_color_viridis() +
      xlab(expression(paste("Pre-ash log"[2]*"(diploid FC)"))) +
      ylab(expression(paste("Post-ash log"[2]*"(diploid FC)"))) +
      labs(color = "SE (pre-ash)")
    
  } else if (beta == 'delta_t'){
    
    pA = ggplot(df %>%
                  filter(is.na(predicted_category) == FALSE), 
                aes(x=pw_mean_log2_delta_t, y=ash_posterior_mean_log2_delta_t, 
                    color=pw_se_log2_delta_t)) +
      geom_abline(intercept = 0, slope = 1, linetype='dashed', color='grey') +
      geom_point(alpha=0.5) + 
      facet_wrap(~predicted_category, nrow=1) +
      theme_light() +
      theme(axis.text = element_text(size=17),
            axis.title = element_text(size=17),
            legend.text = element_text(size=15),
            legend.title = element_text(size=17),
            strip.text = element_text(size=17, color='white'),
            strip.background = element_rect(fill = "black", color = "black")) +
      scale_color_viridis() +
      xlab(expression(paste("Pre-ash log"[2]*" "*Delta*" T"))) +
      ylab(expression(paste("Post-ash log"[2]*" "*Delta*" T"))) +
      labs(color = "SE (pre-ash)")
    
  } 
  
  # Save plot
  png(str_interp('${beta}_ash_by_categories.png'), width=800, height=450)
  plot(pA)
  dev.off()
  
  # Print plot
  print(pA)
  
}

plot_pre_and_post_ash_distributions = function(df, beta){
  
  if (beta == 'diploid'){
    
    pA = ggplot(df %>%
                  filter(is.na(predicted_category) == FALSE), 
                aes(x=pw_mean_log2_diploid)) +
      geom_density() + 
      facet_wrap(~predicted_category, nrow=1, scales='free') +
      theme_light() +
      theme(axis.text = element_text(size=15),
            axis.title = element_text(size=15),
            legend.text = element_text(size=15),
            legend.title = element_text(size=17),
            plot.title = element_text(size=19, hjust=0.5, face='bold'),
            strip.text = element_text(size=17, color='white'),
            strip.background = element_rect(fill = "black", color = "black")) +
      ggtitle("Pre-Ash") +
      xlab(expression(paste("Precision-weighted mean log"[2]*"(diploid FC)")))
    
    pB = ggplot(df %>%
                  filter(is.na(predicted_category) == FALSE), 
                aes(x=ash_posterior_mean_log2_diploid)) +
      geom_density() + 
      facet_wrap(~predicted_category, nrow=1, scales='free') +
      theme_light() +
      theme(axis.text = element_text(size=15),
            axis.title = element_text(size=15),
            legend.text = element_text(size=15),
            legend.title = element_text(size=17),
            plot.title = element_text(size=19, hjust=0.5, face='bold'),
            strip.text = element_text(size=17, color='white'),
            strip.background = element_rect(fill = "black", color = "black")) +
      ggtitle("Post-Ash") +
      xlab(expression(paste("Ash posterior log"[2]*"(diploid FC)")))
    
  } else if (beta == 'delta_t'){
    
    pA = ggplot(df %>%
                  filter(is.na(predicted_category) == FALSE), 
                aes(x=pw_mean_log2_delta_t)) +
      geom_density() + 
      facet_wrap(~predicted_category, nrow=1, scales='free') +
      theme_light() +
      theme(axis.text = element_text(size=15),
            axis.title = element_text(size=15),
            legend.text = element_text(size=15),
            legend.title = element_text(size=17),
            plot.title = element_text(size=19, hjust=0.5, face='bold'),
            strip.text = element_text(size=17, color='white'),
            strip.background = element_rect(fill = "black", color = "black")) +
      ggtitle("Pre-Ash") +
      xlab(expression(paste("Precision-weighted mean log"[2]*" "*Delta*" T")))
    
    pB = ggplot(df %>%
                  filter(is.na(predicted_category) == FALSE), 
                aes(x=ash_posterior_mean_log2_delta_t)) +
      geom_density() + 
      facet_wrap(~predicted_category, nrow=1, scales='free') +
      theme_light() +
      theme(axis.text = element_text(size=15),
            axis.title = element_text(size=15),
            legend.text = element_text(size=15),
            legend.title = element_text(size=17),
            plot.title = element_text(size=19, hjust=0.5, face='bold'),
            strip.text = element_text(size=17, color='white'),
            strip.background = element_rect(fill = "black", color = "black")) +
      ggtitle("Post-Ash") +
      xlab(expression(paste("Ash posterior mean log"[2]*" "*Delta*" T")))
    
  }
  
  # Save plot 
  png(str_interp('${beta}_ash_distributions.png'), 
      width=800, height=450)
  gridExtra::grid.arrange(pA, pB, nrow=2)
  dev.off()
  
  # Print plot
  gridExtra::grid.arrange(pA, pB, nrow=2)
  
}


# Leave-one-out validation
calculate_prec_mean_and_se_3replicates = function(df, beta){
  
  if (beta == 'diploid'){
    
    df = df %>%
      
      rowwise() %>%
      
      # Calculate precision-weighted mean and SE 
      mutate(num = sum((log2_diploid_fc1)/(se_log2_diploid_1^2), 
                       (log2_diploid_fc2)/(se_log2_diploid_2^2), 
                       (log2_diploid_fc3)/(se_log2_diploid_3^2)),
             
             denom = sum(1/(se_log2_diploid_1^2), 
                         1/(se_log2_diploid_2^2), 
                         1/(se_log2_diploid_3^2))) %>%
      
      mutate(pw_mean_3reps = (num/denom),
             pw_se_3reps = sqrt(1/denom)) %>%
      
      mutate(mean_3reps = sum(log2_diploid_fc1 +
                                log2_diploid_fc2 +
                                log2_diploid_fc3)/3) %>%
      
      select(-num, -denom) %>%
      ungroup()
    
  } else if (beta == 'delta_t'){
    
    df = df %>%
      
      rowwise() %>%
      
      # Calculate precision-weighted mean and SE 
      mutate(num = sum((log2_delta_t1)/(se_log2_delta_t1^2), 
                       (log2_delta_t2)/(se_log2_delta_t2^2), 
                       (log2_delta_t3)/(se_log2_delta_t3^2)),
             
             denom = sum(1/(se_log2_delta_t1^2), 
                         1/(se_log2_delta_t2^2), 
                         1/(se_log2_delta_t3^2))) %>%
      
      mutate(pw_mean_3reps = (num/denom),
             pw_se_3reps = sqrt(1/denom)) %>%
      
      mutate(mean_3reps = sum(log2_delta_t1 +
                                log2_delta_t2 +
                                log2_delta_t3)/3) %>%
      
      select(-num, -denom) %>%
      ungroup()
    
  }
  
  return(df)
  
}

leave_one_out_validation = function(df, beta, rep_4_colname){
  
  df_4replicates = df
  
  # Calculate precision-weighted means
  if (beta == 'diploid'){
    df_3replicates = calculate_prec_mean_and_se_3replicates(df, 'diploid')
  } else if (beta == 'delta_t'){
    df_3replicates = calculate_prec_mean_and_se_3replicates(df, 'delta_t')
  }
  
  # Calculate posterior ash means 
  df_3replicates = run_ash(df_3replicates, n_categories=4, 
                           'pw_mean_3reps', 'pw_se_3reps',
                           'ash_posterior_mean_3reps', 
                           'ash_posterior_sd_3reps')
  
  # Merge data 
  df_validate = merge(df_4replicates %>%
                        rename(replicate_4 = rep_4_colname) %>%
                        select(humvar, replicate_4),
                      
                      df_3replicates %>%
                        rename(simple_mean_3reps='mean_3reps') %>%
                        select(humvar, 
                               ash_posterior_mean_3reps, 
                               simple_mean_3reps), by='humvar')
  
  # Calculate means in bins
  df_validate = df_validate %>%
    mutate(bins = ntile(ash_posterior_mean_3reps, 100)) %>%
    
    group_by(bins) %>%
    mutate(mean_ash_posterior_3reps = mean(ash_posterior_mean_3reps),
           mean_replicate_4 = mean(replicate_4),
           mean_3reps = mean(simple_mean_3reps)) %>%
    ungroup()
  
  # Compute MSE
  mse1 = mean((df_validate$mean_ash_posterior_3reps - 
                 df_validate$mean_replicate_4)^2)
  
  mse2 = mean((df_validate$mean_3reps - 
                 df_validate$mean_replicate_4)^2)
  
  fraction = round(mse2/mse1, 1)
  
  print(str_interp("Ash reduces MSE by ${fraction}-fold."))
  
  # Plot 
  range = range(c(df_validate$mean_ash_posterior_3reps, 
                  df_validate$mean_3reps,
                  df_validate$mean_replicate_4))
  
  pA = ggplot(df_validate, aes(x=mean_ash_posterior_3reps, 
                               y=mean_replicate_4)) +
    geom_point() + geom_abline(slope=1) +
    ggtitle(str_interp("MSE = ${round(mse1, 5)}")) +
    theme_light() +
    theme(axis.text = element_text(size=17),
          axis.title = element_text(size=17),
          plot.title = element_text(size=19, face='bold', hjust=0.5)) +
    xlab("Mean ash posterior (3 replicates)") +
    ylab("Mean replicate 4") +
    xlim(range) + ylim(range)
  
  pB = ggplot(df_validate, aes(x=mean_3reps, 
                               y=mean_replicate_4)) +
    geom_point() + geom_abline(slope=1) +
    ggtitle(str_interp("MSE = ${round(mse2, 5)}")) +
    theme_light() +
    theme(axis.text = element_text(size=17),
          axis.title = element_text(size=17),
          plot.title = element_text(size=19, face='bold', hjust=0.5)) +
    xlab("Mean simple mean (3 replicates)") +
    ylab("Mean replicate 4") +
    xlim(range) + ylim(range)
  
  # Save plot 
  png(str_interp('${beta}_mse.png'), width=800, height=450)
  gridExtra::grid.arrange(pA, pB, nrow=1)
  dev.off()
  
  # Print plot
  gridExtra::grid.arrange(pA, pB, nrow=1)
  
}



#------------#
# VARIABLES  #
#------------#
n_bins = 100

coef_matrix = matrix(c(9/16, 1/16, 1/16, 1/16,
                       1/16, 9/16, 1/16, 1/16,
                       1/16, 1/16, 9/16, 1/16,
                       1/16, 1/16, 1/16, 9/16), nrow=4, byrow=TRUE)


#----------------------#
# LOAD & PROCESS DATA  #
#----------------------#

# NaP-TRAP data 
df_raw = read.csv('./data/raw/humvar_5utr_ntrap_v6_pm_reads_counts_041824.csv')
df_spikeins = df_raw %>% filter(grepl("spk", X))

df = df_raw %>%
  separate(X, into=c('gene', 'humvar'), "_")

# Predicted categories 
df_categories_raw = read.csv('./data/processed/annotated_variants_90k_revised.csv')
df_categories = keep_relevant_cols_from_df_categories(df_categories_raw)

# Merge NaP-TRAP and categories
df = merge(df, df_categories, 
           by.x=c('gene','humvar'), 
           by.y=c('gene_id', 'humvar'), all.x=TRUE, all.y=FALSE)


#-------------------#
# MAKE SUMSTATS DF  #
#-------------------#
sumstats_df = data.frame(matrix(nrow=0, ncol=3))
names(sumstats_df) = c('Filter', 'Variants', 'Genes')

sumstats_df = add_sumstat(df, sumstats_df, 'Raw data')


#---------------#
# PROCESS DATA  #
#---------------#

# 1. Remove spike-ins
df = df %>% filter(gene != 'ntrap')
sumstats_df = add_sumstat(df, sumstats_df, 'Remove spike-ins')

# 2. Read filter 
df = read_filter(df, read_filter_spec)
sumstats_df = add_sumstat(df, sumstats_df, 
                          str_interp("Read filter = ${read_filter_spec}"))

# 3. Remove input 0 (if specified)
if (remove_input_0_spec){
  df = remove_reporters_with_any_zeros(df)
  sumstats_df = add_sumstat(df, sumstats_df,
                            str_interp("0 input in any of the 4 replicates"))
}

# 4. Remove genes with median = 0
if (remove_median_0_spec){
  df = remove_genes_with_0_median(df)
  sumstats_df = add_sumstat(df, sumstats_df, "Remove genes with 0 median")
}

# 5. Add epsilon to all input and pulldown 
df = add_epsilon(df, epsilon_spec)

# 6. Normalize (RPM or spikein)
df = normalize_replicates(df, normalization_spec, df_spikeins)

# 7. Calculate translation values as pulldown/input
df = calculate_translation(df)

# 8. Calculate delta translation (based on median or reference)
df = calculate_delta_translation(df, reference_spec)

if (reference_spec == 'median'){
  sumstats_df = add_sumstat(df, sumstats_df, "Remove reference reporters")
} else {
  sumstats_df = add_sumstat(df, sumstats_df, "Remove reference reporters and variants without ref")
}


# 9. Remove variants w/o predicted category 
df = df %>% filter(is.na(predicted_category) == FALSE)
sumstats_df = add_sumstat(df, sumstats_df, 
                          "Remove variants w/o predicted category (NFASC tiling, missing references)")

# 10. Calculate fold changes
df = calculate_fold_changes(df)

# 11. Calculate SE based on bins of expected pulldown
df = calculate_SE_diploid(df)
df = calculate_SE_delta(df)

# QC plots
plot_mean_expected_pulldown_vs_SE(df, 'diploid')
plot_mean_expected_pulldown_vs_SE(df, 'delta_t')

# 12. Calculate precision-weighted mean and SE 
df = calculate_prec_mean_and_se_diploid(df)
df = calculate_prec_mean_and_se_delta_t(df)

# 13. Ash
df = run_ash(df, n_categories=4,
             'pw_mean_log2_diploid', 'pw_se_log2_diploid',
             'ash_posterior_mean_log2_diploid', 
             'ash_posterior_se_log2_diploid')

df = run_ash(df, n_categories=4,
             'pw_mean_log2_delta_t', 'pw_se_log2_delta_t',
             'ash_posterior_mean_log2_delta_t',
             'ash_posterior_se_log2_delta_t')


#-------------#
# SAVE DATA   #
#-------------#
subdir = make_subdirectory_name()
setwd(str_interp("~/Desktop/5prime_utr/output_plots/${subdir}"))

# Save data for gold standard and gold standard - reference 
if (subdir == '0_1_rpm_median_TRUE_TRUE'){
  
  write.csv(df, 
            str_interp('../../data/processed/humvar_5utr_ntrap_v6_ash_gold_standard.csv'), 
            row.names=FALSE)
  
} else if (subdir == '0_1_rpm_reference_TRUE_TRUE'){
  
  write.csv(df, 
            str_interp('../../data/processed/humvar_5utr_ntrap_v6_ash_gold_standard_reference.csv'), 
            row.names=FALSE)
}

write.csv(sumstats_df, 
          str_interp('filtering_sumstats_${subdir}.csv'), 
          row.names=FALSE)


#---------------------------------------#
# ASH PLOTS & LEAVE-ONE-OUT VALIDATION  #
#---------------------------------------#
plot_pre_vs_post_ash_scatter(df, 'diploid')
plot_pre_and_post_ash_distributions(df, 'diploid')
leave_one_out_validation(df, 'diploid', 'log2_diploid_fc4')

plot_pre_vs_post_ash_scatter(df, 'delta_t')
plot_pre_and_post_ash_distributions(df, 'delta_t')
leave_one_out_validation(df, 'delta_t', 'log2_delta_t4')


#-------------------------------------#
# EXTRA ANALYSES FOR JUSTIFICATIONS   #
#-------------------------------------#
analyzing_variants_low_read_count = function(read_filter_spec){
  
  df = df %>%
    rowwise() %>%
    filter(min_reads_input <= abs(read_filter_spec) | 
             min_reads_pulldown <= abs(read_filter_spec)) %>%
    ungroup()
  
  plot_pre_vs_post_ash_scatter(df, 'diploid')
  plot_pre_and_post_ash_distributions(df, 'diploid')
  
  plot_pre_vs_post_ash_scatter(df, 'delta_t')
  plot_pre_and_post_ash_distributions(df, 'delta_t')
  
  # Delta T 
  nrow(df)
  nrow(df %>% filter(pw_se_log2_delta_t > 0.3))
  nrow(df %>% filter(pw_se_log2_delta_t > 0.3))/nrow(df)
  
  nrow(df)
  nrow(df %>% filter(pw_se_log2_delta_t <= 0.3))
  nrow(df %>% filter(pw_se_log2_delta_t <= 0.3))/nrow(df)
  
  # Diploid FC
  nrow(df)
  nrow(df %>% filter(pw_se_log2_diploid > 0.15))
  nrow(df %>% filter(pw_se_log2_diploid > 0.15))/nrow(df)
  
  nrow(df)
  nrow(df %>% filter(pw_se_log2_diploid <= 0.15))
  nrow(df %>% filter(pw_se_log2_diploid <= 0.15))/nrow(df)
  
}
  
zero_input_variants = function(){
  
  # Do variants with zero input also have a zero output?
  pA = ggplot(df %>% filter(input_12h_B1 == 0), 
         aes(x=pulldown_12h_B1)) +
    geom_histogram() + theme_light() +
    ggtitle("Replicate 1") +
    theme(plot.title = element_text(size=14, hjust=0.5, face='bold'))
  
  pB = ggplot(df %>% filter(input_12h_B2 == 0), 
              aes(x=pulldown_12h_B2)) +
    geom_histogram() + theme_light() +
    ggtitle("Replicate 2") +
    theme(plot.title = element_text(size=14, hjust=0.5, face='bold'))
  
  pC = ggplot(df %>% filter(input_12h_B3 == 0), 
              aes(x=pulldown_12h_B3)) +
    geom_histogram() + theme_light() +
    ggtitle("Replicate 3") +
    theme(plot.title = element_text(size=14, hjust=0.5, face='bold'))
  
  pD = ggplot(df %>% filter(input_12h_B4 == 0), 
              aes(x=pulldown_12h_B4)) +
    geom_histogram() + theme_light() +
    ggtitle("Replicate 4") +
    theme(plot.title = element_text(size=14, hjust=0.5, face='bold'))
  
  gridExtra::grid.arrange(pA, pB, pC, pD, nrow=2)
  
  tmp_df = df %>% filter(input_12h_B1 == 0)
  perc_zero1 = round(nrow(tmp_df %>% 
                            filter(pulldown_12h_B1 == 0))/nrow(tmp_df)*100, 1)
  
  tmp_df = df %>% filter(input_12h_B2 == 0)
  perc_zero2 = round(nrow(tmp_df %>% 
                            filter(pulldown_12h_B2 == 0))/nrow(tmp_df)*100, 1)
  
  tmp_df = df %>% filter(input_12h_B3 == 0)
  perc_zero3 = round(nrow(tmp_df %>% 
                            filter(pulldown_12h_B3 == 0))/nrow(tmp_df)*100, 1)
  
  tmp_df = df %>% filter(input_12h_B4 == 0)
  perc_zero4 = round(nrow(tmp_df %>% 
                            filter(pulldown_12h_B4 == 0))/nrow(tmp_df)*100, 1)
  
  perc_zero1
  perc_zero2
  perc_zero3
  perc_zero4
  
}

comparing_medians_and_refs = function(){
  
  df_median = df 
  df_refs = df
  
  df_merged = merge(df_median %>% 
                      select(humvar, pw_mean_log2_diploid, 
                             ash_posterior_mean_log2_diploid,
                             median1:median4) %>%
                      rename(pre_ash_median = 'pw_mean_log2_diploid',
                             ash_median = 'ash_posterior_mean_log2_diploid'),
                    
                    df_refs %>%
                      select(humvar, pw_mean_log2_diploid, 
                             ash_posterior_mean_log2_diploid,
                             ref1:ref4) %>%
                      rename(pre_ash_ref = 'pw_mean_log2_diploid',
                             ash_ref = 'ash_posterior_mean_log2_diploid'))
  
  pA = ggplot(df_merged, aes(x=pre_ash_median, y=pre_ash_ref)) +
    geom_point() +
    geom_smooth(method='lm') +
    ggtitle("Pre-ash log2 diploid") +
    theme(plot.title = element_text(size=19, hjust=0.5, face='bold'))
  
  pB = ggplot(df_merged, aes(x=ash_median, y=ash_ref)) +
    geom_point() +
    geom_smooth(method='lm') +
    ggtitle("Post-ash log2 diploid") +
    theme(plot.title = element_text(size=19, hjust=0.5, face='bold'))
  
  gridExtra::grid.arrange(pA, pB, nrow=1)
  
  range = range(c(df_merged$median1, df_merged$ref1))
  pA = ggplot(df_merged, aes(x=median1, y=ref1)) +
    geom_point() +
    geom_smooth(method='lm') + 
    xlim(range) + ylim(range)
  
  range = range(c(df_merged$median2, df_merged$ref2))
  pB = ggplot(df_merged, aes(x=median2, y=ref2)) +
    geom_point() +
    geom_smooth(method='lm') + 
    xlim(range) + ylim(range)
  
  range = range(c(df_merged$median3, df_merged$ref3))
  pC = ggplot(df_merged, aes(x=median3, y=ref3)) +
    geom_point() +
    geom_smooth(method='lm') + 
    xlim(range) + ylim(range)
  
  range = range(c(df_merged$median4, df_merged$ref4))
  pD = ggplot(df_merged, aes(x=median4, y=ref4)) +
    geom_point() +
    geom_smooth(method='lm') + 
    xlim(range) + ylim(range)
  
  gridExtra::grid.arrange(pA, pB, pC, pD, nrow=1)
  
  # How frequently are all references 0's and the median is different?
  tmp = df_merged %>%
    rowwise() %>%
    mutate(mean_ref = mean(ref1, ref2, ref3, ref4),
           mean_median = mean(median1, median2, median3, median4))
  
  range = range(c(tmp$mean_median, tmp$mean_ref))
  ggplot(tmp, aes(x=mean_ref, y=mean_median)) + 
    geom_point() + xlim(range) + ylim(range) +
    geom_smooth(method='lm') #+ xlim(c(0, 0.5))
  
}

looking_into_zero_medians = function(){
  
  zero_median_genes = df %>%
    filter(median1 == 0 | median2 == 0 | median3 == 0 | median4 == 0)
  
  n_genes = length(unique(zero_median_df$gene))
  
  n_vars_per_gene = zero_median_df %>%
    group_by(gene) %>%
    mutate(n_variants = n()) 
  
  range(n_vars_per_gene$n_variants)
  
  n_vars_per_gene = n_vars_per_gene %>%
    
    filter(humvar == 'ref') %>%
    rowwise() %>%
    mutate(n_zeros_in_ref = 
             (translation1 == 0) + 
             (translation2 == 0) +
             (translation3 == 0) + 
             (translation4 == 0),
           any_zeros_in_ref = ifelse(n_zeros_in_ref > 0, 'yes', 'no')) %>%
    
    select(gene, n_variants, n_zeros_in_ref, any_zeros_in_ref)
  
  table(n_vars_per_gene$any_zeros_in_ref)
  
  BCR_df = zero_median_df %>%
    filter(gene == 'BCR')
  
  pA = ggplot(BCR_df, aes(x=translation1)) + geom_histogram()
  pB = ggplot(BCR_df, aes(x=translation2)) + geom_histogram()
  pC = ggplot(BCR_df, aes(x=translation3)) + geom_histogram()
  pD = ggplot(BCR_df, aes(x=translation4)) + geom_histogram()
  
  gridExtra::grid.arrange(pA, pB, pC, pD, nrow=2)
  
  write.csv(zero_median_genes, '~/Downloads/genes_with_median_0.csv',
            row.names = FALSE)
  
}

comparing_median_and_ref_post_ash = function(){
  setwd('~/Desktop/5prime_utr/')
  
  #------------#
  # LOAD DATA  #
  #------------#
  df_median = read.csv('./data/processed/humvar_5utr_ntrap_v6_ash_gold_standard.csv')
  df_ref = read.csv('./data/processed/humvar_5utr_ntrap_v6_ash_gold_standard_reference.csv')
  
  
  #-------------#
  # MERGE DATA  #
  #-------------#
  df = merge(df_median %>%
               select(humvar, 
                      ash_posterior_mean_log2_delta_t,
                      ash_posterior_se_log2_delta_t,
                      
                      ash_posterior_mean_log2_diploid,
                      ash_posterior_se_log2_diploid) %>%
               
               rename(ash_log2_delta_t_median = 'ash_posterior_mean_log2_delta_t',
                      ash_se_log2_delta_t_median = 'ash_posterior_se_log2_delta_t',
                      
                      ash_log2_diploid_median = 'ash_posterior_mean_log2_diploid',
                      ash_se_log2_diploid_median = 'ash_posterior_se_log2_diploid'),
             
             df_ref %>%
               select(humvar, 
                      ash_posterior_mean_log2_delta_t,
                      ash_posterior_se_log2_delta_t,
                      
                      ash_posterior_mean_log2_diploid,
                      ash_posterior_se_log2_diploid) %>%
               
               rename(ash_log2_delta_t_ref = 'ash_posterior_mean_log2_delta_t',
                      ash_se_log2_delta_t_ref = 'ash_posterior_se_log2_delta_t',
                      
                      ash_log2_diploid_ref = 'ash_posterior_mean_log2_diploid',
                      ash_se_log2_diploid_ref = 'ash_posterior_se_log2_diploid'), by='humvar')
  
  #--------#
  # PLOTS  #
  #--------#
  ggplot(df, aes(x=ash_log2_delta_t_median, y=ash_log2_delta_t_ref)) +
    geom_point() + geom_smooth(method='lm')
  
  ggplot(df, aes(x=ash_log2_diploid_median, y=ash_log2_diploid_ref)) +
    geom_point() + geom_smooth(method='lm')
  
  ggplot(df, aes(x=ash_se_log2_delta_t_median, y=ash_se_log2_delta_t_ref)) +
    geom_point() + geom_smooth(method='lm')
  
  ggplot(df, aes(x=ash_se_log2_diploid_median, y=ash_se_log2_diploid_ref)) +
    geom_point() + geom_smooth(method='lm')
}
