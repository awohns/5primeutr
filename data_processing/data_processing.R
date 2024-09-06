rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(stringr)

library(plotrix)
library(devtools)
library(ashr)
library(nnls)
library(viridis)
library(DescTools)

setwd('/Users/awohns/Documents/5primeutr/error_modelling')


#-----------------#
# SPECIFICATIONS  #
#-----------------#
read_filter_spec = 0
epsilon_spec = 1
normalization_spec = 'rpm'
reference_spec = 'median'
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
    rename(predicted_category = 'predicted_category') %>%
    select(gene_id, humvar, 
           variant_type,
           ref_seq, insert_seq, predicted_category)
  
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

calculate_expected_pulldown = function(df){
  # Bin variants based on mean expected pulldown  
  df = df %>%
    
    # Calculate expected pulldown
    mutate(expected_pulldown_1 = input_12h_B1_norm*median1,
           expected_pulldown_2 = input_12h_B2_norm*median2,
           expected_pulldown_3 = input_12h_B3_norm*median3,
           expected_pulldown_4 = input_12h_B4_norm*median4) %>%
    
    ungroup()
  
  return(df)
}

# Calculate SEs 
get_design_matrix <- function(x, num_bins) {
  n <- nrow(x)
  p <- ncol(x)
  
  # Compute the bin edges based on quantiles
  x_bin <- quantile(array(unlist(x)),
                    probs = seq(0, 1,
                                length.out = num_bins + 1
                    )
  )[-(num_bins + 1)]
  
  # Digitize x based on the computed bin edges
  bins <- matrix(findInterval(array(unlist(x)), x_bin, left.open = TRUE), n, p)
  
  # Initialize the one-hot encoding vector, v
  v <- array(0, dim = c(n, p, num_bins))
  for (i in 1:n) {
    for (j in 1:p) {
      v[i, j, 1:bins[i, j]] <- 1
    }
  }
  
  # Compute the coefficients for the diagonal and off-diagonal elements
  diag_coef <- ((p - 1) / p)^2
  offdiag_coef <- (1 / p)^2
  
  # Calculate u: made up of one hot encoding and coefficent matrix
  u <- array(0, dim = c(n, p, num_bins))
  for (i in 1:n) {
    for (k in 1:p) {
      u[i, k, ] <- (diag_coef - offdiag_coef) *
        v[i, k, ] + offdiag_coef * colSums(v[i, , ])
    }
  }
  
  return(list(u = u, v = v))
}

# Function to estimate standard error
estimate_standard_error <- function(y, x, num_bins) {
  x <- 1 / x
  n <- nrow(x)
  p <- ncol(x)
  
  # Get design matrix
  design_matrices <- get_design_matrix(x, num_bins)
  A <- design_matrices$u
  V <- design_matrices$v
  
  residuals <- y - matrix(
    rep(apply(y, 1, mean), ncol(y)),
    nrow = nrow(y), ncol = ncol(y)
  )
  
  A_transposed <- aperm(A, c(3, 2, 1))
  dims <- dim(A)
  nrow <- prod(dims[1:2]) # Number of rows to reshape to
  ncol <- dims[3]
  A_flat_row_major <- as.vector(A_transposed)
  A <- matrix(A_flat_row_major, nrow = nrow, ncol = ncol, byrow = TRUE)
  V_transposed <- aperm(V, c(3, 2, 1))
  V_flat_row_major <- as.vector(V_transposed)
  V <- matrix(V_flat_row_major, nrow = nrow, ncol = ncol, byrow = TRUE)
  residuals <- as.vector(t(residuals))
  
  # Solve for theta using NNLS
  nnls_result <- nnls(A, residuals^2)
  theta <- coef(nnls_result)
  
  return(matrix(V %*% theta, nrow = n, ncol = p, byrow = TRUE))
}

# Adds estimated standard error to dataframe
calculate_se <- function(df, num_bins) {
  x <- as.matrix(
    df[c(
      "expected_pulldown_1", "expected_pulldown_2",
      "expected_pulldown_3", "expected_pulldown_4"
    )]
  )
  y <- as.matrix(
    df[c(
      "log2_diploid_fc1", "log2_diploid_fc2",
      "log2_diploid_fc3", "log2_diploid_fc4"
    )]
  )
  estimated_error <- sqrt(estimate_standard_error(y, x, num_bins))
  stopifnot(all(estimated_error >= 0))
  df_combined <- cbind(df, estimated_error)
  colnames(df_combined)[
    (ncol(df) + 1):ncol(df_combined)
  ] <- c(
    "se_log2_diploid_1", "se_log2_diploid_2",
    "se_log2_diploid_3", "se_log2_diploid_4"
  )
  df <- df_combined
  
  return(df)
}

# SE QC plots
plot_mean_expected_pulldown_vs_SE = function(df, beta_variable){

  if (beta_variable == 'diploid'){
  
    pA = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_diploid_1)) + 
      geom_point() + theme_light()
    pB = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_diploid_2)) + 
      geom_point() + theme_light()
    pC = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_diploid_3)) +
      geom_point() + theme_light()
    pD = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_diploid_4)) +
      geom_point() + theme_light()
    
    print(gridExtra::grid.arrange(pA, pB, pC, pD, nrow=2))
    
  } else if (beta_variable == 'delta_t'){
    
    pA = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_delta_t1)) +
      geom_point() + theme_light()
    pB = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_delta_t2)) +
      geom_point() + theme_light()
    pC = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_delta_t3)) +
      geom_point() + theme_light()
    pD = ggplot(df, aes(x=mean_expected_pulldown, y=se_log2_delta_t4)) +
      geom_point() + theme_light()
    
    print(gridExtra::grid.arrange(pA, pB, pC, pD, nrow=2))
    
  }
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
              mixcompdist = 'normal', pi_thresh=0, pointmass=TRUE,
              mode='estimate')
        
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
      theme(axis.text = element_text(size=14),
            axis.title = element_text(size=14),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            strip.text = element_text(size=14, color='white'),
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
      theme(axis.text = element_text(size=14),
            axis.title = element_text(size=14),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            strip.text = element_text(size=14, color='white'),
            strip.background = element_rect(fill = "black", color = "black")) +
      scale_color_viridis() +
      xlab(expression(paste("Pre-ash log"[2]*" "*Delta*" T"))) +
      ylab(expression(paste("Post-ash log"[2]*" "*Delta*" T"))) +
      labs(color = "SE (pre-ash)")
    
    range = range(df$pw_mean_log2_delta_t)
    pB = pA + ylim(range)
    
    print(pB)
  } 
  
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
      theme(axis.text = element_text(size=12),
            axis.title = element_text(size=12),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            plot.title = element_text(size=16, hjust=0.5, face='bold'),
            strip.text = element_text(size=14, color='white'),
            strip.background = element_rect(fill = "black", color = "black")) +
      ggtitle("Pre-Ash") +
      xlab(expression(paste("Precision-weighted mean log"[2]*"(diploid FC)")))
    
    pB = ggplot(df %>%
                  filter(is.na(predicted_category) == FALSE), 
                aes(x=ash_posterior_mean_log2_diploid)) +
      geom_density() + 
      facet_wrap(~predicted_category, nrow=1, scales='free') +
      theme_light() +
      theme(axis.text = element_text(size=12),
            axis.title = element_text(size=12),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            plot.title = element_text(size=16, hjust=0.5, face='bold'),
            strip.text = element_text(size=14, color='white'),
            strip.background = element_rect(fill = "black", color = "black")) +
      ggtitle("Post-Ash") +
      xlab(expression(paste("Ash posterior log"[2]*"(diploid FC)")))
    
    gridExtra::grid.arrange(pA, pB, nrow=2)
    
  } else if (beta == 'delta_t'){
    
    pA = ggplot(df %>%
                  filter(is.na(predicted_category) == FALSE), 
                aes(x=pw_mean_log2_delta_t)) +
      geom_density() + 
      facet_wrap(~predicted_category, nrow=1, scales='free') +
      theme_light() +
      theme(axis.text = element_text(size=12),
            axis.title = element_text(size=12),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            plot.title = element_text(size=16, hjust=0.5, face='bold'),
            strip.text = element_text(size=14, color='white'),
            strip.background = element_rect(fill = "black", color = "black")) +
      ggtitle("Pre-Ash") +
      xlab(expression(paste("Precision-weighted mean log"[2]*" "*Delta*" T")))
    
    pB = ggplot(df %>%
                  filter(is.na(predicted_category) == FALSE), 
                aes(x=ash_posterior_mean_log2_delta_t)) +
      geom_density() + 
      facet_wrap(~predicted_category, nrow=1, scales='free') +
      theme_light() +
      theme(axis.text = element_text(size=12),
            axis.title = element_text(size=12),
            legend.text = element_text(size=12),
            legend.title = element_text(size=14),
            plot.title = element_text(size=16, hjust=0.5, face='bold'),
            strip.text = element_text(size=14, color='white'),
            strip.background = element_rect(fill = "black", color = "black")) +
      ggtitle("Post-Ash") +
      xlab(expression(paste("Ash posterior mean log"[2]*" "*Delta*" T")))
    
    gridExtra::grid.arrange(pA, pB, nrow=2)
    
  }
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
  
  df_3replicates = calculate_prec_mean_and_se_3replicates(df, 'diploid')
  df_3replicates = run_ash(df_3replicates, n_categories=4, 
                           'pw_mean_3reps', 'pw_se_3reps',
                           'ash_posterior_mean_3reps', 'ash_posterior_sd_3reps')
  
  if (beta == 'delta_t'){
   
    df_3replicates = transform_diploid_to_haploid(df_3replicates,
                                                  'ash_posterior_mean_3reps',
                                                  'ash_posterior_sd_3reps', all_reps=FALSE)
  }
  
  df_validate = merge(df_4replicates %>%
                        rename(replicate_4 = rep_4_colname) %>%
                        select(humvar, replicate_4),
                      
                      df_3replicates %>%
                        rename(simple_mean_3reps='mean_3reps') %>%
                        select(humvar, 
                               ash_posterior_mean_3reps, 
                               simple_mean_3reps), by='humvar')
  
  df_validate = df_validate %>%
    mutate(bins = ntile(ash_posterior_mean_3reps, 100)) %>%
    
    group_by(bins) %>%
    mutate(mean_ash_posterior_3reps = 2^mean(ash_posterior_mean_3reps),
           mean_replicate_4 = 2^mean(replicate_4),
           mean_3reps = 2^mean(simple_mean_3reps)) %>%
    ungroup()
  
  mse1 = mean((df_validate$mean_ash_posterior_3reps - 
                 df_validate$mean_replicate_4)^2)
  
  mse2 = mean((df_validate$mean_3reps - 
                 df_validate$mean_replicate_4)^2)
  
  range = range(c(df_validate$mean_ash_posterior_3reps, 
                  df_validate$mean_3reps,
                  df_validate$mean_replicate_4))

  pA = ggplot(df_validate, aes(x=mean_ash_posterior_3reps, y=mean_replicate_4)) +
    geom_point() + geom_abline(slope=1) +
    ggtitle(str_interp("MSE = ${round(mse1, 4)}")) +
    theme_light() +
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14),
          plot.title = element_text(size=16, face='bold', hjust=0.5)) +
    xlab("Mean ash posterior (3 replicates)") +
    ylab("Mean replicate 4") +
    xlim(range) + ylim(range)
  
  pB = ggplot(df_validate, aes(x=mean_3reps, y=mean_replicate_4)) +
    geom_point() + geom_abline(slope=1) +
    ggtitle(str_interp("MSE = ${round(mse2, 4)}")) +
    theme_light() +
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14),
          plot.title = element_text(size=16, face='bold', hjust=0.5)) +
    xlab("Mean simple mean (3 replicates)") +
    ylab("Mean replicate 4") +
    xlim(range) + ylim(range)
  
  gridExtra::grid.arrange(pA, pB, nrow=1)
  
}


#------------#
# VARIABLES  #
#------------#
n_bins = 100

#----------------------#
# LOAD & PROCESS DATA  #
#----------------------#

# NaP-TRAP data 
df_raw = read.csv('~/Documents/5primeutr/data/humvar_5utr_ntrap_v6_pm_reads_counts_041824.csv')
df_spikeins = df_raw %>% filter(grepl("spk", X))

df = df_raw %>%
  separate(X, into=c('gene', 'humvar'), "_")

# Predicted categories 
df_categories_raw = read.csv('~/Downloads/annotated_variants_90k_revised.csv')
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

# 3. Remove input 0 (if specified )
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
sumstats_df = add_sumstat(df, sumstats_df, "Remove reference reporters")


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
df = calculate_expected_pulldown(df)
df = calculate_se(df, n_bins)

# QC plots
# plot_mean_expected_pulldown_vs_SE(df, 'diploid')
# plot_mean_expected_pulldown_vs_SE(df, 'delta_t')

# 12. Calculate precision-weighted mean and SE 
df = calculate_prec_mean_and_se_diploid(df)

# 13. Ash
df = run_ash(df, n_categories=4, 
             'pw_mean_log2_diploid', 'pw_se_log2_diploid',
             'ash_posterior_mean_log2_diploid', 'ash_posterior_se_log2_diploid')

# df = run_ash(df, n_categories=4,
#             'pw_mean_log2_delta_t', 'pw_se_log2_delta_t',
#             'ash_posterior_mean_log2_delta_t', 'ash_posterior_se_log2_delta_t')

sumstats_df = add_sumstat(df, sumstats_df, "Remove variants with missing predicted category")

# 14. Transform diploid back to log2 delta
transform_diploid_to_haploid = function(df, 
                                        ash_diploid_mean_colname, 
                                        ash_diploid_se_colname, all_reps){
  
  df = df %>%
    
    rename(diploid=ash_diploid_mean_colname,
           diploid_se=ash_diploid_se_colname) %>%

    rowwise() %>% 
    
    # Naive transformation of diploid to haploid  
    mutate(naive_haploid = log2((2^diploid-1)*2 + 1),
           
           # Second derivative of naive transformation 
           d2_dx2_haploid = - (2^(diploid+1)*log(2)) / ((2^(diploid+1)-1)^2),
           
           # Transformation
           haploid = naive_haploid + 0.5*d2_dx2_haploid*(diploid^2 + diploid_se))
  
  if (all_reps){
    
    # Rename back to ash_posterior_mean_log2_delta_t so it's compatible with 
    # downstream functions
    df = df %>%
      rename(ash_posterior_mean_log2_delta_t='haploid',
             
             # Rename diploid colnames back 
             ash_posterior_mean_log2_diploid='diploid',
             ash_posterior_se_log2_diploid='diploid_se')
    
  } else {
    
    df = df %>%
      rename(ash_posterior_mean_3reps='haploid')
  }
           
  return(df)
  
  
}

df = transform_diploid_to_haploid(df, 
                                  'ash_posterior_mean_log2_diploid',
                                  'ash_posterior_se_log2_diploid', all_reps=TRUE)

# 15. Save data
write.csv(df, '../data/humvar_5utr_ntrap_v6_ash_epsilon_1_08_07_24.csv', row.names=FALSE)


#-------------------------------------#
# ASH PLOTS — COMPARISONS & OUTLIERS  #
#-------------------------------------#
plot_pre_vs_post_ash_scatter(df, 'diploid')
plot_pre_vs_post_ash_scatter(df, 'delta_t')

plot_pre_and_post_ash_distributions(df, 'diploid')
plot_pre_and_post_ash_distributions(df, 'delta_t')


#---------------------------#
# LEAVE ONE OUT VALIDATION  #
#---------------------------#
leave_one_out_validation(df, 'diploid', 'log2_diploid_fc4')
leave_one_out_validation(df, 'delta_t', 'log2_delta_t4')
