rm(list=ls())

library(ggplot2)
library(dplyr)

setwd('~/Desktop/5prime_utr/')

df = read.csv('./data/humvar_5utr_ntrap_ash_v7.csv')


#-------------------#
# OUTLIER ANALYSIS  #
#-------------------#
ggplot(df %>%
         filter(is.na(predicted_category) == FALSE), 
       aes(x=pw_mean_log2_diploid, y=ash_posterior_mean_log2_diploid, 
           color=pw_se_log2_diploid)) +
  geom_abline(intercept = 0, slope = 1, linetype='dashed', color='grey') +
  geom_point(alpha=0.5) + 
  facet_wrap(~predicted_category, nrow=1) +
  theme_light() +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=14),
        strip.text = element_text(size=14, color='white'),
        strip.background = element_rect(fill = "black", color = "black")) +
  scale_color_viridis() +  
  xlab("Raw beta") + ylab("Posterior beta") +
  labs(color = "SE (raw beta)")


# NEGATIVE outliers --------------------------------------#
df_neg_outliers = df %>%
  filter(predicted_category == 'Negative') %>%
  filter(ash_posterior_mean_log2_diploid > 0.58) %>%
  select(humvar, gene, 
         orf_ann_ref, orf_ann_alt,
         variant_type, variant_type_adjusted,
         ash_posterior_mean_delta,
         ref_seq, insert_seq)

tmp = data.frame(table(df_neg_outliers$variant_type_adjusted))
View(tmp)


# POSITIVE outliers --------------------------------------#
df_pos_outliers = df %>%
  filter(predicted_category == 'Positive') %>%
  filter(ash_posterior_mean_log2_diploid < -0.5) %>%
  select(humvar, gene, 
         orf_ann_ref, orf_ann_alt,
         variant_type, variant_type_adjusted,
         ash_posterior_mean_delta,
         ref_seq, insert_seq)

tmp = data.frame(table(df_pos_outliers$variant_type_adjusted))
View(tmp)


# NULL outliers (negative) -------------------------------#
df_null_negative_outliers = df %>%
  filter(predicted_category == 'Null') %>%
  filter(ash_posterior_mean_log2_diploid < -0.5) %>%
  select(humvar, gene, 
         orf_ann_ref, orf_ann_alt,
         variant_type, variant_type_adjusted,
         ash_posterior_mean_delta)

tmp = data.frame(table(df_null_negative_outliers$variant_type_adjusted))
View(tmp)

# NULL outliers (positive) -------------------------------#
df_null_positive_outliers = df %>%
  filter(predicted_category == 'Null') %>%
  filter(ash_posterior_mean_log2_diploid > 0.58) %>%
  select(humvar, gene, 
         orf_ann_ref, orf_ann_alt,
         variant_type, variant_type_adjusted,
         ash_posterior_mean_delta)

tmp = data.frame(table(df_null_positive_outliers$variant_type_adjusted))
View(tmp)

# TOTAL LOSS outliers (positive) --------------------------#
df_total_loss_outliers = df %>%
  filter(predicted_category == 'Total loss') %>%
  filter(ash_posterior_mean_log2_diploid > 0) %>%
  select(humvar, gene, 
         orf_ann_ref, orf_ann_alt,
         variant_type, variant_type_adjusted,
         ash_posterior_mean_delta)

tmp = data.frame(table(df_total_loss_outliers$variant_type_adjusted))
View(tmp)
