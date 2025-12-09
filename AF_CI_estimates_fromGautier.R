################################################################################
# Simple Pool-seq Confidence Intervals
# Accounts for: (1) pooling error, (2) finite sample size, (3) sequencing depth
################################################################################

library(dplyr)

# Load your data
# Expected columns: library_id, snp_id, replicate (if applicable), 
#                   alt_count, ref_count, pool_size
data <- read.csv("poolseqdataorsomethinglikethat.csv")

################################################################################
# STEP 1: Estimate epsilon from technical replicates
################################################################################

################################################################################
# Estimate epsilon from genome-wide SNPs
################################################################################

# Use a subset of genome-wide SNPs with good coverage to estimate epsilon
snp_stats_genomewide <- rep_data %>%
  filter(total_depth >= 5, p >= 0.05, p <= 0.95) %>%  # Filter for depth and MAF
  group_by(library_id, snp_id) %>%
  summarise(
    p_mean = mean(p),
    p_var = var(p),
    cov_mean = mean(total_depth),
    .groups = "drop"
  ) %>%
  mutate(
    var_sequencing = p_mean * (1 - p_mean) / cov_mean,
    var_pooling = pmax(0, p_var - var_sequencing),
    epsilon_sq = var_pooling / (p_mean * (1 - p_mean))
  ) %>%
  filter(epsilon_sq >= 0, epsilon_sq < 1)

# Global epsilon from all SNPs
epsilon <- sqrt(median(snp_stats_genomewide$epsilon_sq))

cat(sprintf("Genome-wide epsilon: %.3f (from %d SNPs)\n", 
            epsilon, nrow(snp_stats_genomewide)))


################################################################################
# Apply to resistance SNPs specifically
################################################################################

# Filter to just resistance SNPs
resistance_snps <- c("ALS_574", "ALS_653", "PPO_210", ...)  # resistance SNP IDs

resistance_data <- data_avg %>%
  filter(snp_id %in% resistance_snps)

# Calculate CIs using genome-wide epsilon
resistance_data <- resistance_data %>%
  mutate(
    total_depth = alt_count + ref_count,
    p = alt_count / total_depth,
    n_eff = pool_size / (1 + epsilon^2),
    n_eff_chrom = 2 * n_eff,
    variance = case_when(
      total_depth < n_eff_chrom ~ 
        (p * (1 - p) / total_depth) * ((n_eff_chrom - total_depth) / (n_eff_chrom - 1)),
      TRUE ~ p * (1 - p) / n_eff_chrom
    ),
    SE = sqrt(variance),
    CI_lower = pmax(0, p - 1.96 * SE),
    CI_upper = pmin(1, p + 1.96 * SE)
  )

# Summary for resistance SNPs
cat("\nResistance SNP summary:\n")
cat(sprintf("  N resistance SNPs: %d\n", nrow(resistance_data)))
cat(sprintf("  Mean frequency: %.3f\n", mean(resistance_data$p)))
cat(sprintf("  Mean SE: %.4f\n", mean(resistance_data$SE)))
