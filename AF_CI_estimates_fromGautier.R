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


#this is actually a rough approximation of gautier because it doesn't account for the fact that increasing depth above pool size marginally increases the precision of AF estimation.
#can directly use equation 4 to estimate effective pool size 

################################################################################
# STEP 2REFINED: Calculate CI using Gautier Equation 4
################################################################################

calc_poolseq_ci_gautier <- function(alt_count, ref_count, pool_size, 
                                    epsilon, s_lambda = 0) {
  #' @param alt_count: alternate allele read count
  #' @param ref_count: reference allele read count  
  #' @param pool_size: number of diploid individuals
  #' @param epsilon: experimental error (pooling CV)
  #' @param s_lambda: overdispersion parameter (default 0 = Poisson)
  #' @return named vector with p, SE, CI_lower, CI_upper
  
  depth <- alt_count + ref_count
  if (depth == 0) return(c(p = NA, SE = NA, CI_lower = NA, CI_upper = NA))
  
  # Allele frequency estimate
  p <- alt_count / depth
  
  # Effective pool size (accounting for unequal contributions)
  n_e <- pool_size / (1 + epsilon^2)
  n_e_chrom <- 2 * n_e  # haploid chromosome sample size
  
  # Gautier Equation 4: V(p̂) = [p(1-p)/np] × [1 + (λp + sλp)/λp² × (np - 1)]
  ####################        \_________/     \____________________________/
  ######################         Part A                   Part B
  #######################    Base variance          Correction factor
  # Base variance from pool size
  var_base <- p * (1 - p) / n_e_chrom
  
  # Correction factor for sequencing
  if (depth > 0) {
    correction_factor <- 1 + ((depth + s_lambda) / depth^2) * (n_e_chrom - 1)
  } else {
    correction_factor <- 1
  }
  
  # Total variance
  variance <- var_base * correction_factor
  
  # Standard error and 95% CI
  se <- sqrt(variance)
  ci_lower <- max(0, p - 1.96 * se)
  ci_upper <- min(1, p + 1.96 * se)
  
  return(c(p = p, SE = se, CI_lower = ci_lower, CI_upper = ci_upper))
}

################################################################################
# STEP 3: Apply to resistance SNPs, as above
################################################################################

add_poolseq_cis_resistance <- function(data, epsilon, s_lambda = 0) {
  #' Add CIs to resistance SNPs using Gautier equation
  #' 
  #' @param data: data.frame with columns:
  #'   - gene (resistance gene name)
  #'   - mutation (amino acid change, e.g., "Ser264Gly")
  #'   - library (sample ID)
  #'   - alt_count (resistance allele reads)
  #'   - ref_count (susceptible allele reads)
  #'   - pool_size (number of diploid individuals)
  #' @param epsilon: experimental error from genome-wide estimate
  #' @param s_lambda: overdispersion parameter (0 = Poisson)
  #' @return data.frame with added CI columns
  
  cat(sprintf("\n=== Calculating CIs for Resistance SNPs ===\n"))
  cat(sprintf("Using Gautier et al. 2013 Equation 4\n"))
  cat(sprintf("  Epsilon (genome-wide estimate): %.3f\n", epsilon))
  cat(sprintf("  Overdispersion parameter: %.1f\n", s_lambda))
  cat(sprintf("  Number of resistance SNPs: %d\n", nrow(data)))
  
  # Check required columns
  required_cols <- c("alt_count", "ref_count", "pool_size")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  
  # Optional columns for nice output
  id_cols <- intersect(c("gene", "mutation", "library"), colnames(data))
  
  # Calculate CIs
  results <- data %>%
    rowwise() %>%
    mutate(
      ci_result = list(calc_poolseq_ci_gautier( #from revised step 2! e.g. equation 4
        alt_count, 
        ref_count, 
        pool_size, 
        epsilon, 
        s_lambda
      ))
    ) %>%
    ungroup() %>%
    mutate(
      p = sapply(ci_result, function(x) x["p"]),
      SE = sapply(ci_result, function(x) x["SE"]),
      CI_lower = sapply(ci_result, function(x) x["CI_lower"]),
      CI_upper = sapply(ci_result, function(x) x["CI_upper"]),
      CI_width = CI_upper - CI_lower
    ) %>%
    select(-ci_result)
