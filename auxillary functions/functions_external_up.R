

# auxillary functions -----------------------------------------------------

dynamic_select <- function(data, ref_col, factor_levels) {
  # Ensure the levels are present in the dataset
  existing_levels <- intersect(factor_levels, colnames(data))
  
  # Select only the reference column and existing levels
  data %>% select(all_of(c(ref_col, existing_levels)))
}


# Survival models ---------------------------------------------------------

check.surv <- function(nma){

# Perform node-splitting to compare direct and indirect estimates
netsplit_results <- netsplit(nma)

# Extract direct and indirect treatment effect estimates
direct_vs_indirect <- data.frame(
  Comparison = netsplit_results$random$comparison,
  Direct = as.numeric(netsplit_results$direct.random$TE),
  Direct_Lower = as.numeric(netsplit_results$direct.random$lower),
  Direct_Upper = as.numeric(netsplit_results$direct.random$upper),
  Indirect = as.numeric(netsplit_results$indirect.random$TE),
  Indirect_Lower = as.numeric(netsplit_results$indirect.random$lower),
  Indirect_Upper = as.numeric(netsplit_results$indirect.random$upper)
)


# Create the forest plot
plot1 = ggplot(direct_vs_indirect, aes(x = Comparison)) +
  geom_point(aes(y = Direct, color = "Direct Estimate"), size = 3, alpha=0.7) +
  geom_errorbar(aes(ymin = Direct_Lower, ymax = Direct_Upper, color = "Direct Estimate"), width = 0.2) +
  geom_point(aes(y = Indirect, color = "Indirect Estimate"), size = 3, shape = 17, alpha=0.7) +
  geom_errorbar(aes(ymin = Indirect_Lower, ymax = Indirect_Upper, color = "Indirect Estimate"), width = 0.2, linetype = "dashed") +
  scale_color_manual(values = c("Direct Estimate" = "blue", "Indirect Estimate" = "red")) +
  theme_minimal() +
  labs(
    y = "log hazard ratio",
    x = "Treatment Comparison",
    color = "Estimate Type"
  ) +
  coord_flip()+
  theme(
    legend.position = "top"
  )


return(plot1)

}


# Reductions --------------------------------------------------------------


check.reduc.bet <- function(nma){
  # Perform node-splitting to compare direct and indirect estimates
  netsplit_results <- netsplit(nma)
  
  #day 0
  allcom = expand.grid(var1=levels(nec$arm),var2=levels(nec$arm))
  allcom = allcom[allcom$var1!=allcom$var2,]
  subcom = paste0(allcom$var1,":",allcom$var2)
  subcom = as.vector(subcom)
  suby = match(subcom, netsplit_results$random$comparison, nomatch = 0)
  
  # Extract direct and indirect treatment effect estimates
  direct_vs_indirect <- data.frame(
    Comparison = netsplit_results$random$comparison[suby],
    Direct = as.numeric(netsplit_results$direct.random$TE[suby]),
    Direct_se = as.numeric(netsplit_results$direct.random$seTE[suby]),
    Indirect = as.numeric(netsplit_results$indirect.random$TE[suby]),
    Indirect_se = as.numeric(netsplit_results$indirect.random$seTE[suby])
  ) %>%
    mutate(Comparison=factor(Comparison),
           Direct_Lower = Direct - 1.96*Direct_se,
           Direct_Upper = Direct + 1.96*Direct_se,
           Indirect_Lower = Indirect - 1.96*Indirect_se,
           Indirect_Upper = Indirect + 1.96*Indirect_se)
  
  # Create the forest plot
  plot2 = ggplot(direct_vs_indirect, aes(x = Comparison)) +
    geom_point(aes(y = Direct, color = "Direct Estimate"), size = 3, alpha=0.7) +
    geom_errorbar(aes(ymin = Direct_Lower, ymax = Direct_Upper, color = "Direct Estimate"), width = 0.2) +
    geom_point(aes(y = Indirect, color = "Indirect Estimate"), size = 3, shape = 17, alpha=0.7) +
    geom_errorbar(aes(ymin = Indirect_Lower, ymax = Indirect_Upper, color = "Indirect Estimate"), width = 0.2, linetype = "dashed") +
    scale_color_manual(values = c("Direct Estimate" = "blue", "Indirect Estimate" = "red")) +
    theme_minimal() +
    labs(
      y = "differences in outcome",
      x = "Treatment Comparison",
      color = "Estimate Type"
    ) +
    coord_flip(ylim=c(-100,100))+
    theme(
      legend.position = "top"
    )
  
  
  return(plot2)
  
}



