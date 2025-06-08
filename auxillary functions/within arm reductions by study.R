

get_reductions = function(model){

newd=modelbased::estimate_contrasts(model, contrast=c("visit"), by=c("study","arm"))
newd = newd %>% select(study,arm,Level1,Level2,Difference,SE)
colnames(newd) = c("study","arm","compare","ref","coef","se")

# Step 1: Keep only rows where Day 0 is in compare or ref
newd_filtered <- newd %>%
  filter(grepl("day 0", compare) | grepl("day 0", ref))

# Step 2: Swap if Day 0 is in compare, and flip coef sign
newd_aligned <- newd_filtered %>%
  mutate(
    swap_flag = grepl("day 0", compare),
    compare_final = ifelse(swap_flag, ref, compare),
    ref_final = ifelse(swap_flag, compare, ref),
    coef_final = ifelse(swap_flag, -coef, coef),
    se_final = se
  )

# Step 3: Calculate relative reductions and 95% CIs
newd_results <- newd_aligned %>%
  mutate(
    RR = (1 - exp(coef_final)) * 100,
    ci_lower = (1 - exp(coef_final + 1.96 * se_final)) * 100,
    ci_upper = (1 - exp(coef_final - 1.96 * se_final)) * 100,
    
    # Step 4: Calculate p-value using Wald test
    z_value = coef_final / se_final,
    p_value = 2 * (1 - pnorm(abs(z_value))),
    se=se_final
  ) %>%
  select(study,arm, visit = compare_final, RR, se, ci_lower, ci_upper, p_value)

return(newd_results)

}
