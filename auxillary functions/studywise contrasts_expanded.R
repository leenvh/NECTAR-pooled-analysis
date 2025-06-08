library(netmeta)
library(gridGraphics)
library(igraph)
library(ggraph)
library(ggtext)


get_studywise_contrasts_bam <- function(model, data, visitval, arms) {
  
  results.allc = results.all %>% mutate(arm = as.character(arm)) %>% arrange(arm, visit) %>% na.omit()
  
  # Ensure visit and arm are factors with correct levels
  visit_levels <- c("day 0", visitval)
  arm_levels <- levels(data$arm)
  
  valid_studyid_str <- data$studyid_str[1]  # Use first valid value for random effect
  
  # Get actual visit*arm combinations per study
  valid_combos <- data %>%
    mutate(visit_arm = paste(visit, arm, sep = "___")) %>%
    distinct(study, visit, arm, visit_arm)
  
  # Create all valid pairwise combinations within each study
  pairwise_contrasts <- valid_combos %>%
    group_by(study) %>%
    do({
      dat <- .
      expand.grid(arm_combo = dat$visit_arm, ref_combo = dat$visit_arm, stringsAsFactors = FALSE) %>%
        filter(arm_combo != ref_combo) %>%
        mutate(study = unique(dat$study))  # Fix scoping for mutate
    }) %>%
    ungroup() %>%
    separate(arm_combo, into = c("visit_arm_visit", "visit_arm_arm"), sep = "___", extra = "merge") %>%
    separate(ref_combo, into = c("ref_arm_visit", "ref_arm_arm"), sep = "___", extra = "merge") %>%
    mutate(
      visit_arm_visit_key = paste(visit_arm_visit, visit_arm_arm, sep = "___"),
      ref_arm_visit_key = paste(ref_arm_visit, ref_arm_arm, sep = "___")
    )
  
  # Prepare prediction dataset with valid factor levels and studyid_str
  pred_data <- valid_combos %>%
    mutate(
      visit = factor(visit, levels = visit_levels),
      arm = factor(arm, levels = arm_levels),
      studyid_str = valid_studyid_str
    ) %>% na.omit()
  
  # Predict fitted values for each visit-arm-study combo
  preds <- predict(model, newdata = pred_data, se.fit = TRUE, exclude = "s(studyid_str)")
  
  pred_data$fit <- preds$fit
  pred_data$se_fit <- preds$se.fit
  pred_data$visit_arm <- paste(pred_data$visit, pred_data$arm, sep = "___")
  
  # Merge predictions for arm and ref correctly
  results <- pairwise_contrasts %>%
    left_join(pred_data %>% select(study, visit_arm, fit, se_fit) %>%
                rename(fit_arm = fit, se_fit_arm = se_fit),
              by = c("study", "visit_arm_visit_key" = "visit_arm")) %>%
    left_join(pred_data %>% select(study, visit_arm, fit, se_fit) %>%
                rename(fit_ref = fit, se_fit_ref = se_fit),
              by = c("study", "ref_arm_visit_key" = "visit_arm")) %>%
    mutate(
      arm = paste(visit_arm_visit, visit_arm_arm),
      ref = paste(ref_arm_visit, ref_arm_arm),
      coef = fit_arm - fit_ref,
      se = sqrt(se_fit_arm^2 + se_fit_ref^2)
    ) %>%
    select(study, arm, ref, coef, se)
  
  
  # Filter to keep only comparisons involving Day 0 and Day 2
results_day0_day2 <- results %>%
  filter(
    (grepl("day 0", arm) | grepl(visitval, arm)) &
      (grepl("day 0", ref) | grepl(visitval, ref))
  )

results_day0_day2_aligned <- results_day0_day2 %>%
  mutate(
    swap_flag = grepl("day 0", arm) & grepl(visitval, ref),
    
    arm_final = ifelse(swap_flag, ref, arm),
    ref_final = ifelse(swap_flag, arm, ref),
    
    coef_final = ifelse(swap_flag, -coef, coef),
    se_final = se
  )



results_day0_day2_aligned <- results_day0_day2_aligned %>%
  mutate(
    arm_clean = gsub("day [0-9]+ ", "", arm_final),
    ref_clean = gsub("day [0-9]+ ", "", ref_final)
  )


# Filter to keep only within-treatment comparisons: Day 2 vs Day 0 of the same treatment
results_day0_day2_final <- results_day0_day2_aligned %>%
  filter(
    grepl(visitval, arm_final) & grepl("day 0", ref_final),  # Ensure correct days
    arm_clean == ref_clean  # Ensure same treatment
  )


# Step 1: Remove visitval and "day 0"
results_clean <- results_day0_day2_final %>%
  mutate(
    arm = gsub("day [0-9]+ ", "", arm_final),
    ref = gsub("day [0-9]+ ", "", ref_final),
    
    # Step 2: Compute relative reduction (1 - exp(beta)) * 100%
    rel_reduction = (1 - exp(coef_final)) * 100,
    
    # Standard error of exp(beta) using Delta method: SE(exp(beta)) = exp(beta) * SE(beta)
    se_rel_reduction = exp(coef_final) * se_final * 100
  ) %>% 
  select(study, arm_clean, rel_reduction, se_rel_reduction) %>%
  rename(arm=arm_clean, coef = rel_reduction, se= se_rel_reduction)

results_clean = results_clean[!duplicated(results_clean),]



# Step 3
# Create all pairwise comparisons within each study
pairwise_reductions <- results_clean %>%
  group_by(study) %>%
  do({
    data = .
    expand.grid(ref = data$arm, arm = data$arm, stringsAsFactors = FALSE) %>%
      filter(ref != arm) %>%
      left_join(data, by = c("ref" = "arm")) %>%
      rename(ref_rel_red = coef, ref_se_rel_red = se) %>%
      left_join(data, by = c("arm" = "arm")) %>%
      rename(arm_rel_red = coef, arm_se_rel_red = se)
  }) %>%
  ungroup() %>%
  mutate(
    # Calculate arm - ref (trust the labels!)
    coef = arm_rel_red - ref_rel_red,
    
    # SE calculation remains standard
    se = sqrt(ref_se_rel_red^2 + arm_se_rel_red^2)
  ) %>%
  select(study, arm, ref, coef, se)  %>%
  mutate(
    treat_min = pmin(arm, ref),
    treat_max = pmax(arm, ref),
    
    # If we flipped to min/max, adjust coef accordingly
    coef_adj = ifelse(arm == treat_max, coef, -coef)
  ) %>%
  group_by(study, treat_min, treat_max) %>%
  slice(1) %>%
  ungroup() %>%
  select(study, arm = treat_max, ref = treat_min, coef = coef_adj, se)



nma = netmeta(TE = pairwise_reductions$coef, seTE=pairwise_reductions$se, treat1=pairwise_reductions$arm, treat2=pairwise_reductions$ref, studlab=pairwise_reductions$study, sm="MD", reference="AL")
sumnma=summary(nma)$random

results.allc = na.omit(results.allc %>% filter(visit %in% visitval))

#day 2 reductions
sumnmadays = list()
sumnmadays$TE = results.allc$RR[results.allc$visit==visitval]
sumnmadays$seTE = results.allc$se[results.allc$visit==visitval]
sumnmadays$lower = results.allc$ci_lower[results.allc$visit==visitval]
sumnmadays$upper = results.allc$ci_upper[results.allc$visit==visitval]
sumnmadays$pval = results.allc$p_value[results.allc$visit==visitval]
outresdays = matrix(NA, nrow=length(arms),ncol=length(arms))


colnames(outresdays)=results.allc$arm[results.allc$visit==visitval]
row.names(outresdays)=results.allc$arm[results.allc$visit==visitval]
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmadays$TE)) {
  for (j in seq_along(sumnmadays$TE)) {
    if (i != j) {
      # Calculate z-score
      
      # Calculate two-sided p-value
      p <- t(sumnma$p[i,j])
      
      p_value <- paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
      
      difference <- format(round((t(sumnma$TE)[i,j]), 2), nsmall = 2)
      lower <- format(round((t(sumnma$lower)[i,j]), 2), nsmall = 2)
      upper <- format(round((t(sumnma$upper)[i,j]), 2), nsmall = 2)
      
      outresdays[i, j] = paste0(difference, "% (", lower, "%",", ",upper,"%), ", p_value)
    }
  }
}



diag(outresdays) = sapply(seq_along(sumnmadays$TE), function(i) {
  
  # Get the two-sided p-value
  p <- sumnmadays$pval[i]
  
  # Format the percentage reduction and confidence interval
  reduction <- format(round((sumnmadays$TE[i]), 2), nsmall = 2)
  lower <- format(round((sumnmadays$lower[i]), 2), nsmall = 2)
  upper <- format(round((sumnmadays$upper[i]), 2), nsmall = 2)
  
  # Format the p-value
  p_value <- ifelse(p < 0.0001, "<0.0001", ifelse(p > 0.9999, ">0.9999", paste0("=", format(round(p, 4), scientific=FALSE, nsmall = 4))))
  
  # Combine the reduction, confidence interval, and p-value
  paste0(reduction, "% (", lower, "%",", ",upper,"%), p", p_value)
})


result_matrixdays <- as.data.frame(outresdays)
result_matrixdays$Reference = factor(row.names(outresdays),levels=arms)
result_matrixdays = result_matrixdays %>% select("Reference",all_of(arms)) %>% arrange(Reference)

ranks <- netrank(nma, small.values = "bad")
rdata <- data.frame(arm=names(ranks$ranking.random), pscore =ranks$ranking.random) %>%
  arrange(pscore) 
rdata$arm = factor(rdata$arm, levels = rdata$arm)

pscoreplot = ggplot(data=rdata, aes(x=visitval,y=arm, fill=pscore))+
  geom_tile()+
  theme_prism()+
  theme(axis.title = element_blank(),
        legend.title = element_text(),
        legend.position = "top")+
  geom_label(aes(label=format(round(pscore,2),nsmall=2)),col="black")+
  scale_fill_gradient2(
    low = "#d9e4e6", high = "#2b6777",
    limits = c(0, 1),
    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    name = "P-score\n"
  )

outcome_var <- all.vars(formula(model))[1]

summation = data %>% 
  filter(visit %in% visitval & !is.na(outcome_var)) %>%
  group_by(arm) %>%
  summarise(patients = n()) %>%
  mutate(arm = as.character(arm)) %>% 
  rename(name=arm) %>%
  ungroup()

study_df <- pairwise_reductions %>%
  mutate(pair = paste(pmin(arm, ref), pmax(arm, ref), sep = " vs ")) %>%
  count(pair, treat1 = pmin(arm, ref), treat2 = pmax(arm, ref)) %>%
  rename(comparisons = n) %>%
  select(-pair)


# Build igraph object
g <- graph_from_data_frame(d = study_df, directed = FALSE, vertices = summation)

# Plot with ggraph
networkplot = ggraph(g, layout = 'fr') +
  geom_edge_link(aes(width = comparisons), color = "black", alpha = 0.7) +
  geom_node_point(aes(size = patients), color = "#2b6777", alpha = 0.9) +
  geom_node_text(aes(label = name), vjust = 1.5, size = 5) +
  scale_size_continuous(limits=c(0,200),range = c(3, 12), name = "Number of Patients") +
  scale_edge_width(range = c(0.5, 4),breaks=c(1,2), name = "Number of Comparisons") +
  theme_void()+
  theme(legend.position = "right",
        plot.margin = unit(c(1, 1, 1, 1), "cm"))+
ggtitle(visitval)




my_base_plot <- function() forest(nma, xlab=visitval)
# Run and capture as a grob
my_base_plot()
grob_object <- grid.grab()





outlist= list(check.reduc.bet(nma),result_matrixdays, pscoreplot, grob_object, networkplot, rdata)

  return(outlist)
}
