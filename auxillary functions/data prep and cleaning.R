#setwd("~/Google Drive/My Drive/PhD_LSHTM/Projects/NECTAR/Meta-analysis")
#Data cleaning

#First we load the packages we need, read in the data, and we relabel the treatment arms in NECTAR2 which were incorrectly labelled, and convert it into a factor variable.

pacman::p_load(
  haven,          # import/export
  rio,            # import/export
  here,           # file pathways
  flextable,      # make HTML tables 
  officer,        # helper functions for tables
  tidyverse,      # data management, summary, and visualization
  patchwork)
library(ggbreak)
library(scales)
library(grid)
library(mgcv)
library(ggprism)
library(survival)
library(survminer)
library(pammtools)
library(lme4)
library(lmerTest)
library(modelbased)
library(netmeta)
library(blme)
library(emmeans)
library(BlandAltmanLeh)

source("functions_external_up.R")

nec <- read_dta("NECTAR_1234_Merged_Edit.dta")

nec$N1_2 = ifelse(is.na(nec$N1_2),0,nec$N1_2)

for(i in 1:nrow(nec)){
  if(nec$N1_2[i]==1 & nec$trt_arm[i]==2){
    nec$arm[i] = "DHA-PPQ + TQ (0.42 mg/kg)"
  } else {
    if(nec$N1_2[i]==1 & nec$trt_arm[i]==3){
      nec$arm[i] = "DHA-PPQ + TQ (0.83 mg/kg)"
    } else {
      if(nec$N1_2[i]==1 & nec$trt_arm[i]==4){
        nec$arm[i] = "DHA-PPQ + TQ (1.66 mg/kg)"
      } } }
  j=i
}


nec$tmpd = nec$mic_pf_asc_ul+nec$mic_pf_gct_ul
nec <- nec %>%
  mutate(arm = as.character(arm)) %>% # Convert arm to character
  mutate(arm = replace(arm, arm == "DHA-PPQ + PQ", "DHA-PPQ + PQ (0.25 mg/kg)"))
nec_merge<-nec[, c("studyid", "studyid_str", "study","studyvisit_num","uniqueid","visitdate","arm","age","sex_bin","temp","tmpd","mic_pf_gct_ul","qrtpcr_gct_ul","qrtpcr_femgct_ul","mosq_pos_prop","mic_pf_asc_ul","person_pos","year","mosq_total","mosq_pos","month","mean_ooc_all","mean_ooc_pos")]
nec_merge <- nec_merge %>%
  mutate(
    year = as.factor(year),
    arm = case_when(
      arm == "AL + PQ" ~ "AL + PQ (0.25 mg/kg)",
      arm == "AS-AQ + PQ" ~ "AS-AQ + PQ (0.25 mg/kg)",
      arm == "AL-AQ + PQ" ~ "AL-AQ + PQ (0.25 mg/kg)",
      arm == "PY-AS + PQ" ~ "PY-AS + PQ (0.25 mg/kg)",
      TRUE ~ arm
    ))




## Data cleaning PQ studies
pq01 <- read_dta("Primaquine efficacy trial _Wwarn dataset_June2016.dta")
pq01_asex <- read_dta("Mali PQ data wide_Teun gametocyte_june2016.dta")
pq01_gams <- read.csv("Mali raw data Dec 2014_hbinfectivity.csv")
pq03 <- read_dta("May17_2017-PQ03-DatabaseWide_EPI + QTPCR.dta")
pq03_gams <- read.csv("Repeat_PQ03_qPCR.csv")

pq01 <- pq01 %>%
  mutate(
    mosq_pos_prop = propinfected * 100,
    qrtpcr_femgct_ul = exp(qRTPCR_logconc),
    qrtpcr_prevalance = qRTPCR_prevalence,
    studyid_str = sprintf("PQ-01-%03d", id),
    uniqueid = paste0(sprintf("PQ-01-%03d", id), "_", round(visit, 1)),
    person_pos = ifelse(propinfected > 0, 1, 0)
  ) %>%
  rename(
    studyid = id,
    studyvisit_num = visit
  ) %>%
  mutate(
    arm = case_when(
      dose == 0.000 ~ "DHA-PPQ",
      dose == 0.125 ~ "DHA-PPQ + PQ (0.125 mg/kg)",
      dose == 0.500 ~ "DHA-PPQ + PQ (0.5 mg/kg)",
      dose == 0.250 ~ "DHA-PPQ + PQ (0.25 mg/kg)",
      dose == 0.0625 ~ "DHA-PPQ + PQ (0.0625 mg/kg)"
    ),
    is_baseline = if_else(studyvisit_num == 0, 1, 0)
  ) %>%
  group_by(studyid) %>%
  mutate(
    arm = if_else(is_baseline == 1, arm, NA_character_)
  ) %>%
  fill(arm, .direction = "downup") %>%
  mutate(
    age = if_else(is_baseline == 1, age, NA_real_),
    sex_bin = if_else(is_baseline == 1, gender, NA_real_)
  ) %>%
  fill(age, sex_bin, .direction = "downup") %>%
  ungroup() %>%
  select(-is_baseline) %>%
  mutate(
    arm = as.factor(arm),
    study = "PQ01",
    year = as.factor("2013-2014")
  )%>%
  filter(studyvisit_num %in% c(0,1,2,3,7,14,28))

#Merge mic_pf_gct_ul, mosq_total and mosq_pos
pq01_merge1 <- pq01_gams %>%
  filter(visit %in% c(0,1,2,3,7,14,28)) %>%
  mutate(studyid_str = sprintf("PQ-01-%03d", id),
         uniqueid = paste0(sprintf("PQ-01-%03d", id), "_", round(visit, 1)),) %>%
  select(uniqueid, numgams,nummostest,nummosqpos) %>%
  left_join(pq01, ., by = "uniqueid") %>%
  mutate(
    numgams = as.numeric(numgams) * 16,
    mic_pf_gct_ul = numgams,
    mosq_total = as.numeric(nummostest),
    mosq_pos = as.numeric(nummosqpos)
  )

#Merge asexual PC (only available at baseline)
pq01_merge2 <- pq01_asex %>%
  filter(studyvisit == "H0") %>%
  mutate(uniqueid = sprintf("PQ-01-%03d_0", id)) %>%
  select(uniqueid, asex_ul_d0, asex_prev_d0) %>%
  left_join(pq01_merge1, by = "uniqueid") %>%
  rename(mic_pf_asc_ul = asex_ul_d0)

#Merge again with PQ01_merge1, to include all other timepoints after baseline
pq01_merge3<-merge(pq01_merge1, pq01_merge2, by = c("studyid", "studyvisit_num", "propinfected", "hb", "age", "dose", "pqdosegrp",
                                                    "AE1", "AE2", "AE3", "AE4", "AE5", "AE6", "AE7", "AE1s", "AE2s", "AE3s", "AE4s", 
                                                    "AE5s", "AE6s", "AE7s", "anySAE", "gender", "qRTPCR_logconc", "qRTPCR_prevalence", 
                                                    "_merge", "mosq_pos_prop", "qrtpcr_femgct_ul", "qrtpcr_prevalance", "studyid_str", 
                                                    "uniqueid", "person_pos", "arm", "sex_bin", "study", "year","numgams","nummostest","nummosqpos","mic_pf_gct_ul","mosq_total","mosq_pos"), all.x = TRUE)

pq01_merge3<-pq01_merge3[, c("studyid", "studyid_str", "study","studyvisit_num","uniqueid","arm","age","sex_bin","mosq_pos_prop","person_pos","mic_pf_asc_ul","mic_pf_gct_ul","qrtpcr_femgct_ul","year","mosq_total","mosq_pos")]

pq01_merge3$qrtpcr_gct_ul = pq01_merge3$qrtpcr_femgct_ul

pq03_long <- pq03 %>%
  select(-starts_with("spo2")) %>%
  # Convert all columns (except the first 9 which are identifiers and should not be converted) to character
  mutate(across(-c(participantid, scrno, txgroup, anyae, status, agetoday, weight, height, group), as.character)) %>%
  pivot_longer(
    cols = -c(participantid, scrno, txgroup, anyae, status, agetoday, weight, height, group), # Specify columns to pivot
    names_to = c(".value", "timepoint"), # Split column names
    names_pattern = "^(.*?)(\\d+)$", # Pattern to extract variable name and timepoint
    values_to = "value" # Column to store original values
  ) %>% 
  rename(
    studyid = participantid,
    arm = txgroup,
    age = agetoday,
    mosq_pos_prop = propinfect,
    studyvisit_num = timepoint,
    mic_pf_asc_ul = pfasexual300wbc,
    mic_pf_gct_ul = pfgameto500wbc
  ) %>%
  mutate(
    studyid_str = sprintf("PQ-03-%03d", studyid),
    sex_bin = 1,
    person_pos = ifelse(mosq_pos_prop > 0, 1, 0),
    studyvisit_num = as.numeric(studyvisit_num),
    mosq_pos_prop = as.numeric(mosq_pos_prop) * 100,
    mic_pf_asc_ul = as.numeric(mic_pf_asc_ul) * 26.66666,
    mic_pf_gct_ul = as.numeric(mic_pf_gct_ul) * 16,
    temp = as.numeric(temp),
    mosq_total = as.numeric(nummosqsurv),
    mosq_pos = as.numeric(numoocystpres),
    uniqueid = paste0(studyid_str, "_", studyvisit_num),
    visitdate = as.Date(visitdate, format = "%Y-%m-%d"),
    study = "PQ03",
    year = as.factor("2016"),
    arm = case_when(
      arm == "DP" ~ "DHA-PPQ",
      arm == "DP +MB" ~ "DHA-PPQ + MB",
      arm == "SP-AQ +PQ" ~ "SP-AQ + PQ (0.25 mg/kg)",
      arm == "SP-AQ" ~ "SP-AQ"
    ),
    visitdate = as.Date(visitdate), 
    month = month(visitdate)
  )

pq03_reorder<-pq03_long[, c("studyid", "studyid_str","visitdate","sex_bin", "study","studyvisit_num","uniqueid","arm","age","temp","mosq_pos_prop","mic_pf_gct_ul","mic_pf_asc_ul","person_pos","year","mosq_total","mosq_pos","month")]


pq03_gams <- pq03_gams %>%
  mutate(qrtpcr_malegct_uL=malegam_mL/1000,
         qrtpcr_femgct_ul=femalegam_mL/1000,
         qrtpcr_gct_ul=totalgam_mL/1000)
pq03_gams_selected <- pq03_gams[, c("uniqueid", "qrtpcr_malegct_uL", "qrtpcr_femgct_ul", "qrtpcr_gct_ul", "remarks")]

#pq03_gams_selected <- pq03_gams_selected %>%
#mutate(qrtpcr_malegct_uL = ifelse(remarks != "", NA, qrtpcr_malegct_uL),
#qrtpcr_femgct_ul = ifelse(remarks != "", NA, qrtpcr_femgct_ul),
#qrtpcr_gct_ul = ifelse(remarks != "", NA, qrtpcr_gct_ul)
#)

pq03_merge <- merge(pq03_reorder, pq03_gams_selected, by = "uniqueid", all.x = TRUE)



nec <- bind_rows(nec_merge, pq01_merge3,pq03_merge)
nec <- nec %>%
  mutate(arm = factor(arm),
         tmpd = mic_pf_asc_ul + mic_pf_gct_ul,
         year = factor(year, levels = c("2013-2014", "2016", "2018", "2019", "2020", "2021", "2022")))


nec$asx = nec$mic_pf_asc_ul
nec$mictot = nec$mic_pf_asc_ul+nec$mic_pf_gct_ul
nec$gam.p = nec$qrtpcr_gct_ul
nec$gam.m = nec$mic_pf_gct_ul
nec$days = nec$studyvisit_num
nec$studyid_str = as.factor(nec$studyid_str)
nec$study = as.factor(nec$study)

nec$arm = recode(nec$arm,
                 "DHA-PPQ + PQ (0.25mg/kg)" = "ACT-PQ",
                 "PY-AS + PQ (0.25 mg/kg)" = "ACT-PQ",
                 "DHA-PPQ + TQ (0.42 mg/kg)" = "drop",
                 "DHA-PPQ + TQ (0.83 mg/kg)" = "ACT-TQ",
                 "DHA-PPQ + TQ (1.66 mg/kg)" = "ACT-TQ",
                 "AL + PQ (0.25 mg/kg)" = "ACT-PQ",
                 "SP-AQ + TQ 1.66mg/kg" = "drop",
                 "AL-AQ" = "AL",
                 "AL-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                 "AS-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                 "DHA-PPQ + PQ (0.0625 mg/kg)" = "drop",
                 "DHA-PPQ + PQ (0.125 mg/kg)" = "drop",
                 "DHA-PPQ + PQ (0.5 mg/kg)" = "ACT-PQ",
                 "DHA-PPQ + MB" = "drop",
                 "SP-AQ + PQ (0.25 mg/kg)" = "Non-ACT-PQ",
                 "DHA-PPQ + PQ (0.25 mg/kg)" = "ACT-PQ"
) 



nec = nec %>%
  filter(arm != "drop")

nec$arm = as.character(nec$arm)

#nec$arm = factor(nec$arm, levels=sort(unique(nec$arm)))

nec$arm = factor(nec$arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

#Exclude 5 participants without baseline infectivity measures
nec <- nec %>%
  filter(!studyid_str %in% c("PQ-01-001","PQ-01-004","PQ-01-005","PQ-01-006","PQ-01-008","PQ-01-052"))


nec = nec %>% arrange(studyid_str, studyvisit_num)


nec0 = nec[nec$studyvisit_num==0,]