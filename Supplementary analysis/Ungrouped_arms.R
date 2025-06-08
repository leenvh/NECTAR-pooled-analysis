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
         "DHA-PPQ + TQ (0.42 mg/kg)" = "drop",
         "SP-AQ + TQ 1.66mg/kg" = "drop",
         "DHA-PPQ + PQ (0.0625 mg/kg)" = "drop",
         "DHA-PPQ + PQ (0.125 mg/kg)" = "drop",
         "DHA-PPQ + MB" = "drop"
  ) 



nec = nec %>%
  filter(arm != "drop")

nec$arm2 = as.character(nec$arm)
nec$arm2 = as.factor(nec$arm2)

nec$arm2 = recode(nec$arm2,
                  "DHA-PPQ + PQ (0.25mg/kg)" = "ACT-PQ",
                  "PY-AS + PQ (0.25 mg/kg)" = "ACT-PQ",
                  "DHA-PPQ + TQ (0.83 mg/kg)" = "ACT-TQ",
                  "DHA-PPQ + TQ (1.66 mg/kg)" = "ACT-TQ",
                  "AL + PQ (0.25 mg/kg)" = "ACT-PQ",
                  "AL-AQ" = "AL",
                  "AL-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                  "AS-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                  "DHA-PPQ + PQ (0.5 mg/kg)" = "ACT-PQ",
                  "SP-AQ + PQ (0.25 mg/kg)" = "Non-ACT-PQ",
                  "DHA-PPQ + PQ (0.25 mg/kg)" = "ACT-PQ"
) 




nec$arm2 = factor(nec$arm2, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

nec$arm = as.character(nec$arm)

nec$arm = factor(nec$arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", sort(unique(nec$arm[nec$arm2 %in% "AL"])), sort(unique(nec$arm[nec$arm2 %in% "Non-ACT-PQ"])), sort(unique(nec$arm[nec$arm2 %in% "ACT-PQ"])), sort(unique(nec$arm[nec$arm2 %in% "ACT-TQ"]))))

armlevels = levels(nec$arm)

#Exclude 5 participants without baseline infectivity measures
nec <- nec %>%
  filter(!studyid_str %in% c("PQ-01-001","PQ-01-004","PQ-01-005","PQ-01-006","PQ-01-008","PQ-01-052"))


nec = nec %>% arrange(studyid_str, studyvisit_num)

armlevs = levels(nec$arm)




# Functions required ------------------------------------------------------

source("within arm reductions.R")
source("studywise contrasts_expanded.R")


# Descriptive results -----------------------------------------------------



#We create a dataset for the baseline data [i.e. at visit 0] and make a summarising function for the descriptives.

nec0 = nec[nec$studyvisit_num==0,]

# Treatment arm differences -----------------------------------------------


# Time to clearance -------------------------------------------------------

#Here we will look at the time to clearance of gametocytes and the time to clearance of infectivity. We can estimate the Kaplan-Meier survival curves and perhaps compare them. We can also look at the Cox model for the treatment effects hazard ratios.

####gametocytes by pcr####

necgam = nec %>% 
  select(study, studyid_str, arm, days, gam.m, gam.p) %>% 
  filter(days>0) %>%
  left_join(
    nec %>% 
      select(studyid_str,days, gam.m, gam.p) %>% 
      filter(days==0) %>%
      mutate(gam0.m = gam.m, gam0.p=gam.p) %>%
      select(studyid_str, gam0.m, gam0.p),
    by="studyid_str"
  )


nectime5a = nec %>% 
  #mutate(gam.p = ifelse(is.na(gam.p),gam.m, gam.p)) %>%
  filter(!is.na(gam.p)) %>%
  select(studyid_str, study, days, gam.p, gam.m,mosq_total, mosq_pos,arm) %>%
  group_by(studyid_str, study) %>%
  mutate(firstzero = cumsum(I(gam.p==0)), status=ifelse(firstzero!=1,0,1), anyevent=sum(status), lastv = n():1, lastv=ifelse(lastv==1 & sum(status)==0,1,0)) %>%
  select(studyid_str, study, days, gam.p, gam.m,mosq_total, mosq_pos, status, lastv,arm) %>%
  filter((status==1 & lastv==0) | (status==0 & lastv==1)) %>%
  filter(days>0) %>% 
  mutate(arm = as.character(arm))

nectime5a$arm = factor(nectime5a$arm, levels=armlevels)


nectime5a = nectime5a %>% left_join(
  necgam %>% 
    select(studyid_str, study, gam0.m, gam0.p) %>%
    distinct(), by=c("study","studyid_str")
)

studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
currentstudy = studies[i]

subdata = nectime5a[nectime5a$study %in% currentstudy,]
subdata$arm = as.factor(as.character(subdata$arm))
mod5a = coxph(Surv(days,status)~arm+I(log(gam0.p+0.01)), data= subdata)

vcov_mat <- vcov(mod5a)
vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
mod5a$var = vcov_mat

sumdat=as.data.frame(summary(mod5a)$coefficients)
sumdat = sumdat[-length(levels(subdata$arm)),]

newd = expand.grid(study=currentstudy, arm=levels(subdata$arm)[-1], ref=levels(subdata$arm)[1])
newd$arm = as.character(newd$arm)
newd$ref = as.character(newd$ref)
newd$coef = sumdat$coef
newd$se = sumdat$`se(coef)`

if(length(levels(subdata$arm))>2){
pairs = combn(levels(subdata$arm)[-1],2)

for(j in 1:ncol(pairs)){
  cpair = pairs[,j]
  newd[nrow(newd)+1,]=NA
  newd[nrow(newd),"arm"] = cpair[1]
  newd[nrow(newd),"ref"] = cpair[2]
  newd[nrow(newd),"study"] = currentstudy
  newd[nrow(newd),"coef"] = mod5a$coefficients[paste0("arm",cpair[1])] - mod5a$coefficients[paste0("arm",cpair[2])]
  newd[nrow(newd),"se"] = sqrt(vcov(mod5a)[paste0("arm",cpair[1]),paste0("arm",cpair[1])]-2*vcov(mod5a)[paste0("arm",cpair[1]),paste0("arm",cpair[2])]+vcov(mod5a)[paste0("arm",cpair[2]),paste0("arm",cpair[2])])
  }
}

results = results %>% rbind(newd) %>% na.omit()

}


results$arm = factor(results$arm, levels=armlevels)

library(netmeta)
nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="HR", reference="DHA-PPQ")

sumnma=summary(nma)$random

outres = matrix(paste0(format(round(exp(sumnma$TE),2),nsmall=2)," (",format(round(exp(sumnma$lower),2), nsmall=2)," - ",format(round(exp(sumnma$upper),2), nsmall=2),") p",ifelse(sumnma$p<0.0001,"<0.0001",ifelse(sumnma$p>0.9999,">0.9999",paste0("=",format(round(sumnma$p,4), nsmall=4))))), nrow=nrow(sumnma$TE), ncol=ncol(sumnma$TE))

colnames(outres)=colnames(sumnma$TE)
row.names(outres)=row.names(sumnma$TE)
diag(outres)=""
outres = t(outres)


netgraph(
  nma,
  number = FALSE,                    # Add node numbers for easier reference
  thickness = "se.fixed",             # Use line thickness inversely proportional to standard errors
  multiarm = F,                  # Adjust for multi-arm studies
  cex = 1.2,                        # Scale text for better readability
  cex.points = 1.5,                 # Scale node sizes
  lty = c("direct" = 1, "indirect" = 6),  # Solid lines for direct, dashed for indirect evidence
  points = TRUE                     # Show treatment nodes as points
)


#pfor=list()

# for(arm in levels(nectime5a$arm)){
#   pfor[[arm]] = forest(nma, pooled="random", reference.group = arm)
# }

comp=netsplit(nma)
indirect=comp$comparison[comp$prop.random<=0.02]
# Separate the 'indirect' comparisons into rows and columns
indirect_pairs <- as.data.frame(do.call(rbind, strsplit(indirect, ":")), stringsAsFactors = FALSE)
colnames(indirect_pairs) <- c("var1", "var2")  # Rename columns for clarity
indirect_pairs = rbind(indirect_pairs, indirect_pairs %>% rename(var1=var2,var2=var1) %>% select(var1, var2))

for (i in 1:nrow(indirect_pairs)){
outres[indirect_pairs$var1[i],indirect_pairs$var2[i]] = paste0("\U01D43",outres[indirect_pairs$var1[i],indirect_pairs$var2[i]])
}


result_matrix5a = as.data.frame(outres)
result_matrix5a$Reference = factor(row.names(outres),levels=armlevels)
result_matrix5a = result_matrix5a %>% dynamic_select("Reference",armlevels) %>% arrange(Reference)


flextable(result_matrix5a) %>% autofit()  %>% save_as_docx( path = "OUTPUT/long_drop_some_up/survival_gampcr_NMA_all_some_drop.docx")



####propinfected adjusted baseline gam by pcr####

nectime5b = nec %>% 
  mutate(prop=mosq_pos/mosq_total) %>%
  filter(!is.na(prop) & days>0) %>%
  select(studyid_str, study, days, gam.p, gam.m,mosq_total, mosq_pos,prop,arm) %>%
  group_by(studyid_str, study) %>%
  mutate(firstzero = cumsum(I(prop==0)), status=ifelse(firstzero!=1,0,1), anyevent=sum(status), lastv = n():1, lastv=ifelse(lastv==1 & sum(status)==0,1,0)) %>%
  select(studyid_str, study, days, gam.p, gam.m,mosq_total, mosq_pos,prop, status, lastv,arm) %>%
  filter((status==1 & lastv==0) | (status==0 & lastv==1))%>%
  filter(days>0) %>% 
  mutate(arm = as.character(arm))

nectime5b$arm = factor(nectime5b$arm, levels=armlevels)

nectime5b = nectime5b %>% left_join(
  necgam %>% 
    select(studyid_str, study, gam0.m, gam0.p) %>%
    distinct(), by=c("study","studyid_str")
)


studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = nectime5b[nectime5b$study %in% currentstudy,]
  subdata$arm = as.factor(as.character(subdata$arm))
  mod5b = coxph(Surv(days,status)~arm+I(log(gam0.p+0.01)), data= subdata)
  vcov_mat <- vcov(mod5b)
  vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
  vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
  mod5b$var = vcov_mat
  sumdat=as.data.frame(summary(mod5b)$coefficients)
  sumdat = sumdat[-length(levels(subdata$arm)),]
  
  newd = expand.grid(study=currentstudy, arm=levels(subdata$arm)[-1], ref=levels(subdata$arm)[1])
  newd$arm = as.character(newd$arm)
  newd$ref = as.character(newd$ref)
  newd$coef = sumdat$coef
  newd$se = sumdat$`se(coef)`
  
  if(length(levels(subdata$arm))>2){
    pairs = combn(levels(subdata$arm)[-1],2)
    
    for(j in 1:ncol(pairs)){
      cpair = pairs[,j]
      newd[nrow(newd)+1,]=NA
      newd[nrow(newd),"arm"] = cpair[1]
      newd[nrow(newd),"ref"] = cpair[2]
      newd[nrow(newd),"study"] = currentstudy
      newd[nrow(newd),"coef"] = mod5b$coefficients[paste0("arm",cpair[1])] - mod5b$coefficients[paste0("arm",cpair[2])]
      newd[nrow(newd),"se"] = sqrt(vcov(mod5b)[paste0("arm",cpair[1]),paste0("arm",cpair[1])]-2*vcov(mod5b)[paste0("arm",cpair[1]),paste0("arm",cpair[2])]+vcov(mod5b)[paste0("arm",cpair[2]),paste0("arm",cpair[2])])
    }
  }
  
  results = results %>% rbind(newd) %>% na.omit()
  
}


results$arm = factor(results$arm, levels=armlevels)

library(netmeta)
nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="HR", reference="DHA-PPQ")

sumnma=summary(nma)$random

outres = matrix(paste0(format(round(exp(sumnma$TE),2),nsmall=2)," (",format(round(exp(sumnma$lower),2), nsmall=2)," - ",format(round(exp(sumnma$upper),2), nsmall=2),") p",ifelse(sumnma$p<0.0001,"<0.0001",ifelse(sumnma$p>0.9999,">0.9999",paste0("=",format(round(sumnma$p,4), nsmall=4))))), nrow=nrow(sumnma$TE), ncol=ncol(sumnma$TE))

colnames(outres)=colnames(sumnma$TE)
row.names(outres)=row.names(sumnma$TE)
diag(outres)=""
outres = t(outres)


netgraph(
  nma,
  number = FALSE,                    # Add node numbers for easier reference
  thickness = "se.fixed",             # Use line thickness inversely proportional to standard errors
  multiarm = F,                  # Adjust for multi-arm studies
  cex = 1.2,                        # Scale text for better readability
  cex.points = 1.5,                 # Scale node sizes
  lty = c("direct" = 1, "indirect" = 6),  # Solid lines for direct, dashed for indirect evidence
  points = TRUE                     # Show treatment nodes as points
)


#pfor=list()

# for(arm in levels(nectime5b$arm)){
#   pfor[[arm]] = forest(nma, pooled="random", reference.group = arm)
# }

comp=netsplit(nma)
indirect=comp$comparison[comp$prop.random<=0.02]
# Separate the 'indirect' comparisons into rows and columns
indirect_pairs <- as.data.frame(do.call(rbind, strsplit(indirect, ":")), stringsAsFactors = FALSE)
colnames(indirect_pairs) <- c("var1", "var2")  # Rename columns for clarity
indirect_pairs = rbind(indirect_pairs, indirect_pairs %>% rename(var1=var2,var2=var1) %>% select(var1, var2))

for (i in 1:nrow(indirect_pairs)){
  outres[indirect_pairs$var1[i],indirect_pairs$var2[i]] = paste0("\U01D43",outres[indirect_pairs$var1[i],indirect_pairs$var2[i]])
}


result_matrix5b = as.data.frame(outres)
result_matrix5b$Reference = factor(row.names(outres),levels=armlevels)
result_matrix5b = result_matrix5b %>% dynamic_select("Reference",armlevels) %>% arrange(Reference)


flextable(result_matrix5b) %>% autofit()  %>% save_as_docx( path = "OUTPUT/long_drop_some_up/survival_infectivity_adjpcr_NMA_all_some_drop.docx")



####gametocytes by mic####

nectime5c = nec %>% 
  #mutate(gam.m = ifelse(is.na(gam.m),gam.p, gam.m)) %>%
  filter(!is.na(gam.m)) %>%
  select(studyid_str, study, days, gam.m, gam.p,mosq_total, mosq_pos,arm) %>%
  group_by(studyid_str, study) %>%
  mutate(firstzero = cumsum(I(gam.m==0)), status=ifelse(firstzero!=1,0,1), anyevent=sum(status), lastv = n():1, lastv=ifelse(lastv==1 & sum(status)==0,1,0)) %>%
  select(studyid_str, study, days, gam.m, gam.p,mosq_total, mosq_pos, status, lastv,arm) %>%
  filter((status==1 & lastv==0) | (status==0 & lastv==1))%>%
  filter(days>0) %>% 
  mutate(arm = as.character(arm))

nectime5c$arm = factor(nectime5c$arm, levels=armlevels)



nectime5c = nectime5c %>% left_join(
  necgam %>% 
    select(studyid_str, study, gam0.p, gam0.m) %>%
    distinct(), by=c("studyid_str","study")
)


studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = nectime5c[nectime5c$study %in% currentstudy,]
  
  if(nrow(subdata)>0){
  
  subdata$arm = as.factor(as.character(subdata$arm))
  mod5c = coxph(Surv(days,status)~arm+I(log(gam0.m+0.1)), data= subdata)
  vcov_mat <- vcov(mod5c)
  vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
  vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
  mod5c$var = vcov_mat
  sumdat=as.data.frame(summary(mod5c)$coefficients)
  sumdat = sumdat[-length(levels(subdata$arm)),]
  
  newd = expand.grid(study=currentstudy, arm=levels(subdata$arm)[-1], ref=levels(subdata$arm)[1])
  newd$arm = as.character(newd$arm)
  newd$ref = as.character(newd$ref)
  newd$coef = sumdat$coef
  newd$se = sumdat$`se(coef)`
  
  if(length(levels(subdata$arm))>2){
    pairs = combn(levels(subdata$arm)[-1],2)
    
    for(j in 1:ncol(pairs)){
      cpair = pairs[,j]
      newd[nrow(newd)+1,]=NA
      newd[nrow(newd),"arm"] = cpair[1]
      newd[nrow(newd),"ref"] = cpair[2]
      newd[nrow(newd),"study"] = currentstudy
      newd[nrow(newd),"coef"] = mod5c$coefficients[paste0("arm",cpair[1])] - mod5c$coefficients[paste0("arm",cpair[2])]
      newd[nrow(newd),"se"] = sqrt(vcov(mod5c)[paste0("arm",cpair[1]),paste0("arm",cpair[1])]-2*vcov(mod5c)[paste0("arm",cpair[1]),paste0("arm",cpair[2])]+vcov(mod5c)[paste0("arm",cpair[2]),paste0("arm",cpair[2])])
    }
  }
  
  results = results %>% rbind(newd) %>% na.omit()
  
}

}

results$arm = factor(results$arm, levels=armlevels)

library(netmeta)
nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="HR", reference="DHA-PPQ")

sumnma=summary(nma)$random

outres = matrix(paste0(format(round(exp(sumnma$TE),2),nsmall=2)," (",format(round(exp(sumnma$lower),2), nsmall=2)," - ",format(round(exp(sumnma$upper),2), nsmall=2),") p",ifelse(sumnma$p<0.0001,"<0.0001",ifelse(sumnma$p>0.9999,">0.9999",paste0("=",format(round(sumnma$p,4), nsmall=4))))), nrow=nrow(sumnma$TE), ncol=ncol(sumnma$TE))

colnames(outres)=colnames(sumnma$TE)
row.names(outres)=row.names(sumnma$TE)
diag(outres)=""
outres = t(outres)


netgraph(
  nma,
  thickness = "number.of.studies",             # Use line thickness inversely proportional to standard errors
  cex = 1.2,                        # Scale text for better readability
  cex.points = 1.5,                 # Scale node sizes
  lty = c("direct" = 1, "indirect" = 6),  # Solid lines for direct, dashed for indirect evidence
  points = TRUE,                     # Show treatment nodes as points, start = "circle", iterate = TRUE
  multiarm = TRUE, col.multiarm=colorspace::rainbow_hcl(n=5)
)


#pfor=list()

# for(arm in levels(nectime5c$arm)){
#   pfor[[arm]] = forest(nma, pooled="random", reference.group = arm)
# }

comp=netsplit(nma)
indirect=comp$comparison[comp$prop.random<=0.02]
# Separate the 'indirect' comparisons into rows and columns
indirect_pairs <- as.data.frame(do.call(rbind, strsplit(indirect, ":")), stringsAsFactors = FALSE)
colnames(indirect_pairs) <- c("var1", "var2")  # Rename columns for clarity
indirect_pairs = rbind(indirect_pairs, indirect_pairs %>% rename(var1=var2,var2=var1) %>% select(var1, var2))

for (i in 1:nrow(indirect_pairs)){
  outres[indirect_pairs$var1[i],indirect_pairs$var2[i]] = paste0("\U01D43",outres[indirect_pairs$var1[i],indirect_pairs$var2[i]])
}


result_matrix5c = as.data.frame(outres)
result_matrix5c$Reference = factor(row.names(outres),levels=armlevels)
result_matrix5c = result_matrix5c %>% dynamic_select("Reference",armlevels) %>% arrange(Reference)

flextable(result_matrix5c) %>% autofit()  %>% save_as_docx( path = "OUTPUT/long_drop_some_up/survival_gammic_NMA_all_some_drop.docx")


####propinfected adjusted baseline gam by mic####

nectime5d = nec %>% 
  mutate(prop=mosq_pos/mosq_total) %>%
  filter(!is.na(prop) & days>0) %>%
  select(studyid_str, study, days, gam.p, gam.m,mosq_total, mosq_pos,prop,arm) %>%
  group_by(studyid_str, study) %>%
  mutate(firstzero = cumsum(I(prop==0)), status=ifelse(firstzero!=1,0,1), anyevent=sum(status), lastv = n():1, lastv=ifelse(lastv==1 & sum(status)==0,1,0)) %>%
  select(studyid_str, study, days, gam.p, gam.m,mosq_total, mosq_pos,prop, status, lastv,arm) %>%
  filter((status==1 & lastv==0) | (status==0 & lastv==1))%>%
  filter(days>0) %>% 
  mutate(arm = as.character(arm))

nectime5d$arm = factor(nectime5d$arm, levels=armlevels)


nectime5d = nectime5d %>% left_join(
  necgam %>% 
    select(studyid_str, study, gam0.m, gam0.p) %>%
    distinct(), by=c("study","studyid_str")
)

studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = nectime5d[nectime5d$study %in% currentstudy,]
  
  if(nrow(subdata)>0){
    
    subdata$arm = as.factor(as.character(subdata$arm))
    mod5d = coxph(Surv(days,status)~arm+I(log(gam0.m+0.1)), data= subdata)
    vcov_mat <- vcov(mod5d)
    vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
    vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
    mod5d$var = vcov_mat
    sumdat=as.data.frame(summary(mod5d)$coefficients)
    sumdat = sumdat[-length(levels(subdata$arm)),]
    
    newd = expand.grid(study=currentstudy, arm=levels(subdata$arm)[-1], ref=levels(subdata$arm)[1])
    newd$arm = as.character(newd$arm)
    newd$ref = as.character(newd$ref)
    newd$coef = sumdat$coef
    newd$se = sumdat$`se(coef)`
    
    if(length(levels(subdata$arm))>2){
      pairs = combn(levels(subdata$arm)[-1],2)
      
      for(j in 1:ncol(pairs)){
        cpair = pairs[,j]
        newd[nrow(newd)+1,]=NA
        newd[nrow(newd),"arm"] = cpair[1]
        newd[nrow(newd),"ref"] = cpair[2]
        newd[nrow(newd),"study"] = currentstudy
        newd[nrow(newd),"coef"] = mod5d$coefficients[paste0("arm",cpair[1])] - mod5d$coefficients[paste0("arm",cpair[2])]
        newd[nrow(newd),"se"] = sqrt(vcov(mod5d)[paste0("arm",cpair[1]),paste0("arm",cpair[1])]-2*vcov(mod5d)[paste0("arm",cpair[1]),paste0("arm",cpair[2])]+vcov(mod5d)[paste0("arm",cpair[2]),paste0("arm",cpair[2])])
      }
    }
    
    results = results %>% rbind(newd) %>% na.omit()
    
  }
  
}

results$arm = factor(results$arm, levels=armlevels)

nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="HR", reference="DHA-PPQ")

sumnma=summary(nma)$random

outres = matrix(paste0(format(round(exp(sumnma$TE),2),nsmall=2)," (",format(round(exp(sumnma$lower),2), nsmall=2)," - ",format(round(exp(sumnma$upper),2), nsmall=2),") p",ifelse(sumnma$p<0.0001,"<0.0001",ifelse(sumnma$p>0.9999,">0.9999",paste0("=",format(round(sumnma$p,4), nsmall=4))))), nrow=nrow(sumnma$TE), ncol=ncol(sumnma$TE))

colnames(outres)=colnames(sumnma$TE)
row.names(outres)=row.names(sumnma$TE)
diag(outres)=""
outres = t(outres)


netgraph(
  nma,
  number = FALSE,                    # Add node numbers for easier reference
  thickness = "se.fixed",             # Use line thickness inversely proportional to standard errors
  multiarm = F,                  # Adjust for multi-arm studies
  cex = 1.2,                        # Scale text for better readability
  cex.points = 1.5,                 # Scale node sizes
  lty = c("direct" = 1, "indirect" = 6),  # Solid lines for direct, dashed for indirect evidence
  points = TRUE                     # Show treatment nodes as points
)


#pfor=list()

# for(arm in levels(nectime5d$arm)){
#   pfor[[arm]] = forest(nma, pooled="random", reference.group = arm)
# }

comp=netsplit(nma)
indirect=comp$comparison[comp$prop.random<=0.02]
# Separate the 'indirect' comparisons into rows and columns
indirect_pairs <- as.data.frame(do.call(rbind, strsplit(indirect, ":")), stringsAsFactors = FALSE)
colnames(indirect_pairs) <- c("var1", "var2")  # Rename columns for clarity
indirect_pairs = rbind(indirect_pairs, indirect_pairs %>% rename(var1=var2,var2=var1) %>% select(var1, var2))

for (i in 1:nrow(indirect_pairs)){
  outres[indirect_pairs$var1[i],indirect_pairs$var2[i]] = paste0("\U01D43",outres[indirect_pairs$var1[i],indirect_pairs$var2[i]])
}


result_matrix5d = as.data.frame(outres)
result_matrix5d$Reference = factor(row.names(outres),levels=armlevels)
result_matrix5d = result_matrix5d %>% dynamic_select("Reference",armlevels) %>% arrange(Reference)

flextable(result_matrix5d) %>% autofit()  %>% save_as_docx( path = "OUTPUT/long_drop_some_up/survival_infectivity_adjmic_NMA_all_some_drop.docx")




# A simpler summary -------------------------------------------------------

#### RRG: Relative reduction in gametocytes####

#Similar to how we look at transmission reducing activity where we consider the decline in oocysts. Now we look at the day 2 and day 7 decline in gametocytes (both microscopy and pcr). 

##### Microscopy #####

necsub = nec %>% 
  filter((days==0 | days==2 | days==7| days==14) & !is.na(gam.m)) %>%
  select(studyid_str, study,days,arm,gam.m) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))


necsub$loggam = log(necsub$gam.m+0.1)

medgams = left_join(rbind(
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 2") %>% summarise(med=exp(mean(loggam,na.rm=T))),
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 7") %>% summarise(med=exp(mean(loggam,na.rm=T))),
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 14") %>% summarise(med=exp(mean(loggam,na.rm=T)))), necsub %>% group_by(arm) %>% filter(visit=="day 0") %>% summarise(med0=exp(mean(loggam,na.rm=T))), by="arm") %>% mutate(RR=100-(med/med0)*100)

mod_main = bam(loggam ~ arm*visit +s(studyid_str, bs="re"), data = necsub, discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

results.all = results.all %>% mutate(arm = factor(arm, levels=armlevs), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(arm, visit)

results.all$arm2 = recode(results.all$arm,
                 "DHA-PPQ + PQ (0.25mg/kg)" = "ACT-PQ",
                 "PY-AS + PQ (0.25 mg/kg)" = "ACT-PQ",
                 "DHA-PPQ + TQ (0.83 mg/kg)" = "ACT-TQ",
                 "DHA-PPQ + TQ (1.66 mg/kg)" = "ACT-TQ",
                 "AL + PQ (0.25 mg/kg)" = "ACT-PQ",
                 "AL-AQ" = "AL",
                 "AL-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                 "AS-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                 "DHA-PPQ + PQ (0.5 mg/kg)" = "ACT-PQ",
                 "SP-AQ + PQ (0.25 mg/kg)" = "Non-ACT-PQ",
                 "DHA-PPQ + PQ (0.25 mg/kg)" = "ACT-PQ"
)

results.all$arm2 = factor(results.all$arm2, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

wrapped_levels <- str_wrap(levels(results.all$arm), width = 10)
results.all$arm = factor(str_wrap(results.all$arm, width=10), levels=wrapped_levels)

ggplot(data=results.all, aes(x=arm, y=RR/100, fill=arm2, group=visit))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(aes(alpha=visit),col="grey90",stat="identity", position = position_dodge(preserve = "single"),, width=0.75)+
  scale_alpha_manual(values=c(0.4,0.6,0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_discrete()+
  ylab("Relative reduction (%) in \ngametocyte densities by microscopy")+
  coord_cartesian(ylim = c(-0.5, NA))+
  scale_y_continuous(labels=scales::percent)


ggsave("OUTPUT/long_drop_some_up/Gametocyte_density_mic_NMA_all_some_drop.png",device="png", width=17, height=8, dpi=600)


##### pCR #####


necsub = nec %>% 
  filter((days==0 | days==2 | days==7| days==14) & !is.na(gam.p)) %>%
  select(studyid_str, study,days,arm,gam.p) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))


necsub$loggam = log(necsub$gam.p+0.1)

medgams = left_join(rbind(
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 2") %>% summarise(med=exp(mean(loggam,na.rm=T))),
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 7") %>% summarise(med=exp(mean(loggam,na.rm=T))),
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 14") %>% summarise(med=exp(mean(loggam,na.rm=T)))), necsub %>% group_by(arm) %>% filter(visit=="day 0") %>% summarise(med0=exp(mean(loggam,na.rm=T))), by="arm") %>% mutate(RR=100-(med/med0)*100)

mod_main = bam(loggam ~ arm*visit +s(studyid_str, bs="re"), data = necsub, discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

results.all = results.all %>% mutate(arm = factor(arm, levels=armlevs), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(arm, visit)

results.all$arm2 = recode(results.all$arm,
                          "DHA-PPQ + PQ (0.25mg/kg)" = "ACT-PQ",
                          "PY-AS + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "DHA-PPQ + TQ (0.83 mg/kg)" = "ACT-TQ",
                          "DHA-PPQ + TQ (1.66 mg/kg)" = "ACT-TQ",
                          "AL + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "AL-AQ" = "AL",
                          "AL-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "AS-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "DHA-PPQ + PQ (0.5 mg/kg)" = "ACT-PQ",
                          "SP-AQ + PQ (0.25 mg/kg)" = "Non-ACT-PQ",
                          "DHA-PPQ + PQ (0.25 mg/kg)" = "ACT-PQ"
)

results.all$arm2 = factor(results.all$arm2, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

wrapped_levels <- str_wrap(levels(results.all$arm), width = 10)
results.all$arm = factor(str_wrap(results.all$arm, width=10), levels=wrapped_levels)

ggplot(data=results.all, aes(x=arm, y=RR/100, fill=arm2, group=visit))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(aes(alpha=visit),col="grey90",stat="identity", position = position_dodge(preserve = "single"),, width=0.75)+
  scale_alpha_manual(values=c(0.4,0.6,0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_discrete()+
  ylab("Relative reduction (%) in \ngametocyte densities by PCR")+
  coord_cartesian(ylim = c(-0.5, NA))+
  scale_y_continuous(labels=scales::percent)


ggsave("OUTPUT/long_drop_some_up/Gametocyte_density_PCR_NMA_all_some_drop.png",device="png", width=17, height=8, dpi=600)


####RRGP: Relative reduction in gametocyte prevalence####

#####RRGP by microscopy (adjusted for baseline microscopy)#####


necsub = nec %>% 
  filter((days==0 | days==2 | days==7| days==14) & !is.na(gam.m)) %>%
  select(studyid_str, study,days,arm,gam.m) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))


necsub$loggam = I(necsub$gam.m>0)


mod_main = bam(loggam ~ arm*visit +s(studyid_str, bs="re"), data = necsub, family=poisson(link="log"), discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

results.all = results.all %>% mutate(arm = factor(arm, levels=armlevs), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(arm, visit)

results.all$arm2 = recode(results.all$arm,
                          "DHA-PPQ + PQ (0.25mg/kg)" = "ACT-PQ",
                          "PY-AS + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "DHA-PPQ + TQ (0.83 mg/kg)" = "ACT-TQ",
                          "DHA-PPQ + TQ (1.66 mg/kg)" = "ACT-TQ",
                          "AL + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "AL-AQ" = "AL",
                          "AL-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "AS-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "DHA-PPQ + PQ (0.5 mg/kg)" = "ACT-PQ",
                          "SP-AQ + PQ (0.25 mg/kg)" = "Non-ACT-PQ",
                          "DHA-PPQ + PQ (0.25 mg/kg)" = "ACT-PQ"
)

results.all$arm2 = factor(results.all$arm2, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

wrapped_levels <- str_wrap(levels(results.all$arm), width = 10)
results.all$arm = factor(str_wrap(results.all$arm, width=10), levels=wrapped_levels)

ggplot(data=results.all, aes(x=arm, y=RR/100, fill=arm2, group=visit))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(aes(alpha=visit),col="grey90",stat="identity", position = position_dodge(preserve = "single"),, width=0.75)+
  scale_alpha_manual(values=c(0.4,0.6,0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_discrete()+
  ylab("Relative reduction (%) in \ngametocyte prevalence by microscopy")+
  coord_cartesian(ylim = c(-0.5, NA))+
  scale_y_continuous(labels=scales::percent)


ggsave("OUTPUT/long_drop_some_up/Gametocyte_prevalence_mic_NMA_all_some_drop.png",device="png", width=17, height=8, dpi=600)


##### pCR #####


necsub = nec %>% 
  filter((days==0 | days==2 | days==7| days==14) & !is.na(gam.p)) %>%
  select(studyid_str, study,days,arm,gam.p) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))


necsub$loggam = I(necsub$gam.p>0)

mod_main = bam(loggam ~ arm*visit +s(studyid_str, bs="re"), data = necsub, family=poisson(link="log"), discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

results.all = results.all %>% mutate(arm = factor(arm, levels=armlevs), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(arm, visit)

results.all$arm2 = recode(results.all$arm,
                          "DHA-PPQ + PQ (0.25mg/kg)" = "ACT-PQ",
                          "PY-AS + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "DHA-PPQ + TQ (0.83 mg/kg)" = "ACT-TQ",
                          "DHA-PPQ + TQ (1.66 mg/kg)" = "ACT-TQ",
                          "AL + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "AL-AQ" = "AL",
                          "AL-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "AS-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "DHA-PPQ + PQ (0.5 mg/kg)" = "ACT-PQ",
                          "SP-AQ + PQ (0.25 mg/kg)" = "Non-ACT-PQ",
                          "DHA-PPQ + PQ (0.25 mg/kg)" = "ACT-PQ"
)

results.all$arm2 = factor(results.all$arm2, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

wrapped_levels <- str_wrap(levels(results.all$arm), width = 10)
results.all$arm = factor(str_wrap(results.all$arm, width=10), levels=wrapped_levels)

ggplot(data=results.all, aes(x=arm, y=RR/100, fill=arm2, group=visit))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(aes(alpha=visit),col="grey90",stat="identity", position = position_dodge(preserve = "single"),, width=0.75)+
  scale_alpha_manual(values=c(0.4,0.6,0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_discrete()+
  ylab("Relative reduction (%) in \ngametocyte prevalence by PCR")+
  coord_cartesian(ylim = c(-0.5, NA))+
  scale_y_continuous(labels=scales::percent)


ggsave("OUTPUT/long_drop_some_up/Gametocyte_prevalence_PCR_NMA_all_some_drop.png",device="png", width=17, height=8, dpi=600)


#### TBA: transmission blocking activity ####


nec2 = nec %>% 
  group_by(studyid_str) %>%
  arrange(studyid_str, studyvisit_num) %>%
  mutate(adjust = ifelse(is.na(mosq_pos_prop) & studyvisit_num==14  & lag(studyvisit_num)==7 & lag(mosq_pos_prop)==0 , T, F),
         mosq_pos_prop = ifelse(adjust, 0, mosq_pos_prop),
         mosq_total = ifelse(adjust, 30, mosq_total),
         mosq_pos = ifelse(adjust, 0, mosq_pos),
         mean_ooc_pos = ifelse(adjust, 0, mean_ooc_pos),
         mean_ooc_all = ifelse(adjust, 0, mean_ooc_all)
  ) %>%
  ungroup()


necfeed =nec2 %>% filter(!is.na(mosq_total) & !(mosq_total==0)) 
necfeed$mosq_total = round(necfeed$mosq_total)
necfeed$prop = necfeed$mosq_pos/necfeed$mosq_total


necfeed = necfeed %>% 
  left_join(
    nec %>% 
      select(studyid_str,study,days, gam.m, gam.p) %>% 
      filter(days==0) %>%
      mutate(gam0.m = gam.m, gam0.p=gam.p) %>%
      select(studyid_str,study, gam0.m, gam0.p),
    by=c("studyid_str","study")
  ) 

necfeedsub  = necfeed %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7", "day 14")), arm=as.factor(as.character(arm)))

necfeedsub$loggam = necfeedsub$prop

necfeed = necfeedsub

necfeedsub<- necfeedsub %>%
  rowwise() %>%
  do(data.frame(
    visit = .$visit,
    arm = .$arm,
    studyid_str = .$studyid_str,
    study = .$study,
    prev = c(rep(1, (.$mosq_pos)), rep(0, .$mosq_total - .$mosq_pos))
  )) %>%
  ungroup()


necfeedsub$loggam = necfeedsub$prev

necfeedsub <- necfeedsub %>%
  mutate(row_id = row_number()) %>%
  group_by(study,arm, visit) %>%
  mutate(mosq_pos = sum(prev)) %>%
  ungroup() %>%
  group_by(study,arm, visit) %>%
  mutate(to_adjust = mosq_pos == 0) %>%
  ungroup()

# Identify the rows to be modified
rows_to_uresults.alle <- necfeedsub %>%
  filter(mosq_pos == 0) %>%
  group_by(study,arm, visit) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  mutate(new_prev = 0.5) %>%
  select(row_id, new_prev)

necfeedsub <- necfeedsub %>%
  left_join(rows_to_uresults.alle, by = "row_id") %>%
  mutate(prev = ifelse(!is.na(new_prev), new_prev, prev)) %>%
  select(-row_id, -new_prev)




mod_main = bam(prev ~ arm*visit+s(studyid_str, bs="re"), family=poisson(link="log"), data = necfeedsub, discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

results.all = results.all %>% mutate(arm = factor(arm, levels=armlevs), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(arm, visit)

results.all = results.all %>%
  mutate(RR = ifelse(visit=="day 14" & arm%in%c("SP-AQ + PQ (0.25 mg/kg)","DHA-PPQ + PQ (0.25 mg/kg)","PY-AS + PQ (0.25 mg/kg)"),NA,RR),
         ci_upper = ifelse(visit=="day 14" & arm%in%c("SP-AQ + PQ (0.25 mg/kg)","DHA-PPQ + PQ (0.25 mg/kg)","PY-AS + PQ (0.25 mg/kg)"),NA,ci_upper),
         ci_lower = ifelse(visit=="day 14" & arm%in%c("SP-AQ + PQ (0.25 mg/kg)","DHA-PPQ + PQ (0.25 mg/kg)","PY-AS + PQ (0.25 mg/kg)"),NA,ci_lower))

results.all$arm2 = recode(results.all$arm,
                          "DHA-PPQ + PQ (0.25mg/kg)" = "ACT-PQ",
                          "PY-AS + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "DHA-PPQ + TQ (0.83 mg/kg)" = "ACT-TQ",
                          "DHA-PPQ + TQ (1.66 mg/kg)" = "ACT-TQ",
                          "AL + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "AL-AQ" = "AL",
                          "AL-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "AS-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                          "DHA-PPQ + PQ (0.5 mg/kg)" = "ACT-PQ",
                          "SP-AQ + PQ (0.25 mg/kg)" = "Non-ACT-PQ",
                          "DHA-PPQ + PQ (0.25 mg/kg)" = "ACT-PQ"
)

results.all$arm2 = factor(results.all$arm2, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

wrapped_levels <- str_wrap(levels(results.all$arm), width = 10)
results.all$arm = factor(str_wrap(results.all$arm, width=10), levels=wrapped_levels)


ggplot(data=results.all, aes(x=arm, y=RR/100, fill=arm2, group=visit))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(aes(alpha=visit),col="grey90",stat="identity", position = position_dodge(preserve = "single"),, width=0.75)+
  scale_alpha_manual(values=c(0.4,0.6,0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_discrete()+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \nproportion infected mosquitos")+
  coord_cartesian(ylim = c(-0.5, NA))+
  geom_text(data=data.frame(arm="SP-AQ +\nPQ (0.25\nmg/kg)" ,arm2="Non-ACT-PQ", visit="Day 14", TRA=0),
            aes(x=arm, y=0.5, label="No Data", group=visit), # Use the same group aesthetic
            position=position_dodge(preserve="single", width=-0.5), # Align with bar dodge
            angle=90, size=5, , col=hue_pal()(length(levels(results.all$arm2)))[which(levels(results.all$arm2)=="Non-ACT-PQ")], inherit.aes=FALSE)+
  geom_text(data=data.frame(arm="DHA-PPQ +\nPQ (0.25\nmg/kg)" ,arm2="ACT-PQ", visit="Day 14", TRA=0),
            aes(x=arm, y=0.5, label="No Data", group=visit), # Use the same group aesthetic
            position=position_dodge(preserve="single", width=-0.5), # Align with bar dodge
            angle=90, size=5, , col=hue_pal()(length(levels(results.all$arm2)))[which(levels(results.all$arm2)=="ACT-PQ")], inherit.aes=FALSE)+
  geom_text(data=data.frame(arm="PY-AS +\nPQ (0.25\nmg/kg)" ,arm2="ACT-PQ", visit="Day 14", TRA=0),
            aes(x=arm, y=0.5, label="No Data", group=visit), # Use the same group aesthetic
            position=position_dodge(preserve="single", width=-0.5), # Align with bar dodge
            angle=90, size=5, , col=hue_pal()(length(levels(results.all$arm2)))[which(levels(results.all$arm2)=="ACT-PQ")], inherit.aes=FALSE)


ggsave("OUTPUT/long_drop_some_up/TBA_NMA_all_some_drop.png",device="png", width=17, height=8, dpi=600)




# Oocyst density ----------------------------------------------------------

#####unadjusted ######
necfeed = nec[!is.na(nec$mosq_total) & !(nec$mosq_total==0),]
necfeed$mosq_total = round(necfeed$mosq_total)
necfeed$prop = necfeed$mosq_pos/necfeed$mosq_total


necfeed = necfeed %>% 
  left_join(
    nec %>% 
      select(studyid_str,study,days, gam.m, gam.p) %>% 
      filter(days==0) %>%
      mutate(gam0.m = gam.m, gam0.p=gam.p) %>%
      select(studyid_str,study, gam0.m, gam0.p),
    by=c("studyid_str","study")
  ) 

necfeedsub  = necfeed %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop) & !is.na(mean_ooc_all) & !(study %in% c("PQ01","PQ03"))) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7", "day 14")), arm=as.factor(as.character(arm)))

necfeedsub$loggam = necfeedsub$mean_ooc_all


necfeedsub <- necfeedsub %>%
  mutate(row_id = row_number()) %>%
  group_by(study,arm, visit) %>%
  mutate(sumooc = sum(mean_ooc_all)) %>%
  ungroup() %>%
  group_by(study,arm, visit) %>%
  mutate(to_adjust = sumooc == 0) %>%
  ungroup()

rows_to_update <- necfeedsub %>%
  filter(sumooc == 0) %>%
  group_by(study,arm, visit) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  mutate(new_mean_ooc_all = 0.01) %>%
  select(row_id, new_mean_ooc_all)

necfeedsub <- necfeedsub %>%
  left_join(rows_to_update, by = "row_id") %>%
  mutate(mean_ooc_all = ifelse(!is.na(new_mean_ooc_all), new_mean_ooc_all, mean_ooc_all)) %>%
  select(-row_id, -new_mean_ooc_all)


mod_main = bam(mean_ooc_all ~ arm*visit+s(studyid_str, bs="re"), family=poisson(link="log"), weights = mosq_total, data = necfeedsub, discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

results.all = results.all %>% mutate(arm = factor(arm, levels=armlevs), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(arm, visit)


missdat= data.frame(arm=c("SP-AQ + PQ (0.25 mg/kg)","SP-AQ + PQ (0.25 mg/kg)","SP-AQ + PQ (0.25 mg/kg)","DHA-PPQ + PQ (0.25 mg/kg)","DHA-PPQ + PQ (0.5 mg/kg)","DHA-PPQ + PQ (0.5 mg/kg)","DHA-PPQ + PQ (0.5 mg/kg)", "PY-AS + PQ (0.25 mg/kg)"), visit=c("day 2","day 7","day 14","Day 14","day 2","day 7","day 14","day 14"), RR=NA, se=NA, ci_lower=NA, ci_upper=NA, p_value=NA) %>% mutate(arm = factor(arm, levels=armlevs), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) 

results.all = rbind(results.all, missdat)

results.all$arm2 = recode(results.all$arm,
                   "DHA-PPQ + PQ (0.25mg/kg)" = "ACT-PQ",
                   "PY-AS + PQ (0.25 mg/kg)" = "ACT-PQ",
                   "DHA-PPQ + TQ (0.83 mg/kg)" = "ACT-TQ",
                   "DHA-PPQ + TQ (1.66 mg/kg)" = "ACT-TQ",
                   "AL + PQ (0.25 mg/kg)" = "ACT-PQ",
                   "AL-AQ" = "AL",
                   "AL-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                   "AS-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                   "DHA-PPQ + PQ (0.5 mg/kg)" = "ACT-PQ",
                   "SP-AQ + PQ (0.25 mg/kg)" = "Non-ACT-PQ",
                   "DHA-PPQ + PQ (0.25 mg/kg)" = "ACT-PQ"
)

results.all$arm2 = factor(results.all$arm2, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

wrapped_levels <- str_wrap(levels(results.all$arm), width = 10)
results.all$arm = factor(str_wrap(results.all$arm, width=10), levels=wrapped_levels)


ggplot(data=results.all, aes(x=arm, y=RR/100, fill=arm2, group=visit))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(aes(alpha=visit),col="grey90",stat="identity", position = position_dodge(preserve = "single"),, width=0.75)+
  scale_alpha_manual(values=c(0.4,0.6,0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_discrete(drop = FALSE)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in oocysts\nper dissected mosquito")+
  coord_cartesian(ylim = c(-0.5, NA))+
  geom_text(data=data.frame(arm="SP-AQ +\nPQ (0.25\nmg/kg)" ,arm2="Non-ACT-PQ", TRA=0),
            aes(x=arm, y=0.5, label="No Data"), # Use the same group aesthetic
            position=position_dodge(preserve="single", width=0), # Align with bar dodge
            angle=90, size=5, col=hue_pal()(length(levels(results.all$arm2)))[which(levels(results.all$arm2)=="Non-ACT-PQ")], inherit.aes=FALSE)+
  geom_text(data=data.frame(arm="DHA-PPQ\n+ PQ (0.5\nmg/kg)" ,arm2="ACT-PQ", TRA=0),
            aes(x=arm, y=0.5, label="No Data"), # Use the same group aesthetic
            position=position_dodge(preserve="single", width=0), # Align with bar dodge
            angle=90, size=5, col=hue_pal()(length(levels(results.all$arm2)))[which(levels(results.all$arm2)=="ACT-PQ")], inherit.aes=FALSE)+
  geom_text(data=data.frame(arm="DHA-PPQ +\nPQ (0.25\nmg/kg)" ,arm2="ACT-PQ", visit="Day 14", TRA=0),
            aes(x=arm, y=0.5, label="No Data", group=visit), # Use the same group aesthetic
            position=position_dodge(preserve="single", width=-0.5), # Align with bar dodge
            angle=90, size=5,  col=hue_pal()(length(levels(results.all$arm2)))[which(levels(results.all$arm2)=="ACT-PQ")], inherit.aes=FALSE)+
  geom_text(data=data.frame(arm="PY-AS +\nPQ (0.25\nmg/kg)" ,arm2="ACT-PQ", visit="Day 14", TRA=0),
            aes(x=arm, y=0.5, label="No Data", group=visit), # Use the same group aesthetic
            position=position_dodge(preserve="single", width=-0.5), # Align with bar dodge
            angle=90, size=5, col=hue_pal()(length(levels(results.all$arm2)))[which(levels(results.all$arm2)=="ACT-PQ")], inherit.aes=FALSE)



ggsave("OUTPUT/long_drop_some_up/TRA_NMA_all_some_drop.png",device="png", width=17, height=8, dpi=600)



# Infectiousness ----------------------------------------------------------


necfeed = nec[!is.na(nec$mosq_total) & !(nec$mosq_total==0),]
necfeed$mosq_total = round(necfeed$mosq_total)
necfeed$prop = necfeed$mosq_pos/necfeed$mosq_total


necfeed = necfeed %>% 
  left_join(
    nec %>% 
      select(studyid_str,study,days, gam.m, gam.p) %>% 
      filter(days==0) %>%
      mutate(gam0.m = gam.m, gam0.p=gam.p) %>%
      select(studyid_str,study, gam0.m, gam0.p),
    by=c("studyid_str","study")
  ) 

necfeedsub  = necfeed %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop) & !is.na(mosq_pos)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7", "day 14")), arm=as.factor(as.character(arm)), loggam = as.numeric(mosq_pos>0), prev=loggam)


necfeedsub <- necfeedsub %>%
  mutate(row_id = row_number()) %>%
  group_by(study,arm, visit) %>%
  mutate(mosq_pos = sum(prev)) %>%
  ungroup() %>%
  group_by(study,arm, visit) %>%
  mutate(to_adjust = mosq_pos == 0) %>%
  ungroup()

# Identify the rows to be modified
rows_to_update <- necfeedsub %>%
  filter(mosq_pos == 0) %>%
  group_by(study,arm, visit) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  mutate(new_prev = 0.3) %>%
  select(row_id, new_prev)

necfeedsub <- necfeedsub %>%
  left_join(rows_to_update, by = "row_id") %>%
  mutate(prev = ifelse(!is.na(new_prev), new_prev, prev)) %>%
  select(-row_id, -new_prev)


mod_main = bam(prev ~ arm*visit, family=poisson(link="log"), data = necfeedsub, discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

results.all = results.all %>% mutate(arm = factor(arm, levels=armlevs), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(arm, visit)

results.all$arm2 = recode(results.all$arm,
                   "DHA-PPQ + PQ (0.25mg/kg)" = "ACT-PQ",
                   "PY-AS + PQ (0.25 mg/kg)" = "ACT-PQ",
                   "DHA-PPQ + TQ (0.83 mg/kg)" = "ACT-TQ",
                   "DHA-PPQ + TQ (1.66 mg/kg)" = "ACT-TQ",
                   "AL + PQ (0.25 mg/kg)" = "ACT-PQ",
                   "AL-AQ" = "AL",
                   "AL-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                   "AS-AQ + PQ (0.25 mg/kg)" = "ACT-PQ",
                   "DHA-PPQ + PQ (0.5 mg/kg)" = "ACT-PQ",
                   "SP-AQ + PQ (0.25 mg/kg)" = "Non-ACT-PQ",
                   "DHA-PPQ + PQ (0.25 mg/kg)" = "ACT-PQ"
)

results.all$arm2 = factor(results.all$arm2, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))


results.all = results.all %>%
  mutate(RR = ifelse(visit=="day 14" &  arm%in%c("SP-AQ + PQ (0.25 mg/kg)","DHA-PPQ + PQ (0.25 mg/kg)","PY-AS + PQ (0.25 mg/kg)"),NA,RR),
         ci_upper = ifelse(visit=="day 14" & arm%in%c("SP-AQ + PQ (0.25 mg/kg)","DHA-PPQ + PQ (0.25 mg/kg)","PY-AS + PQ (0.25 mg/kg)"),NA,ci_upper),
         ci_lower = ifelse(visit=="day 14" &  arm%in%c("SP-AQ + PQ (0.25 mg/kg)","DHA-PPQ + PQ (0.25 mg/kg)","PY-AS + PQ (0.25 mg/kg)"),NA,ci_lower))


wrapped_levels <- str_wrap(levels(results.all$arm), width = 10)
results.all$arm = factor(str_wrap(results.all$arm, width=10), levels=wrapped_levels)

ggplot(data=results.all, aes(x=arm, y=RR/100, fill=arm2, group=visit))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(aes(alpha=visit),col="grey90",stat="identity", position = position_dodge(preserve = "single"),, width=0.75)+
  scale_alpha_manual(values=c(0.4,0.6,0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_discrete()+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \nindividual's infection likelihood")+
  coord_cartesian(ylim = c(-0.5, NA))+
  geom_text(data=data.frame(arm="SP-AQ +\nPQ (0.25\nmg/kg)" ,arm2="Non-ACT-PQ", visit="Day 14", RR=0),
            aes(x=arm, y=0.5, label="No Data", group=visit), # Use the same group aesthetic
            position=position_dodge(preserve="single", width=0), # Align with bar dodge
            angle=90, size=5, , col=hue_pal()(length(levels(results.all$arm2)))[which(levels(results.all$arm2)=="Non-ACT-PQ")], inherit.aes=FALSE)+
  geom_text(data=data.frame(arm="DHA-PPQ +\nPQ (0.25\nmg/kg)" ,arm2="ACT-PQ", visit="Day 14", RR=0),
            aes(x=arm, y=0.5, label="No Data", group=visit), # Use the same group aesthetic
            position=position_dodge(preserve="single", width=-0.5), # Align with bar dodge
            angle=90, size=5, , col=hue_pal()(length(levels(results.all$arm2)))[which(levels(results.all$arm2)=="ACT-PQ")], inherit.aes=FALSE)+
  geom_text(data=data.frame(arm="PY-AS +\nPQ (0.25\nmg/kg)" ,arm2="ACT-PQ", visit="Day 14", RR=0),
            aes(x=arm, y=0.5, label="No Data", group=visit), # Use the same group aesthetic
            position=position_dodge(preserve="single", width=-0.5), # Align with bar dodge
            angle=90, size=5, , col=hue_pal()(length(levels(results.all$arm2)))[which(levels(results.all$arm2)=="ACT-PQ")], inherit.aes=FALSE)


ggsave("OUTPUT/long_drop_some_up/infectiousness_NMA_all_some_drop.png",device="png", width=17, height=8, dpi=600)

