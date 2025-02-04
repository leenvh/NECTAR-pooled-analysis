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


nec_merge[nec_merge$arm %in% c("DHA-PPQ + PQ (0.25 mg/kg)", "PY-AS + PQ (0.25 mg/kg)") & nec_merge$study == "NECTAR1" & nec_merge$studyvisit_num==14,]$mosq_pos_prop = 0

nec_merge[nec_merge$arm %in% c("DHA-PPQ + PQ (0.25 mg/kg)", "PY-AS + PQ (0.25 mg/kg)") & nec_merge$study == "NECTAR1" & nec_merge$studyvisit_num==14,]$mosq_total = 30

nec_merge[nec_merge$arm %in% c("DHA-PPQ + PQ (0.25 mg/kg)", "PY-AS + PQ (0.25 mg/kg)") & nec_merge$study == "NECTAR1" & nec_merge$studyvisit_num==14,]$mosq_pos = 0

nec_merge[nec_merge$arm %in% c("DHA-PPQ + PQ (0.25 mg/kg)", "PY-AS + PQ (0.25 mg/kg)") & nec_merge$study == "NECTAR1" & nec_merge$studyvisit_num==14,]$mean_ooc_pos = 0

nec_merge[nec_merge$arm %in% c("DHA-PPQ + PQ (0.25 mg/kg)", "PY-AS + PQ (0.25 mg/kg)") & nec_merge$study == "NECTAR1" & nec_merge$studyvisit_num==14,]$mean_ooc_all = 0



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


pq01_merge3[pq01_merge3$arm %in% c("DHA-PPQ + PQ (0.25 mg/kg)", "PY-AS + PQ (0.25 mg/kg)") & pq01_merge3$studyvisit_num==14,]$mosq_pos_prop = 0

pq01_merge3[pq01_merge3$arm %in% c("DHA-PPQ + PQ (0.25 mg/kg)", "PY-AS + PQ (0.25 mg/kg)") & pq01_merge3$studyvisit_num==14,]$mosq_total = 30

pq01_merge3[pq01_merge3$arm %in% c("DHA-PPQ + PQ (0.25 mg/kg)", "PY-AS + PQ (0.25 mg/kg)") & pq01_merge3$studyvisit_num==14,]$mosq_pos = 0

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

nectime5a$arm = factor(nectime5a$arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))


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


results$arm = factor(results$arm,, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

library(netmeta)
nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="HR", reference="ACT-PQ")

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
result_matrix5a$Reference = factor(row.names(outres),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrix5a = result_matrix5a %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)


flextable(result_matrix5a) %>% autofit()  %>% save_as_docx( path = "OUTPUT/NMA/survival_gampcr_NMA.docx")



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

nectime5b$arm = factor(nectime5b$arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

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


results$arm = factor(results$arm,, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

library(netmeta)
nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="HR", reference="ACT-PQ")

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
result_matrix5b$Reference = factor(row.names(outres),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrix5b = result_matrix5b %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)


flextable(result_matrix5b) %>% autofit()  %>% save_as_docx( path = "OUTPUT/NMA/survival_infectivity_adjpcr_NMA.docx")



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

nectime5c$arm = factor(nectime5c$arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))



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

results$arm = factor(results$arm,, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

library(netmeta)
nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="HR", reference="ACT-PQ")

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
  multiarm = TRUE, col.multiarm=colorspace::rainbow_hcl(n=4)
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
result_matrix5c$Reference = factor(row.names(outres),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrix5c = result_matrix5c %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrix5c) %>% autofit()  %>% save_as_docx( path = "OUTPUT/NMA/survival_gammic_NMA.docx")


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

nectime5d$arm = factor(nectime5d$arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))


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

results$arm = factor(results$arm,, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="HR", reference="ACT-PQ")

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
result_matrix5d$Reference = factor(row.names(outres),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrix5d = result_matrix5d %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrix5d) %>% autofit()  %>% save_as_docx( path = "OUTPUT/NMA/survival_infectivity_adjmic_NMA.docx")



# A simpler summary -------------------------------------------------------

#### RRG: Relative reduction in gametocytes####

#Similar to how we look at transmission reducing activity where we consider the decline in oocysts. Now we look at the day 2 and day 7 decline in gametocytes (both microscopy and pcr). 

##### Microscopy #####

necsub = nec %>% 
  filter((days==0 | days==2 | days==7| days==14) & !is.na(gam.m)) %>%
  select(studyid_str, study,days,arm,gam.m) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))


studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = necsub[necsub$study %in% currentstudy,]
  
  if(nrow(subdata)>0 & length(unique(subdata$visit))>1){
    
    subdata$arm = as.factor(as.character(subdata$arm))
    mod_main = gam(I(gam.m+0.1) ~ visit*arm + s(studyid_str, bs="re"),
                   data = subdata,
                   family = Gamma(link="log"))
    
    vcov_mat <- vcov(mod_main)
    vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
    vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
    mod_main$Vp = vcov_mat 
    
    newd=estimate_contrasts(mod_main, contrast=c("visit","arm"))
    newd$study = currentstudy
    newd = newd %>% select(study,Level1,Level2,Difference,SE)
    colnames(newd) = c("study","arm","ref","coef","se")
    
    results = results %>% rbind(newd) %>% na.omit()
    
  }
  
}

# Remove outer brackets from `arm` and `ref`
results$arm <- gsub("^\\(|\\)$", "", results$arm)  # Remove leading and trailing parentheses
results$ref <- gsub("^\\(|\\)$", "", results$ref)  # Remove leading and trailing parentheses

library(netmeta)
nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="RR", reference="day 0 ACT-PQ")

sumnma=summary(nma)$random

arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

#day 2 reductions
sumnmaday2 = list()
sumnmaday2$TE = diag(sumnma$TE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$seTE = diag(sumnma$seTE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$lower = sumnmaday2$TE-1.96*sumnmaday2$seTE
sumnmaday2$upper = sumnmaday2$TE+1.96*sumnmaday2$seTE
outresday2 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday2)=arms
row.names(outresday2)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday2$TE)) {
  for (j in seq_along(sumnmaday2$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday2$TE[i] - sumnmaday2$TE[j]) / sqrt(sumnmaday2$seTE[i]^2 + sumnmaday2$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday2[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday2) = paste0(format(round((1-exp(sumnmaday2$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday2$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday2$lower))*100,2), nsmall=2),"%)")


result_matrixday2 = as.data.frame(outresday2)
result_matrixday2$Reference = factor(row.names(outresday2),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday2 = result_matrixday2 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday2) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/Gam_relative_reduction_mic_day2_NMA.docx")

#day 7 reductions
sumnmaday7 = list()
sumnmaday7$TE = diag(sumnma$TE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$seTE = diag(sumnma$seTE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$lower = sumnmaday7$TE-1.96*sumnmaday7$seTE
sumnmaday7$upper = sumnmaday7$TE+1.96*sumnmaday7$seTE
outresday7 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday7)=arms
row.names(outresday7)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday7$TE)) {
  for (j in seq_along(sumnmaday7$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday7$TE[i] - sumnmaday7$TE[j]) / sqrt(sumnmaday7$seTE[i]^2 + sumnmaday7$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday7[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday7) = paste0(format(round((1-exp(sumnmaday7$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday7$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday7$lower))*100,2), nsmall=2),"%)")


result_matrixday7 = as.data.frame(outresday7)
result_matrixday7$Reference = factor(row.names(outresday7),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday7 = result_matrixday7 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday7) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/Gam_relative_reduction_mic_day7_NMA.docx")

#day 14 reductions
sumnmaday14 = list()
sumnmaday14$TE = diag(sumnma$TE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$seTE = diag(sumnma$seTE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$lower = sumnmaday14$TE-1.96*sumnmaday14$seTE
sumnmaday14$upper = sumnmaday14$TE+1.96*sumnmaday14$seTE
outresday14 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday14)=arms
row.names(outresday14)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday14$TE)) {
  for (j in seq_along(sumnmaday14$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday14$TE[i] - sumnmaday14$TE[j]) / sqrt(sumnmaday14$seTE[i]^2 + sumnmaday14$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday14[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday14) = paste0(format(round((1-exp(sumnmaday14$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday14$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday14$lower))*100,2), nsmall=2),"%)")


result_matrixday14 = as.data.frame(outresday14)
result_matrixday14$Reference = factor(row.names(outresday14),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday14 = result_matrixday14 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday14) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/Gam_relative_reduction_mic_day14_NMA.docx")


pdat= rbind(expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 2") %>% mutate(TRA=1-exp(sumnmaday2$TE), ci_lower=1-exp(sumnmaday2$TE+sumnmaday2$seTE), ci_upper=1-exp(sumnmaday2$TE-sumnmaday2$seTE)),expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 7") %>% mutate(TRA=1-exp(sumnmaday7$TE), ci_lower=1-exp(sumnmaday7$TE+sumnmaday7$seTE), ci_upper=1-exp(sumnmaday7$TE-sumnmaday7$seTE)), expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 14") %>% mutate(TRA=1-exp(sumnmaday14$TE), ci_lower=1-exp(sumnmaday14$TE+sumnmaday14$seTE), ci_upper=1-exp(sumnmaday14$TE-sumnmaday14$seTE)))  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 


ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \ngametocyte densities by microscopy")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA/Gametocyte_density_mic_NMA.png",device="png", width=10, height=7, dpi=600)


##### pCR #####

necsub = nec %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(gam.p)) %>%
  select(studyid_str,study,days,arm,gam.p) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))


studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = necsub[necsub$study %in% currentstudy,]
  
  if(nrow(subdata)>0 & length(unique(subdata$visit))>1){
    
    subdata$arm = as.factor(as.character(subdata$arm))
    mod_main = gam(I(gam.p+0.01) ~ visit*arm + s(studyid_str, bs="re"),
                   data = subdata,
                   family = Gamma(link="log"))
    vcov_mat <- vcov(mod_main)
    vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
    vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
    mod_main$Vp = vcov_mat 
    
    newd=estimate_contrasts(mod_main, contrast=c("visit","arm"))
    newd$study = currentstudy
    newd = newd %>% select(study,Level1,Level2,Difference,SE)
    colnames(newd) = c("study","arm","ref","coef","se")
    
    results = results %>% rbind(newd) %>% na.omit()
    
  }
  
}

# Remove outer brackets from `arm` and `ref`
results$arm <- gsub("^\\(|\\)$", "", results$arm)  # Remove leading and trailing parentheses
results$ref <- gsub("^\\(|\\)$", "", results$ref)  # Remove leading and trailing parentheses

library(netmeta)
nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="RR", reference="day 0 ACT-PQ")

sumnma=summary(nma)$random

arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

#day 2 reductions
sumnmaday2 = list()
sumnmaday2$TE = diag(sumnma$TE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$seTE = diag(sumnma$seTE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$lower = sumnmaday2$TE-1.96*sumnmaday2$seTE
sumnmaday2$upper = sumnmaday2$TE+1.96*sumnmaday2$seTE
outresday2 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday2)=arms
row.names(outresday2)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday2$TE)) {
  for (j in seq_along(sumnmaday2$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday2$TE[i] - sumnmaday2$TE[j]) / sqrt(sumnmaday2$seTE[i]^2 + sumnmaday2$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday2[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday2) = paste0(format(round((1-exp(sumnmaday2$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday2$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday2$lower))*100,2), nsmall=2),"%)")


result_matrixday2 = as.data.frame(outresday2)
result_matrixday2$Reference = factor(row.names(outresday2),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday2 = result_matrixday2 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday2) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/Gam_relative_reduction_PCR_day2_NMA.docx")

#day 7 reductions
sumnmaday7 = list()
sumnmaday7$TE = diag(sumnma$TE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$seTE = diag(sumnma$seTE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$lower = sumnmaday7$TE-1.96*sumnmaday7$seTE
sumnmaday7$upper = sumnmaday7$TE+1.96*sumnmaday7$seTE
outresday7 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday7)=arms
row.names(outresday7)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday7$TE)) {
  for (j in seq_along(sumnmaday7$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday7$TE[i] - sumnmaday7$TE[j]) / sqrt(sumnmaday7$seTE[i]^2 + sumnmaday7$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday7[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday7) = paste0(format(round((1-exp(sumnmaday7$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday7$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday7$lower))*100,2), nsmall=2),"%)")


result_matrixday7 = as.data.frame(outresday7)
result_matrixday7$Reference = factor(row.names(outresday7),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday7 = result_matrixday7 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday7) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/Gam_relative_reduction_PCR_day7_NMA.docx")

#day 14 reductions
sumnmaday14 = list()
sumnmaday14$TE = diag(sumnma$TE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$seTE = diag(sumnma$seTE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$lower = sumnmaday14$TE-1.96*sumnmaday14$seTE
sumnmaday14$upper = sumnmaday14$TE+1.96*sumnmaday14$seTE
outresday14 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday14)=arms
row.names(outresday14)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday14$TE)) {
  for (j in seq_along(sumnmaday14$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday14$TE[i] - sumnmaday14$TE[j]) / sqrt(sumnmaday14$seTE[i]^2 + sumnmaday14$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday14[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday14) = paste0(format(round((1-exp(sumnmaday14$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday14$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday14$lower))*100,2), nsmall=2),"%)")


result_matrixday14 = as.data.frame(outresday14)
result_matrixday14$Reference = factor(row.names(outresday14),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday14 = result_matrixday14 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday14) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/Gam_relative_reduction_PCR_day14_NMA.docx")


pdat= rbind(expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 2") %>% mutate(TRA=1-exp(sumnmaday2$TE), ci_lower=1-exp(sumnmaday2$TE+sumnmaday2$seTE), ci_upper=1-exp(sumnmaday2$TE-sumnmaday2$seTE)),expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 7") %>% mutate(TRA=1-exp(sumnmaday7$TE), ci_lower=1-exp(sumnmaday7$TE+sumnmaday7$seTE), ci_upper=1-exp(sumnmaday7$TE-sumnmaday7$seTE)), expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 14") %>% mutate(TRA=1-exp(sumnmaday14$TE), ci_lower=1-exp(sumnmaday14$TE+sumnmaday14$seTE), ci_upper=1-exp(sumnmaday14$TE-sumnmaday14$seTE)))  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 


ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \ngametocyte densities by PCR")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA/Gametocyte_density_PCR_NMA.png",device="png", width=10, height=7, dpi=600)


####RRGP: Relative reduction in gametocyte prevalence####

#####RRGP by microscopy (adjusted for baseline microscopy)#####


necsub = nec %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(gam.m)) %>%
  select(studyid_str,study,days,arm,gam.m) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm))) %>%
  group_by(studyid_str) %>%
  mutate(gam0.m=ifelse(first(days)!=0,NA,first(gam.m))) %>%
  ungroup()

studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = necsub[necsub$study %in% currentstudy,]
  
  if(nrow(subdata)>0 & length(unique(subdata$visit))>1){
    
    subdata$arm = as.factor(as.character(subdata$arm))
    subdata$log = log10(subdata$gam0.m+0.1)
    mod_main = bam(I(gam.m>0) ~ visit:arm-1 + log + s(studyid_str, bs="re"),
                   data = subdata,
                   family = poisson(link="log"))

    vcov_mat <- vcov(mod_main)
    vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
    vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
    mod_main$Vp = vcov_mat 
    
    newd=estimate_contrasts(mod_main, contrast=c("visit","arm"))
    newd$study = currentstudy
    newd = newd %>% select(study,Level1,Level2,Difference,SE)
    colnames(newd) = c("study","arm","ref","coef","se")
    
    results = results %>% rbind(newd) %>% na.omit()
    
  }
  
}

# Remove outer brackets from `arm` and `ref`
results$arm <- gsub("^\\(|\\)$", "", results$arm)  # Remove leading and trailing parentheses
results$ref <- gsub("^\\(|\\)$", "", results$ref)  # Remove leading and trailing parentheses


nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="RR", reference="day 0 AL")

sumnma=summary(nma)$random

arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

#day 2 reductions
sumnmaday2 = list()
sumnmaday2$TE = diag(sumnma$TE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$seTE = diag(sumnma$seTE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$lower = sumnmaday2$TE-1.96*sumnmaday2$seTE
sumnmaday2$upper = sumnmaday2$TE+1.96*sumnmaday2$seTE
outresday2 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday2)=arms
row.names(outresday2)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday2$TE)) {
  for (j in seq_along(sumnmaday2$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday2$TE[i] - sumnmaday2$TE[j]) / sqrt(sumnmaday2$seTE[i]^2 + sumnmaday2$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday2[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday2) = paste0(format(round((1-exp(sumnmaday2$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday2$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday2$lower))*100,2), nsmall=2),"%)")


result_matrixday2 = as.data.frame(outresday2)
result_matrixday2$Reference = factor(row.names(outresday2),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday2 = result_matrixday2 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday2) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/Gam_relative_reduction_prev_mic_day2_NMA.docx")

#day 7 reductions
sumnmaday7 = list()
sumnmaday7$TE = diag(sumnma$TE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$seTE = diag(sumnma$seTE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$lower = sumnmaday7$TE-1.96*sumnmaday7$seTE
sumnmaday7$upper = sumnmaday7$TE+1.96*sumnmaday7$seTE
outresday7 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday7)=arms
row.names(outresday7)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday7$TE)) {
  for (j in seq_along(sumnmaday7$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday7$TE[i] - sumnmaday7$TE[j]) / sqrt(sumnmaday7$seTE[i]^2 + sumnmaday7$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday7[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday7) = paste0(format(round((1-exp(sumnmaday7$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday7$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday7$lower))*100,2), nsmall=2),"%)")


result_matrixday7 = as.data.frame(outresday7)
result_matrixday7$Reference = factor(row.names(outresday7),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday7 = result_matrixday7 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday7) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/Gam_relative_reduction_prev_mic_day7_NMA.docx")

#day 14 reductions
sumnmaday14 = list()
sumnmaday14$TE = diag(sumnma$TE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$seTE = diag(sumnma$seTE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$lower = sumnmaday14$TE-1.96*sumnmaday14$seTE
sumnmaday14$upper = sumnmaday14$TE+1.96*sumnmaday14$seTE
outresday14 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday14)=arms
row.names(outresday14)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday14$TE)) {
  for (j in seq_along(sumnmaday14$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday14$TE[i] - sumnmaday14$TE[j]) / sqrt(sumnmaday14$seTE[i]^2 + sumnmaday14$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday14[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday14) = paste0(format(round((1-exp(sumnmaday14$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday14$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday14$lower))*100,2), nsmall=2),"%)")


result_matrixday14 = as.data.frame(outresday14)
result_matrixday14$Reference = factor(row.names(outresday14),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday14 = result_matrixday14 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday14) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/Gam_relative_reduction_prev_mic_day14_NMA.docx")


pdat= rbind(expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 2") %>% mutate(TRA=1-exp(sumnmaday2$TE), ci_lower=1-exp(sumnmaday2$TE+sumnmaday2$seTE), ci_upper=1-exp(sumnmaday2$TE-sumnmaday2$seTE)),expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 7") %>% mutate(TRA=1-exp(sumnmaday7$TE), ci_lower=1-exp(sumnmaday7$TE+sumnmaday7$seTE), ci_upper=1-exp(sumnmaday7$TE-sumnmaday7$seTE)), expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 14") %>% mutate(TRA=1-exp(sumnmaday14$TE), ci_lower=1-exp(sumnmaday14$TE+sumnmaday14$seTE), ci_upper=1-exp(sumnmaday14$TE-sumnmaday14$seTE)))  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 


ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \ngametocyte prevalence by microscopy")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA/Gametocyte_density_prev_mic_NMA.png",device="png", width=10, height=7, dpi=600)



#####RRGP by PCR (adjusted for baseline PCR)#####

necsub = nec %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(gam.p)) %>%
  select(studyid_str,study,days,arm,gam.p) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))%>%
  group_by(studyid_str) %>%
  mutate(gam0.p=ifelse(first(days)!=0,NA,first(gam.p))) %>%
  ungroup()




studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = necsub[necsub$study %in% currentstudy,]
  
  if(nrow(subdata)>0 & length(unique(subdata$visit))>1){
    
    subdata$arm = as.factor(as.character(subdata$arm))
    subdata$log = log10(subdata$gam0.p+0.01)
    mod_main = bam(I(gam.p>0) ~ visit:arm-1 + log + s(studyid_str, bs="re"),
                   data = subdata,
                   family = poisson(link="log"))
    
    vcov_mat <- vcov(mod_main)
    vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
    vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
    mod_main$Vp = vcov_mat 
    
    newd=estimate_contrasts(mod_main, contrast=c("visit","arm"))
    newd$study = currentstudy
    newd = newd %>% select(study,Level1,Level2,Difference,SE)
    colnames(newd) = c("study","arm","ref","coef","se")
    
    results = results %>% rbind(newd) %>% na.omit()
    
  }
  
}

# Remove outer brackets from `arm` and `ref`
results$arm <- gsub("^\\(|\\)$", "", results$arm)  # Remove leading and trailing parentheses
results$ref <- gsub("^\\(|\\)$", "", results$ref)  # Remove leading and trailing parentheses


nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="RR", reference="day 0 AL")

sumnma=summary(nma)$random

arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

#day 2 reductions
sumnmaday2 = list()
sumnmaday2$TE = diag(sumnma$TE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$seTE = diag(sumnma$seTE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$lower = sumnmaday2$TE-1.96*sumnmaday2$seTE
sumnmaday2$upper = sumnmaday2$TE+1.96*sumnmaday2$seTE
outresday2 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday2)=arms
row.names(outresday2)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday2$TE)) {
  for (j in seq_along(sumnmaday2$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday2$TE[i] - sumnmaday2$TE[j]) / sqrt(sumnmaday2$seTE[i]^2 + sumnmaday2$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday2[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday2) = paste0(format(round((1-exp(sumnmaday2$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday2$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday2$lower))*100,2), nsmall=2),"%)")


result_matrixday2 = as.data.frame(outresday2)
result_matrixday2$Reference = factor(row.names(outresday2),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday2 = result_matrixday2 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday2) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/Gam_relative_reduction_prev_PCR_day2_NMA.docx")

#day 7 reductions
sumnmaday7 = list()
sumnmaday7$TE = diag(sumnma$TE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$seTE = diag(sumnma$seTE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$lower = sumnmaday7$TE-1.96*sumnmaday7$seTE
sumnmaday7$upper = sumnmaday7$TE+1.96*sumnmaday7$seTE
outresday7 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday7)=arms
row.names(outresday7)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday7$TE)) {
  for (j in seq_along(sumnmaday7$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday7$TE[i] - sumnmaday7$TE[j]) / sqrt(sumnmaday7$seTE[i]^2 + sumnmaday7$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday7[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday7) = paste0(format(round((1-exp(sumnmaday7$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday7$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday7$lower))*100,2), nsmall=2),"%)")


result_matrixday7 = as.data.frame(outresday7)
result_matrixday7$Reference = factor(row.names(outresday7),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday7 = result_matrixday7 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday7) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/Gam_relative_reduction_prev_PCR_day7_NMA.docx")

#day 14 reductions
sumnmaday14 = list()
sumnmaday14$TE = diag(sumnma$TE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$seTE = diag(sumnma$seTE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$lower = sumnmaday14$TE-1.96*sumnmaday14$seTE
sumnmaday14$upper = sumnmaday14$TE+1.96*sumnmaday14$seTE
outresday14 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday14)=arms
row.names(outresday14)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday14$TE)) {
  for (j in seq_along(sumnmaday14$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday14$TE[i] - sumnmaday14$TE[j]) / sqrt(sumnmaday14$seTE[i]^2 + sumnmaday14$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday14[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday14) = paste0(format(round((1-exp(sumnmaday14$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday14$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday14$lower))*100,2), nsmall=2),"%)")


result_matrixday14 = as.data.frame(outresday14)
result_matrixday14$Reference = factor(row.names(outresday14),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday14 = result_matrixday14 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday14) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/Gam_relative_reduction_prev_PCR_day14_NMA.docx")


pdat= rbind(expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 2") %>% mutate(TRA=1-exp(sumnmaday2$TE), ci_lower=1-exp(sumnmaday2$TE+sumnmaday2$seTE), ci_upper=1-exp(sumnmaday2$TE-sumnmaday2$seTE)),expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 7") %>% mutate(TRA=1-exp(sumnmaday7$TE), ci_lower=1-exp(sumnmaday7$TE+sumnmaday7$seTE), ci_upper=1-exp(sumnmaday7$TE-sumnmaday7$seTE)), expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 14") %>% mutate(TRA=1-exp(sumnmaday14$TE), ci_lower=1-exp(sumnmaday14$TE+sumnmaday14$seTE), ci_upper=1-exp(sumnmaday14$TE-sumnmaday14$seTE)))  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 


ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \ngametocyte prevalence by PCR")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA/Gametocyte_density_prev_PCR_NMA.png",device="png", width=10, height=7, dpi=600)



#### TBA: transmission blocking activity ####

##### Unadjusted #####

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
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7", "day 14")), arm=as.factor(as.character(arm)))

studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = necfeedsub[necfeedsub$study %in% currentstudy,]
  
  subdata <- subdata %>%
    rowwise() %>%
    do(data.frame(
      visit = .$visit,
      arm = .$arm,
      studyid_str = .$studyid_str,
      prev = c(rep(1, .$mosq_pos), rep(0, .$mosq_total - .$mosq_pos))
    )) %>%
    ungroup()
  
  
  if(nrow(subdata)>0 & length(unique(subdata$visit))>1){
    

    mod_main = bam(prev ~ visit:arm-1 +s(studyid_str, bs="re"),
                      family=poisson(link="log"), 
                      data=subdata)
    
    vcov_mat <- vcov(mod_main)
    vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
    vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
    mod_main$Vp = vcov_mat 
     
    newd=estimate_contrasts(mod_main, contrast=c("visit","arm"))
    newd$study = currentstudy
    newd = newd %>% select(study,Level1,Level2,Difference,SE)
    colnames(newd) = c("study","arm","ref","coef","se")
    
    results = results %>% rbind(newd) %>% na.omit()
    
  }
  
}

# Remove outer brackets from `arm` and `ref`
results$arm <- gsub("^\\(|\\)$", "", results$arm)  # Remove leading and trailing parentheses
results$ref <- gsub("^\\(|\\)$", "", results$ref)  # Remove leading and trailing parentheses

results$se[results$se==0] = 0.001

nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="RR", reference="day 0 AL")

sumnma=summary(nma)$random

arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

#day 2 reductions
sumnmaday2 = list()
sumnmaday2$TE = diag(sumnma$TE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$seTE = diag(sumnma$seTE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$lower = sumnmaday2$TE-1.96*sumnmaday2$seTE
sumnmaday2$upper = sumnmaday2$TE+1.96*sumnmaday2$seTE
outresday2 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday2)=arms
row.names(outresday2)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday2$TE)) {
  for (j in seq_along(sumnmaday2$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday2$TE[i] - sumnmaday2$TE[j]) / sqrt(sumnmaday2$seTE[i]^2 + sumnmaday2$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday2[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday2) = paste0(format(round((1-exp(sumnmaday2$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday2$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday2$lower))*100,2), nsmall=2),"%)")


result_matrixday2 = as.data.frame(outresday2)
result_matrixday2$Reference = factor(row.names(outresday2),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday2 = result_matrixday2 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday2) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TBA_day2_unadjusted_NMA.docx")

#day 7 reductions
sumnmaday7 = list()
sumnmaday7$TE = diag(sumnma$TE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$seTE = diag(sumnma$seTE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$lower = sumnmaday7$TE-1.96*sumnmaday7$seTE
sumnmaday7$upper = sumnmaday7$TE+1.96*sumnmaday7$seTE
outresday7 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday7)=arms
row.names(outresday7)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday7$TE)) {
  for (j in seq_along(sumnmaday7$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday7$TE[i] - sumnmaday7$TE[j]) / sqrt(sumnmaday7$seTE[i]^2 + sumnmaday7$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday7[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday7) = paste0(format(round((1-exp(sumnmaday7$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday7$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday7$lower))*100,2), nsmall=2),"%)")


result_matrixday7 = as.data.frame(outresday7)
result_matrixday7$Reference = factor(row.names(outresday7),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday7 = result_matrixday7 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday7) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TBA_day7_unadjusted_NMA.docx")

#day 14 reductions

#No feed data for Non-ACT-PQ

arms=arms[!(arms %in% "Non-ACT-PQ")]

sumnmaday14 = list()
sumnmaday14$TE = diag(sumnma$TE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$seTE = diag(sumnma$seTE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$lower = sumnmaday14$TE-1.96*sumnmaday14$seTE
sumnmaday14$upper = sumnmaday14$TE+1.96*sumnmaday14$seTE
outresday14 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday14)=arms
row.names(outresday14)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday14$TE)) {
  for (j in seq_along(sumnmaday14$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday14$TE[i] - sumnmaday14$TE[j]) / sqrt(sumnmaday14$seTE[i]^2 + sumnmaday14$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday14[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday14) = paste0(format(round((1-exp(sumnmaday14$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday14$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday14$lower))*100,2), nsmall=2),"%)")


result_matrixday14 = as.data.frame(outresday14)
result_matrixday14$Reference = factor(row.names(outresday14),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday14 = result_matrixday14 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL",  "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday14) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TBA_day14_unadjusted_NMA.docx")


pdat= rbind(expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 2") %>% mutate(TRA=1-exp(sumnmaday2$TE), ci_lower=1-exp(sumnmaday2$TE+sumnmaday2$seTE), ci_upper=1-exp(sumnmaday2$TE-sumnmaday2$seTE)),expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 7") %>% mutate(TRA=1-exp(sumnmaday7$TE), ci_lower=1-exp(sumnmaday7$TE+sumnmaday7$seTE), ci_upper=1-exp(sumnmaday7$TE-sumnmaday7$seTE)), expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 14") %>% mutate(TRA=1-exp(sumnmaday14$TE), ci_lower=1-exp(sumnmaday14$TE+sumnmaday14$seTE), ci_upper=1-exp(sumnmaday14$TE-sumnmaday14$seTE)))  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 


ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \nproportion infected mosquitos")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA/TBA_unadjusted_NMA.png",device="png", width=10, height=7, dpi=600)



# 
# #Heatmap for day2
# result_matrixday2 <- as.data.frame(result_matrixday2)
# rownames(result_matrixday2) <- result_matrixday2[, 1] 
# result_matrixday2 <- result_matrixday2[, -1]
# 
# desired_order <- c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")
# 
# result_long <- result_matrixday2 %>%
#   rownames_to_column("treatment") %>%
#   pivot_longer(cols = -treatment, names_to = "comparison", values_to = "value") %>%
#   mutate(
#     p_value_numeric = ifelse(grepl("p<", value), 0.0001,  # Explicit handling of values "<0.0001"
#                              ifelse(grepl("p=", value), 
#                                     as.numeric(gsub(".*p=([0-9.]+).*", "\\1", value)),
#                                     NA_real_)),
#     label = ifelse(treatment == comparison,
#                    gsub(" \\(", "\n(", value) %>% gsub(", ", ",\n", .),  # Add newlines for diagonal
#                    ifelse(grepl("p<0.0001", value), "p<0.0001", sprintf("p=%.4f", p_value_numeric))),
#     p_value_log = ifelse(is.na(p_value_numeric), NA, -log10(p_value_numeric)),  # Avoid log10 of NA
#     border = ifelse(treatment == comparison, "black", NA)  # Black border for diagonal
#   ) %>%
#   drop_na(p_value_numeric) %>%
#   mutate(
#     treatment = factor(treatment, levels = desired_order),
#     comparison = factor(comparison, levels = desired_order)
#   ) %>%
#   filter(as.numeric(treatment) <= as.numeric(comparison))
# 
# 
# # Function to extract the numeric value from the activity string
# extract_activity <- function(value) {
#   matches <- regmatches(value, regexpr("(-?[0-9]+\\.?[0-9]*)%", value))
#   if (length(matches) > 0) {
#     return(as.numeric(sub("%", "", matches)))
#   }
#   return(NA)
# }
# 
# # Extract the numeric activity values
# result_long <- result_long %>%
#   mutate(Activity_numeric = sapply(value, extract_activity))
# 
# # Initialize Activity_diff column
# result_long$Activity_diff <- NA
# 
# # Calculate the absolute difference in activity
# for (i in 1:nrow(result_long)) {
#   treatment <- result_long$treatment[i]
#   comparison <- result_long$comparison[i]
#   
#   # Find the activity values for the treatment and comparison
#   activity_treatment <- result_long$Activity_numeric[
#     result_long$treatment == treatment & result_long$comparison == treatment]
#   activity_comparison <- result_long$Activity_numeric[
#     result_long$treatment == comparison & result_long$comparison == comparison]
#   
#   if (length(activity_treatment) > 0 && length(activity_comparison) > 0) {
#     # Calculate the absolute difference
#     result_long$Activity_diff[i] <- abs(activity_treatment - activity_comparison)
#   }
# }
# 
# result_long <- result_long %>%
#   mutate(is_diagonal = ifelse(treatment == comparison, TRUE, FALSE))
# 
# # Plotting the heatmap
# ggplot(result_long, aes(x = comparison, y = treatment, fill = Activity_diff)) +
#   geom_tile(color = "white", size = 0.1) +  # Default borders for non-diagonal
#   geom_tile(data = filter(result_long, !is.na(border)), color = "white", size = 1, fill = NA) + 
#   geom_tile(data = filter(result_long, is_diagonal), color = "black", size = 0.1, fill = "white") + 
#   scale_fill_gradient2(low="#bfbbbb",high = "#508ea1", name = "Absolute difference in relative reduction",limits=c(0,120),breaks = c(0, 20, 40, 60, 80,100,120), labels = c("0%", "20%", "40%", "60%", "80%","100%","120%")) +
#   geom_text(
#     aes(label = label),  # Use the formatted label
#     size = 3, 
#     angle = 45  # Rotating text to 45 degrees
#   ) +
#   coord_fixed() +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.text.y = element_text(angle = 45, hjust = 1),
#     axis.title = element_blank(),
#     panel.grid = element_blank(),
#     axis.ticks = element_blank(),
#     legend.position = "right"
#   ) +
#   labs(title = "P-value Heatmap", x = "Comparison Group", y = "Treatment Group")
# 
# ggsave("OUTPUT/NMA/Heatmap_TBA_Day2.png",device="png", width=12, height=8, dpi=600)
# 
# 
# #Heatmap for day7
# result_matrixday7 <- as.data.frame(result_matrixday7)
# rownames(result_matrixday7) <- result_matrixday7[, 1] 
# result_matrixday7 <- result_matrixday7[, -1]
# 
# desired_order <- c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")
# 
# result_long <- result_matrixday7 %>%
#   rownames_to_column("treatment") %>%
#   pivot_longer(cols = -treatment, names_to = "comparison", values_to = "value") %>%
#   mutate(
#     p_value_numeric = ifelse(grepl("p<", value), 0.0001,  # Explicit handling of values "<0.0001"
#                              ifelse(grepl("p=", value), 
#                                     as.numeric(gsub(".*p=([0-9.]+).*", "\\1", value)),
#                                     NA_real_)),
#     label = ifelse(treatment == comparison,
#                    gsub(" \\(", "\n(", value) %>% gsub(", ", ",\n", .),  # Add newlines for diagonal
#                    ifelse(grepl("p<0.0001", value), "p<0.0001", sprintf("p=%.4f", p_value_numeric))),
#     p_value_log = ifelse(is.na(p_value_numeric), NA, -log10(p_value_numeric)),  # Avoid log10 of NA
#     border = ifelse(treatment == comparison, "black", NA)  # Black border for diagonal
#   ) %>%
#   drop_na(p_value_numeric) %>%
#   mutate(
#     treatment = factor(treatment, levels = desired_order),
#     comparison = factor(comparison, levels = desired_order)
#   ) %>%
#   filter(as.numeric(treatment) <= as.numeric(comparison))
# 
# 
# # Function to extract the numeric value from the activity string
# extract_activity <- function(value) {
#   matches <- regmatches(value, regexpr("(-?[0-9]+\\.?[0-9]*)%", value))
#   if (length(matches) > 0) {
#     return(as.numeric(sub("%", "", matches)))
#   }
#   return(NA)
# }
# 
# # Extract the numeric activity values
# result_long <- result_long %>%
#   mutate(Activity_numeric = sapply(value, extract_activity))
# 
# # Initialize Activity_diff column
# result_long$Activity_diff <- NA
# 
# # Calculate the absolute difference in activity
# for (i in 1:nrow(result_long)) {
#   treatment <- result_long$treatment[i]
#   comparison <- result_long$comparison[i]
#   
#   # Find the activity values for the treatment and comparison
#   activity_treatment <- result_long$Activity_numeric[
#     result_long$treatment == treatment & result_long$comparison == treatment]
#   activity_comparison <- result_long$Activity_numeric[
#     result_long$treatment == comparison & result_long$comparison == comparison]
#   
#   if (length(activity_treatment) > 0 && length(activity_comparison) > 0) {
#     # Calculate the absolute difference
#     result_long$Activity_diff[i] <- abs(activity_treatment - activity_comparison)
#   }
# }
# 
# result_long <- result_long %>%
#   mutate(is_diagonal = ifelse(treatment == comparison, TRUE, FALSE))
# 
# # Plotting the heatmap
# ggplot(result_long, aes(x = comparison, y = treatment, fill = Activity_diff)) +
#   geom_tile(color = "white", size = 0.1) +  # Default borders for non-diagonal
#   geom_tile(data = filter(result_long, !is.na(border)), color = "white", size = 1, fill = NA) + 
#   geom_tile(data = filter(result_long, is_diagonal), color = "black", size = 0.1, fill = "white") + 
#   scale_fill_gradient2(low="#bfbbbb",high = "#508ea1", name = "Absolute difference in relative reduction",limits=c(0,120),breaks = c(0, 20, 40, 60, 80,100,120), labels = c("0%", "20%", "40%", "60%", "80%","100%","120%")) +
#   geom_text(
#     aes(label = label),  # Use the formatted label
#     size = 3, 
#     angle = 45  # Rotating text to 45 degrees
#   ) +
#   coord_fixed() +
#   theme_minimal() +
#   theme(
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     axis.text.y = element_text(angle = 45, hjust = 1),
#     axis.title = element_blank(),
#     panel.grid = element_blank(),
#     axis.ticks = element_blank(),
#     legend.position = "right"
#   ) +
#   labs(title = "P-value Heatmap", x = "Comparison Group", y = "Treatment Group")
# 
# ggsave("OUTPUT/NMA/Heatmap_TBA_Day7.png",device="png", width=12, height=8, dpi=600)
# 
# 




##### Adjusted by baseline gametocyte density (pcr) #####


library(modelbased)


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
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7", "day 14")), arm=as.factor(as.character(arm)))

necfeedsub$log = I(log10(necfeedsub$gam0.p+0.01))

studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = necfeedsub[necfeedsub$study %in% currentstudy,]
  
  subdata <- subdata %>%
    rowwise() %>%
    do(data.frame(
      visit = .$visit,
      arm = .$arm,
      studyid_str = .$studyid_str,
      log=.$log,
      prev = c(rep(1, .$mosq_pos), rep(0, .$mosq_total - .$mosq_pos))
    )) %>%
    ungroup()
  
  
  if(nrow(subdata)>0 & length(unique(subdata$visit))>1){
    
    
    mod_main = bam(prev ~ visit:arm-1 + log +s(studyid_str, bs="re"),
                   family=poisson(link="log"), 
                   data=subdata)
    
    vcov_mat <- vcov(mod_main)
    vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
    vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
    mod_main$Vp = vcov_mat 
    
    newd=estimate_contrasts(mod_main, contrast=c("visit","arm"))
    newd$study = currentstudy
    newd = newd %>% select(study,Level1,Level2,Difference,SE)
    colnames(newd) = c("study","arm","ref","coef","se")
    
    results = results %>% rbind(newd) %>% na.omit()
    
  }
  
}

# Remove outer brackets from `arm` and `ref`
results$arm <- gsub("^\\(|\\)$", "", results$arm)  # Remove leading and trailing parentheses
results$ref <- gsub("^\\(|\\)$", "", results$ref)  # Remove leading and trailing parentheses

results$se[results$se==0] = 0.001

nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="RR", reference="day 0 AL")

sumnma=summary(nma)$random

arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

#day 2 reductions
sumnmaday2 = list()
sumnmaday2$TE = diag(sumnma$TE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$seTE = diag(sumnma$seTE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$lower = sumnmaday2$TE-1.96*sumnmaday2$seTE
sumnmaday2$upper = sumnmaday2$TE+1.96*sumnmaday2$seTE
outresday2 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday2)=arms
row.names(outresday2)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday2$TE)) {
  for (j in seq_along(sumnmaday2$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday2$TE[i] - sumnmaday2$TE[j]) / sqrt(sumnmaday2$seTE[i]^2 + sumnmaday2$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday2[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday2) = paste0(format(round((1-exp(sumnmaday2$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday2$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday2$lower))*100,2), nsmall=2),"%)")


result_matrixday2 = as.data.frame(outresday2)
result_matrixday2$Reference = factor(row.names(outresday2),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday2 = result_matrixday2 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday2) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TBA_day2_adjusted_pcr_NMA.docx")

#day 7 reductions
sumnmaday7 = list()
sumnmaday7$TE = diag(sumnma$TE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$seTE = diag(sumnma$seTE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$lower = sumnmaday7$TE-1.96*sumnmaday7$seTE
sumnmaday7$upper = sumnmaday7$TE+1.96*sumnmaday7$seTE
outresday7 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday7)=arms
row.names(outresday7)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday7$TE)) {
  for (j in seq_along(sumnmaday7$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday7$TE[i] - sumnmaday7$TE[j]) / sqrt(sumnmaday7$seTE[i]^2 + sumnmaday7$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday7[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday7) = paste0(format(round((1-exp(sumnmaday7$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday7$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday7$lower))*100,2), nsmall=2),"%)")


result_matrixday7 = as.data.frame(outresday7)
result_matrixday7$Reference = factor(row.names(outresday7),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday7 = result_matrixday7 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday7) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TBA_day7_adjusted_pcr_NMA.docx")

#day 14 reductions

#No feed data for Non-ACT-PQ

arms=arms[!(arms %in% "Non-ACT-PQ")]

sumnmaday14 = list()
sumnmaday14$TE = diag(sumnma$TE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$seTE = diag(sumnma$seTE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$lower = sumnmaday14$TE-1.96*sumnmaday14$seTE
sumnmaday14$upper = sumnmaday14$TE+1.96*sumnmaday14$seTE
outresday14 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday14)=arms
row.names(outresday14)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday14$TE)) {
  for (j in seq_along(sumnmaday14$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday14$TE[i] - sumnmaday14$TE[j]) / sqrt(sumnmaday14$seTE[i]^2 + sumnmaday14$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday14[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday14) = paste0(format(round((1-exp(sumnmaday14$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday14$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday14$lower))*100,2), nsmall=2),"%)")


result_matrixday14 = as.data.frame(outresday14)
result_matrixday14$Reference = factor(row.names(outresday14),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday14 = result_matrixday14 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL",  "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday14) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TBA_day14_adjusted_pcr_NMA.docx")


pdat= rbind(expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 2") %>% mutate(TRA=1-exp(sumnmaday2$TE), ci_lower=1-exp(sumnmaday2$TE+sumnmaday2$seTE), ci_upper=1-exp(sumnmaday2$TE-sumnmaday2$seTE)),expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 7") %>% mutate(TRA=1-exp(sumnmaday7$TE), ci_lower=1-exp(sumnmaday7$TE+sumnmaday7$seTE), ci_upper=1-exp(sumnmaday7$TE-sumnmaday7$seTE)), expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 14") %>% mutate(TRA=1-exp(sumnmaday14$TE), ci_lower=1-exp(sumnmaday14$TE+sumnmaday14$seTE), ci_upper=1-exp(sumnmaday14$TE-sumnmaday14$seTE)))  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 


ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \nproportion infected mosquitos")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA/TBA_adjusted_pcr_NMA.png",device="png", width=10, height=7, dpi=600)



##### Adjusted by baseline gametocyte density (mic) #####


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
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7", "day 14")), arm=as.factor(as.character(arm)))

necfeedsub$log = I(log10(necfeedsub$gam0.m+0.1))

studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = necfeedsub[necfeedsub$study %in% currentstudy,]
  
  subdata <- subdata %>%
    rowwise() %>%
    do(data.frame(
      visit = .$visit,
      arm = .$arm,
      studyid_str = .$studyid_str,
      log=.$log,
      prev = c(rep(1, .$mosq_pos), rep(0, .$mosq_total - .$mosq_pos))
    )) %>%
    ungroup()
  
  
  if(nrow(subdata)>0 & length(unique(subdata$visit))>1){
    
    
    mod_main = bam(prev ~ visit:arm-1 + log +s(studyid_str, bs="re"),
                   family=poisson(link="log"), 
                   data=subdata)
    
    vcov_mat <- vcov(mod_main)
    vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
    vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
    mod_main$Vp = vcov_mat 
    
    
    newd=estimate_contrasts(mod_main, contrast=c("visit","arm"))
    newd$study = currentstudy
    newd = newd %>% select(study,Level1,Level2,Difference,SE)
    colnames(newd) = c("study","arm","ref","coef","se")
    
    results = results %>% rbind(newd) %>% na.omit()
    
  }
  
}

# Remove outer brackets from `arm` and `ref`
results$arm <- gsub("^\\(|\\)$", "", results$arm)  # Remove leading and trailing parentheses
results$ref <- gsub("^\\(|\\)$", "", results$ref)  # Remove leading and trailing parentheses

results$se[results$se==0] = 0.001

nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="RR", reference="day 0 AL")

sumnma=summary(nma)$random

arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

#day 2 reductions
sumnmaday2 = list()
sumnmaday2$TE = diag(sumnma$TE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$seTE = diag(sumnma$seTE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$lower = sumnmaday2$TE-1.96*sumnmaday2$seTE
sumnmaday2$upper = sumnmaday2$TE+1.96*sumnmaday2$seTE
outresday2 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday2)=arms
row.names(outresday2)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday2$TE)) {
  for (j in seq_along(sumnmaday2$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday2$TE[i] - sumnmaday2$TE[j]) / sqrt(sumnmaday2$seTE[i]^2 + sumnmaday2$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday2[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday2) = paste0(format(round((1-exp(sumnmaday2$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday2$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday2$lower))*100,2), nsmall=2),"%)")


result_matrixday2 = as.data.frame(outresday2)
result_matrixday2$Reference = factor(row.names(outresday2),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday2 = result_matrixday2 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday2) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TBA_day2_adjusted_mic_NMA.docx")

#day 7 reductions
sumnmaday7 = list()
sumnmaday7$TE = diag(sumnma$TE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$seTE = diag(sumnma$seTE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$lower = sumnmaday7$TE-1.96*sumnmaday7$seTE
sumnmaday7$upper = sumnmaday7$TE+1.96*sumnmaday7$seTE
outresday7 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday7)=arms
row.names(outresday7)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday7$TE)) {
  for (j in seq_along(sumnmaday7$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday7$TE[i] - sumnmaday7$TE[j]) / sqrt(sumnmaday7$seTE[i]^2 + sumnmaday7$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday7[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday7) = paste0(format(round((1-exp(sumnmaday7$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday7$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday7$lower))*100,2), nsmall=2),"%)")


result_matrixday7 = as.data.frame(outresday7)
result_matrixday7$Reference = factor(row.names(outresday7),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday7 = result_matrixday7 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday7) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TBA_day7_adjusted_mic_NMA.docx")

#day 14 reductions

#No feed data for Non-ACT-PQ

arms=arms[!(arms %in% "Non-ACT-PQ")]

sumnmaday14 = list()
sumnmaday14$TE = diag(sumnma$TE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$seTE = diag(sumnma$seTE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$lower = sumnmaday14$TE-1.96*sumnmaday14$seTE
sumnmaday14$upper = sumnmaday14$TE+1.96*sumnmaday14$seTE
outresday14 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday14)=arms
row.names(outresday14)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday14$TE)) {
  for (j in seq_along(sumnmaday14$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday14$TE[i] - sumnmaday14$TE[j]) / sqrt(sumnmaday14$seTE[i]^2 + sumnmaday14$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday14[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday14) = paste0(format(round((1-exp(sumnmaday14$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday14$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday14$lower))*100,2), nsmall=2),"%)")


result_matrixday14 = as.data.frame(outresday14)
result_matrixday14$Reference = factor(row.names(outresday14),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday14 = result_matrixday14 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL",  "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday14) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TBA_day14_adjusted_mic_NMA.docx")


pdat= rbind(expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 2") %>% mutate(TRA=1-exp(sumnmaday2$TE), ci_lower=1-exp(sumnmaday2$TE+sumnmaday2$seTE), ci_upper=1-exp(sumnmaday2$TE-sumnmaday2$seTE)),expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 7") %>% mutate(TRA=1-exp(sumnmaday7$TE), ci_lower=1-exp(sumnmaday7$TE+sumnmaday7$seTE), ci_upper=1-exp(sumnmaday7$TE-sumnmaday7$seTE)), expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 14") %>% mutate(TRA=1-exp(sumnmaday14$TE), ci_lower=1-exp(sumnmaday14$TE+sumnmaday14$seTE), ci_upper=1-exp(sumnmaday14$TE-sumnmaday14$seTE)))  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 


ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \nproportion infected mosquitos")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA/TBA_adjusted_mic_NMA.png",device="png", width=10, height=7, dpi=600)






# Oocyst density ----------------------------------------------------------

#####unadjusted#####

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
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7", "day 14")), arm=as.factor(as.character(arm)))

necfeedsub$log = I(log10(necfeedsub$gam0.m+0.1))

studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = necfeedsub[necfeedsub$study %in% currentstudy,] %>% select(mean_ooc_all,mosq_total, visit, arm, studyid_str) %>% na.omit() 
  
  if(nrow(subdata)>0 & length(unique(subdata$visit))>1){
    
    
    mod_main = bam(mean_ooc_all ~ visit:arm-1 +s(studyid_str, bs="re"),
                   weights = mosq_total,
                   family=poisson(link="log"), 
                   data=subdata)
    
    vcov_mat <- vcov(mod_main)
    vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
    vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
    mod_main$Vp = vcov_mat 
    
    newd=estimate_contrasts(mod_main, contrast=c("visit","arm"))
    newd$study = currentstudy
    newd = newd %>% select(study,Level1,Level2,Difference,SE)
    colnames(newd) = c("study","arm","ref","coef","se")
    
    results = results %>% rbind(newd) %>% na.omit()
    
  }
  
}

# Remove outer brackets from `arm` and `ref`
results$arm <- gsub("^\\(|\\)$", "", results$arm)  # Remove leading and trailing parentheses
results$ref <- gsub("^\\(|\\)$", "", results$ref)  # Remove leading and trailing parentheses

results$se[results$se==0] = 0.001

nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="RR", reference="day 0 AL")

sumnma=summary(nma)$random

arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ")

#day 2 reductions
sumnmaday2 = list()
sumnmaday2$TE = diag(sumnma$TE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$seTE = diag(sumnma$seTE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$lower = sumnmaday2$TE-1.96*sumnmaday2$seTE
sumnmaday2$upper = sumnmaday2$TE+1.96*sumnmaday2$seTE
outresday2 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday2)=arms
row.names(outresday2)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday2$TE)) {
  for (j in seq_along(sumnmaday2$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday2$TE[i] - sumnmaday2$TE[j]) / sqrt(sumnmaday2$seTE[i]^2 + sumnmaday2$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday2[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday2) = paste0(format(round((1-exp(sumnmaday2$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday2$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday2$lower))*100,2), nsmall=2),"%)")


result_matrixday2 = as.data.frame(outresday2)
result_matrixday2$Reference = factor(row.names(outresday2),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"))
result_matrixday2 = result_matrixday2 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday2) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TRA_day2_unadjusted_NMA.docx")

#day 7 reductions
sumnmaday7 = list()
sumnmaday7$TE = diag(sumnma$TE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$seTE = diag(sumnma$seTE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$lower = sumnmaday7$TE-1.96*sumnmaday7$seTE
sumnmaday7$upper = sumnmaday7$TE+1.96*sumnmaday7$seTE
outresday7 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday7)=arms
row.names(outresday7)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday7$TE)) {
  for (j in seq_along(sumnmaday7$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday7$TE[i] - sumnmaday7$TE[j]) / sqrt(sumnmaday7$seTE[i]^2 + sumnmaday7$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday7[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday7) = paste0(format(round((1-exp(sumnmaday7$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday7$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday7$lower))*100,2), nsmall=2),"%)")


result_matrixday7 = as.data.frame(outresday7)
result_matrixday7$Reference = factor(row.names(outresday7),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"))
result_matrixday7 = result_matrixday7 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday7) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TRA_day7_unadjusted_NMA.docx")

#day 14 reductions

#No feed data for Non-ACT-PQ


sumnmaday14 = list()
sumnmaday14$TE = diag(sumnma$TE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$seTE = diag(sumnma$seTE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$lower = sumnmaday14$TE-1.96*sumnmaday14$seTE
sumnmaday14$upper = sumnmaday14$TE+1.96*sumnmaday14$seTE
outresday14 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday14)=arms
row.names(outresday14)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday14$TE)) {
  for (j in seq_along(sumnmaday14$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday14$TE[i] - sumnmaday14$TE[j]) / sqrt(sumnmaday14$seTE[i]^2 + sumnmaday14$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday14[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday14) = paste0(format(round((1-exp(sumnmaday14$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday14$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday14$lower))*100,2), nsmall=2),"%)")


result_matrixday14 = as.data.frame(outresday14)
result_matrixday14$Reference = factor(row.names(outresday14),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"))
result_matrixday14 = result_matrixday14 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL",  "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday14) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TRA_day14_unadjusted_NMA.docx")


pdat= rbind(expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 2") %>% mutate(TRA=1-exp(sumnmaday2$TE), ci_lower=1-exp(sumnmaday2$TE+sumnmaday2$seTE), ci_upper=1-exp(sumnmaday2$TE-sumnmaday2$seTE)),expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 7") %>% mutate(TRA=1-exp(sumnmaday7$TE), ci_lower=1-exp(sumnmaday7$TE+sumnmaday7$seTE), ci_upper=1-exp(sumnmaday7$TE-sumnmaday7$seTE)), expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 14") %>% mutate(TRA=1-exp(sumnmaday14$TE), ci_lower=1-exp(sumnmaday14$TE+sumnmaday14$seTE), ci_upper=1-exp(sumnmaday14$TE-sumnmaday14$seTE)))  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL","Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 

pdat[22:24,]=NA
pdat$arm[22:24]="Non-ACT-PQ"
pdat$TRA[22:24]=NA
pdat$ci_lower[22:24]=NA
pdat$ci_upper[22:24]=NA
pdat$day[22:24]=c("Day 2","Day 7","Day 14")


ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in oocysts\nper dissected mosquito")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA/TRA_unadjusted_NMA.png",device="png", width=10, height=7, dpi=600)


#####adjusted_mic######


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
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7", "day 14")), arm=as.factor(as.character(arm)))

necfeedsub$log = I(log10(necfeedsub$gam0.m+0.1))

studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = necfeedsub[necfeedsub$study %in% currentstudy,] %>% select(mean_ooc_all,mosq_total, visit, arm, studyid_str, log) %>% na.omit() 
  
  if(nrow(subdata)>0 & length(unique(subdata$visit))>1){
    
    
    mod_main = bam(mean_ooc_all ~ visit:arm-1 +log+s(studyid_str, bs="re"),
                   weights = mosq_total,
                   family=poisson(link="log"), 
                   data=subdata)
    
    vcov_mat <- vcov(mod_main)
    vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
    vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
    mod_main$Vp = vcov_mat 
    
    newd=estimate_contrasts(mod_main, contrast=c("visit","arm"))
    newd$study = currentstudy
    newd = newd %>% select(study,Level1,Level2,Difference,SE)
    colnames(newd) = c("study","arm","ref","coef","se")
    
    results = results %>% rbind(newd) %>% na.omit()
    
  }
  
}

# Remove outer brackets from `arm` and `ref`
results$arm <- gsub("^\\(|\\)$", "", results$arm)  # Remove leading and trailing parentheses
results$ref <- gsub("^\\(|\\)$", "", results$ref)  # Remove leading and trailing parentheses

results$se[results$se==0] = 0.001

nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="RR", reference="day 0 AL")

sumnma=summary(nma)$random

arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ")

#day 2 reductions
sumnmaday2 = list()
sumnmaday2$TE = diag(sumnma$TE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$seTE = diag(sumnma$seTE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$lower = sumnmaday2$TE-1.96*sumnmaday2$seTE
sumnmaday2$upper = sumnmaday2$TE+1.96*sumnmaday2$seTE
outresday2 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday2)=arms
row.names(outresday2)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday2$TE)) {
  for (j in seq_along(sumnmaday2$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday2$TE[i] - sumnmaday2$TE[j]) / sqrt(sumnmaday2$seTE[i]^2 + sumnmaday2$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday2[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday2) = paste0(format(round((1-exp(sumnmaday2$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday2$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday2$lower))*100,2), nsmall=2),"%)")


result_matrixday2 = as.data.frame(outresday2)
result_matrixday2$Reference = factor(row.names(outresday2),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"))
result_matrixday2 = result_matrixday2 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday2) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TRA_day2_adjusted_mic_NMA.docx")

#day 7 reductions
sumnmaday7 = list()
sumnmaday7$TE = diag(sumnma$TE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$seTE = diag(sumnma$seTE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$lower = sumnmaday7$TE-1.96*sumnmaday7$seTE
sumnmaday7$upper = sumnmaday7$TE+1.96*sumnmaday7$seTE
outresday7 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday7)=arms
row.names(outresday7)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday7$TE)) {
  for (j in seq_along(sumnmaday7$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday7$TE[i] - sumnmaday7$TE[j]) / sqrt(sumnmaday7$seTE[i]^2 + sumnmaday7$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday7[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday7) = paste0(format(round((1-exp(sumnmaday7$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday7$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday7$lower))*100,2), nsmall=2),"%)")


result_matrixday7 = as.data.frame(outresday7)
result_matrixday7$Reference = factor(row.names(outresday7),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"))
result_matrixday7 = result_matrixday7 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday7) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TRA_day7_adjusted_mic_NMA.docx")

#day 14 reductions

#No feed data for Non-ACT-PQ


sumnmaday14 = list()
sumnmaday14$TE = diag(sumnma$TE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$seTE = diag(sumnma$seTE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$lower = sumnmaday14$TE-1.96*sumnmaday14$seTE
sumnmaday14$upper = sumnmaday14$TE+1.96*sumnmaday14$seTE
outresday14 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday14)=arms
row.names(outresday14)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday14$TE)) {
  for (j in seq_along(sumnmaday14$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday14$TE[i] - sumnmaday14$TE[j]) / sqrt(sumnmaday14$seTE[i]^2 + sumnmaday14$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday14[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday14) = paste0(format(round((1-exp(sumnmaday14$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday14$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday14$lower))*100,2), nsmall=2),"%)")


result_matrixday14 = as.data.frame(outresday14)
result_matrixday14$Reference = factor(row.names(outresday14),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"))
result_matrixday14 = result_matrixday14 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL",  "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday14) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TRA_day14_adjusted_mic_NMA.docx")


pdat= rbind(expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 2") %>% mutate(TRA=1-exp(sumnmaday2$TE), ci_lower=1-exp(sumnmaday2$TE+sumnmaday2$seTE), ci_upper=1-exp(sumnmaday2$TE-sumnmaday2$seTE)),expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 7") %>% mutate(TRA=1-exp(sumnmaday7$TE), ci_lower=1-exp(sumnmaday7$TE+sumnmaday7$seTE), ci_upper=1-exp(sumnmaday7$TE-sumnmaday7$seTE)), expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 14") %>% mutate(TRA=1-exp(sumnmaday14$TE), ci_lower=1-exp(sumnmaday14$TE+sumnmaday14$seTE), ci_upper=1-exp(sumnmaday14$TE-sumnmaday14$seTE)))  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL","Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 

pdat[22:24,]=NA
pdat$arm[22:24]="Non-ACT-PQ"
pdat$TRA[22:24]=NA
pdat$ci_lower[22:24]=NA
pdat$ci_upper[22:24]=NA
pdat$day[22:24]=c("Day 2","Day 7","Day 14")


ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in oocysts\nper dissected mosquito")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA/TRA_adjusted_mic_NMA.png",device="png", width=10, height=7, dpi=600)



#####adjusted_PCR######


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
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7", "day 14")), arm=as.factor(as.character(arm)))

necfeedsub$log = I(log10(necfeedsub$gam0.p+0.01))

studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = necfeedsub[necfeedsub$study %in% currentstudy,] %>% select(mean_ooc_all,mosq_total, visit, arm, studyid_str, log) %>% na.omit() 
  
  if(nrow(subdata)>0 & length(unique(subdata$visit))>1){
    
    
    mod_main = bam(mean_ooc_all ~ visit:arm-1 +log+s(studyid_str, bs="re"),
                   weights = mosq_total,
                   family=poisson(link="log"), 
                   data=subdata)
    
    vcov_mat <- vcov(mod_main)
    vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
    vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
    mod_main$Vp = vcov_mat 
    
    newd=estimate_contrasts(mod_main, contrast=c("visit","arm"))
    newd$study = currentstudy
    newd = newd %>% select(study,Level1,Level2,Difference,SE)
    colnames(newd) = c("study","arm","ref","coef","se")
    
    results = results %>% rbind(newd) %>% na.omit()
    
  }
  
}

# Remove outer brackets from `arm` and `ref`
results$arm <- gsub("^\\(|\\)$", "", results$arm)  # Remove leading and trailing parentheses
results$ref <- gsub("^\\(|\\)$", "", results$ref)  # Remove leading and trailing parentheses

results$se[results$se==0] = 0.001

nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="RR", reference="day 0 AL")

sumnma=summary(nma)$random

arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ")

#day 2 reductions
sumnmaday2 = list()
sumnmaday2$TE = diag(sumnma$TE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$seTE = diag(sumnma$seTE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$lower = sumnmaday2$TE-1.96*sumnmaday2$seTE
sumnmaday2$upper = sumnmaday2$TE+1.96*sumnmaday2$seTE
outresday2 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday2)=arms
row.names(outresday2)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday2$TE)) {
  for (j in seq_along(sumnmaday2$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday2$TE[i] - sumnmaday2$TE[j]) / sqrt(sumnmaday2$seTE[i]^2 + sumnmaday2$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday2[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday2) = paste0(format(round((1-exp(sumnmaday2$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday2$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday2$lower))*100,2), nsmall=2),"%)")


result_matrixday2 = as.data.frame(outresday2)
result_matrixday2$Reference = factor(row.names(outresday2),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"))
result_matrixday2 = result_matrixday2 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday2) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TRA_day2_adjusted_PCR_NMA.docx")

#day 7 reductions
sumnmaday7 = list()
sumnmaday7$TE = diag(sumnma$TE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$seTE = diag(sumnma$seTE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$lower = sumnmaday7$TE-1.96*sumnmaday7$seTE
sumnmaday7$upper = sumnmaday7$TE+1.96*sumnmaday7$seTE
outresday7 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday7)=arms
row.names(outresday7)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday7$TE)) {
  for (j in seq_along(sumnmaday7$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday7$TE[i] - sumnmaday7$TE[j]) / sqrt(sumnmaday7$seTE[i]^2 + sumnmaday7$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday7[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday7) = paste0(format(round((1-exp(sumnmaday7$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday7$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday7$lower))*100,2), nsmall=2),"%)")


result_matrixday7 = as.data.frame(outresday7)
result_matrixday7$Reference = factor(row.names(outresday7),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"))
result_matrixday7 = result_matrixday7 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday7) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TRA_day7_adjusted_PCR_NMA.docx")

#day 14 reductions

#No feed data for Non-ACT-PQ


sumnmaday14 = list()
sumnmaday14$TE = diag(sumnma$TE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$seTE = diag(sumnma$seTE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$lower = sumnmaday14$TE-1.96*sumnmaday14$seTE
sumnmaday14$upper = sumnmaday14$TE+1.96*sumnmaday14$seTE
outresday14 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday14)=arms
row.names(outresday14)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday14$TE)) {
  for (j in seq_along(sumnmaday14$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday14$TE[i] - sumnmaday14$TE[j]) / sqrt(sumnmaday14$seTE[i]^2 + sumnmaday14$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday14[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday14) = paste0(format(round((1-exp(sumnmaday14$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday14$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday14$lower))*100,2), nsmall=2),"%)")


result_matrixday14 = as.data.frame(outresday14)
result_matrixday14$Reference = factor(row.names(outresday14),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"))
result_matrixday14 = result_matrixday14 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL",  "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday14) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/TRA_day14_adjusted_PCR_NMA.docx")


pdat= rbind(expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 2") %>% mutate(TRA=1-exp(sumnmaday2$TE), ci_lower=1-exp(sumnmaday2$TE+sumnmaday2$seTE), ci_upper=1-exp(sumnmaday2$TE-sumnmaday2$seTE)),expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 7") %>% mutate(TRA=1-exp(sumnmaday7$TE), ci_lower=1-exp(sumnmaday7$TE+sumnmaday7$seTE), ci_upper=1-exp(sumnmaday7$TE-sumnmaday7$seTE)), expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 14") %>% mutate(TRA=1-exp(sumnmaday14$TE), ci_lower=1-exp(sumnmaday14$TE+sumnmaday14$seTE), ci_upper=1-exp(sumnmaday14$TE-sumnmaday14$seTE)))  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL","Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 

pdat[22:24,]=NA
pdat$arm[22:24]="Non-ACT-PQ"
pdat$TRA[22:24]=NA
pdat$ci_lower[22:24]=NA
pdat$ci_upper[22:24]=NA
pdat$day[22:24]=c("Day 2","Day 7","Day 14")


ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in oocysts\nper dissected mosquito")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA/TRA_adjusted_PCR_NMA.png",device="png", width=10, height=7, dpi=600)



# Infectiousness ----------------------------------------------------------

#####unadjusted#####

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
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7", "day 14")), arm=as.factor(as.character(arm)))

necfeedsub$log = I(log10(necfeedsub$gam0.m+0.1))

studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = necfeedsub[necfeedsub$study %in% currentstudy,] %>% select(mosq_pos,mosq_total, visit, arm, studyid_str) %>% na.omit() %>% mutate(prev=as.numeric(mosq_pos>0))
  
  if(nrow(subdata)>0 & length(unique(subdata$visit))>1){
    
    
    mod_main = bam(prev ~ visit:arm-1 +s(studyid_str, bs="re"),
                   weights = mosq_total,
                   family=poisson(link="log"), 
                   data=subdata)
    
    vcov_mat <- vcov(mod_main)
    vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
    vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
    mod_main$Vp = vcov_mat 
    
    newd=estimate_contrasts(mod_main, contrast=c("visit","arm"))
    newd$study = currentstudy
    newd = newd %>% select(study,Level1,Level2,Difference,SE)
    colnames(newd) = c("study","arm","ref","coef","se")
    
    results = results %>% rbind(newd) %>% na.omit()
    
  }
  
}

# Remove outer brackets from `arm` and `ref`
results$arm <- gsub("^\\(|\\)$", "", results$arm)  # Remove leading and infectiousnessiling parentheses
results$ref <- gsub("^\\(|\\)$", "", results$ref)  # Remove leading and infectiousnessiling parentheses

results$se[results$se==0] = 0.001

nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="RR", reference="day 0 AL")

sumnma=summary(nma)$random

arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

#day 2 reductions
sumnmaday2 = list()
sumnmaday2$TE = diag(sumnma$TE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$seTE = diag(sumnma$seTE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$lower = sumnmaday2$TE-1.96*sumnmaday2$seTE
sumnmaday2$upper = sumnmaday2$TE+1.96*sumnmaday2$seTE
outresday2 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday2)=arms
row.names(outresday2)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday2$TE)) {
  for (j in seq_along(sumnmaday2$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday2$TE[i] - sumnmaday2$TE[j]) / sqrt(sumnmaday2$seTE[i]^2 + sumnmaday2$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday2[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday2) = paste0(format(round((1-exp(sumnmaday2$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday2$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday2$lower))*100,2), nsmall=2),"%)")


result_matrixday2 = as.data.frame(outresday2)
result_matrixday2$Reference = factor(row.names(outresday2),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday2 = result_matrixday2 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday2) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/infectiousness_day2_unadjusted_NMA.docx")

#day 7 reductions
sumnmaday7 = list()
sumnmaday7$TE = diag(sumnma$TE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$seTE = diag(sumnma$seTE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$lower = sumnmaday7$TE-1.96*sumnmaday7$seTE
sumnmaday7$upper = sumnmaday7$TE+1.96*sumnmaday7$seTE
outresday7 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday7)=arms
row.names(outresday7)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday7$TE)) {
  for (j in seq_along(sumnmaday7$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday7$TE[i] - sumnmaday7$TE[j]) / sqrt(sumnmaday7$seTE[i]^2 + sumnmaday7$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday7[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday7) = paste0(format(round((1-exp(sumnmaday7$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday7$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday7$lower))*100,2), nsmall=2),"%)")


result_matrixday7 = as.data.frame(outresday7)
result_matrixday7$Reference = factor(row.names(outresday7),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday7 = result_matrixday7 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday7) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/infectiousness_day7_unadjusted_NMA.docx")

#day 14 reductions

#No feed data for Non-ACT-PQ

arms=arms[!(arms %in% "Non-ACT-PQ")]

sumnmaday14 = list()
sumnmaday14$TE = diag(sumnma$TE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$seTE = diag(sumnma$seTE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$lower = sumnmaday14$TE-1.96*sumnmaday14$seTE
sumnmaday14$upper = sumnmaday14$TE+1.96*sumnmaday14$seTE
outresday14 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday14)=arms
row.names(outresday14)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday14$TE)) {
  for (j in seq_along(sumnmaday14$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday14$TE[i] - sumnmaday14$TE[j]) / sqrt(sumnmaday14$seTE[i]^2 + sumnmaday14$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday14[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday14) = paste0(format(round((1-exp(sumnmaday14$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday14$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday14$lower))*100,2), nsmall=2),"%)")


result_matrixday14 = as.data.frame(outresday14)
result_matrixday14$Reference = factor(row.names(outresday14),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday14 = result_matrixday14 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL",  "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday14) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/infectiousness_day14_unadjusted_NMA.docx")


pdat= rbind(expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 2") %>% mutate(infectiousness=1-exp(sumnmaday2$TE), ci_lower=1-exp(sumnmaday2$TE+sumnmaday2$seTE), ci_upper=1-exp(sumnmaday2$TE-sumnmaday2$seTE)),expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 7") %>% mutate(infectiousness=1-exp(sumnmaday7$TE), ci_lower=1-exp(sumnmaday7$TE+sumnmaday7$seTE), ci_upper=1-exp(sumnmaday7$TE-sumnmaday7$seTE)), expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 14") %>% mutate(infectiousness=1-exp(sumnmaday14$TE), ci_lower=1-exp(sumnmaday14$TE+sumnmaday14$seTE), ci_upper=1-exp(sumnmaday14$TE-sumnmaday14$seTE)))  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL","Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 


ggplot(data=pdat, aes(x=arm, y=infectiousness, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \nindividual's infection likelihood")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA/infectiousness_unadjusted_NMA.png",device="png", width=10, height=7, dpi=600)




#####adjusted_mic#####

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
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7", "day 14")), arm=as.factor(as.character(arm)))

necfeedsub$log = I(log10(necfeedsub$gam0.m+0.1))

studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = necfeedsub[necfeedsub$study %in% currentstudy,] %>% select(mosq_pos,mosq_total, visit, arm, studyid_str) %>% na.omit() %>% mutate(prev=as.numeric(mosq_pos>0))
  
  if(nrow(subdata)>0 & length(unique(subdata$visit))>1){
    
    
    mod_main = bam(prev ~ visit:arm-1 +s(studyid_str, bs="re"),
                   weights = mosq_total,
                   family=poisson(link="log"), 
                   data=subdata)
    
    vcov_mat <- vcov(mod_main)
    vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
    vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
    mod_main$Vp = vcov_mat 
    
    newd=estimate_contrasts(mod_main, contrast=c("visit","arm"))
    newd$study = currentstudy
    newd = newd %>% select(study,Level1,Level2,Difference,SE)
    colnames(newd) = c("study","arm","ref","coef","se")
    
    results = results %>% rbind(newd) %>% na.omit()
    
  }
  
}

# Remove outer brackets from `arm` and `ref`
results$arm <- gsub("^\\(|\\)$", "", results$arm)  # Remove leading and infectiousnessiling parentheses
results$ref <- gsub("^\\(|\\)$", "", results$ref)  # Remove leading and infectiousnessiling parentheses

results$se[results$se==0] = 0.001

nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="RR", reference="day 0 AL")

sumnma=summary(nma)$random

arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

#day 2 reductions
sumnmaday2 = list()
sumnmaday2$TE = diag(sumnma$TE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$seTE = diag(sumnma$seTE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$lower = sumnmaday2$TE-1.96*sumnmaday2$seTE
sumnmaday2$upper = sumnmaday2$TE+1.96*sumnmaday2$seTE
outresday2 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday2)=arms
row.names(outresday2)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday2$TE)) {
  for (j in seq_along(sumnmaday2$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday2$TE[i] - sumnmaday2$TE[j]) / sqrt(sumnmaday2$seTE[i]^2 + sumnmaday2$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday2[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday2) = paste0(format(round((1-exp(sumnmaday2$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday2$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday2$lower))*100,2), nsmall=2),"%)")


result_matrixday2 = as.data.frame(outresday2)
result_matrixday2$Reference = factor(row.names(outresday2),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday2 = result_matrixday2 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday2) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/infectiousness_day2_adjusted_mic_NMA.docx")

#day 7 reductions
sumnmaday7 = list()
sumnmaday7$TE = diag(sumnma$TE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$seTE = diag(sumnma$seTE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$lower = sumnmaday7$TE-1.96*sumnmaday7$seTE
sumnmaday7$upper = sumnmaday7$TE+1.96*sumnmaday7$seTE
outresday7 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday7)=arms
row.names(outresday7)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday7$TE)) {
  for (j in seq_along(sumnmaday7$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday7$TE[i] - sumnmaday7$TE[j]) / sqrt(sumnmaday7$seTE[i]^2 + sumnmaday7$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday7[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday7) = paste0(format(round((1-exp(sumnmaday7$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday7$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday7$lower))*100,2), nsmall=2),"%)")


result_matrixday7 = as.data.frame(outresday7)
result_matrixday7$Reference = factor(row.names(outresday7),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday7 = result_matrixday7 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday7) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/infectiousness_day7_adjusted_mic_NMA.docx")

#day 14 reductions

#No feed data for Non-ACT-PQ

arms=arms[!(arms %in% "Non-ACT-PQ")]

sumnmaday14 = list()
sumnmaday14$TE = diag(sumnma$TE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$seTE = diag(sumnma$seTE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$lower = sumnmaday14$TE-1.96*sumnmaday14$seTE
sumnmaday14$upper = sumnmaday14$TE+1.96*sumnmaday14$seTE
outresday14 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday14)=arms
row.names(outresday14)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday14$TE)) {
  for (j in seq_along(sumnmaday14$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday14$TE[i] - sumnmaday14$TE[j]) / sqrt(sumnmaday14$seTE[i]^2 + sumnmaday14$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday14[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday14) = paste0(format(round((1-exp(sumnmaday14$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday14$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday14$lower))*100,2), nsmall=2),"%)")


result_matrixday14 = as.data.frame(outresday14)
result_matrixday14$Reference = factor(row.names(outresday14),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday14 = result_matrixday14 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL",  "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday14) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/infectiousness_day14_adjusted_mic_NMA.docx")


pdat= rbind(expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 2") %>% mutate(infectiousness=1-exp(sumnmaday2$TE), ci_lower=1-exp(sumnmaday2$TE+sumnmaday2$seTE), ci_upper=1-exp(sumnmaday2$TE-sumnmaday2$seTE)),expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 7") %>% mutate(infectiousness=1-exp(sumnmaday7$TE), ci_lower=1-exp(sumnmaday7$TE+sumnmaday7$seTE), ci_upper=1-exp(sumnmaday7$TE-sumnmaday7$seTE)), expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 14") %>% mutate(infectiousness=1-exp(sumnmaday14$TE), ci_lower=1-exp(sumnmaday14$TE+sumnmaday14$seTE), ci_upper=1-exp(sumnmaday14$TE-sumnmaday14$seTE)))  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL","Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 


ggplot(data=pdat, aes(x=arm, y=infectiousness, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \nindividual's infection likelihood")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA/infectiousness_adjusted_mic_NMA.png",device="png", width=10, height=7, dpi=600)





#####adjusted_PCR#####

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
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7", "day 14")), arm=as.factor(as.character(arm)))

necfeedsub$log = I(log10(necfeedsub$gam0.p+0.01))

studies = c("NECTAR1", "NECTAR2", "NECTAR3", "NECTAR4", "PQ01", "PQ03")

results = data.frame(study=NA,arm=NA, ref=NA,coef=NA,se=NA)

for (i in 1:length(studies)){
  currentstudy = studies[i]
  
  subdata = necfeedsub[necfeedsub$study %in% currentstudy,] %>% select(mosq_pos,mosq_total, visit, arm, studyid_str) %>% na.omit() %>% mutate(prev=as.numeric(mosq_pos>0))
  
  if(nrow(subdata)>0 & length(unique(subdata$visit))>1){
    
    
    mod_main = bam(prev ~ visit:arm-1 +s(studyid_str, bs="re"),
                   weights = mosq_total,
                   family=poisson(link="log"), 
                   data=subdata)
    
    vcov_mat <- vcov(mod_main)
    vcov_mat[abs(vcov_mat) > 5] <- sign(vcov_mat[abs(vcov_mat) > 5]) * 5
    vcov_mat[abs(vcov_mat) == 0] <- sign(vcov_mat[abs(vcov_mat) == 0]) + 0.001^2
    mod_main$Vp = vcov_mat 
    
    newd=estimate_contrasts(mod_main, contrast=c("visit","arm"))
    newd$study = currentstudy
    newd = newd %>% select(study,Level1,Level2,Difference,SE)
    colnames(newd) = c("study","arm","ref","coef","se")
    
    results = results %>% rbind(newd) %>% na.omit()
    
  }
  
}

# Remove outer brackets from `arm` and `ref`
results$arm <- gsub("^\\(|\\)$", "", results$arm)  # Remove leading and infectiousnessiling parentheses
results$ref <- gsub("^\\(|\\)$", "", results$ref)  # Remove leading and infectiousnessiling parentheses

results$se[results$se==0] = 0.001

nma = netmeta(TE = results$coef, seTE=results$se, treat1=results$arm, treat2=results$ref, studlab=results$study, sm="RR", reference="day 0 AL")

sumnma=summary(nma)$random

arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

#day 2 reductions
sumnmaday2 = list()
sumnmaday2$TE = diag(sumnma$TE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$seTE = diag(sumnma$seTE[paste("day 2",arms),paste("day 0",arms)])
sumnmaday2$lower = sumnmaday2$TE-1.96*sumnmaday2$seTE
sumnmaday2$upper = sumnmaday2$TE+1.96*sumnmaday2$seTE
outresday2 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday2)=arms
row.names(outresday2)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday2$TE)) {
  for (j in seq_along(sumnmaday2$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday2$TE[i] - sumnmaday2$TE[j]) / sqrt(sumnmaday2$seTE[i]^2 + sumnmaday2$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday2[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday2) = paste0(format(round((1-exp(sumnmaday2$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday2$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday2$lower))*100,2), nsmall=2),"%)")


result_matrixday2 = as.data.frame(outresday2)
result_matrixday2$Reference = factor(row.names(outresday2),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday2 = result_matrixday2 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday2) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/infectiousness_day2_adjusted_PCR_NMA.docx")

#day 7 reductions
sumnmaday7 = list()
sumnmaday7$TE = diag(sumnma$TE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$seTE = diag(sumnma$seTE[paste("day 7",arms),paste("day 0",arms)])
sumnmaday7$lower = sumnmaday7$TE-1.96*sumnmaday7$seTE
sumnmaday7$upper = sumnmaday7$TE+1.96*sumnmaday7$seTE
outresday7 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday7)=arms
row.names(outresday7)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday7$TE)) {
  for (j in seq_along(sumnmaday7$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday7$TE[i] - sumnmaday7$TE[j]) / sqrt(sumnmaday7$seTE[i]^2 + sumnmaday7$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday7[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday7) = paste0(format(round((1-exp(sumnmaday7$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday7$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday7$lower))*100,2), nsmall=2),"%)")


result_matrixday7 = as.data.frame(outresday7)
result_matrixday7$Reference = factor(row.names(outresday7),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday7 = result_matrixday7 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday7) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/infectiousness_day7_adjusted_PCR_NMA.docx")

#day 14 reductions

#No feed data for Non-ACT-PQ

arms=arms[!(arms %in% "Non-ACT-PQ")]

sumnmaday14 = list()
sumnmaday14$TE = diag(sumnma$TE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$seTE = diag(sumnma$seTE[paste("day 14",arms),paste("day 0",arms)])
sumnmaday14$lower = sumnmaday14$TE-1.96*sumnmaday14$seTE
sumnmaday14$upper = sumnmaday14$TE+1.96*sumnmaday14$seTE
outresday14 = matrix(NA, nrow=length(arms),ncol=length(arms))
colnames(outresday14)=arms
row.names(outresday14)=arms
# Compute p-values for all pairwise differences
for (i in seq_along(sumnmaday14$TE)) {
  for (j in seq_along(sumnmaday14$TE)) {
    if (i != j) {
      # Calculate z-score
      z <- (sumnmaday14$TE[i] - sumnmaday14$TE[j]) / sqrt(sumnmaday14$seTE[i]^2 + sumnmaday14$seTE[j]^2)
      
      # Calculate two-sided p-value
      p <- 2 * (1 - pnorm(abs(z)))
      
      outresday14[i, j] = paste0("p",ifelse(p<0.0001,"<0.0001",ifelse(p>0.9999,">0.9999",paste0("=",format(round(p,4), nsmall=4)))))
    }
  }
}

diag(outresday14) = paste0(format(round((1-exp(sumnmaday14$TE))*100,2),nsmall=2),"% (",format(round((1-exp(sumnmaday14$upper))*100,2), nsmall=2),"% - ",format(round((1-exp(sumnmaday14$lower))*100,2), nsmall=2),"%)")


result_matrixday14 = as.data.frame(outresday14)
result_matrixday14$Reference = factor(row.names(outresday14),levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
result_matrixday14 = result_matrixday14 %>% select(c("Reference","DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL",  "ACT-PQ", "ACT-TQ")) %>% arrange(Reference)

flextable(result_matrixday14) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA/infectiousness_day14_adjusted_PCR_NMA.docx")


pdat= rbind(expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 2") %>% mutate(infectiousness=1-exp(sumnmaday2$TE), ci_lower=1-exp(sumnmaday2$TE+sumnmaday2$seTE), ci_upper=1-exp(sumnmaday2$TE-sumnmaday2$seTE)),expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), day="Day 7") %>% mutate(infectiousness=1-exp(sumnmaday7$TE), ci_lower=1-exp(sumnmaday7$TE+sumnmaday7$seTE), ci_upper=1-exp(sumnmaday7$TE-sumnmaday7$seTE)), expand.grid(arm=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"), day="Day 14") %>% mutate(infectiousness=1-exp(sumnmaday14$TE), ci_lower=1-exp(sumnmaday14$TE+sumnmaday14$seTE), ci_upper=1-exp(sumnmaday14$TE-sumnmaday14$seTE)))  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL","Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 


ggplot(data=pdat, aes(x=arm, y=infectiousness, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \nindividual's infection likelihood")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA/infectiousness_adjusted_PCR_NMA.png",device="png", width=10, height=7, dpi=600)

