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

#Exclude 6 participants without baseline infectivity or gametocyte measures
nec <- nec %>%
  filter(!studyid_str %in% c("PQ-01-001","PQ-01-004","PQ-01-005","PQ-01-006","PQ-01-008","PQ-01-052"))


nec = nec %>% arrange(studyid_str, studyvisit_num)


# Descriptive results -----------------------------------------------------



#We create a dataset for the baseline data [i.e. at visit 0] and make a summarising function for the descriptives.

nec0 = nec[nec$studyvisit_num==0,]
mese = function(x,y=1){
  paste0(format(round(mean(x, na.rm=TRUE),y),nsmall=y), " (",format(round(sd(x, na.rm=TRUE),y),nsmall=y),")")
}

mesel = function(x,y=1){
  if(sum(x==0, na.rm=TRUE)>0){x[x==0]=0.1} else {}
  
  paste0(format(round(mean(log(x), na.rm=TRUE),y),nsmall=y), " (",format(round(sd(log(x),na.rm=TRUE),y),nsmall=y),")")
}

mediqr = function(x,y=1){
  paste0(format(round(median(x, na.rm=TRUE),y),nsmall=y), " (",format(round(quantile(x,probs = c(0.25),na.rm=TRUE)[[1]],y),nsmall=y)," - ",format(round(quantile(x,probs = c(0.75),na.rm=TRUE)[[1]],y),nsmall=y),")")
}

percy = function(x){
  paste0(table(x)," (",format(round(prop.table(table(x))*100,1),nsmall=1),"%)")
}
percyt = function(x,y){
  paste0(table(nec0$arm,x==y)[,2]," (",format(round(prop.table(table(nec0$arm,x==y),1)[,2]*100,1),nsmall=1),"%)")
}

percyg = function(x,y){
  paste0(table(nec0$arm,x==y)[,1]," (",format(round(prop.table(table(nec0$arm,x==y),1)[,1]*100,1),nsmall=1),"%)")
}

nec0$parpos = as.numeric(nec0$mic_pf_asc_ul>0)
nec0$gampos =as.numeric(nec0$mic_pf_gct_ul>0)

#Next we need to create a dataframe table for our treatment counts and percentages.

out1 = data.frame(
  treat = levels(nec0$arm),
  counts = percy(nec0$arm))
colnames(out1) = c("Treatment","N (%)")

#Further, we want to see the distribution of several other variables over the different treatments.

out2  = nec0[order(nec$arm),] %>% 
  group_by(arm) %>% 
  summarise(age = mese(age,1),
            temp = mese(temp,1),
            logpar = mediqr(tmpd,1),
            loggam1 = mediqr(mic_pf_gct_ul,1),
            loggam2 = mediqr(qrtpcr_gct_ul,1),
            loggamf = mediqr(qrtpcr_femgct_ul,1),
            proppos = mediqr(mosq_pos_prop,1)) %>% 
  dplyr::select(arm, age, temp, logpar, loggam1, loggam2, loggamf, proppos)

out2=na.omit(out2)
colnames(out2)= c("Treatment","Age", "Temperature", "Total parasites/uL", "Gametocytes/uL (microscopy)", "Gametocytes/uL (pcr)", "Female gametocytes/uL (pcr)", "Proportion infected")

out3 = data.frame(
  treat = levels(nec0$arm),
  counts = percyt(nec0$sex_bin,1))
colnames(out3) = c("Treatment","Males")

out4 = data.frame(
  treat = levels(nec0$arm),
  counts = percyt(nec0$gampos,1))
colnames(out4) = c("Treatment","Gametocyte prevalence")

out5 = data.frame(
  treat = levels(nec0$arm),
  counts = percyt(nec0$parpos,1))
colnames(out5) = c("Treatment","Asexual parasite prevalence")

out6 = data.frame(
  treat = levels(nec0$arm),
  counts = percyt(nec0$person_pos,1))
colnames(out6) = c("Treatment","Proportion participants infectious")

out.c = out1 %>% left_join(
  out3 %>% left_join(out2, by="Treatment"),
  by="Treatment"
) %>% left_join(out4 %>% left_join(out5 %>% left_join(out6, by="Treatment"), by="Treatment"), by="Treatment")

flextable(out.c) %>% autofit()  %>% save_as_docx( path = "OUTPUT/descriptives.docx")



#We create a dataset for the baseline data [i.e. at visit 2] and make a summarising function for the descriptives.

nec2 = nec[nec$studyvisit_num==2,]
mese = function(x,y=1){
  paste0(format(round(mean(x, na.rm=TRUE),y),nsmall=y), " (",format(round(sd(x, na.rm=TRUE),y),nsmall=y),")")
}

mesel = function(x,y=1){
  if(sum(x==0, na.rm=TRUE)>0){x[x==0]=0.1} else {}
  
  paste0(format(round(mean(log(x), na.rm=TRUE),y),nsmall=y), " (",format(round(sd(log(x),na.rm=TRUE),y),nsmall=y),")")
}

mediqr = function(x,y=1){
  paste0(format(round(median(x, na.rm=TRUE),y),nsmall=y), " (",format(round(quantile(x,probs = c(0.25),na.rm=TRUE)[[1]],y),nsmall=y)," - ",format(round(quantile(x,probs = c(0.75),na.rm=TRUE)[[1]],y),nsmall=y),")")
}

percy = function(x){
  paste0(table(x)," (",format(round(prop.table(table(x))*100,1),nsmall=1),"%)")
}
percyt = function(x,y){
  paste0(table(nec2$arm,x==y)[,2]," (",format(round(prop.table(table(nec2$arm,x==y),1)[,2]*100,1),nsmall=1),"%)")
}

percyg = function(x,y){
  paste0(table(nec2$arm,x==y)[,1]," (",format(round(prop.table(table(nec2$arm,x==y),1)[,1]*100,1),nsmall=1),"%)")
}

nec2$parpos = as.numeric(nec2$mic_pf_asc_ul>0)
nec2$gampos =as.numeric(nec2$mic_pf_gct_ul>0)

#Next we need to create a dataframe table for our treatment counts and percentages.

out1 = data.frame(
  treat = levels(nec2$arm),
  counts = percy(nec2$arm))
colnames(out1) = c("Treatment","N (%)")

#Further, we want to see the distribution of several other variables over the different treatments.

out2  = nec2[order(nec$arm),] %>% 
  group_by(arm) %>% 
  summarise(age = mese(age,1),
            temp = mese(temp,1),
            logpar = mediqr(tmpd,1),
            loggam1 = mediqr(mic_pf_gct_ul,1),
            loggam2 = mediqr(qrtpcr_gct_ul,1),
            loggamf = mediqr(qrtpcr_femgct_ul,1),
            proppos = mediqr(mosq_pos_prop,1)) %>% 
  dplyr::select(arm, age, temp, logpar, loggam1, loggam2, loggamf, proppos)

out2=na.omit(out2)
colnames(out2)= c("Treatment","Age", "Temperature", "Total parasites/uL", "Gametocytes/uL (microscopy)", "Gametocytes/uL (pcr)", "Female gametocytes/uL (pcr)", "Proportion infected")

out3 = data.frame(
  treat = levels(nec2$arm),
  counts = percyt(nec2$sex_bin,1))
colnames(out3) = c("Treatment","Males")

out4 = data.frame(
  treat = levels(nec2$arm),
  counts = percyt(nec2$gampos,1))
colnames(out4) = c("Treatment","Gametocyte prevalence")

out5 = data.frame(
  treat = levels(nec2$arm),
  counts = percyt(nec2$parpos,1))
colnames(out5) = c("Treatment","Asexual parasite prevalence")

out6 = data.frame(
  treat = levels(nec2$arm),
  counts = percyt(nec2$person_pos,1))
colnames(out6) = c("Treatment","Proportion participants infectious")

out.c = out1 %>% left_join(
  out3 %>% left_join(out2, by="Treatment"),
  by="Treatment"
) %>% left_join(out4 %>% left_join(out5 %>% left_join(out6, by="Treatment"), by="Treatment"), by="Treatment")

flextable(out.c) %>% autofit()  %>% save_as_docx( path = "OUTPUT/descriptives_day2.docx")



# Baseline figures by year ------------------------------------------------



library(ggplot2)
library(ggprism)
library(ggpubr)
library(epitools)

nec0s = nec0 %>% 
  group_by(year) %>% 
  summarise(n=n(),count=sum(parpos),prev = mean(parpos)) %>% 
  mutate(lwr = binom.exact(count, n, conf.level = 0.95)[[4]], upr=binom.exact(count, n, conf.level = 0.95)[[5]])


ggplot(data=nec0s, aes(x=as.factor(year), y=prev,fill=as.factor(year)))+
  geom_bar(stat="identity",alpha=0.5)+
  theme_prism()+
  ylab("Asexual parasite prevalence")+
  xlab("Cohort year")+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,1,by=0.2), limits=c(0,1))+
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.1, linewidth=1.2)+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_1.png",device="png", width=12, height=8, dpi=600)

options(scipen=999)
ggplot(data=nec0,aes(x=as.factor(year),y=mic_pf_asc_ul,fill=as.factor(year)))+
  geom_violin(alpha=0.5)+
  geom_jitter(alpha=0.1, width=0.2)+
  scale_y_log10()+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  ylab("Asexual parasite density/µL")+
  xlab("Cohort year")+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))  

ggsave("OUTPUT/Fig_2.png",device="png", width=9, height=6, dpi=600)


#gametocyte positives 
nec0g = nec0 %>% 
  group_by(year) %>% 
  summarise(n=sum(!is.na(gampos)),count=sum(gampos,na.rm=TRUE),prev = mean(gampos,na.rm=TRUE)) %>% 
  mutate(lwr = binom.exact(count, n, conf.level = 0.95)[[4]], upr=binom.exact(count, n, conf.level = 0.95)[[5]])


ggplot(data=nec0g, aes(x=as.factor(year), y=prev,fill=as.factor(year)))+
  geom_bar(stat="identity",alpha=0.5)+
  theme_prism()+
  ylab("Gametocyte prevalence")+
  xlab("Cohort year")+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,1,by=0.2), limits=c(0,1))+
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.1, linewidth=1.2)+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_3.png",device="png", width=12, height=8, dpi=600)


options(scipen=999)
ggplot(data=nec0,aes(x=as.factor(year),y=mic_pf_gct_ul,fill=as.factor(year)))+
  geom_violin(alpha=0.5)+
  geom_jitter(alpha=0.1, width=0.2)+
  scale_y_log10()+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  ylab("Gametocyte density/µL (microscopy)")+
  xlab("Cohort year")+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))  

ggsave("OUTPUT/Fig_4.png",device="png", width=9, height=6, dpi=600)


ggplot(data=nec0,aes(x=as.factor(year),y=qrtpcr_gct_ul,fill=as.factor(year)))+
  geom_violin(alpha=0.5)+
  geom_jitter(alpha=0.1, width=0.2)+
  scale_y_log10()+
  theme_prism()+
  ylab("Gametocyte density pcr/µL")+
  xlab("Cohort year")+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_5.png",device="png", width=12, height=8, dpi=600)


#mean infectivity over the years

library(mgcv)

modi = gam(mosq_pos_prop/100 ~ year, family=binomial(), data=nec0)
summary(modi)


nec0i = expand.grid(year=sort(unique((nec0$year))))
ilink <- family(modi)$linkinv
## add fit and se.fit on the **link** scale
nec0i <- bind_cols(nec0i, setNames(as_tibble(predict(modi, nec0i, se.fit = TRUE)[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i <- mutate(nec0i,
                prev  = ilink(fit_link),
                upr = ilink(fit_link + (2 * se_link)),
                lwr = ilink(fit_link - (2 * se_link)))


ggplot(data=nec0i, aes(x=as.factor(year), y=prev,fill=as.factor(year)))+
  geom_bar(stat="identity",alpha=0.5)+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  ylab("Proportion infected mosquitoes")+
  xlab("Cohort year")+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,0.3,by=0.05), limits=c(0,0.35))+
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.1, linewidth=0.6)+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))  


ggsave("OUTPUT/Fig_6.png",device="png", width=9, height=6, dpi=600)



#baseline infectivity vs gametocyte density (pcr) 

modi = gam(mosq_pos_prop/100 ~ I(log(gam.p))+s(year, bs="re"), family=binomial(), data=nec0[nec0$year!="2016",])
summary(modi)

nec0i = expand.grid(gam.p=seq(0.01,max(nec0$gam.p,na.rm=T), length.out=10000), year=nec0$year[1])
ilink <- family(modi)$linkinv
## add fit and se.fit on the **link** scale
nec0i <- bind_cols(nec0i, setNames(as_tibble(predict(modi, nec0i, se.fit = TRUE,exclude = c("s(year)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i <- mutate(nec0i,
                prev  = ilink(fit_link),
                upr = ilink(fit_link + (2 * se_link)),
                lwr = ilink(fit_link - (2 * se_link)))


modi2 = gam(mosq_pos_prop/100 ~ I(log(gam.p))*year, family=binomial(), data=nec0)
summary(modi2)

nec0i2 = expand.grid(gam.p=seq(0.01,max(nec0$gam.p,na.rm=T), length.out=10000),year=sort(unique((nec0$year))))
ilink <- family(modi2)$linkinv
## add fit and se.fit on the **link** scale
nec0i2 <- bind_cols(nec0i2, setNames(as_tibble(predict(modi2, nec0i2, se.fit = TRUE)[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i2 <- mutate(nec0i2,
                prev  = ilink(fit_link),
                upr = ilink(fit_link + (2 * se_link)),
                lwr = ilink(fit_link - (2 * se_link)))


ggplot(data=nec0i, aes(x=gam.p, y=prev))+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2)+
  geom_line(data=nec0i2, aes(x=gam.p, y=prev, col=year),linewidth=0.8,alpha=0.7)+
  geom_line(linewidth=1.2)+
  theme_prism()+
  ylab("Proportion infected mosquitoes")+
  xlab("Gametocyte density/µL (pcr)")+
  scale_x_log10()+
  scale_color_viridis_d()+
  theme(legend.position = "top")+
  scale_y_continuous(labels = scales::percent)

ggsave("OUTPUT/Fig_6b.png",device="png", width=12, height=8, dpi=600)



#baseline infectivity vs gametocyte density (microscopy) 

modi = gam(mosq_pos_prop/100 ~ I(log(gam.m+0.1))+s(year, bs="re"), family=binomial(), data=nec0)
summary(modi)

nec0i = expand.grid(gam.m=seq(0,max(nec0$gam.m,na.rm=T), length.out=10000), year=nec0$year[1])
ilink <- family(modi)$linkinv
## add fit and se.fit on the **link** scale
nec0i <- bind_cols(nec0i, setNames(as_tibble(predict(modi, nec0i, se.fit = TRUE,exclude = c("s(year)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i <- mutate(nec0i,
                prev  = ilink(fit_link),
                upr = ilink(fit_link + (2 * se_link)),
                lwr = ilink(fit_link - (2 * se_link)))


modi2 = gam(mosq_pos_prop/100 ~ I(log(gam.m+0.1))*year, family=binomial(), data=nec0)
summary(modi2)

nec0i2 = expand.grid(gam.m=seq(0,max(nec0$gam.m,na.rm=T), length.out=10000),year=sort(unique((nec0$year))))
ilink <- family(modi2)$linkinv
## add fit and se.fit on the **link** scale
nec0i2 <- bind_cols(nec0i2, setNames(as_tibble(predict(modi2, nec0i2, se.fit = TRUE)[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i2 <- mutate(nec0i2,
                 prev  = ilink(fit_link),
                 upr = ilink(fit_link + (2 * se_link)),
                 lwr = ilink(fit_link - (2 * se_link)))


ggplot(data=nec0i, aes(x=gam.m, y=prev))+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2)+
  geom_line(data=nec0i2, aes(x=gam.m, y=prev, col=year),linewidth=1.05,alpha=0.7)+
  geom_line(linewidth=1.2)+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  ylab("Proportion infected mosquitoes")+
  xlab("Gametocyte density/µL (microscopy)")+
  scale_x_log10()+
  scale_color_viridis_d()+
  theme(legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))+
  scale_y_continuous(labels = scales::percent)

ggsave("OUTPUT/Fig_6c.png",device="png", width=8, height=6, dpi=600)


#correspondence between gametocyte densities by microscopy and by PCR

modi = gam(I(log10(gam.p+0.01)) ~ I(log10(gam.m+0.1))+s(year, bs="re"), data=nec0)
summary(modi)

nec0i = expand.grid(gam.m=seq(min(nec0$gam.m[nec0$gam.m>0],na.rm=T),max(nec0$gam.m,na.rm=T), length.out=10000), year=nec0$year[1])
ilink <- family(modi)$linkinv
## add fit and se.fit on the **link** scale
nec0i <- bind_cols(nec0i, setNames(as_tibble(predict(modi, nec0i, se.fit = TRUE,exclude = c("s(year)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i <- mutate(nec0i,
                gam.p  = 10^(fit_link),
                upr = 10^(fit_link + (2 * se_link)),
                lwr = 10^(fit_link - (2 * se_link)))



modi2 = gam(I(log10(gam.p+0.01)) ~ I(log10(gam.m+0.1))*year, data=nec0)
summary(modi2)
smi2 = summary(modi2)

new_labels = paste0(levels(nec0$year)," (\u03B2=",round(as.vector(c(smi2$p.coeff["I(log10(gam.m + 0.1))"], smi2$p.coeff["I(log10(gam.m + 0.1))"]+smi2$p.coeff[c("I(log10(gam.m + 0.1)):year2016","I(log10(gam.m + 0.1)):year2019","I(log10(gam.m + 0.1)):year2020","I(log10(gam.m + 0.1)):year2021","I(log10(gam.m + 0.1)):year2022")])),2),")")

nec0i2 = expand.grid(gam.m=seq(min(nec0$gam.m[nec0$gam.m>0],na.rm=T),max(nec0$gam.m,na.rm=T), length.out=10000), study=nec0$study[1], year=unique(sort(nec0$year)))
ilink <- family(modi2)$linkinv
## add fit and se.fit on the **link** scale
nec0i2 <- bind_cols(nec0i2, setNames(as_tibble(predict(modi2, nec0i2, se.fit = TRUE,exclude = c("s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i2 <- mutate(nec0i2,
                 gam.p  = 10^(fit_link),
                 upr = 10^(fit_link + (2 * se_link)),
                 lwr = 10^(fit_link - (2 * se_link)))


ggplot(data=nec0i, aes(x=gam.m, y=gam.p))+
  geom_point(data=nec0, aes(x=gam.m, y=gam.p, col=year),alpha=0.4)+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2)+
  geom_line(data=nec0i2, aes(x=gam.m,y=gam.p,col=year),linewidth=0.8, alpha=0.7)+
  geom_line(linewidth=1.2)+
  theme_prism(base_size = 14, base_line_size = 0.1)+
  ylab("Gametocyte density/µL (pcr)")+
  xlab("Gametocyte density/µL (microscopy)")+
  scale_x_log10(breaks=c(1,10,100,1000,10000,100000),limits=c(1,100000))+
  scale_y_log10(breaks=c(1,10,100,1000,10000,100000),limits=c(1,100000))+
  scale_color_viridis_d(labels = new_labels)+
  theme(legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  geom_abline(slope=1, intercept=0, linetype="dashed")

ggsave("OUTPUT/Fig_6d.png",device="png", width=9, height=6, dpi=600)


#variation in proportion infected mosquitoes + contribution of submicroscopic gametocytes
#BASELINE

# Round up mosq_pos_prop to the nearest integer
nec0_plot<-nec0
nec0_plot$mosq_pos_prop_bin <- ceiling(nec0_plot$mosq_pos_prop)

# Create a new variable for color based on mic_pf_gct_ul
nec0_plot$color <- ifelse(nec0_plot$mic_pf_gct_ul > 0, "Positive", "Negative")
nec0_plot$color <- factor(nec0_plot$color, levels = c("Negative", "Positive"))
nec0_subset<-nec0_plot[!is.na(nec0_plot$mosq_pos_prop_bin) & !is.na(nec0_plot$color),]
nec0_subset %>%
  group_by(study) %>%
  summarise(unique_studyids = n_distinct(studyid_str))
nec0_subset<-nec0_subset[nec0_subset$mosq_pos_prop_bin>1,]
p<-ggplot(nec0_subset, aes(x = mosq_pos_prop_bin, fill = color)) +
  geom_bar(stat = "count", position = position_stack(reverse = TRUE),  aes(y = after_stat(count) + 0.1)) +
  labs(x = "Mosq Pos Prop", y = "Day 0", fill = "Mic_pf_gct_ul") +
  scale_fill_manual(values = c("Negative" = "#e28743", "Positive" = "#94bbc6"),
                    labels = c("Negative" = "Microscopy-, PCR+ ", "Positive" = "Microscopy+")) + 
  theme_prism(base_size = 14, base_line_size = 0.2) +
  scale_y_log10() +
  guides(fill = guide_legend(title = "Gametocytes")) +
  scale_x_continuous(breaks = c(1, 20, 40, 60, 80,100),limits=c(1,100))+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))




#Day2
# Round up mosq_pos_prop to the nearest integer
nec2_plot<-nec2
nec2_plot$mosq_pos_prop_bin <- ceiling(nec2_plot$mosq_pos_prop)

# Create a new variable for color based on mic_pf_gct_ul
nec2_plot$color <- ifelse(nec2_plot$mic_pf_gct_ul > 0, "Positive", "Negative")
nec2_plot$color <- factor(nec2_plot$color, levels = c("Negative", "Positive"))
nec2_subset<-nec2_plot[!is.na(nec2_plot$mosq_pos_prop_bin) & !is.na(nec2_plot$color),]
nec2_subset %>%
  group_by(study) %>%
  summarise(unique_studyids = n_distinct(studyid_str))
nec2_subset<-nec2_subset[nec2_subset$mosq_pos_prop_bin>1,]
q<-ggplot(nec2_subset, aes(x = mosq_pos_prop_bin, fill = color)) +
  geom_bar(stat = "count", position = position_stack(reverse = TRUE),  aes(y = after_stat(count) + 0.1)) +
  labs(x = "Mosq Pos Prop", y = "Day 2", fill = "Mic_pf_gct_ul") +
  scale_fill_manual(values = c("Negative" = "#e28743","Positive" = "#94bbc6")) + 
  theme_prism(base_size = 14, base_line_size = 0.2) +
  scale_y_log10() + 
  scale_x_continuous(breaks = c(1, 20, 40, 60, 80,100),limits=c(1,100))+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))


#Day7
nec7 = nec[nec$studyvisit_num==7,]
nec7$mosq_pos_prop_bin <- ceiling(nec7$mosq_pos_prop)
# Create a new variable for color based on mic_pf_gct_ul
nec7$color <- ifelse(nec7$mic_pf_gct_ul > 0, "Positive", "Negative")
nec7$color <- factor(nec7$color, levels = c("Negative", "Positive"))
nec7_subset<-nec7[!is.na(nec7$mosq_pos_prop_bin) & !is.na(nec7$color),]
nec7_subset %>%
  group_by(study) %>%
  summarise(unique_studyids = n_distinct(studyid_str))
nec7_subset<-nec7_subset[nec7_subset$mosq_pos_prop_bin>1,]
r<-ggplot(nec7_subset, aes(x = mosq_pos_prop_bin, fill = color)) +
  geom_bar(stat = "count", position = position_stack(reverse = TRUE),  aes(y = after_stat(count) + 0.1)) +
  labs(x = "Mosq Pos Prop", y = "Day 7", fill = "Mic_pf_gct_ul") +
  scale_fill_manual(values = c("Negative" = "#e28743","Positive" = "#94bbc6")) + 
  theme_prism(base_size = 14, base_line_size = 0.2) +
  scale_y_log10() +
  scale_x_continuous(breaks = c(1, 20, 40, 60, 80,100),limits=c(1,100))+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))



#Day14
nec14 = nec[nec$studyvisit_num==14,]
nec14$mosq_pos_prop_bin <- ceiling(nec14$mosq_pos_prop)
# Create a new variable for color based on mic_pf_gct_ul
nec14$color <- ifelse(nec14$mic_pf_gct_ul > 0, "Positive", "Negative")
nec14$color <- factor(nec14$color, levels = c("Negative", "Positive"))
nec14_subset<-nec14[!is.na(nec14$mosq_pos_prop_bin) & !is.na(nec14$color),]
nec14_subset %>%
  group_by(study) %>%
  summarise(unique_studyids = n_distinct(studyid_str))
nec14_subset<-nec14_subset[nec14_subset$mosq_pos_prop_bin>1,]
s<-ggplot(nec14_subset, aes(x = mosq_pos_prop_bin, fill = color)) +
  geom_bar(stat = "count", position = position_stack(reverse = TRUE),  aes(y = after_stat(count) + 0.05)) +
  labs(x = "Proportion infected mosquitoes (%)", y = "Day 14", fill = "Mic_pf_gct_ul") +
  scale_fill_manual(values = c("Negative" = "#e28743","Positive" = "#94bbc6")) + 
  theme_prism(base_size = 14, base_line_size = 0.2) +
  scale_y_log10() +
  scale_x_continuous(breaks = c(1, 20, 40, 60, 80,100),limits=c(1,100))+
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))


MIRvariation <- p / q / r / s


ggsave("OUTPUT/Fig_6e.png",device="png", width=12, height=8, dpi=600)



#by year and by month (are earlier infections more infectious?)

modi = gam(mosq_pos_prop/100 ~ year*as.factor(month), family=binomial(), data=nec0)
summary(modi)


nec0i = na.omit(nec0 %>% select(year, month) %>% distinct() %>% arrange(year, month))
ilink <- family(modi)$linkinv
## add fit and se.fit on the **link** scale
nec0i <- bind_cols(nec0i, setNames(as_tibble(predict(modi, nec0i, se.fit = TRUE)[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i <- mutate(nec0i,
                prev  = ilink(fit_link),
                upr = ilink(fit_link + (2 * se_link)),
                lwr = ilink(fit_link - (2 * se_link)))



ggplot(data=nec0i, aes(x=as.factor(year), y=prev,fill=as.factor(month)))+
  geom_bar(stat="identity",alpha=0.5, position="dodge")+
  theme_prism()+
  ylab("Proportion infected mosquitos")+
  xlab("Cohort year")+
  scale_fill_viridis_d(option = "G", labels=c("Jun","July","Aug","Sep","Oct","Nov","Dec"))+
  scale_color_viridis_d(option = "G", labels=c("Jun","July","Aug","Sep","Oct","Nov","Dec"))+
  theme(legend.position = "top")+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,1,by=0.05), limits=c(0,1))+
  geom_errorbar(aes(ymin=lwr, ymax=upr,col=as.factor(month)), width=0.1, linewidth=1.2, position=position_dodge(0.9))+
  coord_cartesian(ylim=c(0,0.45))

ggsave("OUTPUT/Fig_7.png",device="png", width=12, height=8, dpi=600)



#intensity of infection amongst all dissected mosquitos
ggplot(data=nec0,aes(x=as.factor(year),y=mean_ooc_all,fill=as.factor(year)))+
  geom_violin(alpha=0.5)+
  geom_jitter(alpha=0.1, width=0.2)+
  scale_y_log10()+
  theme_prism()+
  ylab("Oocysts/dissected mosquito")+
  xlab("Cohort year")+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_8.png",device="png", width=12, height=8, dpi=600)



#intensity of infections amongst infected mosquitos

ggplot(data=nec0,aes(x=as.factor(year),y=mean_ooc_pos,fill=as.factor(year)))+
  geom_violin(alpha=0.5)+
  geom_jitter(alpha=0.1, width=0.2)+
  scale_y_log10()+
  theme_prism()+
  ylab("Oocysts/infected mosquito")+
  xlab("Cohort year")+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_9.png",device="png", width=12, height=8, dpi=600)


#relationship between oocysts prevalence and intensity
#baseline infectivity vs gametocyte density (microscopy) 

modi = gam(mosq_pos_prop/100 ~ I(log(mean_ooc_all+0.1))+s(year, bs="re"), family=binomial(), data=nec0)
summary(modi)

nec0i = expand.grid(mean_ooc_all=seq(0,max(nec0$mean_ooc_all,na.rm=T), length.out=10000), year=nec0$year[1])
ilink <- family(modi)$linkinv
## add fit and se.fit on the **link** scale
nec0i <- bind_cols(nec0i, setNames(as_tibble(predict(modi, nec0i, se.fit = TRUE,exclude = c("s(year)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i <- mutate(nec0i,
                prev  = ilink(fit_link),
                upr = ilink(fit_link + (2 * se_link)),
                lwr = ilink(fit_link - (2 * se_link)))


modi2 = gam(mosq_pos_prop/100 ~ I(log(mean_ooc_all+0.1))*year, family=binomial(), data=nec0)
summary(modi2)

nec0i2 = expand.grid(mean_ooc_all=seq(0,max(nec0$mean_ooc_all,na.rm=T), length.out=10000),year=sort(unique((nec0$year[nec0$year %in% c("2019","2020","2021","2022")]))))
ilink <- family(modi2)$linkinv
## add fit and se.fit on the **link** scale
nec0i2 <- bind_cols(nec0i2, setNames(as_tibble(predict(modi2, nec0i2, se.fit = TRUE)[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i2 <- mutate(nec0i2,
                 prev  = ilink(fit_link),
                 upr = ilink(fit_link + (2 * se_link)),
                 lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=nec0i, aes(x=mean_ooc_all, y=prev))+
  geom_point(data=nec0[nec0$year %in% c("2019","2020","2021","2022"),],aes(x=mean_ooc_all,y=mosq_pos_prop/100,color=as.factor(year),alpha=0.5))+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2)+
  geom_line(data=nec0i2, aes(x=mean_ooc_all, y=prev, col=year),linewidth=1.05,alpha=0.7)+
  geom_line(linewidth=1.2)+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  xlab("Oocyst density")+
  ylab("Oocyst prevalence")+
  scale_x_log10()+
  scale_color_viridis_d()+
  theme(legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))+
  scale_y_continuous(labels = scales::percent)



ggsave("OUTPUT/Fig_9b.png",device="png", width=12, height=8, dpi=600)


# Baseline figures by arm ------------------------------------------------



library(ggplot2)
library(ggprism)
library(ggpubr)
library(epitools)

nec0s = nec0 %>% 
  group_by(arm) %>% 
  summarise(n=n(),count=sum(parpos),prev = mean(parpos)) %>% 
  mutate(lwr = binom.exact(count, n, conf.level = 0.95)[[4]], upr=binom.exact(count, n, conf.level = 0.95)[[5]])


ggplot(data=nec0s, aes(x=as.factor(arm), y=prev,fill=as.factor(arm)))+
  geom_bar(stat="identity",alpha=0.5)+
  theme_prism()+
  ylab("Asexual parasite prevalence")+
  theme(axis.title.x = element_blank())+
  scale_fill_viridis_d()+
  theme(legend.position = "none")+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,1,by=0.2), limits=c(0,1))+
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.1, linewidth=1.2)

ggsave("OUTPUT/Fig_1_arm.png",device="png", width=12, height=8, dpi=600)

options(scipen=999)
ggplot(data=nec0,aes(x=as.factor(arm),y=mic_pf_asc_ul,fill=as.factor(arm)))+
  geom_violin(alpha=0.5)+
  geom_jitter(alpha=0.1, width=0.2)+
  scale_y_log10()+
  theme_prism()+
  ylab("Asexual parasite density/µL")+
  xlab("Cohort arm")+
  scale_fill_viridis_d()+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_2_arm.png",device="png", width=12, height=8, dpi=600)


#gametocyte positives 
nec0g = nec0 %>% 
  group_by(arm) %>% 
  summarise(n=sum(!is.na(gampos)),count=sum(gampos,na.rm=TRUE),prev = mean(gampos,na.rm=TRUE)) %>% 
  mutate(lwr = binom.exact(count, n, conf.level = 0.95)[[4]], upr=binom.exact(count, n, conf.level = 0.95)[[5]])


ggplot(data=nec0g, aes(x=as.factor(arm), y=prev,fill=as.factor(arm)))+
  geom_bar(stat="identity",alpha=0.5)+
  theme_prism()+
  ylab("Gametocyte prevalence")+
  xlab("Cohort arm")+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,1,by=0.2), limits=c(0,1))+
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.1, linewidth=1.2)+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_3_arm.png",device="png", width=12, height=8, dpi=600)


options(scipen=999)
ggplot(data=nec0,aes(x=as.factor(arm),y=mic_pf_gct_ul,fill=as.factor(arm)))+
  geom_violin(alpha=0.5)+
  geom_jitter(alpha=0.1, width=0.2)+
  scale_y_log10()+
  theme_prism()+
  ylab("Gametocyte density mic/µL")+
  xlab("Cohort arm")+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_4_arm.png",device="png", width=12, height=8, dpi=600)


ggplot(data=nec0,aes(x=as.factor(arm),y=qrtpcr_gct_ul,fill=as.factor(arm)))+
  geom_violin(alpha=0.5)+
  geom_jitter(alpha=0.1, width=0.2)+
  scale_y_log10()+
  theme_prism()+
  ylab("Gametocyte density pcr/µL")+
  xlab("Cohort arm")+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_5_arm.png",device="png", width=12, height=8, dpi=600)


#mean infectivity over the arms



modi = gam(mosq_pos_prop/100 ~ arm, family=binomial(), data=nec0)
summary(modi)


nec0i = expand.grid(arm=sort(unique((nec0$arm))))
ilink <- family(modi)$linkinv
## add fit and se.fit on the **link** scale
nec0i <- bind_cols(nec0i, setNames(as_tibble(predict(modi, nec0i, se.fit = TRUE)[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i <- mutate(nec0i,
                prev  = ilink(fit_link),
                upr = ilink(fit_link + (2 * se_link)),
                lwr = ilink(fit_link - (2 * se_link)))


ggplot(data=nec0i, aes(x=as.factor(arm), y=prev,fill=as.factor(arm)))+
  geom_bar(stat="identity",alpha=0.5)+
  theme_prism()+
  ylab("Proportion infected mosquitos")+
  xlab("Cohort arm")+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,0.3,by=0.05), limits=c(0,0.35))+
  geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.1, linewidth=1.2)+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_6_arm.png",device="png", width=12, height=8, dpi=600)



#baseline infectivity vs gametocyte density (pcr) 

modi = gam(mosq_pos_prop/100 ~ I(log(gam.p))+s(study, bs="re"), family=binomial(), data=nec0)
summary(modi)

nec0i = expand.grid(gam.p=seq(0.01,max(nec0$gam.p,na.rm=T), length.out=10000), study=nec0$study[1])
ilink <- family(modi)$linkinv
## add fit and se.fit on the **link** scale
nec0i <- bind_cols(nec0i, setNames(as_tibble(predict(modi, nec0i, se.fit = TRUE,exclude = c("s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i <- mutate(nec0i,
                prev  = ilink(fit_link),
                upr = ilink(fit_link + (2 * se_link)),
                lwr = ilink(fit_link - (2 * se_link)))


modi2 = gam(mosq_pos_prop/100 ~ I(log(gam.p))*arm +s(study, bs="re"), family=binomial(), data=nec0)
summary(modi2)

nec0i2 = expand.grid(gam.p=seq(0.01,max(nec0$gam.p,na.rm=T), length.out=10000),arm=sort(unique((nec0$arm))), study = nec0$study[1])
ilink <- family(modi2)$linkinv
## add fit and se.fit on the **link** scale
nec0i2 <- bind_cols(nec0i2, setNames(as_tibble(predict(modi2, nec0i2, se.fit = TRUE, exclude=c("s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i2 <- mutate(nec0i2,
                 prev  = ilink(fit_link),
                 upr = ilink(fit_link + (2 * se_link)),
                 lwr = ilink(fit_link - (2 * se_link)))


ggplot(data=nec0i, aes(x=gam.p, y=prev))+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2)+
  geom_line(data=nec0i2, aes(x=gam.p, y=prev, col=arm),linewidth=1.05,alpha=0.7)+
  geom_line(linewidth=1.2)+
  theme_prism()+
  ylab("Proportion infected mosquitos")+
  xlab("Gametocyte density/µL (pcr)")+
  scale_x_log10()+
  scale_color_viridis_d()+
  theme(legend.position = "top")+
  scale_y_continuous(labels = scales::percent)

ggsave("OUTPUT/Fig_6b_arm.png",device="png", width=12, height=8, dpi=600)



#baseline infectivity vs gametocyte density (microscopy) 

modi = gam(mosq_pos_prop/100 ~ I(log(gam.m+0.1))+s(study, bs="re"), family=binomial(), data=nec0)
summary(modi)

nec0i = expand.grid(gam.m=seq(0,max(nec0$gam.m,na.rm=T), length.out=10000), study=nec0$study[1])
ilink <- family(modi)$linkinv
## add fit and se.fit on the **link** scale
nec0i <- bind_cols(nec0i, setNames(as_tibble(predict(modi, nec0i, se.fit = TRUE,exclude = c("s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i <- mutate(nec0i,
                prev  = ilink(fit_link),
                upr = ilink(fit_link + (2 * se_link)),
                lwr = ilink(fit_link - (2 * se_link)))


modi2 = gam(mosq_pos_prop/100 ~ I(log(gam.m+0.1))*arm +s(study, bs="re"), family=binomial(), data=nec0)
summary(modi2)

nec0i2 = expand.grid(gam.m=seq(0,max(nec0$gam.m,na.rm=T), length.out=10000),arm=sort(unique((nec0$arm))), study=nec0$study[1])
ilink <- family(modi2)$linkinv
## add fit and se.fit on the **link** scale
nec0i2 <- bind_cols(nec0i2, setNames(as_tibble(predict(modi2, nec0i2, se.fit = TRUE, exclude=c("s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i2 <- mutate(nec0i2,
                 prev  = ilink(fit_link),
                 upr = ilink(fit_link + (2 * se_link)),
                 lwr = ilink(fit_link - (2 * se_link)))


ggplot(data=nec0i, aes(x=gam.m, y=prev))+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2)+
  geom_line(data=nec0i2, aes(x=gam.m, y=prev, col=arm),linewidth=1.05,alpha=0.7)+
  geom_line(linewidth=1.2)+
  theme_prism()+
  ylab("Proportion infected mosquitos")+
  xlab("Gametocyte density/µL (microscopy)")+
  scale_x_log10()+
  scale_color_viridis_d()+
  theme(legend.position = "top")+
  scale_y_continuous(labels = scales::percent)

ggsave("OUTPUT/Fig_6c_arm.png",device="png", width=12, height=8, dpi=600)


#correspondence between gametocyte densities by microscopy and by PCR

modi = gam(I(log10(gam.p+0.01)) ~ I(log10(gam.m+0.1))+s(study, bs="re"), data=nec0)
summary(modi)

nec0i = expand.grid(gam.m=seq(min(nec0$gam.m[nec0$gam.m>0],na.rm=T),max(nec0$gam.m,na.rm=T), length.out=10000), study=nec0$study[1])
ilink <- family(modi)$linkinv
## add fit and se.fit on the **link** scale
nec0i <- bind_cols(nec0i, setNames(as_tibble(predict(modi, nec0i, se.fit = TRUE,exclude = c("s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i <- mutate(nec0i,
                gam.p  = 10^(fit_link),
                upr = 10^(fit_link + (2 * se_link)),
                lwr = 10^(fit_link - (2 * se_link)))

modi2 = gam(I(log10(gam.p+0.01)) ~ I(log10(gam.m+0.1))*arm+s(study, bs="re"), data=nec0)
summary(modi2)

nec0i2 = expand.grid(gam.m=seq(min(nec0$gam.m[nec0$gam.m>0],na.rm=T),max(nec0$gam.m,na.rm=T), length.out=10000), study=nec0$study[1], arm=unique(sort(nec0$arm)))
ilink <- family(modi2)$linkinv
## add fit and se.fit on the **link** scale
nec0i2 <- bind_cols(nec0i2, setNames(as_tibble(predict(modi2, nec0i2, se.fit = TRUE,exclude = c("s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i2 <- mutate(nec0i2,
                gam.p  = 10^(fit_link),
                upr = 10^(fit_link + (2 * se_link)),
                lwr = 10^(fit_link - (2 * se_link)))


ggplot(data=nec0i, aes(x=gam.m, y=gam.p))+
  geom_point(data=nec0, aes(x=gam.m, y=gam.p, col=arm),alpha=0.4)+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2)+
  geom_line(data=nec0i2, aes(x=gam.m,y=gam.p,col=arm),linewidth=1.05, alpha=0.7)+
  geom_line(linewidth=1.2)+
  theme_prism()+
  ylab("Gametocyte density/µL (pcr)")+
  xlab("Gametocyte density/µL (microscopy)")+
  scale_x_log10()+
  scale_y_log10()+
  scale_color_viridis_d()+
  theme(legend.position = "top")

ggsave("OUTPUT/Fig_6d_arm.png",device="png", width=12, height=8, dpi=600)




#intensity of infection amongst all dissected mosquitos
ggplot(data=nec0,aes(x=as.factor(arm),y=mean_ooc_all,fill=as.factor(arm)))+
  geom_violin(alpha=0.5)+
  geom_jitter(alpha=0.1, width=0.2)+
  scale_y_log10()+
  theme_prism()+
  ylab("Oocysts/dissected mosquito")+
  xlab("Cohort arm")+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_8_arm.png",device="png", width=12, height=8, dpi=600)



#intensity of infections amongst infected mosquitos

ggplot(data=nec0,aes(x=as.factor(arm),y=mean_ooc_pos,fill=as.factor(arm)))+
  geom_violin(alpha=0.5)+
  geom_jitter(alpha=0.1, width=0.2)+
  scale_y_log10()+
  theme_prism()+
  ylab("Oocysts/infected mosquito")+
  xlab("Cohort arm")+
  scale_fill_viridis_d()+
  theme(legend.position = "top")+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_9_arm.png",device="png", width=12, height=8, dpi=600)






# Asexual parasite density (mic) over time -----------------------------

####Overall####


library(mgcv)

mod1a = gam(I(asx+0.1) ~ s(days)+s(studyid_str, bs="re"), family=Gamma(link="log"), data=nec)
summary(mod1a)
##plot(mod1a)

newd = expand.grid(days = seq(0,max(nec$days),by=0.1), study = nec$study[1], studyid_str = nec$studyid_str[1])
newd$prop=predict(mod1a,newd, type="response", exclude=c("s(studyid_str)"))
fam <- family(mod1a)
ilink <- family(mod1a)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod1a, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),
                                  c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=days, y= fit_resp))+
  geom_line(col="darkblue", linewidth=1.2)+
  geom_ribbon(aes(ymin=right_lwr, ymax=right_upr), fill="darkblue", alpha=0.3)+
  theme_prism()+
  coord_cartesian(ylim=c(0,600),xlim=c(0,5))+
  ylab("Asexual parasite density/µL")+
  xlab("Days after treatment")


ggsave("OUTPUT/Fig_010.png",device="png", width=12, height=8, dpi=600)


####by cohort prev on day 2####
#mod1b = gam(I(asx>0) ~ study, family=binomial(), data=nec[nec$days==2,])
#summary(mod1b)

#newd = expand.grid(study = sort(unique((nec[nec$days==2,]$study))), studyid_str = nec$studyid_str[1])
#newd$prop=predict(mod1b,newd, type="response", exclude=c("s(studyid_str)"))
#ilink <- family(mod1b)$linkinv
## add fit and se.fit on the **link** scale
#ndata <- bind_cols(newd, setNames(as_tibble(predict(mod1b, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
#ndata <- mutate(ndata,
#fit_resp  = ilink(fit_link),
#right_upr = ilink(fit_link + (2 * se_link)),
#right_lwr = ilink(fit_link - (2 * se_link)))

#ggplot(data=ndata, aes(x=study, y= fit_resp))+
#geom_bar(aes(fill=study),linewidth=1.2,stat="identity",position = position_dodge())+
#geom_errorbar(aes(ymin=right_lwr, ymax=right_upr, group=study),width=0.25, linewidth=0.7)+
#theme_prism()+
#ylab("Asexual parasite prevalence (day 2)")+
#xlab("Cohort")+
#scale_y_continuous(labels=scales::percent)


#ggsave("OUTPUT/Fig_011.png",device="png", width=12, height=8, dpi=600)


####by ARM prev on day 2####
#mod1b = gam(I(asx>0) ~ arm, family=binomial(), data=nec[nec$days==2,])
#summary(mod1b)

#newd = expand.grid(arm = sort(unique((nec[nec$days==2,]$arm))), studyid_str = nec$studyid_str[1])
#newd$prop=predict(mod1b,newd, type="response", exclude=c("s(studyid_str)"))
#ilink <- family(mod1b)$linkinv
## add fit and se.fit on the **link** scale
#ndata <- bind_cols(newd, setNames(as_tibble(predict(mod1b, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
#ndata <- mutate(ndata,
#fit_resp  = ilink(fit_link),
#right_upr = ilink(fit_link + (2 * se_link)),
#right_lwr = ilink(fit_link - (2 * se_link)))

#ggplot(data=ndata, aes(x=arm, y= fit_resp))+
#geom_bar(aes(fill=arm),linewidth=1.2,stat="identity",position = position_dodge())+
#geom_errorbar(aes(ymin=right_lwr, ymax=right_upr, group=arm),width=0.25, linewidth=0.7)+
#theme_prism()+
#ylab("Asexual parasite prevalence (day 2)")+
#xlab("Cohort")+
#scale_y_continuous(labels=scales::percent)+
#theme(axis.title.x = element_blank())+
#theme(legend.position = "none")



#ggsave("OUTPUT/Fig_011_arm.png",device="png", width=12, height=8, dpi=600)



# Gametocytes over time ---------------------------------------------------


####microscopy overall####
mod2a = gam(I(gam.m+0.01) ~  s(days)+s(studyid_str, bs="re"), family=Gamma(link="log"), data=nec)
summary(mod2a)

newd = expand.grid(days = seq(0,max(nec$days),by=0.1), study = nec$study[1], studyid_str = nec$studyid_str[1])
newd$prop=predict(mod2a,newd, type="response", exclude=c("s(studyid_str)"))
ilink <- family(mod2a)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod2a, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),
                                  c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=days, y= fit_resp))+
  geom_line(col="darkgreen", linewidth=1.2)+
  geom_ribbon(aes(ymin=right_lwr, ymax=right_upr), fill="darkgreen", alpha=0.3)+
  theme_prism()+
  coord_cartesian(ylim=c(0,150))+
  ylab("Gametocyte density/µL (microscopy)")+
  xlab("Days after treatment")

ggsave("OUTPUT/Fig_012.png",device="png", width=12, height=8, dpi=600)



####microscopy over cohort####

mod2b = gam(I(gam.m+0.1) ~  study + s(days, k=5) + s(days, by = as.ordered(study), k=5)+s(studyid_str, bs="re"), family=Gamma(link="log"), data=nec[!is.na(nec$gam.m),])
summary(mod2b)


newd = expand.grid(days = seq(0,max(nec$days),by=0.1), study = sort(unique((nec[!is.na(nec$gam.m),]$study))), studyid_str = nec[!is.na(nec$gam.m),]$studyid_str[1])
newd$prop=predict(mod2b,newd, type="response", exclude=c("s(studyid_str)"))
ilink <- family(mod2b)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod2b, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),
                                  c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=days, y= fit_resp))+
  geom_line(aes(col=study), linewidth=1.2)+
  geom_ribbon(aes(ymin=right_lwr, ymax=right_upr, fill=study), alpha=0.25)+
  theme_prism()+
  coord_cartesian(ylim=c(0,150), xlim=c(0,30))+
  ylab("Gametocyte density/µL (microscopy)")+
  xlab("Days after treatment")


ggsave("OUTPUT/Fig_013.png",device="png", width=12, height=8, dpi=600)


####microscopy over arm####

mod2b = gam(I(gam.m+0.1) ~  arm + s(days, k=5) + s(days, by = as.ordered(arm), k=5)+s(studyid_str, bs="re"), family=Gamma(link="log"), data=nec[!is.na(nec$gam.m),])
summary(mod2b)


newd = expand.grid(days = seq(0,max(nec$days),by=0.1), arm = sort(unique((nec[!is.na(nec$gam.m),]$arm))), studyid_str = nec[!is.na(nec$gam.m),]$studyid_str[1])
newd$prop=predict(mod2b,newd, type="response", exclude=c("s(studyid_str)"))
ilink <- family(mod2b)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod2b, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),
                                  c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=days, y= fit_resp))+
  geom_line(aes(col=arm), linewidth=1.2)+
  geom_ribbon(aes(ymin=right_lwr, ymax=right_upr, fill=arm), alpha=0.25)+
  theme_prism()+
  coord_cartesian(ylim=c(0,150), xlim=c(0,30))+
  ylab("Gametocyte density/µL (microscopy)")+
  xlab("Days after treatment")


ggsave("OUTPUT/Fig_013_arm.png",device="png", width=12, height=8, dpi=600)




####Microscopy prevalence####

mod2c = gam(I(gam.m>0) ~  study*I(log10(days+1))+s(studyid_str, bs="re"), family=binomial(), data=nec[!is.na(nec$gam.m),])
summary(mod2c)


newd = expand.grid(days = seq(0,max(nec$days),by=0.1), study = sort(unique((nec[!is.na(nec$gam.m),]$study))), studyid_str = nec[!is.na(nec$gam.m),]$studyid_str[1])
newd$prop=predict(mod2c,newd, type="response", exclude=c("s(studyid_str)"))
ilink <- family(mod2c)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod2c, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),
                                  c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=days, y= fit_resp))+
  geom_line(aes(col=study), linewidth=1.2, alpha=0.9)+
  geom_ribbon(aes(ymin=right_lwr, ymax=right_upr, fill=study), alpha=0.2)+
  theme_prism()+
  #coord_cartesian(ylim=c(0,150), xlim=c(0,30))+
  ylab("Gametocyte (mic) prevalence")+
  xlab("Days after treatment")+
  scale_y_continuous(labels=scales::percent)

ggsave("OUTPUT/Fig_014.png",device="png", width=12, height=8, dpi=600)

##by ARM##
mod2c = gam(I(gam.m>0) ~  arm*I(log10(days+1))+s(studyid_str, bs="re"), family=binomial(), data=nec[!is.na(nec$gam.m),])
summary(mod2c)


newd = expand.grid(days = seq(0,max(nec$days),by=0.1), arm = sort(unique((nec[!is.na(nec$gam.m),]$arm))), studyid_str = nec[!is.na(nec$gam.m),]$studyid_str[1])
newd$prop=predict(mod2c,newd, type="response", exclude=c("s(studyid_str)"))
ilink <- family(mod2c)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod2c, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),
                                  c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=days, y= fit_resp))+
  geom_line(aes(col=arm), linewidth=1.2, alpha=0.9)+
  geom_ribbon(aes(ymin=right_lwr, ymax=right_upr, fill=arm), alpha=0.2)+
  theme_prism()+
  #coord_cartesian(ylim=c(0,150), xlim=c(0,30))+
  ylab("Gametocyte (mic) prevalence")+
  xlab("Days after treatment")+
  scale_y_continuous(labels=scales::percent)

ggsave("OUTPUT/Fig_014_arm.png",device="png", width=12, height=8, dpi=600)


####pcr overall####

mod3a = gam(I(gam.p+0.01) ~  s(days)+s(studyid_str, bs="re"), family=Gamma(link="log"), data=nec[!is.na(nec$gam.p),])
summary(mod3a)


newd = expand.grid(days = seq(0,max(nec$days),by=0.1),  studyid_str = nec[!is.na(nec$gam.p),]$studyid_str[1])

newd$prop=predict(mod3a,newd, type="response", exclude=c("s(studyid_str)"))
ilink <- family(mod3a)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod3a, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=days, y= fit_resp))+
  geom_line(col="darkgreen", linewidth=1.2)+
  geom_ribbon(aes(ymin=right_lwr, ymax=right_upr), fill="darkgreen", alpha=0.3)+
  theme_prism()+
  coord_cartesian(ylim=c(0,150))+
  ylab("Gametocyte density/µL (pcr)")+
  xlab("Days after treatment")

ggsave("OUTPUT/Fig_015.png",device="png", width=12, height=8, dpi=600)


####pcr by cohort####

mod3b = gam(I(gam.p+0.01) ~  study + s(days, k=5) + s(days, by = as.ordered(study), k=5)+s(studyid_str, bs="re"), family=Gamma(link="log"), data=nec[!is.na(nec$gam.p),])
summary(mod3b)


newd = expand.grid(days = seq(0,max(nec$days),by=0.1), study = sort(unique((nec[!is.na(nec$gam.p),]$study))), studyid_str = nec[!is.na(nec$gam.p),]$studyid_str[1])
newd$prop=predict(mod3b,newd, type="response", exclude=c("s(study)", "s(studyid_str)"))
ilink <- family(mod3b)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod3b, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),
                                  c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=days, y= fit_resp))+
  geom_line(aes(col=study), linewidth=1.2)+
  geom_ribbon(aes(ymin=right_lwr, ymax=right_upr, fill=study), alpha=0.1)+
  theme_prism()+
  coord_cartesian(ylim=c(0,150), xlim=c(0,30))+
  ylab("Gametocyte density/µL (pcr)")+
  xlab("Days after treatment")

ggsave("OUTPUT/Fig_016.png",device="png", width=12, height=8, dpi=600)


####pcr by arm####

mod3b = gam(I(gam.p+0.01) ~  arm + s(days, k=5) + s(days, by = as.ordered(arm), k=5)+s(studyid_str, bs="re"), family=Gamma(link="log"), data=nec[!is.na(nec$gam.p),])
summary(mod3b)


newd = expand.grid(days = seq(0,max(nec$days),by=0.1), arm = sort(unique((nec[!is.na(nec$gam.p),]$arm))), studyid_str = nec[!is.na(nec$gam.p),]$studyid_str[1])
newd$prop=predict(mod3b,newd, type="response", exclude=c("s(arm)", "s(studyid_str)"))
ilink <- family(mod3b)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod3b, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),
                                  c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=days, y= fit_resp))+
  geom_line(aes(col=arm), linewidth=1.2)+
  geom_ribbon(aes(ymin=right_lwr, ymax=right_upr, fill=arm), alpha=0.1)+
  theme_prism()+
  coord_cartesian(ylim=c(0,150), xlim=c(0,30))+
  ylab("Gametocyte density/µL (pcr)")+
  xlab("Days after treatment")

ggsave("OUTPUT/Fig_016_arm.png",device="png", width=12, height=8, dpi=600)



####pcr prevalence####

mod3c = gam(I(gam.p>0) ~  study*I(log10(days+1))+s(studyid_str, bs="re"), family=binomial(), data=nec[!is.na(nec$gam.p),])
summary(mod3c)


newd = expand.grid(days = seq(0,max(nec$days),by=0.1), study = sort(unique((nec[!is.na(nec$gam.p),]$study))), studyid_str = nec[!is.na(nec$gam.p),]$studyid_str[1])
newd$prop=predict(mod3c,newd, type="response", exclude=c("s(studyid_str)"))
ilink <- family(mod3c)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod3c, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),
                                  c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=days, y= fit_resp))+
  geom_line(aes(col=study), linewidth=1.2, alpha=0.9)+
  geom_ribbon(aes(ymin=right_lwr, ymax=right_upr, fill=study), alpha=0.2)+
  theme_prism()+
  #coord_cartesian(ylim=c(0,150), xlim=c(0,30))+
  ylab("Gametocyte (pcr) prevalence")+
  xlab("Days after treatment")+
  scale_y_continuous(labels=scales::percent)


ggsave("OUTPUT/Fig_017.png",device="png", width=12, height=8, dpi=600)


##by ARM##

mod3c = gam(I(gam.p>0) ~  arm*I(log10(days+1))+s(studyid_str, bs="re"), family=binomial(), data=nec[!is.na(nec$gam.p),])
summary(mod3c)


newd = expand.grid(days = seq(0,max(nec$days),by=0.1), arm = sort(unique((nec[!is.na(nec$gam.p),]$arm))), studyid_str = nec[!is.na(nec$gam.p),]$studyid_str[1])
newd$prop=predict(mod3c,newd, type="response", exclude=c("s(studyid_str)"))
ilink <- family(mod3c)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod3c, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),
                                  c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=days, y= fit_resp))+
  geom_line(aes(col=arm), linewidth=1.2, alpha=0.9)+
  geom_ribbon(aes(ymin=right_lwr, ymax=right_upr, fill=arm), alpha=0.2)+
  theme_prism()+
  #coord_cartesian(ylim=c(0,150), xlim=c(0,30))+
  ylab("Gametocyte (pcr) prevalence")+
  xlab("Days after treatment")+
  scale_y_continuous(labels=scales::percent)


ggsave("OUTPUT/Fig_017_arm.png",device="png", width=12, height=8, dpi=600)



# Proportion infected mosquitos over time ---------------------------------

####overall####

necfeed = nec[!is.na(nec$mosq_total) & !(nec$mosq_total==0),]
necfeed$prop = necfeed$mosq_pos/necfeed$mosq_total
necfeed$prop = ifelse(necfeed$prop==0,0.0001,ifelse(necfeed$prop==1,0.9999,necfeed$prop))

mod4a = gam(prop ~ log10(days+1)+s(studyid_str, bs="re"), family=binomial(), weights = mosq_total, data=necfeed)
summary(mod4a)
#plot(mod4a)


newd = expand.grid(days = seq(0,max(nec$days),by=0.1), study = nec$study[1], studyid_str = nec$studyid_str[1])
newd$prop=predict(mod4a,newd, type="response", exclude=c("s(studyid_str)"))
ilink <- family(mod4a)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod4a, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),
                                  c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=days, y= fit_resp))+
  geom_line(col="darkred", linewidth=1.2)+
  geom_ribbon(aes(ymin=right_upr, ymax=right_lwr), fill="darkred", alpha=0.3)+
  theme_prism()+
  ylab("Proportion infected mosquitos")+
  xlab("Days after treatment")+
  scale_y_continuous(labels=scales::percent)+
  coord_cartesian(ylim=c(0,0.05))


ggsave("OUTPUT/Fig_018.png",device="png", width=12, height=8, dpi=600)



mod4b = gam(prop ~ study*I(log10(days+1))+s(studyid_str, bs="re"), family=binomial(), weights = mosq_total, necfeed[!is.na(necfeed$prop),])
summary(mod4b)


newd = expand.grid(days = seq(0,max(necfeed$days),by=0.1), study = sort(unique((necfeed[!is.na(necfeed$prop),]$study))), studyid_str = necfeed[!is.na(necfeed$prop),]$studyid_str[1])
newd$prop=predict(mod4b,newd, type="response", exclude=c("s(studyid_str)"))
ilink <- family(mod4b)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod4b, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=days, y= fit_resp))+
  geom_line(aes(col=study), linewidth=1.2)+
  geom_ribbon(aes(ymin=right_upr, ymax=right_lwr, fill=study), alpha=0.1)+
  theme_prism()+
  ylab("Proportion infected mosquitos")+
  xlab("Days after treatment")+
  scale_y_continuous(labels=scales::percent)+
  coord_cartesian(ylim=c(0,0.05))


ggsave("OUTPUT/Fig_019.png",device="png", width=12, height=8, dpi=600)



##by ARM##

mod4b = gam(prop ~ arm*I(log10(days+1))+s(studyid_str, bs="re"), family=binomial(), weights = mosq_total, necfeed[!is.na(necfeed$prop),])
summary(mod4b)


newd = expand.grid(days = seq(0,max(necfeed$days),by=0.1), arm = sort(unique((necfeed[!is.na(necfeed$prop),]$arm))), studyid_str = necfeed[!is.na(necfeed$prop),]$studyid_str[1])
newd$prop=predict(mod4b,newd, type="response", exclude=c("s(studyid_str)"))
ilink <- family(mod4b)$linkinv
## add fit and se.fit on the **link** scale
ndata <- bind_cols(newd, setNames(as_tibble(predict(mod4b, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
ndata <- mutate(ndata,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=ndata, aes(x=days, y= fit_resp))+
  geom_line(aes(col=arm), linewidth=1.2)+
  geom_ribbon(aes(ymin=right_upr, ymax=right_lwr, fill=arm), alpha=0.1)+
  theme_prism()+
  ylab("Proportion infected mosquitos")+
  xlab("Days after treatment")+
  scale_y_continuous(labels=scales::percent)+
  coord_cartesian(ylim=c(0,0.05))


ggsave("OUTPUT/Fig_019_arm.png",device="png", width=12, height=8, dpi=600)




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


library(survival)
library(survminer)
library(pammtools)

km5a = survfit(Surv(days,status)~1, data= nectime5a)
survout5a = survminer::ggsurvplot(km5a)
survout5a$plot+
  xlab("Time to clearance of gametocytes (pcr)")+
  scale_y_continuous(labels=scales::percent)+
  theme_prism()+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_020.png",device="png", width=12, height=8, dpi=600)


nectime5a = nectime5a %>% left_join(
  necgam %>% 
    select(studyid_str, study, gam0.m, gam0.p) %>%
    distinct(), by=c("study","studyid_str")
)

mod5a = coxph(Surv(days,status)~arm+I(log(gam0.p+0.01)) + cluster(study), data= nectime5a)
sumdat=summary(mod5a)

ref = levels(nectime5a$arm)

newd0 = nectime5a %>% make_newdata(
  arm=levels(arm)
) %>% select(arm) %>%
  mutate(gmean="")

current = ref[1]
mod5a = coxph(Surv(days,status)~relevel(arm, ref=current) + I(log10(gam0.p+0.01)) + cluster(study), data= nectime5a)
sumdat = summary(mod5a)
newd = data.frame(arm = levels(relevel(nectime5a$arm, ref=current))) %>% 
  mutate(hazard = c(1,sumdat[["conf.int"]][1:(length(ref)-1),1]),
         lwr = c(1,sumdat[["conf.int"]][1:(length(ref)-1),3]),
         upr =c(1,sumdat[["conf.int"]][1:(length(ref)-1),4]),
         pval = c(1,sumdat[["coefficients"]][1:(length(ref)-1),6]),
         reference = current,result= ifelse(arm==current, newd0$gmean[newd0$arm==current][[1]]  ,paste0(ifelse(hazard>999,Inf,format(round(hazard,2), nsmall = 2)), " (",format(round(lwr,2), nsmall = 2),"-",ifelse(upr>999,Inf,format(round(upr,2), nsmall = 2)),")\n",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))) %>%
  dplyr::select(arm, reference, result)


for (i in 2:length(ref)){
  current = ref[i]
  mod5a = coxph(Surv(days,status)~relevel(arm, ref=current) + I(log10(gam0.p+0.01)) + cluster(study), data= nectime5a)
  sumdat = summary(mod5a)
  newd1 = data.frame(arm = levels(relevel(nectime5a$arm, ref=current))) %>% 
    mutate(hazard = c(1,sumdat[["conf.int"]][1:(length(ref)-1),1]),
           lwr = c(1,sumdat[["conf.int"]][1:(length(ref)-1),3]),
           upr =c(1,sumdat[["conf.int"]][1:(length(ref)-1),4]),
           pval = c(1,sumdat[["coefficients"]][1:(length(ref)-1),6]),
           reference = current,result= ifelse(arm==current, newd0$gmean[newd0$arm==current][[1]]  ,paste0(ifelse(hazard>999,Inf,format(round(hazard,2), nsmall = 2)), " (",format(round(lwr,2), nsmall = 2),"-",ifelse(upr>999,Inf,format(round(upr,2), nsmall = 2)),")\n",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))) %>%
    dplyr::select(arm, reference, result)
  
  newd = rbind(newd, newd1)
  
  print(i)
  
}

newd = unique(as.data.frame(newd))

result_matrix5a <- pivot_wider(newd, names_from = arm, values_from = result)

flextable(result_matrix5a) %>% autofit()  %>% save_as_docx( path = "OUTPUT/survival_gampcr_.docx")



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


library(survival)
library(survminer)

km5b = survfit(Surv(days,status)~1, data= nectime5b)
survout5b = survminer::ggsurvplot(km5b)
survout5b$plot+
  xlab("Time to clearance of infectivity")+
  scale_y_continuous(labels=scales::percent)+
  theme_prism()+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_021.png",device="png", width=12, height=8, dpi=600)


nectime5b = nectime5b %>% left_join(
  necgam %>% 
    select(studyid_str, study, gam0.m, gam0.p) %>%
    distinct(), by=c("study","studyid_str")
)

mod5b = coxph(Surv(days,status)~arm+I(log(gam0.p+0.01)) + cluster(study), data= nectime5b)
sumdat=summary(mod5b)

ref = levels(nectime5b$arm)

newd0 = nectime5b %>% make_newdata(
  arm=levels(arm)
) %>% select(arm) %>%
  mutate(gmean="")

current = ref[1]
mod5b = coxph(Surv(days,status)~relevel(arm, ref=current) + I(log10(gam0.p+0.01)) + cluster(study), data= nectime5b)
sumdat = summary(mod5b)
newd = data.frame(arm = levels(relevel(nectime5b$arm, ref=current))) %>% 
  mutate(hazard = c(1,sumdat[["conf.int"]][1:(length(ref)-1),1]),
         lwr = c(1,sumdat[["conf.int"]][1:(length(ref)-1),3]),
         upr =c(1,sumdat[["conf.int"]][1:(length(ref)-1),4]),
         pval = c(1,sumdat[["coefficients"]][1:(length(ref)-1),6]),
         reference = current,result= ifelse(arm==current, newd0$gmean[newd0$arm==current][[1]]  ,paste0(ifelse(hazard>999,Inf,format(round(hazard,2), nsmall = 2)), " (",format(round(lwr,2), nsmall = 2),"-",ifelse(upr>999,Inf,format(round(upr,2), nsmall = 2)),")\n",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))) %>%
  dplyr::select(arm, reference, result)


for (i in 2:length(ref)){
  current = ref[i]
  mod5b = coxph(Surv(days,status)~relevel(arm, ref=current) + I(log10(gam0.p+0.01)) + cluster(study), data= nectime5b)
  sumdat = summary(mod5b)
  newd1 = data.frame(arm = levels(relevel(nectime5b$arm, ref=current))) %>% 
    mutate(hazard = c(1,sumdat[["conf.int"]][1:(length(ref)-1),1]),
           lwr = c(1,sumdat[["conf.int"]][1:(length(ref)-1),3]),
           upr =c(1,sumdat[["conf.int"]][1:(length(ref)-1),4]),
           pval = c(1,sumdat[["coefficients"]][1:(length(ref)-1),6]),
           reference = current,result= ifelse(arm==current, newd0$gmean[newd0$arm==current][[1]]  ,paste0(ifelse(hazard>999,Inf,format(round(hazard,2), nsmall = 2)), " (",format(round(lwr,2), nsmall = 2),"-",ifelse(upr>999,Inf,format(round(upr,2), nsmall = 2)),")\n",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))) %>%
    dplyr::select(arm, reference, result)
  
  newd = rbind(newd, newd1)
  
  print(i)
  
}

newd = unique(as.data.frame(newd))

result_matrix5b <- pivot_wider(newd, names_from = arm, values_from = result)

flextable(result_matrix5b) %>% autofit()  %>% save_as_docx( path = "OUTPUT/survival_infectivity_adjpcr_.docx")



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


library(survival)
library(survminer)

km5c = survfit(Surv(days,status)~1, data= nectime5c)
survout5c = survminer::ggsurvplot(km5c)
survout5c$plot+
  xlab("Time to clearance of gametocytes (microscopy)")+
  scale_y_continuous(labels=scales::percent)+
  theme_prism()+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_022.png",device="png", width=12, height=8, dpi=600)

nectime5c = nectime5c %>% left_join(
  necgam %>% 
    select(studyid_str, study, gam0.p, gam0.m) %>%
    distinct(), by=c("studyid_str","study")
)

mod5c = coxph(Surv(days,status)~arm+I(log(gam0.m+0.1)) + cluster(study), data= nectime5c)
sumdat=summary(mod5c)

ref = levels(nectime5c$arm)

newd0 = nectime5c %>% make_newdata(
  arm=levels(arm)
) %>% select(arm) %>%
  mutate(gmean="")

current = ref[1]
mod5c = coxph(Surv(days,status)~relevel(arm, ref=current) + I(log10(gam0.m+0.1)) + cluster(study), data= nectime5c)
sumdat = summary(mod5c)
newd = data.frame(arm = levels(relevel(nectime5c$arm, ref=current))) %>% 
  mutate(hazard = c(1,sumdat[["conf.int"]][1:(length(ref)-1),1]),
         lwr = c(1,sumdat[["conf.int"]][1:(length(ref)-1),3]),
         upr =c(1,sumdat[["conf.int"]][1:(length(ref)-1),4]),
         pval = c(1,sumdat[["coefficients"]][1:(length(ref)-1),6]),
         reference = current,result= ifelse(arm==current, newd0$gmean[newd0$arm==current][[1]]  ,paste0(ifelse(hazard>999,Inf,format(round(hazard,2), nsmall = 2)), " (",format(round(lwr,2), nsmall = 2),"-",ifelse(upr>999,Inf,format(round(upr,2), nsmall = 2)),")\n",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))) %>%
  dplyr::select(arm, reference, result)


for (i in 2:length(ref)){
  current = ref[i]
  mod5c = coxph(Surv(days,status)~relevel(arm, ref=current) + I(log10(gam0.m+0.1)) + cluster(study), data= nectime5c)
  sumdat = summary(mod5c)
  newd1 = data.frame(arm = levels(relevel(nectime5c$arm, ref=current))) %>% 
    mutate(hazard = c(1,sumdat[["conf.int"]][1:(length(ref)-1),1]),
           lwr = c(1,sumdat[["conf.int"]][1:(length(ref)-1),3]),
           upr =c(1,sumdat[["conf.int"]][1:(length(ref)-1),4]),
           pval = c(1,sumdat[["coefficients"]][1:(length(ref)-1),6]),
           reference = current,result= ifelse(arm==current, newd0$gmean[newd0$arm==current][[1]]  ,paste0(ifelse(hazard>999,Inf,format(round(hazard,2), nsmall = 2)), " (",format(round(lwr,2), nsmall = 2),"-",ifelse(upr>999,Inf,format(round(upr,2), nsmall = 2)),")\n",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))) %>%
    dplyr::select(arm, reference, result)
  
  newd = rbind(newd, newd1)
  
  print(i)
  
}

newd = unique(as.data.frame(newd))

result_matrix5c <- pivot_wider(newd, names_from = arm, values_from = result)

flextable(result_matrix5c) %>% autofit()  %>% save_as_docx( path = "OUTPUT/survival_gammic_.docx")


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


library(survival)
library(survminer)

km5d = survfit(Surv(days,status)~1, data= nectime5d)
survout5d = survminer::ggsurvplot(km5d)
survout5d$plot+
  xlab("Time to clearance of infectivity")+
  scale_y_continuous(labels=scales::percent)+
  theme_prism()+
  theme(legend.position = "none")

ggsave("OUTPUT/Fig_022.png",device="png", width=12, height=8, dpi=600)


nectime5d = nectime5d %>% left_join(
  necgam %>% 
    select(studyid_str, study, gam0.m, gam0.p) %>%
    distinct(), by=c("study","studyid_str")
)

mod5d = coxph(Surv(days,status)~arm+I(log(gam0.m+0.1)) + cluster(study), data= nectime5d)
sumdat=summary(mod5d)

ref = levels(nectime5d$arm)

newd0 = nectime5d %>% make_newdata(
  arm=levels(arm)
) %>% select(arm) %>%
  mutate(gmean="")

current = ref[1]
mod5d = coxph(Surv(days,status)~relevel(arm, ref=current) + I(log10(gam0.m+0.1)) + cluster(study), data= nectime5d)
sumdat = summary(mod5d)
newd = data.frame(arm = levels(relevel(nectime5d$arm, ref=current))) %>% 
  mutate(hazard = c(1,sumdat[["conf.int"]][1:(length(ref)-1),1]),
         lwr = c(1,sumdat[["conf.int"]][1:(length(ref)-1),3]),
         upr =c(1,sumdat[["conf.int"]][1:(length(ref)-1),4]),
         pval = c(1,sumdat[["coefficients"]][1:(length(ref)-1),6]),
         reference = current,result= ifelse(arm==current, newd0$gmean[newd0$arm==current][[1]]  ,paste0(ifelse(hazard>999,Inf,format(round(hazard,2), nsmall = 2)), " (",format(round(lwr,2), nsmall = 2),"-",ifelse(upr>999,Inf,format(round(upr,2), nsmall = 2)),")\n",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))) %>%
  dplyr::select(arm, reference, result)


for (i in 2:length(ref)){
  current = ref[i]
  mod5d = coxph(Surv(days,status)~relevel(arm, ref=current) + I(log10(gam0.m+0.1)) + cluster(study), data= nectime5d)
  sumdat = summary(mod5d)
  newd1 = data.frame(arm = levels(relevel(nectime5d$arm, ref=current))) %>% 
    mutate(hazard = c(1,sumdat[["conf.int"]][1:(length(ref)-1),1]),
           lwr = c(1,sumdat[["conf.int"]][1:(length(ref)-1),3]),
           upr =c(1,sumdat[["conf.int"]][1:(length(ref)-1),4]),
           pval = c(1,sumdat[["coefficients"]][1:(length(ref)-1),6]),
           reference = current,result= ifelse(arm==current, newd0$gmean[newd0$arm==current][[1]]  ,paste0(ifelse(hazard>999,Inf,format(round(hazard,2), nsmall = 2)), " (",format(round(lwr,2), nsmall = 2),"-",ifelse(upr>999,Inf,format(round(upr,2), nsmall = 2)),")\n",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))) %>%
    dplyr::select(arm, reference, result)
  
  newd = rbind(newd, newd1)
  
  print(i)
  
}

newd = unique(as.data.frame(newd))

result_matrix5d <- pivot_wider(newd, names_from = arm, values_from = result)

flextable(result_matrix5d) %>% autofit()  %>% save_as_docx( path = "OUTPUT/survival_infectivity_adjmic_.docx")






####combined log rank - using pcr adjustments####

nectime5a$type = "gametocytes (pcr)"
nectime5b$type = "infectivity"
nectime5c$type = "gametocytes (microscopy)"

nectime5 = rbind(nectime5a, nectime5b, nectime5c)

library(survival)
library(survminer)

km5 = survfit(Surv(days,status)~type, data= nectime5)
survout5 = survminer::ggsurvplot(km5,data=nectime5, surv.median.line = "hv", legend.labs=c("gametocytes (microscopy)","gametocytes (pcr)","infectivity"),pval=TRUE, conf.int = TRUE,pval.coord=c(30,0.9))

survout5$plot+
  xlab("Time (days) to clearance of...")+
  scale_y_continuous(labels=scales::percent)+
  theme_prism()+
  scale_color_viridis_d(begin=0.2, end=0.9)+
  scale_fill_viridis_d(begin=0.2, end=0.9)+
  theme(legend.position = "bottom")


ggsave("OUTPUT/Fig_024.png",device="png", width=12, height=8, dpi=600)


km6a = survfit(Surv(days,status)~arm, data= nectime5a)
survout6a = survminer::ggsurvplot(km6a)
dt6a = survout6a$data.survplot

km6b = survfit(Surv(days,status)~arm, data= nectime5b)
survout6b = survminer::ggsurvplot(km6b)
dt6b = survout6b$data.survplot


km6c = survfit(Surv(days,status)~arm, data= nectime5c)
survout6c = survminer::ggsurvplot(km6c)
dt6c = survout6c$data.survplot


dt6a$type = "gametocytes (pcr)"
dt6b$type = "infectivity"
dt6c$type = "gametocytes (microscopy)"


dt6 = rbind(dt6a, dt6b, dt6c) %>%
  select(time, surv, upper, lower, arm, type)

dt6sub = expand.grid(time=0, surv=1, upper=1, lower=1, arm=sort(unique(dt6$arm)),type=sort(unique(dt6$type)))

dt6 = unique(rbind(dt6,dt6sub))

dt6$arm = factor(dt6$arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

ggplot(data=dt6, aes(x=time, y=surv))+
  geom_step(aes(col=type))+
  geom_stepribbon(aes(ymin=lower, ymax=upper, fill=type),alpha=0.2)+
  xlab("Time (days) to clearance of ...")+
  scale_y_continuous(labels=scales::percent)+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  scale_color_viridis_d(begin=0.2, end=0.9)+
  scale_fill_viridis_d(begin=0.2, end=0.9)+
  theme(legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill = NA, size = 1))+  
  facet_wrap(~arm, nrow=2)+
  ylab("Probability of remaining uncleared")


ggsave("OUTPUT/Fig_024_arms.png",device="png", width=15, height=8, dpi=600)



# A simpler summary -------------------------------------------------------

#### RRG: Relative reduction in gametocytes####

#Similar to how we look at transmission reducing activity where we consider the decline in oocysts. Now we look at the day 2 and day 7 decline in gametocytes (both microscopy and pcr). 

##### Microscopy #####

necsub = nec %>% 
  filter((days==0 | days==2 | days==7| days==14) & !is.na(gam.m)) %>%
  select(studyid_str,days,arm,gam.m) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))


ref = levels(necsub$arm)

necsub$arm2 = necsub$arm
necsub$arm = factor(necsub$arm, levels=levels(necsub$arm2))

mod2f = gam(I(gam.m+0.1) ~ visit*arm + s(studyid_str, bs="re"),
            data = necsub,
            family = Gamma(link="log")
)

sb2f = summary(mod2f)


newd0 = data.frame(arm=levels(necsub$arm)) 

newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(necsub$arm)[1],"visitday 2", paste0("visitday 2:arm", arm)),
    beta = sb2f$p.coeff[[values]],
    se = sb2f$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)"),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))

day7newd0 = data.frame(arm=levels(necsub$arm)) 

day7newd0 = day7newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(necsub$arm)[1],"visitday 7", paste0("visitday 7:arm", arm)),
    beta = sb2f$p.coeff[[values]],
    se = sb2f$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)"),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))

day14newd0 = data.frame(arm=levels(necsub$arm)) 

day14newd0 = day14newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(necsub$arm)[1],"visitday 14", paste0("visitday 14:arm", arm)),
    beta = sb2f$p.coeff[[values]],
    se = sb2f$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)"),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))


pdat = rbind(
  newd0 %>% mutate(day="Day 2"),
  day7newd0 %>% mutate(day="Day 7"),
  day14newd0 %>% mutate(day="Day 14")
)  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")),
              day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 


ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.25, linewidth=0.5,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \ngametocyte densities by micrscopy")+
  coord_cartesian(ylim = c(-0.5, NA))


#annotation_custom(grob = linesGrob(gp = gpar(col = "black", lwd = 1)), ymin = Inf, ymax = Inf, xmin = -Inf, xmax = Inf)

ggsave("OUTPUT/Gametocyte_density_mic_day14.png",device="png", width=10, height=7, dpi=600)

newd0 = newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day7newd0 = day7newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day14newd0 = day14newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

current = ref[1]

newd =  data.frame(arm=levels(necsub$arm)) 

newd = newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 2:arm", arm), reference = current,
    pval = sb2f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)

day7newd =  data.frame(arm=levels(necsub$arm)) 

day7newd = day7newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 7:arm", arm), reference = current,
    pval = sb2f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)

day14newd =  data.frame(arm=levels(necsub$arm)) 

day14newd = day14newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 14:arm", arm), reference = current,
    pval = sb2f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


for (i in 2:length(ref)){
  
  current = ref[i]
  necsub$arm = factor(necsub$arm, levels=levels(necsub$arm2))
  necsub$arm = relevel(necsub$arm,current)
  mod2f = gam(I(gam.m+0.1) ~ visit*arm + s(studyid_str, bs="re"),
              data = necsub,
              family = Gamma(link="log")
  )
  
  sb2f = summary(mod2f)
  
  
  newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  newd1 = newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 2:arm", arm), reference = current,
      pval = sb2f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse( pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  newd = rbind(newd, newd1)
  
  day7newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  day7newd1 = day7newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 7:arm", arm), reference = current,
      pval = sb2f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse( pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day7newd = rbind(day7newd, day7newd1)
  
  day14newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  day14newd1 = day14newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 14:arm", arm), reference = current,
      pval = sb2f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse( pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day14newd = rbind(day14newd, day14newd1)
  
  print(i)
}


newd = rbind(newd,newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

newd = unique(as.data.frame(newd))

result_matrix2f <- pivot_wider(newd, names_from = arm, values_from = result)

result_matrix2f = result_matrix2f %>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix2f) %>% autofit() %>% save_as_docx( path = "OUTPUT/Gam_relative_reduction_mic_day2.docx")

day7newd = rbind(day7newd,day7newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day7newd = unique(as.data.frame(day7newd))

result_matrix2g <- pivot_wider(day7newd, names_from = arm, values_from = result)

result_matrix2g = result_matrix2g %>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix2g) %>% autofit() %>% save_as_docx( path = "OUTPUT/Gam_relative_reduction_mic_day7.docx")

day14newd = rbind(day14newd,day14newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day14newd = unique(as.data.frame(day14newd))

result_matrix2g <- pivot_wider(day14newd, names_from = arm, values_from = result)

result_matrix2g = result_matrix2g %>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix2g) %>% autofit() %>% save_as_docx( path = "OUTPUT/Gam_relative_reduction_mic_day14.docx")

##### pCR #####

necsub = nec %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(gam.p)) %>%
  select(studyid_str,days,arm,gam.p) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))


ref = levels(necsub$arm)

necsub$arm2 = necsub$arm
necsub$arm = factor(necsub$arm, levels=levels(necsub$arm2))

mod3f = gam(I(gam.p+0.01) ~ visit*arm + s(studyid_str, bs="re"),
            data = necsub,
            family = Gamma(link="log")
)

sb3f = summary(mod3f)


newd0 = data.frame(arm=levels(necsub$arm)) 

newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(necsub$arm)[1],"visitday 2", paste0("visitday 2:arm", arm)),
    beta = sb3f$p.coeff[[values]],
    se = sb3f$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)"),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))

day7newd0 = data.frame(arm=levels(necsub$arm)) 

day7newd0 = day7newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(necsub$arm)[1],"visitday 7", paste0("visitday 7:arm", arm)),
    beta = sb3f$p.coeff[[values]],
    se = sb3f$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)"),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))

day14newd0 = data.frame(arm=levels(necsub$arm)) 

day14newd0 = day14newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(necsub$arm)[1],"visitday 14", paste0("visitday 14:arm", arm)),
    beta = sb3f$p.coeff[[values]],
    se = sb3f$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)"),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))

pdat = rbind(
  newd0 %>% mutate(day="Day 2"),
  day7newd0 %>% mutate(day="Day 7"),
  day14newd0 %>% mutate(day="Day 14")
)  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")),
              day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 


ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.25, linewidth=0.5,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \ngametocyte densities by PCR")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/Gametocyte_density_PCR_day14.png",device="png", width=10, height=7, dpi=600)

newd0 = newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day7newd0 = day7newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day14newd0 = day14newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

current = ref[1]

newd =  data.frame(arm=levels(necsub$arm)) 

newd = newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 2:arm", arm), reference = current,
    pval = sb3f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)

day7newd =  data.frame(arm=levels(necsub$arm)) 

day7newd = day7newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 7:arm", arm), reference = current,
    pval = sb3f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)

day14newd =  data.frame(arm=levels(necsub$arm)) 

day14newd = day14newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 14:arm", arm), reference = current,
    pval = sb3f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


for (i in 2:length(ref)){
  
  current = ref[i]
  necsub$arm = factor(necsub$arm, levels=levels(necsub$arm2))
  necsub$arm = relevel(necsub$arm,current)
  mod3f = gam(I(gam.p+0.01) ~ visit*arm + s(studyid_str, bs="re"),
              data = necsub,
              family = Gamma(link="log")
  )
  
  sb3f = summary(mod3f)
  
  
  newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  newd1 = newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 2:arm", arm), reference = current,
      pval = sb3f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse( pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  newd = rbind(newd, newd1)
  
  day7newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  day7newd1 = day7newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 7:arm", arm), reference = current,
      pval = sb3f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse( pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day7newd = rbind(day7newd, day7newd1)
  
  day14newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  day14newd1 = day14newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 14:arm", arm), reference = current,
      pval = sb3f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse( pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day14newd = rbind(day14newd, day14newd1)
  
  print(i)
}


newd = rbind(newd,newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

newd = unique(as.data.frame(newd))

result_matrix3f <- pivot_wider(newd, names_from = arm, values_from = result)

result_matrix3f = result_matrix3f %>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix3f) %>% autofit() %>% save_as_docx( path = "OUTPUT/Gam_relative_reduction_pcr_day2.docx")

day7newd = rbind(day7newd,day7newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day7newd = unique(as.data.frame(day7newd))

result_matrix3g <- pivot_wider(day7newd, names_from = arm, values_from = result) %>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix3g) %>% autofit() %>% save_as_docx( path = "OUTPUT/Gam_relative_reduction_pcr_day7.docx")

day14newd = rbind(day14newd,day14newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day14newd = unique(as.data.frame(day14newd))

result_matrix3g <- pivot_wider(day14newd, names_from = arm, values_from = result) %>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix3g) %>% autofit() %>% save_as_docx( path = "OUTPUT/Gam_relative_reduction_pcr_day14.docx")


####RRGP: Relative reduction in gametocyte prevalence####

#####RRGP by microscopy (adjusted for baseline microscopy)#####


necsub = nec %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(gam.m)) %>%
  select(studyid_str,days,arm,gam.m) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm))) %>%
  group_by(studyid_str) %>%
  mutate(gam0.m=ifelse(first(days)!=0,NA,first(gam.m))) %>%
  ungroup()

necsub$arm = factor(necsub$arm, levels = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

ref = levels(necsub$arm)

necsub$arm2 = necsub$arm
necsub$arm = factor(necsub$arm, levels=levels(necsub$arm2))

mod2h = gam(I(gam.m>0) ~ visit*arm + I(log10(gam0.m+0.1)) +s(studyid_str, bs="re"), family=poisson(link="log"), data=necsub)
sb2h = summary(mod2h)

newd0 = data.frame(arm=levels(necsub$arm)) 
newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 2", paste0("visitday 2:arm", arm)),
    beta = sb2h$p.coeff[[values]],
    se = sb2h$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)"),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))


day7newd0 = data.frame(arm=levels(necsub$arm)) 
day7newd0 = day7newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 7", paste0("visitday 7:arm", arm)),
    beta = sb2h$p.coeff[[values]],
    se = sb2h$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)"),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))


day14newd0 = data.frame(arm=levels(necsub$arm)) 
day14newd0 = day14newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 14", paste0("visitday 14:arm", arm)),
    beta = sb2h$p.coeff[[values]],
    se = sb2h$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)"),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))

pdat = rbind(
  newd0 %>% mutate(day="Day 2"),
  day7newd0 %>% mutate(day="Day 7"),
  day14newd0 %>% mutate(day="Day 14")
)  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")),
              day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 


ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.25, linewidth=0.5,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-1.1,-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \ngametocyte prevalence by microscopy")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/Gametocyte_prev_mic_day14.png",device="png", width=10, height=7, dpi=600)


newd0 = newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day7newd0 = day7newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day14newd0 = day14newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

current = ref[1]

newd =  data.frame(arm=levels(necsub$arm)) 

newd = newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 2:arm", arm), reference = current,
    pval = sb2h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


day7newd =  data.frame(arm=levels(necsub$arm)) 

day7newd = day7newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 7:arm", arm), reference = current,
    pval = sb2h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)

day14newd =  data.frame(arm=levels(necsub$arm)) 

day14newd = day14newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 14:arm", arm), reference = current,
    pval = sb2h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)

for (i in 2:length(ref)){
  
  current = ref[i]
  
  necsub$arm = relevel(necsub$arm2,current)
  mod2h = gam(I(gam.m>0) ~ visit*arm + I(log10(gam0.m+0.1)) +s(studyid_str, bs="re"), family=poisson(link="log"), data=necsub)
  
  sb2h = summary(mod2h)
  
  
  newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  newd1 = newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 2:arm", arm), reference = current,
      pval = sb2h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  newd = rbind(newd, newd1)
  
  
  day7newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  day7newd1 = day7newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 7:arm", arm), reference = current,
      pval = sb2h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day7newd = rbind(day7newd, day7newd1)
  
  
  day14newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  day14newd1 = day14newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 14:arm", arm), reference = current,
      pval = sb2h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day14newd = rbind(day14newd, day14newd1)
}


newd = rbind(newd,
             newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

newd = unique(as.data.frame(newd))

result_matrix2h <- pivot_wider(newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)


day7newd = rbind(day7newd,
                 day7newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day7newd = unique(as.data.frame(day7newd))

result_matrix2h2 <- pivot_wider(day7newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)


day14newd = rbind(day14newd,
                 day14newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day14newd = unique(as.data.frame(day14newd))

result_matrix2h3 <- pivot_wider(day14newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)


flextable(result_matrix2h) %>% autofit() %>% save_as_docx( path = "OUTPUT/Gam_relative_reduction_prev_mic_day2.docx")
flextable(result_matrix2h2) %>% autofit() %>% save_as_docx( path = "OUTPUT/Gam_relative_reduction_prev_mic_day7.docx")
flextable(result_matrix2h3) %>% autofit() %>% save_as_docx( path = "OUTPUT/Gam_relative_reduction_prev_mic_day14.docx")



#####RRGP by PCR (adjusted for baseline PCR)#####

necsub = nec %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(gam.p)) %>%
  select(studyid_str,days,arm,gam.p) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))%>%
  group_by(studyid_str) %>%
  mutate(gam0.p=ifelse(first(days)!=0,NA,first(gam.p))) %>%
  ungroup()

necsub$arm = factor(necsub$arm, levels = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

ref = levels(necsub$arm)

necsub$arm2 = necsub$arm
necsub$arm = factor(necsub$arm, levels=levels(necsub$arm2))

mod3h = gam(I(gam.p>0) ~ visit*arm + I(log10(gam0.p+0.01)) +s(studyid_str, bs="re"), family=poisson(link="log"), data=necsub)

sb3h = summary(mod3h)

newd0 = data.frame(arm=levels(necsub$arm)) 
newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 2", paste0("visitday 2:arm", arm)),
    beta = sb3h$p.coeff[[values]],
    se = sb3h$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)"),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))

day7newd0 = data.frame(arm=levels(necsub$arm)) 
day7newd0 = day7newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 7", paste0("visitday 7:arm", arm)),
    beta = sb3h$p.coeff[[values]],
    se = sb3h$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)"),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))

day14newd0 = data.frame(arm=levels(necsub$arm)) 
day14newd0 = day14newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 14", paste0("visitday 14:arm", arm)),
    beta = sb3h$p.coeff[[values]],
    se = sb3h$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)"),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))


pdat = rbind(
  newd0 %>% mutate(day="Day 2"),
  day7newd0 %>% mutate(day="Day 7"),
  day14newd0 %>% mutate(day="Day 14")
)  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")),
              day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 


P<-ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.25, linewidth=0.5,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-1,-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \ngametocyte prevalence by PCR")+
  annotation_custom(grob = linesGrob(gp = gpar(col = "black", lwd = 1)), ymin = Inf, ymax = Inf, xmin = -Inf, xmax = Inf)


P+scale_y_break(c(-1, -0.5))


ggsave("OUTPUT/Gametocyte_prev_PCR_day14.png",device="png", width=12, height=8, dpi=600)

newd0 = newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day7newd0 = day7newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day14newd0 = day14newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

current = ref[1]

newd =  data.frame(arm=levels(necsub$arm)) 

newd = newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 2:arm", arm), reference = current,
    pval = sb3h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


day7newd =  data.frame(arm=levels(necsub$arm)) 

day7newd = day7newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 7:arm", arm), reference = current,
    pval = sb3h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


day14newd =  data.frame(arm=levels(necsub$arm)) 

day14newd = day14newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 14:arm", arm), reference = current,
    pval = sb3h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


for (i in 2:length(ref)){
  
  current = ref[i]
  
  necsub$arm = relevel(necsub$arm2,current)
  mod3h = gam(I(gam.p>0) ~ visit*arm + I(log10(gam0.p+0.01)) +s(studyid_str, bs="re"), family=poisson(link="log"), data=necsub)
  
  sb3h = summary(mod3h)
  
  
  newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  newd1 = newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 2:arm", arm), reference = current,
      pval = sb3h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  newd = rbind(newd, newd1)
  
  
  day7newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  day7newd1 = day7newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 7:arm", arm), reference = current,
      pval = sb3h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day7newd = rbind(day7newd, day7newd1)
  
  
  day14newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  day14newd1 = day14newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 14:arm", arm), reference = current,
      pval = sb3h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day14newd = rbind(day14newd, day14newd1)
}


newd = rbind(newd,
             newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

newd = unique(as.data.frame(newd))

result_matrix3h <- pivot_wider(newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)


day7newd = rbind(day7newd,
                 day7newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day7newd = unique(as.data.frame(day7newd))

result_matrix3h2 <- pivot_wider(day7newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)


day14newd = rbind(day14newd,
                 day14newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day14newd = unique(as.data.frame(day14newd))

result_matrix3h3 <- pivot_wider(day14newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)


flextable(result_matrix3h) %>% autofit() %>% save_as_docx( path = "OUTPUT/Gam_relative_reduction_prev_pcr_day2.docx")
flextable(result_matrix3h2) %>% autofit() %>% save_as_docx( path = "OUTPUT/Gam_relative_reduction_prev_pcr_day7.docx")
flextable(result_matrix3h3) %>% autofit() %>% save_as_docx( path = "OUTPUT/Gam_relative_reduction_prev_pcr_day14.docx")




#### TBA: transmission blocking activity ####

##### Unadjusted #####

necfeed = necfeed %>% 
  left_join(
    nec %>% 
      select(studyid_str,days, gam.m, gam.p) %>% 
      filter(days==0) %>%
      mutate(gam0.m = gam.m, gam0.p=gam.p) %>%
      select(studyid_str, gam0.m, gam0.p),
    by="studyid_str"
  ) 

necfeedsub  = necfeed %>% 
  filter((days==0 | days==2 | days==7) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7")), arm=as.factor(as.character(arm)))

necfeedsub$arm = factor(necfeedsub$arm, levels = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))


necfeedsub$prop = ifelse(necfeedsub$prop==0, 0.0001, ifelse(necfeedsub$prop==1, 0.9999, necfeedsub$prop))

ref = levels(necfeedsub$arm)
necfeedsub$arm2 = necfeedsub$arm 
necfeedsub$arm = relevel(necfeedsub$arm,ref[1])


mod4f = gam(prop ~ visit*arm +s(studyid_str, bs="re"), family=binomial(link="log"), weights = mosq_total, data=necfeedsub)

sb4f = summary(mod4f)

newd0 = data.frame(arm=levels(necfeedsub$arm)) 
newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(necfeedsub$arm)[1],"visitday 2", paste0("visitday 2:arm", arm)),
    beta = sb4f$p.coeff[[values]],
    se = sb4f$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))


day7newd0 = data.frame(arm=levels(necfeedsub$arm)) 
day7newd0 = day7newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(necfeedsub$arm)[1],"visitday 7", paste0("visitday 7:arm", arm)),
    beta = sb4f$p.coeff[[values]],
    se = sb4f$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))

pdat = rbind(
  newd0 %>% mutate(day="Day 2"),
  day7newd0 %>% mutate(day="Day 7")
)  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) 

ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  geom_point(position = position_dodge(0.8))+
  geom_bar(stat="identity", position = position_dodge(), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(0.8), width=0.25, linewidth=1.1)+
  theme_prism()+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent)+
  ylab("Relative reduction (%) in \nproportion infected mosquitos")+
  annotation_custom(grob = linesGrob(gp = gpar(col = "black", lwd = 1)), ymin = Inf, ymax = Inf, xmin = -Inf, xmax = Inf)


ggsave("OUTPUT/TBA_unadjusted.png",device="png", width=12, height=8, dpi=600)



newd0 = newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day7newd0 = day7newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)


current = ref[1]

newd =  data.frame(arm=levels(necfeedsub$arm)) 

newd = newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 2:arm", arm), reference = current,
    pval = sb4f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


day7newd =  data.frame(arm=levels(necfeedsub$arm)) 

day7newd = day7newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 7:arm", arm), reference = current,
    pval = sb4f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


for (i in 2:length(ref)){
  
  current = ref[i]
  
  necfeedsub$arm = relevel(necfeedsub$arm2,current)
  mod4f = gam(prop ~ visit*arm +s(studyid_str, bs="re"), family=binomial(link="log"), weights = mosq_total, data=necfeedsub)
  
  sb4f = summary(mod4f)
  
  
  newd1 =  data.frame(arm=levels(necfeedsub$arm)) 
  
  newd1 = newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 2:arm", arm), reference = current,
      pval = sb4f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  newd = rbind(newd, newd1)
  
  
  day7newd1 =  data.frame(arm=levels(necfeedsub$arm)) 
  
  day7newd1 = day7newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 7:arm", arm), reference = current,
      pval = sb4f$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day7newd = rbind(day7newd, day7newd1)
  
  print(i)
}


newd = rbind(newd,
             newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

newd = unique(as.data.frame(newd))

result_matrix4f <- pivot_wider(newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix4f) %>% autofit() %>% save_as_docx( path = "OUTPUT/TBA_day2_unadjusted.docx")



day7newd = rbind(day7newd,
             day7newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day7newd = unique(as.data.frame(day7newd))

result_matrix4f2 <- pivot_wider(day7newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix4f2) %>% autofit() %>% save_as_docx( path = "OUTPUT/TBA_day7_unadjusted.docx")


##### Adjusted by baseline gametocyte density (pcr) #####



necfeedsub  = necfeed %>% 
  filter((days==0 | days==2 | days==7) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop) & !is.na(gam0.p)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7")), arm=as.factor(as.character(arm)))

necfeedsub$prop = ifelse(necfeedsub$prop==0, 0.0001, ifelse(necfeedsub$prop==1, 0.9999, necfeedsub$prop))
necfeedsub$arm = factor(necfeedsub$arm, levels = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
ref = levels(necfeedsub$arm)
necfeedsub$arm2 = necfeedsub$arm 
necfeedsub$arm = relevel(necfeedsub$arm,ref[1])


mod4g = gam(prop ~ visit*arm + I(log10(gam0.p+0.01)) +s(studyid_str, bs="re"), family=binomial(link="log"), weights = mosq_total, data=necfeedsub)

sb4g = summary(mod4g)

newd0 = data.frame(arm=levels(necfeedsub$arm)) 
newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 2", paste0("visitday 2:arm", arm)),
    beta = sb4g$p.coeff[[values]],
    se = sb4g$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
  ) %>%
  dplyr::select(arm, RRG)


day7newd0 = data.frame(arm=levels(necfeedsub$arm)) 
day7newd0 = day7newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 7", paste0("visitday 7:arm", arm)),
    beta = sb4g$p.coeff[[values]],
    se = sb4g$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
  ) %>%
  dplyr::select(arm, RRG)



current = ref[1]

newd =  data.frame(arm=levels(necfeedsub$arm)) 

newd = newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 2:arm", arm), reference = current,
    pval = sb4g$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


day7newd =  data.frame(arm=levels(necfeedsub$arm)) 

day7newd = day7newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 7:arm", arm), reference = current,
    pval = sb4g$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


for (i in 2:length(ref)){
  
  current = ref[i]
  
  necfeedsub$arm = relevel(necfeedsub$arm2,current)
  mod4g = gam(prop ~ visit*arm + I(log10(gam0.p+0.01)) +s(studyid_str, bs="re"), family=binomial(link="log"), weights = mosq_total, data=necfeedsub)
  
  sb4g = summary(mod4g)
  
  
  newd1 =  data.frame(arm=levels(necfeedsub$arm)) 
  
  newd1 = newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 2:arm", arm), reference = current,
      pval = sb4g$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  newd = rbind(newd, newd1)
  
  
  day7newd1 =  data.frame(arm=levels(necfeedsub$arm)) 
  
  day7newd1 = day7newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 7:arm", arm), reference = current,
      pval = sb4g$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day7newd = rbind(day7newd, day7newd1)
  
  print(i)
}


newd = rbind(newd,
             newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

newd = unique(as.data.frame(newd))

result_matrix4g <- pivot_wider(newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix4g) %>% autofit() %>% save_as_docx( path = "OUTPUT/TBA_day2_adjustedpcr.docx")



day7newd = rbind(day7newd,
                 day7newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day7newd = unique(as.data.frame(day7newd))

result_matrix4g2 <- pivot_wider(day7newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix4g2) %>% autofit() %>% save_as_docx( path = "OUTPUT/TBA_day7_adjustedpcr.docx")



##### Adjusted by baseline gametocyte density (mic) #####


necfeedsub  = necfeed %>% 
  filter((days==0 | days==2 | days==7) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop) & !is.na(gam0.m)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7")), arm=as.factor(as.character(arm)))

necfeedsub$prop = ifelse(necfeedsub$prop==0, 0.0001, ifelse(necfeedsub$prop==1, 0.9999, necfeedsub$prop))
necfeedsub$arm = factor(necfeedsub$arm, levels = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
ref = levels(necfeedsub$arm)
necfeedsub$arm2 = necfeedsub$arm 
necfeedsub$arm = relevel(necfeedsub$arm,ref[1])


mod4h = gam(prop ~ visit*arm + I(log10(gam0.m+0.1)) +s(studyid_str, bs="re"), family=binomial(link="log"), weights = mosq_total, data=necfeedsub)

sb4h = summary(mod4h)

newd0 = data.frame(arm=levels(necfeedsub$arm)) 
newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 2", paste0("visitday 2:arm", arm)),
    beta = sb4h$p.coeff[[values]],
    se = sb4h$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
  ) %>%
  dplyr::select(arm, RRG)


day7newd0 = data.frame(arm=levels(necfeedsub$arm)) 
day7newd0 = day7newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 7", paste0("visitday 7:arm", arm)),
    beta = sb4h$p.coeff[[values]],
    se = sb4h$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
  ) %>%
  dplyr::select(arm, RRG)



current = ref[1]

newd =  data.frame(arm=levels(necfeedsub$arm)) 

newd = newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 2:arm", arm), reference = current,
    pval = sb4h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


day7newd =  data.frame(arm=levels(necfeedsub$arm)) 

day7newd = day7newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 7:arm", arm), reference = current,
    pval = sb4h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


for (i in 2:length(ref)){
  
  current = ref[i]
  
  necfeedsub$arm = relevel(necfeedsub$arm2,current)
  mod4h = gam(prop ~ visit*arm + I(log10(gam0.m+0.1)) +s(studyid_str, bs="re"), family=binomial(link="log"), weights = mosq_total, data=necfeedsub)
  
  sb4h = summary(mod4h)
  
  
  newd1 =  data.frame(arm=levels(necfeedsub$arm)) 
  
  newd1 = newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 2:arm", arm), reference = current,
      pval = sb4h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  newd = rbind(newd, newd1)
  
  
  day7newd1 =  data.frame(arm=levels(necfeedsub$arm)) 
  
  day7newd1 = day7newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 7:arm", arm), reference = current,
      pval = sb4h$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day7newd = rbind(day7newd, day7newd1)
  
  print(i)
}


newd = rbind(newd,
             newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

newd = unique(as.data.frame(newd))

result_matrix4h <- pivot_wider(newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix4h) %>% autofit() %>% save_as_docx( path = "OUTPUT/TBA_day2_adjustedmic.docx")



day7newd = rbind(day7newd,
                 day7newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day7newd = unique(as.data.frame(day7newd))

result_matrix4h2 <- pivot_wider(day7newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix4h2) %>% autofit() %>% save_as_docx( path = "OUTPUT/TBA_day7_adjustedmic.docx")



##### With day 14 #####

necfeedsub  = necfeed %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7", "day 14")), arm=as.factor(as.character(arm)))

necfeedsub$arm = factor(necfeedsub$arm, levels = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))


necfeedsub$prop = ifelse(necfeedsub$prop==0, 0.0001, ifelse(necfeedsub$prop==1, 0.9999, necfeedsub$prop))

ref = levels(necfeedsub$arm)
necfeedsub$arm2 = necfeedsub$arm 
necfeedsub$arm = relevel(necfeedsub$arm,ref[1])


mod4i = gam(prop ~ visit*arm +s(studyid_str, bs="re"), family=binomial(link="log"), weights = mosq_total, data=necfeedsub)

sb4i = summary(mod4i)

newd0 = data.frame(arm=levels(necfeedsub$arm)) 
newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 2", paste0("visitday 2:arm", arm)),
    beta = sb4i$p.coeff[[values]],
    se = sb4i$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))


day7newd0 = data.frame(arm=levels(necfeedsub$arm)) 
day7newd0 = day7newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 7", paste0("visitday 7:arm", arm)),
    beta = sb4i$p.coeff[[values]],
    se = sb4i$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))

day14newd0 = data.frame(arm=levels(necfeedsub$arm)) 
day14newd0 = day14newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 14", paste0("visitday 14:arm", arm)),
    beta = sb4i$p.coeff[[values]],
    se = sb4i$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))

pdat = rbind(
  newd0 %>% mutate(day="Day 2"),
  day7newd0 %>% mutate(day="Day 7"),
  day14newd0 %>% mutate(day="Day 14")
)  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")),
              day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) %>%
  filter(arm!="Non-ACT-PQ" | day!="Day 14")



ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.25, linewidth=0.5,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent)+
  ylab("Relative reduction (%) in \nproportion infected mosquitoes")


ggsave("OUTPUT/TBA_unadjusted_day14.png",device="png", width=10, height=7, dpi=600)



newd0 = newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day7newd0 = day7newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day14newd0 = day14newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)


current = ref[1]

newd =  data.frame(arm=levels(necfeedsub$arm)) 

newd = newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 2:arm", arm), reference = current,
    pval = sb4i$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


day7newd =  data.frame(arm=levels(necfeedsub$arm)) 

day7newd = day7newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 7:arm", arm), reference = current,
    pval = sb4i$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)



day14newd =  data.frame(arm=levels(necfeedsub$arm)) 

day14newd = day14newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 14:arm", arm), reference = current,
    pval = sb4i$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)



for (i in 2:length(ref)){
  
  current = ref[i]
  
  necfeedsub$arm = relevel(necfeedsub$arm2,current)
  mod4j = gam(prop ~ visit*arm +s(studyid_str, bs="re"), family=binomial(link="log"), weights = mosq_total, data=necfeedsub)
  
  sb4j = summary(mod4j)
  
  
  newd1 =  data.frame(arm=levels(necfeedsub$arm)) 
  
  newd1 = newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 2:arm", arm), reference = current,
      pval = sb4j$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  newd = rbind(newd, newd1)
  
  
  day7newd1 =  data.frame(arm=levels(necfeedsub$arm)) 
  
  day7newd1 = day7newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 7:arm", arm), reference = current,
      pval = sb4j$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day7newd = rbind(day7newd, day7newd1)
  
  
  day14newd1 =  data.frame(arm=levels(necfeedsub$arm)) 
  
  day14newd1 = day14newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 14:arm", arm), reference = current,
      pval = sb4j$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day14newd = rbind(day14newd, day14newd1)
  
  
  print(i)
}


newd = rbind(newd,
             newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

newd = unique(as.data.frame(newd))

result_matrix4j <- pivot_wider(newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix4j) %>% autofit() %>% save_as_docx( path = "OUTPUT/TBA_day2_unadjusted_withday14.docx")



day7newd = rbind(day7newd,
                 day7newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day7newd = unique(as.data.frame(day7newd))

result_matrix4j2 <- pivot_wider(day7newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix4j2) %>% autofit() %>% save_as_docx( path = "OUTPUT/TBA_day7_unadjusted_withday14.docx")


day14newd = rbind(day14newd,
                 day14newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day14newd = unique(as.data.frame(day14newd))

result_matrix4j3 <- pivot_wider(day14newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix4j3) %>% autofit() %>% save_as_docx( path = "OUTPUT/TBA_day14_unadjusted_withday14.docx")


#Heatmap for day2
result_matrix4j <- as.data.frame(result_matrix4j)
rownames(result_matrix4j) <- result_matrix4j[, 1] 
result_matrix4j <- result_matrix4j[, -1]

desired_order <- c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

result_long <- result_matrix4j %>%
  rownames_to_column("treatment") %>%
  pivot_longer(cols = -treatment, names_to = "comparison", values_to = "value") %>%
  mutate(
    p_value_numeric = ifelse(grepl("p<", value), 0.0001,  # Explicit handling of values "<0.0001"
                             ifelse(grepl("p=", value), 
                                    as.numeric(gsub(".*p=([0-9.]+).*", "\\1", value)),
                                    NA_real_)),
    label = ifelse(treatment == comparison,
                   gsub(" \\(", "\n(", value) %>% gsub(", ", ",\n", .),  # Add newlines for diagonal
                   ifelse(grepl("p<0.0001", value), "p<0.0001", sprintf("p=%.4f", p_value_numeric))),
    p_value_log = ifelse(is.na(p_value_numeric), NA, -log10(p_value_numeric)),  # Avoid log10 of NA
    border = ifelse(treatment == comparison, "black", NA)  # Black border for diagonal
  ) %>%
  drop_na(p_value_numeric) %>%
  mutate(
    treatment = factor(treatment, levels = desired_order),
    comparison = factor(comparison, levels = desired_order)
  ) %>%
  filter(as.numeric(treatment) <= as.numeric(comparison))


# Function to extract the numeric value from the activity string
extract_activity <- function(value) {
  matches <- regmatches(value, regexpr("(-?[0-9]+\\.?[0-9]*)%", value))
  if (length(matches) > 0) {
    return(as.numeric(sub("%", "", matches)))
  }
  return(NA)
}

# Extract the numeric activity values
result_long <- result_long %>%
  mutate(Activity_numeric = sapply(value, extract_activity))

# Initialize Activity_diff column
result_long$Activity_diff <- NA

# Calculate the absolute difference in activity
for (i in 1:nrow(result_long)) {
  treatment <- result_long$treatment[i]
  comparison <- result_long$comparison[i]
  
  # Find the activity values for the treatment and comparison
  activity_treatment <- result_long$Activity_numeric[
    result_long$treatment == treatment & result_long$comparison == treatment]
  activity_comparison <- result_long$Activity_numeric[
    result_long$treatment == comparison & result_long$comparison == comparison]
  
  if (length(activity_treatment) > 0 && length(activity_comparison) > 0) {
    # Calculate the absolute difference
    result_long$Activity_diff[i] <- abs(activity_treatment - activity_comparison)
  }
}

result_long <- result_long %>%
  mutate(is_diagonal = ifelse(treatment == comparison, TRUE, FALSE))

# Plotting the heatmap
ggplot(result_long, aes(x = comparison, y = treatment, fill = Activity_diff)) +
  geom_tile(color = "white", size = 0.1) +  # Default borders for non-diagonal
  geom_tile(data = filter(result_long, !is.na(border)), color = "white", size = 1, fill = NA) + 
  geom_tile(data = filter(result_long, is_diagonal), color = "black", size = 0.1, fill = "white") + 
  scale_fill_gradient2(low="#bfbbbb",high = "#508ea1", name = "Absolute difference in relative reduction",limits=c(0,120),breaks = c(0, 20, 40, 60, 80,100,120), labels = c("0%", "20%", "40%", "60%", "80%","100%","120%")) +
  geom_text(
    aes(label = label),  # Use the formatted label
    size = 3, 
    angle = 45  # Rotating text to 45 degrees
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "P-value Heatmap", x = "Comparison Group", y = "Treatment Group")

ggsave("OUTPUT/Heatmap_TBA_Day2.png",device="png", width=12, height=8, dpi=600)


#Heatmap for day7
result_matrix4j2 <- as.data.frame(result_matrix4j2)
rownames(result_matrix4j2) <- result_matrix4j2[, 1] 
result_matrix4j2 <- result_matrix4j2[, -1]

desired_order <- c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

result_long <- result_matrix4j2 %>%
  rownames_to_column("treatment") %>%
  pivot_longer(cols = -treatment, names_to = "comparison", values_to = "value") %>%
  mutate(
    p_value_numeric = ifelse(grepl("p<", value), 0.0001,  # Explicit handling of values "<0.0001"
                             ifelse(grepl("p=", value), 
                                    as.numeric(gsub(".*p=([0-9.]+).*", "\\1", value)),
                                    NA_real_)),
    label = ifelse(treatment == comparison,
                   gsub(" \\(", "\n(", value) %>% gsub(", ", ",\n", .),  # Add newlines for diagonal
                   ifelse(grepl("p<0.0001", value), "p<0.0001", sprintf("p=%.4f", p_value_numeric))),
    p_value_log = ifelse(is.na(p_value_numeric), NA, -log10(p_value_numeric)),  # Avoid log10 of NA
    border = ifelse(treatment == comparison, "black", NA)  # Black border for diagonal
  ) %>%
  drop_na(p_value_numeric) %>%
  mutate(
    treatment = factor(treatment, levels = desired_order),
    comparison = factor(comparison, levels = desired_order)
  ) %>%
  filter(as.numeric(treatment) <= as.numeric(comparison))


# Function to extract the numeric value from the activity string
extract_activity <- function(value) {
  matches <- regmatches(value, regexpr("(-?[0-9]+\\.?[0-9]*)%", value))
  if (length(matches) > 0) {
    return(as.numeric(sub("%", "", matches)))
  }
  return(NA)
}

# Extract the numeric activity values
result_long <- result_long %>%
  mutate(Activity_numeric = sapply(value, extract_activity))

# Initialize Activity_diff column
result_long$Activity_diff <- NA

# Calculate the absolute difference in activity
for (i in 1:nrow(result_long)) {
  treatment <- result_long$treatment[i]
  comparison <- result_long$comparison[i]
  
  # Find the activity values for the treatment and comparison
  activity_treatment <- result_long$Activity_numeric[
    result_long$treatment == treatment & result_long$comparison == treatment]
  activity_comparison <- result_long$Activity_numeric[
    result_long$treatment == comparison & result_long$comparison == comparison]
  
  if (length(activity_treatment) > 0 && length(activity_comparison) > 0) {
    # Calculate the absolute difference
    result_long$Activity_diff[i] <- abs(activity_treatment - activity_comparison)
  }
}

result_long <- result_long %>%
  mutate(is_diagonal = ifelse(treatment == comparison, TRUE, FALSE))

# Plotting the heatmap
ggplot(result_long, aes(x = comparison, y = treatment, fill = Activity_diff)) +
  geom_tile(color = "white", size = 0.1) +  # Default borders for non-diagonal
  geom_tile(data = filter(result_long, !is.na(border)), color = "white", size = 1, fill = NA) + 
  geom_tile(data = filter(result_long, is_diagonal), color = "black", size = 0.1, fill = "white") + 
  scale_fill_gradient2(low="#bfbbbb",high = "#508ea1", name = "Absolute difference in relative reduction",limits=c(0,120),breaks = c(0, 20, 40, 60, 80,100,120), labels = c("0%", "20%", "40%", "60%", "80%","100%","120%")) +
  geom_text(
    aes(label = label),  # Use the formatted label
    size = 3, 
    angle = 45  # Rotating text to 45 degrees
  ) +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "P-value Heatmap", x = "Comparison Group", y = "Treatment Group")

ggsave("OUTPUT/Heatmap_TBA_Day7.png",device="png", width=12, height=8, dpi=600)




####### Combining figures for reductions in gam prevalence, gam density and TBA



# 3d analysis: Gam densities, time and prop infected mosquitos ------------

mod6a = gam(prop ~ ti(I(log(days+1)))+ ti(I(log(gam.p+0.01))) + ti(I(log(days+1)),I(log(gam.p+0.01))) +s(studyid_str, bs="re"), family=binomial(link="logit"), weights = mosq_total, data=necfeed)

summary(mod6a)


newd = expand.grid(days = seq(0,max(necfeed$days,na.rm=T),length.out=1000),
                   gam.p=exp(seq(-10, log(quantile(necfeed$gam.p,probs=1,na.rm=T)[1]), length.out=1000)),
                   studyid_str = necfeed$studyid_str[1])

newd$prop=predict(mod6a,newd, type="response", exclude=c("s(studyid_str)"))

ilink <- family(mod6a)$linkinv
## add fit and se.fit on the **link** scale
newd <- bind_cols(newd, setNames(as_tibble(predict(mod6a, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
newd <- mutate(newd,
                fit_resp  = ilink(fit_link),
                right_upr = ilink(fit_link + (2 * se_link)),
                right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=newd, aes(x=days, y=gam.p, fill= fit_resp*100))+
  geom_tile()+
  theme_prism()+
  ylab("Gametocyte density/µL")+
  xlab("Days after treatment")+
  scale_y_log10()+
  #scale_x_log10()+
  scale_fill_viridis_c(name="Proportion\ninfected\nmosquitos")+
  theme(legend.position = "top",
        legend.title = element_text())

ggsave("OUTPUT/Fig_025.png",device="png", width=12, height=8, dpi=600)



# Infection rates per gametocytes -----------------------------------------

mod6b = gam(I(prop*100/(gam.p+0.01)) ~ I(log10(days+1))+s(studyid_str, bs="re"), family=gaussian(link="log"), data=necfeed)
summary(mod6b)


newd = expand.grid(days = seq(0,max(necfeed$days,na.rm=T),length.out=1000),
                   studyid_str = necfeed$studyid_str[1])

newd$prop=predict(mod6b,newd, type="response", exclude=c("s(studyid_str)"))

ilink <- family(mod6b)$linkinv
## add fit and se.fit on the **link** scale
newd <- bind_cols(newd, setNames(as_tibble(predict(mod6b, newd, se.fit = TRUE, exclude=c("s(studyid_str)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
newd <- mutate(newd,
               fit_resp  = ilink(fit_link),
               right_upr = ilink(fit_link + (2 * se_link)),
               right_lwr = ilink(fit_link - (2 * se_link)))

ggplot(data=newd, aes(x=days, y=fit_resp))+
  geom_line() +
  geom_ribbon(aes(ymin=right_lwr,ymax=right_upr), alpha=0.2)+
  theme_prism()+
  ylab("Infectivity per gametocyte density (%/µL)")+
  xlab("Days after treatment")+
  scale_y_log10()+
  scale_fill_viridis_c()

ggsave("OUTPUT/Fig_026.png",device="png", width=12, height=8, dpi=600)


####RRIG: Relative reduction in infectivity per gametocyte####

#by visualization first (pcr)

##with day 7

necfeedsub  = necfeed %>% 
  filter((days==0 | days==2 | days==7) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop) & !is.na(gam.p)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7")), arm=as.factor(as.character(arm)))

necfeedsub$prop = ifelse(necfeedsub$prop==0, 0.0001, ifelse(necfeedsub$prop==1, 0.9999, necfeedsub$prop))

modi = gam(prop ~ I(log(gam.p+0.01))+s(arm, bs="re")+s(study, bs="re"), family=binomial(), weights = mosq_total, data=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),])
summary(modi)

nec0i = expand.grid(gam.p=seq(0.01,max(necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),]$gam.p,na.rm=T), length.out=10000), study=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),]$study[1], arm=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),]$arm[1])
ilink <- family(modi)$linkinv
## add fit and se.fit on the **link** scale
nec0i <- bind_cols(nec0i, setNames(as_tibble(predict(modi, nec0i, se.fit = TRUE,exclude = c("s(arm)","s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i <- mutate(nec0i,
                prev  = ilink(fit_link),
                upr = ilink(fit_link + (2 * se_link)),
                lwr = ilink(fit_link - (2 * se_link)))


modi2 = gam(prop ~ I(log(gam.p+0.01))+s(arm, bs="re")+s(study, bs="re"), family=binomial(), weights = mosq_total, data=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),])
summary(modi2)

nec0i2 = expand.grid(gam.p=seq(0.01,max(necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),]$gam.p,na.rm=T), length.out=10000), study=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),]$study[1], arm=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),]$arm[1])
ilink <- family(modi2)$linkinv
## add fit and se.fit on the **link** scale
nec0i2 <- bind_cols(nec0i2, setNames(as_tibble(predict(modi2, nec0i2, se.fit = TRUE,exclude = c("s(arm)","s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i2 <- mutate(nec0i2,
                prev  = ilink(fit_link),
                upr = ilink(fit_link + (2 * se_link)),
                lwr = ilink(fit_link - (2 * se_link)))
nec0.full = rbind(nec0i %>% mutate(group="group1"),nec0i2 %>% mutate(group="group2")) %>% mutate(group=as.factor(group))

necfeedsub$arm = factor(necfeedsub$arm, levels = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

necfeedsub$group = as.factor(ifelse(necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), "group2", "group1"))

ggplot(data=nec0.full, aes(x=gam.p, y=prev))+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2)+
  geom_point(data=necfeedsub, aes(x=gam.p, y=prop, col=arm, shape=visit),alpha=0.4,size=2)+
  geom_line(linewidth=1.2)+
  theme_prism()+
  ylab("Proportion infected mosquitos")+
  xlab("Gametocyte density/µL (pcr)")+
  scale_x_log10()+
  scale_color_viridis_d()+
  theme(legend.position = "top")+
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(~group)

ggsave("OUTPUT/Fig_027_with_day7.png",device="png", width=12, height=8, dpi=600)


##withOUT day 7

necfeedsub  = necfeed %>% 
  filter((days==0 | days==2) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop) & !is.na(gam.p)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2")), arm=as.factor(as.character(arm)))

necfeedsub$prop = ifelse(necfeedsub$prop==0, 0.0001, ifelse(necfeedsub$prop==1, 0.9999, necfeedsub$prop))

modi = gam(prop ~ I(log(gam.p+0.01))+s(arm, bs="re")+s(study, bs="re"), family=binomial(), weights = mosq_total, data=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),])
summary(modi)

nec0i = expand.grid(gam.p=seq(0.01,max(necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),]$gam.p,na.rm=T), length.out=10000), study=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),]$study[1], arm=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),]$arm[1])
ilink <- family(modi)$linkinv
## add fit and se.fit on the **link** scale
nec0i <- bind_cols(nec0i, setNames(as_tibble(predict(modi, nec0i, se.fit = TRUE,exclude = c("s(arm)","s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i <- mutate(nec0i,
                prev  = ilink(fit_link),
                upr = ilink(fit_link + (2 * se_link)),
                lwr = ilink(fit_link - (2 * se_link)))


modi2 = gam(prop ~ I(log(gam.p+0.01))+s(arm, bs="re")+s(study, bs="re"), family=binomial(), weights = mosq_total, data=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),])
summary(modi2)

nec0i2 = expand.grid(gam.p=seq(0.01,max(necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),]$gam.p,na.rm=T), length.out=10000), study=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),]$study[1], arm=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),]$arm[1])
ilink <- family(modi2)$linkinv
## add fit and se.fit on the **link** scale
nec0i2 <- bind_cols(nec0i2, setNames(as_tibble(predict(modi2, nec0i2, se.fit = TRUE,exclude = c("s(arm)","s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i2 <- mutate(nec0i2,
                 prev  = ilink(fit_link),
                 upr = ilink(fit_link + (2 * se_link)),
                 lwr = ilink(fit_link - (2 * se_link)))
nec0.full = rbind(nec0i %>% mutate(group="group1"),nec0i2 %>% mutate(group="group2")) %>% mutate(group=as.factor(group))

necfeedsub$arm = factor(necfeedsub$arm, levels = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

necfeedsub$group = as.factor(ifelse(necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), "group2", "group1"))

ggplot(data=nec0.full, aes(x=gam.p, y=prev))+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2)+
  geom_point(data=necfeedsub, aes(x=gam.p, y=prop, col=arm, shape=visit),alpha=0.4,size=2)+
  geom_line(linewidth=1.2)+
  theme_prism()+
  ylab("Proportion infected mosquitos")+
  xlab("Gametocyte density/µL (pcr)")+
  scale_x_log10()+
  scale_color_viridis_d()+
  theme(legend.position = "top")+
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(~group)

ggsave("OUTPUT/Fig_027_withOUT_day7.png",device="png", width=12, height=8, dpi=600)

#by visualization first (microscopy)

##with day 7

necfeedsub  = necfeed %>% 
  filter((days==0 | days==2 | days==7) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop) & !is.na(gam.m)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7")), arm=as.factor(as.character(arm)))

necfeedsub$prop = ifelse(necfeedsub$prop==0, 0.0001, ifelse(necfeedsub$prop==1, 0.9999, necfeedsub$prop))

modi = gam(prop ~ I(log(gam.m+0.1))+s(arm, bs="re")+s(study, bs="re"), family=binomial(), weights = mosq_total, data=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),])
summary(modi)

nec0i = expand.grid(gam.m=seq(0.01,max(necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),]$gam.m,na.rm=T), length.out=10000), study=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),]$study[1], arm=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),]$arm[1])
ilink <- family(modi)$linkinv
## add fit and se.fit on the **link** scale
nec0i <- bind_cols(nec0i, setNames(as_tibble(predict(modi, nec0i, se.fit = TRUE,exclude = c("s(arm)","s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i <- mutate(nec0i,
                prev  = ilink(fit_link),
                upr = ilink(fit_link + (2 * se_link)),
                lwr = ilink(fit_link - (2 * se_link)))


modi2 = gam(prop ~ I(log(gam.m+0.1))+s(arm, bs="re")+s(study, bs="re"), family=binomial(), weights = mosq_total, data=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),])
summary(modi2)

nec0i2 = expand.grid(gam.m=seq(0.01,max(necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),]$gam.m,na.rm=T), length.out=10000), study=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),]$study[1], arm=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),]$arm[1])
ilink <- family(modi2)$linkinv
## add fit and se.fit on the **link** scale
nec0i2 <- bind_cols(nec0i2, setNames(as_tibble(predict(modi2, nec0i2, se.fit = TRUE,exclude = c("s(arm)","s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i2 <- mutate(nec0i2,
                 prev  = ilink(fit_link),
                 upr = ilink(fit_link + (2 * se_link)),
                 lwr = ilink(fit_link - (2 * se_link)))
nec0.full = rbind(nec0i %>% mutate(group="group1"),nec0i2 %>% mutate(group="group2")) %>% mutate(group=as.factor(group))

necfeedsub$arm = factor(necfeedsub$arm, levels = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

necfeedsub$group = as.factor(ifelse(necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), "group2", "group1"))

ggplot(data=nec0.full, aes(x=gam.m, y=prev))+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2)+
  geom_point(data=necfeedsub, aes(x=gam.m, y=prop, col=arm, shape=visit),alpha=0.4,size=2)+
  geom_line(linewidth=1.2)+
  theme_prism()+
  ylab("Proportion infected mosquitos")+
  xlab("Gametocyte density/µL (microscopy)")+
  scale_x_log10()+
  scale_color_viridis_d()+
  theme(legend.position = "top")+
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(~group)

ggsave("OUTPUT/Fig_028_with_day7.png",device="png", width=12, height=8, dpi=600)


##withOUT day 7

necfeedsub  = necfeed %>% 
  filter((days==0 | days==2) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop) & !is.na(gam.m)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2")), arm=as.factor(as.character(arm)))

necfeedsub$prop = ifelse(necfeedsub$prop==0, 0.0001, ifelse(necfeedsub$prop==1, 0.9999, necfeedsub$prop))

modi = gam(prop ~ I(log(gam.m+0.1))+s(arm, bs="re")+s(study, bs="re"), family=binomial(), weights = mosq_total, data=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),])
summary(modi)

nec0i = expand.grid(gam.m=seq(0.01,max(necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),]$gam.m,na.rm=T), length.out=10000), study=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),]$study[1], arm=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ"),]$arm[1])
ilink <- family(modi)$linkinv
## add fit and se.fit on the **link** scale
nec0i <- bind_cols(nec0i, setNames(as_tibble(predict(modi, nec0i, se.fit = TRUE,exclude = c("s(arm)","s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i <- mutate(nec0i,
                prev  = ilink(fit_link),
                upr = ilink(fit_link + (2 * se_link)),
                lwr = ilink(fit_link - (2 * se_link)))


modi2 = gam(prop ~ I(log(gam.m+0.1))+s(arm, bs="re")+s(study, bs="re"), family=binomial(), weights = mosq_total, data=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),])
summary(modi2)

nec0i2 = expand.grid(gam.m=seq(0.01,max(necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),]$gam.m,na.rm=T), length.out=10000), study=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),]$study[1], arm=necfeedsub[necfeedsub$days==0 & necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"),]$arm[1])
ilink <- family(modi2)$linkinv
## add fit and se.fit on the **link** scale
nec0i2 <- bind_cols(nec0i2, setNames(as_tibble(predict(modi2, nec0i2, se.fit = TRUE,exclude = c("s(arm)","s(study)"))[1:2]),c('fit_link','se_link')))
## create the interval and backtransform
nec0i2 <- mutate(nec0i2,
                 prev  = ilink(fit_link),
                 upr = ilink(fit_link + (2 * se_link)),
                 lwr = ilink(fit_link - (2 * se_link)))
nec0.full = rbind(nec0i %>% mutate(group="group1"),nec0i2 %>% mutate(group="group2")) %>% mutate(group=as.factor(group))

necfeedsub$arm = factor(necfeedsub$arm, levels = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))

necfeedsub$group = as.factor(ifelse(necfeedsub$arm %in% c("AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"), "group2", "group1"))

ggplot(data=nec0.full, aes(x=gam.m, y=prev))+
  geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.2)+
  geom_point(data=necfeedsub, aes(x=gam.m, y=prop, col=arm, shape=visit),alpha=0.4,size=2)+
  geom_line(linewidth=1.2)+
  theme_prism()+
  ylab("Proportion infected mosquitos")+
  xlab("Gametocyte density/µL (microscopy)")+
  scale_x_log10()+
  scale_color_viridis_d()+
  theme(legend.position = "top")+
  scale_y_continuous(labels = scales::percent)+
  facet_wrap(~group)

ggsave("OUTPUT/Fig_028_withOUT_day7.png",device="png", width=12, height=8, dpi=600)



#modelling

necfeedsub  = necfeed %>% 
  filter((days==0 | days==2 | days==7) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop) & !is.na(gam.p)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7")), arm=as.factor(as.character(arm)))

necfeedsub$prop = ifelse(necfeedsub$prop==0, 0.0001, ifelse(necfeedsub$prop==1, 0.9999, necfeedsub$prop))
necfeedsub$arm = factor(necfeedsub$arm, levels = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))
ref = levels(necfeedsub$arm)
necfeedsub$arm2 = necfeedsub$arm 
necfeedsub$arm = relevel(necfeedsub$arm,ref[1])


mod6c = gam(I(prop*100/(gam.p+0.01)) ~ visit*arm +s(studyid_str, bs="re"), family=gaussian(link="log"), data=necfeedsub)

sb6c = summary(mod6c)

newd0 = data.frame(arm=levels(necfeedsub$arm)) 
newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 2", paste0("visitday 2:arm", arm)),
    beta = sb6c$p.coeff[[values]],
    se = sb6c$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)")
  ) %>%
  dplyr::select(arm, RRG)


day7newd0 = data.frame(arm=levels(necfeedsub$arm)) 
day7newd0 = day7newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 7", paste0("visitday 7:arm", arm)),
    beta = sb6c$p.coeff[[values]],
    se = sb6c$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)")
  ) %>%
  dplyr::select(arm, RRG)



current = ref[1]

newd =  data.frame(arm=levels(necfeedsub$arm)) 

newd = newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 2:arm", arm), reference = current,
    pval = sb6c$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


day7newd =  data.frame(arm=levels(necfeedsub$arm)) 

day7newd = day7newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 7:arm", arm), reference = current,
    pval = sb6c$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


for (i in 2:length(ref)){
  
  current = ref[i]
  
  necfeedsub$arm = relevel(necfeedsub$arm2,current)
  mod6c = gam(I(prop*100/(gam.p+0.01)) ~ visit*arm +s(studyid_str, bs="re"), family=gaussian(link="log"), data=necfeedsub)
  
  sb6c = summary(mod6c)
  
  
  newd1 =  data.frame(arm=levels(necfeedsub$arm)) 
  
  newd1 = newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 2:arm", arm), reference = current,
      pval = sb6c$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  newd = rbind(newd, newd1)
  
  
  day7newd1 =  data.frame(arm=levels(necfeedsub$arm)) 
  
  day7newd1 = day7newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 7:arm", arm), reference = current,
      pval = sb6c$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day7newd = rbind(day7newd, day7newd1)
  
  print(i)
}


newd = rbind(newd,
             newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

newd = unique(as.data.frame(newd))

result_matrix6c <- pivot_wider(newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix6c) %>% autofit() %>% save_as_docx( path = "OUTPUT/RRIG_day2_adjustedpcr.docx")



day7newd = rbind(day7newd,
                 day7newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day7newd = unique(as.data.frame(day7newd))

result_matrix6c2 <- pivot_wider(day7newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix6c2) %>% autofit() %>% save_as_docx( path = "OUTPUT/RRIG_day7_adjustedpcr.docx")




# Oocyst density ----------------------------------------------------------


necsub = nec %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mean_ooc_all)) %>%
  select(studyid_str,days,arm,mean_ooc_all) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm))) %>%
  group_by(studyid_str) %>%
  ungroup()

necsub$arm = factor(necsub$arm, levels = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"))
#no oocyst densities for Non-ACT-PQ studies


ref = levels(necsub$arm)

necsub$arm2 = necsub$arm
necsub$arm = factor(necsub$arm, levels=levels(necsub$arm2))

necsub$mean_ooc_all[necsub$arm == "AL" & necsub$days==7][1]=0.01
necsub$mean_ooc_all[necsub$arm == "AL" & necsub$days==14][1]=0.01
necsub$mean_ooc_all[necsub$arm == "ACT-PQ" & necsub$days==14][1]=0.01
necsub$mean_ooc_all[necsub$arm == "ACT-TQ" & necsub$days==14][1]=0.01


k=100
necsub$Dissected = k
necsub$oocysts_t = necsub$Dissected*necsub$mean_ooc_all
mod2j = gam(oocysts_t ~ visit*arm +s(studyid_str, bs="re"),offset=necsub$Dissected, family=nb(), data=necsub)
sb2j = summary(mod2j)
sb2j


newd0 = data.frame(arm=levels(necsub$arm)) 
newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(necsub$arm)[1],"visitday 2", paste0("visitday 2:arm", arm)),
    beta = sb2j$p.coeff[[values]],
    se = sb2j$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)")
  )


day7newd0 = data.frame(arm=levels(necsub$arm)) 
day7newd0 = day7newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 7", paste0("visitday 7:arm", arm)),
    beta = sb2j$p.coeff[[values]],
    se = sb2j$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)")
  ) 


day14newd0 = data.frame(arm=levels(necsub$arm)) 
day14newd0 = day14newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 14", paste0("visitday 14:arm", arm)),
    beta = sb2j$p.coeff[[values]],
    se = sb2j$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%)")
  ) 




pdat = rbind(
  newd0 %>% mutate(day="Day 2"),
  day7newd0 %>% mutate(day="Day 7"),
  day14newd0 %>% mutate(day="Day 14")
)   %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")),
               day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) 

pdat[22:24,]=NA
pdat$arm[22:24]="Non-ACT-PQ"
pdat$TRA[22:24]=NA
pdat$ci_lower[22:24]=NA
pdat$ci_upper[22:24]=NA
pdat$day[22:24]=c("Day 2","Day 7","Day 14")


ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.25, linewidth=1.1)+
  theme_prism()+
  theme(axis.title.x = element_blank(),
        legend.position = "top")+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent)+
  ylab("Relative reduction (%) in oocysts\n per dissected mosquito")+
  coord_cartesian(ylim=c(-1,1))

ggsave("OUTPUT/Oocyst_density_reduction.png",device="png", width=12, height=8, dpi=600)


newd0 = newd0  %>%
  dplyr::select(arm, RRG)

day7newd0 = day7newd0  %>%
  dplyr::select(arm, RRG)

day14newd0 = day14newd0  %>%
  dplyr::select(arm, RRG)




current = ref[1]

newd =  data.frame(arm=levels(necsub$arm)) 

newd = newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 2:arm", arm), reference = current,
    pval = sb2j$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


day7newd =  data.frame(arm=levels(necsub$arm)) 

day7newd = day7newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 7:arm", arm), reference = current,
    pval = sb2j$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)

day14newd =  data.frame(arm=levels(necsub$arm)) 

day14newd = day14newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 14:arm", arm), reference = current,
    pval = sb2j$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)

for (i in 2:length(ref)){
  
  current = ref[i]
  
  necsub$arm = relevel(necsub$arm2,current)
  mod2j = gam(oocysts_t ~ visit*arm +s(studyid_str, bs="re"),offset=necsub$Dissected, family=nb(), data=necsub)
  sb2j = summary(mod2j)
  
  
  newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  newd1 = newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 2:arm", arm), reference = current,
      pval = sb2j$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  newd = rbind(newd, newd1)
  
  
  day7newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  day7newd1 = day7newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 7:arm", arm), reference = current,
      pval = sb2j$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day7newd = rbind(day7newd, day7newd1)
  
  
  day14newd1 =  data.frame(arm=levels(necsub$arm)) 
  
  day14newd1 = day14newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 14:arm", arm), reference = current,
      pval = sb2j$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day14newd = rbind(day14newd, day14newd1)
}


newd = rbind(newd,
             newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

newd = unique(as.data.frame(newd))

result_matrix2j <- pivot_wider(newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)


day7newd = rbind(day7newd,
                 day7newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day7newd = unique(as.data.frame(day7newd))

result_matrix2j2 <- pivot_wider(day7newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)


day14newd = rbind(day14newd,
                  day14newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day14newd = unique(as.data.frame(day14newd))

result_matrix2j3 <- pivot_wider(day14newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)


flextable(result_matrix2j) %>% autofit() %>% save_as_docx( path = "OUTPUT/oocyst_relative_reduction_day2.docx")
flextable(result_matrix2j2) %>% autofit() %>% save_as_docx( path = "OUTPUT/oocyst_relative_reduction_day7.docx")
flextable(result_matrix2j3) %>% autofit() %>% save_as_docx( path = "OUTPUT/oocyst_relative_reduction_day14.docx")



# Infectiousness ----------------------------------------------------------

necfeedsub  = necfeed %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop)) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2", "day 7", "day 14")), arm=as.factor(as.character(arm)))

necfeedsub$arm = factor(necfeedsub$arm, levels = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))


ref = levels(necfeedsub$arm)
necfeedsub$arm2 = necfeedsub$arm 
necfeedsub$arm = relevel(necfeedsub$arm,ref[1])

necfeedsub$prop[necfeedsub$arm == "AL" & necfeedsub$days==7][1]=1
necfeedsub$prop[necfeedsub$arm == "AL" & necfeedsub$days==14][1]=1
necfeedsub$prop[necfeedsub$arm == "Non-ACT-PQ" & necfeedsub$days==7][1]=1
necfeedsub$prop[necfeedsub$arm == "ACT-PQ" & necfeedsub$days==14][1]=1
necfeedsub$prop[necfeedsub$arm == "ACT-TQ" & necfeedsub$days==14][1]=1



mod4i = gam(I(prop>0.0001) ~ visit*arm +s(studyid_str, bs="re"), family=poisson(link="log"), data=necfeedsub)
sb4i = summary(mod4i)

newd0 = data.frame(arm=levels(necfeedsub$arm)) 
newd0 = newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 2", paste0("visitday 2:arm", arm)),
    beta = sb4i$p.coeff[[values]],
    se = sb4i$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))


day7newd0 = data.frame(arm=levels(necfeedsub$arm)) 
day7newd0 = day7newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 7", paste0("visitday 7:arm", arm)),
    beta = sb4i$p.coeff[[values]],
    se = sb4i$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))

day14newd0 = data.frame(arm=levels(necfeedsub$arm)) 
day14newd0 = day14newd0 %>% 
  group_by(arm) %>%
  mutate(
    values = ifelse(arm==levels(nec$arm)[1],"visitday 14", paste0("visitday 14:arm", arm)),
    beta = sb4i$p.coeff[[values]],
    se = sb4i$se[[values]]) %>%
  ungroup() %>%
  mutate(n = 1:n(), 
         beta=ifelse(n==1,beta,first(beta)+beta), 
         se=ifelse(n==1,se,sqrt(first(se)^2+se^2)),
         TRA = 1-exp(beta),
         ci_lower = 1-exp(beta+1.96*se),
         ci_upper = 1-exp(beta-1.96*se),
         p.pv = pnorm(abs(beta/se), lower.tail = F)*2,
         pval = ifelse(p.pv<0.0001,"<0.0001",ifelse(p.pv>0.9999,">0.9999", paste0("=",format(round(p.pv,4),nsmall=4)))))

pdat = rbind(
  newd0 %>% mutate(day="Day 2"),
  day7newd0 %>% mutate(day="Day 7"),
  day14newd0 %>% mutate(day="Day 14")
)  %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")),
              day = factor(day, levels=c("Day 2", "Day 7", "Day 14"))) %>%
  filter(arm!="Non-ACT-PQ" | day!="Day 14")

ggplot(data=pdat, aes(x=arm, y=TRA, fill=day, group=day))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper), position = position_dodge(preserve = "single", width=0.8), width=0.25, linewidth=1.1)+
  theme_prism()+
  theme(axis.title.x = element_blank(),
        legend.position = "top")+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent)+
  ylab("Relative reduction (%) in probability\nof infecting at least 1 mosquito")

ggsave("OUTPUT/infectiousness.png",device="png", width=12, height=8, dpi=600)


newd0 = newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day7newd0 = day7newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)

day14newd0 = day14newd0 %>% mutate(RRG = paste0(format((round(TRA*100,2)),2),"% (",format((round(ifelse(ci_lower*100 < -10000, -Inf, ci_lower*100),2)), nsmall = 2),"% -",format((round(ci_upper*100,2)), nsmall = 2),"%, p",pval,")")
) %>%
  dplyr::select(arm, RRG)


current = ref[1]

newd =  data.frame(arm=levels(necfeedsub$arm)) 

newd = newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 2:arm", arm), reference = current,
    pval = sb4i$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)


day7newd =  data.frame(arm=levels(necfeedsub$arm)) 

day7newd = day7newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 7:arm", arm), reference = current,
    pval = sb4i$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)



day14newd =  data.frame(arm=levels(necfeedsub$arm)) 

day14newd = day14newd %>% filter(!(arm %in% current)) %>% 
  group_by(arm) %>%
  mutate(
    values = paste0("visitday 14:arm", arm), reference = current,
    pval = sb4i$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
  dplyr::select(arm, reference, result)



for (i in 2:length(ref)){
  
  current = ref[i]
  
  necfeedsub$arm = relevel(necfeedsub$arm2,current)
  mod4j = gam(I(prop>0.0001) ~ visit*arm +s(studyid_str, bs="re"), family=poisson(link="log"), data=necfeedsub)
  
  sb4j = summary(mod4j)
  
  
  newd1 =  data.frame(arm=levels(necfeedsub$arm)) 
  
  newd1 = newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 2:arm", arm), reference = current,
      pval = sb4j$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  newd = rbind(newd, newd1)
  
  
  day7newd1 =  data.frame(arm=levels(necfeedsub$arm)) 
  
  day7newd1 = day7newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 7:arm", arm), reference = current,
      pval = sb4j$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day7newd = rbind(day7newd, day7newd1)
  
  
  day14newd1 =  data.frame(arm=levels(necfeedsub$arm)) 
  
  day14newd1 = day14newd1 %>% filter(!(arm %in% current)) %>% 
    group_by(arm) %>%
    mutate(
      values = paste0("visitday 14:arm", arm), reference = current,
      pval = sb4j$p.pv[[values]], result=ifelse(is.nan(pval)|pval>0.9999,"p>0.9999",ifelse(pval<0.0001,paste0("p<0.0001"),paste0("p=",format(round(pval,4), nsmall = 4)))))%>%
    dplyr::select(arm, reference, result)
  
  
  day14newd = rbind(day14newd, day14newd1)
  
  
  print(i)
}


newd = rbind(newd,
             newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

newd = unique(as.data.frame(newd))

result_matrix4j <- pivot_wider(newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix4j) %>% autofit() %>% save_as_docx( path = "OUTPUT/prob_infecting_atleast_one_mosq_day2.docx")



day7newd = rbind(day7newd,
                 day7newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day7newd = unique(as.data.frame(day7newd))

result_matrix4j2 <- pivot_wider(day7newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix4j2) %>% autofit() %>% save_as_docx( path = "OUTPUT/prob_infecting_atleast_one_mosq_day7.docx")


day14newd = rbind(day14newd,
                  day14newd0 %>% mutate(reference = arm, result=RRG) %>% select(arm, reference, result)) %>% arrange(reference, arm)

day14newd = unique(as.data.frame(day14newd))

result_matrix4j3 <- pivot_wider(day14newd, names_from = arm, values_from = result)%>% select(reference,`DHA-PPQ`, `SP-AQ`, `PY-AS`, `AS-AQ`, AL, `Non-ACT-PQ`, `ACT-PQ`, `ACT-TQ`) %>% mutate(reference = factor(reference, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ"))) %>% arrange(reference)

flextable(result_matrix4j3) %>% autofit() %>% save_as_docx( path = "OUTPUT/prob_infecting_atleast_one_mosq_day14.docx")







######Calculate numbers to mention in text/legend ######
summary(nec0$arm)
summary(nec0$study)
summary(nec2$study)
summary(nec7$study)
summary(nec0$gampos)
table(nec0$gampos)
summary(nec0$mic_pf_gct_ul[nec0$mic_pf_gct_ul>0]) #Mean 55.0 (among gametocyte carriers)
exp(mean(log(nec0$mic_pf_gct_ul[nec0$mic_pf_gct_ul>0]), na.rm = TRUE)) #Geometric mean 63.7
table(nec0$person_pos)
summary(nec0$mosq_pos_prop[nec0$mosq_pos_prop>0])

#Numbers for Figure 3 legend
nec2_sub = nec2%>%
  filter(!is.na(gam.p) & !is.na(gam.m))
summary(nec2_sub$arm)  

nec7_sub = nec7%>%
  filter(!is.na(gam.p) & !is.na(gam.m))
summary(nec7_sub$arm)  

nec14_sub = nec14%>%
  filter(!is.na(gam.p) & !is.na(gam.m))
summary(nec14_sub$arm)  

#Numbers for Figure 4 legend
nec2_sub = necfeedsub%>%
  filter((days==2) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop))
summary(nec2_sub$arm)  

nec7_sub = necfeedsub%>%
  filter((days==7) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop))
summary(nec7_sub$arm)  

nec14_sub = necfeedsub%>%
  filter((days==14) & !is.na(mosq_total) & mosq_total>0 & !is.na(prop))
summary(nec14_sub$arm)  


#Numbers for Figure 6 legend
#infectivity
nectime5b %>%
  group_by(arm) %>%
  summarise(unique_studyids = n_distinct(studyid_str))      

#microscopy
nectime5c %>%
  group_by(arm) %>%
  summarise(unique_studyids = n_distinct(studyid_str))       

#PCR 
nectime5a %>%
  group_by(arm) %>%
  summarise(unique_studyids = n_distinct(studyid_str)) 

         
nec0%>% 
  summarise(age = mese(age,1),
            temp = mese(temp,1),
            logpar = mediqr(tmpd,1),
            loggam1 = mediqr(mic_pf_gct_ul,1),
            loggam2 = mediqr(qrtpcr_gct_ul,1),
            loggamf = mediqr(qrtpcr_femgct_ul,1),
            proppos = mediqr(mosq_pos_prop,1))

nec0 %>%
  group_by(study) %>%
  summarise(
    Count = n(),                 
    Median = median(mosq_total, na.rm = TRUE),
    Q25 = quantile(mosq_total, 0.25, na.rm = TRUE), 
    Q75 = quantile(mosq_total, 0.75, na.rm = TRUE),
    Min = min(mosq_total, na.rm = TRUE),
    Max = max(mosq_total, na.rm = TRUE)
  )

summary(nec0$gam.p[nec0$study=="PQ03"])
summary(nec0$gam.m[nec0$study=="PQ03"])

table(nec0$gam.p>0&nec0$gam.m==0)
filtered_data <- subset(nec0, gam.p > 0 & gam.m == 0)
study_ids <- filtered_data$studyid_str
filtered_data <- subset(nec, studyid_str == "PQ-03-001" & gam.m > 0)
visit_numbers <- filtered_data$studyvisit_num
filtered_data <- subset(nec, studyid_str == "PQ-03-023" & gam.m > 0)
visit_numbers <- filtered_data$studyvisit_num

#Number of infectious individuals that have gametocytes by microscopy
table(nec0$person_pos==1)
nec0_inf <- nec0[!is.na(nec0$person_pos) & !is.na(nec0$gam.m) & nec0$person_pos == 1,]
sum(nec0_inf$gam.m > 0, na.rm = TRUE)
sum(nec0_inf$gam.m == 0, na.rm = TRUE)

#% of all infectious individuals has a submicroscopic gametocyte infection at days 2, 7, 14
nec2_inf<-nec2[nec2$person_pos>0,]
sum(nec2_inf$gam.m == 0 & nec2_inf$gam.p > 0, na.rm = TRUE)

nec7_inf<-nec7[nec7$person_pos>0,]
sum(nec7_inf$gam.m == 0 & nec7_inf$gam.p > 0, na.rm = TRUE)

nec14_inf<-nec14[nec14$person_pos>0,]
sum(nec14_inf$gam.m == 0 & nec14_inf$gam.p > 0, na.rm = TRUE)
length(na.omit(nec14_inf$studyid_str))

#%individuals became infectious to mosquitoes after treatment while initially not infecting mosquitoes

results<-nec %>%
  filter(studyvisit_num %in% c(0, 1, 2, 3, 5, 7, 10, 14, 21, 28, 35, 42, 49)) %>%
  group_by(studyid_str,arm) %>%
  summarise(
    meets_condition_0 = any(person_pos == 0 & studyvisit_num == 0),
    meets_condition_following = any(person_pos == 1 & studyvisit_num %in% c(1, 2, 3, 5, 7, 10, 14, 21, 28, 35, 42, 49)),
    meets_both = meets_condition_0 & meets_condition_following
  ) %>%
  filter(meets_both) %>%
  select(studyid_str,arm)

n_distinct(results$studyid_str)
summary(results$arm)
summary(nec0$arm)



#% of all gametocyte infections that were submicroscopic

intermediate_results <- nec %>%
  filter(!is.na(gam.p))
denominator <- nrow(intermediate_results)
final_results <- intermediate_results %>%
  filter(gam.p > 0, gam.m == 0)
numerator <- nrow(final_results)




