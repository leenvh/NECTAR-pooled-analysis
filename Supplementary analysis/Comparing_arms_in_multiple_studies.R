
source("data prep and cleaning.R")
source("within arm reductions by study.R")


# A simpler summary -------------------------------------------------------

#### RRG: Relative reduction in gametocytes####

#Similar to how we look at transmission reducing activity where we consider the decline in oocysts. Now we look at the day 2 and day 7 decline in gametocytes (both microscopy and pcr). 

##### Microscopy #####

necsub = nec %>% 
  filter((days==0 | days==2 | days==7| days==14) & !is.na(gam.m)) %>%
  select(studyid_str, study,days,arm,gam.m) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))

necsub$loggam = log(necsub$gam.m+0.1)

mod_main = bam(loggam ~ arm*visit*study +s(studyid_str, bs="re"), data = necsub, discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

studylevs = levels(nec$study)

results.all = results.all %>% mutate(study=factor(study, levels=studylevs),arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(study, arm, visit)

results.all = left_join(necsub %>% select(study, arm, visit) %>% filter(visit != "day 0") %>% unique(),results.all, by=c("study","arm", "visit")) %>% mutate(arm2= interaction(study,arm) )  %>% arrange(study, arm, visit) %>%    group_by(arm, visit) %>%   mutate(n=n()) %>%   filter(n!=1) %>%   ungroup()



ggplot(data=results.all, aes(x=arm2, y=RR/100, fill=study, group=visit))+
  geom_bar(aes(alpha=visit),col="grey90",stat="identity", position = position_dodge(preserve = "single"), width=0.75)+
  scale_alpha_manual(values=c(0.4,0.6,0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_discrete()+
  ylab("Relative reduction (%) in \ngametocyte densities by microscopy")+
  coord_cartesian(ylim = c(-0.5, NA))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels=(results.all %>% filter(visit=="day 2") %>% select(arm, arm2) %>% arrange(arm2))$arm)


ggsave("OUTPUT/NMA_bystudy/Gametocyte_density_mic_bystudy.png",device="png", width=20, height=10, dpi=600)



##### pCR #####


necsub = nec %>% 
  filter((days==0 | days==2 | days==7| days==14) & !is.na(gam.p)) %>%
  select(studyid_str, study,days,arm,gam.p) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))

necsub$loggam = log(necsub$gam.p+0.01)


mod_main = bam(loggam ~ arm*visit*study +s(studyid_str, bs="re"), data = necsub, discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

studylevs = levels(nec$study)

results.all = results.all %>% mutate(study=factor(study, levels=studylevs),arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(study, arm, visit)

results.all = left_join(necsub %>% select(study, arm, visit) %>% filter(visit != "day 0") %>% unique(),results.all, by=c("study","arm", "visit")) %>% mutate(arm2= interaction(study,arm) )  %>% arrange(study, arm, visit) %>%    group_by(arm, visit) %>%   mutate(n=n()) %>%   filter(n!=1) %>%   ungroup()


ggplot(data=results.all, aes(x=arm2, y=RR/100, fill=study, group=visit))+
  geom_bar(aes(alpha=visit),col="grey90",stat="identity", position = position_dodge(preserve = "single"), width=0.75)+
  scale_alpha_manual(values=c(0.4,0.6,0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_discrete()+
  ylab("Relative reduction (%) in \ngametocyte densities by PCR")+
  coord_cartesian(ylim = c(-0.5, NA))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels=(results.all %>% filter(visit=="day 2") %>% select(arm, arm2) %>% arrange(arm2))$arm)


ggsave("OUTPUT/NMA_bystudy/Gametocyte_density_PCR_bystudy.png",device="png", width=20, height=10, dpi=600)


####RRGP: Relative reduction in gametocyte prevalence####

#####RRGP by microscopy #####


necsub = nec %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(gam.m)) %>%
  select(studyid_str,study,days,arm,gam.m) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm))) %>%
  group_by(studyid_str) %>%
  mutate(gam0.m=ifelse(first(days)!=0,NA,first(gam.m))) %>%
  ungroup()

necsub$loggam = as.numeric(I(necsub$gam.m>0))

mod_main = bam(loggam ~ arm*visit*study +s(studyid_str, bs="re"), data = necsub, discrete=T, family=poisson(link="log"), nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

studylevs = levels(nec$study)

results.all = results.all %>% mutate(study=factor(study, levels=studylevs),arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(study, arm, visit)

results.all = left_join(necsub %>% select(study, arm, visit) %>% filter(visit != "day 0") %>% unique(),results.all, by=c("study","arm", "visit")) %>% mutate(arm2= interaction(study,arm) )  %>% arrange(study, arm, visit) %>%    group_by(arm, visit) %>%   mutate(n=n()) %>%   filter(n!=1) %>%   ungroup()


ggplot(data=results.all, aes(x=arm2, y=RR/100, fill=study, group=visit))+
  geom_bar(aes(alpha=visit),col="grey90",stat="identity", position = position_dodge(preserve = "single"), width=0.75)+
  scale_alpha_manual(values=c(0.4,0.6,0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_discrete()+
  ylab("Relative reduction (%) in \ngametocyte prevalence by microscopy")+
  coord_cartesian(ylim = c(-0.5, NA))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels=(results.all %>% filter(visit=="day 2") %>% select(arm, arm2) %>% arrange(arm2))$arm)


ggsave("OUTPUT/NMA_bystudy/Gametocyte_prevalence_mic_bystudy.png",device="png", width=20, height=10, dpi=600)


#####RRGP by PCR#####

necsub = nec %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(gam.p)) %>%
  select(studyid_str,study,days,arm,gam.p) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))%>%
  group_by(studyid_str) %>%
  mutate(gam0.p=ifelse(first(days)!=0,NA,first(gam.p))) %>%
  ungroup()


necsub$loggam = as.numeric(I(necsub$gam.p>0))

mod_main = bam(loggam ~ arm*visit*study +s(studyid_str, bs="re"), data = necsub, discrete=T, family=poisson(link="log"), nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

studylevs = levels(nec$study)

results.all = results.all %>% mutate(study=factor(study, levels=studylevs),arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(study, arm, visit)

results.all = left_join(necsub %>% select(study, arm, visit) %>% filter(visit != "day 0") %>% unique(),results.all, by=c("study","arm", "visit")) %>% mutate(arm2= interaction(study,arm) )  %>% arrange(study, arm, visit) %>%    group_by(arm, visit) %>%   mutate(n=n()) %>%   filter(n!=1) %>%   ungroup()


ggplot(data=results.all, aes(x=arm2, y=RR/100, fill=study, group=visit))+
  geom_bar(aes(alpha=visit),col="grey90",stat="identity", position = position_dodge(preserve = "single"), width=0.75)+
  scale_alpha_manual(values=c(0.4,0.6,0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_discrete()+
  ylab("Relative reduction (%) in \ngametocyte prevalence by PCR")+
  coord_cartesian(ylim = c(-0.5, NA))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels=(results.all %>% filter(visit=="day 2") %>% select(arm, arm2) %>% arrange(arm2))$arm)


ggsave("OUTPUT/NMA_bystudy/Gametocyte_prevalence_PCR_bystudy.png",device="png", width=20, height=10, dpi=600)



#### TBA: transmission blocking activity ####

nec2 = nec %>% 
  group_by(studyid_str) %>%
  arrange(studyid_str, studyvisit_num) %>%
  mutate(adjust = ifelse(is.na(mosq_pos_prop) & studyvisit_num==14  & lag(studyvisit_num)==7 & lag(mosq_pos_prop)==0 & arm!="Non-ACT-PQ", T, F),
         mosq_pos_prop = ifelse(adjust, 0, mosq_pos_prop),
         mosq_total = ifelse(adjust, 5, mosq_total),
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
  rows_to_update <- necfeedsub %>%
    filter(mosq_pos == 0) %>%
    group_by(study,arm, visit) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    mutate(new_prev = 0.5) %>%
    select(row_id, new_prev)
  
  necfeedsub <- necfeedsub %>%
    left_join(rows_to_update, by = "row_id") %>%
    mutate(prev = ifelse(!is.na(new_prev), new_prev, prev)) %>%
    select(-row_id, -new_prev)
  



mod_main = bam(prev ~ arm*visit*study+s(studyid_str, bs="re"), family=poisson(link="log"), data = necfeedsub, discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

studylevs = levels(nec$study)

results.all = results.all %>% mutate(study=factor(study, levels=studylevs),arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(study, arm, visit)

results.all = left_join(necsub %>% select(study, arm, visit) %>% filter(visit != "day 0") %>% unique(),results.all, by=c("study","arm", "visit")) %>% mutate(arm2= interaction(study,arm) )  %>% arrange(study, arm, visit) %>%    group_by(arm, visit) %>%   mutate(n=n()) %>%   filter(n!=1) %>%   ungroup()


results.all = results.all %>%
  mutate(RR = ifelse(visit=="day 14" & arm=="Non-ACT-PQ",NA,RR),
         ci_upper = ifelse(visit=="day 14" & arm=="Non-ACT-PQ",NA,ci_upper),
         ci_lower = ifelse(visit=="day 14" & arm=="Non-ACT-PQ",NA,ci_lower))


ggplot(data=results.all, aes(x=arm2, y=RR/100, fill=study, group=visit))+
  geom_bar(aes(alpha=visit),col="grey90",stat="identity", position = position_dodge(preserve = "single"), width=0.75)+
  scale_alpha_manual(values=c(0.4,0.6,0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_discrete()+
  ylab("Relative reduction (%) in \nproportion infected mosquitos")+
  coord_cartesian(ylim = c(-0.5, NA))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels=(results.all %>% filter(visit=="day 2") %>% select(arm, arm2) %>% arrange(arm2))$arm)


ggsave("OUTPUT/NMA_bystudy/TBA_bystudy.png",device="png", width=20, height=10, dpi=600)


# Oocyst density ----------------------------------------------


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


mod_main = bam(mean_ooc_all ~ arm*visit*study+s(studyid_str, bs="re"), family=poisson(link="log"), weights = mosq_total, data = necfeedsub, discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

studylevs = as.character(sort(unique(necfeedsub$study[!is.na(necfeedsub$mean_ooc_all)])))


results.all = results.all %>% mutate(study=factor(study, levels=studylevs),arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(study, arm, visit)

results.all = left_join(necsub %>% select(study, arm, visit) %>% filter(visit != "day 0") %>% unique(),results.all, by=c("study","arm", "visit")) %>% mutate(arm2= interaction(study,arm) )  %>% arrange(study, arm, visit) %>%    group_by(arm, visit) %>%   mutate(n=n()) %>%   filter(n!=1) %>%   ungroup() %>% na.omit()



ggplot(data=results.all, aes(x=arm2, y=RR/100, fill=study, group=visit))+
  geom_bar(aes(alpha=visit),col="grey90",stat="identity", position = position_dodge(preserve = "single"), width=0.75)+
  scale_alpha_manual(values=c(0.4,0.6,0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_discrete()+
  ylab("Relative reduction (%) in oocysts\nper dissected mosquito")+
  coord_cartesian(ylim = c(-0.5, NA))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels=(results.all %>% filter(visit=="day 2") %>% select(arm, arm2) %>% arrange(arm2))$arm)+
  geom_text(data = data.frame(arm2 = "NECTAR1.ACT-PQ", arm="ACT-PQ",study="NECTAR1",RR=50, visit="day 14"), aes(x=arm2, y=RR/100,  group=visit,col=study,label="No data"), angle=90, position=position_dodge(preserve="single", width=-0.5), show.legend = FALSE)


ggsave("OUTPUT/NMA_bystudy/TRA_bystudy.png",device="png", width=20, height=10, dpi=600)





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

mod_main = bam(prev ~ arm*visit*study, family=poisson(link="log"), data = necfeedsub, discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

studylevs = levels(nec$study)

results.all = results.all %>% mutate(study=factor(study, levels=studylevs),arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(study, arm, visit)

results.all = left_join(necsub %>% select(study, arm, visit) %>% filter(visit != "day 0") %>% unique(),results.all, by=c("study","arm", "visit")) %>% mutate(arm2= interaction(study,arm) )  %>% arrange(study, arm, visit) %>%    group_by(arm, visit) %>%   mutate(n=n()) %>%   filter(n!=1) %>%   ungroup()


results.all = results.all %>%
  mutate(RR = ifelse(visit=="day 14" & arm=="Non-ACT-PQ",NA,RR),
         ci_upper = ifelse(visit=="day 14" & arm=="Non-ACT-PQ",NA,ci_upper),
         ci_lower = ifelse(visit=="day 14" & arm=="Non-ACT-PQ",NA,ci_lower))


ggplot(data=results.all, aes(x=arm2, y=RR/100, fill=study, group=visit))+
  geom_bar(aes(alpha=visit),col="grey90",stat="identity", position = position_dodge(preserve = "single"), width=0.75)+
  scale_alpha_manual(values=c(0.4,0.6,0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5))+
  scale_fill_discrete()+
  ylab("Relative reduction (%) in \nindividual's infectious likelihood")+
  coord_cartesian(ylim = c(-0.5, NA))+
  scale_y_continuous(labels = scales::percent)+
  scale_x_discrete(labels=(results.all %>% filter(visit=="day 2") %>% select(arm, arm2) %>% arrange(arm2))$arm)

ggsave("OUTPUT/NMA_bystudy/infectiousness_bystudy.png",device="png", width=20, height=10, dpi=600)
