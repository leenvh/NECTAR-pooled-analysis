
source("data prep and cleaning.R")
source("within arm reductions.R")
source("studywise contrasts.R")


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

nectime5a_summary = nectime5a %>% group_by(study,arm) %>% summarise(n=n());
# Aggregate number of patients per treatment arm across all studies
patient_counts <- nectime5a_summary %>%
  group_by(arm) %>%
  summarise(patients = sum(n)) %>%
  ungroup() %>%
  mutate(arm = as.character(arm))


# Count number of studies comparing each treatment pair
study_df <- results %>%
  mutate(treat1 = as.character(arm), treat2 = as.character(ref)) %>%
  select(treat1, treat2) %>%
  mutate(pair = ifelse(as.character(treat1) < as.character(treat2),
                       paste(treat1, treat2, sep = " vs "),
                       paste(treat2, treat1, sep = " vs "))) %>%
  group_by(pair, treat1, treat2) %>%
  summarise(comparisons = n()) %>%
  ungroup()

patient_counts <- patient_counts %>%
  rename(name = arm)

study_df_fixed <- study_df %>%
  mutate(t1 = pmin(treat1, treat2),
         t2 = pmax(treat1, treat2)) %>%
  group_by(t1, t2) %>%
  summarise(comparisons = sum(comparisons), .groups = "drop") %>%
  rename(treat1 = t1, treat2 = t2)


library(igraph)
library(ggraph)
# Create igraph object
g <- graph_from_data_frame(study_df_fixed, directed = FALSE, vertices = patient_counts)

ggraph(g, layout = 'fr') +
  geom_edge_link(aes(width = comparisons), color = "black", alpha = 0.7) +
  geom_node_point(aes(size = patients), color = "#2b6777", alpha = 0.9) +
  geom_node_text(aes(label = name), vjust = 1.5, size = 5) +
  scale_size_continuous(range = c(3, 12), name = "Number of Patients") +
  scale_edge_width(range = c(0.5, 4),breaks=c(1,2), name = "Number of Comparisons") +
  theme_void()+
  theme(legend.position = "right",
        plot.margin = unit(c(1.5, 1.5, 1.5, 1.5), "cm"))  


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


flextable(result_matrix5a) %>% autofit()  %>% save_as_docx( path = "OUTPUT/NMA_updated/survival_gampcr_NMA.docx")

check.surv(nma)
ggsave("OUTPUT/NMA_updated/CHECK_survival_gampcr_NMA.png", dpi=600, device="png", width=15, height=12, units="cm")


my_base_plot <- function() forest(nma)
# Run and capture as a grob
my_base_plot()
grob_object1 <- grid.grab()

ggarrange(grob_object1)
ggsave("OUTPUT/NMA_updated/Forest_survival_gampcr_NMA.png", dpi=600, device="png", width=15, height=12, units="cm")


ranks <- netrank(nma, small.values = "bad")
rdata <- data.frame(arm=names(ranks$ranking.random), pscore =ranks$ranking.random) %>%
  arrange(pscore) 
rdata$arm = factor(rdata$arm, levels = rdata$arm)

pscoreplot = ggplot(data=rdata, aes(x="",y=arm, fill=pscore))+
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

pscoreplot

ggsave("OUTPUT/NMA_updated/pscore_survival_gampcr_NMA.png", dpi=600, device="png", width=12, height=15, units="cm")


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


flextable(result_matrix5b) %>% autofit()  %>% save_as_docx( path = "OUTPUT/NMA_updated/survival_infectivity_adjpcr_NMA.docx")

check.surv(nma)
ggsave("OUTPUT/NMA_updated/CHECK_survival_infectivity_adjpcr_NMA.png", dpi=600, device="png", width=15, height=12, units="cm")



my_base_plot <- function() forest(nma)
# Run and capture as a grob
my_base_plot()
grob_object1 <- grid.grab()

ggarrange(grob_object1)
ggsave("OUTPUT/NMA_updated/Forest_survival_infectivity_adjpcr_NMA.png", dpi=600, device="png", width=15, height=12, units="cm")


ranks <- netrank(nma, small.values = "bad")
rdata <- data.frame(arm=names(ranks$ranking.random), pscore =ranks$ranking.random) %>%
  arrange(pscore) 
rdata$arm = factor(rdata$arm, levels = rdata$arm)

pscoreplot = ggplot(data=rdata, aes(x="",y=arm, fill=pscore))+
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

pscoreplot

ggsave("OUTPUT/NMA_updated/pscore_survival_infectivity_adjpcr_NMA.png", dpi=600, device="png", width=12, height=15, units="cm")


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

flextable(result_matrix5c) %>% autofit()  %>% save_as_docx( path = "OUTPUT/NMA_updated/survival_gammic_NMA.docx")

check.surv(nma)
ggsave("OUTPUT/NMA_updated/CHECK_survival_gammic_NMA.png", dpi=600, device="png", width=15, height=12, units="cm")



my_base_plot <- function() forest(nma)
# Run and capture as a grob
my_base_plot()
grob_object1 <- grid.grab()

ggarrange(grob_object1)
ggsave("OUTPUT/NMA_updated/Forest_survival_gammic_NMA.png", dpi=600, device="png", width=15, height=12, units="cm")


ranks <- netrank(nma, small.values = "bad")
rdata <- data.frame(arm=names(ranks$ranking.random), pscore =ranks$ranking.random) %>%
  arrange(pscore) 
rdata$arm = factor(rdata$arm, levels = rdata$arm)

pscoreplot = ggplot(data=rdata, aes(x="",y=arm, fill=pscore))+
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

pscoreplot

ggsave("OUTPUT/NMA_updated/pscore_survival_gammic_NMA.png", dpi=600, device="png", width=12, height=15, units="cm")

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

flextable(result_matrix5d) %>% autofit()  %>% save_as_docx( path = "OUTPUT/NMA_updated/survival_infectivity_adjmic_NMA.docx")


check.surv(nma)
ggsave("OUTPUT/NMA_updated/CHECK_survival_infectivity_adjmic_NMA.png", dpi=600, device="png", width=15, height=12, units="cm")


my_base_plot <- function() forest(nma)
# Run and capture as a grob
my_base_plot()
grob_object1 <- grid.grab()

ggarrange(grob_object1)
ggsave("OUTPUT/NMA_updated/Forest_survival_infectivity_adjmic_NMA.png", dpi=600, device="png", width=15, height=12, units="cm")


ranks <- netrank(nma, small.values = "bad")
rdata <- data.frame(arm=names(ranks$ranking.random), pscore =ranks$ranking.random) %>%
  arrange(pscore) 
rdata$arm = factor(rdata$arm, levels = rdata$arm)

pscoreplot = ggplot(data=rdata, aes(x="",y=arm, fill=pscore))+
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

pscoreplot

ggsave("OUTPUT/NMA_updated/pscore_survival_infectivity_adjmic_NMA.png", dpi=600, device="png", width=12, height=15, units="cm")


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

results.all = results.all %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(arm, visit)


ggplot(data=results.all, aes(x=arm, y=RR/100, fill=visit, group=visit))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_point(data=medgams, aes(x=arm, y=RR/100, group=visit), col="red", position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  #geom_errorbar(data=medgams,aes(x=arm,ymin=RR_lower/100, ymax=RR_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="red",alpha=0.2)+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \ngametocyte densities by microscopy")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA_updated/Gametocyte_density_mic_NMA.png",device="png", width=10, height=7, dpi=600)


arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")


mod_main2 = bam(loggam ~ arm*visit*study + s(studyid_str, bs="re"), data = necsub, discrete=T, nthread=4)
summary(mod_main2)

result_day2 <- get_studywise_contrasts_bam(mod_main2, necsub, "day 2", arms)
flextable(result_day2[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/Gam_relative_reduction_mic_day2_NMA.docx")
result_day7 <- get_studywise_contrasts_bam(mod_main2, necsub, "day 7", arms)
flextable(result_day7[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/Gam_relative_reduction_mic_day7_NMA.docx")
result_day14 <- get_studywise_contrasts_bam(mod_main2, necsub, "day 14", arms)
flextable(result_day14[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/Gam_relative_reduction_mic_day14_NMA.docx")


comparisons = result_day2[[1]]$data$Comparison

ggarrange(result_day2[[1]]+
            ggtitle("Day 2")+
            theme(axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_x_discrete(limits = comparisons),
          result_day7[[1]]+
                      theme(axis.text.y = element_blank(),
                            axis.title = element_blank(),
                            panel.border = element_rect(color = "black", fill = NA, size = 0.5)
                      )+
                      ggtitle("Day 7")+
            scale_x_discrete(limits = comparisons),
          result_day14[[1]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            ggtitle("Day 14")+
            scale_x_discrete(limits = comparisons),
          legend="top",
          common.legend = TRUE,
          nrow=1, ncol=3,
          widths = c(1.4, 1, 1)
          )
ggsave("OUTPUT/NMA_updated/CHECK1_Gametocyte_density_mic_NMA.png", dpi=600, device="png", width=28, height=24, units="cm")


png("OUTPUT/NMA_updated/FOREST_Gametocyte_density_mic_NMA.png", width = 18, height = 24, units = "cm", res = 1200)
grid.newpage()
gridExtra::grid.arrange(result_day2[[4]],result_day7[[4]],result_day14[[4]])
dev.off()

# Extract y-axis levels from plot 1 (to use as standard)
p1_data <- result_day2[[6]]
common_levels <- p1_data$arm


ggarrange(result_day2[[3]]+
            theme(axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            ),
          result_day7[[3]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels),
          result_day14[[3]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels),
          legend="top",
          common.legend = TRUE,
          nrow=1, ncol=3,
          widths = c(1.4, 1, 1)
)
ggsave("OUTPUT/NMA_updated/pscore_Gametocyte_density_mic_NMA.png", dpi=600, device="png", width=28, height=14, units="cm")


result_day2[[5]]+theme(title = element_blank())  + coord_cartesian(clip = "off")

ggsave("OUTPUT/NMA_updated/network_Gametocyte_density_mic_NMA.png", dpi=600, device="png", width=14, height=22, units="cm")


##### pCR #####


necsub = nec %>% 
  filter((days==0 | days==2 | days==7| days==14) & !is.na(gam.p)) %>%
  select(studyid_str, study,days,arm,gam.p) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))

necsub$loggam = log(necsub$gam.p+0.01)

medgams = left_join(rbind(
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 2") %>% summarise(med=exp(mean(loggam,na.rm=T))),
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 7") %>% summarise(med=exp(mean(loggam,na.rm=T))),
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 14") %>% summarise(med=exp(mean(loggam,na.rm=T)))), necsub %>% group_by(arm) %>% filter(visit=="day 0") %>% summarise(med0=exp(mean(loggam,na.rm=T))), by="arm") %>% mutate(RR=100-(med/med0)*100)

mod_main = bam(loggam ~ arm*visit + s(studyid_str, bs="re"), data = necsub, discrete=T, nthread=4)
summary(mod_main)



results.all = get_reductions(mod_main)

results.all = results.all %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(arm, visit)


ggplot(data=results.all, aes(x=arm, y=RR/100, fill=visit, group=visit))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_point(data=medgams, aes(x=arm, y=RR/100, group=visit), col="red", position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \ngametocyte densities by pcr")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA_updated/Gametocyte_density_pcr_NMA.png",device="png", width=10, height=7, dpi=600)


arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

mod_main2 = bam(loggam ~ arm*visit*study + s(studyid_str, bs="re"), data = necsub, discrete=T, nthread=4)
summary(mod_main2)

result_day2 <- get_studywise_contrasts_bam(mod_main2, necsub, "day 2", arms)
flextable(result_day2[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/Gam_relative_reduction_pcr_day2_NMA.docx")
result_day7 <- get_studywise_contrasts_bam(mod_main2, necsub, "day 7", arms)
flextable(result_day7[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/Gam_relative_reduction_pcr_day7_NMA.docx")
result_day14 <- get_studywise_contrasts_bam(mod_main2, necsub, "day 14", arms)
flextable(result_day14[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/Gam_relative_reduction_pcr_day14_NMA.docx")



ggarrange(result_day2[[1]]+
            ggtitle("Day 2")+
            theme(axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_x_discrete(limits = comparisons),
          result_day7[[1]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            ggtitle("Day 7")+
            scale_x_discrete(limits = comparisons),
          result_day14[[1]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            ggtitle("Day 14")+
            scale_x_discrete(limits = comparisons),
          legend="top",
          common.legend = TRUE,
          nrow=1, ncol=3,
          widths = c(1.4, 1, 1)
)
ggsave("OUTPUT/NMA_updated/CHECK1_Gametocyte_density_pcr_NMA.png", dpi=600, device="png", width=28, height=24, units="cm")


png("OUTPUT/NMA_updated/FOREST_Gametocyte_density_pcr_NMA.png", width = 18, height = 24, units = "cm", res = 1200)
grid.newpage()
gridExtra::grid.arrange(result_day2[[4]],result_day7[[4]],result_day14[[4]])
dev.off()

# Extract y-axis levels from plot 1 (to use as standard)
p1_data <- result_day2[[6]]
common_levels <- p1_data$arm


ggarrange(result_day2[[3]]+
            theme(axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            ),
          result_day7[[3]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels),
          result_day14[[3]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels),
          legend="top",
          common.legend = TRUE,
          nrow=1, ncol=3,
          widths = c(1.4, 1, 1)
)
ggsave("OUTPUT/NMA_updated/pscore_Gametocyte_density_pcr_NMA.png", dpi=600, device="png", width=28, height=14, units="cm")


result_day7[[5]]+theme(title = element_blank())  + coord_cartesian(clip = "off")

ggsave("OUTPUT/NMA_updated/network_Gametocyte_density_pcr_NMA.png", dpi=600, device="png", width=14, height=22, units="cm")



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

medgams = left_join(rbind(
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 2") %>% summarise(med=mean(loggam,na.rm=T)),
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 7") %>% summarise(med=mean(loggam,na.rm=T)),
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 14") %>% summarise(med=mean(loggam,na.rm=T))), necsub %>% group_by(arm) %>% filter(visit=="day 0") %>% summarise(med0=mean(loggam,na.rm=T)), by="arm") %>% mutate(RR=100-(med/med0)*100)

mod_main = bam(loggam ~ arm*visit+s(studyid_str, bs="re"), family=poisson(link="log"), data = necsub, discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

results.all = results.all %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(arm, visit)


ggplot(data=results.all, aes(x=arm, y=RR/100, fill=visit, group=visit))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_point(data=medgams, aes(x=arm, y=RR/100, group=visit), col="red", position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  #geom_errorbar(data=medgams,aes(x=arm,ymin=RR_lower/100, ymax=RR_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="red",alpha=0.2)+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \ngametocyte prevalence by microscopy")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA_updated/Gametocyte_prevalence_mic_NMA.png",device="png", width=10, height=7, dpi=600)


arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

mod_main2 = bam(loggam ~ arm*visit*study + s(studyid_str, bs="re"), family=poisson(link="log"), data = necsub, discrete=T, nthread=4)
summary(mod_main2)

result_day2 <- get_studywise_contrasts_bam(mod_main2, necsub, "day 2", arms)
flextable(result_day2[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/Gam_relative_reduction_prev_mic_day2_NMA.docx")
result_day7 <- get_studywise_contrasts_bam(mod_main2, necsub, "day 7", arms)
flextable(result_day7[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/Gam_relative_reduction_prev_mic_day7_NMA.docx")
result_day14 <- get_studywise_contrasts_bam(mod_main2, necsub, "day 14", arms)
flextable(result_day14[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/Gam_relative_reduction_prev_mic_day14_NMA.docx")



ggarrange(result_day2[[1]]+
            ggtitle("Day 2")+
            theme(axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_x_discrete(limits = comparisons),
          result_day7[[1]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            ggtitle("Day 7")+
            scale_x_discrete(limits = comparisons),
          result_day14[[1]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            ggtitle("Day 14")+
            scale_x_discrete(limits = comparisons),
          legend="top",
          common.legend = TRUE,
          nrow=1, ncol=3,
          widths = c(1.4, 1, 1)
)
ggsave("OUTPUT/NMA_updated/CHECK1_Gametocyte_prevalence_mic_NMA.png", dpi=600, device="png", width=28, height=24, units="cm")



png("OUTPUT/NMA_updated/FOREST_Gametocyte_prevalence_mic_NMA.png", width = 18, height = 24, units = "cm", res = 1200)
grid.newpage()
gridExtra::grid.arrange(result_day2[[4]],result_day7[[4]],result_day14[[4]])
dev.off()

# Extract y-axis levels from plot 1 (to use as standard)
p1_data <- result_day2[[6]]
common_levels <- p1_data$arm


ggarrange(result_day2[[3]]+
            theme(axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            ),
          result_day7[[3]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels),
          result_day14[[3]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels),
          legend="top",
          common.legend = TRUE,
          nrow=1, ncol=3,
          widths = c(1.4, 1, 1)
)
ggsave("OUTPUT/NMA_updated/pscore_Gametocyte_prevalence_mic_NMA.png", dpi=600, device="png", width=28, height=14, units="cm")


result_day2[[5]]+theme(title = element_blank())  + coord_cartesian(clip = "off")

ggsave("OUTPUT/NMA_updated/network_Gametocyte_prevalence_mic_NMA.png", dpi=600, device="png", width=14, height=22, units="cm")



#####RRGP by PCR#####

necsub = nec %>% 
  filter((days==0 | days==2 | days==7 | days==14) & !is.na(gam.p)) %>%
  select(studyid_str,study,days,arm,gam.p) %>%
  mutate(visit= factor(days, labels=c("day 0", "day 2","day 7","day 14")), arm=as.factor(as.character(arm)))%>%
  group_by(studyid_str) %>%
  mutate(gam0.p=ifelse(first(days)!=0,NA,first(gam.p))) %>%
  ungroup()


necsub$loggam = as.numeric(I(necsub$gam.p>0))

medgams = left_join(rbind(
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 2") %>% summarise(med=mean(loggam,na.rm=T)),
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 7") %>% summarise(med=mean(loggam,na.rm=T)),
  necsub %>% group_by(arm, visit) %>% filter(visit=="day 14") %>% summarise(med=mean(loggam,na.rm=T))), necsub %>% group_by(arm) %>% filter(visit=="day 0") %>% summarise(med0=mean(loggam,na.rm=T)), by="arm") %>% mutate(RR=100-(med/med0)*100)

mod_main = bam(loggam ~ arm*visit+s(studyid_str, bs="re"), family=poisson(link="log"), data = necsub, discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

results.all = results.all %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(arm, visit)


ggplot(data=results.all, aes(x=arm, y=RR/100, fill=visit, group=visit))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  geom_point(data=medgams, aes(x=arm, y=RR/100, group=visit), col="red", position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  #geom_errorbar(data=medgams,aes(x=arm,ymin=RR_lower/100, ymax=RR_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="red",alpha=0.2)+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \ngametocyte prevalence by pcr")+
  coord_cartesian(ylim = c(-0.5, NA))


ggsave("OUTPUT/NMA_updated/Gametocyte_prevalence_pcr_NMA.png",device="png", width=10, height=7, dpi=600)


arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

mod_main2 = bam(loggam ~ arm*visit*study + s(studyid_str, bs="re"), data = necsub, family=poisson(link="log"), discrete=T, nthread=4)
summary(mod_main2)

result_day2 <- get_studywise_contrasts_bam(mod_main2, necsub, "day 2", arms)
flextable(result_day2[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/Gam_relative_reduction_prev_pcr_day2_NMA.docx")
result_day7 <- get_studywise_contrasts_bam(mod_main2, necsub, "day 7", arms)
flextable(result_day7[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/Gam_relative_reduction_prev_pcr_day7_NMA.docx")
result_day14 <- get_studywise_contrasts_bam(mod_main2, necsub, "day 14", arms)
flextable(result_day14[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/Gam_relative_reduction_prev_pcr_day14_NMA.docx")




ggarrange(result_day2[[1]]+
            ggtitle("Day 2")+
            theme(axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_x_discrete(limits = comparisons),
          result_day7[[1]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            ggtitle("Day 7")+
            scale_x_discrete(limits = comparisons),
          result_day14[[1]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            ggtitle("Day 14")+
            scale_x_discrete(limits = comparisons),
          legend="top",
          common.legend = TRUE,
          nrow=1, ncol=3,
          widths = c(1.4, 1, 1)
)
ggsave("OUTPUT/NMA_updated/CHECK1_Gametocyte_prevalence_pcr_NMA.png", dpi=600, device="png", width=28, height=24, units="cm")



png("OUTPUT/NMA_updated/FOREST_Gametocyte_prevalence_pcr_NMA.png", width = 18, height = 24, units = "cm", res = 1200)
grid.newpage()
gridExtra::grid.arrange(result_day2[[4]],result_day7[[4]],result_day14[[4]])
dev.off()

# Extract y-axis levels from plot 1 (to use as standard)
p1_data <- result_day2[[6]]
common_levels <- p1_data$arm


ggarrange(result_day2[[3]]+
            theme(axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            ),
          result_day7[[3]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels),
          result_day14[[3]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels),
          legend="top",
          common.legend = TRUE,
          nrow=1, ncol=3,
          widths = c(1.4, 1, 1)
)
ggsave("OUTPUT/NMA_updated/pscore_Gametocyte_prevalence_pcr_NMA.png", dpi=600, device="png", width=28, height=14, units="cm")


result_day2[[5]]+theme(title = element_blank())  + coord_cartesian(clip = "off")

ggsave("OUTPUT/NMA_updated/network_Gametocyte_prevalence_pcr_NMA.png", dpi=600, device="png", width=14, height=22, units="cm")





#### TBA: transmission blocking activity ####

nec = nec %>% 
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

nec$ageb = cut(nec$age,breaks=c(-1,5,15.99999999,100), labels=c("<5yrs","5-15yrs","16+yrs"))
nec$ageb = factor(nec$ageb, levels=c("16+yrs","5-15yrs","<5yrs"))


necfeed =nec %>% filter(!is.na(mosq_total) & !(mosq_total==0)) 
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

medgams = left_join(rbind(
  necfeedsub %>% group_by(arm, visit) %>% filter(visit=="day 2") %>% summarise(med=mean(loggam,na.rm=T)),
  necfeedsub %>% group_by(arm, visit) %>% filter(visit=="day 7") %>% summarise(med=mean(loggam,na.rm=T)),
  necfeedsub %>% group_by(arm, visit) %>% filter(visit=="day 14") %>% summarise(med=mean(loggam,na.rm=T))), necfeedsub %>% group_by(arm) %>% filter(visit=="day 0") %>% summarise(med0=mean(loggam,na.rm=T)), by="arm") %>% mutate(RR=100-(med/med0)*100)

medgams[24,]=NA
medgams$arm[24]="Non-ACT-PQ"
medgams$RR[24]=100
medgams$visit[24]="day 14"


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
  



mod_main = bam(prev ~ arm*visit+s(studyid_str, bs="re"), family=poisson(link="log"), data = necfeedsub, discrete=T, nthread=4)

summary(mod_main)

results.all = get_reductions(mod_main)

results.all = results.all %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(arm, visit)

results.all = results.all %>%
  mutate(RR = ifelse(visit=="day 14" & arm=="Non-ACT-PQ",NA,RR),
         ci_upper = ifelse(visit=="day 14" & arm=="Non-ACT-PQ",NA,ci_upper),
         ci_lower = ifelse(visit=="day 14" & arm=="Non-ACT-PQ",NA,ci_lower))

ggplot(data=results.all, aes(x=arm, y=RR/100, fill=visit, group=visit))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  #geom_point(data=medgams, aes(x=arm, y=RR/100, group=visit), col="red", position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  #geom_errorbar(data=medgams,aes(x=arm,ymin=RR_lower/100, ymax=RR_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="red",alpha=0.2)+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \nproportion infected mosquitos")+
  coord_cartesian(ylim = c(-0.5, NA))+
  geom_text(data=data.frame(arm="Non-ACT-PQ", day="Day 14", TRA=0),
            aes(x=arm, y=0.5, label="No Data", group=day), # Use the same group aesthetic
            position=position_dodge(preserve="single", width=-0.5), # Align with bar dodge
            angle=90, size=5, color="black", inherit.aes=FALSE)


ggsave("OUTPUT/NMA_updated/TBA_NMA.png",device="png", width=10, height=7, dpi=600)


arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

mod_main2 = bam(prev ~ arm*visit*study + s(studyid_str, bs="re"), data = necfeedsub, family=poisson(link="log"), discrete=T, nthread=4)
summary(mod_main2)

result_day2 <- get_studywise_contrasts_bam(mod_main2, necfeed, "day 2", arms)
flextable(result_day2[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/TBA_day2_NMA.docx")
result_day7 <- get_studywise_contrasts_bam(mod_main2, necfeed, "day 7", arms)
flextable(result_day7[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/TBA_day7_NMA.docx")
arms_up = arms[which(!(arms %in% "Non-ACT-PQ"))]
result_day14 <- get_studywise_contrasts_bam(mod_main2, necfeed, "day 14", arms_up)
flextable(result_day14[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/TBA_day14_NMA.docx")


df_day2 <- data.frame(result_day2[[2]])
df_day7 <- data.frame(result_day7[[2]])

process_data <- function(df) {
  treatment_levels <- df$Reference
  
  df <- df %>% select(-Reference)
  
  df[upper.tri(df,diag=T)] = NA
  
  colnames(df) <- treatment_levels
  df$treatment <- treatment_levels
  
  df_long <- df %>%
    pivot_longer(-treatment, names_to = "comparison", values_to = "text") %>%
    filter(treatment != comparison) %>%
    mutate(
      activity = as.numeric(str_extract(text, "^[-]?[0-9.]+")),
      abs_diff = abs(activity),
      label_wrapped = str_replace_all(text, ", ", ",<br>") %>%
        str_replace(" \\(", "<br>(")
    ) %>%
    drop_na(activity) %>%
    mutate(
      treatment = factor(treatment, levels = treatment_levels),
      comparison = factor(comparison, levels = treatment_levels)
    )
  
  return(df_long)
}


df2_long <- process_data(df_day2)
df7_long <- process_data(df_day7)


diamp2 = ggplot(df2_long, aes(x = comparison, y = treatment, fill = abs_diff)) +
  geom_tile(color = "white", size = 0.3) +
  geom_richtext(aes(label = label_wrapped),
                size = 2.5,
                fill = NA, label.color = NA,
                angle = 0, hjust = 0.5, vjust = 0.5,
                lineheight = 1) +
  scale_fill_gradient2(
    low = "#d9e4e6", high = "#2b6777",
    limits = c(0, 135),
    breaks = c(0, 20, 40, 60, 80, 100, 120),
    labels = c("0%", "20%", "40%", "60%", "80%", "100%", "120%"),
    name = "Absolute difference\nin relative reduction"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = 10),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  scale_x_discrete(position = "top")
diamp2



diamp7 = ggplot(df7_long, aes(x = comparison, y = treatment, fill = abs_diff)) +
  geom_tile(color = "white", size = 0.3) +
  geom_richtext(aes(label = label_wrapped),
                size = 2.5,
                fill = NA, label.color = NA,
                angle = 0, hjust = 0.5, vjust = 0.5,
                lineheight = 1) +
  scale_fill_gradient2(
    low = "#d9e4e6", high = "#2b6777",
    limits = c(0, 135),
    breaks = c(0, 20, 40, 60, 80, 100, 120),
    labels = c("0%", "20%", "40%", "60%", "80%", "100%", "120%"),
    name = "Absolute difference\nin relative reduction"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = 10),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  scale_x_discrete(position = "top")
diamp7



ggarrange(result_day2[[1]]+
            ggtitle("Day 2")+
            theme(axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_x_discrete(limits = comparisons),
          result_day7[[1]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            ggtitle("Day 7")+
            scale_x_discrete(limits = comparisons),
          result_day14[[1]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            ggtitle("Day 14")+
            scale_x_discrete(limits = comparisons),
          legend="top",
          common.legend = TRUE,
          nrow=1, ncol=3,
          widths = c(1.4, 1, 1)
)
ggsave("OUTPUT/NMA_updated/CHECK1_TBA_NMA.png", dpi=600, device="png", width=28, height=24, units="cm")


png("OUTPUT/NMA_updated/FOREST_TBA_NMA.png", width = 18, height = 24, units = "cm", res = 1200)
grid.newpage()
gridExtra::grid.arrange(result_day2[[4]],result_day7[[4]],result_day14[[4]])
dev.off()

# Extract y-axis levels from plot 1 (to use as standard)
p1_data <- result_day2[[6]]
common_levels <- p1_data$arm


pscore_full = ggarrange(result_day2[[3]]+
            theme(axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels)+
              theme(plot.background = element_blank()),
          result_day7[[3]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels)+
            theme(plot.background = element_blank()),
          result_day14[[3]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels)+
            theme(plot.background = element_blank()),
          legend="top",
          common.legend = TRUE,
          nrow=1, ncol=3,
          widths = c(1.6, 1, 1)
)

pscore_full 

ggsave("OUTPUT/NMA_updated/pscore_TBA_NMA.png", dpi=600, device="png", width=28, height=14, units="cm")


netplot = result_day2[[5]]+theme(title = element_blank())  + coord_cartesian(clip = "off") + theme(legend.position = "top")

netplot

ggsave("OUTPUT/NMA_updated/network_TBA_NMA.png", dpi=600, device="png", width=14, height=22, units="cm")

ggarrange(
ggarrange(netplot, pscore_full, nrow=1, ncol=2, widths = c(0.4, 0.6), labels=c("A","B")),
ggarrange(diamp2+theme(legend.position = "bottom", legend.key.width = unit(2,"cm")), diamp7, nrow=1, ncol=2, common.legend = T, legend="bottom", labels=c("C","D")),
nrow=2, ncol=1, heights = c(0.3, 0.7))

ggsave("OUTPUT/NMA_updated/COMBINED_TBA_NMA.png", dpi=600, device="png", width=28, height=26, units="cm")

# Oocyst density ----------------------------------------------------------


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

results.all = results.all %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(arm, visit)

results.all[22:24,]=NA
results.all$arm[22:24]="Non-ACT-PQ"
results.all$RR[22:24]=NA
results.all$ci_lower[22:24]=NA
results.all$ci_upper[22:24]=NA
results.all$visit[22:24]=c("day 2","day 7","day 14")


ggplot(data=results.all, aes(x=arm, y=RR/100, fill=visit, group=visit))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  #geom_point(data=medgams, aes(x=arm, y=RR/100, group=visit), col="red", position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  #geom_errorbar(data=medgams,aes(x=arm,ymin=RR_lower/100, ymax=RR_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="red",alpha=0.2)+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in oocysts\nper dissected mosquito")+
  coord_cartesian(ylim = c(-0.5, NA))+
  geom_text(data=data.frame(arm="Non-ACT-PQ", visit="day 7", RR=0),
            aes(x=arm, y=0.5, label="No Data", group=visit), # Use the same group aesthetic
            position=position_dodge(preserve="single", width=0), # Align with bar dodge
            angle=90, size=5, color="black", inherit.aes=FALSE)


ggsave("OUTPUT/NMA_updated/TRA_NMA.png",device="png", width=10, height=7, dpi=600)


arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "ACT-PQ", "ACT-TQ")

mod_main2 = bam(mean_ooc_all  ~ arm*visit*study + s(studyid_str, bs="re"), weights=mosq_total, data = necfeedsub, family=poisson(link="log"), discrete=T, nthread=4)
summary(mod_main2)

result_day2 <- get_studywise_contrasts_bam(mod_main2, necfeedsub, "day 2", arms)
flextable(result_day2[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/TRA_day2_NMA.docx")
result_day7 <- get_studywise_contrasts_bam(mod_main2, necfeedsub, "day 7", arms)
flextable(result_day7[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/TRA_day7_NMA.docx")
result_day14 <- get_studywise_contrasts_bam(mod_main2, necfeedsub, "day 14", arms_up)
flextable(result_day14[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/TRA_day14_NMA.docx")




ggarrange(result_day2[[1]]+
            ggtitle("Day 2")+
            theme(axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_x_discrete(limits = comparisons),
          result_day7[[1]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            ggtitle("Day 7")+
            scale_x_discrete(limits = comparisons),
          result_day14[[1]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5))+
            ggtitle("Day 14")+
            scale_x_discrete(limits = comparisons),
          legend="top",
          common.legend = TRUE,
          nrow=1, ncol=3,
          widths = c(1.4, 1, 1)
)
ggsave("OUTPUT/NMA_updated/CHECK1_TRA_NMA.png", dpi=600, device="png", width=28, height=24, units="cm")


png("OUTPUT/NMA_updated/FOREST_TRA_NMA.png", width = 18, height = 24, units = "cm", res = 1200)
grid.newpage()
gridExtra::grid.arrange(result_day2[[4]],result_day7[[4]],result_day14[[4]])
dev.off()

# Extract y-axis levels from plot 1 (to use as standard)
p1_data <- result_day2[[6]]
common_levels <- p1_data$arm


ggarrange(result_day2[[3]]+
            theme(axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            ),
          result_day7[[3]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels),
          result_day14[[3]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels),
          legend="top",
          common.legend = TRUE,
          nrow=1, ncol=3,
          widths = c(1.4, 1, 1)
)
ggsave("OUTPUT/NMA_updated/pscore_TRA_NMA.png", dpi=600, device="png", width=28, height=14, units="cm")


result_day2[[5]]+theme(title = element_blank())  + coord_cartesian(clip = "off")

ggsave("OUTPUT/NMA_updated/network_TRA_NMA.png", dpi=600, device="png", width=14, height=22, units="cm")



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


medgams = left_join(rbind(
  necfeedsub %>% group_by(arm, visit) %>% filter(visit=="day 2") %>% summarise(med=mean(loggam,na.rm=T)),
  necfeedsub %>% group_by(arm, visit) %>% filter(visit=="day 7") %>% summarise(med=mean(loggam,na.rm=T)),
  necfeedsub %>% group_by(arm, visit) %>% filter(visit=="day 14") %>% summarise(med=mean(loggam,na.rm=T))), necfeedsub %>% group_by(arm) %>% filter(visit=="day 0") %>% summarise(med0=mean(loggam,na.rm=T)), by="arm") %>% mutate(RR=100-(med/med0)*100)

medgams[24,]=NA
medgams$arm[24]="Non-ACT-PQ"
medgams$RR[24]=100
medgams$visit[24]="day 14"


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

results.all = results.all %>% mutate(arm = factor(arm, levels=c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")), visit = factor(visit, levels=c("day 2", "day 7", "day 14"))) %>% arrange(arm, visit)

results.all = results.all %>%
  mutate(RR = ifelse(visit=="day 14" & arm=="Non-ACT-PQ",NA,RR),
         ci_upper = ifelse(visit=="day 14" & arm=="Non-ACT-PQ",NA,ci_upper),
         ci_lower = ifelse(visit=="day 14" & arm=="Non-ACT-PQ",NA,ci_lower))

ggplot(data=results.all, aes(x=arm, y=RR/100, fill=visit, group=visit))+
  #geom_point(position = position_dodge2(0.8, preserve = "single"))+
  geom_bar(stat="identity", position = position_dodge(preserve = "single"), alpha=0.8, width=0.75)+
  #geom_point(data=medgams, aes(x=arm, y=RR/100, group=visit), col="red", position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=ci_lower/100, ymax=ci_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="grey30")+
  #geom_errorbar(data=medgams,aes(x=arm,ymin=RR_lower/100, ymax=RR_upper/100), position = position_dodge(preserve = "single", width=0.8), width=0.15, linewidth=0.8,color="red",alpha=0.2)+
  theme_prism(base_size = 14, base_line_size = 0.2)+
  theme(axis.title.x = element_blank(),
        legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))+
  scale_fill_viridis_d(begin=0.2, end=0.8)+
  scale_y_continuous(labels=scales::percent,breaks=c(-0.5,0,0.5,1))+
  ylab("Relative reduction (%) in \nindividual's infectious likelihood")+
  coord_cartesian(ylim = c(-0.5, NA))+
  geom_text(data=data.frame(arm="Non-ACT-PQ", day="Day 14", TRA=0),
            aes(x=arm, y=0.5, label="No Data", group=day), # Use the same group aesthetic
            position=position_dodge(preserve="single", width=-0.5), # Align with bar dodge
            angle=90, size=5, color="black", inherit.aes=FALSE)


ggsave("OUTPUT/NMA_updated/infectiousness_NMA.png",device="png", width=10, height=7, dpi=600)


arms = c("DHA-PPQ", "SP-AQ", "PY-AS", "AS-AQ", "AL", "Non-ACT-PQ", "ACT-PQ", "ACT-TQ")

mod_main2 = bam(prev ~ arm*visit*study , data = necfeedsub, family=poisson(link="log"), discrete=T, nthread=4)
summary(mod_main2)

result_day2 <- get_studywise_contrasts_bam(mod_main2, necfeedsub, "day 2", arms)
flextable(result_day2[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/infectiousness_day2_NMA.docx")
result_day7 <- get_studywise_contrasts_bam(mod_main2, necfeedsub, "day 7", arms)
flextable(result_day7[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/infectiousness_day7_NMA.docx")
arms_up = arms[which(!(arms %in% "Non-ACT-PQ"))]
result_day14 <- get_studywise_contrasts_bam(mod_main2, necfeedsub, "day 14", arms_up)
flextable(result_day14[[2]]) %>% autofit() %>% save_as_docx( path = "OUTPUT/NMA_updated/infectiousness_day14_NMA.docx")



ggarrange(result_day2[[1]]+
            ggtitle("Day 2")+
            theme(axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_x_discrete(limits = comparisons)
          ,
          result_day7[[1]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            ggtitle("Day 7")+
            scale_x_discrete(limits = comparisons)
          ,
          result_day14[[1]]+
            theme(axis.text.y = element_blank(),
              axis.title = element_blank(),
              panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            ggtitle("Day 14")+
            scale_x_discrete(limits = comparisons)
          ,
          legend="top",
          common.legend = TRUE,
          nrow=1, ncol=3,
          widths = c(1.4, 1, 1)
)
ggsave("OUTPUT/NMA_updated/CHECK1_infectiousness_NMA.png", dpi=600, device="png", width=28, height=24, units="cm")



png("OUTPUT/NMA_updated/FOREST_infectiousness_NMA.png", width = 18, height = 24, units = "cm", res = 1200)
grid.newpage()
gridExtra::grid.arrange(result_day2[[4]],result_day7[[4]],result_day14[[4]])
dev.off()

# Extract y-axis levels from plot 1 (to use as standard)
p1_data <- result_day2[[6]]
common_levels <- p1_data$arm


ggarrange(result_day2[[3]]+
            theme(axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            ),
          result_day7[[3]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels),
          result_day14[[3]]+
            theme(axis.text.y = element_blank(),
                  axis.title = element_blank(),
                  panel.border = element_rect(color = "black", fill = NA, size = 0.5)
            )+
            scale_y_discrete(limits=common_levels),
          legend="top",
          common.legend = TRUE,
          nrow=1, ncol=3,
          widths = c(1.4, 1, 1)
)
ggsave("OUTPUT/NMA_updated/pscore_infectiousness_NMA.png", dpi=600, device="png", width=28, height=14, units="cm")


result_day2[[5]]+theme(title = element_blank())  + coord_cartesian(clip = "off")

ggsave("OUTPUT/NMA_updated/network_infectiousness_NMA.png", dpi=600, device="png", width=14, height=22, units="cm")



