source("data prep and cleaning.R")

# Descriptive results -----------------------------------------------------

#We create a dataset for the baseline data [i.e. at visit 0] and make a summarising function for the descriptives.

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


#We create a dataset for the baseline data [i.e. at visit 0] and make a summarising function for the descriptives.


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


#correspondence between gametocyte densities by microscopy and by PCR - BLAND ALTMAN



# Subset nec0$year to match the length of ba.stats$means and ba.stats$diffs
nec0_subset <- nec0[!is.na(nec0$gam.m) & !is.na(nec0$gam.p), ]

# Bland-Altman statistics
ba.stats <- bland.altman.stats(nec0_subset$gam.m, nec0_subset$gam.p)

# Create the data frame for plotting
ba.data <- data.frame(
  means = ba.stats$means,
  diffs = ba.stats$diffs,
  year = factor(nec0_subset$year)  # Subset year to match
)

# Generate the plot
ggplot(ba.data, aes(x = means, y = diffs, color = year)) +
  geom_point(size = 3,alpha=0.5) +
  geom_hline(yintercept = ba.stats$lines[2], linetype = "dashed", color = "black", size = 1) +
  geom_hline(yintercept = ba.stats$lines[1], linetype = "dotted", color = "black", size = 1) +
  geom_hline(yintercept = ba.stats$lines[3], linetype = "dotted", color = "black", size = 1) +
  scale_color_viridis_d(name = "Year") +
  xlim(0, 1000) +
  ylim(-750, 750) +
  theme_minimal() +
  labs(
    x = "Mean of Measurements",
    y = "Difference of Measurements",
    title = "Bland-Altman Plot"
  )


### Calculated per group (year)
# Subset nec0 for non-missing values
nec0_subset <- nec0 %>% filter(!is.na(gam.m) & !is.na(gam.p))

# Calculate Bland-Altman stats for each group (year)
ba_list <- nec0_subset %>%
  group_by(year) %>%
  summarise(
    mean_diff = mean(gam.m - gam.p),
    lower_limit = mean(gam.m - gam.p) - 1.96 * sd(gam.m - gam.p),
    upper_limit = mean(gam.m - gam.p) + 1.96 * sd(gam.m - gam.p),
    means = list((gam.m + gam.p) / 2),
    diffs = list(gam.m - gam.p)
  ) %>%
  unnest(cols = c(means, diffs))

# Create a data frame for horizontal lines
hline_data <- nec0_subset %>%
  group_by(year) %>%
  summarise(
    mean_diff = mean(gam.m - gam.p),
    lower_limit = mean(gam.m - gam.p) - 1.96 * sd(gam.m - gam.p),
    upper_limit = mean(gam.m - gam.p) + 1.96 * sd(gam.m - gam.p)
  )

# Plot Bland-Altman points and horizontal lines
ggplot(ba_list, aes(x = means, y = diffs, color = factor(year))) +
  # Black horizontal lines for overall statistics
  geom_hline(yintercept = ba.stats$lines[2], linetype = "solid", color = "black", size = 1) +  # Overall mean difference
  geom_hline(yintercept = ba.stats$lines[1], linetype = "dotted", color = "black", size = 1) +  # Overall lower limit
  geom_hline(yintercept = ba.stats$lines[3], linetype = "dotted", color = "black", size = 1) +  # Overall upper limit
  # Group-specific horizontal lines
  geom_hline(data = hline_data, aes(yintercept = mean_diff, color = factor(year)), linetype = "dashed", size = 0.7) +
  # Points colored by year
  geom_point(size = 2, alpha = 0.5) +
  # Viridis color scale for groups
  scale_color_viridis_d(name = NULL) +
  theme_minimal() +
  xlim(0, 1000) +
  ylim(-1000, 1000) +
  labs(
    x = "Mean of Measurements",
    y = "Difference of Measurements"
  )+
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = element_text(size = 10),
    plot.margin = margin(20, 20, 40, 20),   # Adjust margins for more space
    panel.border = element_rect(color = "black", fill = NA, size = 0.5))+
  guides(color = guide_legend(nrow = 1)) 


ggsave("OUTPUT/Fig_BlandAltman.pdf",device="pdf", width=7, height=5.5, dpi=600)


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


