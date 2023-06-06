rm(list=ls())
data(PBC, package='gmvjoint')

PBC <- na.omit(PBC[,c("id", "survtime", "status", "drug", "sex", "age", "time", 
                      "ascites", "hepatomegaly", "spiders", "serBilir",
                      "albumin", "alkaline", "SGOT", "platelets", "prothrombin")])
PBC$serBilir <- log(PBC$serBilir)
PBC$prothrombin <- (PBC$prothrombin * .1)^ (-4)
PBC$SGOT <- log(PBC$SGOT)

library(ggplot2)
library(dplyr)
library(forcats)
library(tidyr)
source('../gmvjoint/zzz/theme_csda.R')

PBC %>% 
  filter(status == 1) %>% 
  pivot_longer(cols=`serBilir`:`prothrombin`, names_to='biomarker') %>% 
  mutate(biomarker = case_when(
    biomarker == 'serBilir' ~ "log(Serum~bilirubin)",
    biomarker == 'albumin' ~ "Albumin",
    biomarker == 'prothrombin' ~ '(0.1 ~ x ~ Prothrombin~time)^{-4}',
    biomarker == 'SGOT' ~ "log(AST)",
    biomarker == "platelets" ~ "Platelet~count",
    biomarker == 'alkaline' ~ "Alkaline~phosphatase",
    T ~ 'AA'
  ),
  tt = -1 * (survtime-time)
  ) %>% 
  mutate(f.biomarker = factor(biomarker, levels = c('log(Serum~bilirubin)',
                                                    'log(AST)', 'Albumin', '(0.1 ~ x ~ Prothrombin~time)^{-4}',
                                                    'Platelet~count', "Alkaline~phosphatase"))) %>% 
  filter(biomarker != 'AA') %>% 
  ggplot(aes(x=tt, y = value, group = id)) + 
  geom_vline(xintercept = 0, colour = 'black', alpha = .25) + 
  geom_line(alpha = .10) + 
  geom_smooth(aes(group=NULL), colour = RColorBrewer::brewer.pal(3,'YlOrRd')[3], method = 'loess', formula = y~x) +
  facet_wrap(~f.biomarker, scales = 'free', strip.position = 'left', labeller = label_parsed) + 
  scale_y_continuous(breaks = scales::pretty_breaks(10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  labs(y = NULL,
       x = 'Time (years) from death (0: time of death)') + 
  theme_csda() + 
  theme(strip.placement = 'outside',
        strip.text = element_text(vjust = 1))
ggsave('~/Downloads/PBCtrajectories.png', width = 190, height = 120, units = 'mm')
tiff('~/Downloads/PBCtrajectories.tiff',
     width = 190, height = 120, units = 'mm',
     res = 1e3, compression = 'lzw')
last_plot()
dev.off()


# Version showing survived/died -------------------------------------------
PBC %>% 
  mutate(prothrombin = (10 * prothrombin^-.25)) %>% 
  pivot_longer(cols=`serBilir`:`prothrombin`, names_to='biomarker') %>% 
  mutate(biomarker = case_when(
    biomarker == 'serBilir' ~ "log(Serum~bilirubin)",
    biomarker == 'albumin' ~ "Albumin",
    # biomarker == 'prothrombin' ~ '(0.1 ~ x ~ Prothrombin~time)^{-4}',
    biomarker == 'prothrombin' ~ 'Prothrombin~time',
    biomarker == 'SGOT' ~ "log(AST)",
    biomarker == "platelets" ~ "Platelet~count",
    biomarker == 'alkaline' ~ "Alkaline~phosphatase",
    T ~ 'AA'
  ),
  tt = -1 * (survtime-time)
  ) %>% 
  mutate(f.biomarker = factor(biomarker, levels = c('log(Serum~bilirubin)',
                                                    'log(AST)', 'Albumin', 'Prothrombin~time',
                                                    'Platelet~count', "Alkaline~phosphatase")),
         f.surv = factor(status, levels = c(0,1),labels=c("Alive", "Dead"))) %>% 
  filter(biomarker != 'AA') %>% 
  ggplot(aes(x=time, y = value, group = id)) + 
  geom_vline(xintercept = 0, colour = 'black', alpha = .25) + 
  geom_line(alpha = .10) + 
  # geom_smooth(aes(colour=f.surv, group = f.surv), method = 'loess', formula = y~x, se = F) +
  geom_smooth(aes(colour=f.surv, group = f.surv), method = 'lm', formula = y~splines::ns(x, knots = c(1,4)), se = F) +
  facet_wrap(~f.biomarker, scales = 'free', strip.position = 'left', labeller = label_parsed) + 
  scale_y_continuous(breaks = scales::pretty_breaks(10)) +
  scale_x_continuous(breaks = scales::pretty_breaks(10)) +
  scale_color_brewer(palette = "Dark2") +
  labs(y = NULL,
       x = 'Time (years)',
       colour = "") + 
  theme_csda() + 
  theme(strip.placement = 'outside', legend.position = 'bottom',
        strip.text = element_text(vjust = 1))
  
ggsave("~/Downloads/aaa.png", width = 190, height = 120, units = 'mm')
