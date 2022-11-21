# PBC -- Making 'base' dataset

rm(list=ls())
library(dplyr)
pbc <- joineRML::pbc2 %>% as_tibble %>% 
  select(-status) %>% 
  rename(survtime = years, time = year, status = status2)
pbc$drug <- as.numeric(pbc$drug) - 1 # 1: active, 0: placebo
pbc$sex <- as.numeric(pbc$sex) - 1 # 1: female, 0: male

survdata <- pbc %>%
  distinct(id, survtime, status, age, drug) %>% 
  mutate(`age` = as.numeric(scale(age)))

pbc <- left_join(
  pbc %>% select(-age),
  survdata %>% select(id, age), 
  'id'
)

# Response clear-up
pbc$ascites <- as.numeric(pbc$ascites) - 1 # 1: Yes
pbc$spiders <- as.numeric(pbc$spiders) - 1 # 1: Yes
pbc$hepatomegaly <- as.numeric(pbc$hepatomegaly) - 1 # 1: Yes
pbc <- as.data.frame(pbc)
save(pbc, file = 'PBC.RData')