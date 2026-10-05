###################################################
##AUTHOR: Sean Galagan 
##DATE STARTED: 18NOV2022
##PURPOSE: Builds Last Step Datasets for Analysis
##NOTES:
##LAST UPDATED: 30JAN2026
###################################################

## ENSURE NECESSARY PACKAGES ARE INSTALLED
packages = c("ezknitr","tidyverse","haven","wesanderson","knitr","ggplot2", "ggthemes", "scales", "data.table", "lubridate", "janitor","stringr", "consort", "anthro", "rmarkdown", "zoo", "conflicted")

## Now load or install & load all
package.check <- lapply( packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
}
)

# build quick function
rm(list=ls())
`%notin%` <- Negate(`%in%`)

## Load necessary data files & split by country
today <- as.Date((paste0(format(Sys.Date(),"%Y%m%d"))),format="%Y%m%d")

# Set preferred conflicts
conflicts_prefer(lubridate::month) # conflicts with data.table
conflicts_prefer(lubridate::year) # conflicts with data.table
conflicts_prefer(dplyr::filter) # conflicts with stats

######LOAD DATA#####
# PE/HUS
pehus <- readRDS("./Data/PE_HUS_Final.Rds") 
cluster_comp <- readRDS("./Data/PE_ClusterCompletion.Rds") %>% 
  mutate(pop_cluster_id=clcom_cluster_id) 
careseek <- readRDS("./Data/PE_HUS_Careseek.Rds") %>% 
  filter(!is.na(hus_careseek_efgh)) %>%
  inner_join(pehus %>% select(pop_hh_id,link_key,pop_country,hus_blood,hus_consent,re_consent), by = c("pop_hh_id","link_key")) 
# Load full cluster datasets
bangladesh <- read.csv("./Data/bg_v3.csv") %>% 
  mutate(country=1) %>% select(country, ClusterID, Area_m3)
kenya <- read.csv("./Data/ke.csv") %>% 
  mutate(country=2) %>% select(country, ClusterID, Area_m3)
malawi <- read.csv("./Data/mw.csv") %>% 
  mutate(country=3) %>%  rename(ClusterID=Cluster_id) %>% select(country, ClusterID, Area_m3)
mali <- read.csv("./Data/ml_v3.csv") %>% 
  mutate(country=4) %>% select(country, ClusterID, Area_m3)
pakistan <- read.csv("./Data/pk_v7.csv") %>% 
  mutate(country=5) %>% select(country, ClusterID, Area_m3)
peru <- read.csv("./Data/pe_v3.csv") %>% 
  mutate(country=6) %>% rename(Area_m3=Area__m3_) %>% select(country, ClusterID, Area_m3)
gambia <- read.csv("./Data/gm_v4.csv") %>% 
  mutate(country=7) %>%  rename(ClusterID=NewCluster) %>% select(country, ClusterID, Area_m3)
clusters <- bind_rows(bangladesh, kenya, malawi, mali, peru, pakistan, gambia) %>%
  rename(clcom_country=country) |>
  distinct()
# DCS screening, enrollment and lab testing - excluded 6302484 and 7104928
prescreening <- readRDS("./Data/DCS_01_prescreening.Rds")  %>% 
  mutate(temp=substr(as.character(serial_id),1,1))
screening <- readRDS("./Data/DCS_02_screening.Rds")
preenrollment <- readRDS("./Data/DCS_03_preenrollment.Rds")
enrollment <- readRDS("./Data/DCS_04_enrollment.Rds") 
anthro <- readRDS("./Data/DCS_anthro.Rds")
rectal_swab_results <- readRDS("./Data/DCS_16_rectal_swab_results.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid")
tac_long <- readRDS("./Data/DCS_tac_long.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid")
tac <- readRDS("./Data/DCS_tac_automated.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid")
severity <- readRDS("./Data/DCS_severity_scores.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid") 
stool_collection <- readRDS("./Data/DCS_13_stool_collection.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid")
# Propensity to seek care weights
prop_care <- readRDS("./Data/enrollment_ps_weights_20250408.Rds") |>
  inner_join(enrollment %>% select(pid), by = "pid")
prop_care_7days <- readRDS("./Data/enrollment_ps_weights_7days_20250408.Rds") |>
  inner_join(enrollment %>% select(pid), by = "pid")

# remove partially enumerated cluster in Pakistan
pehus2 <- pehus |>
  filter(pop_country!=5 | (pop_cluster_id!=7 & pop_country==5))

############################################
######## BUILD DENOMINATOR #################
############################################
# Prep cluster completion dataset 
cluster_comp2 <- cluster_comp %>%
  filter(clcom_households %in% c(1,2)) |> # remove unenumeratd clusters
  select(clcom_country, clcom_cluster_id,clcom_subcluster_id) %>%
  rename(pop_cluster_id=clcom_cluster_id, pop_subcluster_id=clcom_subcluster_id) %>%
  mutate(cluster_complete=1,
         pop_country=clcom_country)

cluster_comp3 <- cluster_comp %>%
  filter(clcom_households %in% c(1,2)) |> # remove unenumeratd clusters
  select(clcom_country, clcom_cluster_id,clcom_subcluster_id) %>%
  rename(pop_cluster_id=clcom_cluster_id, pop_subcluster_id=clcom_subcluster_id) %>%
  mutate(cluster_complete=1,
         pop_country=clcom_country)

# Compute number of children under-five enumerated
raw_denom_site <- pehus2 %>%
  left_join(cluster_comp2, by = c("pop_country", "pop_cluster_id","pop_subcluster_id")) %>% 
  filter(cluster_complete==1) %>% # Only keep records from completed clusters
  filter(pop_age>=6 & pop_age <36) %>% # Keep kids 6-35 months
  group_by(pop_country,pop_cluster_id) %>%
  summarise(n = n()) %>%
  rename(raw_denom=n)

# FOR CM analysis - raw total population per site
raw_denom_adults <- pehus2 %>%
  left_join(cluster_comp2, by = c("pop_country", "pop_cluster_id","pop_subcluster_id")) %>% 
  filter(cluster_complete==1) %>% # Only keep records from completed clusters
  filter(pop_hhold_consent==1) %>% # Keep enumerated households
  select(pop_country,pop_cluster_id, pop_hh_id, pop_num_residents) |>
  unique() |>
  group_by(pop_country,pop_cluster_id) %>%
  summarise(n = sum(pop_num_residents)) %>%
  rename(adult_denom=n)

# Combine
raw_denom_overall <- raw_denom_adults |>
  left_join(raw_denom_site, by = c("pop_country","pop_cluster_id")) |>
  mutate(raw_denom=ifelse(is.na(raw_denom),0,raw_denom))
  
  # SUB_GROUPS
  # SUBGROUP 1 - by age
  age_denom_site1 <- pehus2 %>%
    left_join(cluster_comp2, by = c("pop_country", "pop_cluster_id","pop_subcluster_id")) %>% 
    filter(cluster_complete==1) %>% # Only keep records from completed clusters
    filter(pop_age>=6 & pop_age <36) %>% # Keep kids 6-35 months
    mutate(c_agegroup=ifelse(pop_age>=6  & pop_age <9, 1,
                      ifelse(pop_age>=9  & pop_age <12,2,     
                      ifelse(pop_age>=12 & pop_age <18,3,
                      ifelse(pop_age>=18 & pop_age <24,4,
                      ifelse(pop_age>=24 & pop_age <36,5,NA)))))) %>%
    group_by(pop_country,pop_cluster_id, c_agegroup) %>%
    summarise(n = n()) %>%
    rename(raw_denom=n)
  
raw_denom <- bind_rows(raw_denom_overall,
                       age_denom_site1) 

# Create proportion of households enumerated variable by site
pop_enum <- pehus2 %>%
  left_join(cluster_comp2, by = c("pop_country", "pop_cluster_id")) %>% 
  filter(cluster_complete==1) %>% # Only keep records from completed clusters
  mutate(enumerated=ifelse(pop_hhold_consent==2 | pop_hhold_approach %in% c(3,4),0,1)) %>% # Create variable that marks households as completed
  select(pop_country, pop_hh_id, pop_cluster_id,enumerated, pop_hhold_approach, pop_hhold_consent)

pop_enum_stats_site <- unique(pop_enum) %>% # remove duplicate entries -want one per house
  group_by(pop_country, pop_cluster_id, enumerated) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n),
         enum_weight=1/freq) %>%
  filter(enumerated==1) %>%
  select(pop_country, pop_cluster_id,freq, enum_weight)

# Cluster adjustment
# part 1 - determine total area of site
total_area_site <- clusters %>%
  group_by(clcom_country) %>%
  summarise(total_area=sum(Area_m3),
            total_clusters=n())

total_area_total <- clusters %>%
  mutate(clcom_country=8) %>% #allows groupby to work
  group_by(clcom_country) %>%
  summarise(total_area=sum(Area_m3),
            total_clusters=n())

# part 2: area of cluster completed 
complete_area <- cluster_comp2 %>%
  rename(ClusterID=pop_cluster_id) %>%
  left_join(clusters, by = c("clcom_country", "ClusterID")) 

complete_area_site <- complete_area %>%
  group_by(clcom_country) %>%
  summarise(clusters_completed = n(),
            complete_area=sum(Area_m3)) %>%
  left_join(total_area_site, by = "clcom_country") %>%
  mutate(cluster_freq=complete_area/total_area,
         cluster_weight=1/cluster_freq) %>%
  rename(pop_country=clcom_country)

complete_area_total <- complete_area %>%
  mutate(clcom_country=8) |>
  group_by(clcom_country) %>%
  summarise(clusters_completed = n(),
            complete_area=sum(Area_m3)) %>%
  left_join(total_area_total, by = "clcom_country") %>%
  mutate(cluster_freq=complete_area/total_area,
         cluster_weight=1/cluster_freq) %>%
  rename(pop_country=clcom_country)

# Final denominator
denom_final_merge <- raw_denom |>
  left_join(pop_enum_stats_site, by=c("pop_country","pop_cluster_id")) |>
  left_join(complete_area_site, by = "pop_country") |>
  mutate(adj_denom=raw_denom*enum_weight*cluster_weight, # creates estimated number of kids
         adj_adults=adult_denom*enum_weight*cluster_weight,
         cluster_freq=cluster_freq*100,
         freq=freq*100) %>%
  # Create country sum
  group_by(pop_country,c_agegroup) |>
  summarise(raw_denom=sum(raw_denom),
            raw_adults=sum(raw_denom_adults),
            hh_freq=mean(freq),
            clusters_completed=max(clusters_completed),
            cluster_freq=max(cluster_freq),
            adj_denom=sum(adj_denom),
            adj_adults=sum(adj_adults)) %>%
  rename(enroll_site=pop_country)

denom_final_total <- denom_final_merge %>%
  mutate(enroll_site=8) %>%
  group_by(enroll_site,c_agegroup) %>%
  summarise(raw_denom=sum(raw_denom),
            raw_adults=sum(raw_adults),
            hh_freq=mean(hh_freq),
            adj_denom=sum(adj_denom),
            adj_adults=sum(adj_adults)) |>
  left_join(complete_area_total |> select(pop_country, cluster_freq, clusters_completed) |> 
                                   rename(enroll_site=pop_country) |>
                                   mutate(cluster_freq=cluster_freq*100), by = "enroll_site")
  
denom_final <- denom_final_merge %>%
  bind_rows(denom_final_total) 

# Export key dataset
write_rds(denom_final |> select(-c(adj_adults,raw_adults)), 
          paste0("./Last Step Datasets/denominator.Rds"))

############################################
### BUILD HEALTHCARE SEEKING ADJUSTMENT ####
############################################
# Identify missing careseking information - CURRENTLY NONE (11/25/2024 SRG)
test_careseeking <- pehus2 %>%
  filter((hus_consent %in% c(1,2) | re_consent %in% c(1,2))) %>%
  select(pop_hh_id,pop_country,link_key,hus_consent,re_consent,hus_blood) %>%
  full_join(careseek %>% select(pop_hh_id,link_key,hus_careseek_efgh), by = c("pop_hh_id","link_key")) %>%
  filter(is.na(hus_careseek_efgh))

###### Create A-adjustment
# First figure out total diarrhea cases
## First pull out number of cases by diarrhea type
total_diarrhea <- pehus2 %>%
  filter((hus_consent %in% c(1,2) | re_consent %in% c(1,2)) & hus_blood!=3) 

total_diarrhea_site <- total_diarrhea %>%
  group_by(pop_country, hus_blood) %>%
  summarise(diarrhea=n()) 

total_diarrhea_total <- total_diarrhea %>%
  mutate(pop_country=8) %>%
  group_by(pop_country, hus_blood) %>%
  summarise(diarrhea=n())

## Combine into single dataset
total_diarrhea_merge <- total_diarrhea_site %>%
  bind_rows(total_diarrhea_total)

# now figure out who sought care
## Pull out kids who sought care at PHC facilities
total_care_site <- careseek %>%
  filter(hus_blood %in% c(1,2) & hus_careseek_mad==1) %>%
  group_by(pop_country, hus_blood) %>%
  summarize(sought_care=n())

total_care_total <- careseek %>%
  filter(hus_blood %in% c(1,2) & hus_careseek_mad==1) %>%
  mutate(pop_country=8) %>%
  group_by(pop_country, hus_blood) %>%
  summarize(sought_care=n())

## Pull out kids who sought care at PHC facilities - WITHIN 7 DAYS FOR SENSITIVITY ANALYSIS
total_care_site_7days <- careseek %>%
  filter(hus_blood %in% c(1,2) & hus_careseek_mad_7days==1) %>%
  group_by(pop_country, hus_blood) %>%
  summarize(sought_care_7days=n())

total_care_total_7days <- careseek %>%
  filter(hus_blood %in% c(1,2) & hus_careseek_mad_7days==1) %>%
  mutate(pop_country=8) %>%
  group_by(pop_country, hus_blood) %>%
  summarize(sought_care_7days=n())

## Combine
total_care_merge <- bind_rows(total_care_site,total_care_total) |>
  left_join(bind_rows(total_care_site_7days,total_care_total_7days), by = c("pop_country","hus_blood"))

# Pull out total by site sought care and sought care at EFGH facility to determine ratio for secondary weight
## Any care seeking
total_care_site_overall <- careseek %>%
  filter(hus_careseek_mad==1) %>%
  group_by(pop_country) %>%
  summarize(sought_care_overall=n())

total_care_total_overall <- careseek %>%
  filter(hus_careseek_mad==1) %>%
  mutate(pop_country=8) %>%
  group_by(pop_country) %>%
  summarize(sought_care_overall=n())

## Any care seeking - WITHIN 7 DAYS FOR SENSITIVITY ANALYSIS
total_care_site_overall_7days <- careseek %>%
  filter(hus_careseek_mad_7days==1) %>%
  group_by(pop_country) %>%
  summarize(sought_care_overall_7days=n())

total_care_total_overall_7days <- careseek %>%
  filter(hus_careseek_mad_7days==1) %>%
  mutate(pop_country=8) %>%
  group_by(pop_country) %>%
  summarize(sought_care_overall_7days=n())

## Combine
total_care_merge_overall <- bind_rows(total_care_site_overall,total_care_total_overall) |>
  left_join(bind_rows(total_care_site_overall_7days,total_care_total_overall_7days), by = "pop_country")

# Care seeking at EFGH facilities
## EFGH facilities
total_efgh_site_overall <- careseek %>%
  filter(hus_careseek_efgh==1) %>%
  group_by(pop_country) %>%
  summarize(sought_care_efgh=n())

total_efgh_total_overall <- careseek %>%
  filter(hus_careseek_efgh==1) %>%
  mutate(pop_country=8) %>%
  group_by(pop_country) %>%
  summarize(sought_care_efgh=n())

## EFGH facilities- WITHIN 7 DAYS FOR SENSITIVITY ANALYSIS
total_efgh_site_overall_7days <- careseek %>%
  filter(hus_careseek_efgh_7days==1) %>%
  group_by(pop_country) %>%
  summarize(sought_care_efgh_7days=n())

total_efgh_total_overall_7days <- careseek %>%
  filter(hus_careseek_efgh_7days==1) %>%
  mutate(pop_country=8) %>%
  group_by(pop_country) %>%
  summarize(sought_care_efgh_7days=n())

## Combine
total_efgh_merge_overall <- bind_rows(total_efgh_total_overall,total_efgh_site_overall) |>
  left_join(bind_rows(total_efgh_total_overall_7days,total_efgh_site_overall_7days), by = "pop_country")

# calculate a_factor
a_factor <- total_diarrhea_merge %>%
  left_join(total_care_merge, by = c("pop_country", "hus_blood")) %>%
  left_join(total_efgh_merge_overall, by = "pop_country") %>%
  left_join(total_care_merge_overall, by = "pop_country") %>%
  mutate(# Overall
         a_factor_step1=(sought_care/diarrhea)*100,
         a_factor_step2=(sought_care_efgh/sought_care_overall)*100,
         a_factor_weight=1/((a_factor_step1/100)*(a_factor_step2/100)),
         # sensitivity analysis
         a_factor_step1_sens=(sought_care_7days/diarrhea)*100,
         a_factor_step2_sens=(sought_care_efgh_7days/sought_care_overall_7days)*100,
         a_factor_weight_sens=1/((a_factor_step1_sens/100)*(a_factor_step2_sens/100)),
         b_dysentery=ifelse(hus_blood==2,0,1)) %>%
  rename(enroll_site=pop_country) |>
  select(-hus_blood)

# Export key dataset
write_rds(a_factor, paste0("./Last Step Datasets/a_factor.Rds"))

############################################
####### BUILD ENROLLMENT ADJUSTMENT ########
############################################
# See which facilities have screened and enrolled
screen_facilities <- prescreening %>%
  mutate(fac_id=substr(serial_id,1,2)) %>%
  group_by(pscr_country, fac_id) %>%
  summarise(n=n())

enroll_facilities <- enrollment %>%
  mutate(fac_id=substr(pid,1,2)) %>%
  group_by(enroll_site, fac_id) %>%
  summarise(n=n())

##### PART 1: Screening Adjustment ######
# Create base dataset with definition of diarrhea
enrollment_dataset <- prescreening %>%
  left_join(screening, by="serial_id") %>%
  left_join(preenrollment, by="serial_id") %>%
  left_join(enrollment %>% mutate(enrolled=1), by="pid") %>%
  mutate(enrolled=ifelse(is.na(enrolled),0,1),
         fac_id=substr(serial_id,1,2), # recode Bangladesh additional IDs for icddr,b
         fac_id=ifelse(fac_id %in% c(16,17,18,19,81,82,83),14,fac_id),
         b_dysentery=scr_diar_blood) |>
  mutate(fac_id2=substr(pid,1,2)) 

table(enrollment_dataset$fac_id,useNA="ifany") # all accounted for

# Check for mismatched facility IDs in the screening versus enrollment datasets
test <- enrollment_dataset |>
  filter(fac_id!=fac_id2) |>
  select(serial_id, pid, fac_id, fac_id2)

# LOOK FOR MISMATCHES
test <- enrollment_dataset %>%
  filter(pre_enroll==1 & !is.na(pid)) %>%
  full_join(enrollment, by="pid")  %>%
  select(serial_id,pid,language.y) # 11/25/2024 NO MISMATCHES

# Total enrolled 
total_enrolled_facility <- enrollment_dataset %>%
  filter(pre_enroll==1) %>%
  group_by(pscr_country,b_dysentery,fac_id) %>%
  summarise(enrolled=n())

total_enrolled_site <- enrollment_dataset %>%
  filter(pre_enroll==1) %>%
  group_by(pscr_country,b_dysentery) %>%
  summarise(enrolled=n())

total_enrolled_total <- enrollment_dataset %>%
  filter(pre_enroll==1) %>%
  group_by(b_dysentery) %>%
  summarise(enrolled=n())%>%
  mutate(pscr_country=8)

total_enrolled <- total_enrolled_site %>%
  bind_rows(total_enrolled_total) %>%
  bind_rows(total_enrolled_facility) 

# Total eligible for enrollment
total_eligible <- enrollment_dataset %>%
  filter(scr_elig_diar==1 & scr_elig_area==1 & scr_elig_age==1 & scr_elig_dhours==1 & scr_elig_ddays==1) %>% # Filter to individuals who meet study definition of age/diarrhea/study area
  filter(pre_enroll==1 | pre_consent==0 | pre_enroll==0 | scr_refer==2 | scr_enrollcap==1 | scr_diartime==1 | scr_elig_study==0 | scr_elig_months==0 | scr_elig_consent==0) %>% # Filter to individuals who were enrolled OR would have been able to be enrolled
  select(fac_id,b_dysentery,pscr_country,scr_elig_diar,scr_elig_area,scr_elig_age,scr_elig_dhours,scr_elig_ddays,scr_elig_months,scr_elig_consent,scr_refer,scr_enrollcap,scr_diartime,scr_elig_study,pre_enroll,pre_consent)

total_eligible_facility <- total_eligible %>%
  group_by(pscr_country,b_dysentery,fac_id) %>%
  summarise(eligible=n())

total_eligible_site <- total_eligible %>%
  group_by(pscr_country,b_dysentery) %>%
  summarise(eligible=n())

total_eligible_total <- total_eligible %>%
  group_by(b_dysentery) %>%
  summarise(eligible=n()) %>%
  mutate(pscr_country=8)

total_eligible_final <- total_eligible_site %>%
  bind_rows(total_eligible_total) %>%
  bind_rows(total_eligible_facility)

##### PART 2: Prescreening Adjustment ######
# Total pre-screened
prescreen_facility <- enrollment_dataset %>%
  group_by(pscr_country,fac_id) %>%
  summarise(prescreen=n())

prescreen_site <- enrollment_dataset %>%
  group_by(pscr_country) %>%
  summarise(prescreen=n())

prescreen_total <- enrollment_dataset %>%
  summarise(prescreen=n()) %>%
  mutate(pscr_country=8)

prescreen_final <- prescreen_site %>%
  bind_rows(prescreen_total) %>%
  bind_rows(prescreen_facility)

# Total screened
screen_site <- enrollment_dataset %>%
  subset(scr_consent==1 | scr_consent==2) %>%
  group_by(pscr_country) %>%
  summarise(screen=n())

screen_facility <- enrollment_dataset %>%
  subset(scr_consent==1 | scr_consent==2) %>%
  group_by(pscr_country,fac_id) %>%
  summarise(screen=n())

screen_total <-  enrollment_dataset %>%
  subset(scr_consent==1 | scr_consent==2) %>%
  summarise(screen=n()) %>%
  mutate(pscr_country=8)

screen_final<- screen_site %>%
  bind_rows(screen_total) %>%
  bind_rows(screen_facility)

# Total eligible for screening
screen_eligible_facility <- enrollment_dataset %>%
  filter((pscr_age>=6 & pscr_age<36 & pscr_diar==1 & pscr_already_enr==0) | (scr_consent==1 | scr_consent==2)) %>% # met screening criteria
  filter(is.na(scr_verify_cache2) | scr_verify_cache2==1 | (scr_verify_cache2 %in% c(2,3) & !is.na(scr_careg_yn))) %>% # remove Bangladesh children who were outside of the catchment area
  group_by(pscr_country,fac_id) %>%
  summarise(screen_elig=n())

screen_eligible_site <- enrollment_dataset %>%
  filter((pscr_age>=6 & pscr_age<36 & pscr_diar==1 & pscr_already_enr==0) | (scr_consent==1 | scr_consent==2)) %>% # met screening criteria
  filter(is.na(scr_verify_cache2) | scr_verify_cache2==1 | (scr_verify_cache2 %in% c(2,3) & !is.na(scr_careg_yn))) %>% # remove Bangladesh children who were outside of the catchment area
  group_by(pscr_country) %>%
  summarise(screen_elig=n()) 

screen_eligible_total <- screen_eligible_site %>%
  summarise(screen_elig=sum(screen_elig)) %>%
  mutate(pscr_country=8)

screen_eligible <- screen_eligible_site %>%
  bind_rows(screen_eligible_total) %>%
  bind_rows(screen_eligible_facility)

# add 0s for facility where no one was enrolled
missing_facility= data.frame(pscr_country=c(1,1), fac_id=c("11","11"), b_dysentery=c(1,0), eligible=c(0,0), enrolled=c(0,0))

# Final b_factor dataset
b_factor_p1 <- total_eligible_final %>%
  left_join(total_enrolled, by=c("pscr_country","b_dysentery","fac_id")) %>% # merge total eligible for enrollment and total enrolled
  bind_rows(missing_facility)
            
b_factor <- list(prescreen_final,screen_eligible,screen_final,b_factor_p1) %>% # merge in total prescreened, total eligible for screening and total screened
  reduce(left_join, by = c("pscr_country","fac_id")) %>%
  # create screening adjustments
  mutate(screen_adj=enrolled/eligible, # percent of those eligible for enrollment for non-operational reasons who were enrolled
         screen_adj_weight=1/screen_adj, # screening adjustment weight
         prescreen_adj=screen/screen_elig, # % of those with diarrhea in the right age who were screened
         prescreen_adj_weight=1/prescreen_adj, # prescreening weight
         b_factor=screen_adj*prescreen_adj, # combined % of those who met study diarrhea defiunition who were enrolled
         b_factor_weight=1/b_factor, # final weight
         # Convert to percentages
         screen_adj=screen_adj*100,
         prescreen_adj=prescreen_adj*100,
         b_factor=b_factor*100,
         # fill in 0s for fac_id==11
         screen_adj=case_when(fac_id=="11"  ~  0,
                              is.na(fac_id) ~ screen_adj,
                              TRUE          ~ screen_adj),
         screen_adj_weight=case_when(fac_id=="11"  ~  1,
                                     is.na(fac_id) ~ screen_adj_weight,
                                     TRUE          ~ screen_adj_weight),
         b_factor=case_when(fac_id=="11"  ~  0,
                            is.na(fac_id) ~ b_factor,
                            TRUE          ~ b_factor),
         b_factor_weight=case_when(fac_id=="11"  ~  1,
                                   is.na(fac_id) ~ b_factor_weight,
                                   TRUE          ~ b_factor_weight)) %>%
  rename(enroll_site=pscr_country) 

# Export key dataset
write_rds(b_factor, paste0("./Last Step Datasets/b_factor.Rds"))

############################################
######### BUILD NUMERATOR ##################
############################################
# Prep screening dataset
scr_data <- preenrollment |>
  left_join(screening, by="serial_id") %>% 
  select(scr_diar_blood,pid) |>
  filter(!is.na(pid)) %>%
  rename(b_dysentery=scr_diar_blood)

# Create adjustment for number of individuals with available culture and TAC results 
  # ALL cultures complete - no weight needed
  lab_avail <- enrollment %>%
    full_join(rectal_swab_results %>% select(pid) %>% unique() %>% mutate(lab_avail=1), by= "pid") %>%
    group_by(enroll_site,lab_avail) %>%
    summarise(culture_N = n()) %>%
    mutate(culture_avail = culture_N / sum(culture_N),
           culture_weight=1/culture_avail) %>%
    filter(lab_avail==1) %>%
    select(enroll_site, starts_with("culture"))
  
  # TAC weight computed by study month to account for seasonality (study month) of missing results
  tac_avail <- enrollment %>%
    mutate(pid=as.numeric(pid),
           #year=year(enroll_date),
           month=month(enroll_date)) %>%
    select(enroll_site, pid,month#,year
           ) |>
    full_join(tac_long %>% filter(target=="Shigella EIEC") |> select(pid) %>% mutate(tac_avail=1), by= "pid") %>%
    group_by(enroll_site,#year,
             month,tac_avail) %>%
    summarise(tac_n = n()) %>%
    mutate(tac_N=sum(tac_n),
           tac_avail2 = tac_n / tac_N,
           tac_weight=1/tac_avail2) %>%
    filter(tac_avail==1) 

  # Check that AT LEAST one sample exists in the TAC data for each year and month - TRUE - N=173
  test <- enrollment %>%
    mutate(pid=as.numeric(pid),
           year=year(enroll_date),
           month=month(enroll_date)) %>%
    select(enroll_site, pid,month,year) |>
    full_join(tac_long %>% filter(target=="Shigella EIEC") %>% select(pid) %>% mutate(tac_avail=1), by= "pid") %>%
    group_by(enroll_site,year,month) %>%
    summarise(tac_N = n())
  
  # Check for overall missingness of TAC overall and by site
  tac_check_site <- enrollment %>%
    mutate(pid=as.numeric(pid)) %>%
    left_join(tac_long %>% filter(target=="Shigella EIEC") %>% select(pid) %>% mutate(tac_avail=1), by= "pid") %>%
    group_by(enroll_site,tac_avail) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n),
           N=sum(n),
           tac_weight=1/freq) %>%
    filter(tac_avail==1) %>%
    select(n, N, freq, tac_weight) %>%
    rename(tac_avail=freq)
  
  tac_check_site_rectal <- enrollment %>%
    mutate(pid=as.numeric(pid)) %>%
    left_join(tac_long %>% filter(target=="Shigella EIEC" & sample_type=="Rectal") %>% select(pid) %>% mutate(tac_avail=1), by= "pid") %>%
    group_by(enroll_site,tac_avail) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n),
           N=sum(n),
           tac_weight=1/freq) %>%
    filter(tac_avail==1) %>%
    select(n, N, freq, tac_weight) %>%
    rename(tac_avail=freq)
  
  tac_check_overall <- enrollment %>%
    mutate(pid=as.numeric(pid)) %>%
    left_join(tac_long %>% filter(target=="Shigella EIEC") %>% select(pid) %>% mutate(tac_avail=1), by= "pid") %>%
    group_by(tac_avail) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n),
           N=sum(n),
           tac_weight=1/freq) %>%
    filter(tac_avail==1) %>%
    select(n, N, freq, tac_weight) %>%
    rename(tac_avail=freq)
  
  tac_check_overall_rectal <- enrollment %>%
    mutate(pid=as.numeric(pid)) %>%
    left_join(tac_long %>% filter(target=="Shigella EIEC" & sample_type=="Rectal") %>% select(pid) %>% mutate(tac_avail=1), by= "pid") %>%
    group_by(tac_avail) %>%
    summarise(n = n()) %>%
    mutate(freq = n / sum(n),
           N=sum(n),
           tac_weight=1/freq) %>%
    filter(tac_avail==1) %>%
    select(n, N, freq, tac_weight) %>%
    rename(tac_avail=freq)

# Build culture dataset
culture_results <- rectal_swab_results %>%
  filter(sw_isolate_num>0) %>%
  select(pid, swi_isolate_group, swi_isolate_type,swi_isolate_type_dys) %>%
  mutate(source="culture",
         dys_serotype=case_when(swi_isolate_group==1 & swi_isolate_type_dys=="Type 1" ~ 1,
                                swi_isolate_group==1 & (swi_isolate_type_dys=="Non-type 1") ~ 2,
                                swi_isolate_group==1 & is.na(swi_isolate_type_dys) ~ 3)) %>%
  distinct()

# Build TAC dataset  
tac_results <- tac_long %>%
  filter(target=="Shigella EIEC" & b_attributable==1) |>
  mutate(enroll_site=case_when(groupid == "Bangladesh" ~ 1,
                               groupid == "Kenya"      ~ 2,
                               groupid == "Malawi"     ~ 3,
                               groupid == "Mali"       ~ 4,
                               groupid == "Pakistan"   ~ 5,
                               groupid == "Peru"       ~ 6,
                               groupid == "Gambia"     ~ 7),
         # Species
         swi_isolate_group=case_when(str_detect(serotype,"flexneri")                 ~ 2,
                                     str_detect(serotype,"sonnei")                   ~ 4,
                                     serotype %in% c("not typed","Not typed")        ~ 5,
                                     serotype=="" | is.na(serotype)                  ~ 5), # n=3 (1500506, 1500631, 7100230) to be QCed 
         # Serotype
         swi_isolate_type=case_when(str_detect(serotype,"1a") & swi_isolate_group==2 ~ 1,
                                    str_detect(serotype,"1b") & swi_isolate_group==2 ~ 2,
                                    str_detect(serotype,"1d") & swi_isolate_group==2 ~ 3,
                                    str_detect(serotype,"2a") & swi_isolate_group==2 ~ 4,
                                    str_detect(serotype,"2b") & swi_isolate_group==2 ~ 5,
                                    str_detect(serotype,"3a") & swi_isolate_group==2 ~ 6,
                                    str_detect(serotype,"3b") & swi_isolate_group==2 ~ 7,
                                    str_detect(serotype,"4a") & swi_isolate_group==2 ~ 8,
                                    str_detect(serotype,"4b") & swi_isolate_group==2 ~ 9,
                                    str_detect(serotype,"5a") & swi_isolate_group==2 ~ 10,
                                    str_detect(serotype,"5b") & swi_isolate_group==2 ~ 11,
                                    str_detect(serotype,"6")  & swi_isolate_group==2 ~ 12,
                                    str_detect(serotype,"7a") & swi_isolate_group==2 ~ 17,
                                    str_detect(serotype,"X")  & swi_isolate_group==2 ~ 13,
                                    str_detect(serotype,"Y")  & swi_isolate_group==2 ~ 14,
                                    swi_isolate_group==2                             ~ 16),
         # Hardcode PIDs where species/type order creates a problem
         swi_isolate_group=ifelse(pid %in% c(1504661,4104587,7105978,7112013,6501040,6201142,6301775,
                                             6301782,6201972,6201985,6100845,6201997,6302722,5601798,5602202),4,swi_isolate_group),
         swi_isolate_type=ifelse(pid %in% c(6501040,6201142,6301775,6301782,4105513,5602202,6201972,6100845,6201985,
                                            6302722,6201997,6301571,7100522,6201011,1504661,7105978,4104587,6201030),NA,
                          ifelse(pid %in% c(5601935,6201468,6302189,6400742,6501376,6100393),4,
                          ifelse(pid==5100965,5,
                          ifelse(pid %in% c(7101152,7110848,1504864),6,
                          ifelse(pid==4105463,11,
                          ifelse(pid==5602838,12,
                          ifelse(pid==3102209,13,
                          ifelse(pid %in% c(7112013,5601798),NA,swi_isolate_type))))))))) %>%
  select(pid,ct,sample_type,swi_isolate_group,swi_isolate_type,serotype) %>%
  mutate(source="tac")

  # PULL OUT COINFECTIONS AND APPEND THEM INTO DATASET - HAVE TO CHECK THESE EVERY TIME
  tac_results_coinfections <- tac_long %>%
    filter(target=="Shigella EIEC" & b_attributable==1 & str_detect(serotype,",")) %>%
    mutate(enroll_site=case_when(groupid == "Bangladesh" ~ 1,
                                 groupid == "Kenya"      ~ 2,
                                 groupid == "Malawi"     ~ 3,
                                 groupid == "Mali"       ~ 4,
                                 groupid == "Pakistan"   ~ 5,
                                 groupid == "Peru"       ~ 6,
                                 groupid == "Gambia"     ~ 7),
           swi_isolate_group=case_when(str_detect(serotype,", S. flexneri") | str_detect(serotype,", S.flexneri") ~ 2,
                                       str_detect(serotype,", S. sonnei")  | str_detect(serotype,", S.sonnei")    ~ 4,
                                       serotype %in% c("not typed","Not typed")                                   ~ 5,
                                       TRUE                                                                       ~ 5), 
           swi_isolate_type=case_when(str_detect(str_remove(serotype,".*,"),"1a") & swi_isolate_group==2 ~ 1,
                                      str_detect(str_remove(serotype,".*,"),"1b") & swi_isolate_group==2 ~ 2,
                                      str_detect(str_remove(serotype,".*,"),"1d") & swi_isolate_group==2 ~ 3,
                                      str_detect(str_remove(serotype,".*,"),"2a") & swi_isolate_group==2 ~ 4,
                                      str_detect(str_remove(serotype,".*,"),"2b") & swi_isolate_group==2 ~ 5,
                                      str_detect(str_remove(serotype,".*,"),"3a") & swi_isolate_group==2 ~ 6,
                                      str_detect(str_remove(serotype,".*,"),"3b") & swi_isolate_group==2 ~ 7,
                                      str_detect(str_remove(serotype,".*,"),"4a") & swi_isolate_group==2 ~ 8,
                                      str_detect(str_remove(serotype,".*,"),"4b") & swi_isolate_group==2 ~ 9,
                                      str_detect(str_remove(serotype,".*,"),"5a") & swi_isolate_group==2 ~ 10,
                                      str_detect(str_remove(serotype,".*,"),"5b") & swi_isolate_group==2 ~ 11,
                                      str_detect(str_remove(serotype,".*,"),"6")  & swi_isolate_group==2 ~ 12,
                                      str_detect(str_remove(serotype,".*,"),"7a") & swi_isolate_group==2 ~ 17,
                                      str_detect(str_remove(serotype,".*,"),"X")  & swi_isolate_group==2 ~ 13,
                                      str_detect(str_remove(serotype,".*,"),"Y")  & swi_isolate_group==2 ~ 14,
                                      swi_isolate_group==2                                               ~ 16)) %>%
    select(pid,ct,sample_type,swi_isolate_group,swi_isolate_type,serotype) %>%
    mutate(source="tac")

  # Combine culture and TAC results and merge in covariates
  shigella_results <- culture_results %>% # combine lab results
    bind_rows(tac_results |> mutate(sens2=1) |> rename(tac_serotype=serotype)) %>%
    bind_rows(tac_results_coinfections |> rename(tac_serotype=serotype)) %>%
    # Merge in screening dysentery status
    left_join(scr_data, by = "pid") |>
    # Now merge in covariates
    left_join(enrollment %>% select(pid,enroll_date,enroll_site,enroll_diar_blood), by = "pid") %>%
    left_join(anthro %>% select(pid, enr_age_months) , by = "pid") %>%
    left_join(severity %>% select(pid, dysentery, episode_start_date,hosp_during_episode,gems_msd:maled_classification), by = "pid") %>%
    # Merge in healthcare seeking weights
    left_join(prop_care |> select(wt_trim,pid), by = "pid") %>% # propensity scores
    left_join(prop_care_7days |> select(wt_trim,pid) |> rename(wt_trim_7days=wt_trim), by = "pid") %>% # propensity scores
    left_join(a_factor %>% select(enroll_site,b_dysentery,a_factor_weight,a_factor_step2,a_factor_step2_sens), by = c("enroll_site","b_dysentery")) %>%
    mutate(# create variables for sensitivity analyses
           sens1=ifelse((source=="tac" & ct<29.8 & sample_type=="Rectal") | sample_type=="Whole Stool",1,NA), #CT<29.8 instead of 31.1
             # NOTE - sens2 created during bind_rows step - excludes secondary infections
           sens3=ifelse(source=="tac" & sample_type=="Rectal",1,NA),
           # Convert % of careseeking at EFGH to weight
           a_factor_step2_wt=100/a_factor_step2,
           # Finalize propensity scores
           ps_wt_final=wt_trim*(1/(a_factor_step2/100)),
           # Finalize propensity scores for sensitivity analysis 
           ps_wt_final_sens=wt_trim_7days*(1/(a_factor_step2_sens/100)),
           # Create vaccineindicators of interest
           quad=ifelse(swi_isolate_group==4 | swi_isolate_type %in% c(4,6,12),1,0),  # Create quadrivalent vaccine target indicator(S. flexneri 2a/3a/6 + S. sonnei)
           biv=ifelse(swi_isolate_group==4 | swi_isolate_type %in% c(4),1,0),  # Create bivalent vaccine target indicator (S. flexneri 2a and S. sonnei)
           quad2=ifelse(swi_isolate_group==4 | swi_isolate_type %in% c(2,4,6),1,0), # additional quadrivalent formulation (S. flexneri 1b/2a/3a and S. sonnei)
           # Cross-protection variables
           quad_cross=ifelse(swi_isolate_type %in% c(1,5,7:11,13),1,0), 
           biv_cross=ifelse(swi_isolate_type %in% c(1,5,7,8,10),1,0),
           quad2_cross=ifelse(swi_isolate_type %in% c(1,5,7:11,13),1,0), 
           # Age group
           agegroup=ifelse(enr_age_months>=6  & enr_age_months<9, 1, 
                    ifelse(enr_age_months>=9  & enr_age_months<12,2,
                    ifelse(enr_age_months>=12 & enr_age_months<18,3,
                    ifelse(enr_age_months>=18 & enr_age_months<24,4,
                    ifelse(enr_age_months>=24 & enr_age_months<36,5,99))))),
           fac_id=substr(pid,1,2),
           fac_id=ifelse(pid==1401167,15,
                  ifelse(pid==7100171,72,fac_id)),month=month(enroll_date),
           year=year(enroll_date),
           month=month(enroll_date),
           mon_yr=paste0(month,"-",year),
           # FOR JEN - Shig 6 by month
           mon_yr_sf6=ifelse(swi_isolate_type==12,mon_yr,NA),
           # Fill in missing data for severity indicators
           gems_shig=ifelse(is.na(gems_shig),"Missing",gems_shig),
           mvs=ifelse(is.na(mvs_classification),"Missing",mvs_classification),
           clark=ifelse(is.na(clark_classification),"Missing",clark_classification),
           maled=ifelse(is.na(maled_classification),"Missing",maled_classification),
           mvs_dys=ifelse(dysentery==1 | mvs_classification %in% c("Moderate","Severe"),1,
                   ifelse(mvs_classification=="Missing",99,0)),
           # Create combined indicators for bivalent and quadrivalent vaccine targets and severity indicators
           gems_msd_biv=ifelse(biv==1,gems_msd,NA),
           mvs_dys_biv=ifelse(biv==1,mvs_dys,NA), 
           gems_msd_quad=ifelse(quad==1,gems_msd,NA),
           mvs_dys_quad=ifelse(quad==1,mvs_dys,NA), 
           gems_msd_quad2=ifelse(quad2==1,gems_msd,NA),
           mvs_dys_quad2=ifelse(quad2==1,mvs_dys,NA)) %>%
    rename(serotype=swi_isolate_type,serogroup=swi_isolate_group) |>
    # Finally - merge in weight accouning for missing TAC results and set culture weight to 1
    left_join(tac_avail |> mutate(source ="tac") |> select(enroll_site,source,month,tac_weight),by=c("enroll_site","month","source")) |>
    mutate(tac_weight=ifelse(is.na(tac_weight),1,tac_weight))

# Run numerators
  # Total number of shigella positives 
  num_shigella_overall_facility <- shigella_results %>%
    select(pid,enroll_site,fac_id,b_dysentery,source,ps_wt_final,ps_wt_final_sens,a_factor_step2_wt,tac_weight,a_factor_weight) %>%
    unique() %>%
    group_by(enroll_site,b_dysentery,source,fac_id) %>%
    summarise(shigella_pos=n(), # crude cases
              shigella_est=sum(tac_weight), # crude cases accounting for missing TAC results
              shigella_est2=sum(tac_weight*a_factor_step2_wt), #crude cases accounting for missing TAC results extroplated to other care seeking
              shigella_weight=sum(ps_wt_final*tac_weight),# care seeking adjusted incidence
              shigella_weight_sens=sum(ps_wt_final_sens*tac_weight), # care seeking adjust incidence - SENSITIVITY ANALYSIS
              shigella_weight_old=sum(a_factor_weight*tac_weight)) 
  
  num_shigella_overall_site <- shigella_results %>%
    select(pid,enroll_site,b_dysentery,source,ps_wt_final,ps_wt_final_sens,a_factor_step2_wt,tac_weight,a_factor_weight) %>%
    unique() %>%
    group_by(enroll_site,b_dysentery,source) %>%
    summarise(shigella_pos=n(), # crude cases
              shigella_est=sum(tac_weight), # crude cases accounting for missing TAC results
              shigella_est2=sum(tac_weight*a_factor_step2_wt), #crude cases accounting for missing TAC results extroplated to other care seeking
              shigella_weight=sum(ps_wt_final*tac_weight),# care seeking adjusted incidence
              shigella_weight_sens=sum(ps_wt_final_sens*tac_weight), # care seeking adjust incidence - SENSITIVITY ANALYSIS
              shigella_weight_old=sum(a_factor_weight*tac_weight)) 
  
  num_shigella_overall_total <- shigella_results %>%
    select(pid,b_dysentery,source,ps_wt_final,ps_wt_final_sens,tac_weight,a_factor_step2_wt,a_factor_weight) %>%
    mutate(enroll_site=8) %>%
    unique() %>%
    group_by(enroll_site,b_dysentery,source) %>%
    summarise(shigella_pos=n(), # crude cases
              shigella_est=sum(tac_weight), # crude cases accounting for missing TAC results
              shigella_est2=sum(tac_weight*a_factor_step2_wt), #crude cases accounting for missing TAC results extroplated to other care seeking
              shigella_weight=sum(ps_wt_final*tac_weight),# care seeking adjusted incidence
              shigella_weight_sens=sum(ps_wt_final_sens*tac_weight), # care seeking adjust incidence - SENSITIVITY ANALYSIS
              shigella_weight_old=sum(a_factor_weight*tac_weight)) 
  
  # SUBGROUPS 2: by serotype/serogroup, vaccine strains, diarrhea severity indicators, hospitalization and severity/vaccines
  # Write functions generating the facility, site and overall number of shigella positives by TAC and culture stratified by dysentery
  # HAVE
  facility_fun=function(a) {
    as.data.frame(shigella_results %>%
                    filter(!is.na(UQ(sym(a)))) %>%
                    select(pid, eval(substitute(a)),enroll_site,b_dysentery,source,a_factor_step2_wt,fac_id,ps_wt_final,ps_wt_final_sens,tac_weight,a_factor_weight) |>
                    unique() |>
                    group_by(enroll_site,.data[[a]],b_dysentery,source,fac_id) |>
                    summarise(shigella_pos=n(), # crude cases
                              shigella_est=sum(tac_weight), # crude cases accounting for missing TAC results
                              shigella_est2=sum(tac_weight*a_factor_step2_wt), #crude cases accounting for missing TAC results extroplated to other care seeking
                              shigella_weight=sum(ps_wt_final*tac_weight),# care seeking adjusted incidence
                              shigella_weight_sens=sum(ps_wt_final_sens*tac_weight), # care seeking adjust incidence - SENSITIVITY ANALYSIS
                              shigella_weight_old=sum(a_factor_weight*tac_weight))) 
  }
  site_fun=function(a) {
    as.data.frame(shigella_results %>%
                    filter(!is.na(UQ(sym(a)))) %>%
                    select(pid,eval(substitute(a)),enroll_site,b_dysentery,source,a_factor_step2_wt,ps_wt_final,ps_wt_final_sens,tac_weight,a_factor_weight) %>%
                    unique() %>%
                    group_by(enroll_site,.data[[a]],b_dysentery,source) %>%
                    summarise(shigella_pos=n(), # crude cases
                              shigella_est=sum(tac_weight), # crude cases accounting for missing TAC results
                              shigella_est2=sum(tac_weight*a_factor_step2_wt), #crude cases accounting for missing TAC results extroplated to other care seeking
                              shigella_weight=sum(ps_wt_final*tac_weight),# care seeking adjusted incidence
                              shigella_weight_sens=sum(ps_wt_final_sens*tac_weight), # care seeking adjust incidence - SENSITIVITY ANALYSIS
                              shigella_weight_old=sum(a_factor_weight*tac_weight)))
  }
  total_fun=function(a) {
    as.data.frame(shigella_results %>%
                    filter(!is.na(UQ(sym(a)))) %>%
                    select(pid,eval(substitute(a)),b_dysentery,source,a_factor_step2_wt,ps_wt_final,ps_wt_final_sens,tac_weight,a_factor_weight) %>%
                    mutate(enroll_site=8) %>%
                    unique() %>%
                    group_by(enroll_site,.data[[a]],b_dysentery,source) %>%
                    summarise(shigella_pos=n(), # crude cases
                              shigella_est=sum(tac_weight), # crude cases accounting for missing TAC results
                              shigella_est2=sum(tac_weight*a_factor_step2_wt), #crude cases accounting for missing TAC results extroplated to other care seeking
                              shigella_weight=sum(ps_wt_final*tac_weight),# care seeking adjusted incidence
                              shigella_weight_sens=sum(ps_wt_final_sens*tac_weight), # care seeking adjust incidence - SENSITIVITY ANALYSIS
                              shigella_weight_old=sum(a_factor_weight*tac_weight)))
  }
  
  # wrapper function
  subgroups2_fun = function(a, b) {
    if(b=="facility") {
      facility_fun(a)
    } else if(b=="site") {
      site_fun(a)
    } else if(b=="total") {
      total_fun(a)
    } else {print("Incorrect level selected")}
  }

  # specifying objects needed to run through all level/subgroup combinations
  fun_subgroup = c(# Vaccine subgroups
                    "serogroup", "serotype", "dys_serotype","quad", "quad2","biv", 
                    # Age
                    "agegroup", 
                    # Diarrhea severity
                    "gems_msd", "gems_shig", "mvs", "mvs_dys", "clark","maled","hosp_during_episode", "dysentery",
                    # Seasonality
                    "mon_yr",
                    # Vaccine targets by GEMS/MVS severity
                    "mvs_dys_biv", "mvs_dys_quad", "gems_msd_biv", "gems_msd_quad", "mvs_dys_quad2", "gems_msd_quad2",
                    # sensitivity analyses
                    "sens1","sens2","sens3")
  fun_level = c("facility", "site", "total")
  grid = expand.grid(fun_subgroup, fun_level, stringsAsFactors = FALSE)
  
  # run wrapper function through each combination
  num_subgroup = setNames(do.call(mapply, c("subgroups2_fun", unname(as.list(grid)))), paste0("num_shigella_", grid$Var1, "_", grid$Var2) )
  
  # Rbind list plus overall non-stratified
  numerator = bind_rows(num_subgroup) %>%
    bind_rows(num_shigella_overall_facility, num_shigella_overall_site, num_shigella_overall_total) %>%
    # num_shigella_season_facility, num_shigella_season_site, num_shigella_season_total) %>%
    # re-create month and year variables
    mutate(month = str_split_i(mon_yr, "-", 1),
           year = str_split_i(mon_yr, "-", 2))
  
  # care final numerator dataset
  numerator_merge <- numerator %>%
    # Pull out # positive among those with valid sample results
    select(-shigella_weight,-shigella_weight_old,-shigella_est,-shigella_est2,-shigella_weight_sens) %>%
    pivot_wider(names_from=source, values_from=shigella_pos, values_fill = 0)  %>%
    rename(culture_pos=culture,tac_pos=tac) %>%
    # Pull out # estimated cases among enrolled
    left_join(numerator |> select(-shigella_weight,-shigella_weight_old,-shigella_pos,-shigella_est2,-shigella_weight_sens) |>
                           pivot_wider(names_from=source, values_from=shigella_est, values_fill = 0),
              by = c("enroll_site","fac_id","sens1","sens2","sens3","b_dysentery","dysentery","agegroup","serogroup","dys_serotype", "serotype","quad",
                     "quad2","biv","hosp_during_episode","gems_msd","gems_msd_biv","gems_msd_quad","gems_shig","mvs","mvs_dys","mvs_dys_biv","mvs_dys_quad",
                     "mvs_dys_quad2","gems_msd_quad2","month","year", "clark","maled","mon_yr")) %>%
    rename(culture_est=culture,tac_est=tac) %>%
    # Pull out estimate number of cases that sought care using EFGH weight
    left_join(numerator %>% select(-shigella_pos,-shigella_weight,-shigella_weight_old,-shigella_est,-shigella_weight_sens) %>% 
                pivot_wider(names_from=source, values_from=shigella_est2, values_fill = 0),
              by = c("enroll_site","fac_id","sens1","sens2","sens3","b_dysentery","dysentery","agegroup","serogroup","dys_serotype", "serotype","quad","quad2","biv","hosp_during_episode","gems_msd","gems_msd_biv","gems_msd_quad","gems_shig",
                     "mvs","mvs_dys","mvs_dys_biv","mvs_dys_quad","mvs_dys_quad2","gems_msd_quad2","month","year", "clark","maled","mon_yr")) %>%
    rename(culture_est2=culture,tac_est2=tac) %>%
    # Pull out community diarrhea cases using propensity scores and old weights
    left_join(numerator %>% select(-shigella_pos,-shigella_weight_old,-shigella_est,-shigella_est2,-shigella_weight_sens) %>% 
                            pivot_wider(names_from=source, values_from=shigella_weight, values_fill = 0),
              by = c("enroll_site","fac_id","sens1","sens2","sens3","b_dysentery","dysentery","agegroup","serogroup","dys_serotype", "serotype","quad","quad2","biv","hosp_during_episode","gems_msd","gems_msd_biv","gems_msd_quad","gems_shig",
                     "mvs","mvs_dys","mvs_dys_biv","mvs_dys_quad","mvs_dys_quad2","gems_msd_quad2","month","year", "clark","maled","mon_yr")) %>%
    rename(culture_wt=culture,tac_wt=tac) %>%
    left_join(numerator %>% select(-shigella_pos,-shigella_weight,-shigella_est,-shigella_est2,-shigella_weight_sens) %>% 
                            pivot_wider(names_from=source, values_from=shigella_weight_old, values_fill = 0),
              by = c("enroll_site","fac_id","sens1","sens2","sens3","b_dysentery","dysentery","agegroup","serogroup","dys_serotype", "serotype","quad","quad2","biv","hosp_during_episode","gems_msd","gems_msd_biv","gems_msd_quad","gems_shig",
                     "mvs","mvs_dys","mvs_dys_biv","mvs_dys_quad","mvs_dys_quad2","gems_msd_quad2","month","year", "clark","maled","mon_yr")) %>%
    rename(culture_wt_old=culture,tac_wt_old=tac) %>%
    left_join(numerator %>% select(-shigella_pos,-shigella_weight,-shigella_est,-shigella_est2,-shigella_weight_old) %>% 
                pivot_wider(names_from=source, values_from=shigella_weight_sens, values_fill = 0),
              by = c("enroll_site","fac_id","sens1","sens2","sens3","b_dysentery","dysentery","agegroup","serogroup","dys_serotype", "serotype","quad","quad2","biv","hosp_during_episode","gems_msd","gems_msd_biv","gems_msd_quad","gems_shig",
                     "mvs","mvs_dys","mvs_dys_biv","mvs_dys_quad","mvs_dys_quad2","gems_msd_quad2","month","year", "clark","maled","mon_yr")) %>%
    rename(culture_wt_sens=culture,tac_wt_sens=tac)

write_rds(numerator_merge, paste0("./Last Step Datasets/numerator.Rds")) 

# Export TAC and culture results to look at serotypes
serotype_dist <- shigella_results |>
  select(pid, enroll_site,source, serogroup, serotype,dys_serotype,sens2,quad,quad2,biv,quad_cross,quad2_cross,biv_cross) |>
  distinct() |>
  mutate(serotype=factor(serotype,levels = c(1:17),
                                  labels = c("1a","1b","1d","2a","2b","3a","3b","4a","4b","5a","5b","6","X","Y","Non-typeable","Other","7a")),
         serogroup=factor(serogroup,levels = c(1:5),
                                  labels = c("S. dysenteriae","S. flexneri","S. boydii","S. sonnei","Undetermined")))

write_rds(serotype_dist, paste0("./Last Step Datasets/serotype_dist.Rds")) 

###############################################
######### incidence of all diarrhea ###########
###############################################
# Combine culture and TAC results and merge in covariates
diarrhea_results <- enrollment |>
  select(pid,enroll_date,enroll_site,enroll_diar_blood) |>
  # Merge in screening dysentery status
  left_join(scr_data, by = "pid") |>
  # Now merge in covariates
  left_join(anthro %>% select(pid, enr_age_months) , by = "pid") %>%
  left_join(severity %>% select(pid, dysentery, episode_start_date,hosp_during_episode,gems_msd:maled_classification), by = "pid") %>%
  # Merge in healthcare seeking weights
  left_join(prop_care |> select(wt_trim,pid), by = "pid") %>% # propensity scores
  left_join(prop_care_7days |> select(wt_trim,pid) |> rename(wt_trim_7days=wt_trim), by = "pid") %>% # propensity scores
  left_join(a_factor %>% select(enroll_site,b_dysentery,a_factor_weight,a_factor_step2,a_factor_step2_sens), by = c("enroll_site","b_dysentery")) %>%
  mutate(# create variables for sensitivity analyses
    # Convert % of careseeking at EFGH to weight
    a_factor_step2_wt=100/a_factor_step2,
    # Finalize propensity scores
    ps_wt_final=wt_trim*(1/(a_factor_step2/100)),
    # Finalize propensity scores for sensitivity analysis 
    ps_wt_final_sens=wt_trim_7days*(1/(a_factor_step2_sens/100)),
    # Age group
    agegroup=ifelse(enr_age_months>=6  & enr_age_months<9, 1, 
             ifelse(enr_age_months>=9  & enr_age_months<12,2,
             ifelse(enr_age_months>=12 & enr_age_months<18,3,
             ifelse(enr_age_months>=18 & enr_age_months<24,4,
             ifelse(enr_age_months>=24 & enr_age_months<36,5,99))))),
    # Facility ID
    fac_id=substr(pid,1,2),
    fac_id=ifelse(pid==1401167,15,
           ifelse(pid==7100171,72,fac_id)),month=month(enroll_date),
    # Seasonality / study month
    year=year(enroll_date),
    month=month(enroll_date),
    mon_yr=paste0(month,"-",year),
    # Fill in missing data for severity indicators
    gems_shig=ifelse(is.na(gems_shig),"Missing",gems_shig),
    mvs=ifelse(is.na(mvs_classification),"Missing",mvs_classification),
    clark=ifelse(is.na(clark_classification),"Missing",clark_classification),
    maled=ifelse(is.na(maled_classification),"Missing",maled_classification),
    mvs_dys=ifelse(dysentery==1 | mvs_classification %in% c("Moderate","Severe"),1,
                   ifelse(mvs_classification=="Missing",99,0)))

# Run numerators
# Total number of diarrhea positives 
num_diarrhea_overall_facility <- diarrhea_results %>%
  select(pid,enroll_site,fac_id,b_dysentery,ps_wt_final,a_factor_step2_wt) %>%
  unique() %>%
  group_by(enroll_site,b_dysentery,fac_id) %>%
  summarise(diarrhea_pos=n(), # crude cases
            diarrhea_est=sum(a_factor_step2_wt), #crude cases extrapolated to other care seeking
            diarrhea_weight=sum(ps_wt_final))# care seeking adjusted incidence
            

num_diarrhea_overall_site <- diarrhea_results %>%
  select(pid,enroll_site,b_dysentery,ps_wt_final,a_factor_step2_wt) %>%
  unique() %>%
  group_by(enroll_site,b_dysentery) %>%
  summarise(diarrhea_pos=n(), # crude cases
            diarrhea_est=sum(a_factor_step2_wt), #crude cases extroplated to other care seeking
            diarrhea_weight=sum(ps_wt_final))# care seeking adjusted incidence
           
num_diarrhea_overall_total <- diarrhea_results %>%
  select(pid,b_dysentery,ps_wt_final,a_factor_step2_wt,a_factor_weight) %>%
  mutate(enroll_site=8) %>%
  unique() %>%
  group_by(enroll_site,b_dysentery) %>%
  summarise(diarrhea_pos=n(), # crude cases
            diarrhea_est=sum(a_factor_step2_wt), #crude cases extroplated to other care seeking
            diarrhea_weight=sum(ps_wt_final))# care seeking adjusted incidence
   
# Rbind list plus overall non-stratified
diarrhea_numerator <- bind_rows(num_diarrhea_overall_facility, num_diarrhea_overall_site, num_diarrhea_overall_total)
    
write_rds(diarrhea_numerator, paste0("./Last Step Datasets/diarrhea_numerator.Rds")) 
