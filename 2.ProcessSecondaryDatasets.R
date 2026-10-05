###################################################
##AUTHOR: Sean Galagan 
##DATE STARTED: 18NOV2022
##PURPOSE: Builds Last Step Datasets for Analysis
## NOTES:
## LAST UPDATED: 30JAN2026
###################################################

## ENSURE NECESSARY PACKAGES ARE INSTALLED
packages = c("ezknitr", "tidyverse", "haven", "knitr", "scales", "lubridate", "conflicted",
             "janitor", "consort", "anthro", "rmarkdown", "zoo", "ggrepel","gtsummary", "flextable")

## Now load or install & load all
package.check <- lapply( packages, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
}
)

# build quick function
`%notin%` <- Negate(`%in%`)

## Load necessary data files & split by country
today <- as.Date((paste0(format(Sys.Date(),"%Y%m%d"))),format="%Y%m%d")

# Declare conflict preferences
conflicts_prefer(dplyr::filter)

######LOAD DATA#####
# PE/HUS
pehus <- readRDS("./Data/PE_HUS_Final.Rds") |>
  filter(pop_country!=5 | (pop_cluster_id!=7 & pop_country==5))
cluster_comp <- readRDS("./Data/PE_ClusterCompletion.Rds") |>
  mutate(pop_cluster_id=clcom_cluster_id) 
careseek <- readRDS("./Data/PE_HUS_Careseek.Rds") %>% 
  filter(!is.na(hus_careseek_efgh)) %>%
  inner_join(pehus %>% select(pop_country,pop_hh_id,link_key,pop_country,hus_blood,hus_consent,re_consent, pop_visit_dt, hus_diarrh_start_dt, hus_diarrh_start, hus_diarrh_start_days,starts_with("hus_care")), by = c("pop_hh_id","link_key")) 
# Clusters
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
  rename(clcom_country=country) 
# DCS screening, enrollment and lab testing 
prescreening <- readRDS("./Data/DCS_01_prescreening.Rds")
screening <- readRDS("./Data/DCS_02_screening.Rds")
preenrollment <- readRDS("./Data/DCS_03_preenrollment.Rds")
enrollment <- readRDS("./Data/DCS_04_enrollment.Rds") %>%
  mutate(pid=as.numeric(pid)) 
rectal_swab_results <- readRDS("./Data/DCS_16_rectal_swab_results.Rds") %>%
  mutate(pid=as.numeric(pid))%>%
  inner_join(enrollment %>% select(pid), by = "pid")
whole_stool_results <- readRDS("./Data/DCS_17_whole_stool_results.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid")
specimen_accession <- readRDS("./Data/DCS_15_specimen_accession.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid")
stool_collection <- readRDS("./Data/DCS_13_stool_collection.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid")
tac <- readRDS("./Data/DCS_tac_automated.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid")
tac_long <- readRDS("./Data/DCS_tac_long.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid")
severity <- readRDS("./Data/DCS_severity_scores.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid") 
anthro <- readRDS("./Data/DCS_anthro.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid") 
wealth <- readRDS("./Data/DCS_wealth_index.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid")
derived <- readRDS("./Data/DCS_derived_variables.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid")
visit_tracker <- readRDS("./Data/DCS_visit_tracker.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid")
closeout <- readRDS("./Data/DCS_11_close_out.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid")
mortality <- readRDS("./Data/DCS_10a_mortality.Rds") %>%
  inner_join(enrollment %>% select(pid), by = "pid")
cod <- readRDS("./Data/cause_of_death.Rds") 
shig <- readRDS("./Data/DCS_shigella.Rds") 
  
#################################################
######### CREATE SUMMARY FIGURE DATASET #########
#################################################
# Build dataset for screening and enrollment 
all_screening <- prescreening |>
  left_join(screening, by="serial_id") |>
  left_join(preenrollment, by="serial_id") |>
  mutate(country=case_when(redcap_data_access_group=="bangladesh_team" ~ "Bangladesh",
                           redcap_data_access_group=="kenya_team"      ~ "Kenya",
                           redcap_data_access_group=="malawi_team"     ~ "Malawi",
                           redcap_data_access_group=="mali_team"       ~ "Mali",
                           redcap_data_access_group=="pakistan_team"   ~ "Pakistan",
                           redcap_data_access_group=="peru_team"       ~ "Peru",
                           redcap_data_access_group=="the_gambia_team" ~ "The Gambia"),
         id_check=substr(paste0(serial_id),1,1),
         prescreen_age=ifelse(((pscr_age>=0 & pscr_age<6) | (pscr_age>=36 & pscr_age<=99) | pscr_age_known==0) & (!is.na(pscr_age) | !is.na(pscr_age_known)),0,1),
         prescreen_pres=ifelse(pscr_pres2==2 & !is.na(pscr_pres2),0,1),
         prescreen_alrenr=ifelse(pscr_already_enr==1 & !is.na(pscr_already_enr),0,1),
         screen_consent=ifelse((scr_careg_yn==0 & !is.na(scr_careg_yn)) | (scr_consent==0 & !is.na(scr_consent)),0,1),
         prescreen_catchment=ifelse(scr_verify_cache2 %in% c(0,2) & !is.na(scr_verify_cache2) & scr_consent %notin% c(1,2),0,1),         
         prescreen_noreason=ifelse(scr_consent %notin% c(1,2) & pscr_diar==1 & prescreen_age==1 & prescreen_pres==1 & prescreen_alrenr==1 & screen_consent==1 & prescreen_catchment==1,0,1),
         eligible_4hours=case_when(scr_elig_yn==1 & (pre_enroll!=1 | is.na(pre_enroll)) & scr_diartime==1 ~ 0,
                                   scr_elig_yn==1 ~ 1),
         eligible_enrollcap=case_when(scr_elig_yn==1 & (pre_enroll!=1 | is.na(pre_enroll)) & scr_diartime==0 & scr_enrollcap==1 ~ 0,
                                      scr_elig_yn==1 ~ 1),
         eligible_refer=case_when(scr_elig_yn==1 & (pre_enroll!=1 | is.na(pre_enroll)) & scr_diartime==0 & scr_enrollcap==0 & scr_refer==2 ~ 0,
                                  scr_elig_yn==1 ~ 1),
         eligible_consent=case_when(scr_elig_yn==1 & (pre_enroll!=1 | is.na(pre_enroll)) & scr_diartime==0 & scr_enrollcap==0 & scr_refer %in% c(1,3) & (pre_consent==0 | pre_enroll_whynot==1) ~ 0,
                                    scr_elig_yn==1 ~ 1),
         eligible_other=case_when(scr_elig_yn==1 & (pre_enroll!=1 | is.na(pre_enroll)) & scr_diartime==0 & scr_enrollcap==0 & scr_refer %in% c(1,3) & pre_consent==1 & pre_enroll_whynot %in% c(2,3,77) ~ 0,
                                  scr_elig_yn==1 ~ 1),
         screen_nodiar=case_when(scr_elig_diar==1 | scr_elig_dhours==1 ~ 1,
                              scr_elig_diar==0 & scr_elig_dhours==0 ~ 0),
         # HARDCODE BANGLADESH CORRECTIONS
         prescreen_catchment=ifelse(country=="Bangladesh" & prescreen_noreason==0,0,prescreen_catchment),
         prescreen_noreason=ifelse(country=="Bangladesh" & prescreen_noreason==0,1,prescreen_noreason))

# Build data for retention
retention <- visit_tracker |>
  left_join(closeout |> select(-sc_date), by = c("enroll_site","pid")) |>
  left_join(mortality, by = c("enroll_site","pid")) |>
  mutate(country=case_when(enroll_site==1 ~ "Bangladesh",
                           enroll_site==2 ~ "Kenya",
                           enroll_site==3 ~ "Malawi",
                           enroll_site==4 ~ "Mali",
                           enroll_site==5 ~ "Pakistan",
                           enroll_site==6 ~ "Peru",
                           enroll_site==7 ~ "The Gambia"), 
         noattend_week4=case_when(week4_attended==1                                                                                ~ 0, # attended visit
                                  (week4_attended==0 & month3_attended==1) |          
                                  (week4_attended==0 & month3_attended==0 & sc_reason %in% c(2,3,6) & sc_date>week4_hardclose_date) |
                                  (week4_attended==0 & month3_attended==0 & sc_reason==5 & mort_date>week4_hardclose_date) |
                                  (is.na(sc_reason) & week4_attended==0 & month3_attended==0)                                      ~ 1, # missed visit, includes withdrawls/deaths that occurred after end of visit window as well as missed visits where the participant has yet to complete the month 3 visit
                                  (week4_attended==0 & month3_attended==0 & sc_reason==4) | pid==6501209                           ~ 2, # LTFU (PID is temporary)
                                  week4_attended==0 & month3_attended==0 & sc_reason==5 & mort_date<=week4_hardclose_date          ~ 3, # died
                                  week4_attended==0 & month3_attended==0 & sc_reason==2 & sc_date<=week4_hardclose_date            ~ 4, # Voluntary withdrawl
                                  week4_attended==0 & month3_attended==0 & sc_reason %in% c(3,6) & sc_date<=week4_hardclose_date   ~ 5), # Other reason withdrew (eligibiliy revoked or other reason)
         noattend_month3=case_when(month3_attended==1                                                                              ~ 0, # attended visit
                                   (month3_attended==0 & noattend_week4 %in% c(0,1) & sc_reason==4) | pid %in% c(3104190,5300833)  ~ 1, # LTFU
                                   month3_attended==0 & noattend_week4 %in% c(0,1) & sc_reason==5                                  ~ 2, # Died
                                   month3_attended==0 & noattend_week4 %in% c(0,1) & sc_reason==2                                  ~ 3, # Voluntary withdrawl
                                   month3_attended==0 & noattend_week4 %in% c(0,1) & sc_reason %in% c(3,6)                         ~ 4, # Other reason withdrew (eligibiliy revoked or other reason)
                                   month3_attended==0 & noattend_week4 %in% c(0,1) & is.na(sc_reason)                              ~ 5)) # Visit pending


# prescreened
prescreened_sites <- all_screening %>%
  group_by(country) %>%
  summarize(frequency=n()) %>%
  mutate(indicator="n prescreened")

prescreened <- prescreened_sites |>
  group_by(indicator) |>
  summarize(frequency=sum(frequency)) |>
  mutate(country="Total") %>%
  bind_rows(prescreened_sites)

# Screened
screened_sites <- all_screening |>
  filter(scr_consent %in% c(1,2)) |>
  group_by(country) |>
  summarize(frequency=n()) |>
  mutate(indicator="n screened")

screened <- screened_sites |>
  group_by(indicator) |>
  summarize(frequency=sum(frequency)) |>
  mutate(country="Total") %>%
  bind_rows(screened_sites)

# Reasons not screened
# Write functions generating the facility, site and overall number of shigella positives by TAC and culture stratified by dysentery
prescreen_site=function(a) {
  as.data.frame(all_screening %>%
                  group_by(country,.data[[a]]) %>%
                  summarise(frequency=n()) %>%
                  mutate(percent = (frequency / sum(frequency))*100) %>%
                  filter(UQ(sym(a))==0))
}
prescreen_total=function(a) {
  as.data.frame(all_screening %>%
                  group_by(.data[[a]]) %>%
                  summarise(frequency=n()) %>%
                  mutate(percent = (frequency / sum(frequency))*100) %>%
                  filter(UQ(sym(a))==0))
}

# wrapper function
prescreen_fun = function(a, b) {
  if(b=="site") {
    prescreen_site(a)
  } else if(b=="total") {
    prescreen_total(a)
  } else {print("Incorrect level selected")}
}

# specifying objects needed to run through all level/subgroup combinations
prescreen_fun_subgroup = c("pscr_diar","prescreen_age","prescreen_pres","prescreen_alrenr","screen_consent","prescreen_catchment","prescreen_noreason")
prescreen_fun_level = c("site", "total")
grid = expand.grid(prescreen_fun_subgroup, prescreen_fun_level, stringsAsFactors = FALSE)

# run wrapper function through each combination
prescreen_num_subgroup = setNames(do.call(mapply, c("prescreen_fun", unname(as.list(grid)))), paste0("prescreen_", grid$Var1, "_", grid$Var2))

# Rbind list plus overall non-stratified
noscreenreasons = bind_rows(prescreen_num_subgroup) %>%
  mutate(country=ifelse(is.na(country),"Total",country),
         indicator=case_when(pscr_diar==0 ~ "Not screened: Diarrhea/dysentery/gastroenteritis not a complaint of the child",
                             prescreen_age==0 ~ "Not screened: Child <6 months or 36+ months of age, or age is unknown",
                             prescreen_pres==0 ~ "Not screened: Child no longer present in the facility",
                             prescreen_alrenr==0 ~ "Not screened: Already enrolled in EFGH follow-up",
                             screen_consent==0 ~ "Not screened: Refused consent to screening or adult not available",
                             prescreen_catchment==0 ~ "Not screened: Lives outside of the catchment area",
                             prescreen_noreason==0 ~ "Not screened: Not screened, reason unknown")) |>
  select(country,frequency,percent,indicator)

# Eligible for enrollment
eligible_sites <- all_screening |>
  filter(scr_elig_yn==1) |>
  group_by(country) |>
  summarize(frequency=n()) |>
  mutate(indicator="n eligible")

eligible <- eligible_sites |>
  group_by(indicator) |>
  summarize(frequency=sum(frequency)) |>
  mutate(country="Total") %>%
  bind_rows(eligible_sites)

# wrapper function
# Write functions generating the facility, site and overall number of shigella positives by TAC and culture stratified by dysentery
screen_site=function(a) {
  as.data.frame(all_screening %>%
                  filter(scr_consent %in% c(1,2)) |>
                  group_by(country,.data[[a]]) %>%
                  summarise(frequency=n()) %>%
                  mutate(percent = (frequency / sum(frequency))*100) %>%
                  filter(UQ(sym(a))==0))
}
screen_total=function(a) {
  as.data.frame(all_screening %>%
                  filter(scr_consent %in% c(1,2)) |>
                  group_by(.data[[a]]) %>%
                  summarise(frequency=n()) %>%
                  mutate(percent = (frequency / sum(frequency))*100) %>%
                  filter(UQ(sym(a))==0))
}

# wrapper function
screen_fun = function(a, b) {
  if(b=="site") {
    screen_site(a)
  } else if(b=="total") {
    screen_total(a)
  } else {print("Incorrect level selected")}
}

# specifying objects needed to run through all level/subgroup combinations
screen_fun_subgroup = c("screen_nodiar","scr_elig_consent","scr_elig_area","scr_elig_age","scr_elig_ddays","scr_elig_months","scr_elig_study")
screen_fun_level = c("site", "total")
screen_grid = expand.grid(screen_fun_subgroup, screen_fun_level, stringsAsFactors = FALSE)

# run wrapper function through each combination
screen_num_subgroup = setNames(do.call(mapply, c("screen_fun", unname(as.list(screen_grid)))), paste0("eligible_", screen_grid$Var1, "_", screen_grid$Var2))

# Rbind list plus overall non-stratified
noteligible_reasons = bind_rows(screen_num_subgroup) %>%
  mutate(country=ifelse(is.na(country),"Total",country),
         indicator=case_when(screen_nodiar==0 ~ "Not eligible: no diarrhea",
                             scr_elig_consent==0 ~ "Not eligible: consent not possible",
                             scr_elig_area==0 ~ "Not eligible: not in catchment area",
                             scr_elig_age==0 ~ "Not eligible: out of age range",
                             scr_elig_ddays==0 ~ "Not eligible: 7+ days of diarrhea",
                             scr_elig_months==0 ~ "Not eligible: will leave area within 3 months",
                             scr_elig_study==0 ~ "Not eligible: other study")) |>
  select(country,frequency,percent,indicator)

# Enrolled
enrolled_sites <- all_screening |>
  filter(pre_enroll==1) |>
  group_by(country) |>
  summarize(frequency=n()) |>
  mutate(indicator="n enrolled")

enrolled <- enrolled_sites |>
  group_by(indicator) |>
  summarize(frequency=sum(frequency)) |>
  mutate(country="Total") %>%
  bind_rows(enrolled_sites)

# Reasons not enrolled
  # wrapper function - reasons not enrolled
  eligible_site=function(a) {
    as.data.frame(all_screening %>%
                    filter(scr_elig_yn==1) |>
                    group_by(country,.data[[a]]) %>%
                    summarise(frequency=n()) %>%
                    mutate(percent = (frequency / sum(frequency))*100) %>%
                    filter(UQ(sym(a))==0))
  }
  eligible_total=function(a) {
    as.data.frame(all_screening %>%
                    filter(scr_elig_yn==1) |>
                    group_by(.data[[a]]) %>%
                    summarise(frequency=n()) %>%
                    mutate(percent = (frequency / sum(frequency))*100) %>%
                    filter(UQ(sym(a))==0))
  }
  
  # wrapper function
  eligible_fun = function(a, b) {
    if(b=="site") {
      eligible_site(a)
    } else if(b=="total") {
      eligible_total(a)
    } else {print("Incorrect level selected")}
  }
  # specifying objects needed to run through all level/subgroup combinations
  eligible_fun_subgroup = c("eligible_4hours","eligible_enrollcap","eligible_refer","eligible_consent","eligible_other")
  eligible_fun_level = c("site", "total")
  eligible_grid = expand.grid(eligible_fun_subgroup, eligible_fun_level, stringsAsFactors = FALSE)
  
  # run wrapper function through each combination
  eligible_num_subgroup = setNames(do.call(mapply, c("eligible_fun", unname(as.list(eligible_grid)))), paste0("eligible_", eligible_grid$Var1, "_", eligible_grid$Var2))
  
  # Rbind list plus overall non-stratified
  notenroll_reasons = bind_rows(eligible_num_subgroup) %>%
    mutate(country=ifelse(is.na(country),"Total",country),
           indicator=case_when(eligible_4hours==0 ~ "Not enrolled: 4 hours have passed since child presented to facility",
                               eligible_enrollcap==0 ~ "Not enrolled: enrollment cap met",
                               eligible_refer==0 ~ "Not enrolled: Child referred to non-study facility",
                               eligible_consent==0 ~ "Not enrolled: Caregiver refused or withdrew consent",
                               eligible_other==0 ~ "Not enrolled: Other reason")) |>
    select(country,frequency,percent,indicator)

# attended week four visit
week_four_sites <- retention |>
  group_by(country, noattend_week4) |>
  summarize(frequency=n()) |>
  mutate(percent = (frequency / sum(frequency))*100) 
  
week_four <- retention |>
  group_by(noattend_week4) |>
  summarize(frequency=n()) |>
  mutate(percent= (frequency / sum(frequency))*100,
         country="Total") |>
  bind_rows(week_four_sites) |>
  mutate(indicator=case_when(noattend_week4==0 ~ "Week 4: Attended visit",
                             noattend_week4==1 ~ "Week 4: Missed visit",
                             noattend_week4==2 ~ "Week 4: LTFU",
                             noattend_week4==3 ~ "Week 4: Died",
                             noattend_week4==4 ~ "Week 4: Voluntary withdrawl",
                             noattend_week4==5 ~ "Week 4: Other reason withdrew")) |>
  select(-noattend_week4)

# attended week four visit
month_three_sites <- retention |>
  filter(!is.na(noattend_month3)) |>
  group_by(country, noattend_month3) |>
  summarize(frequency=n()) |>
  mutate(percent = (frequency / sum(frequency))*100) 

month_three <- retention |>
  filter(!is.na(noattend_month3)) |>
  group_by(noattend_month3) |>
  summarize(frequency=n()) |>
  mutate(percent= (frequency / sum(frequency))*100,
         country="Total") |>
  bind_rows(month_three_sites) |>
  mutate(indicator=case_when(noattend_month3==0 ~ "Month 3: Attended visit",
                             noattend_month3==1 ~ "Month 3: LTFU",
                             noattend_month3==2 ~ "Month 3: Died",
                             noattend_month3==3 ~ "Month 3: Voluntary withdrawl",
                             noattend_month3==4 ~ "Month 3: Other reason withdrew",
                             noattend_month3==5 ~ "Month 3: Visit pendiong")) |>
  select(-noattend_month3)

# Combine into a single dataset
summary_final <- bind_rows(prescreened,screened,noscreenreasons,eligible,noteligible_reasons,enrolled,notenroll_reasons, week_four, month_three)      %>%
  mutate(value=ifelse(is.na(percent),paste0(prettyNum(frequency, big.mark=",", scientific = FALSE)),
                                                paste0(prettyNum(frequency, big.mark=",", scientific = FALSE), " (", sprintf("%.1f", round(percent,1)), ")"))) %>%
  select(-c( frequency, percent)) |>
  pivot_wider(names_from = "country", values_from = "value") |>
  mutate(Bangladesh=ifelse(is.na(Bangladesh), "0 (0.0)",Bangladesh),
         Kenya=ifelse(is.na(Kenya), "0 (0.0)",Kenya),
         Malawi=ifelse(is.na(Malawi), "0 (0.0)",Malawi),
         Mali=ifelse(is.na(Mali), "0 (0.0)",Mali),
         Pakistan=ifelse(is.na(Pakistan), "0 (0.0)",Pakistan),
         Peru=ifelse(is.na(Peru), "0 (0.0)",Peru),
         `The Gambia`=ifelse(is.na(`The Gambia`), "0 (0.0)",`The Gambia`))

# Export last step dataset
write_rds(summary_final, paste0("./Last Step Datasets/enrollment_summary.Rds")) 

################################################
############## CREATE DCS TABLE 1 ##############
################################################
# Merge and create dataset
enroll_data <- enrollment %>%
  left_join(severity, by="pid") %>%
  left_join(anthro, by="pid") %>%
  left_join(wealth, by="pid") %>%
  left_join(derived, by="pid") %>%
  left_join(screening %>% select (serial_id, scr_sex, scr_date,scr_diar_when), by = "serial_id") %>%
  mutate(mvs_or_dysentery_ind = ifelse(mvs_or_dysentery=="Mild with no dysentery" & !is.na(mvs_or_dysentery), 0,
                                ifelse(mvs_or_dysentery=="Moderate/Severe or Dysentery" & !is.na(mvs_or_dysentery), 1, NA)),
         country=case_when(enroll_site==1 ~ "Bangladesh",
                           enroll_site==2 ~ "Kenya",
                           enroll_site==3 ~ "Malawi",
                           enroll_site==4 ~ "Mali",
                           enroll_site==5 ~ "Pakistan",
                           enroll_site==6 ~ "Peru",
                           enroll_site==7 ~ "The Gambia"),
         agegroup=factor(case_when(enr_age_months>=5  & enr_age_months <9  ~ 1,
                                   enr_age_months>=9  & enr_age_months <12 ~ 2,
                                   enr_age_months>=12 & enr_age_months <18 ~ 3,
                                   enr_age_months>=18 & enr_age_months <24 ~ 4,
                                   enr_age_months>=24 & enr_age_months <36 ~ 5), levels=1:5, labels=c("6-8","9-11","12-17","18-23","24-35")),
         mat_edu=factor(case_when(enroll_cg_moth_ed %in% c(1,2)   ~ 1,
                                  enroll_cg_moth_ed==7            ~ 2,
                                  enroll_cg_moth_ed %in% c(3,4,5) ~ 3,
                                  enroll_cg_moth_ed %in% c(6,8)   ~ 4), levels=1:4, labels=c("Less than primary school","Koranic school only","Primary school or greater","Unknown or declined")),
         sex=factor(scr_sex, levels=1:2, labels=c("Male", "Female sex: n (%)")),
         sex2=case_when(scr_sex==1 ~ 0,
                        scr_sex==2 ~ 1),
         wealth_quint=as.numeric(final_quintile),
         care_employ=factor(case_when(enroll_cg_work %in% c(1,2,3) ~ 1,
                                      enroll_cg_work %in% c(4,5,6,7,8) ~ 2,
                                      enroll_cg_work %in% c(9,10) ~ 3), levels=1:3, labels=c("Not employed","Employed","Unknown or other")),
         underweight_ordered=factor(case_when(enr_underweight=="None" ~ 1,
                                              enr_underweight=="Moderate" ~ 2,
                                              enr_underweight=="Severe" ~ 3,
                                              is.na(enr_underweight) ~ 4), levels=1:4, labels=c("None","Moderate","Severe","Unknown")),
         underweight_ordered2=factor(ifelse(underweight_ordered!="Unknown",as.numeric(underweight_ordered),NA), levels=1:3, labels=c("None","Moderate","Severe")),
         wasting_ordered=factor(case_when(enr_wasting=="None" ~ 1,
                                          enr_wasting=="Moderate" ~ 2,
                                          enr_wasting=="Severe" ~ 3,
                                          is.na(enr_wasting) ~ 4), levels=1:4, labels=c("None","Moderate","Severe","Unknown")),
         wasting_ordered2=factor(ifelse(wasting_ordered!="Unknown",as.numeric(wasting_ordered),NA), levels=1:3, labels=c("None","Moderate","Severe")),
         stunting_ordered=factor(case_when(enr_stunting=="Not Stunted" ~ 1,
                                           enr_stunting=="Stunted" ~ 2,
                                           is.na(enr_stunting) ~ 3), levels=1:3, labels=c("Not stunted","Stunted","Unknown")),
         stunting_ordered2=factor(ifelse(stunting_ordered!="Unknown",as.numeric(stunting_ordered),NA), levels=1:2, labels=c("Not stunted","Stunted")),
         water_factor=factor(water_drink_level, levels=1:3, labels=c("Improved source","Unimproved source","Surface water")),
         water_factor2=factor(case_when(water_drink_level==1 ~ 1,
                                        water_drink_level %in% c(2,3) ~ 2), levels=1:2, labels=c("Improved source","Unimproved source")),
         sanitation_factor=factor(sanitation_level, levels=1:3, labels=c("Improved source","Unimproved source","Open defecation")),
         sanitation_factor2=factor(case_when(sanitation_level==1 ~ 1,
                                             sanitation_level %in% c(2,3) ~ 2), levels=1:2, labels=c("Improved source","Unimproved source")),
         deyhd_factor=factor(case_when(enroll_cond_dehyd==1 ~ 3,
                                       enroll_cond_dehyd==2 ~ 2,
                                       enroll_cond_dehyd==3 ~ 1), levels=1:3, labels=c("None","Some","Severe")),
         careseek_factor=factor(case_when(enroll_hs_beh_wh___1==0 & 
                                          enroll_hs_beh_wh___2==0 & 
                                          enroll_hs_beh_wh___3==0 & 
                                          enroll_hs_beh_wh___4==0 & 
                                          enroll_hs_beh_wh___5==0 & 
                                          enroll_hs_beh_wh___6==0 & 
                                          enroll_hs_beh_wh___7==0 & 
                                          enroll_hs_beh_wh___8==0 & 
                                          enroll_hs_beh_wh___9==0 ~ 0,
                                          enroll_hs_beh_wh___7==1 | 
                                          enroll_hs_beh_wh___8==1 ~ 1,
                                          enroll_hs_beh_wh___1==1 | 
                                          enroll_hs_beh_wh___2==1 | 
                                          enroll_hs_beh_wh___3==1 | 
                                          enroll_hs_beh_wh___4==1 | 
                                          enroll_hs_beh_wh___5==1 | 
                                          enroll_hs_beh_wh___6==1 | 
                                          enroll_hs_beh_wh___9==1 ~ 2), levels=0:2, labels=c("None", "Inpatient or outpatient health facility", "Other care-seeking")),
         healer=case_when(enroll_hs_beh_wh___1==1 | enroll_hs_beh_wh___2==1 ~ 1,
                          enroll_hs_beh_wh___1==0 & enroll_hs_beh_wh___2==0 ~ 0),
         oedema=case_when(enroll_oedema.x %in% c(2,3) ~ 0,
                          enroll_oedema.x==1          ~ 1),
         lrti=case_when(enroll_cough==1 | enroll_sev_resp==1 | enroll_dif_breath==1 | enroll_chest==1 | enroll_pneum==1 | enroll_cyan==1 | enroll_ox_sat==1 ~ 1,
                        enroll_cough %in% c(2,3) & enroll_sev_resp %in% c(2,3) & enroll_dif_breath %in% c(2,3) & enroll_chest %in% c(2,3) & enroll_pneum %in% c(2,3) & enroll_cyan %in% c(2,3) & enroll_ox_sat %in% c(2,3) ~ 0),
         lrti=ifelse(is.na(lrti),0,lrti),
         stiffneck=case_when(enroll_stiff_neck %in% c(2,3) ~ 0,
                             enroll_stiff_neck==1          ~ 1),
         rash=case_when(enroll_rash %in% c(2,3) | is.na(enroll_rash) ~ 0,
                        enroll_rash==1          ~ 1),
         convulsion=case_when(enroll_conv %in% c(2,3) ~ 0,
                              enroll_conv==1          ~ 1),
         lethargy=case_when(enroll_leth %in% c(2,3) ~ 0,
                            enroll_leth==1          ~ 1),
         pallor=case_when(enroll_palmor %in% c(2,3) ~ 0,
                          enroll_palmor==1          ~ 1),
         total=1,
         days_diarrhea=as.numeric(as.Date(scr_date, format="%Y-%m-%d")-as.Date(scr_diar_when, format="%Y-%m-%d")))

# Export last step dataset
write_rds(enroll_data, paste0("./Last Step Datasets/dcs_summary.Rds")) 

################################################
######### CREATE PEHUS SUMMARY DATASET #########
################################################
# Number of clusters completed
cluster <- clusters %>%
  mutate(cluster_num=floor(ClusterID),
         cluster=1) %>%
  select(clcom_country, ClusterID,cluster) %>%
  distinct() %>%
  group_by(clcom_country,cluster) %>%
  summarise(frequency = n()) %>%
  rename(pop_country=clcom_country)

cluster_total <- cluster %>%
  mutate(pop_country=8) %>%
  group_by(pop_country,cluster) %>%
  summarise(frequency = sum(frequency)) %>%
  bind_rows(cluster)

# area of clusters
cluster_area <- clusters %>%
  mutate(area=1) %>%
  group_by(clcom_country,area) %>%
  summarise(sum_total = sum(Area_m3)/1000000) %>%
  rename(pop_country=clcom_country)

cluster_area_total <- cluster_area %>%
  mutate(pop_country=8) %>%
  group_by(pop_country,area) %>%
  summarise(sum_total = sum(sum_total)) %>%
  bind_rows(cluster_area)

# % clusters completed 
cluster_complete <- cluster_comp %>%
  rename(ClusterID=clcom_cluster_id) %>%
  full_join(clusters, by = c("clcom_country", "ClusterID")) %>%
  mutate(cluster_summary=case_when(!is.na(clcom_households) & clcom_households==1          ~ 1,
                                   !is.na(clcom_households) & clcom_households==2          ~ 2,
                                   !is.na(clcom_households) & clcom_households %in% c(3,4) ~ 3,
                                   is.na(clcom_households)                                 ~ 4)) %>%
  group_by(clcom_country,cluster_summary) %>%
  summarise(frequency = n()) %>%
  mutate(percent = (frequency / sum(frequency))*100) %>%
  rename(pop_country=clcom_country)

cluster_complete_total <- cluster_comp %>%
  rename(ClusterID=clcom_cluster_id) %>%
  full_join(clusters, by = c("clcom_country", "ClusterID")) %>%
  mutate(cluster_summary=case_when(!is.na(clcom_households) & clcom_households==1          ~ 1,
                                   !is.na(clcom_households) & clcom_households==2          ~ 2,
                                   !is.na(clcom_households) & clcom_households %in% c(3,4) ~ 3,
                                   is.na(clcom_households)                                 ~ 4)) %>%
  group_by(cluster_summary) %>%
  summarise(frequency = n()) %>%
  mutate(percent = (frequency / sum(frequency))*100,
         pop_country=8) %>%
  bind_rows(cluster_complete)

# Total individuals enumerated
ind_enum <- pehus %>%
  filter(pop_hhold_consent==1) %>%
  select(pop_hh_id, pop_country, pop_num_residents) %>%
  distinct() %>%
  group_by(pop_country) %>%
  summarise(frequency=sum(pop_num_residents)) %>%
  mutate(ind=1)

ind_enum_total <- ind_enum %>%
  mutate(pop_country=8) %>%
  group_by(pop_country, ind) %>%
  summarise(frequency=sum(frequency)) %>%
  bind_rows(ind_enum)


children_enum_total <- pehus %>%
  filter(pop_u5==1) |>
  distinct() %>%
  group_by(pop_country) %>%
  summarise(frequency=sum(pop_num_residents)) %>%
  mutate(ind=1)

# Children 6-35 months of age
children_enum <- pehus %>%
  filter(pop_age>=6 & pop_age<36) %>%
  group_by(pop_country) %>%
  summarise(frequency=n()) %>%
  mutate(children=1)

children_enum_total <- children_enum %>%
  mutate(pop_country=8) %>%
  group_by(pop_country, children) %>%
  summarise(frequency=sum(frequency)) %>%
  bind_rows(children_enum)

# Children with diarrhea
children_diar <- pehus %>%
  filter(pop_age>=6 & pop_age<36) %>%
  group_by(pop_country,pop_diarrh_2weeks) %>%
  summarise(frequency=n()) %>%
  mutate(percent = (frequency / sum(frequency))*100)

children_diar_total <- pehus %>%
  filter(pop_age>=6 & pop_age<36) %>%
  group_by(pop_diarrh_2weeks) %>%
  summarise(frequency=n()) %>%
  mutate(percent = (frequency / sum(frequency))*100,
         pop_country=8) %>%
  bind_rows(children_diar)

# Consent to the HUS
hus_consent <- pehus %>%
  filter(pop_age>=6 & pop_age<36 & pop_diarrh_2weeks==1) %>%
  mutate(consent=case_when(re_consent %in% c(1,2) | hus_consent %in% c(1,2) ~ 1,
                           re_consent==0 | hus_consent==0 ~ 0,
                           hus_consent==3 & is.na(re_consent) ~ 3)) %>%
  group_by(pop_country,consent) %>%
  summarise(frequency=n()) %>%
  mutate(percent = (frequency / sum(frequency))*100)

hus_consent_total <- pehus %>%
  filter(pop_age>=6 & pop_age<36 & pop_diarrh_2weeks==1) %>%
  mutate(consent=case_when(re_consent %in% c(1,2) | hus_consent %in% c(1,2) ~ 1,
                           re_consent==0 | hus_consent==0 ~ 0,
                           hus_consent==3 & is.na(re_consent) ~ 3)) %>%
  group_by(consent) %>%
  summarise(frequency=n()) %>%
  mutate(percent = (frequency / sum(frequency))*100,
         pop_country=8) %>%
  bind_rows(hus_consent)

table(pehus$hus_consent, pehus$re_consent, useNA = "ifany")

# Sought care
careseek_data <- careseek %>%
  mutate(careseek_none=ifelse(hus_care==0,1,0),
         careseek_efgh=ifelse(hus_careseek_efgh==1,1,0),
         careseek_mad=ifelse(hus_careseek_mad==1 & hus_careseek_efgh==0,1,0),
         careseek_outpost=ifelse(hus_care_oth==6,1,0),
         careseek_pharm=ifelse(hus_care_oth %in% c(3,4),1,0),
         careseek_healer=ifelse(hus_care_oth %in% c(1,2),1,0),
         careseek_oth=ifelse(hus_care_oth %in% c(5,7),1,0),
         careseek_sum=case_when(hus_careseek_efgh==1 | hus_careseek_mad==1 ~ 1,
                                hus_care==0                                ~ 0,
                                hus_care==1                                ~ 2))

# Write functions generating the facility, site and overall number of shigella positives by TAC and culture stratified by dysentery
site_fun=function(a) {
  as.data.frame(careseek_data %>%
                  group_by(pop_country,.data[[a]]) %>%
                  summarise(frequency=n()) %>%
                  mutate(percent = (frequency / sum(frequency))*100) %>%
                  filter(!is.na(UQ(sym(a)))))
}
total_fun=function(a) {
  as.data.frame(careseek_data %>%
                  group_by(.data[[a]]) %>%
                  summarise(frequency=n()) %>%
                  mutate(percent = (frequency / sum(frequency))*100) %>%
                  filter(!is.na(UQ(sym(a)))))
}

# wrapper function
subgroups_fun = function(a, b) {
  if(b=="site") {
    site_fun(a)
  } else if(b=="total") {
    total_fun(a)
  } else {print("Incorrect level selected")}
}

# specifying objects needed to run through all level/subgroup combinations
fun_subgroup = c("careseek_none","careseek_efgh","careseek_sum","careseek_mad","careseek_outpost","careseek_pharm","careseek_healer","careseek_oth")
fun_level = c("site", "total")
grid = expand.grid(fun_subgroup, fun_level, stringsAsFactors = FALSE)

# run wrapper function through each combination
num_subgroup = setNames(do.call(mapply, c("subgroups_fun", unname(as.list(grid)))), paste0("careseek_", grid$Var1, "_", grid$Var2))

# Rbind list plus overall non-stratified
careseek_final = bind_rows(num_subgroup) %>%
  mutate(pop_country=ifelse(is.na(pop_country),8,pop_country))

# days between diarrhea and care seeking
careseeking_days <- careseek %>%
  mutate(start_dt=ifelse(hus_diarrh_start==1,as.Date(hus_diarrh_start_dt, format="%B %d, %Y"), 
                         ifelse(hus_diarrh_start==0,as.Date(pop_visit_dt, format="%m/%d/%Y")-hus_diarrh_start_days+1,NA)), 
         start_dt2=ifelse(hus_diarrh_start==1,as.Date(hus_diarrh_start_dt, format="%m/%d/%Y"),NA),
         start_dt_final=ifelse(is.na(start_dt),start_dt2,start_dt),
         efgh_date=coalesce(as.Date(hus_care_fac_ba_dt, format="%B %d, %Y"),
                            as.Date(hus_care_fac_ke_dt, format="%B %d, %Y"),
                            as.Date(hus_care_fac_mw_dt, format="%B %d, %Y"),
                            as.Date(hus_care_fac_ml_dt, format="%B %d, %Y"),
                            as.Date(hus_care_fac_pa_dt, format="%B %d, %Y"),
                            as.Date(hus_care_fac_pe_dt, format="%B %d, %Y"),
                            as.Date(hus_care_fac_ga_dt, format="%B %d, %Y")),
         outpatient_date=as.Date(hus_care_out_dt, format="%B %d, %Y"),
         inpatient_date=as.Date(hus_care_in_dt, format="%B %d, %Y"),
         diarrhea_days=as.numeric(ifelse(!is.na(efgh_date) & !is.na(start_dt_final) & 
                                           (efgh_date<=inpatient_date | is.na(inpatient_date)) & 
                                           (efgh_date<=outpatient_date | is.na(outpatient_date)),efgh_date-start_dt_final+1,
                                         ifelse(!is.na(inpatient_date) & !is.na(start_dt_final) & 
                                                  (inpatient_date<=outpatient_date | is.na(outpatient_date)),inpatient_date-start_dt_final+1,
                                                ifelse(!is.na(outpatient_date) & !is.na(start_dt_final),outpatient_date-start_dt_final+1,NA)))),
         diarrhea_days=ifelse(diarrhea_days<=0,NA,diarrhea_days)) # Excluded n=36 negative values

careseeking_days_site <- careseeking_days %>%
  filter(!is.na(diarrhea_days)) %>%
  group_by(pop_country) %>%
  summarise(n=n(),
            median = quantile(diarrhea_days, probs = 0.5), 
            quant25= quantile(diarrhea_days, probs = 0.25),
            quant75= quantile(diarrhea_days, probs = 0.75)) %>%
  mutate(careseek_days=1)

careseeking_days_total <- careseeking_days %>%
  filter(!is.na(diarrhea_days)) %>%
  summarise(n=n(),
            median = quantile(diarrhea_days, probs = 0.5), 
            quant25= quantile(diarrhea_days, probs = 0.25),
            quant75= quantile(diarrhea_days, probs = 0.75)) %>%
  mutate(pop_country=8,
         careseek_days=1) %>%
  bind_rows(careseeking_days_site)

# Female sex
female <- pehus %>%
  filter(hus_consent %in% c(1,2) | re_consent %in% c(1,2)) %>%
  group_by(pop_country,pop_sex) %>%
  summarise(frequency=n()) %>%
  mutate(percent = (frequency / sum(frequency))*100)

female_total <- pehus %>%
  filter(hus_consent %in% c(1,2) | re_consent %in% c(1,2)) %>%
  group_by(pop_sex) %>%
  summarise(frequency=n()) %>%
  mutate(percent = (frequency / sum(frequency))*100,
         pop_country=8) %>%
  bind_rows(female) %>%
  filter(pop_sex==1)

# Age category
age_site <- pehus %>%
  filter((hus_consent %in% c(1,2) | re_consent %in% c(1,2))) %>%
  mutate(c_age=case_when(pop_age>=6 & pop_age<9   ~ 1,
                         pop_age>=9 & pop_age<12  ~ 2,
                         pop_age>=12 & pop_age<18 ~ 3,
                         pop_age>=18 & pop_age<24 ~ 4,
                         pop_age>=24 & pop_age<36 ~ 5)) %>%
  group_by(pop_country,c_age) %>%
  summarise(frequency=n()) %>%
  mutate(percent = (frequency / sum(frequency))*100)

age_total <- pehus %>%
  filter((hus_consent %in% c(1,2) | re_consent %in% c(1,2))) %>%
  mutate(c_age=case_when(pop_age>=6 & pop_age<9   ~ 1,
                         pop_age>=9 & pop_age<12  ~ 2,
                         pop_age>=12 & pop_age<18 ~ 3,
                         pop_age>=18 & pop_age<24 ~ 4,
                         pop_age>=24 & pop_age<36 ~ 5)) %>%
  group_by(c_age) %>%
  summarise(frequency=n()) %>%
  mutate(percent = (frequency / sum(frequency))*100,
         pop_country=8) %>%
  bind_rows(age_site) 

# Continuous age
agecont_site <- pehus %>%
  filter((hus_consent %in% c(1,2) | re_consent %in% c(1,2))) %>%
  mutate(age_months=floor(pop_age)) %>%
  group_by(pop_country) %>%
  summarise(median = quantile(age_months, probs = 0.5), 
            quant25= quantile(age_months, probs = 0.25),
            quant75= quantile(age_months, probs = 0.75)) %>%
  mutate(agecont=1)

agecont_total <- pehus %>%
  filter((hus_consent %in% c(1,2) | re_consent %in% c(1,2))) %>%
  mutate(age_months=floor(pop_age)) %>%
  summarise(median = quantile(age_months, probs = 0.5), 
            quant25= quantile(age_months, probs = 0.25),
            quant75= quantile(age_months, probs = 0.75)) %>%
  mutate(pop_country=8,
         agecont=1) %>%
  bind_rows(agecont_site) 

# symptoms
# Write functions generating site and overall 
site_fun2=function(a) {
  as.data.frame(pehus %>%
                  filter((hus_consent %in% c(1,2) | re_consent %in% c(1,2))) %>%
                  group_by(pop_country,.data[[a]]) %>%
                  summarise(frequency=n()) %>%
                  mutate(percent = (frequency / sum(frequency))*100) %>%
                  filter(.data[[a]]==1))
}
total_fun2=function(a) {
  as.data.frame(pehus %>%
                  filter((hus_consent %in% c(1,2) | re_consent %in% c(1,2))) %>%
                  group_by(.data[[a]]) %>%
                  summarise(frequency=n()) %>%
                  mutate(percent = (frequency / sum(frequency))*100)%>%
                  filter(.data[[a]]==1))
}

# wrapper function
subgroups_fun2 = function(a, b) {
  if(b=="site") {
    site_fun2(a)
  } else if(b=="total") {
    total_fun2(a)
  } else {print("Incorrect level selected")}
}

# specifying objects needed to run through all level/subgroup combinations
fun_subgroup2 = c("hus_blood","hus_irritable","hus_thirsty","hus_drinks","hus_wrinkled_skin","hus_nodrink","hus_sunken_eyes","hus_lethargy")
fun_level2 = c("site", "total")
grid2 = expand.grid(fun_subgroup2, fun_level2, stringsAsFactors = FALSE)

# run wrapper function through each combination
num_subgroup2 = setNames(do.call(mapply, c("subgroups_fun2", unname(as.list(grid2)))), paste0("hus_", grid2$Var1, "_", grid2$Var2))

# Rbind list plus overall non-stratified
symptoms = bind_rows(num_subgroup2) %>%
  mutate(pop_country=ifelse(is.na(pop_country),8,pop_country))

# Combine datasets
pehus_final=bind_rows(cluster_total, cluster_area_total, cluster_complete_total,ind_enum_total, children_enum_total, children_diar_total,  
                      hus_consent_total, careseek_final,careseeking_days_total,female_total, age_total, agecont_total,symptoms) %>%
  mutate(country=case_when(pop_country==1 ~ "Bangladesh",
                           pop_country==2 ~ "Kenya",
                           pop_country==3 ~ "Malawi",
                           pop_country==4 ~ "Mali",
                           pop_country==5 ~ "Pakistan",
                           pop_country==6 ~ "Peru",
                           pop_country==7 ~ "The Gambia",
                           pop_country==8 ~ "Total"))

# Export last step dataset
write_rds(pehus_final, paste0("./Last Step Datasets/pehus_table1_.Rds")) 

#########################
#### Antibiotic data ####
#########################
# First create a dataset showing graph of important species
antibiotic_data <- enrollment %>%
  left_join(rectal_swab_results, by = "pid") %>%
  filter(!is.na(swi_isolate_group)) %>%
  #filter(!duplicated(pid, swi_isolate_group, swi_isolate_type)) %>%
  select(pid, ends_with("_zone"), ends_with("_int"),swi_ast_az_mic,swi_isolate_group,swi_isolate_type) |>
  group_by(pid, swi_isolate_group, swi_isolate_type) |>
  summarise(swi_ast_amp_zone=  min(swi_ast_amp_zone),
            swi_ast_az_zone=   min(swi_ast_az_zone),
            swi_ast_ceft_zone= min(swi_ast_ceft_zone),
            swi_ast_cipro_zone=min(swi_ast_cipro_zone),
            swi_ast_nal_zone=  min(swi_ast_nal_zone),
            swi_ast_piv_zone=  min(swi_ast_piv_zone),
            swi_ast_tri_zone=  min(swi_ast_tri_zone),
            swi_ast_amp_int=   max(swi_ast_amp_int),
            swi_ast_az_int=    max(swi_ast_az_int),
            swi_ast_ceft_int=  max(swi_ast_ceft_int),
            swi_ast_cipro_int= max(swi_ast_cipro_int),
            swi_ast_nal_int=   max(swi_ast_nal_int),
            swi_ast_piv_int=   max(swi_ast_piv_int),
            swi_ast_tri_int=   max(swi_ast_tri_int),
            swi_ast_az_mic=    max(swi_ast_az_mic))
  
antibiotic_data_zone <- antibiotic_data %>%
  pivot_longer(!c(pid,swi_isolate_group,swi_isolate_type), names_to = "antibiotic", values_to = "zone") %>%
  mutate(antibiotic=ifelse(antibiotic=="swi_ast_amp_zone","Ampicillin",
                           ifelse(antibiotic=="swi_ast_az_zone","Azithromycin",
                                  ifelse(antibiotic=="swi_ast_ceft_zone","Ceftriaxone",
                                         ifelse(antibiotic=="swi_ast_cipro_zone","Ciprofloxacin",
                                                ifelse(antibiotic=="swi_ast_nal_zone","Nalidixic acid",
                                                       ifelse(antibiotic=="swi_ast_piv_zone","Pivemicellinam",
                                                              ifelse(antibiotic=="swi_ast_tri_zone","Trimethoprim/sulfamethoxazole", NA)))))))) %>%
  filter(!is.na(antibiotic))

antibiotic_data_int <- antibiotic_data %>%
  pivot_longer(!c(pid,swi_isolate_group,swi_isolate_type), names_to = "antibiotic", values_to = "int") %>%
  mutate(antibiotic=ifelse(antibiotic=="swi_ast_amp_int","Ampicillin",
                           ifelse(antibiotic=="swi_ast_az_int","Azithromycin",
                                  ifelse(antibiotic=="swi_ast_ceft_int","Ceftriaxone",
                                         ifelse(antibiotic=="swi_ast_cipro_int","Ciprofloxacin",
                                                ifelse(antibiotic=="swi_ast_nal_int","Nalidixic acid",
                                                       ifelse(antibiotic=="swi_ast_piv_int","Pivemicellinam",
                                                              ifelse(antibiotic=="swi_ast_tri_int","Trimethoprim/sulfamethoxazole", NA)))))))) %>%
  filter(!is.na(antibiotic))

antibiotic_data_mic <- antibiotic_data %>%
  pivot_longer(!c(pid,swi_isolate_group,swi_isolate_type), names_to = "antibiotic", values_to = "mic") %>%
  mutate(antibiotic=ifelse(antibiotic=="swi_ast_az_mic","Azithromycin",NA)) %>%
  filter(!is.na(antibiotic))  

antibiotic_data2 <- antibiotic_data_zone %>%
  left_join(antibiotic_data_int, by = c("pid","swi_isolate_group","swi_isolate_type","antibiotic")) %>%
  left_join(antibiotic_data_mic, by = c("pid","swi_isolate_group","swi_isolate_type","antibiotic")) %>%
  mutate(int_final=ifelse((antibiotic=="Ampicillin" & zone>16 & !is.na(zone)) |
                          (antibiotic=="Azithromycin" & ((zone>15 & !is.na(zone)) | (int==999 & mic<16))) |
                          (antibiotic=="Ceftriaxone" & zone>22 & !is.na(zone)) |
                          (antibiotic=="Ciprofloxacin" & zone>25 & !is.na(zone)) |
                          (antibiotic=="Nalidixic acid" & zone>18 & !is.na(zone)) |
                          (antibiotic=="Pivemicellinam" & zone>14 & !is.na(zone)) |
                          (antibiotic=="Trimethoprim/sulfamethoxazole" & zone>15 & !is.na(zone)),1,
                   ifelse((antibiotic=="Ampicillin" & zone %in% 14:16) |
                          (antibiotic=="Azithromycin" & ((zone %in% 11:15) | (int==999 & mic==16))) |
                          (antibiotic=="Ceftriaxone" & zone %in% 20:22) |
                          (antibiotic=="Ciprofloxacin" & zone %in% 22:25) |
                          (antibiotic=="Nalidixic acid" & zone %in% 14:18) |
                          (antibiotic=="Pivemicellinam" & zone %in% 12:14) |
                          (antibiotic=="Trimethoprim/sulfamethoxazole" & zone %in% 11:15),2,
                   ifelse((antibiotic=="Ampicillin" & zone %in% 0:13) |
                          (antibiotic=="Azithromycin" & ((zone %in% 0:10) | (int==999 & mic >16))) |
                          (antibiotic=="Ceftriaxone" & zone %in% 0:19) |
                          (antibiotic=="Ciprofloxacin" & zone %in% 0:21) |
                          (antibiotic=="Nalidixic acid" & zone %in% 0:13) |
                          (antibiotic=="Pivemicellinam" & zone %in% 0:11) |
                          (antibiotic=="Trimethoprim/sulfamethoxazole" & zone %in% 0:10),3,int))),
         test=ifelse(int_final!=int,1,0),
         b_nonsuscept=ifelse(int_final %in% c(2,3),1,
                      ifelse(int_final==1,0,NA)),
         b_resistant=ifelse(int_final==3,1,
                     ifelse(int_final %in% c(1,2),0,NA)),,
         Resistance=factor(int_final,levels=c(1,2,3),
                           labels = c("Susceptible","Intermediate","Resistant")),
         group=factor(swi_isolate_group,levels=c(1:5), labels=c("S. dysenteriae","S. flexneri","S. boydii","S. sonnei","Undetermined")))

# All antibiotic results
antibiotic_all <- antibiotic_data2 %>%
  filter(!is.na(Resistance) & !is.na(antibiotic)) %>%
  left_join(enrollment %>% select(pid, enroll_site), by = "pid")

###### BUILD DATA FOR FIGURES #######
# Total antibiotic resistance
abx_n_all <- antibiotic_all %>%
  group_by(antibiotic,Resistance) %>%
  summarise(`Number of Samples`=n())

abx_tot_all <- antibiotic_all %>%
  group_by(antibiotic) %>%
  summarise(tot=n())

antibiotic_tot_all <- abx_n_all %>%
  left_join(abx_tot_all,by="antibiotic") %>%
  mutate(Percent=paste0(sprintf("%.1f", round(`Number of Samples`/tot,3)*100),"%"),
         Antibiotic=antibiotic) %>%
  mutate(type='All isolates')

# Total flexneri
abx_n_flex <- antibiotic_all %>%
  filter(group=='S. flexneri') %>%
  group_by(antibiotic,Resistance) %>%
  summarise(`Number of Samples`=n()) 
abx_tot_flex <- antibiotic_all %>%
  filter(group=='S. flexneri') %>%
  group_by(antibiotic) %>%
  summarise(tot=n())

antibiotic_flex <- abx_n_flex %>%
  left_join(abx_tot_flex,by="antibiotic") %>%
  mutate(Percent=paste0(sprintf("%.1f", round(`Number of Samples`/tot,3)*100),"%"),
         Antibiotic=antibiotic) %>%
  mutate(type='S. flexneri')

# Total sonnei
abx_n_sonnei <- antibiotic_all %>%
  filter(group=='S. sonnei') %>%
  group_by(antibiotic,Resistance) %>%
  summarise(`Number of Samples`=n()) 

abx_tot_sonnei <- antibiotic_all %>%
  filter(group=='S. sonnei') %>%
  group_by(antibiotic) %>%
  summarise(tot=n())

antibiotic_sonnei <- abx_n_sonnei %>%
  left_join(abx_tot_sonnei,by="antibiotic") %>%
  mutate(Percent=paste0(sprintf("%.1f", round(`Number of Samples`/tot,3)*100),"%"),
         Antibiotic=antibiotic) %>%
  mutate(type='S. sonnei')

# Total boydii
abx_n_boydii <- antibiotic_all %>%
  filter(group=='S. boydii') %>%
  group_by(antibiotic,Resistance) %>%
  summarise(`Number of Samples`=n()) 

abx_tot_boydii <- antibiotic_all %>%
  filter(group=='S. boydii') %>%
  group_by(antibiotic) %>%
  summarise(tot=n())

antibiotic_boydii <- abx_n_boydii %>%
  left_join(abx_tot_boydii,by="antibiotic") %>%
  mutate(Percent=paste0(sprintf("%.1f", round(`Number of Samples`/tot,3)*100),"%"),
         Antibiotic=antibiotic) %>%
  mutate(type='S. boydii')

# Total dysentariae
abx_n_dysenteriae <- antibiotic_all %>%
  filter(group=='S. dysenteriae') %>%
  group_by(antibiotic,Resistance) %>%
  summarise(`Number of Samples`=n()) 

abx_tot_dysenteriae <- antibiotic_all %>%
  filter(group=='S. dysenteriae') %>%
  group_by(antibiotic) %>%
  summarise(tot=n())

antibiotic_dysenteriae <- abx_n_dysenteriae %>%
  left_join(abx_tot_dysenteriae,by="antibiotic") %>%
  mutate(Percent=paste0(sprintf("%.1f", round(`Number of Samples`/tot,3)*100),"%"),
         Antibiotic=antibiotic) %>%
  mutate(type='S. dysenteriae')

# total in South America
abx_n_peru <- antibiotic_all %>%
  filter(enroll_site==6) %>%
  group_by(antibiotic,Resistance) %>%
  summarise(`Number of Samples`=n())

abx_tot_peru <- antibiotic_all %>%
  filter(enroll_site==6) %>%
  group_by(antibiotic) %>%
  summarise(tot=n())

antibiotic_tot_peru <- abx_n_peru %>%
  left_join(abx_tot_peru,by="antibiotic") %>%
  mutate(Percent=paste0(sprintf("%.1f", round(`Number of Samples`/tot,3)*100),"%"),
         Antibiotic=antibiotic) %>%
  mutate(type='All isolates',
         area="South America")

# total in East Africa
abx_n_east_africa <- antibiotic_all %>%
  filter(enroll_site %in% c(2,3)) %>%
  group_by(antibiotic,Resistance) %>%
  summarise(`Number of Samples`=n())

abx_tot_east_africa <- antibiotic_all %>%
  filter(enroll_site %in% c(2,3)) %>%
  group_by(antibiotic) %>%
  summarise(tot=n())

antibiotic_tot_east_africa <- abx_n_east_africa %>%
  left_join(abx_tot_east_africa,by="antibiotic") %>%
  mutate(Percent=paste0(sprintf("%.1f", round(`Number of Samples`/tot,3)*100),"%"),
         Antibiotic=antibiotic) %>%
  mutate(type='All isolates',
         area="East Africa")

# total in West Africa
abx_n_west_africa <- antibiotic_all %>%
  filter(enroll_site %in% c(4,7)) %>%
  group_by(antibiotic,Resistance) %>%
  summarise(`Number of Samples`=n())

abx_tot_west_africa <- antibiotic_all %>%
  filter(enroll_site %in% c(4,7)) %>%
  group_by(antibiotic) %>%
  summarise(tot=n())

antibiotic_tot_west_africa <- abx_n_west_africa %>%
  left_join(abx_tot_west_africa,by="antibiotic") %>%
  mutate(Percent=paste0(sprintf("%.1f", round(`Number of Samples`/tot,3)*100),"%"),
         Antibiotic=antibiotic) %>%
  mutate(type='All isolates',
         area="West Africa")

# total in Asia
abx_n_asia <- antibiotic_all %>%
  filter(enroll_site %in% c(1,5)) %>%
  group_by(antibiotic,Resistance) %>%
  summarise(`Number of Samples`=n())

abx_tot_asia <- antibiotic_all %>%
  filter(enroll_site %in% c(1,5)) %>%
  group_by(antibiotic) %>%
  summarise(tot=n())

antibiotic_tot_asia <- abx_n_asia %>%
  left_join(abx_tot_asia,by="antibiotic") %>%
  mutate(Percent=paste0(sprintf("%.1f", round(`Number of Samples`/tot,3)*100),"%"),
         Antibiotic=antibiotic) %>%
  mutate(type='All isolates',
         area="Asia")

# Total by site
# Total antibiotic resistance
abx_n_site <- antibiotic_all %>%
  group_by(antibiotic,enroll_site,Resistance) %>%
  summarise(`Number of Samples`=n())

abx_tot_site <- antibiotic_all %>%
  group_by(antibiotic,enroll_site) %>%
  summarise(tot=n())

antibiotic_tot_site <- abx_n_site %>%
  left_join(abx_tot_site,by=c("antibiotic","enroll_site")) %>%
  mutate(Percent=paste0(sprintf("%.1f", round(`Number of Samples`/tot,3)*100),"%"),
         Antibiotic=antibiotic) %>%
  mutate(type='Site - all isolates')

antibiotic_graphs <- antibiotic_tot_all %>%
  bind_rows(antibiotic_tot_site,
            antibiotic_flex,           antibiotic_sonnei,         antibiotic_boydii,  antibiotic_dysenteriae,
            antibiotic_tot_east_africa,antibiotic_tot_west_africa,antibiotic_tot_asia,antibiotic_tot_peru)

# Export last step dataset
write_rds(antibiotic_graphs, paste0("./Last Step Datasets/abx_fig.Rds")) 

### Create table
# Reshape results BACK into long form
antibiotic_wide <- antibiotic_all |>
  ungroup() |>
  select(enroll_site,pid,group,swi_isolate_type,antibiotic,b_nonsuscept) %>%
  mutate(antibiotic=ifelse(antibiotic=="Trimethoprim/sulfamethoxazole","Trimethoprim",antibiotic)) |>
  pivot_wider(id_cols = c("pid","group","swi_isolate_type"), names_from = "antibiotic", values_from = "b_nonsuscept") 

# Combine - one will be used for total, and the merge in antibiotic data and create indicators
antibiotic_wide2 <- enrollment %>%
  select(pid, enroll_site) |>
  bind_rows(enrollment |> select(pid, enroll_site) %>% mutate(enroll_site=9)) %>%
  left_join(rectal_swab_results |> select(pid,sw_isolate_num) |> distinct(), by="pid") |>
  left_join(antibiotic_wide, by = "pid") %>%
  mutate(available=ifelse(pid %in% rectal_swab_results$pid,1,0),
         shigella=ifelse(sw_isolate_num>0,1,0),
         who=ifelse(Azithromycin==1 & Ciprofloxacin==1 & Ceftriaxone==1,1,0),
         mdr=ifelse(Ampicillin+Azithromycin+Ceftriaxone+Ciprofloxacin+`Nalidixic acid`+Pivemicellinam+Trimethoprim>=3,1,0),
         xdr=ifelse(Ampicillin==1 & Azithromycin==1 & Ceftriaxone==1 & Ciprofloxacin==1 & Trimethoprim==1,1,0))

### Create table
# Reshape results BACK into long form
antibiotic_wide_temp <- antibiotic_all |>
  ungroup() |>
  select(enroll_site,pid,group,swi_isolate_type,antibiotic,b_resistant) %>%
  mutate(antibiotic=ifelse(antibiotic=="Trimethoprim/sulfamethoxazole","Trimethoprim",antibiotic)) |>
  pivot_wider(id_cols = c("pid","group","swi_isolate_type"), names_from = "antibiotic", values_from = "b_resistant") 

# Combine - one will be used for total, and the merge in antibiotic data and create indicators
antibiotic_wide_temp2 <- enrollment %>%
  select(pid, enroll_site) |>
  bind_rows(enrollment |> select(pid, enroll_site) %>% mutate(enroll_site=9)) %>%
  left_join(rectal_swab_results |> select(pid,sw_isolate_num) |> distinct(), by="pid") |>
  left_join(antibiotic_wide_temp, by = "pid") %>%
  mutate(available=ifelse(pid %in% rectal_swab_results$pid,1,0),
         shigella=ifelse(sw_isolate_num>0,1,0),
         who=ifelse(Azithromycin==1 & Ciprofloxacin==1 & Ceftriaxone==1,1,0),
         mdr=ifelse(Ampicillin+Azithromycin+Ceftriaxone+Ciprofloxacin+`Nalidixic acid`+Pivemicellinam+Trimethoprim>=3,1,0),
         xdr=ifelse(Ampicillin==1 & Azithromycin==1 & Ceftriaxone==1 & Ciprofloxacin==1 & Trimethoprim==1,1,0)) |>
  group_by(enroll_site,mdr) |>
  summarise(n=n()) |>
  mutate(N=sum(n))

# Export last step dataset
write_rds(antibiotic_wide2, paste0("./Last Step Datasets/abx_tab.Rds")) 

### ADDITIONAL DATA ANALYSIS FOR FIGURE ####
# Combine culture data with individual data to run indicators
shig_final <- antibiotic_data2 |>
  left_join(enrollment |> select(pid, enroll_site), by = "pid") |>
  ungroup() |>
  select(pid,enroll_site,group,swi_isolate_type,antibiotic,b_nonsuscept) |>
  pivot_wider(names_from = "antibiotic", values_from =  "b_nonsuscept") |>
  mutate(who=ifelse(Azithromycin==1 & Ciprofloxacin==1 & Ceftriaxone==1,1,0),
         who_any=ifelse(Azithromycin==1 | Ciprofloxacin==1 | Ceftriaxone==1,1,0),
         mdr=ifelse(Ampicillin+Azithromycin+Ceftriaxone+Ciprofloxacin+`Nalidixic acid`+Pivemicellinam+`Trimethoprim/sulfamethoxazole`>=3,1,0),
         xdr=ifelse(Ampicillin==1 & Azithromycin==1 & Ceftriaxone==1 & Ciprofloxacin==1 & `Trimethoprim/sulfamethoxazole`==1,1,0)) |>
  select(pid,enroll_site,group,swi_isolate_type,who,who_any,mdr,xdr) |>
  pivot_longer(!c(pid,enroll_site,group,swi_isolate_type), names_to = "antibiotic", values_to = "b_nonsuscept") |>
  bind_rows(antibiotic_data2 |> left_join(enrollment |> ungroup() |> select(pid, enroll_site), by = "pid") |>
                                select(pid,enroll_site,group,swi_isolate_type,antibiotic,b_nonsuscept)) |>
  # Make factors
  mutate(abx=factor(case_when(antibiotic=="Ampicillin" ~ 1,
                              antibiotic=="Azithromycin" ~ 2,
                              antibiotic=="Ciprofloxacin" ~ 3,
                              antibiotic=="Ceftriaxone" ~ 4,
                              antibiotic=="Nalidixic acid" ~ 5,
                              antibiotic=="Pivemicellinam" ~ 6,
                              antibiotic=="Trimethoprim/sulfamethoxazole" ~ 7,
                              antibiotic=="mdr" ~ 8,
                              antibiotic=="xdr" ~ 9,
                              antibiotic=="who_any" ~ 10,
                              antibiotic=="who" ~ 11), levels=1:11, labels=c("AMP","AZI","CIP","CEF","NA","PIV","TRI","MDR","XDR","Any WHO","All WHO")),
         country=factor(case_when(enroll_site==2 ~ 1,
                                  enroll_site==3 ~ 2,
                                  enroll_site==6 ~ 3,
                                  enroll_site==1 ~ 4,
                                  enroll_site==5 ~ 5,
                                  enroll_site==4 ~ 6,
                                  enroll_site==7 ~ 7), levels=1:7, labels=c("Kenya","Malawi","Peru","Bangladesh","Pakistan","Mali","The Gambia")),
         resistance=factor(b_nonsuscept, levels=0:1, labels=c("No","Yes")),
         region=factor(case_when(enroll_site %in% c(2,3) ~ 1,
                                 enroll_site==6 ~ 2,
                                 enroll_site %in% c(1,5) ~ 3,
                                 enroll_site %in% c(4,7) ~ 4), levels=1:4, labels=c("East Africa","South America","South Asia","West Africa")))

# Add in data again for total
shig_final2 <- bind_rows(shig_final,
                         shig_final |> mutate(country=factor(8,levels=8,labels="Total"),
                                              region=NA))
  
# All antibiotic results - by country and serotype
abx_n_country_sero <- shig_final2  %>%
  group_by(abx,country,group,resistance) %>%
  summarise(n=n())

abx_tot_country_sero <- shig_final2  %>%
  group_by(abx,country,group) %>%
  summarise(tot=n())

antibiotic_tot_country_sero <- abx_n_country_sero  %>%
  left_join(abx_tot_country_sero,by=c("abx","country","group")) %>%
  mutate(Percent=(n/tot)*100) 

# All antibiotic results - by country 
abx_n_country <- shig_final2 %>%
  group_by(abx,country,resistance) %>%
  summarise(n=n())

abx_tot_country <- shig_final2 %>%
  group_by(abx,country) %>%
  summarise(tot=n())

antibiotic_tot_country <- abx_n_country %>%
  left_join(abx_tot_country,by=c("abx","country")) %>%
  mutate(Percent=(n/tot)*100,
         group=factor(0,levels=0,labels="All Isolates"))

# Combine and merge in totals AGAIN to make 0 rows
abx_country_sero_merge <- bind_rows(antibiotic_tot_country,antibiotic_tot_country_sero) |> 
  filter(resistance=="Yes") |>
  filter(group %notin% c("Undetermined"))

abx_missing <- with(abx_country_sero_merge, expand.grid(group = levels(group), country = levels(country), abx = levels(abx)))

abx_country_sero_final <- abx_missing |>
  left_join(abx_country_sero_merge, by=c("country","abx","group")) |>
  filter(group %notin% c("Undetermined")) |>
  # fill in missing groups
  mutate(Percent=ifelse(is.na(Percent),0,Percent))

write_rds(abx_country_sero_final, paste0("./Last Step Datasets/abx.Rds")) 

############################
#### CO-OTHER PATHOGENS ####
############################
# Make any enteric pathogen variable (except Shigella)
any_enteric <- tac_long |>
  filter(target!="Shigella EIEC" & !is.na(b_attributable)) |>
  group_by(groupid,pid) |>
  summarize(b_attributable=max(b_attributable)) |> mutate(target="Any enteric pathogen")

# Pull out species of interest
oth_path <- tac_long %>%
  # Only keep attributable pathogens
  filter(!is.na(b_attributable)) |>
  # Add any enteric pathogen among shigella
  bind_rows(any_enteric) |>
  # merge in data needed for subgroups
  left_join(severity |> select(mvs_or_dysentery,gems_msd,pid), by = "pid") |> #MVS severity
  left_join(rectal_swab_results |> filter(sw_isolate_num>0) |>       # Culture positivity for shigella
                                   select(pid,sw_isolate_num) |>
                                   distinct(), by = "pid") |>
  left_join(tac_long |> filter(target=="Shigella EIEC" & b_attributable==1) |> 
                        mutate(tac_pos=1) |>
                        select(pid,tac_pos), by = "pid") |>
  left_join(anthro |> select(pid,enr_wasting), by="pid") |>
   mutate(# for culture comparison
          culture_pos=ifelse(!is.na(sw_isolate_num),1,0),
          tac_pos=ifelse(is.na(tac_pos),0,tac_pos),
          b_wasting=factor(case_when(enr_wasting %in% c("Moderate","Severe") ~ 1,
                                     enr_wasting=="None" ~ 2), levels=1:2,labels=c("Moderate/severe wasting","No wasting")))

####### FIRST PROPORTIONS!
# Run site-specific and overall proportions
oth_path_site <- oth_path |>
  filter(target!="Any enteric pathogen") |>
  group_by(groupid,target,b_attributable) |>
  summarise(n=n()) |>
  mutate(percent = (n / sum(n))*100,
         cat=1) |>
  filter(b_attributable==1)

oth_path_total <- oth_path |>
  filter(target!="Any enteric pathogen") |>
   group_by(target,b_attributable)|>
   summarise(n=n()) |>
   mutate(percent = (n / sum(n))*100,
          cat=1) |>
  filter(b_attributable==1)

# Run site-specific and overall proportions by MVS mod/heavy and or dysentery
oth_path_site_mvsdys <- oth_path |>
                  filter(mvs_or_dysentery=="Moderate/Severe or Dysentery" & target!="Any enteric pathogen") |>
                  group_by(groupid,target,b_attributable) |>
                  summarise(n=n()) %>%
                  mutate(percent = (n / sum(n))*100,
                         cat=2)|>
  filter(b_attributable==1)

oth_path_total_mvsdys <- oth_path |>
                  filter(mvs_or_dysentery=="Moderate/Severe or Dysentery" & target!="Any enteric pathogen") |>
                  group_by(target,b_attributable) %>% 
                  summarise(n=n()) %>%
                  mutate(percent = (n / sum(n))*100,
                         cat=2) |>
  filter(b_attributable==1)

# Run site-specific and overall proportions among shigella infections by TAC
oth_path_site_shig <- oth_path |>
                  filter(tac_pos==1 & target!="Shigella EIEC") %>%
                  group_by(groupid,target,b_attributable) %>%
                  summarise(n=n()) %>%
                  mutate(percent = (n / sum(n))*100,
                         cat=3) |>
  filter(b_attributable==1)

oth_path_total_shig <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC") %>%
  group_by(target,b_attributable) %>%
  summarise(n=n()) %>%
  mutate(percent = (n / sum(n))*100,
         cat=3) |>
  filter(b_attributable==1)

# Run site-specific and overall proportions among MVS mod/heavy +/- dysentery shigella-attribuntable infections
oth_path_site_shig_mvs <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC" & mvs_or_dysentery=="Moderate/Severe or Dysentery") %>%
  group_by(groupid,target,b_attributable) %>%
  summarise(n=n()) %>%
  mutate(percent = (n / sum(n))*100,
         cat=4) |>
  filter(b_attributable==1)

oth_path_total_shig_mvs <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC" & mvs_or_dysentery=="Moderate/Severe or Dysentery") %>%
  group_by(target,b_attributable) %>%
  summarise(n=n()) %>%
  mutate(percent = (n / sum(n))*100,
         cat=4) |>
  filter(b_attributable==1)

# Run site-specific and overall proportions among GEMS MSD shigella-attributable infections
oth_path_site_shig_gems <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC" & gems_msd=="Moderate-to-severe") %>%
  group_by(groupid,target,b_attributable) %>%
  summarise(n=n()) %>%
  mutate(percent = (n / sum(n))*100,
         cat=5) |>
  filter(b_attributable==1)

oth_path_total_shig_gems <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC" & gems_msd=="Moderate-to-severe") %>%
  group_by(target,b_attributable) %>%
  summarise(n=n()) %>%
  mutate(percent = (n / sum(n))*100,
         cat=5) |>
  filter(b_attributable==1)

# Run site-specific and overall proportions among shigella infections by culture
oth_path_site_shig2_mvs <- oth_path |>
  filter(culture_pos==1 & target!="Shigella EIEC") %>%
  group_by(groupid,target,b_attributable) %>%
  summarise(n=n()) %>%
  mutate(percent = (n / sum(n))*100,
         cat=6) |>
  filter(b_attributable==1)

oth_path_total_shig2_mvs <- oth_path |>
  filter(culture_pos==1 & target!="Shigella EIEC") %>%
  group_by(target,b_attributable) %>%
  summarise(n=n()) %>%
  mutate(percent = (n / sum(n))*100,
         cat=6) |>
  filter(b_attributable==1)

# Run site-specific and overall proportions among MVS moderate/severe or dysentery shigella infections by culture
oth_path_site_shig2 <- oth_path |>
  filter(culture_pos==1 & target!="Shigella EIEC" & mvs_or_dysentery=="Moderate/Severe or Dysentery") %>%
  group_by(groupid,target,b_attributable) %>%
  summarise(n=n()) %>%
  mutate(percent = (n / sum(n))*100,
         cat=7) |>
  filter(b_attributable==1)

oth_path_total_shig2 <- oth_path |>
  filter(culture_pos==1 & target!="Shigella EIEC" & mvs_or_dysentery=="Moderate/Severe or Dysentery") %>%
  group_by(target,b_attributable) %>%
  summarise(n=n()) %>%
  mutate(percent = (n / sum(n))*100,
         cat=7) |>
  filter(b_attributable==1)

# Run site-specific and overall proportions among GEMS MSD shigella infections by culture
oth_path_site_shig2_gems <- oth_path |>
  filter(culture_pos==1 & target!="Shigella EIEC" & gems_msd=="Moderate-to-severe") %>%
  group_by(groupid,target,b_attributable) %>%
  summarise(n=n()) %>%
  mutate(percent = (n / sum(n))*100,
         cat=8) |>
  filter(b_attributable==1)

oth_path_total_shig2_gems <- oth_path |>
  filter(culture_pos==1 & target!="Shigella EIEC" & gems_msd=="Moderate-to-severe") %>%
  group_by(target,b_attributable) %>%
  summarise(n=n()) %>%
  mutate(percent = (n / sum(n))*100,
         cat=8) |>
  filter(b_attributable==1)

# Run site-specific and overall proportions among wasted shigella-attributable infections
oth_path_site_shig_wasting <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC" & !is.na(b_wasting)) %>%
  group_by(groupid,target,b_wasting,b_attributable) %>%
  summarise(n=n()) %>%
  mutate(percent = (n / sum(n))*100,
         cat=9) |>
  filter(b_attributable==1)

oth_path_total_shig_wasting <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC" & !is.na(b_wasting)) %>%
  group_by(target,b_wasting,b_attributable) %>%
  summarise(n=n()) %>%
  mutate(percent = (n / sum(n))*100,
         cat=9) |>
  filter(b_attributable==1)

####### NEXT TOTALS
# Run site-specific and overall
oth_path_site_N <- oth_path |>
  filter(target!="Any enteric pathogen") |>
  group_by(groupid,target) |>
  summarise(N=n()) |>
  mutate(cat=1)

oth_path_total_N <- oth_path |>
  filter(target!="Any enteric pathogen") |>
  group_by(target)|>
  summarise(N=n()) |>
  mutate(cat=1) 

# Run site-specific and overall MVS mod/heavy and or dysentery
oth_path_site_mvsdys_N <- oth_path |>
  filter(mvs_or_dysentery=="Moderate/Severe or Dysentery" & target!="Any enteric pathogen") |>
  group_by(groupid,target) |>
  summarise(N=n()) %>%
  mutate(cat=2)

oth_path_total_mvsdys_N <- oth_path |>
  filter(mvs_or_dysentery=="Moderate/Severe or Dysentery" & target!="Any enteric pathogen") |>
  group_by(target) %>% 
  summarise(N=n()) %>%
  mutate(cat=2) 

# Run site-specific and overall proportions among shigella infections by TAC
oth_path_site_shig_N <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC") %>%
  group_by(groupid,target) %>%
  summarise(N=n()) %>%
  mutate(cat=3) 

oth_path_total_shig_N <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC") %>%
  group_by(target) %>%
  summarise(N=n()) %>%
  mutate(cat=3) 

# Run site-specific and overall proportions among MVS mod/severe or dysentery shigella infections by TAC
oth_path_site_shig_mvs_N <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC" & mvs_or_dysentery=="Moderate/Severe or Dysentery") %>%
  group_by(groupid,target) %>%
  summarise(N=n()) %>%
  mutate(cat=4) 

oth_path_total_shig_mvs_N <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC" & mvs_or_dysentery=="Moderate/Severe or Dysentery") %>%
  group_by(target) %>%
  summarise(N=n()) %>%
  mutate(cat=4) 

# Run site-specific and overall proportions among GEMS MSD shigella infections by TAC
oth_path_site_shig_gems_N <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC" & gems_msd=="Moderate-to-severe") %>%
  group_by(groupid,target) %>%
  summarise(N=n()) %>%
  mutate(cat=5) 

oth_path_total_shig_gems_N <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC" & gems_msd=="Moderate-to-severe") %>%
  group_by(target) %>%
  summarise(N=n()) %>%
  mutate(cat=5) 

# Run site-specific and overall proportions among shigella infections by culture
oth_path_site_shig2_N <- oth_path |>
  filter(culture_pos==1 & target!="Shigella EIEC") %>%
  group_by(groupid,target) %>%
  summarise(N=n()) %>%
  mutate(cat=6) 

oth_path_total_shig2_N <- oth_path |>
  filter(culture_pos==1 & target!="Shigella EIEC") %>%
  group_by(target) %>%
  summarise(N=n()) %>%
  mutate(cat=6) 

# Run site-specific and overall proportions among MVS mod/severe or dysentery shigella infections by culture
oth_path_site_shig2_mvs_N <- oth_path |>
  filter(culture_pos==1 & target!="Shigella EIEC" & mvs_or_dysentery=="Moderate/Severe or Dysentery") %>%
  group_by(groupid,target) %>%
  summarise(N=n()) %>%
  mutate(cat=7) 

oth_path_total_shig2_mvs_N <- oth_path |>
  filter(culture_pos==1 & target!="Shigella EIEC" & mvs_or_dysentery=="Moderate/Severe or Dysentery") %>%
  group_by(target) %>%
  summarise(N=n()) %>%
  mutate(cat=7) 

# Run site-specific and overall proportions among GEMS MSD shigella infections by culture
oth_path_site_shig2_gems_N <- oth_path |>
  filter(culture_pos==1 & target!="Shigella EIEC" & gems_msd=="Moderate-to-severe") %>%
  group_by(groupid,target) %>%
  summarise(N=n()) %>%
  mutate(cat=8) 

oth_path_total_shig2_gems_N <- oth_path |>
  filter(culture_pos==1 & target!="Shigella EIEC" & gems_msd=="Moderate-to-severe") %>%
  group_by(target) %>%
  summarise(N=n()) %>%
  mutate(cat=8) 

# Run site-specific and overall proportions among GEMS MSD shigella infections by TAC
oth_path_site_shig_wasting_N <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC" & !is.na(b_wasting)) %>%
  group_by(groupid,target,b_wasting) %>%
  summarise(N=n()) %>%
  mutate(cat=9) 

oth_path_total_shig_wasting_N <- oth_path |>
  filter(tac_pos==1 & target!="Shigella EIEC" & !is.na(b_wasting)) %>%
  group_by(target,b_wasting) %>%
  summarise(N=n()) %>%
  mutate(cat=9) 

# Rbind estimates and add missing groups
oth_path_N <- bind_rows(oth_path_site_N,
                        oth_path_total_N,
                        oth_path_site_mvsdys_N,
                        oth_path_total_mvsdys_N,
                        oth_path_site_shig_N,
                        oth_path_total_shig_N,
                        oth_path_site_shig_mvs_N,
                        oth_path_total_shig_mvs_N,
                        oth_path_site_shig_gems_N,
                        oth_path_total_shig_gems_N,
                        oth_path_site_shig2_N,
                        oth_path_total_shig2_N,
                        oth_path_site_shig2_mvs_N,
                        oth_path_total_shig2_mvs_N,
                        oth_path_site_shig2_gems_N,
                        oth_path_total_shig2_gems_N,
                        oth_path_site_shig_wasting_N,
                        oth_path_total_shig_wasting_N,) |>
  mutate(groupid=ifelse(is.na(groupid),"Total",groupid)) 

oth_path_final <- bind_rows(oth_path_site,
                            oth_path_total,
                            oth_path_site_mvsdys,
                            oth_path_total_mvsdys,
                            oth_path_site_shig,
                            oth_path_total_shig,
                            oth_path_site_shig_mvs,
                            oth_path_total_shig_mvs,
                            oth_path_site_shig_gems,
                            oth_path_total_shig_gems,
                            oth_path_site_shig2,
                            oth_path_total_shig2,
                            oth_path_site_shig2_mvs,
                            oth_path_total_shig2_mvs,
                            oth_path_site_shig2_gems,
                            oth_path_total_shig2_gems,
                            oth_path_site_shig_wasting,
                            oth_path_total_shig_wasting) |>
  mutate(groupid=ifelse(is.na(groupid),"Total",groupid)) |>
  full_join(oth_path_N, by=c("groupid","target","cat","b_wasting")) |>
  # fill in missing categories (basically 0s)
  mutate(n=ifelse(is.na(n),0,n),
         percent=ifelse(is.na(percent),0,percent)) |>
  select(-b_attributable)

# Export last step dataset
write_rds(oth_path_final, paste0("./Last Step Datasets/other_pathogens.Rds")) 

###########################################
######### CAUSE OF DEATH ANALYSIS #########
###########################################
cause_of_death <- cod |>
  left_join(enrollment |> select(pid,enroll_site,enroll_date), by="pid") |> 
  left_join(anthro |> select(pid,dob), by="pid") |>
  left_join(mortality |> select(pid,mort_cert_dt,mort_no_cert_dt), by="pid") |>
  left_join(shig |> select(pid,positive_tac_or_culture), by = "pid") |>
  mutate(country=case_when(enroll_site==1 ~ "Bangladesh",
                           enroll_site==2 ~ "Kenya",
                           enroll_site==3 ~ "Malawi",
                           enroll_site==4 ~ "Mali",
                           enroll_site==5 ~ "Pakistan",
                           enroll_site==6 ~ "Peru",
                           enroll_site==7 ~ "The Gambia",
                           enroll_site==9 ~ "Overall",
                           TRUE ~ NA),
         date_of_death=ifelse(is.na(mort_cert_dt),mort_no_cert_dt,mort_cert_dt),
         enroll_to_death=round(difftime(date_of_death,enroll_date,units="days")),
         age_at_death=interval(ymd(dob), ymd(date_of_death)) %/% months(1),
         shigella=ifelse(positive_tac_or_culture==0,"No Shigella",
                         ifelse(positive_tac_or_culture==1,"Shigella",NA)),
         cause_1a=ifelse(cause_1a_code=="","",paste0(cause_1a_code,": ",cause_1a_desc)),
         cause_1b=ifelse(cause_1b_code=="","",paste0(cause_1b_code,": ",cause_1b_desc)),
         cause_1c=ifelse(cause_1c_code=="","",paste0(cause_1c_code,": ",cause_1c_desc)),
         cause_2=ifelse(cause_2_code=="","",paste0(cause_2_code,": ",cause_2_desc)),
         cause_2_oth=ifelse(cause_2_oth1_code=="","",paste0(cause_2_oth1_code,": ",cause_2_oth1_desc)),
         cause_2_final=ifelse(cause_2_oth=="",cause_2,paste0(cause_2,", ",cause_2_oth))) |>
  select(country:cause_1c,cause_2_final)

write_rds(cause_of_death, paste0("./Last Step Datasets/cause_of_death.Rds")) 
