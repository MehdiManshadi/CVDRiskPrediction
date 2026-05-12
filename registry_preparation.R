# ==============================================================================
# Project : CVD Risk Prediction
# Script  : Drug Exposure Matrix Construction
# Author  : Mehdi Dehghan Manshadi
# Date    : 2026-05-11
# Description: Builds pre-diagnosis drug exposure matrices (binary, frequency,
#              first-use, last-use) for sporadic breast cancer patients in the
#              Uppsala/Örebro region using the Swedish drug registry.
# ==============================================================================

library(haven)
library(data.table)
# Load necessary libraries
library(dplyr)    # For data manipulation
library(ggplot2)  # For visualizations
library(openxlsx)
library(lubridate)
library(tidyr)
library(forcats)
library(survival)
library(survminer)
library(lubridate)
library(cmprsk)
library(mstate)
# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
PROJECT_DIR   <- "/Users/mehman/Projects/CVDRiskPrediction/"
REF_DATA_DIR  <- "../0_Reference_Data/0_otherRegions/"
OUTPUT_DIR    <- "../0_Reference_Data/0_otherRegions/Processed"

setwd(PROJECT_DIR)

# ------------------------------------------------------------------------------
# Load and filter sporadic breast cancer cohort
# ------------------------------------------------------------------------------
load(file.path(REF_DATA_DIR, "dt_fall.RData"))
Sporadic <- dt_fall

Sporadic$Date_of_BC <- as.Date(
  Sporadic$a_diag_dat,
  format = "%m/%d/%Y %I:%M:%S %p"
)

Sporadic <- Sporadic[
  Sporadic$Date_of_BC >= as.Date("2008-01-01") &
    Sporadic$Date_of_BC <  as.Date("2020-01-01"),
]

Sporadic <- Sporadic[Sporadic$REGION_NAMN == "Region Uppsala/Örebro", ]


Settings <- list(
  exclude_na_stage = FALSE,
  exclude_stage_50 = TRUE,
  exclude_tb_na = TRUE,
  exclude_n2_n3 = FALSE,
  exclude_t4_spo = TRUE,
  exclude_tis = FALSE,
  exclude_m1 = TRUE,
  exclude_chemo_na = TRUE,
  exclude_relapse_met = TRUE,
  exclude_rt_na = TRUE,
  exclude_chemo_missing = TRUE,
  min_followup_years = 0,
  max_end_age = 120,
  followup_start = as.Date("2008-01-01"),
  followup_end = as.Date("2020-01-01"),
  target_event = "first_MACE_after"
)
md <- Process_Raw_Data(Settings, BRCA, Sporadic)


# #########################

Process_Raw_Data <- function(Settings, Sporadic){
  # Add a new column to each data set to indicate the group (BRCA = 1, Sporadic = 0)
  load("../0_Reference_Data/0_otherRegions/Processed/AnyCVDCauseOfDeath.RData")
  # Add a Patient ID column to BRCA
  # lopnr
  
  # ---===---===---===---===---===---===---===---===---===---===---===---===---===
  # ---===---===---===---===---===---===---===---===---===---===---===---===---===
  # ---===---===---===---===---===---===---===---===---===---===---===---===---===
  # Settings and parameters: #####################################################  ******** Exclude all NA Stage, T and N patients ********
  
  MaxEndAge <- Settings$max_end_age
  FUstartDate <- Settings$followup_start
  FUendDate <- Settings$followup_end
  
  # BC diagnosis date, 2008-2019 (end of): #######################################  ******** end is 2020 not 2021 (correction is needed) ********    
  Sporadic$EXC <- 0
  Sporadic$Date_of_BC <- as.Date(Sporadic$a_diag_dat, format = "%m/%d/%Y %I:%M:%S %p")
  # Marking Excludes
  Sporadic$EXC <- ifelse(Sporadic$Date_of_BC < as.Date("2008-01-01") | 
                           Sporadic$Date_of_BC >= as.Date("2020-01-01"), 100 + Sporadic$EXC, Sporadic$EXC)
  # ---===---===---===---===---===---===---===---===---===---===---===---===---===
  d <- Sporadic #                ************ EXLCUDING out of range DATEs **************    
  d$Age <- d$a_pat_alder
  # ---===---===---===---===---===---===---===---===---===---===---===---===---===
  # ---===---===---===---===---===---===---===---===---===---===---===---===---===
  # Refining data, reading from the csv file #####################################
  # Side: ####
  # Excluding 9 and null from "side" for BRCA
  d$Side <- d$a_pat_sida_Varde
  d$EXC <- ifelse(d$Side == "Null" | is.na(d$Side) | d$Side == 9, 1 + d$EXC, d$EXC) #             ****** Excluding Null Sides, No case in BRCA******
  
  # Chemotherapy :####
  # ***first try. for the chemotherapy, lots of missing values
  # Excluding Nulls from "post_kemo_Varde" for BRCA
  # BRCA$PreChemo <- factor(ifelse(BRCA$post_kemo_Varde == "Null", NA, BRCA$post_kemo_Varde))
  # BRCA$EXC <- ifelse(BRCA$post_kemo_Varde == "Null", 10+BRCA$EXC, BRCA$EXC)
  # # Excluding Nulls from post_kemo_Varde for Sporadic
  # Sporadic$PreChemo <- factor(ifelse(Sporadic$post_kemo_Varde == "Null", NA, Sporadic$post_kemo_Varde))
  # Sporadic$EXC <- ifelse(Sporadic$post_kemo_Varde == "Null", 10+Sporadic$EXC, Sporadic$EXC)
  # ---===---===---===---===---===---===---===---===---===---===---===---===---===
  
  # Trying to fill the missing values using drug registry:                        
  BRCA_chem_data <- readRDS("/Users/mehman/1-KI-Projects/1_Matching/3-Extracting data/BRCA_chem_data.rds")
  Sporadic_chem_data <- readRDS("/Users/mehman/1-KI-Projects/1_Matching/3-Extracting data/Sporadic_chem_data.rds")
  chem_data <- rbind(Sporadic_chem_data, BRCA_chem_data)
  
  d$chemotherapy_adj <- "NA"
  d$chemotherapy_neoadj <- "NA" #                                                   *** Some of the FALSE Chemo values are based on ^no drug record^ ***
  d$chemotherapy_status <- "NA"
  
  for(indv in d$lopnr){
    ID <- d$lopnr == indv
    if(indv%in%chem_data$lopnr){
      if(!is.na(chem_data$op_diagnos[chem_data$lopnr == indv]) && chem_data$op_diagnos[chem_data$lopnr == indv] < 0){
        d$EXC[ID] <- 10 + d$EXC[ID]
      }else if(Settings$exclude_chemo_na&is.na(chem_data$chemotherapy_status[chem_data$lopnr == indv])){
        d$EXC[ID] <- 10 + d$EXC[ID]
      }else{
        d$chemotherapy_adj[ID] <- chem_data$chemotherapy_adj[chem_data$lopnr == indv]
        d$chemotherapy_neoadj[ID] <- chem_data$chemotherapy_neoadj[chem_data$lopnr == indv]
        d$chemotherapy_status[ID] <- chem_data$chemotherapy_status[chem_data$lopnr == indv]
      }
    }else{
      if(Settings$exclude_chemo_missing){
        d$EXC[ID] <- 10 + d$EXC[ID]
      }
    }
  }
  d$chemotherapy_adj[is.na(d$chemotherapy_adj)] <- "NA"
  d$chemotherapy_neoadj[is.na(d$chemotherapy_neoadj)] <- "NA"
  d$chemotherapy_status[is.na(d$chemotherapy_status)] <- "NA"
  d$chemotherapy_adj <- factor(d$chemotherapy_adj, levels = c("FALSE", "NA", "TRUE"))
  d$chemotherapy_neoadj <- factor(d$chemotherapy_neoadj, levels = c("FALSE", "NA", "TRUE"))
  d$chemotherapy_status <- factor(d$chemotherapy_status, levels = c("FALSE", "NA", "TRUE"))
  # ---===---===---===---===---===---===---===---===---===---===---===---===---===
  # Considering Risk Factors according to their date: ############################     ***** For RFs, Nulls are considered as FALSE *****
  RiskFactors <- c("Hypertension", "Diabetes")
  d$RF <- FALSE
  for (Risk in RiskFactors){
    Var <- paste("Date_",Risk, sep = "")
    d[[Var]] <- as.Date(d[[Var]], format = "%m/%d/%Y %H:%M")
    d$RF <- ifelse(!is.na(d[[Var]]) & d[[Var]] <= d$Date_of_BC, TRUE, d$RF) #       ***** For RFs, No Excludes *****
  }
  d$RF <- as.factor(d$RF)
  
  # Considering distance metastases: #############################################
  d$first_dist_metas_date <-as.Date(d$r_met_dat, format = "%m/%d/%Y %H:%M")
  d$dist_metas <- ifelse(d$r_met == "1" & !is.na(d$r_met), "TRUE", "FALSE")
  d$dist_metas <- as.factor(d$dist_metas)
  
  # Considering relapse: #########################################################
  d$relapse_date <-as.Date(d$r_rec_brostdat, format = "%m/%d/%Y %H:%M")
  d$Relapse <- ifelse(d$r_rec_brost == "1" & !is.na(d$r_rec_brost), "TRUE", "FALSE")
  d$Relapse <- as.factor(d$Relapse)

  
  # Considering CVDs according to their date (Before and After): #################
  MACE = c(
    "STEMI",
    "NSTEMI",
    "subMI",
    "HeartFailure",
    "Stroke"
  )
  
  d$MACE_ov_before <- FALSE
  d$MACE_sv_before <- FALSE
  
  d$MACE_ov_after <- FALSE
  d$MACE_sv_after <- FALSE
  d$first_MACE_ov_after <- NA
  d$first_MACE_sv_after <- NA
  
  for(CVD in MACE){
    Var_ov = paste(CVD, "_ov_dtc", sep = "")
    Var_sv = paste(CVD, "_sv_dtc", sep = "")
    
    Date = paste(CVD, "_date", sep = "")
    d[[Var_ov]] <- as.Date(d[[Var_ov]], format = "%m/%d/%Y %H:%M")
    d[[Var_sv]] <- as.Date(d[[Var_sv]], format = "%m/%d/%Y %H:%M")
    
    
    d$MACE_ov_before <- ifelse(!is.na(d[[Var_ov]]) & d[[Var_ov]] <= 
                                 d$Date_of_BC, TRUE, d$MACE_ov_before)
    d$MACE_sv_before <- ifelse(!is.na(d[[Var_sv]]) & d[[Var_sv]] <= 
                                 d$Date_of_BC, TRUE, d$MACE_sv_before)
    d$MACE_ov_after <- ifelse(!is.na(d[[Var_ov]]) & d[[Var_ov]] > 
                                d$Date_of_BC, TRUE, d$MACE_ov_after)
    d$MACE_sv_after <- ifelse(!is.na(d[[Var_sv]]) & d[[Var_sv]] > 
                                d$Date_of_BC, TRUE, d$MACE_sv_after)
    
    firstCondov <- !is.na(d[[Var_ov]]) & d[[Var_ov]] > d$Date_of_BC
    d$first_MACE_ov_after <- as.Date(ifelse(firstCondov&!is.na(d[[Var_ov]])&is.na(d$first_MACE_ov_after), d[[Var_ov]], 
                                            ifelse(firstCondov&!is.na(d[[Var_ov]])&!is.na(d$first_MACE_ov_after)&(d[[Var_ov]] <
                                                                                                                    d$first_MACE_ov_after), d[[Var_ov]], d$first_MACE_ov_after)))
    
    firstCondsv <- !is.na(d[[Var_sv]]) & d[[Var_sv]] > d$Date_of_BC
    d$first_MACE_sv_after <- as.Date(ifelse(firstCondsv&!is.na(d[[Var_sv]])&is.na(d$first_MACE_sv_after), d[[Var_sv]], 
                                            ifelse(firstCondsv&!is.na(d[[Var_sv]])&!is.na(d$first_MACE_sv_after)& (d[[Var_sv]] <
                                                                                                                     d$first_MACE_sv_after), d[[Var_sv]], d$first_MACE_sv_after)))
  }
  
  d$MACE_ov_before <- as.factor(d$MACE_ov_before)
  d$MACE_sv_before <- as.factor(d$MACE_sv_before)
  d$MACE_ov_after <- as.factor(d$MACE_ov_after)
  d$MACE_sv_after <- as.factor(d$MACE_sv_after)
  # Creating MACE_after for 15 CVDs: #############################################
  d$MACE_before <- ifelse(d$MACE_ov_before == "TRUE" | d$MACE_sv_before == "TRUE", TRUE, FALSE)
  d$MACE_before <- as.factor(d$MACE_before)
  
  d$MACE_after <- ifelse(d$MACE_ov_after == "TRUE" | d$MACE_sv_after == "TRUE", TRUE, FALSE)
  d$MACE_after <- as.factor(d$MACE_after)
  
  d$first_MACE_after <- as.Date(ifelse(is.na(d$first_MACE_sv_after), d$first_MACE_ov_after,
                                       ifelse(is.na(d$first_MACE_ov_after), d$first_MACE_sv_after,
                                              ifelse(d$first_MACE_sv_after <
                                                       d$first_MACE_ov_after, d$first_MACE_sv_after, d$first_MACE_ov_after))))
  
  # Considering CVDs according to their date (Before and After): #################
  CVDs = c(
    "STEMI",
    "NSTEMI",
    "subMI",
    "AnginaPectoris",
    "AtrialFibrilFlutter",
    "AtrioventBlock",
    "Venttachycardia",
    "VentFibrilFlutter",
    "HeartFailure",
    "Stroke",
    "Claudication",
    "CoronaryStent",
    "CABG",
    "Pacemaker",
    "CRT_ICD"
  )
  
  d$CVDs_ov_before <- FALSE
  d$CVDs_sv_before <- FALSE
  
  d$CVDs_ov_after <- FALSE
  d$CVDs_sv_after <- FALSE
  d$first_CVDs_ov_after <- NA
  d$first_CVDs_sv_after <- NA
  
  for(CVD in CVDs){
    Var_ov = paste(CVD, "_ov_dtc", sep = "")
    Var_sv = paste(CVD, "_sv_dtc", sep = "")
    
    Date = paste(CVD, "_date", sep = "")
    d[[Var_ov]] <- as.Date(d[[Var_ov]], format = "%m/%d/%Y %H:%M")
    d[[Var_sv]] <- as.Date(d[[Var_sv]], format = "%m/%d/%Y %H:%M")
    
    
    d$CVDs_ov_before <- ifelse(!is.na(d[[Var_ov]]) & d[[Var_ov]] <= 
                                 d$Date_of_BC, TRUE, d$CVDs_ov_before)
    d$CVDs_sv_before <- ifelse(!is.na(d[[Var_sv]]) & d[[Var_sv]] <= 
                                 d$Date_of_BC, TRUE, d$CVDs_sv_before)
    d$CVDs_ov_after <- ifelse(!is.na(d[[Var_ov]]) & d[[Var_ov]] > 
                                d$Date_of_BC, TRUE, d$CVDs_ov_after)
    d$CVDs_sv_after <- ifelse(!is.na(d[[Var_sv]]) & d[[Var_sv]] > 
                                d$Date_of_BC, TRUE, d$CVDs_sv_after)
    
    firstCondov <- !is.na(d[[Var_ov]]) & d[[Var_ov]] > d$Date_of_BC
    d$first_CVDs_ov_after <- as.Date(ifelse(firstCondov&!is.na(d[[Var_ov]])&is.na(d$first_CVDs_ov_after), d[[Var_ov]], 
                                            ifelse(firstCondov&!is.na(d[[Var_ov]])&!is.na(d$first_CVDs_ov_after)&(d[[Var_ov]] <
                                                                                                                    d$first_CVDs_ov_after), d[[Var_ov]], d$first_CVDs_ov_after)))
    
    firstCondsv <- !is.na(d[[Var_sv]]) & d[[Var_sv]] > d$Date_of_BC
    d$first_CVDs_sv_after <- as.Date(ifelse(firstCondsv&!is.na(d[[Var_sv]])&is.na(d$first_CVDs_sv_after), d[[Var_sv]], 
                                            ifelse(firstCondsv&!is.na(d[[Var_sv]])&!is.na(d$first_CVDs_sv_after)& (d[[Var_sv]] <
                                                                                                                     d$first_CVDs_sv_after), d[[Var_sv]], d$first_CVDs_sv_after)))
  }
  
  d$CVDs_ov_before <- as.factor(d$CVDs_ov_before)
  d$CVDs_sv_before <- as.factor(d$CVDs_sv_before)
  d$CVDs_ov_after <- as.factor(d$CVDs_ov_after)
  d$CVDs_sv_after <- as.factor(d$CVDs_sv_after)
  
  # Creating CVDs_after for 15 CVDs: #############################################
  d$CVDs_before <- ifelse(d$CVDs_ov_before == "TRUE" | d$CVDs_sv_before == "TRUE", TRUE, FALSE)
  d$CVDs_before <- as.factor(d$CVDs_before)
  
  d$CVDs_after <- ifelse(d$CVDs_ov_after == "TRUE" | d$CVDs_sv_after == "TRUE", TRUE, FALSE)
  d$CVDs_after <- as.factor(d$CVDs_after)
  
  d$first_CVDs_after <- as.Date(ifelse(is.na(d$first_CVDs_sv_after), d$first_CVDs_ov_after,
                                       ifelse(is.na(d$first_CVDs_ov_after), d$first_CVDs_sv_after,
                                              ifelse(d$first_CVDs_sv_after <
                                                       d$first_CVDs_ov_after, d$first_CVDs_sv_after, d$first_CVDs_ov_after))))
  
  # Considering ***All*** CVDs according to their date (Before): #################
  otherCVDs <- c(
    "Cardiac_Arrest_date",  
    "date_of_Cardiomyopati_I42_43_",  
    "Arteries_diseases_date",  
    "Atherosclerosis_date", ####
    "Pulmonary_embolism_date",
    "Venous_embolism_and_thrombosis_date",
    "date_of_VOC_I34_37_",  
    "Cerebrovascular_diseases_date",
    "Date_of_Arrytmi_I44_45_I47_49_",  ####
    "Ischemic_heart_disease_date" ####
  )
  
  d$otherCVDs_before <- FALSE
  d$otherCVDs_after <- FALSE
  d$first_otherCVDs_after <- NA
  
  for(CVD in otherCVDs){
    
    d[[CVD]] <- as.Date(d[[CVD]], format = "%m/%d/%Y %H:%M")
    
    d$otherCVDs_before <- ifelse(!is.na(d[[CVD]]) & d[[CVD]] 
                                 <= d$Date_of_BC, TRUE, d$otherCVDs_before)
    d$otherCVDs_after <- ifelse(!is.na(d[[CVD]]) & d[[CVD]]
                                > d$Date_of_BC, TRUE, d$otherCVDs_after)
    firstCond <- !is.na(d[[CVD]]) & d[[CVD]] > d$Date_of_BC
    d$first_otherCVDs_after <- as.Date(ifelse(firstCond&!is.na(d[[CVD]])&is.na(d$first_otherCVDs_after), d[[CVD]], 
                                              ifelse(firstCond&!is.na(d[[CVD]])&!is.na(d$first_otherCVDs_after)& (d[[CVD]] <
                                                                                                                    d$first_otherCVDs_after), d[[CVD]], d$first_otherCVDs_after)))
  }
  
  d$otherCVDs_before <- as.factor(d$otherCVDs_before)
  d$otherCVDs_after <- as.factor(d$otherCVDs_after)
  
  # Any CVDs: ####################################################################
  d$anyCVDs_before <- ifelse(d$CVDs_ov_before=="TRUE" | d$CVDs_sv_before=="TRUE"
                             | d$otherCVDs_before == "TRUE", TRUE, FALSE)
  d$anyCVDs_after <- ifelse(d$CVDs_ov_after=="TRUE" | d$CVDs_sv_after=="TRUE"
                            | d$otherCVDs_after == "TRUE", TRUE, FALSE)
  
  d$first_anyCVDs_after <- as.Date(ifelse(is.na(d$first_otherCVDs_after), d$first_CVDs_after,
                                          ifelse(is.na(d$first_CVDs_after), d$first_otherCVDs_after,
                                                 ifelse(d$first_otherCVDs_after <
                                                          d$first_CVDs_after, d$first_otherCVDs_after, d$first_CVDs_after))))
  
  
  d$anyCVDs_before <- as.factor(d$anyCVDs_before)
  d$anyCVDs_after <- as.factor(d$anyCVDs_after)
  
  # Considering the relapse and dist_metas: ######################################
  if (Settings$exclude_relapse_met){
    d$anyCVDs_after <- ifelse(d$anyCVDs_after == "TRUE" & (!is.na(d$relapse_date)&(d$first_anyCVDs_after > d$relapse_date) 
                                                           | !is.na(d$first_dist_metas_date)&(d$first_anyCVDs_after > d$first_dist_metas_date)), "FALSE", as.character(d$anyCVDs_after))
    d$anyCVDs_after <- as.factor(d$anyCVDs_after)
  }
  # Tumor Biology (for matching): ################################################
  # Stage:
  d$Stage <- NA
  
  # ---===---===---===---===---===---===---===---===---===---===---===---===---===
  # Tumor Biology, T class:                                                       ********** Tx ??? ************
  d$Tclass <- NA
  d$Tclass <- sub(" .*", "", d$a_tnm_tklass_Beskrivning)
  d$Tclass[d$Tclass=="T1a" | d$Tclass=="T1b" | d$Tclass=="T1c"] <- "T1"
  d$Tclass[d$Tclass=="T4a" | d$Tclass=="T4b" | d$Tclass=="T4c" | d$Tclass=="T4d"] <- "T4"
  
  d$Tclass <- ifelse(d$a_tnm_tklass_Beskrivning == "", NA, d$Tclass)  
  d <- d[!is.na(d$Tclass), ]
  if(Settings$exclude_tb_na){
    d$EXC[is.na(d$Tclass)] <- 1000+d$EXC[is.na(d$Tclass)]             
  }
  
  if(Settings$exclude_tis){
    d$EXC[!is.na(d$Tclass) & d$Tclass == "Tis"] <- 
      1000+d$EXC[!is.na(d$Tclass) & d$Tclass == "Tis"]             
  }
  
  if(Settings$exclude_t4_spo){
    d$EXC[d$Tclass == "T4" & !is.na(d$Tclass)] <-
      1000+d$EXC[d$Tclass == "T4" & !is.na(d$Tclass)]       # ********** Because of errors caused by different factors BRCA/Sporadic ************
  }
  
  d$Tclass <- factor(d$Tclass)
  summary(d$Tclass, useNA = "always")
  # ---===---===---===---===---===---===---===---===---===---===---===---===---===
  
  # Tumor Biology, N class:                                                       ********** Nx ??? ************
  d$Nclass <- NA
  d$Nclass <- d$a_tnm_nklass_Beskrivning
  d$Nclass <- sub(" .*", "", d$a_tnm_nklass_Beskrivning)
  
  d$Nclass <- ifelse(d$a_tnm_nklass_Beskrivning == "", NA, d$Nclass)  
  d <- d[!is.na(d$Nclass), ]
  
  d$Nclass[d$Nclass=="N1a" | d$Nclass=="N1b"] <- "N1"
  d$Nclass[d$Nclass=="N2a"] <- "N2"                                             # ***** Because of errors caused by different factors BRCA/Sporadic *****
  
  if (Settings$exclude_n2_n3){
    d$EXC[!is.na(d$Nclass) & d$Nclass == "N2"] <-                               # ****** N2 and N3 was excluded peacuse of very low occurrence *******
      1000+d$EXC[!is.na(d$Nclass) & d$Nclass == "N2"]
    d$EXC[!is.na(d$Nclass) & d$Nclass == "N3"] <- 
      1000+d$EXC[!is.na(d$Nclass) & d$Nclass == "N3"]
  }
  
  if(Settings$exclude_tb_na){
    d$EXC[is.na(d$Nclass)] <- 1000+d$EXC[is.na(d$Nclass)]
  }
  
  d$Nclass <- factor(d$Nclass)
  summary(d$Nclass)
  # ---===---===---===---===---===---===---===---===---===---===---===---===---===
  
  # Tumor Biology, M class:                                                       ********** Mx and M1 Excluded, M is Not Used for analysis ************
  d$Mclass <- NA
  d$Mclass <- d$a_tnm_mklass_Beskrivning
  d$Mclass <- sub(" .*", "", d$a_tnm_mklass_Beskrivning)
  
  d$Mclass <- ifelse(d$a_tnm_mklass_Beskrivning == "", NA, d$Mclass)  
  d <- d[!is.na(d$Mclass), ]
  
  if(Settings$exclude_tb_na){
    d$EXC[is.na(d$Mclass)] <- 1000+d$EXC[is.na(d$Mclass)]
  }
  
  if(Settings$exclude_m1){
    d$EXC[!is.na(d$Mclass) & (d$M_Tumorutbredning_fjarrmetastaser_Cancer_type_1 == "M1")] <- 
      1000+d$EXC[!is.na(d$Mclass) & (d$M_Tumorutbredning_fjarrmetastaser_Cancer_type_1 == "M1")] 
  }
  
  d$Mclass <- factor(d$Mclass)
  summary(d$Mclass)
  
  # Radio Therapy (RT): ##########################################################
  d$RT <- NA
  d$RT <- ifelse(d$post_rt_Varde == 1 | d$pre_rt_Varde == 1, "TRUE",
                 ifelse(is.na(d$post_rt_Varde) & d$pre_rt_Varde == 0, "FALSE", 
                        ifelse(is.na(d$post_rt_Varde) & is.na(d$pre_rt_Varde), NA,
                               ifelse(d$post_rt_Varde == 0, "FALSE", d$post_rt_Varde))))
  
  d$RT <- ifelse(is.na(d$RT) & (d$a_planbeh_rt_Varde == 1 | d$op_planbeh_rt_Varde == 1), "TRUE", 
                 ifelse(is.na(d$RT) & d$a_planbeh_rt_Varde == 0 | d$op_planbeh_rt_Varde == 0, "FALSE", d$RT))
  if(Settings$exclude_rt_na){
    d$EXC[is.na(d$RT)] <- 10000 + d$EXC[is.na(d$RT)]  #   ********* Excluding NA RTs *****************
    
  }
  d$RT <- as.factor(d$RT)
  summary(d$RT)
  
  
  # Anthracycline Therapy: #######################################################     ********* Should be reconsidered for 0 **********
  # d$Anth <- ifelse(d$pre_kemo_antra == 1 | d$post_kemo_antra == 1, TRUE, 
  #                  ifelse(d$pre_kemo_antra == 0 | d$post_kemo_antra == 0, FALSE, "NA"))
  # d$Anth[d$chemotherapy_status == "FALSE"] <- FALSE
  # d$Anth <- factor(d$Anth)
  # summary(d$Anth)
  
  # HER2 Therapy (Antikropp): ####################################################    ********* Should be reconsidered for 0 **********
  d$Antikropp <- ifelse(!is.na(d$pre_antikropp_Varde) & d$pre_antikropp_Varde == 1 | !is.na(d$post_antikropp_Varde) & d$post_antikropp_Varde == 1, TRUE,
                        ifelse(!is.na(d$pre_antikropp_Varde) & d$pre_antikropp_Varde == 0 | !is.na(d$post_antikropp_Varde) & d$post_antikropp_Varde == 0, FALSE, "NA"))
  d$Antikropp <- factor(d$Antikropp)
  summary(d$Antikropp)
  
  # Ki67: ########################################################################    ********* Lots of missing values, excluded **********
  # d$Ki67_corrected <- NA
  # d$Ki67_corrected <- ifelse(
  #   d$Ki67 == "Null" | d$Ki67 == "97" | d$Ki67 == "98",
  #   NA,  # If the condition is true, assign 1
  #   d$Ki67   # If the condition is false (NA), assign 0
  # )
  # d$Ki67_corrected <- factor(d$Ki67_corrected)
  # summary(d$Ki67_corrected)
  
  # ErposYN:  ####################################################################
  # d$Erpos <- NA
  # d$Erpos <- ifelse(d$a_pad_er_Varde == 1, TRUE, ifelse(
  #   d$ErposYN == 0, FALSE, "NA"))
  # d$Erpos <- factor(d$Erpos)
  # summary(d$Erpos)
  # 
  # PrposYN: #####################################################################
  # d$Prpos <- NA
  # d$Prpos <- ifelse(d$PrposYN == 1, TRUE, ifelse(
  #   d$PrposYN == 0, FALSE, "NA"))
  # d$Prpos <- factor(d$Prpos)
  # summary(d$Prpos)
  
  # Her2YN: ######################################################################    ******** Missing values *******
  # d$Her2 <- NA
  # d$Her2 <- ifelse(d$Her2YN == 1, TRUE, ifelse(
  #   d$Her2YN == 0, FALSE, "NA"))
  # # d$Her2[d$Antikropp == FALSE] <- FALSE
  # d$Her2 <- factor(d$Her2)
  # summary(d$Her2)
  # 
  # Follow-up duration: ##########################################################
  d$Date_of_death <- as.Date(d$AVLIDDAT , format = "%m/%d/%Y")
  d$Date_of_Birth <- as.Date(d$Fodelsedatum, format = "%m/%d/%Y")
  d$Age_at_death_decimal <- time_length(interval(d$Date_of_Birth, d$Date_of_death), "years")
  # d$fup <- as.numeric(difftime(FUendDate, d$Date_of_Birth))
  d$fup <- as.numeric(difftime(FUendDate, d$Date_of_BC, units = "days"))
  
  
  d$FU <- NA
  Event <- d[[Settings$target_event]]
  
  # FU_death <- ifelse(!is.na(d$Date_of_death), difftime(d$Date_of_death, d$Date_of_Birth, units = "days"), d$fup)
  FU_death <- ifelse(!is.na(d$Date_of_death), difftime(d$Date_of_death, d$Date_of_BC, units = "days"), d$fup)
  
  # FU_relapse <- ifelse(!is.na(d$relapse_date), difftime(d$relapse_date, d$Date_of_BC, units = "days"), d$fup)
  # FU_meta <- ifelse(!is.na(d$first_dist_metas_date), difftime(d$first_dist_metas_date, d$Date_of_BC, units = "days"), d$fup)
  # FU_CVD <- ifelse(!is.na(Event), difftime(Event, d$Date_of_Birth, units = "days"), d$fup)
  FU_CVD <- ifelse(!is.na(Event), difftime(Event, d$Date_of_BC, units = "days"), d$fup)
  
  d$FU <- pmin(FU_death, FU_CVD)
  # d$FU <- pmin(FU_death, FU_CVD, FU_meta)
  
  # Considering Cause of death 
  d$CVDdeath <- FALSE
  for(i in 1:length(d$lopnr)){
    # d$CVDdeath[i] <- any(strsplit(d$Cause_of_death_ICD_code[i],"")[[1]]=="I")
    if (any(CauseOfDeath$lopnr%in%d$lopnr[i]) ) {
      d$CVDdeath[i] <- any(CauseOfDeath$CVDdeath[which(CauseOfDeath$lopnr%in%d$lopnr[i])])
    }
  }
  
  d$FU_outcome <- ifelse(d$FU == FU_CVD & !is.na(Event), 1, 0)
  d$FU_outcome[d$CVDdeath&d$FU_outcome == 0 & FU_death==d$FU] <- 1
  d$CompetingOutcome <- ifelse(!is.na(d$Date_of_death) & d$FU_outcome==0, 2, d$FU_outcome)
  # d$CompetingOutcome <- ifelse(d$dist_metas=="TRUE" & d$FU_outcome==0 & d$FU == FU_meta, 2, d$CompetingOutcome)
  
  
  d$Age_at_event <- time_length(interval(d$Date_of_Birth, Event[d$lopnr%in%d$lopnr]), "years")
  d$Age_decimal <- time_length(interval(d$Date_of_Birth, d$Date_of_BC), "years")
  d$Age_at_event_corrected <- ifelse(d$FU_outcome == 1, d$Age_at_event, NA)
  # d$Age_at_end_of_FU <- ifelse(!is.na(d$Age_at_event_corrected), d$Age_at_event_corrected,
  #                               ifelse(!is.na(d$Date_of_death), time_length(interval(d$Date_of_Birth, d$Date_of_death), "years"),
  #                                      time_length(interval(d$Date_of_Birth, FUendDate), "years")))
  d$Age_at_end_of_FU <- ifelse(!is.na(d$Age_at_event_corrected), d$Age_at_event_corrected,
                               # ifelse(!is.na(d$first_dist_metas_date), time_length(interval(d$Date_of_Birth, d$first_dist_metas_date), "years"),
                               ifelse(!is.na(d$Date_of_death), time_length(interval(d$Date_of_Birth, d$Date_of_death), "years"),
                                      time_length(interval(d$Date_of_Birth, FUendDate), "years")))
  d$FU_years <- d$Age_at_end_of_FU - d$Age_decimal
  d$Age_at_end_rounded <- round(d$Age_at_end_of_FU,digits = 0)
  # d$Stage[is.na(d$Stage)] <- 50
  d$Tclass[is.na(d$Tclass)] <- "Tx"
  d$Nclass[is.na(d$Nclass)] <- "Nx"
  d$Mclass[is.na(d$Mclass)] <- "Mx"
  d$Stage <- relevel(d$Stage, ref = "1")
  
  # Categorizing continuous variables #############################################
  d$Age_cat <- cut(d$Age,
                   breaks = c(0, 35, 60, Inf),
                   labels = c("<35", "35-60", "60+"),
                   right = FALSE)
  
  d$Age_at_end_cat <- cut(d$Age_at_end_of_FU,
                          breaks = c(0, 35, 60, Inf),
                          labels = c("<35", "35-60", "60+"),
                          right = FALSE)
  d$FU_cat <- cut(d$FU,
                  breaks = c(summary(d$FU)[1], summary(d$FU)[3], summary(d$FU)[6]),
                  labels = c("1st", "2nd"),
                  right = FALSE)
  
  

  md <- d[d$group=="Sporadic", ]
 
  
  return(md)
}