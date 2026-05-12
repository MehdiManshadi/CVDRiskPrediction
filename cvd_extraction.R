library(dplyr)
library(haven)

# ---- function to load RData safely ----
load_rdata <- function(path) {
  env <- new.env()
  load(path, envir = env)
  as.list(env)
}

# ---- load data ----
data_list <- load_rdata("/Users/mehman/Projects/0_Reference_Data/0_otherRegions/dt_fall.RData")
db <- data_list$dt_fall
date_of_diag <- data.frame(LopNr = db$LopNr)
date_of_diag$Date_of_BC <- as.Date(db$a_diag_dat, format = "%m/%d/%Y %H:%M")
date_of_diag$REGION_NAMN <- db$REGION_NAMN
date_of_diag <- date_of_diag %>% filter(REGION_NAMN == "Region Uppsala/Örebro")

date_of_diag <- date_of_diag[
  date_of_diag$Date_of_BC >= as.Date("2008-01-01") &
    date_of_diag$Date_of_BC <  as.Date("2020-01-01"),
]

rm(data_list)
rm(db)


data_list <- load_rdata("/Users/mehman/Projects/0_Reference_Data/0_otherRegions/patregoppen_fall.RData")
ov <- data_list$patregoppen_fall
rm(data_list)

data_list <- load_rdata("/Users/mehman/Projects/0_Reference_Data/0_otherRegions/patregsluten_2008_2019_fall.RData")
sv <- data_list$patregsluten_2008_2019_fall
rm(data_list)
ov$INDATUMA <- as.Date(ifelse(is.na(ov$INDATUMA), as.Date(paste0(ov$AR, "-07-01")), 
                              as.Date(ov$INDATUMA, format = "%Y%m%d")))


sv$INDATUMA <- as.Date(sv$INDATUMA, format = "%Y%m%d")
sv$UTDATUMA <- as.Date(sv$UTDATUMA, format = "%Y%m%d")


sv <- sv %>% filter(LopNr %in% data)

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

MACE = c(
  "STEMI",
  "NSTEMI",
  "subMI",
  "HeartFailure",
  "Stroke"
)