library(dplyr)

# ---- Constants ----
DATA_DIR   <- "/Users/mehman/Projects/0_Reference_Data/0_otherRegions"
DATE_START <- as.Date("2008-01-01")
DATE_END   <- as.Date("2020-01-01")
FAR_FUTURE <- as.Date("2100-12-31")
REGION     <- "Region Uppsala/Örebro"

CVDs = c(
  "I21", #"STEMI" & "NSTEMI"
  "I22", #"subMI"
  "I20", #"AnginaPectoris"
  "I48", #"AtrialFibrilFlutter"
  "I44", #"AtrioventBlock"
  "I47", #"Venttachycardia"
  "I49", #"VentFibrilFlutter"
  "I50", #"HeartFailure"
  "I61", #"Stroke"
  "I70", #"Claudication"
  "Z95" #"CoronaryStent", "CABG", "Pacemaker", "CRT_ICD"
)

otherCVDs <- c(
  "I46", #"Cardiac_Arrest_date"
  "I42", "I43",  #"date_of_Cardiomyopati_I42_43_"
  "I71", "I72", "I73", "I74", "I75",  #"Arteries_diseases_date"
  "I70", #"Atherosclerosis_date"
  "I26", #"Pulmonary_embolism_date"
  "I81", "I82", #"Venous_embolism_and_thrombosis_date"
  "I34", "I35", "I36", "I37", #"date_of_VOC_I34_37_",  
  "I60", "I61", "I62", "I63", "I64", "I65", "I66", "I67", # "Cerebrovascular_diseases_date"
  "I44", "I45", "I47", "I48", "I49", #"Date_of_Arrytmi_I44_45_I47_49_",  ####
  "I20", "I21", "I22", "I23", "I24", "I25" #  "Ischemic_heart_disease_date" ####
)

any_cvds <- unique(c(CVDs, other_cvds))

MACE <- c("I21", "I22", "I50", "I61")

# ---- Helper functions ----
load_rdata <- function(file, obj_name) {
  env <- new.env()
  load(file.path(DATA_DIR, file), envir = env)
  env[[obj_name]]
}

parse_indate <- function(df) {
  df$INDATUMA <- dplyr::if_else(
    is.na(df$INDATUMA),
    as.Date(paste0(df$AR, "-07-01")),
    as.Date(df$INDATUMA, format = "%Y%m%d")
  )
  df
}

has_cvd_code <- function(codes, code_set) {
  any(sapply(code_set, function(x) any(startsWith(codes, x))))
}

earliest_mace_date <- function(df, diag_col, date_col, after_date) {
  after    <- df[df[[date_col]] > after_date, ]
  hit_rows <- which(vapply(after[[diag_col]], function(x) any(startsWith(x, MACE)), logical(1)))
  
  if (length(hit_rows) == 0) return(list(flag = 0L, date = FAR_FUTURE))
  list(flag = 1L, date = min(after[[date_col]][hit_rows]))
}

# ---- Load and filter diagnosis index ----
date_of_diag <- load_rdata("dt_fall.RData", "dt_fall") %>%
  transmute(
    LopNr,
    Date_of_BC  = as.Date(a_diag_dat, format = "%m/%d/%Y %H:%M"),
    REGION_NAMN
  ) %>%
  filter(
    REGION_NAMN == REGION,
    Date_of_BC  >= DATE_START,
    Date_of_BC  <  DATE_END
  )

cohort_ids <- date_of_diag$LopNr

# ---- Load, parse, and filter register data ----
ov <- load_rdata("patregoppen_fall.RData", "patregoppen_fall") %>%
  parse_indate() %>%
  filter(LopNr %in% cohort_ids)

sv <- load_rdata("patregsluten_2008_2019_fall.RData", "patregsluten_2008_2019_fall") %>%
  parse_indate() %>%
  mutate(UTDATUMA = as.Date(UTDATUMA, format = "%Y%m%d")) %>%
  filter(LopNr %in% cohort_ids)

sv_pre2008 <- load_rdata("patregsluten_1987_2007_fall.RData", "patregsluten_1987_2007_fall") %>%
  parse_indate() %>%
  mutate(UTDATUMA = as.Date(UTDATUMA, format = "%Y%m%d")) %>%
  filter(LopNr %in% cohort_ids)

# ---- Derive CVD status per patient ----
CVDstatus <- data.frame(
  LopNr          = cohort_ids,
  any_cvd_before = 0L,
  MACE           = 0L,
  MACE_date      = as.Date(NA)
)

for (pat in cohort_ids) {
  print(pat)
  bc_date <- min(date_of_diag$Date_of_BC[date_of_diag$LopNr == pat])
  
  # Inpatient 2008–2019 (single diagnosis column)
  sv_pat    <- sv[sv$LopNr == pat, c("hdia", "INDATUMA")] %>% filter(hdia != "")
  sv_before <- filter(sv_pat, INDATUMA <= bc_date)
  sv_mace   <- earliest_mace_date(sv_pat, "hdia", "INDATUMA", bc_date)
  
  # Outpatient (wide diagnosis columns)
  ov_pat       <- ov[ov$LopNr == pat, c(1, 7:37)]
  ov_before    <- filter(ov_pat, INDATUMA <= bc_date)
  ov_mace      <- earliest_mace_date(ov_pat, "hdia", "INDATUMA", bc_date)
  
  # Inpatient pre-2008 (wide diagnosis columns)
  sp_pat       <- sv_pre2008[sv_pre2008$LopNr == pat, c(1, 7, 10:39)]
  sp_before    <- filter(sp_pat, INDATUMA <= bc_date)
  sp_mace      <- earliest_mace_date(sp_pat, "hdia", "INDATUMA", bc_date)
  
  # Any CVD before breast cancer diagnosis
  CVDstatus$any_cvd_before[CVDstatus$LopNr == pat] <- as.integer(
    has_cvd_code(sv_before$hdia,  any_cvds) |
      has_cvd_code(ov_before$hdia,  any_cvds) |
      has_cvd_code(sp_before$hdia,  any_cvds)
  )
  
  # MACE after breast cancer diagnosis
  if (any(sv_mace$flag, ov_mace$flag, sp_mace$flag)) {
    CVDstatus$MACE[CVDstatus$LopNr == pat]      <- 1L
    CVDstatus$MACE_date[CVDstatus$LopNr == pat] <- as.Date(min(sv_mace$date, ov_mace$date, sp_mace$date))
  }
}

