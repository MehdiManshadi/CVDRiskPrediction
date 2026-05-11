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

# ------------------------------------------------------------------------------
# Load and prepare drug registry data
# ------------------------------------------------------------------------------
load(file.path(REF_DATA_DIR, "lakmedreg_fall.RData"))
drugData <- lakmedreg_fall[, c("LopNr", "ATC", "FDATUM", "EDATUM")]

drugData$FDATUM <- as.Date(drugData$FDATUM)
drugData$EDATUM <- as.Date(drugData$EDATUM)

# ------------------------------------------------------------------------------
# Initialize drug exposure matrices
# ------------------------------------------------------------------------------
uniqueATC    <- unique(drugData$ATC)
uniqueLopnr  <- unique(Sporadic$LopNr)
n_indv       <- length(uniqueLopnr)
n_atc        <- length(uniqueATC)

.init_matrix <- function(lopnr_vec, atc_vec) {
  mat           <- as.data.frame(matrix(nrow = length(lopnr_vec), ncol = length(atc_vec)))
  colnames(mat) <- atc_vec
  rownames(mat) <- lopnr_vec
  mat$lopnr     <- lopnr_vec
  mat$anyData   <- FALSE
  mat
}

medicineData <- .init_matrix(uniqueLopnr, uniqueATC)  # Binary exposure flag
medicineFreq <- .init_matrix(uniqueLopnr, uniqueATC)  # Prescription frequency
medicineFT   <- .init_matrix(uniqueLopnr, uniqueATC)  # Days from first use to diagnosis
medicineLT   <- .init_matrix(uniqueLopnr, uniqueATC)  # Days from last use to diagnosis

# ------------------------------------------------------------------------------
# Populate matrices: loop over individuals
# ------------------------------------------------------------------------------
for (indv in uniqueLopnr) {
  if (!indv %in% drugData$LopNr) next

  diagDate <- as.Date(Sporadic$Date_of_BC[Sporadic$LopNr == indv])


  # Set row defaults to 0 and flag as having data
  for (mat in list(medicineData, medicineFreq, medicineFT, medicineLT)) {
    mat[mat$lopnr == indv, uniqueATC] <- 0
    mat$anyData[mat$lopnr == indv]    <- TRUE
  }

  # Pre-diagnosis drug records
  preDiagDrugs <- drugData[drugData$LopNr == indv & drugData$FDATUM < diagDate, ]
  drugCounts   <- table(preDiagDrugs$ATC)

  if (length(drugCounts) == 0) {
    cat(sprintf("[%5.1f%%] Individual %s: no pre-diagnosis records\n",
                which(indv == uniqueLopnr) / n_indv * 100, indv))
    next
  }

  atcNames  <- names(drugCounts)
  firstUse  <- numeric(length(atcNames))
  lastUse   <- numeric(length(atcNames))

  for (i in seq_along(atcNames)) {
    dates        <- drugData$FDATUM[drugData$LopNr == indv & drugData$ATC == atcNames[i]]
    preDiagDates <- dates[dates < diagDate]
    firstUse[i]  <- as.numeric(difftime(diagDate, min(dates),         units = "days"))
    lastUse[i]   <- as.numeric(difftime(diagDate, max(preDiagDates),  units = "days"))
  }

  medicineData[medicineData$lopnr == indv, atcNames] <- 1
  medicineFreq[medicineFreq$lopnr == indv, atcNames] <- as.integer(drugCounts)
  medicineFT  [medicineFT$lopnr   == indv, atcNames] <- firstUse
  medicineLT  [medicineLT$lopnr   == indv, atcNames] <- lastUse

  cat(sprintf("[%5.1f%%] Individual %s processed\n",
              which(indv == uniqueLopnr) / n_indv * 100, indv))
}

# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------
save(medicineData, file = file.path(OUTPUT_DIR, "medicineData.RData"))
save(medicineLT,   file = file.path(OUTPUT_DIR, "medicineLT.RData"))
save(medicineFreq, file = file.path(OUTPUT_DIR, "medicineFreq.RData"))
save(medicineFT,   file = file.path(OUTPUT_DIR, "medicineFT.RData"))

message("All drug exposure matrices saved successfully.")