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

# ---- sanity checks ----
stopifnot(all(c("LopNr", "dors_AR", "dors_DODSDAT", "dors_alder") %in% names(db)))

# ---- inspect (optional, for debugging only) ----
print(names(db))
print(head(db$dors_AR, 10))
print(head(db$dors_DODSDAT, 10))

# ---- build output cleanly ----
cause_of_death <- db |>
  transmute(
    lopnr    = LopNr,
    AR       = dors_AR,
    DODSDAT  = dors_DODSDAT,
    alder    = dors_alder
  )

db_morsak <- db |> select(contains("MORSAK"))


mask <- apply(db_morsak, 1:2, function(x) grepl("^I", x, ignore.case = FALSE))
mask[is.na(mask)] <- FALSE  # optional: treat NAs as no-match

cause_of_death$CVDdeath <- apply(mask,MARGIN = 1, any)
save(cause_of_death, file = "../0_Reference_Data/0_otherRegions/Processed/AnyCVDCauseOfDeath.RData")
