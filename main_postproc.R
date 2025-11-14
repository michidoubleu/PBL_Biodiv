###############################################################
# Main Script for PBL_Biodiv Pipeline
# -------------------------------------------------------------
# Processes G4M → DS biodiversity link outputs:
#   * Reads scenario mapping
#   * Loads link-output RData files for selected scenario
#   * Applies two-stage land-cover (LC) mapping
#   * Aggregates link outputs
#   * Generates CSVs for netCDF conversion
#   * Optionally writes netCDF
###############################################################

library(dplyr)
library(tidyr)
library(data.table)

###############################################################
# User parameters
###############################################################

create.masked <- TRUE
create.global <- TRUE
write.nc <- FALSE    # netCDF writing controlled at the end

###############################################################
# Path Settings
# NOTE: This is a network location used within PBL
###############################################################

define.path <- "P:/globiom/Projects/SSPs/ScenarioMIP7/Results/lookup_table_results_withSSPV3/lookup_table_v5_16Jun2025/raw_results/Biodiversity_Link"
template.path <- "P:/globiom/Projects/PBL_BIODIV_2025/Postprocessing_ncdf/template"

###############################################################
# Load scenario mapping
###############################################################

scenarios <- read.csv(paste0(define.path, "/scenario_mapping_BetterLookup2SSP2a_20250613_2SCENS.csv"))
scen.setting <- NULL

###############################################################
# Utility Function: Filter files by numeric suffix
###############################################################

filter_paths <- function(paths, keep, digits = 6) {
  ids <- sub(".*(?:\\.|_)(\\d+)\\.RData$", "\\1", paths, perl = TRUE)
  ids[grepl("^[0-9]+$", ids) == FALSE] <- NA
  ids_num <- as.numeric(ids)
  mod <- 10^digits
  last_digits <- ids_num %% mod
  idx <- !is.na(last_digits) & (last_digits %in% keep)
  paths[idx]
}

###############################################################
# Select scenario to run
###############################################################

curr.scen <- 3     # <-- loop index (single-run for now)

curr.loops <- scenarios$ScenNr[scenarios$ScenLoop == curr.scen]

# Make scenario name + record scenario settings
sc_row <- unique(subset(scenarios, ScenLoop == curr.scen,
                        select = c("SCEN1", "SCEN2", "SCEN3")))
Scen <- paste(sc_row, collapse = "-")

curr.scen <- data.frame(Scen = Scen, sc_row, row.names = NULL)
scen.setting <- scen.setting %>% bind_rows(curr.scen)

###############################################################
# Read RData results to be processed
###############################################################

all.files <- list.files(define.path, full.names = TRUE)
filtered.files <- filter_paths(all.files, curr.loops)

# keep only files from the latest timestamp
filtered.files <- filtered.files[
  substr(filtered.files, nchar(filtered.files)-17, nchar(filtered.files)-13) ==
    as.character(max(as.numeric(substr(filtered.files,
                                       nchar(filtered.files)-17,
                                       nchar(filtered.files)-13))))
]

###############################################################
# Load link result list objects and merge
###############################################################

df.link.results <- NULL

for (jjj in seq_along(filtered.files)) {
  temp <- readRDS(filtered.files[jjj])
  df.link.results <- df.link.results %>% bind_rows(temp[[5]])
}

###############################################################
# Land-Cover Mapping (two-stage)
###############################################################

mapping_LC_names_1 <- readRDS(paste0(template.path,"/mapping_for_G4MDSlink.RData"))[[3]]
mapping_LC_names_2 <- readRDS(paste0(template.path,"/mapping_for_G4MDSlink.RData"))[[4]]

setDT(df.link.results)
mapping1 <- as.data.table(mapping_LC_names_1)
mapping2 <- as.data.table(mapping_LC_names_2)

results2 <- df.link.results

### Stage 1 — mapping1

# map lu.from
results2 <- mapping1[results2, on = .(lu.linkoutput = lu.from)]
results2[, lu.from := lu.new]
results2[, lu.new := NULL]

# map lu.to
results2 <- mapping1[results2, on = .(lu.linkoutput = lu.to)]
results2[, lu.to := lu.new]
results2[, lu.new := NULL]

### Stage 2 — protected-area mapping2

# map lu.from
results2 <- mapping2[results2, on = .(lu.new = lu.from)]
results2[, lu.from := lu.final]
results2[, lu.final := NULL]

# map lu.to
results2 <- mapping2[results2, on = .(lu.new = lu.to)]
results2[, lu.to := lu.final]
results2[, lu.final := NULL]

### Stage 3 — aggregate
results3 <- results2[
  , .(value = sum(value) * 1000),
  by = .(REGION, times, ns, lu.to, lu.from)
]

results3 <- as.data.frame(results3)

###############################################################
# Prepare for CSV Output (following Prep_CSV_for_netcdf structure)
###############################################################

full_simu_map <- read.csv(paste0(template.path,"/full_simu_map_biodiv.csv"), stringsAsFactors = FALSE)

mapping_simuID <- full_simu_map %>%
  dplyr::select(SimUID, country, REGION_37, colrowID) %>%
  rename(COUNTRY = country, REGION = REGION_37) %>%
  mutate(SimUID = as.numeric(as.character(SimUID)))

REGION_AG_Array <- c(
  "ArgentinaReg", "AustraliaReg", "BrazilReg", "CanadaReg", "ChinaReg",
  "CongoBasin", "EU_Baltic", "EU_CentralEast", "EU_MidWest", "EU_North",
  "EU_South", "Former_USSR", "IndiaReg", "IndonesiaReg", "JapanReg",
  "MalaysiaReg", "MexicoReg", "MiddleEast", "NewZealandReg", "NorthernAf",
  "Pacific_Islands", "RCAM", "RCEU", "ROWE", "RSAM", "RSAS", "RSEA_OPA",
  "RSEA_PAC", "RussiaReg", "SouthAfrReg", "SouthKorea", "EasternAf",
  "SouthernAf", "WesternAf", "TurkeyReg", "UkraineReg", "USAReg"
)

linking_out <- results3

###############################################################
# Compile LC CSV
###############################################################

# LC in base year
linking.resultLC.2000 <- linking_out %>%
  subset(times %in% c("2010")) %>%
  group_by(ns, times, lu.from) %>%
  summarize(value = sum(value), .groups = "keep") %>%
  mutate(times = "2000",
         times = as.integer(times)) %>%
  rename(LC = lu.from, SimUID = ns, Year = times)

# LC transitions
linking.resultLC <- linking_out %>%
  group_by(ns, times, lu.to) %>%
  summarize(value = sum(value), .groups = "keep") %>%
  rename(LC = lu.to, SimUID = ns, Year = times) %>%
  bind_rows(linking.resultLC.2000) %>%
  spread(key = LC, value = value) %>%
  mutate(SimUID = as.numeric(as.character(SimUID))) %>%
  mutate(
    urban = 0,
    priforest = ifelse(is.na(priforest), 0, priforest),
    mngforest = ifelse(is.na(mngforest), 0, mngforest),
    cropland = ifelse(is.na(cropland), 0, cropland),
    grassland = ifelse(is.na(grassland), 0, grassland),
    other = ifelse(is.na(other), 0, other),
    protected_other = 0,
    protected_priforest = 0
  )

linking.resultLC <- linking.resultLC %>% mutate(restored = 0)
if (!"SRP" %in% names(linking.resultLC)) linking.resultLC$SRP <- 0

linking.resultLC$Area <- rowSums(linking.resultLC[, -c(1:2)])

# Add colrow + region information
linking.resultLC.final <- linking.resultLC %>%
  left_join(mapping_simuID %>% dplyr::select(-REGION), by = "SimUID") %>%
  rename(Colrow = colrowID) %>%
  dplyr::select(
    SimUID, Area, Year, COUNTRY, Colrow,
    cropland, grassland, priforest, mngforest,
    SRP, restored, other,
    protected_priforest, protected_other, urban
  ) %>%
  arrange(Year, SimUID)

###############################################################
# LUC transitions
###############################################################

Array_fullLUC_list <- c(
  "cropland.cropland", "cropland.grassland", "cropland.mngforest",
  "cropland.SRP", "cropland.restored", "cropland.other",
  "grassland.grassland", "grassland.cropland", "grassland.mngforest",
  "grassland.SRP", "grassland.restored", "grassland.other",
  "priforest.priforest", "priforest.cropland", "priforest.grassland",
  "priforest.mngforest", "mngforest.mngforest", "SRP.SRP",
  "SRP.mngforest", "SRP.other", "restored.restored", "other.other",
  "other.cropland", "other.grassland", "other.mngforest", "other.SRP",
  "protected_priforest.protected_priforest",
  "protected_other.protected_other", "urban.urban"
)

# prep
ds.resultLUC0 <- linking_out %>%
  mutate(LUC = paste0(lu.from, ".", lu.to)) %>%
  dplyr::select(-lu.from, -lu.to, -REGION) %>%
  left_join(
    mapping_simuID[, c("SimUID", "REGION")] %>% rename(ns = SimUID) %>% mutate(ns = as.character(ns)),
    by = "ns"
  ) %>%
  spread(key = LUC, value = value) %>%
  mutate(urban.urban = 0) %>%
  ungroup()

# ensure all LUC columns exist
for (k in seq_along(Array_fullLUC_list)) {
  if (!(Array_fullLUC_list[k] %in% colnames(ds.resultLUC0))) {
    ds.resultLUC0[[Array_fullLUC_list[k]]] <- 0
  }
}

ds.resultLUC1 <- ds.resultLUC0 %>%
  rename(SimUID = ns, Year = times) %>%
  mutate(SimUID = as.numeric(as.character(SimUID))) %>%
  dplyr::select(SimUID, Year, all_of(Array_fullLUC_list))

linking.resultLCLUC.final <- linking.resultLC.final %>%
  left_join(ds.resultLUC1, by = c("SimUID", "Year"))

linking.resultLCLUC.final[is.na(linking.resultLCLUC.final)] <- 0

###############################################################
# Write CSV Outputs
###############################################################

write.csv(
  linking.resultLCLUC.final %>% arrange(Year, SimUID),
  file = paste0("./output/", Scen, "_", Sys.Date(), "_LULUC.csv"),
  row.names = FALSE
)

write.csv(
  scen.setting,
  file = paste0("./scen_setting_", Sys.Date(), ".csv"),
  row.names = FALSE
)

###############################################################
# Optional: netCDF writing (external script)
###############################################################

if (write.nc) {
  source("codes/results2netcdf_halfdegree_MW.R")
}

###############################################################
# End of main script
###############################################################
