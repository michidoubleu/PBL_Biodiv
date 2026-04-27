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
source("./codes/helper_functions.R")

###############################################################
# User parameters
###############################################################
date.tag <- "2025-12-10"
correct.G4M <- TRUE
create.masked <- TRUE
create.global <- TRUE
write.nc <- TRUE    # netCDF writing controlled at the end
amanda.files <- FALSE

###############################################################
# Path Settings
# NOTE: This is a network location used within PBL
###############################################################


# #### SCENMIP7 Lookup
# define.path <- "P:/globiom/Projects/SSPs/ScenarioMIP7/Results/lookup_table_results_withSSPV3/lookup_table_v5_16Jun2025/raw_results/Biodiversity_Link"
# G4M.path <- "P:/globiom/Projects/SSPs/ScenarioMIP7/Results/lookup_table_results_withSSPV3/lookup_table_v5_16Jun2025/raw_results/G4M/out/BetterLookup2SSP2a_20250613"
# proj.name <- "BetterLookup2SSP2a_20250613"
# all.scenarios <- read.csv(paste0(define.path, "/scenario_mapping_BetterLookup2SSP2a_20250613_selectedScens.csv"))
# scen.selection <- c("GHG000_BIO03", "GHG000_BIO06", "GHG100_BIO03", "GHG100_BIO06", "GHG400_BIO03", "GHG400_BIO06")
# choose.scen <- all.scenarios %>% filter(SCEN2%in%scen.selection)

#### ScenarioMIP7 - SSPV3
define.path <- "P:/globiom/Projects/SSPs/ScenarioMIP7/Results/scenario_results_with_SSPV3/2025v5/final_version/raw_results/Biodiversity_Link/Output/Trunk5266_MSGfdbkOct2025_final"
G4M.path <- "P:/globiom/Projects/SSPs/ScenarioMIP7/Results/scenario_results_with_SSPV3/2025v5/final_version/raw_results/G4M/out/Trunk5266_MSGfdbkOct2025_final_20251024"
proj.name <- "Trunk5266_MSGfdbkOct2025_final_20251024"
all.scenarios <- read.csv(paste0(define.path, "/scenario_mapping_Trunk5266_MSGfdbkOct2025_final_20251024.csv"))
scen.selection <- c("SSP2_VL" ,"SSP2_H", "SSP2_M", "SSP2_LOS", "LED_VL")
choose.scen <- all.scenarios %>% filter(SCEN3%in%scen.selection)

#### UK scaled
# define.path <- "P:/globiom/Projects/PBL_BIODIV_2025/UK_runs/Biodiversity_Link/Output"
# G4M.path <- "P:/globiom/Projects/PBL_BIODIV_2025/UK_runs/G4M/out/DESNZ_22082025"
# proj.name <- "DESNZ_22082025"
# all.scenarios <- read.csv(paste0(define.path, "/../scenario_mapping.csv"))
# scen.selection <- c("GHG000_LINEAR", "GHG100_LINEAR", "GHG400_LINEAR")
# choose.scen <- all.scenarios %>% filter(SCEN3=="SCENRCP4P5", SCEN2%in%scen.selection)



template.path <- "P:/globiom/Projects/PBL_BIODIV_2025/Postprocessing_ncdf/template"


mapping.G4M <- readRDS(file = 'input/G4M_mapping.RData')[[1]]
mapping.G4M <- apply(mapping.G4M, 2, as.character)
mapping.G4M <- data.frame(mapping.G4M)
mapping.G4M <- mapping.G4M %>% rename("g4m_id" = "g4m_05_id", "ns" = "SimUID")

###############################################################
# Select scenario to run 
###############################################################

curr.scen <- 1     # <-- loop index (single-run for now)
for(curr.scen in unique(choose.scen$ScenLoop)){

scen.setting <- NULL
curr.loops <- all.scenarios$ScenNr[all.scenarios$ScenLoop == curr.scen]

# Make scenario name + record scenario settings
sc_row <- unique(subset(all.scenarios, ScenLoop == curr.scen,
                        select = c("SCEN1", "SCEN2", "SCEN3")))
Scen <- paste(sc_row, collapse = "-")

curr.scen <- data.frame(Scen = Scen, sc_row, row.names = NULL)
scen.setting <- scen.setting %>% bind_rows(curr.scen)

G4M.path.scen <- paste0(G4M.path, "/area_harvest_map_",proj.name,"_",curr.scen$SCEN1,"_",curr.scen$SCEN3,"_",curr.scen$SCEN2,".csv")

###############################################################
# Read RData results to be processed
###############################################################

all.files <- list.files(define.path, full.names = TRUE)
filtered.files <- filter_paths(all.files, curr.loops)

if(amanda.files){
# Extract the 4-digit number before the last 6-digit number
four_digit <- sub(".*_(\\d{4})\\.\\d{6}\\.RData$", "\\1", filtered.files)
# Convert to numeric (invalid extractions become NA)
four_digit_num <- as.numeric(four_digit)
# Find the maximum 4-digit number
max_val <- max(four_digit_num, na.rm = TRUE)
# Filter files that match the maximum
filtered.files <- filtered.files[four_digit_num == max_val]
}

# # keep only files from the latest timestamp
# filtered.files <- filtered.files[
#   substr(filtered.files, nchar(filtered.files)-17, nchar(filtered.files)-13) ==
#     as.character(max(as.numeric(substr(filtered.files,
#                                        nchar(filtered.files)-17,
#                                        nchar(filtered.files)-13))))
# ]

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
results2 <- results2[
  , .(value = sum(value) * 1000),
  by = .(REGION, times, ns, lu.to, lu.from)
]

results2 <- as.data.frame(results2)

###### read and apply the forest is used G4M results to get closer to Biodiversity-Link res
if(correct.G4M){
  G4M.res <- read.csv(G4M.path.scen) %>% dplyr::select(g4m_id, year, used)
  
  
  ### with 2000 change. potentially wrong as lu.from is not mng for
  # results2 <- results2 %>% 
  #   left_join(mapping.G4M) %>% 
  #   left_join(G4M.res %>% 
  #       rename(times = year) %>%
  #       mutate(g4m_id = as.character(g4m_id))
  #   ) %>%
  #   mutate(used = ifelse(is.na(used), 0, used)) %>%
  #   mutate(
  #     mask = lu.to %in% c("priforest", "mngforest"),
  #     lu.to = if_else(mask, if_else(used == 1, "mngforest", "priforest"), lu.to)
  #   ) %>% dplyr::select(-mask) %>% 
  #   left_join(G4M.res %>% 
  #                 filter(year==2000) %>% 
  #                 mutate(times=2010, used_2000=used, g4m_id = as.character(g4m_id)) %>%
  #                 dplyr::select(-year,-used)
  #   ) %>%
  #   mutate(used_2000 = ifelse(is.na(used_2000), 0, used_2000)) %>%
  #   mutate(
  #     mask = lu.from %in% c("priforest", "mngforest"),
  #     lu.from = if_else(mask, if_else(used_2000 == 1, "mngforest", "priforest"), lu.from)
  #   ) %>% group_by(REGION, times, ns, lu.to, lu.from) %>%
  #   summarise(value=sum(value))
  
  
  results2 <- results2 %>% 
    left_join(mapping.G4M) %>% 
    left_join(G4M.res %>% 
                rename(times = year) %>%
                mutate(g4m_id = as.character(g4m_id))
    ) %>%
    mutate(used = ifelse(is.na(used), 0, used)) %>%
    mutate(
      mask = lu.to %in% c("priforest", "mngforest"),
      lu.to = if_else(mask, if_else(used == 1, "mngforest", "priforest"), lu.to)
    ) 
  
  
  setDT(results2)  
  results2 <- results2[
    , .(value = sum(value) * 1000),
    by = .(REGION, times, ns, lu.to, lu.from)
  ]
  results2 <- as.data.frame(results2)
  
}






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

linking_out <- results2

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
  rename(LC = lu.to, SimUID = ns, Year = times) 



linking.resultLC <- linking.resultLC %>%
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

linking.resultLC$Area <- rowSums(linking.resultLC[, -c(1:2)], na.rm = T)

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
  file = paste0("./output/", Scen, "_", date.tag, "_LULUC.csv"),
  row.names = FALSE
)

write.csv(
  scen.setting,
  file = paste0("./scen_setting_", date.tag, ".csv"),
  row.names = FALSE
)

###############################################################
# Optional: netCDF writing (external script)
###############################################################

if (write.nc) {
  source("codes/results2netcdf_halfdegree_MW.R")
}
}
###############################################################
# End of main script
###############################################################
