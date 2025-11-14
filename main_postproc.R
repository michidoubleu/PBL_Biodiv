####### main script for PBL_Biodiv pipeline
library(dplyr)
library(tidyr)
library(data.table)

#### initial timings:
### approx 10 min to write the csv for 1 scenario
### approx XX min to write netcdf from 1 scenario (startetd at 15.20)


##### read out files to be processed

filter_paths <- function(paths, keep, digits = 6) {
  # extract the digits immediately before .RData (after last dot or underscore)
  ids <- sub(".*(?:\\.|_)(\\d+)\\.RData$", "\\1", paths, perl = TRUE)
  
  # if extraction failed for some paths, set to NA
  ids[grepl("^[0-9]+$", ids) == FALSE] <- NA
  
  # numeric value
  ids_num <- as.numeric(ids)
  
  # compute remainder for last `digits` digits
  mod <- 10^digits
  last_digits <- ids_num %% mod
  
  # logical index: keep when last_digits in keep (and not NA)
  idx <- !is.na(last_digits) & (last_digits %in% keep)
  
  paths[idx]
}
create.masked <- TRUE
create.global <- TRUE
write.nc <- FALSE


#define.path <- "P:/globiom/Projects/PBL_BIODIV_2025/SSP_runs/raw_results/Biodiversity_Link/output/"
#define.path <- "P:/globiom/Projects/PBL_BIODIV_2025/UK_runs/raw_results/Biodiversity_Link/output/"
define.path <- "P:/globiom/Projects/SSPs/ScenarioMIP7/Results/lookup_table_results_withSSPV3/lookup_table_v5_16Jun2025/raw_results/Biodiversity_Link"

#load scenario setup 
#scenarios <- read.csv(paste0(define.path,"../../scenario_mapping_Trunk5264nutrient_MSGfdbk22Jul2025_20250724.csv"))
scenarios <- read.csv(paste0(define.path,"/scenario_mapping_BetterLookup2SSP2a_20250613_2SCENS.csv"))
scen.setting <- NULL



curr.scen <- 3
##### start loop here





curr.loops <- scenarios$ScenNr[scenarios$ScenLoop==curr.scen]
# extract the unique scenario row
sc_row <- unique(subset(scenarios, ScenLoop == curr.scen, select = c("SCEN1", "SCEN2", "SCEN3")))
# construct full scenario name
Scen <- paste(sc_row, collapse = "-")

# build final data.frame
curr.scen <- data.frame(
  Scen = Scen,
  sc_row,
  row.names = NULL
)
scen.setting <- scen.setting %>% bind_rows(curr.scen)

all.files <- list.files(define.path, full.names = T)
filtered.files <- filter_paths(all.files, curr.loops)
filtered.files <- filtered.files[substr(filtered.files, nchar(filtered.files)-17, nchar(filtered.files)-13)==as.character(max(as.numeric(substr(filtered.files, nchar(filtered.files)-17, nchar(filtered.files)-13))))]

####### COUNTRY selection
## default, GLOBAL
## masked, topup

####### SCENARIO selection


df.link.results <- NULL
jjj <- 1
for(jjj in 1:length(filtered.files)){
#for(jjj in 1){
temp <- readRDS(filtered.files[jjj])
df.link.results <- df.link.results %>% bind_rows(temp[[5]])
}

########### (3)	Refer to the last part of the “get_G4M_DS_Link.R” in attachment to map the LC naming twice
#### updated by Michael to use data.table
mapping_LC_names_1 <- readRDS(file="./codes/mapping_for_G4MDSlink.RData")[[3]]
mapping_LC_names_2 <- readRDS(file="./codes/mapping_for_G4MDSlink.RData")[[4]]

# convert inputs
setDT(df.link.results)
mapping1 <- as.data.table(mapping_LC_names_1)
mapping2 <- as.data.table(mapping_LC_names_2)

# start from original table
results2 <- df.link.results

#### --- Stage 1: first LC mapping (mapping1) ---

# map lu.from
results2 <- mapping1[results2, on = .(lu.linkoutput = lu.from)]
results2[, lu.from := lu.new]
results2[, lu.new := NULL]

# map lu.to
results2 <- mapping1[results2, on = .(lu.linkoutput = lu.to)]
results2[, lu.to := lu.new]
results2[, lu.new := NULL]


#### --- Stage 2: protected area mapping (mapping2) ---

# map lu.from -> lu.final
results2 <- mapping2[results2, on = .(lu.new = lu.from)]
results2[, lu.from := lu.final]
results2[, lu.final := NULL]

# map lu.to -> lu.final
results2 <- mapping2[results2, on = .(lu.new = lu.to)]
results2[, lu.to := lu.final]
results2[, lu.final := NULL]


#### --- Stage 3: aggregate ---

results3 <- results2[
  , .(value = sum(value)),
  by = .(REGION, times, ns, lu.to, lu.from)
]

# optional: back to data.frame
results3 <- as.data.frame(results3)

########### (4)	Continue with the “Prep_CSV_for_netcdf_2025again.R” to write out the dataframe (after the two LC name mappings) into CSV
#source("codes/Prep_CSV_for_netcdf_2025again.R")

full_simu_map = read.csv("./template/full_simu_map_biodiv.csv",stringsAsFactors = FALSE)

mapping_simuID <- full_simu_map %>% 
  dplyr::select(SimUID,country,REGION_37,colrowID) %>% 
  rename(COUNTRY=country) %>% 
  rename(REGION=REGION_37) %>% 
  mutate(SimUID=as.character(SimUID)) %>% 
  mutate(SimUID=as.numeric(as.character(SimUID))) 




REGION_AG_Array <- c("ArgentinaReg","AustraliaReg","BrazilReg", "CanadaReg", "ChinaReg","CongoBasin", "EU_Baltic","EU_CentralEast", "EU_MidWest",  "EU_North","EU_South","Former_USSR",
                     "IndiaReg","IndonesiaReg","JapanReg","MalaysiaReg", "MexicoReg", "MiddleEast",
                     "NewZealandReg","NorthernAf", "Pacific_Islands", "RCAM","RCEU","ROWE",
                     "RSAM","RSAS","RSEA_OPA","RSEA_PAC","RussiaReg","SouthAfrReg",
                     "SouthKorea", "EasternAf", "SouthernAf",  "WesternAf",  "TurkeyReg", "UkraineReg",
                     "USAReg")


# linking_out <- readList2nd(path(str_glue( "../G4M_DS_Link/Output/G4M_DS_Link_result_",cluster_nr_G4MDSlink, ".", sprintf("%06d",i),".RData")))
linking_out <- results3

# linking_out_2 <- readList2nd(path(str_glue( "../G4M_DS_Link/Output/G4M_DS_Link_result_",cluster_nr_G4MDSlink, ".", sprintf("%06d",2),".RData")))

## Compile LC from Link results-----
linking.resultLC.2000 <- linking_out %>% subset(times %in% c("2010")) %>% 
  group_by(ns,times,lu.from) %>%
  summarize(value = sum(value),.groups = "keep") %>%
  # mutate(times=recode(times,"2010","2000")) %>% 
  mutate(times="2000") %>% 
  mutate(times=as.integer(times)) %>% 
  rename(LC=lu.from,SimUID=ns,Year=times)

linking.resultLC <- linking_out %>% group_by(ns,times,lu.to) %>%
  summarize(value = sum(value),.groups = "keep") %>%
  rename(LC=lu.to,SimUID=ns,Year=times) %>% 
  bind_rows(linking.resultLC.2000) %>% 
  spread(key=LC,value = value) %>%
  mutate(SimUID=as.numeric(as.character(SimUID))) %>% 
  mutate(urban=0) %>% 
  mutate(priforest=ifelse(is.na(priforest),0,priforest)) %>% 
  mutate(mngforest=ifelse(is.na(mngforest),0,mngforest)) %>% 
  mutate(cropland=ifelse(is.na(cropland),0,cropland)) %>% 
  mutate(grassland=ifelse(is.na(grassland),0,grassland)) %>% 
  mutate(other=ifelse(is.na(other),0,other)) %>% 
  # mutate(SRP=ifelse(is.na(SRP),0,SRP)) %>% 
  # mutate(mutate(across(any_of(SRP), ifelse(is.na(SRP),0,SRP)))) %>% 
  mutate(protected_other=0) %>% 
  mutate(protected_priforest=0) 

  linking.resultLC <- linking.resultLC %>% 
    mutate(restored=0) 

if("SRP" %in% names(linking.resultLC)){
  linking.resultLC <- linking.resultLC %>% mutate(SRP=ifelse(is.na(SRP),0,SRP))
}else{
  linking.resultLC$SRP=0 }

## Calculate SimU Area------
linking.resultLC$Area=rowSums(linking.resultLC[,-c(1:2)])

## Add colrow, REGION info, and compute final------
linking.resultLC.final <- linking.resultLC %>% 
  left_join(mapping_simuID %>% dplyr::select(-REGION),by=c("SimUID")) %>% 
  rename(Colrow=colrowID) %>% 
  dplyr::select(c(SimUID,Area,Year,COUNTRY,Colrow,cropland,grassland,priforest,mngforest,SRP,restored,other,protected_priforest,protected_other,urban)) %>% 
  arrange(Year,SimUID)



# linking_out <- NULL
# linking.resultLC.final <- NULL

Array_fullLUC_list <- c("cropland.cropland","cropland.grassland","cropland.mngforest","cropland.SRP","cropland.restored","cropland.other","grassland.grassland","grassland.cropland","grassland.mngforest","grassland.SRP","grassland.restored","grassland.other","priforest.priforest","priforest.cropland","priforest.grassland","priforest.mngforest","mngforest.mngforest","SRP.SRP","SRP.mngforest",  'SRP.other',"restored.restored",	"other.other",	'other.cropland',	"other.grassland",	"other.mngforest",	'other.SRP',	'protected_priforest.protected_priforest',	'protected_other.protected_other',	"urban.urban")



# 2.2 Extracting LUC results ------

ds.resultLUC0 <- linking_out %>%
  mutate(LUC=paste0(lu.from,".",lu.to)) %>% 
  dplyr::select(-c(lu.from,lu.to)) %>%
  # group_by(REGION,times,ns,LUC) %>% mutate(row_number() == 1) %>% 
  # subset(`row_number() == 1`==TRUE) %>% 
  dplyr::select(-c(REGION)) %>%
  left_join(mapping_simuID[,c("SimUID","REGION")]%>%rename(ns=SimUID) %>% mutate(ns=as.character(ns)),by=c("ns")) %>%  # YW20240305: Addressing duplicated cells in REGION37
  spread(key=LUC,value = value) %>%
  mutate(urban.urban=0) %>% 
  ungroup()

for(k in 1:length(Array_fullLUC_list)){
  if(!(Array_fullLUC_list[k] %in% colnames(ds.resultLUC0))){
    ds.resultLUC0$NewCol <- 0
    names(ds.resultLUC0)[names(ds.resultLUC0)=="NewCol"] <- Array_fullLUC_list[k]
  }
}

ds.resultLUC1 <- ds.resultLUC0 %>% 
  rename(SimUID=ns, Year=times) %>% 
  mutate(SimUID=as.numeric(as.character(SimUID))) %>% 
  dplyr::select(SimUID,Year,all_of(Array_fullLUC_list))


linking.resultLCLUC.final <- linking.resultLC.final %>% 
  left_join(ds.resultLUC1,by=c("SimUID","Year"))

linking.resultLCLUC.final[is.na(linking.resultLCLUC.final)] <- 0



linking.resultLC.final.all <- linking.resultLC.final %>% 
  arrange(Year,SimUID)
linking.resultLCLUC.final.all <- linking.resultLCLUC.final %>% 
  arrange(Year,SimUID)



#====================Output results================================
# SCEN3 <- "CLIM_CONS"
# write.csv(linking.resultLC.final.all, 
#           file = paste0("./output/",Scen,"_",Sys.Date(),"_LC.csv"),
#           row.names = FALSE)

write.csv(linking.resultLCLUC.final.all, 
          file = paste0("./output/",Scen,"_",Sys.Date(),"_LULUC.csv"),
          row.names = FALSE)

### end loop here

write.csv(scen.setting, 
          file = paste0("./scen_setting_",Sys.Date(),".csv"),
          row.names = FALSE)


########### (5)	Continue with the “results2netcdf_halfdegree_aug19_YZtestMOEJ_2025again.R” (shared in the last email) to convert the CSV into netCDF.

if(write.nc){
source("codes/results2netcdf_halfdegree_MW.R")
}




