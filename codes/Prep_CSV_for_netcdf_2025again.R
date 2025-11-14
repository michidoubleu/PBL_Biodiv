# Prep_CSV_for_netcdf
## Transforming the MOEJ scenario results into a csv in the same format as \downscale_biodiversity_G4M\resJan19_SPA0_RCPref_SSP2_NOBIOD_SSP2_lowp2cor3dev3_6jun19_mod1__luc.csv
## based on Michi's script "G4M_DS_to_link_final_limpopo.R", as here we also need to do similar: to merge GLOBIOM and G4M results
## Input: Downscaled GLOBIOM results; G4M results
## NB: in the future, it is worth to double check with Tamas's script of csv writing and merging
## By: Yazhen, 2023/09/15

rm(list=ls())

# Library -----------
library(gdxrrw)
suppressWarnings(library(tidyverse))
library(stringr) # part of the tidyverse, but not attached by default
library(fs)

# Define function ----------
# loading function
rfunc <- function(x) readRDS(x)[[3]]$out.res
readList1st <- function(x) readRDS(x)[[1]]
readList2nd <- function(x) readRDS(x)[[2]]
readList3rd <- function(x) readRDS(x)[[3]]

CD <- "e:\\IIASA\\Model\\GLOBIOM-G4M-link\\"
# CD <- "f:\\Model\\GLOBIOM-G4M-link-REG37(OldTrunk4318)\\"
# WD_DOWNSCALING <- "e:\\IIASA\\Model\\GLOBIOM-G4M-link\\DownScale\\"
WD_DOWNSCALING <- "DownScale"
setwd(CD)
setwd(WD_DOWNSCALING)

# source("../Biodiversity_Link/G4M_DS_to_simu_link_final_limpopo.R")
source("G4M_DS_link_function_YZtest.R")

# config <- readRDS(path("../Biodiversity_Link/Input","config.RData"))
# scen_map <- config[[3]]


# Define task -----

# DEFINE clusterNr, functions
# project <- "ScenTest5"
# lab <- "27082023"
# cluster_nr_downscaling <- "1940"
# cluster_nr_G4MDSlink <- "2197"
# cluster_nr_G4MDSlink <- "2247

# project <- "ScenTest10FP"
# lab <- "04112023"
# cluster_nr_globiom <- "3107"
# cluster_nr_downscaling <- "3114"
# cluster_nr_G4MDSlink <- "3247" #ScenTestFP10

# project <- "OldTrunkDefault"
# lab <- "13082024"
# cluster_nr_globiom <- "6613"
# cluster_nr_downscaling <- "6616"
# cluster_nr_G4MDSlink <- "6719"

project <- "ScenDriver5new"
lab <- "13012025"
#ScenDriver5new
cluster_nr_globiom <- "7147" 
cluster_nr_downscaling <- "7158"
cluster_nr_G4MDSlink <- "7180"


PROJECT <- project
DATE_LABEL <- lab
source(paste0(CD,"\\R\\run_submission.R"))
source(paste0(CD,"\\R\\helper_functions.R"))
source(paste0(CD,"\\R/configuration/custom_",project,".R"), local=TRUE, echo=FALSE)
scenario_mapping <- get_mapping()


# Load supporting data ------
# Supprting information
full_simu_map = read.csv("full_simu_map_biodiv.csv",stringsAsFactors = FALSE)

mapping_simuID <- full_simu_map %>% 
  dplyr::select(SimUID,country,REGION_37,colrowID) %>% 
  rename(COUNTRY=country) %>% 
  rename(REGION=REGION_37) %>% 
  mutate(SimUID=as.character(SimUID)) %>% 
  mutate(SimUID=as.numeric(as.character(SimUID))) 


# Define scenarios ------
# SCEN3_ARRAY <- c("BASE","CLIM", "CLIM_CONS_SDGL")
# SCEN3_ARRAY <- c("BASE","CLIM","CLIM_CONS","CLIM_CONS_SDGL")
SCEN3_ARRAY <- list(
  c("SSP2","SCEN","BASE"),
  c("SSP2","SCEN","CLIM"),
  c("SSP2","SCEN","CONS"),
  c("SSP2","SCEN","CLIM_CONS"),
  c("SSP2","SCEN","CLIM_CONS_SDGL"))

# SCEN3_ARRAY <- list(
#   c("SSP1","SCEN","SCENRCPREF"),
#   c("SSP1","SCEN","SCENRCP6P0"),
#   c("SSP1","SCEN","SCENRCP4P5"),
#   c("SSP1","SCEN","SCENRCP3P7"),
#   c("SSP1","SCEN","SCENRCP2P6"),
#   c("SSP1","SCEN","SCENRCP2P0"),
#   c("SSP1","SCEN","SCENRCP1P9"),
#   c("SSP1","SCEN","SCENRCP1P8"),
#   c("SSP2","SCEN","SCENRCPREF"),
#   c("SSP2","SCEN","SCENRCP6P0"),
#   c("SSP2","SCEN","SCENRCP4P5"),
#   c("SSP2","SCEN","SCENRCP3P7"),
#   c("SSP2","SCEN","SCENRCP2P6"),
#   c("SSP2","SCEN","SCENRCP2P0"),
#   c("SSP2","SCEN","SCENRCP1P9"),
#   c("SSP2","SCEN","SCENRCP1P8"),
#   c("SSP3","SCEN","SCENRCPREF"),
#   c("SSP3","SCEN","SCENRCP6P0"),
#   c("SSP3","SCEN","SCENRCP4P5"),
#   c("SSP3","SCEN","SCENRCP3P7"),
#   c("SSP3","SCEN","SCENRCP2P6"),
#   c("SSP4","SCEN","SCENRCPREF"),
#   c("SSP4","SCEN","SCENRCP6P0"),
#   c("SSP4","SCEN","SCENRCP4P5"),
#   c("SSP4","SCEN","SCENRCP3P7"),
#   c("SSP4","SCEN","SCENRCP2P6"),
#   c("SSP4","SCEN","SCENRCP2P0"),
#   c("SSP4","SCEN","SCENRCP1P9"),
#   c("SSP4","SCEN","SCENRCP1P8"),
#   c("SSP5","SCEN","SCENRCPREF"),
#   c("SSP5","SCEN","SCENRCP6P0"),
#   c("SSP5","SCEN","SCENRCP4P5"),
#   c("SSP5","SCEN","SCENRCP3P7"),
#   c("SSP5","SCEN","SCENRCP2P6")
#   )

for(sc3 in 1:length(SCEN3_ARRAY) ){
# SCEN1 <- "SSP2"
# SCEN2 <- "SCEN"
# SCEN3 <- SCEN3_ARRAY[sc3]
SCEN1 <- SCEN3_ARRAY[[sc3]][1]
SCEN2 <- SCEN3_ARRAY[[sc3]][2]
SCEN3 <- SCEN3_ARRAY[[sc3]][3]
# SCEN3 <- "BASE"
# SCEN3 <- "CLIM_CONS"
# SCEN3 <- "CLIM_CONS_SDGL"
# "BASE","CLIM", "CONS","CLIM_CONS","CLIM_CONS_SDGL"

ScenNr_currentScen <- scenario_mapping$ScenNr[scenario_mapping$SCEN1==SCEN1&scenario_mapping$SCEN2==SCEN2&scenario_mapping$SCEN3==SCEN3]
i_start <- min(ScenNr_currentScen)
i_end <- max(ScenNr_currentScen)


#====================2. New method: first do the G4M_DS_link, then merge=========================

cat(paste0("Start processing scenario: ", SCEN1,"-",SCEN2,"-",SCEN3,"\n"))
cat("Read downscaling results: \n")

# i=0

linking.resultLCLUC.final.all <- NULL
linking.resultLC.final.all <- NULL

REGION_AG_Array <- c("ArgentinaReg","AustraliaReg","BrazilReg", "CanadaReg", "ChinaReg","CongoBasin", "EU_Baltic","EU_CentralEast", "EU_MidWest",  "EU_North","EU_South","Former_USSR",
                     "IndiaReg","IndonesiaReg","JapanReg","MalaysiaReg", "MexicoReg", "MiddleEast",
                     "NewZealandReg","NorthernAf", "Pacific_Islands", "RCAM","RCEU","ROWE",
                     "RSAM","RSAS","RSEA_OPA","RSEA_PAC","RussiaReg","SouthAfrReg",
                     "SouthKorea", "EasternAf", "SouthernAf",  "WesternAf",  "TurkeyReg", "UkraineReg",
                     "USAReg")
i_count=1

# for(i in 0:36){ # CONS
# for(i in 37:73){ # BASE
  # for(i in 74:110){ # CLIM
  # for(i in 111:147){ # CLIM_CONS
  # for(i in 185:221){ # CLIM_CONS_SDG
  for(i in i_start:i_end){ # now flexible, for any scenario
  
  ## 2.1 Extract LC result ------
  cat(i,"- ")
  
  # current_file <- path(str_glue( "../G4M_DS_Link/Output/G4M_DS_Link_result_",cluster_nr_G4MDSlink, ".", sprintf("%06d",i),".RData"))
  # current_file <- path(str_glue( "../G4M_DS_Link/Output/G4M_DS_Link_results_3247/G4M_DS_Link_result_",cluster_nr_G4MDSlink, ".", sprintf("%06d",i),".RData"))
  current_file <- path(str_glue( "../G4M_DS_Link/Output/G4M_DS_Link_result_",cluster_nr_G4MDSlink, ".", sprintf("%06d",i),".RData"))
  
  if(file.exists(current_file)){
  
  # linking_out <- readList2nd(path(str_glue( "../G4M_DS_Link/Output/G4M_DS_Link_result_",cluster_nr_G4MDSlink, ".", sprintf("%06d",i),".RData")))
  linking_out <- readList3rd(current_file)
  
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
  
  if(grepl(pattern = "CONS",x = SCEN3)){
    linking.resultLC <- linking.resultLC %>% 
      mutate(restored=ifelse(is.na(restored),0,restored))
  }else{
    linking.resultLC <- linking.resultLC %>% 
      mutate(restored=0) 
    }
  ##Yazhen: temp code for restored and protected lands, need to be adapted in the next
  
  if("SRP" %in% names(linking.resultLC)){
    linking.resultLC <- linking.resultLC %>% mutate(SRP=ifelse(is.na(SRP),0,SRP))
  }else{
    linking.resultLC$SRP=0 }
  
  ## Calculate SimU Area------
  linking.resultLC$Area=rowSums(linking.resultLC[,-c(1:2)])
  
  ## Add colrow, REGION info, and compute final------
  linking.resultLC.final <- linking.resultLC %>% 
    left_join(mapping_simuID,by=c("SimUID")) %>% 
    rename(Colrow=colrowID) %>% 
    dplyr::select(c(SimUID,Area,Year,REGION,COUNTRY,Colrow,cropland,grassland,priforest,mngforest,SRP,restored,other,protected_priforest,protected_other,urban)) %>% 
    arrange(Year,SimUID) %>% 
    filter((REGION %in% REGION_AG_Array[i_count])) ## YW 20240305: Addressing duplicated cells in REGION37
  
  
  ## Merge the current region into the full matrix------
  # if(is.null(linking.resultLC.final.all)){
  if(i==i_start){
    linking.resultLC.final.all <- linking.resultLC.final
  }else{
    linking.resultLC.final.all <- rbind(linking.resultLC.final.all,linking.resultLC.final)
  }
  
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
    subset((REGION %in% REGION_AG_Array[i_count])) %>%  ## YW 20240305: Addressing duplicated cells in REGION37
    dplyr::select(SimUID,Year,all_of(Array_fullLUC_list))

  
  linking.resultLCLUC.final <- linking.resultLC.final %>% 
    left_join(ds.resultLUC1,by=c("SimUID","Year"))
  
  linking.resultLCLUC.final[is.na(linking.resultLCLUC.final)] <- 0
  

  ## Merge the current region into the full matrix------
  # if(i==0){
  if(i==i_start){
  # if(is.null(linking.resultLCLUC.final.all)){
    linking.resultLCLUC.final.all <- linking.resultLCLUC.final
  }else{
    linking.resultLCLUC.final.all <- rbind(linking.resultLCLUC.final.all,linking.resultLCLUC.final)
  } ## end merging results

  } # if current_file exist
  
  i_count = i_count + 1
  
} # loop across 37 regions to merge the linking results
  
cat("\n")

linking.resultLC.final.all <- linking.resultLC.final.all %>% 
  arrange(Year,SimUID)
linking.resultLCLUC.final.all <- linking.resultLCLUC.final.all %>% 
  arrange(Year,SimUID)



#====================Output results================================
# SCEN3 <- "CLIM_CONS"
write.csv(linking.resultLC.final.all, 
          # file = paste0("e:\\IIASA\\Model\\downscale_biodiversity_G4M\\Yazhen_test\\csv_G4M_DS_link_output\\",project,"_",lab,"_",SCEN1,"_",SCEN2,"_",SCEN3,"_2024Mar06.csv"),
          file = paste0("e:\\IIASA\\Model\\downscale_biodiversity_G4M\\Yazhen_test\\csv_G4M_DS_link_output\\",project,"_",lab,"_",SCEN1,"_",SCEN2,"_",SCEN3,"_2025Sept01.csv"),
          row.names = FALSE)

write.csv(linking.resultLCLUC.final.all,
          # file = paste0("e:\\IIASA\\Model\\downscale_biodiversity_G4M\\Yazhen_test\\csv_G4M_DS_link_output\\",project,"_",lab,"_",SCEN1,"_",SCEN2,"_",SCEN3,"_test0922.csv"),
          # file = paste0("e:\\IIASA\\Model\\downscale_biodiversity_G4M\\Yazhen_test\\csv_G4M_DS_link_output\\",project,"_",lab,"_",SCEN1,"_",SCEN2,"_",SCEN3,"_LULUC_2024Mar06.csv"),
          file = paste0("e:\\IIASA\\Model\\downscale_biodiversity_G4M\\Yazhen_test\\csv_G4M_DS_link_output\\",project,"_",lab,"_",SCEN1,"_",SCEN2,"_",SCEN3,"_LULUC_2025Sept01.csv"),
          row.names = FALSE)

} #loop for sc3(SCEN3)

#====================Checking================================
if(FALSE){
#1. Check DS result: whether DS was successfully performed
ds.resultLC.World <- DS.res.G4M %>% group_by(times,lu.to) %>%
  summarize(downscale.value = sum(value),.groups = "keep") %>%
  spread(key=times,value = downscale.value) %>%
  rename(LC=lu.to)
}

#====================Draft Area================================
if(FALSE){
downscalr_out1 <- downscalr_out %>%
  # mutate(lu.to=recode(lu.to,"forest_new_ha"="managed_forest","forest_old_ha"="unmanaged_forest"))
  mutate(lu.to=recode(lu.to,"Forest"="managed_forest"))

downscalr_out1 <- downscalr_out %>%
  mutate(lu.from=recode("forest_new_ha"="managed_forest","forest_old_ha"="unmanaged_forest"))
}




if(FALSE){
  ds.resultLUC <- ds.resultLUC %>% 
    select(cropland.cropland,cropland.grassland,cropland.mngforest,cropland.SRP,cropland.restored,cropland.other,grassland.grassland,grassland.cropland,grassland.mngforest,grassland.SRP,grassland.restored,grassland.other,priforest.priforest,priforest.cropland,priforest.grassland,priforest.mngforest,mngforest.mngforest,SRP.SRP,SRP.mngforest,
           SRP.other,restored.restored,	other.other,	other.cropland,	other.grassland,	other.mngforest,	other.SRP,	protected_priforest.protected_priforest,	protected_other.protected_other,	urban.urban
    )
}


#===================BACKUP=================
#====================1.Old method: first merge DS result================================

## Load DS result ------
# i=0
# for(i in 0:(nr.output-1)){
cat("Read downscaling results: ")

downscalr_out_all <- NULL

# for(i in 0:36){ # CONS
for(i in 37:73){ # BASE
  # for(i in 74:110){ # CLIM
  # for(i in 111:147){ # CLIM_CONS
  # for(i in 185:221){ # CLIM_CONS
  # file_DS <- paste0("gdx/output_",cluster_nr_downscaling,".",as.character(formatC(i, width = 6, format = "d", flag = "0")),".RData")
  # downscalr_out <- rfunc(file_DS)
  cat(i,"- ")
  downscalr_out <- rfunc(path(str_glue( "gdx/output_",cluster_nr_downscaling, ".",
                                        sprintf("%06d",i),".RData")))
  
  if(is.null(downscalr_out_all)){
    downscalr_out_all <- downscalr_out
  }else{
    downscalr_out_all <- rbind(downscalr_out_all,downscalr_out)
  }
  
  downscalr_out <- NULL
  
}

# downscalr_out_all.0 <- downscalr_out_all


## Do the G4M_DS_link processing------
# Method 1: call the slightly modified g4mid_to_simuid_YZtest function and send the full-region downscaled results to process, and get the merged (GLOBIOM and G4M) results 
df.link.results <- g4mid_to_simuid_YZtest(downscalr_out_all %>% mutate(value=value*1000),lab,project,SCEN1,SCEN2,SCEN3)
write.csv(df.link.results,file="results_GLOBIOM_G4M_Link_ScenTest5_CONS.csv",row.names = FALSE)

# Method 2: manually run the g4mid_to_simuid_YZtest function line-by-line
DS.res.G4M = downscalr_out_all %>% mutate(value=value*1000)
lab=lab;project=project;curr.SCEN1=SCEN1;curr.SCEN2=SCEN2;curr.SCEN3=SCEN3
# (return DS.res.G4M)
DS.res.G4M.0 <- DS.res.G4M


## Recode LC_Type names (for managed forest)------
## Method1
DS.res.G4M <- df.link.results %>% 
  mutate(lu.from=recode(lu.from,"Forest"="priforest","CrpLnd"="cropland","Grass"="grassland","PltFor"="SRP","RstLnd"="restored","OthNatLnd"="other")) %>% 
  mutate(lu.to=recode(lu.to,"Forest"="priforest","CrpLnd"="cropland","Grass"="grassland","PltFor"="SRP","RstLnd"="restored","OthNatLnd"="other")) 

## Method2
# DS.res.G4M <- DS.res.G4M %>%
# mutate(lu.from=recode("forest_new_ha"="managed_forest","forest_old_ha"="priforest","CrpLnd"="cropland","GrsLnd"="grassland","PltFor"="SRP","RstLnd"="restored","OthNatLnd"="other")) %>% 
# mutate(lu.to=recode("forest_new_ha"="managed_forest","forest_old_ha"="priforest")) 
DS.res.G4M <- downscalr_out_all %>%
  mutate(lu.from=recode(lu.from,"Forest"="priforest","CrpLnd"="cropland","Grass"="grassland","PltFor"="SRP","RstLnd"="restored","OthNatLnd"="other")) %>% 
  mutate(lu.to=recode(lu.to,"Forest"="priforest","CrpLnd"="cropland","Grass"="grassland","PltFor"="SRP","RstLnd"="restored","OthNatLnd"="other")) 


## Compile LC from Link results-----
ds.resultLC <- DS.res.G4M %>% group_by(ns,times,lu.to) %>%
  summarize(value = sum(value),.groups = "keep") %>%
  rename(LC=lu.to) %>% 
  spread(key=LC,value = value) %>%
  mutate(urban=0)

## Calculate SimU Area------
ds.resultLC <- ds.resultLC %>% 
  rename(grassland=Grass,SimUID=ns)
ds.resultLC$Area=rowSums(ds.resultLC[,-c(1:2)])


# Mapping: add colrows and region information ----
na_sum = function(x) {return(sum(x[!is.na(x)]))}
loncr = function(x) {89.75 - (x-1) * .5}
latcr = function(x) {-179.75 + (x-1) * .5}
loncr_inv = function(x) { (x + 89.75) / .5 + 1}
latcr_inv = function(x) { (x + 179.75) / .5 + 1}

full_simu_map = read.csv("full_simu_map_biodiv.csv",stringsAsFactors = FALSE)

mapping_simuID <- full_simu_map %>% 
  select(SimUID,country,REGION_37,colrowID) %>% 
  rename(COUNTRY=country) %>% 
  rename(REGION=REGION_37) %>% 
  mutate(SimUID=as.character(SimUID))

ds.resultLC.final <- ds.resultLC %>% 
  left_join(mapping_simuID,by=c("SimUID")) %>% 
  rename(Year=times) %>% rename(Colrow=colrowID) %>% 
  # select(c(SimUID,Area,Year,REGION,COUNTRY,Colrow,cropland,priforest,mngforest,SRP,restored,other,protected_priforest,protected_other,urban
  select(c(SimUID,Area,Year,REGION,COUNTRY,Colrow,cropland,priforest,SRP,restored,other,protected_priforest,protected_other,urban)) %>% mutate(mngforest=0) %>% relocate(mngforest,.after=priforest)



### Compile LUC from Link results
# ds.resultLUC <- DS.res.G4M %>% group_by(ns,times,lu.from,lu.to) %>%
# ds.resultLUC <- df.link.results %>% group_by(ns,times,lu.from,lu.to) %>%
# summarize(value = sum(value),.groups = "keep") %>%
ds.resultLUC <- df.link.results %>%
  mutate(LUC=paste0(lu.from,".",lu.to)) %>% 
  select(-c(lu.from,lu.to)) %>%
  spread(key=LUC,value = value) %>%
  mutate(urban.urban=0) 

ds.resultLUC <- ds.resultLUC %>% 
  select(cropland.cropland,cropland.grassland,cropland.mngforest,cropland.SRP,cropland.restored,cropland.other,grassland.grassland,grassland.cropland,grassland.mngforest,grassland.SRP,grassland.restored,grassland.other,priforest.priforest,priforest.cropland,priforest.grassland,priforest.mngforest,mngforest.mngforest,SRP.SRP,SRP.mngforest,
         SRP.other,restored.restored,	other.other,	other.cropland,	other.grassland,	other.mngforest,	other.SRP,	protected_priforest.protected_priforest,	protected_other.protected_other,	urban.urban
  )

ds.resultLC.final <- ds.resultLC.final %>% 
  left_join(ds.resultLUC,by=c("SimUID"))
