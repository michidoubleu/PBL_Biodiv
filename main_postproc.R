####### main script for PBL_Biodiv pipeline
library(dplyr)
library(tidyr)
##### read out files to be processed

define.path <- "P:/globiom/Projects/PBL_BIODIV_2025/SSP_runs/raw_results/Biodiversity_Link/output/"

#load scenario setup 
scenarios <- read.csv(paste0(define.path,"../../scenario_mapping_Trunk5264nutrient_MSGfdbk22Jul2025_20250724.csv"))

curr.scen <- 0
curr.loops <- scenarios$ScenNr[scenarios$ScenLoop==curr.scen]

all.files <- list.files(define.path, full.names = T)
ids <- sub(".*([0-9]{6})\\.RData$", "\\1", all.files)
filtered.files <- all.files[as.numeric(ids) %in% curr.loops]



df.link.results <- NULL
jjj <- 1
#for(jjj in 1:length(filtered.files)){
for(jjj in 1){
temp <- readRDS(filtered.files[jjj])
df.link.results <- full.res %>% bind_rows(temp[[5]])
}

########### (3)	Refer to the last part of the “get_G4M_DS_Link.R” in attachment to map the LC naming twice

full_simu_map = read.csv("./input/full_simu_map_biodiv.csv",stringsAsFactors = FALSE)

mapping_simuID <- full_simu_map %>% 
  dplyr::select(SimUID,country,REGION_37,colrowID) %>% 
  rename(COUNTRY=country) %>% 
  rename(REGION=REGION_37) %>% 
  mutate(SimUID=as.character(SimUID)) %>% 
  mutate(SimUID=as.numeric(as.character(SimUID))) 




## Compile LC from Link results-----
linking.resultLC.2000 <- df.link.results %>% subset(times %in% c("2010")) %>% 
  group_by(ns,times,lu.from) %>%
  summarize(value = sum(value),.groups = "keep") %>%
  # mutate(times=recode(times,"2010","2000")) %>% 
  mutate(times="2000") %>% 
  mutate(times=as.integer(times)) %>% 
  rename(LC=lu.from,SimUID=ns,Year=times)

linking.resultLC <- df.link.results %>% group_by(ns,times,lu.to) %>%
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

} # loop across 37 regions to




########### (4)	Continue with the “Prep_CSV_for_netcdf_2025again.R” to write out the dataframe (after the two LC name mappings) into CSV





########### (5)	Continue with the “results2netcdf_halfdegree_aug19_YZtestMOEJ_2025again.R” (shared in the last email) to convert the CSV into netCDF.





