args <- commandArgs(trailing = TRUE)


library(tidyverse)
library(fs)


# Get scenario------
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} else {
  # default output file
  scen <- as.integer(args[1])
}


# Load G4M output conversion function------
  # source("G4M_DS_to_simu_link_final_limpopo.R")
  source("G4M_DS_to_simu_link_YW1.R")


# Define cluster number of downscaling post-processing------
config <- readRDS(path("Input","config.RData"))

project <- config[[1]]
lab <- config[[2]]
scen_map <- config[[3]]
cluster_nr <- config[[4]]
get_bii <- config[[5]]
get_csar <- config[[6]]


# Define loading function------
rfunc <- function(x) readRDS(x)[[3]]$out.res
# rfunc_newRstLndrep <- function(x) readRDS(x)[[7]]


# Prep Mapping data ------
cat("==> Read mapping data","\n")

# mapping: ns to g4m 0.5degree grids
mapping <- readRDS(file='G4m_mapping.RData')[[1]]
mapping <- data.frame(apply(mapping, 2, as.numeric))

# mapping between: x-y-g4m_id (mapping lat-lon data to g4m 0.5degree grids)
# mapping_g4mid_xy_new <- read.csv(file=path("oecd_ssp2_decc2016_pbl_r31_16032022_0.csv")) %>%
  # dplyr::select(x,y,simuid) %>% rename("g4m_05_id"="simuid")
mapping_g4mid_xy_new <- readRDS(file="mapping_for_G4MDSlink.RData")[[1]]


# maplayer_isforest_new: isForest or not at g4m 0.5degree grids
# maplayer_isforest0 <- read.csv(file="forest_nonforest_LUH2_halfdeg.csv")
# maplayer_isforest0 <- read.csv(file=path("forest_nonforest_LUH2_halfdeg.csv"))
# maplayer_isforest_new <- maplayer_isforest0 %>%
#   rename("IsForest"=colnames(maplayer_isforest0)[3]) %>%
#   left_join(mapping_g4mid_xy_new,by=c("x","y")) %>% na.omit() %>%
#   dplyr::select(g4m_05_id,IsForest) %>%
#   mutate(IsForest=as.numeric(IsForest)) %>%
#   mutate(g4m_05_id=as.character(g4m_05_id))
maplayer_isforest_new <- readRDS(file="mapping_for_G4MDSlink.RData")[[2]]


# cat(paste0("Checking RstLnd for: ",scen1,"_",scen2,"_",scen3,", ScenNr=",scen),"\n")


# Define scenarios------
cat("==> Define scenarios","\n")
current_scen <- scen_map %>% filter(ScenNr == scen)

scen1 <- current_scen$SCEN1
scen2 <- current_scen$SCEN2
scen3 <- current_scen$SCEN3

# Load downscalr output------
cat("==> Load downscalr output","\n")
# in local PC it takes about 10 sec
# downscalr_out <- rfunc(path(str_glue( "output_", cluster_nr, ".",
                                      # sprintf("%06d",current_scen$ScenNr),".RData")))

downscalr_out0 <- rfunc(path(str_glue( "output_", cluster_nr, ".",
                                      sprintf("%06d",current_scen$ScenNr),".RData")))

#  group_by(times) %>% group_split()


# Process RstLnd in res------
cat("==> Process RstLnd in res","\n")
res0 <- downscalr_out0
# res0$ns <- as.numeric(res0$ns)
mapping$SimUID <- as.character(mapping$SimUID)
mapping$g4m_05_id <- as.character(mapping$g4m_05_id)

# 1) mapping dat to g4mid
res1 <- res0 %>% left_join(mapping, by=c("ns"="SimUID"))

# Assign RstLnd to afforestable or non-afforestable
res2 <- res1 %>% left_join(maplayer_isforest_new) %>%
  mutate(IsForest=ifelse(is.na(IsForest),0,IsForest)) %>%
  mutate(lu.from=ifelse(IsForest==1,
                        recode(lu.from,"RstLnd"="RstLnd_YesAffor"),
                        recode(lu.from,"RstLnd"="RstLnd_NoAffor")))%>%
  mutate(lu.to=ifelse(IsForest==1,
                      recode(lu.to,"RstLnd"="RstLnd_YesAffor"),
                      recode(lu.to,"RstLnd"="RstLnd_NoAffor")))

# Reallocate Rstlnd_YesAffor to OtherNatLnd ------

res3 <- res2 %>%
  mutate(lu.from=recode(lu.from,"RstLnd_YesAffor"="OthNatLnd")) %>%
  mutate(lu.to=recode(lu.to,"RstLnd_YesAffor"="OthNatLnd"))

res4 <- res3 %>%
  # group_by(REGION,times,ns,lu.to,value,lu.from,g4m_05_id) %>% summarise(value=sum(value)) %>%
  group_by(REGION,times,ns,lu.to,lu.from,g4m_05_id) %>% summarise(value=sum(value)) %>%
  select(-c(g4m_05_id)) %>% ungroup()

downscalr_out <- res4


## Yazhen: temp code for adding region
if(!("REGION" %in% colnames(downscalr_out))){
  REGION_AG_Array <- c("ArgentinaReg","AustraliaReg","BrazilReg", "CanadaReg", "ChinaReg","CongoBasin", "EU_Baltic","EU_CentralEast", "EU_MidWest",  "EU_North","EU_South","Former_USSR",
                       "IndiaReg","IndonesiaReg","JapanReg","MalaysiaReg", "MexicoReg", "MiddleEast",
                       "NewZealandReg","NorthernAf", "Pacific_Islands", "RCAM","RCEU","ROWE",
                       "RSAM","RSAS","RSEA_OPA","RSEA_PAC","RussiaReg","SouthAfrReg",
                       "SouthKorea", "EasternAf", "SouthernAf",  "WesternAf",  "TurkeyReg", "UkraineReg",
                       "USAReg")

  downscalr_out <- bind_cols(REGION=REGION_AG_Array[current_scen$ScenNr+1],downscalr_out)

}


# Get output------
cat("==> Get G4M_DS_Link output","\n")

# YW note: results can be used for get_biodiversity.R for computing BIs; result2/3 can be used for Prep_csv_for_netcdf.R for netCDF output

# results <- g4mid_to_simuid(downscalr_out,lab,project,scen1,scen2,scen3)
resultsAll <- g4mid_to_simuid_YW1(downscalr_out %>% mutate(value=value*1000),lab,project,scen1,scen2,scen3) #YW: adjust DS output values (originally in 1000ha, x1000 to be ha) to be consistent with G4M results (in ha)
      # Time track: starting 1:55pm, ending 3:22
results1 <- resultsAll[[1]]

#results <- downscalr_out %>% map_df(~g4mid_to_simuid(.,lab,project,scen1,scen2,scen3)) %>% rbind


# Rename LCs------
cat("==> Rename LCs (to get results2,3)","\n")


## Rename RstLnd and other LCs to be in line with output matrix for netCDF
### Old method: using "recode" function (a bit slow)

# results2 <- results %>%
  # mutate(lu.from=recode(
    # lu.from,"forest_old_ha"="priforest","forest_new_ha"="mngforest","CrpLnd"="cropland","Grass"="grassland","PltFor"="SRP","RstLnd"="restored","RstLnd_NoAffor"="restored","OthNatLnd"="other")) %>%
  # mutate(lu.to=recode(
    # lu.to,"forest_old_ha"="priforest","forest_new_ha"="mngforest","CrpLnd"="cropland","Grass"="grassland","PltFor"="SRP","RstLnd"="restored","RstLnd_NoAffor"="restored","OthNatLnd"="other"))

### New method: using mapping and left_join() (is potentially a bit faster)
# mapping_LC_names_1 <- read.csv(file=path("mapping_LC_names_1.csv"))
mapping_LC_names_1 <- readRDS(file="mapping_for_G4MDSlink.RData")[[3]]

results2 <- results1 %>%
  left_join(mapping_LC_names_1,by=c("lu.from"="lu.linkoutput")) %>%
  select(-c(lu.from)) %>%
  rename(lu.from=lu.new) %>%
  left_join(mapping_LC_names_1,by=c("lu.to"="lu.linkoutput")) %>%
  select(-c(lu.to)) %>%
  rename(lu.to=lu.new) %>%
  select(REGION,times,ns,lu.to,value,lu.from)


## Merge protected lands into corrsponding LC types
# mapping_LC_names_2 <- read.csv(file=path("mapping_LC_names_2.csv"))
mapping_LC_names_2 <- readRDS(file="mapping_for_G4MDSlink.RData")[[4]]

results3 <- results2 %>%
  left_join(mapping_LC_names_2,by=c("lu.from"="lu.new")) %>%
  select(-c(lu.from)) %>%
  rename(lu.from=lu.final) %>%
  left_join(mapping_LC_names_2,by=c("lu.to"="lu.new")) %>%
  select(-c(lu.to)) %>%
  rename(lu.to=lu.final) %>%
  select(REGION,times,ns,lu.to,value,lu.from) %>%
  group_by(REGION,times,ns,lu.to,lu.from) %>%
  summarise(value=sum(value)) %>%
  ungroup()


### Of course should check here later, whether the merging is correct (similar as checking before and after allocating RstLnd_Affor to OthNatLnd)


## Prepare checking results
return.check1.init.g4mland <- resultsAll[[2]]
return.check1.excess.forest <- resultsAll[[3]]
checking_G4M_DS_results <- list(return.check1.init.g4mland,return.check1.excess.forest)

# Save results to .RData------
cat("==> Write G4M_DS_Link output to RData","\n")

# saveRDS(list(bii_simu,bii,csar_simu,csar),"Output/biodiversity.RData")
saveRDS(list(results1,results2,results3,checking_G4M_DS_results),"Output/G4M_DS_Link_result.RData")



# temp1 <- downscalr_out %>%
#   mutate(lu.from=recode(
#     lu.from,"forest_old_ha"="priforest"))

# test1 <- downscalr_out %>%
#   mutate(lu.from=recode(lu.from,"CrpLnd"="cropland"),
#          lu.to=recode(lu.to,"CrpLnd"="cropland"))



## ------Backup: Checking------
## 1. compare IsForest mapping by using diff G4M csv input for x-y-g4mid information
if(FALSE){
# map: ns to g4m 0.5degree grids
mapping <- readRDS(file='G4m_mapping.RData')[[1]]
mapping <- data.frame(apply(mapping, 2, as.numeric))

# map: isForest or not at g4m 0.5degree grids
mapping_g4mid_xy <- read.csv(file="oecd_ssp2_decc2016_reg31.csv") %>%
  dplyr::select(x,y,simuid) %>% rename("g4m_05_id"="simuid")
mapping_g4mid_xy_new <- read.csv(file="oecd_ssp2_decc2016_pbl_r31_16032022_0.csv") %>%
  dplyr::select(x,y,simuid) %>% rename("g4m_05_id"="simuid")

maplayer_isforest0 <- read.csv(file="forest_nonforest_LUH2_halfdeg.csv")
maplayer_isforest <- maplayer_isforest0 %>%
  rename("IsForest"=colnames(maplayer_isforest0)[3]) %>%
  left_join(mapping_g4mid_xy,by=c("x","y")) %>% na.omit() %>%
  dplyr::select(g4m_05_id,IsForest)
maplayer_isforest$IsForest <- as.numeric(maplayer_isforest$IsForest)

maplayer_isforest_new <- maplayer_isforest0 %>%
  rename("IsForest"=colnames(maplayer_isforest0)[3]) %>%
  left_join(mapping_g4mid_xy_new,by=c("x","y")) %>% na.omit() %>%
  dplyr::select(g4m_05_id,IsForest) %>%
  mutate(IsForest=as.numeric(IsForest))

compare_maplayer <- maplayer_isforest %>%
  left_join(maplayer_isforest_new,by=c("g4m_05_id")) %>%
  mutate(diff=IsForest.y-IsForest.x)
}

## 2. Check RstLnd for all 47 regions in CONS scenario
if(FALSE){
  #YW temp: checking the reallocation of RstLnd after DS & before G4M_DS_Link, for CONS scen (Nr:0-36)

  for(scen in 0:36){

    cat(paste0("Checking RstLnd for: ",scen1,"_",scen2,"_",scen3,", ScenNr=",scen),"\n")

    # Define scenarios
    current_scen <- scen_map %>% filter(ScenNr == scen)

    scen1 <- current_scen$SCEN1
    scen2 <- current_scen$SCEN2
    scen3 <- current_scen$SCEN3

    # Load downscalr output
    # in local PC it takes about 10 sec
    # downscalr_out <- rfunc(path(str_glue( "output_", cluster_nr, ".",
    # sprintf("%06d",current_scen$ScenNr),".RData")))

    downscalr_out0 <- rfunc(path(str_glue( "output_", cluster_nr, ".",
                                           sprintf("%06d",current_scen$ScenNr),".RData")))

    # downscalr_out0 <- rfunc(path(str_glue( "output_", cluster_nr, ".",
    # sprintf("%06d",1),".RData")))
    #  group_by(times) %>% group_split()


    # Process RstLnd in res
    res0 <- downscalr_out0
    res0$ns <- as.numeric(res0$ns)

    # 1) mapping dat to g4mid
    res1 <- res0 %>% left_join(mapping, by=c("ns"="SimUID"))# %>%
    # group_by(g4m_05_id, lu.from, times) %>%
    # summarise(value=sum(value)) %>%
    # rename("ScenYear"="times")

    # res1a <- res1 %>% left_join(maplayer_isforest_new) %>%
    # mutate(IsForest=ifelse(is.na(IsForest),0,IsForest))

    # Assign RstLnd to afforestable or non-afforestable
    res2 <- res1 %>% left_join(maplayer_isforest_new) %>%
      mutate(IsForest=ifelse(is.na(IsForest),0,IsForest)) %>%
      mutate(lu.from=ifelse(IsForest==1,
                            recode(lu.from,"RstLnd"="RstLnd_YesAffor"),
                            recode(lu.from,"RstLnd"="RstLnd_NoAffor")))%>%
      mutate(lu.to=ifelse(IsForest==1,
                          recode(lu.to,"RstLnd"="RstLnd_YesAffor"),
                          recode(lu.to,"RstLnd"="RstLnd_NoAffor")))
    # res2view <- res2 %>%
    # subset((lu.from %in% c("OthNatLnd","RstLnd_YesAffor")) & (lu.to %in% c("OthNatLnd","RstLnd_YesAffor"))) %>%
    # arrange(times,ns)

    downscalr_out <- res2


    # Reallocate Rstlnd_YesAffor to OtherNatLnd

    res3 <- res2 %>%
      mutate(lu.from=recode(lu.from,"RstLnd_YesAffor"="OthNatLnd")) %>%
      mutate(lu.to=recode(lu.to,"RstLnd_YesAffor"="OthNatLnd"))

    res4 <- res3 %>%
      # group_by(REGION,times,ns,lu.to,value,lu.from,g4m_05_id) %>% summarise(value=sum(value))%>% ungroup()
      group_by(REGION,times,ns,lu.to,lu.from,g4m_05_id) %>% summarise(value=sum(value))%>% ungroup()

    ## (Check whether it works)
    ### Check0: compare RstLnd_YesAffor and RstLnd_NoAffor
    Check0_CompRstIsForest <- res2 %>% group_by(REGION,times,lu.to) %>%
      summarise(value=sum(value)) %>%
      subset(lu.to %in% c("RstLnd_YesAffor","RstLnd_NoAffor")) %>% ungroup()

    Check1_b4AllocRstYes <- res2 %>% group_by(REGION,times,lu.from,lu.to) %>% summarise(value=sum(value))%>% ungroup()

    # the results in Check2 and Check3 for "OthNatLnd->OthNatLnd" is the same (from res3-right after rename and res4-aggregate OthNatLnd and the renamed RstLnd_YesAffor)
    Check2_afterAllocRstYes <- res3 %>% group_by(REGION,times,lu.from,lu.to) %>% summarise(value=sum(value))%>% ungroup()

    Check3_afterAllocRstYes <- res4 %>% group_by(REGION,times,lu.from,lu.to) %>% summarise(value=sum(value))%>% ungroup()


    list_CheckRstLnd <- list()
    list_CheckRstLnd[[1]] <- res0
    list_CheckRstLnd[[2]] <- res4
    list_CheckRstLnd[[3]] <- Check0_CompRstIsForest
    list_CheckRstLnd[[4]] <- Check1_b4AllocRstYes
    list_CheckRstLnd[[5]] <- Check2_afterAllocRstYes
    list_CheckRstLnd[[6]] <- Check3_afterAllocRstYes

    res_Region <- unique(Check0_CompRstIsForest$REGION)[1]

    write.csv(Check0_CompRstIsForest,file=paste0("Check_RstLnd_IsForest",scen1,"_",scen2,"_",scen3,"_",scen,"_",res_Region,".csv"))
    saveRDS(list_CheckRstLnd,paste0("Check_RstLnd_",scen1,"_",scen2,"_",scen3,"_",scen,".RData"))

  } # loop through regions to check (amount and the correctness of reallocation of) RstLnd

}


# 3. Check protected land area
if(FALSE){
  #YW temp: checking the reallocation of RstLnd after DS & before G4M_DS_Link, for CONS scen (Nr:0-36)

  for(scen in 0:36){

    cat(paste0("Checking Protected land for: ",scen1,"_",scen2,"_",scen3,", ScenNr=",scen),"\n")

    # Define scenarios
    current_scen <- scen_map %>% filter(ScenNr == scen)

    scen1 <- current_scen$SCEN1
    scen2 <- current_scen$SCEN2
    scen3 <- current_scen$SCEN3


    downscalr_out0 <- rfunc(path(str_glue( "output_", cluster_nr, ".",
                                           sprintf("%06d",current_scen$ScenNr),".RData")))


    # Process RstLnd in res

    Check0_AreaProtect <- downscalr_out0 %>% group_by(REGION,times,lu.to) %>%
      summarise(value=sum(value)) %>%
      subset(lu.to %in% c("protected_priforest","protected_other")) %>% ungroup()

    res_Region <- unique(Check0_AreaProtect$REGION)[1]

    write.csv(Check0_AreaProtect,file=paste0("Check_ProtArea",scen1,"_",scen2,"_",scen3,"_",scen,"_",res_Region,".csv"),row.names = FALSE)

  } # loop through regions to check protected land

}


# 4. Save mappings to RData so that it can be used directly when running on limpopo
if(FALSE){
  list_mapping <- list()
  list_mapping[[1]] <- mapping_g4mid_xy_new
  list_mapping[[2]] <- maplayer_isforest_new
  list_mapping[[3]] <- mapping_LC_names_1
  list_mapping[[4]] <- mapping_LC_names_2
  list_mapping[[5]] <- mapping

  saveRDS(list_mapping,"mapping_for_G4MDSlink.RData")

}  ## Section: saving mappings to .RData file for easy sourcing

# End of file ------
