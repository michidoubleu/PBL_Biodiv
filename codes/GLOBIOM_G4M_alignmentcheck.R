library(dplyr)
library(tidyr)
library(data.table)
source("./codes/helper_functions.R")

define.path <- "P:/globiom/Projects/SSPs/ScenarioMIP7/Results/lookup_table_results_withSSPV3/lookup_table_v5_16Jun2025/raw_results"
all.scenarios <- read.csv(paste0(define.path, "/Biodiversity_Link/scenario_mapping_BetterLookup2SSP2a_20250613_selectedScens.csv"))
choose.scen <- all.scenarios %>% filter(ScenLoop=="0", RegionName%in%c("IndonesiaReg", "BrazilReg"))
cluster.nr <- 2821

mapping <- readRDS(file = 'input/G4M_mapping.RData')[[1]]
mapping <- apply(mapping, 2, as.character)
mapping <- data.frame(mapping)
mapping <- mapping %>% rename("g4m_id" = "g4m_05_id", "ns" = "SimUID")

full.diff <- NULL
jjj <- 1
for(jjj in 1:nrow(choose.scen)){
curr.reg <- choose.scen$RegionName[jjj]
curr.SCEN1 <- choose.scen$SCEN1[jjj]
curr.SCEN2 <- choose.scen$SCEN2[jjj]
curr.SCEN3 <- choose.scen$SCEN3[jjj]
project <- "BetterLookup2SSP2a"
lab <- "20250613"
label <- choose.scen$ScenNr[jjj]



DS.res.full <- list.files(paste0(define.path, "/DownScale"), full.names = TRUE)
DS.res.full <- filter_paths(DS.res.full, label)
DS.res.full <- readRDS(DS.res.full)
DS.res <- DS.res.full[[3]]$out.res

G4M.res <- read.csv(paste0(define.path, "/G4M/out/area_harvest_map_",project,"_",lab,"_",curr.SCEN1,"_",curr.SCEN3,"_",curr.SCEN2,".csv"))
#G4M.res <- read.csv("P:/globiom/Projects/SSPs/ScenarioMIP7/Results/scenario_results_with_SSPV3/2025v5/final_version/raw_results/G4M/out/SSP2_H/area_harvest_map_Trunk5266_MSGfdbkOct2025_final_20251024_SSP2_SSP2_H_SCEN.csv")

G4M.res <-
  G4M.res %>% dplyr::select(g4m_id, year, forest_old_ha, forest_new_ha,forest_tot_ha) %>%
  group_by(g4m_id) %>% mutate(forest_old_ha=forest_old_ha/1000,
                              forest_new_ha=forest_new_ha/1000,forest_tot_ha=forest_tot_ha/1000) %>%
  mutate(
    deforestation = forest_old_ha - lag(forest_old_ha, n = 2),
    afforestation = forest_new_ha - lag(forest_new_ha, n = 2),
    g4m_id = as.character(g4m_id)
  ) %>%
  filter(year %% 10 == 0) %>% na.omit()

DS.res.g4mland <-
  DS.res %>% mutate(
    lu.from = recode(
      lu.from,
      "OthNatLnd" = "G4Mland",
      "Forest" = "G4Mland",
      "PriFor" = "G4Mland",
      "MngFor" = "G4Mland",
      "protected_other" =
        "G4Mland",
      "protected_priforest" =
        "G4Mland",
    ),
    lu.to = recode(
      lu.to,
      "OthNatLnd" = "G4Mland",
      "Forest" = "G4Mland",
      "PriFor" = "G4Mland",
      "MngFor" = "G4Mland",
      "protected_other" = "G4Mland",
      "protected_priforest" =
        "G4Mland",
    )
  ) %>%
  group_by(REGION, times, ns, lu.from, lu.to) %>% summarise(value = mysum(value))

DS.init.g4mland <-
  DS.res.g4mland %>%
  filter(times == 2010) %>%
  group_by(REGION, ns, lu.from, times) %>%
  summarise(value = mysum(value)) %>%
  subset(lu.from == 'G4Mland') %>%
  mutate(times = times - 10)


forest.weighting <- mapping %>%
  left_join((DS.init.g4mland %>% ungroup() %>% dplyr::select(ns, value))) %>%
  na.omit() %>% group_by(g4m_id) %>% mutate(weight = value / mysum(value)) %>% ungroup() %>% dplyr::select(g4m_id, ns, weight)


simu.forest <- forest.weighting %>%
  left_join(
    G4M.res %>% dplyr::select(g4m_id, year, forest_old_ha, forest_new_ha, forest_tot_ha) %>%
      subset(year == 2000)
  ) %>%
  mutate(forest_old_ha = forest_old_ha * weight,
         forest_new_ha = forest_new_ha * weight,
         forest_tot_ha = forest_tot_ha * weight) %>%
  na.omit() %>%
  dplyr::select(ns, forest_old_ha, forest_new_ha, forest_tot_ha)


DS.init.g4mland.forest <-
  DS.init.g4mland %>%
  left_join(simu.forest %>%
              mutate(ns = as.character(ns))) %>%
  mutate(across(starts_with("forest"), ~ifelse(is.na(.), 0, .))) %>%
  mutate(OthNatLnd = value - forest_old_ha - forest_new_ha)

save.diff <- DS.init.g4mland.forest %>% group_by(REGION) %>% summarise(GLOBIOM_max=sum(value,na.rm = T),G4M=sum(forest_tot_ha,na.rm = T),diff=sum(OthNatLnd,na.rm = T)) %>%mutate(SCEN1=curr.SCEN1,SCEN2=curr.SCEN2,SCEN3=curr.SCEN3)


full.diff <- full.diff %>% bind_rows(save.diff)




if (!all(DS.init.g4mland.forest$OthNatLnd > 0)) {
  warning(
    "For some simus there is more G4M forest than G4Mland \n
            (sum of PriFor, MngFor, OthNatland + protected! \n
            Forest reduced to G4Mland levels"
  )
  diff <- sum(DS.init.g4mland.forest$OthNatLnd[DS.init.g4mland.forest$OthNatLnd < 0])
  DS.init.g4mland.forest <- DS.init.g4mland.forest %>%
    mutate(
      forest_old_ha = ifelse(
        OthNatLnd < 0,
        forest_old_ha + (forest_old_ha / (forest_old_ha + forest_new_ha)) * OthNatLnd,
        forest_old_ha
      ),
      forest_new_ha = ifelse(
        OthNatLnd < 0,
        forest_new_ha + (forest_new_ha / (forest_old_ha + forest_new_ha)) * OthNatLnd,
        forest_new_ha
      ),
      OthNatLnd = ifelse(OthNatLnd < 0, 0, OthNatLnd)
    )
}

final.forest.alloc <-
  DS.init.g4mland.forest %>% ungroup() %>% dplyr::select(REGION, times, ns, forest_old_ha, forest_new_ha, OthNatLnd) %>%
  pivot_longer(
    cols = c(forest_old_ha, forest_new_ha, OthNatLnd),
    names_to = "lu.from",
    values_to = "value"
  )

sum(final.forest.alloc$value)



DS.res.after <- list.files(paste0(define.path, "/Biodiversity_Link"), full.names = TRUE)
DS.res.after <- filter_paths(DS.res.after, label)
DS.res.after <- readRDS(DS.res.after)
DS.res.new <- DS.res.after[[6]]



DS.res.new <-
  DS.res.new %>%
  filter(times == 2010) %>%
  group_by(REGION, ns, lu.from, times) %>%
  summarise(value = mysum(value)) %>%
  subset(lu.from %in% c("forest_old_ha", "forest_new_ha")) %>%
  mutate(times = times - 10)

sum(DS.res.new$value)

}


G4M.map <- read.csv('P:/globiom/Projects/uk-desnz/DESNZ_22082025/merge/DESNZ_22082025/tmp/glc_glfm.csv')
G4M.map.indo <- G4M.map %>% filter(iso=="IDN")
G4M_IDS <- unique(G4M.map.indo$g4m_05_id)

G4M.res <- read.csv("P:/globiom/Projects/SSPs/ScenarioMIP7/Results/scenario_results_with_SSPV3/2025v5/final_version/raw_results/G4M/out/SSP2_H/area_harvest_map_Trunk5266_MSGfdbkOct2025_final_20251024_SSP2_SSP2_H_SCEN.csv")
G4M.res <- G4M.res %>% filter(year==2000, g4m_id %in% G4M_IDS)
sum(G4M.res$forest_tot_ha)


dim(G4M.map[which(G4M.map$iso=='IDN'),])
dim(G4M.map[which(G4M.map$iso=='IDN' & !is.na(G4M.map$g4m_05_id)),])
























path_to_gridded <- 'P:/globiom/Projects/SSPs/ScenarioMIP7/Results/scenario_results_with_SSPV3/2025v5/final_version/raw_results/G4M/out/SSP2_H/area_harvest_map_Trunk5266_MSGfdbkOct2025_final_20251024_SSP2_SSP2_H_SCEN.csv'
gridded_SCENMIP7_SSP2_H <- read.csv(path_to_gridded)
path_to_map <- 'P:/globiom/Projects/uk-desnz/DESNZ_22082025/merge/DESNZ_22082025/tmp/glc_glfm.csv'
map_G4M_id <- read.csv(path_to_map) |>
  dplyr::select(x,y,code,iso,cell_area_fm_kha)
ag_map <- gridded_SCENMIP7_SSP2_H |> 
  left_join(map_G4M_id,by = join_by(x, y)) |>
  group_by(iso,year) |>
  summarise(area_forest_old_ha=sum(forest_old_ha),
            area_forest_new_ha=sum(forest_new_ha),
            area_forest_tot_ha=sum(forest_tot_ha),.groups = 'drop')
ag_map.extract <- ag_map |>
  dplyr::filter(iso%in%c('BRA','IDN'),year==2000) |>
  dplyr::mutate(Country=case_when(iso=='BRA'~ 'Brazil', iso=='IDN'~'Indonesia'),
                SCEN='SSP2_SSP2_H_SCEN') |>
  dplyr::select(Country,SCEN,area_forest_new_ha,area_forest_old_ha,area_forest_tot_ha) |>
  dplyr::mutate(source='gridded')



test <- gridded_SCENMIP7_SSP2_H |> 
  left_join(map_G4M_id,by = join_by(x, y)) |>
  dplyr::filter(iso%in%c('IDN'),year==2000)
sum(test$forest_tot_ha)



G4M.map <- read.csv('P:/globiom/Projects/uk-desnz/DESNZ_22082025/merge/DESNZ_22082025/tmp/glc_glfm.csv')
G4M.map.indo <- G4M.map %>% filter(iso=="IDN")
G4M_IDS <- unique(G4M.map.indo$g4m_05_id)

G4M.res <- read.csv("P:/globiom/Projects/SSPs/ScenarioMIP7/Results/scenario_results_with_SSPV3/2025v5/final_version/raw_results/G4M/out/SSP2_H/area_harvest_map_Trunk5266_MSGfdbkOct2025_final_20251024_SSP2_SSP2_H_SCEN.csv")
G4M.res <- G4M.res %>% filter(year==2000, g4m_id %in% G4M_IDS)
sum(G4M.res$forest_tot_ha)



setdiff(G4M_IDS, test$g4m_id) 

