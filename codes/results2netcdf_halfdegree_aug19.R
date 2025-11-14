# save results to netcdf file
rm(list=ls())

#require(ncdf)
require(raster)
# LIBRARIES ---------------------------------------------------------------
library(ncdf4) 
library(ncdf4.helpers)
library(PCICt)
library(lattice)
library(ggplot2)
require(reshape2)
require(stringr)

p4s = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

#source_path = "./output/"
#target_path = "../../../../large_output/"
#target_path = "./output/"
#source_path = "P:/globiom/Projects/FABLE/BendingTheCurve/Results/downscaled/"
#target_path = "P:/globiom/Projects/FABLE/BendingTheCurve/Results/downscaled/"
#date_tag = "low2p_1feb18"
# load from file
settings = read.csv("_settings.csv",stringsAsFactors = FALSE)
source_path = settings$source_path
#source_path = "./"
target_path = settings$target_path
date_tag = settings$date_tag
#source_path = "P:/globiom/Projects/FABLE/BendingTheCurve/Results/downscaled_comparison/"
#target_path = settings$final_target_path
#target_path = "P:/globiom/Projects/FABLE/BendingTheCurve/Results/downscaled_4pierro/"

 # scenarios = data.frame(
 #   RCPs = c('RCPref_SSP2_NOBIOD'),
 #   SPAs = c('SPA0'),
 #   SSPs = c('SSP2'),
 #   stringsAsFactors = FALSE
 # )
# load from file
scenarios = read.csv("_scenarios.csv",stringsAsFactors = FALSE)

full_simu_map = read.csv("full_simu_map_biodiv.csv",stringsAsFactors = FALSE)
na_sum = function(x) {return(sum(x[!is.na(x)]))}

area_open <- nc_open('template/BendingTheCurveFT-LCproj-MODEL-SCEN-DATE.nc')
#print(area_open)
area_variable_Mha <- ncvar_get(area_open, "pixel_area")

#add_file <- nc_open('to_be_added_to_other.nc')
#print(area_open)
#fix_toAddToOther <- ncvar_get(add_file, "fix_toAddToOther")
add_file = read.csv("half_degree_share_pixel_area_underSimU0.csv",stringsAsFactors = FALSE)
## coordinates not fully correct, set right
add_file$POINT_X = round(add_file$POINT_X,2)
add_file$POINT_Y = round(add_file$POINT_Y,2)
# values not completelz shares, correct
add_file = add_file[-which(add_file$share_pixel_under_SimU0 == 0),]
add_file$share_pixel_under_SimU0[add_file$share_pixel_under_SimU0 > 1] = 1
area_simu0_share = rasterFromXYZ(add_file[,c("POINT_X","POINT_Y","share_pixel_under_SimU0")],res = 0.5,crs = p4s)
area_simu0_share =  extend(area_simu0_share,extent(-180,180,-90,90))
area_simu0_share = t(getValues(area_simu0_share,format = "matrix"))

#### create colrow raster
colrows = unique(full_simu_map$colrowID)
colrows = data.frame(
  colrows = colrows, 
  x = as.double(substr(colrows,3,5)), y = as.double(substr(colrows,6,8)),
  id = 1:length(colrows)
)
loncr = function(x) {89.75 - (x-1) * .5}
latcr = function(x) {-179.75 + (x-1) * .5}
loncr_inv = function(x) { (x + 89.75) / .5 + 1}
latcr_inv = function(x) { (x + 179.75) / .5 + 1}
colrows = data.frame(
  colrows,
  long = loncr(colrows$y),
  lat = latcr(colrows$x)
)
ttt1 = paste0("CR",
              sprintf("%03d",latcr_inv(add_file$POINT_X)),
              sprintf("%03d",loncr_inv(add_file$POINT_Y)))
add_file = cbind(add_file,colrow = ttt1)
add_file_colrow = aggregate(add_file$share_pixel_under_SimU0 * add_file$Shape_Area,
                           by = list(colrow = add_file$colrow),
                           FUN = sum, na.rm = T)


#GLOBIOM's simulation units; we have to convert to this
geosims <- raster("./shape_GLOBIOM/raster_simU/w001001.adf")
proj4string(geosims)<-p4s #apply projection manually
save_geovals <- data.frame(getValues(geosims)) #just the SimU IDs
colnames(save_geovals)<-"SimUID"

geosims_area <- raster::area(geosims)
ttemp = getValues(geosims_area)
ttemp[is.na(getValues(geosims))] = NA
geosims_area = setValues(geosims_area,ttemp)
### calculate SimU percents of areas

### protected areas
## static raster in BIOD scenarios
protectedSimU = read.csv("forbidBII_10Dec2018.csv",stringsAsFactors = FALSE)
protectedSimU = rowSums(protectedSimU == 0)
protectedSimU[protectedSimU > 0] = 1
protectedCR = aggregate(list(value = protectedSimU),by = list(colrow = full_simu_map$colrowID),FUN = sum)
protectedCR = protectedCR[match(colrows$colrows,protectedCR$colrow),]
protectedCR$value[protectedCR$value > 0] = 1
ddata = cbind(x = colrows$lat,y = colrows$long,z=protectedCR$value)
protectedMap = rasterFromXYZ(ddata,crs = p4s)
protectedMap =  extend(protectedMap,extent(-180,180,-90,90))
not_protectedMap = protectedMap
data = getValues(not_protectedMap)
data[data == 1] = 0
not_protectedMap = setValues(not_protectedMap,data)

## protected pct of pixel from 2000 onwards
LU2000 = read.csv("SimU_LU_biodiv_G4M_jan19.csv",stringsAsFactors = FALSE)
protected_othfor_CR = aggregate(list(
  value = (LU2000$protected_priforest + LU2000$protected_other)/10^6
),by = list(colrow = full_simu_map$colrowID),FUN = sum)
protected_othfor_CR = protected_othfor_CR[match(colrows$colrows,protected_othfor_CR$colrow),]
ddata = cbind(x = colrows$lat,y = colrows$long,z=protected_othfor_CR$value)
protected_othfor_Map = rasterFromXYZ(ddata,crs = p4s)
protected_othfor_Map =  extend(protected_othfor_Map,extent(-180,180,-90,90))
data = t(as.matrix(protected_othfor_Map)); dataNA = is.na(data)
data = data / area_variable_Mha
data[is.na(data)] = 0; data[data > 1] = 1
data[dataNA] = NA
protected_othfor_Map = setValues(protected_othfor_Map,t(data))

####
load_date_tag = settings$date_tag

to00 = read.csv("to00.csv",row.names = 1)
to10 = read.csv("to10.csv",row.names = 1)

### LU & LUC classes that are extracted drom the downscaled output 
from_file2000_LU_name = c("cropland","grassland","priforest","mngforest",
                           "SRP","restored","other","protected_priforest","protected_other",
                           "urban")
from_file_LUC_name = c("cropland.grassland","cropland.mngforest","cropland.SRP","cropland.restored","cropland.other",                         
              "grassland.cropland","grassland.mngforest","grassland.SRP","grassland.restored","grassland.other",                        
              "priforest.cropland","priforest.grassland","priforest.mngforest",
              "SRP.mngforest","SRP.other",
              "other.cropland","other.grassland","other.mngforest","other.SRP")

### processing names
# without the specific year tracking for abandoned flows
process_LU_name_short = c("cropland","grassland","priforest","mngforest",
                          "SRP","restored","other","urban",
                          "cropland.other","SRP.other","grassland.other","mngforest.other")
# with year tracking
process_LU_name = c("cropland","grassland","priforest","mngforest",
                         "SRP","restored","other","urban",
                         paste0("cropland.other_",seq(2010,2100,by = 10)),
                         paste0("SRP.other_",seq(2010,2100,by = 10)),
                         paste0("grassland.other_",seq(2010,2100,by = 10)),
                         paste0("mngforest.other_",seq(2010,2100,by = 10)))
process_LUC_names_df = str_split_fixed(from_file_LUC_name,"\\.",2)
process_LUC_names_df[from_file_LUC_name == "cropland.other",2] = "cropland.other"
process_LUC_names_df[from_file_LUC_name == "SRP.other",2] = "SRP.other"
process_LUC_names_df[from_file_LUC_name == "grassland.other",2] = "grassland.other"



#just the LU classes that we export
to_file_EXP1to3_LU_name = c("cropland","SRP","grassland","priforest","mngforest",
                         "restored","other","urban")
to_file_EXP4to5_LU_name = c("cropland","SRP","grassland","priforest","mngforest",
                         "restored","other","urban",
                         "cropland.other","SRP.other","grassland.other","mngforest.other")


scenarios = scenarios[-1,]


sss =1
# sss = 2; sss = 7
for (sss in 1:nrow(scenarios)) {
  RCP = scenarios$RCPs[sss]
  SPA = scenarios$SPAs[sss]
  SSP = scenarios$SSPs[sss]
  cat(RCP,SPA,SSP,"\n")
  
  T = 10
  MIN_T = 0
  MAX_T = 90
  
  ftarget = paste(source_path,"resJan19_",SPA,"_",RCP,"_",SSP,"_",
                  as.character(load_date_tag),"__luc.csv",sep="")
  res_file <-read.csv(ftarget ,stringsAsFactors=FALSE)
  
  ### calculate and adjust initial starting values
  ### just do this for the first scenario
  if (sss == 1) {
    ## LU in 2000 - starting value for output
    res_file2000 = res_file[res_file$Year == 2000,]
    res_file_LU2000 = aggregate(res_file2000[,c("Area",from_file2000_LU_name)] ,
                            by=list(colrow=res_file2000$Colrow),FUN=sum)
    res_file_LU2000$priforest = res_file_LU2000$priforest + res_file_LU2000$protected_priforest
    res_file_LU2000$other = res_file_LU2000$other + res_file_LU2000$protected_other
    res_file_LU2000 = cbind(res_file_LU2000,cropland.other = 0,grassland.other = 0,mngforest.other = 0,SRP.other=0)
    # reorder classes to match template
    res_file_LU2000 = res_file_LU2000[,c("colrow","Area",to_file_EXP4to5_LU_name)]
    res_file_LU2000 = melt(res_file_LU2000,id.vars = c(1:2),variable.name = "LU")
    # convert to Mha
    res_file_LU2000$Area = res_file_LU2000$Area / 10^6
    res_file_LU2000$value = res_file_LU2000$value / 10^6
    
    ### add/check with missing other
    # match to colrow for plotting
    # Define some straightforward dimensions
    lon <- seq(-179.75, 179.75, by = 0.5)
    lat <- seq(89.75, -89.75, by = -0.5)
    time <- seq(2010,2100,by=10)
    CLASSES = 8
    lc_class <- seq(1,CLASSES,1)
    dim_lon <- ncdf4::ncdim_def("lon", "degrees_east",lon)
    dim_lat <- ncdf4::ncdim_def("lat", "degrees_north",lat)
    dim_time <- ncdf4::ncdim_def("time", "years", calendar="standard",time,unlim=TRUE)
    o = match(res_file_LU2000$colrow,colrows$colrows)
    LandCover_pixshare_array <- array(NA, dim = c(length(lon), length(lat), length(lc_class)),
                                      dimnames = list(lon, lat, lc_class))
    CLASSES_NAME = to_file_EXP1to3_LU_name
    for (c in 1:CLASSES) {
      to_plot = res_file_LU2000$value[res_file_LU2000$LU == CLASSES_NAME[c]]
      
      ddata = cbind(x = colrows$lat[o],y = colrows$long[o],z=to_plot)
      temp_map = rasterFromXYZ(ddata,crs = p4s)
      temp_map =  extend(temp_map,extent(-180,180,-90,90))
      data <- t(as.matrix(temp_map)) / area_variable_Mha
      
      data[is.na(data) & !is.na(area_simu0_share)] = 0
      if (CLASSES_NAME[c] == "other") {
       ttemp = area_simu0_share
       ttemp[is.na(ttemp) & !is.na(data)] = 0
       data = data + ttemp
      } 
      LandCover_pixshare_array[,,c] = data
      
    }
    pix_totals = apply(LandCover_pixshare_array,c(1,2),sum) # sum over all LU - classes
    # take away first from area_simu_share0, then from other
    problem_pix_totals = which(pix_totals>1,arr.ind = T)
    problem_pix_totals = data.frame(problem_pix_totals,
                               share = pix_totals[problem_pix_totals],
                               pix_area_Mha = 0,
                               err_share = 0,
                               simu0_share = 0)
    problem_pix_totals$err_share = problem_pix_totals$share - 1
    # match with SimU0
    for (ii in 1:nrow(problem_pix_totals)) {
      problem_pix_totals$pix_area_Mha[ii] = area_variable_Mha[problem_pix_totals$row[ii],
                                                            problem_pix_totals$col[ii]]
      problem_pix_totals$simu0_share[ii] = area_simu0_share[problem_pix_totals$row[ii],
                                                            problem_pix_totals$col[ii]]
    }
    problem_pix_totals$simu0_share[is.na(problem_pix_totals$simu0_share)] = 0
    problem_pix_totals = cbind(problem_pix_totals,
                               colrow = paste0("CR",
                                             sprintf("%03d",problem_pix_totals[,1]),
                                             sprintf("%03d",problem_pix_totals[,2])),
                               err_Mha = problem_pix_totals$err_share * problem_pix_totals$pix_area_Mha,
                               simu0_Mha = problem_pix_totals$simu0_share * problem_pix_totals$pix_area_Mha)
    ### allocate error to simu0
    ttemp = problem_pix_totals$simu0_Mha - problem_pix_totals$err_Mha
    simu0_Mha_adj = ttemp; simu0_Mha_adj[simu0_Mha_adj<0] = 0
    err_Mha_non_simu0 = ttemp; err_Mha_non_simu0[err_Mha_non_simu0>0] = 0
    err_Mha_non_simu0 = abs(err_Mha_non_simu0)
    problem_pix_totals = cbind(problem_pix_totals,
                               simu0_Mha_adj = simu0_Mha_adj, err_Mha_non_simu0 = err_Mha_non_simu0)
    # adjust simu0
    for (ii in 1:nrow(problem_pix_totals)) {
      area_simu0_share[problem_pix_totals$row[ii],
                       problem_pix_totals$col[ii]] = problem_pix_totals$simu0_Mha_adj[ii] / problem_pix_totals$pix_area_Mha[ii]
    }
    ### allocate to other
    other_colrow_index = which(res_file_LU2000$LU == "other" & res_file_LU2000$colrow %in% problem_pix_totals$colrow)
    problem_pix_totals = problem_pix_totals[
      match(res_file_LU2000$colrow[other_colrow_index],problem_pix_totals$colrow),]
    ttemp = res_file_LU2000[other_colrow_index,]
    ttemp$value = ttemp$value - problem_pix_totals$err_Mha_non_simu0
    # some values could not be allocated, substract these from cropland
    problem_colrow = ttemp$colrow[which(ttemp$value<0)]
    problem_values = abs(ttemp$value[which(ttemp$value<0)])
    res_file_LU2000$value[res_file_LU2000$colrow %in% problem_colrow & res_file_LU2000$LU=="cropland"] =
      res_file_LU2000$value[res_file_LU2000$colrow %in% problem_colrow & res_file_LU2000$LU=="cropland"] - problem_values
    # put ttemp others back into res_file_LU2000
    ttemp$value[ttemp$value<0] = 0
    res_file_LU2000$value[other_colrow_index] = ttemp$value
  }
  
  ## aggregate LUC res_file
  res_file_LUC = aggregate(res_file[,c("Area",from_file_LUC_name)] ,
                          by=list(colrow=res_file$Colrow,Year = res_file$Year),FUN=sum)
  res_file_LUC = melt(res_file_LUC,id.vars = c(1:3),variable.name = "LU")
  # convert to Mha
  res_file_LUC$Area = res_file_LUC$Area / 10^6
  res_file_LUC$value = res_file_LUC$value / 10^6
  
  #### set up files for storage 
  # current land-use
  curr_LU = dcast(res_file_LU2000,colrow + Area ~ LU)
  curr_LU = curr_LU[,c("colrow","Area",
                       colnames(curr_LU)[colnames(curr_LU) %in% process_LU_name])]
  add_colnames = process_LU_name[!process_LU_name %in% colnames(curr_LU)]
  add_cols = matrix(0,nrow(curr_LU),length(add_colnames)); colnames(add_cols) = add_colnames
  curr_LU = cbind(curr_LU,add_cols)  
  curr_LU = melt(curr_LU,id.vars = c(1:2),variable.name = "LU")
  
  # save LUC as process file
  process_LUC = res_file_LUC
  
  # list of current land-uses
  LUs_2010to2100 = list()
  
  ##### allocate land use changes
  yyears = seq(2010,2100,by = 10)
  
  # loop through years and allocate land-use change
  # flows to other (abandoned land) are tracked specifically per year
  #  if not enough land in other.crop/grass/mng/SRP --> take from abandoned flows from 30,20,10,00 years ago
  #  all other not enough flows are ignored
  t = 1
  for (t in 1:T) {
    cat(yyears[t],"..")
    
    #### substract from LUs through loop
    ###   check whether we have to lower LUC flows (because there is not enough LU present)
    ###   adjust flows if necessary
    ccc = 1
    for (ccc in 1:length(process_LU_name_short)) {
      currClass = process_LU_name_short[ccc]
      
      # first adjust LUC
      outflow_id =  which(process_LUC_names_df[,1] == currClass)
      if (length(outflow_id) > 0) {
        from_land = process_LUC[process_LUC$Year == yyears[t] & 
                                  process_LUC$LU %in% from_file_LUC_name[outflow_id],]
        if (nrow(from_land)>0) {
          from_land = dcast(from_land,colrow ~ LU)
          # check values
          if (ncol(from_land) > 2) {
            rs_from_land = rowSums(from_land[,-1])
          } else {
            rs_from_land = from_land[,-1]
          }
          err = curr_LU$value[curr_LU$LU == currClass] - rs_from_land
          err[err>0] = 0; err = abs(err)
          if (any(err > 0)) {
            ## for other attempt to allocate missing LUC from previous abandoned land
            if (currClass == "other") {
              ## allocate from past abandoned
              abn_years = seq(max(yyears[t]-30,2020),max(yyears[t],2020),by = 10)
              abn_years = abn_years[-length(abn_years)]
              if (length(abn_years) > 0) {
                for (yyy in 1:length(abn_years)) {
                  if (any(err>0)) {
                    abn_labels = paste0(c("cropland.other","SRP.other","grassland.other"),"_",abn_years[yyy])
                    curr_abn = curr_LU[curr_LU$LU %in% abn_labels,]
                    curr_abn = dcast(curr_abn,colrow ~ LU)
                    shares_from_abn = curr_abn[,-1] / rowSums(curr_abn[,-1])
                    shares_from_abn[is.na(shares_from_abn)] = 0
                    ttemp = curr_abn[,-1] - shares_from_abn*err
                    ttemp[ttemp<0] = 0
                    allocated_err = rowSums(curr_abn[,-1]) - rowSums(ttemp)
                    err = err - allocated_err
                    # substract
                    curr_LU$value[curr_LU$LU %in% abn_labels] = 
                      melt(cbind(curr_abn$colrow,ttemp),c(1),variable.name = "LU")$value
                    # add
                    curr_LU$value[curr_LU$LU == currClass] = 
                      curr_LU$value[curr_LU$LU == currClass] + allocated_err
                  }
                }
              }
              # # if any residuals still present -- allocate from restored
              # if (any(err>0)) {
              #   curr_restored = curr_LU[curr_LU$LU == "restored",]
              #   ttemp = curr_restored$value - err
              #   ttemp[ttemp<0] = 0
              #   allocated_err = curr_restored$value - ttemp
              #   err = err - allocated_err
              #   # substract
              #   curr_LU$value[curr_LU$LU == "restored"] = ttemp
              #   # add
              #   curr_LU$value[curr_LU$LU == currClass] = 
              #     curr_LU$value[curr_LU$LU == currClass] + allocated_err
              # }
            }
            # otherwise set error values to zero and ignore
            # if any residuals still present -- ignore flow values
            if (any(err>0)) {
              share_from_land = from_land[,-1] / rs_from_land
              share_from_land[is.na(share_from_land)] = 0
              ttemp = from_land[,-1] - share_from_land*err
              ttemp[ttemp<0] = 0
              allocated_err = rowSums(from_land[,-1]) - rowSums(ttemp)
              err = err - allocated_err
              process_LUC$value[process_LUC$Year == yyears[t] & 
                                  process_LUC$LU %in% from_file_LUC_name[outflow_id]] =
                melt(cbind(from_land$colrow,ttemp),c(1),variable.name = "LU")$value
            }
          }
        }
      }
    }
    
    #### substract corresponding values from LUs
    ccc = 1
    for (ccc in 1:length(process_LU_name_short)) {
      currClass = process_LU_name_short[ccc]
      
      # first adjust LUC
      outflow_id =  which(process_LUC_names_df[,1] == currClass)
      if (length(outflow_id) > 0) {
        from_land = process_LUC[process_LUC$Year == yyears[t] & 
                                  process_LUC$LU %in% from_file_LUC_name[outflow_id],]
        if (nrow(from_land)>0) {
          from_land = dcast(from_land,colrow ~ LU)
          # check values
          if (ncol(from_land) > 2) {
            rs_from_land = rowSums(from_land[,-1])
          } else {
            rs_from_land = from_land[,-1]
          }
          curr_LU$value[curr_LU$LU == currClass] =
            curr_LU$value[curr_LU$LU == currClass] - rs_from_land
          
          #if (any(curr_LU$value < 0)) {
          #  stop("Error in LUC allocations, negative land!")
          #}
        }
      }
    }
    if (any(curr_LU$value <  0)) {
      curr_LU$value[curr_LU$value<0] = 0
    }
    
    #### add corresponding values to LUs
    ccc = 1
    for (ccc in 1:length(process_LU_name_short)) {
      currClass = process_LU_name_short[ccc]
      if (currClass %in% c("cropland.other","SRP.other","grassland.other","mngforest.other")) {
        out_currClass = paste0(currClass,"_",yyears[t])
      } else {
        out_currClass = currClass
      }
      
      inflow_id =  which(process_LUC_names_df[,2] == currClass)
      if (length(inflow_id) > 0) {
        to_land = process_LUC[process_LUC$Year == yyears[t] & 
                                process_LUC$LU %in% from_file_LUC_name[inflow_id],]
        if (nrow(to_land)>0) {
          to_land = dcast(to_land,colrow ~ LU,fun.aggregate = sum)
          
          if (ncol(to_land) > 2) {
            rs_to_land = rowSums(to_land[,-1])
          } else {
            rs_to_land = to_land[,-1]
          }
          curr_LU$value[curr_LU$LU == out_currClass] =
            curr_LU$value[curr_LU$LU == out_currClass] + rs_to_land 
        }
      }
    }
    # save LU
    LUs_2010to2100[[yyears[t]]] = curr_LU
  }
  
  ###### plot to netCDF - one for each EXP
  david_versions = c("EXP1","EXP2","EXP3","EXP4","EXP5")
  
  #### create netCDF storage variables
  output_netCDFs = list()
  david_version = "EXP1"
  for (david_version in david_versions) {
    if (david_version %in% c("EXP4","EXP5")) {
      CLASSES = length(to_file_EXP4to5_LU_name)
      CLASSES_NAME = to_file_EXP4to5_LU_name
    } else {
      CLASSES = length(to_file_EXP1to3_LU_name)
      CLASSES_NAME = to_file_EXP1to3_LU_name
    }
    lc_class <- seq(1,CLASSES,1)
    
    LandCover_pixshare_array <- array(NA, dim = c(length(lon), length(lat), length(lc_class), length(time)),dimnames = list(lon, lat, lc_class, time))
    output_netCDFs[[david_version]] = LandCover_pixshare_array
  }
  
  ### loop through time and david_versions and calculate scenarios
  ### also check for consistency
  # pre-calculate matching
  o = match(res_file_LU2000$colrow[res_file_LU2000$LU == "cropland"],colrows$colrows)
  t = 1
  for (t in 1:T) {
    cat(yyears[t],"-")
    # load current LU
    curr_LU = LUs_2010to2100[[yyears[t]]]
    
    # calculate time lags
    yrs_till_now = seq(2010,yyears[t],by = 10)[-1]
    if (yyears[t] == 2010) {
      last30yrs = c()
    } else {
      last30yrs = seq(max(yyears[t]-30,2020),max(yyears[t],2020),by = 10)
    }
    upto_30yrs_ago = yrs_till_now[!yrs_till_now %in% c(last30yrs,yyears[t])]
    
    david_version = "EXP1"
    for (david_version in david_versions) {
      cat(david_version,",")
      if (david_version %in% c("EXP4","EXP5")) {
        CLASSES = length(to_file_EXP4to5_LU_name)
        CLASSES_NAME = to_file_EXP4to5_LU_name
      } else {
        CLASSES = length(to_file_EXP1to3_LU_name)
        CLASSES_NAME = to_file_EXP1to3_LU_name
      }
      LandCover_pixshare_array = output_netCDFs[[david_version]]
      
      ccc = 1
      for (ccc in 1:CLASSES) {
        currClass = CLASSES_NAME[ccc]
        ### summing up rules for abandoned land
        if (currClass %in% c("cropland","grassland","SRP")) {
          to_plot = curr_LU$value[curr_LU$LU == currClass]
          ## EXP3 - add also last 30 years of abandoned [not 2010]
          if (david_version %in% c("EXP3")) {
            if (length(last30yrs) > 0) {
              add_classes = paste0(currClass,".other_",last30yrs)
              ttemp = curr_LU[curr_LU$LU %in% add_classes,]
              ttemp = aggregate(ttemp$value,by = list(ttemp$colrow),FUN = sum)
              to_plot = to_plot + ttemp$x
            }
          }
        } else if (currClass %in% c("cropland.other","SRP.other",      
                             "grassland.other","mngforest.other")) {
          to_plot = rep(0,length(o))  #EXP4 all times, EXP5:2010
          # EXP5 sum up over abandoned last 30 years [not 2010]
          if (david_version == "EXP5") {
            if (length(last30yrs) > 0) {
              add_classes = paste0(currClass,"_",last30yrs)
              ttemp = curr_LU[curr_LU$LU %in% add_classes,]
              ttemp = aggregate(ttemp$value,by = list(ttemp$colrow),FUN = sum)
              to_plot = to_plot + ttemp$x
            }
          }
        } else if (currClass == "restored") {
          to_plot = curr_LU$value[curr_LU$LU == currClass]
          ## EXP2/EXP4 - add sum over all abandoned till now [not 2010]
          if (david_version %in% c("EXP2","EXP4")) {
            if (length(yrs_till_now) > 0) {
              add_classes = apply(expand.grid(c("cropland.other","SRP.other",
                                                "grassland.other","mngforest.other"),
                                              yrs_till_now),c(1),paste,collapse="_")
              ttemp = curr_LU[curr_LU$LU %in% add_classes,]
              ttemp = aggregate(ttemp$value,by = list(ttemp$colrow),FUN = sum)
              to_plot = to_plot + ttemp$x
            } 
          }
          ## EXP3/5 - add sum over all abandoned from over 30 years ago [not 2010]
          if (david_version %in% c("EXP3","EXP5")) {
            if (length(upto_30yrs_ago) > 0) {
              add_classes = apply(expand.grid(c("cropland.other","SRP.other",
                                                "grassland.other","mngforest.other"),
                                              upto_30yrs_ago),c(1),paste,collapse="_")
              ttemp = curr_LU[curr_LU$LU %in% add_classes,]
              ttemp = aggregate(ttemp$value,by = list(ttemp$colrow),FUN = sum)
              to_plot = to_plot + ttemp$x
            } 
          }
        } else if (currClass == "other") {
          to_plot = curr_LU$value[curr_LU$LU == currClass]
          # EXP1 - add sum over all abandoned till now
          if (david_version %in% c("EXP1")) {
            add_classes = apply(expand.grid(c("cropland.other","SRP.other",
                                                "grassland.other","mngforest.other"),
                                              c(2010,yrs_till_now)),c(1),paste,collapse="_")
            ttemp = curr_LU[curr_LU$LU %in% add_classes,]
            ttemp = aggregate(ttemp$value,by = list(ttemp$colrow),FUN = sum)
            to_plot = to_plot + ttemp$x 
          }
          # EXP2-5 - add sum over all abandoned in 2010
          if (david_version %in% c("EXP2","EXP3","EXP4","EXP5")) {
            add_classes = paste0(c("cropland.other","SRP.other",
                                                "grassland.other","mngforest.other"),
                                              "_2010")
              ttemp = curr_LU[curr_LU$LU %in% add_classes,]
              ttemp = aggregate(ttemp$value,by = list(ttemp$colrow),FUN = sum)
              to_plot = to_plot + ttemp$x 
          }
        } else {
          to_plot = curr_LU$value[curr_LU$LU == currClass]
        }
        
        ddata = cbind(x = colrows$lat[o],y = colrows$long[o],z=to_plot)
        temp_map = rasterFromXYZ(ddata,crs = p4s)
        temp_map =  extend(temp_map,extent(-180,180,-90,90))
        data <- t(as.matrix(temp_map)) / area_variable_Mha
        
        data[is.na(data) & !is.na(area_simu0_share)] = 0
        if (currClass == "other") {
         ttemp = area_simu0_share
         ttemp[is.na(ttemp) & !is.na(data)] = 0
         data = data + ttemp
        }
        
        LandCover_pixshare_array[,,ccc,t] = data
      }
      # test for range of LU & sum of LU
      TOL_VAL = 10^-7
      rrange1 = range(LandCover_pixshare_array,na.rm=T)
      if (any(rrange1< 0-TOL_VAL) || any(rrange1> 1 + TOL_VAL)) {
        stop("NetCDF LU share error!")
      }
      pix_totals = apply(LandCover_pixshare_array[,,,t],c(1,2),sum) # sum over all LU - classes
      rrange2 = range(pix_totals,na.rm = T)
      if (any(rrange2< 0-TOL_VAL) || any(rrange2> 1 + TOL_VAL)) {
        stop("NetCDF LU share error!")
      }
      pix_ranges = apply(LandCover_pixshare_array[,,,t],c(3),range,na.rm = T)
      
      output_netCDFs[[david_version]] = LandCover_pixshare_array
    }
    ### test here for consistency, that is EXP3[crp] - EXP2[crp] = EXP5[crp.abn]
    exp2 = output_netCDFs[["EXP2"]][,,,t]
    exp3 = output_netCDFs[["EXP3"]][,,,t]
    exp5 = output_netCDFs[["EXP5"]][,,,t]
    lcc = "cropland"
    for (lcc in c("cropland","grassland","SRP")) {
      tt1 = exp3[,,which(to_file_EXP1to3_LU_name == lcc)] - 
        exp2[,,which(to_file_EXP1to3_LU_name == lcc)]
      tt2 = exp5[,,which(to_file_EXP4to5_LU_name == paste0(lcc,".other"))] - tt1
      if ( any(abs(range(tt2,na.rm=T))> TOL_VAL ) ) {
        stop(lcc,": EXP1-5 test failed!")
      }
    }
    cat("test.\n")
  }
  exp1 = output_netCDFs[["EXP1"]]
  exp2 = output_netCDFs[["EXP2"]]
  exp3 = output_netCDFs[["EXP3"]]
  exp4 = output_netCDFs[["EXP4"]]
  exp5 = output_netCDFs[["EXP5"]]
  
  # extract_data = rbind(
  #   c("group2",233,265),
  #   c("group2",233,266),
  #   c("group1",399,096),
  #   c("group1",400,096),
  #   c("group1",400,095)
  # )
  # extract_data = data.frame(extract_data,stringsAsFactors = F)
  # extract_data[,2] = as.double(extract_data[,2])
  # extract_data[,3] = as.double(extract_data[,3])
  # extract_data = rbind(
  #   cbind(extract_data,year = 1),
  #   cbind(extract_data,year = 2)
  # )
  # extract_data[,4] = as.double(extract_data[,4])
  # extract_data = cbind(extract_data,matrix(0,nrow(extract_data),12))
  # colnames(extract_data)[-c(1:4)] = to_file_EXP4to5_LU_name
  # for (i in 1:nrow(extract_data)) {
  #   extract_data[i,-c(1:4)] = exp5[extract_data[i,2],extract_data[i,3],,extract_data[i,4]]
  # }
  
  
  ###### export netCDFs
  david_version = "EXP1"
  for (david_version in david_versions) {
    if (david_version %in% c("EXP4","EXP5")) {
      CLASSES = length(to_file_EXP4to5_LU_name)
      CLASSES_NAME = to_file_EXP4to5_LU_name
    } else {
      CLASSES = length(to_file_EXP1to3_LU_name)
      CLASSES_NAME = to_file_EXP1to3_LU_name
    }
    lc_class <- seq(1,CLASSES,1)
    if (david_version %in% c("EXP4","EXP5")) {
      dim_lc_class <- ncdf4::ncdim_def("lc_class", 
                                       paste0("1=cropland_other/2=cropland_2Gbioen/3=grassland/",
                                              "4=forest_unmanaged/5=forest_managed/6=restored/",
                                              "7=other/8=built-up/9=abn_cropland_other/",
                                              "10=abn_cropland_2Gbioen/11=abn_grassland/",
                                              "12=abn_forest_managed"),
                                       lc_class)
    } else {
      dim_lc_class <- ncdf4::ncdim_def("lc_class", 
                                       paste0("1=cropland_other/2=cropland_2Gbioen/3=grassland/",
                                              "4=forest_unmanaged/5=forest_managed/6=restored/",
                                              "7=other/8=built-up"),
                                       lc_class)
    }
    fillvalue <- NA
    LandCover_pixshare_long <- 'share of pixel occupied by various land covers'
    LandCover_pixshare <- ncvar_def(name = "LC_area_share",units = "share of pixel area",longname = LandCover_pixshare_long,
                                    dim = list(dim_lon,dim_lat,dim_lc_class,dim_time), missval=fillvalue,prec = "double",compression = 9)
    PixelArea_Mha_long <- 'total area of the pixel'
    PixelArea_Mha <- ncvar_def(name = "pixel_area",units = "million ha",  longname = PixelArea_Mha_long,
                               dim = list(dim_lon,dim_lat), missval=fillvalue,prec = "double",compression = 9)
    
    LandCover_pixshare_array = output_netCDFs[[david_version]]
    
    save_date_tag = paste0(settings$date_tag,"_",david_version)
    time_date_tag = "26Aug19"
    fsave = paste0(target_path,"BendingTheCurveFT-LCproj-GLOBIOM-",RCP,"-",save_date_tag,"-",time_date_tag,".nc")
    ncnew <- nc_create(fsave, list(LandCover_pixshare,PixelArea_Mha),force_v4=TRUE)
    ncvar_put(ncnew,PixelArea_Mha,area_variable_Mha)
    ncvar_put(ncnew,LandCover_pixshare,LandCover_pixshare_array)
    
    # add attributes
    ncatt_put(ncnew,0,"title",'template for IAM land cover projections')
    ncatt_put(ncnew,0,"institution",'IIASA')
    history <- paste("T. Krisztin, krisztin@iiasa.ac.at", date(), sep=", ")
    ncatt_put(ncnew,0,"history",history)
    
    nc_close(ncnew)
  }
  
}

# 
# #### test protected areas
# 
# # load protected areas
# pa_30by30_cr = read.csv("Export_for_Tamas_30by30_12Feb2020.csv")
# ddata = cbind(x = pa_30by30_cr$x,y = pa_30by30_cr$y,z=pa_30by30_cr$PA2020)
# temp_map = rasterFromXYZ(ddata,crs = p4s)
# temp_map =  extend(temp_map,extent(-180,180,-90,90))
# pa_data <- t(as.matrix(temp_map))
# 
# # load netCDF
# fname = c("P:/globiom/Projects/FABLE/BendingTheCurve/Results/downscaled/BendingTheCurveFT-LCproj-GLOBIOM-RCPref_SSP2_NOBIOD-lowp2cor3dev3_6jun19_mod1_EXP5-26Aug19.nc")
# ncin = nc_open(fname)
# tmp.array = ncvar_get(ncin,"LC_area_share") 
# exp5_2020 = tmp.array[,,,2]
# chck_exp5 =apply(exp5_2020[,,c(1,2,3,8)],c(1,2),sum)
# 
# sum(chck_exp5 > 1 - pa_data,na.rm = T)
# a1 = chck_exp5 - (1 - pa_data)
# range(a1,na.rm=T)
#
# compare baseline data
chck_LU = res_file[,1:16]
chck_LU = chck_LU[chck_LU$Year == "2020",]
chck_CR = aggregate(chck_LU[,-c(1:6)],by = list(Colrow = chck_LU$Colrow),FUN = sum)
chck_CR1 = rowSums(chck_CR[,c("cropland","grassland","SRP","other")]) / rowSums(chck_CR[,-c(1)])
pa_30by30_cr2 = pa_30by30_cr[match(chck_CR$Colrow,pa_30by30_cr$ColRow30),]
sum(chck_CR1 > 1 - pa_30by30_cr2$PA2020)


