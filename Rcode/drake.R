
pkgconfig::set_config("drake::strings_in_dots" = "literals")
library(readxl)
library(tidyr)
library(lubridate)
library(kableExtra)
library(piecewiseSEM)
library(dplyr)   # needs to come after MASS above so select() isn't masked
library(raster)
library(drake)  # 6.1.0
library(ncdf4)
library(cosore)

OUTPUT_DIR		<- "outputs"
DATA_DIR <- 'data'
if(!file.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

#*****************************************************************************************************************
# functions 
#*****************************************************************************************************************
# Source all needed functions
source('Rcode/functions.R')
source('Rcode/functions_for_cosore.R')
source('Rcode/functions_dgrsd_q10.R')
source('Rcode/COSORE_Q10.R')

# load and processing srdb
load_srdb <- function(){
  srdb_v4 <- read.csv('srdbv4/srdbv4.csv', stringsAsFactors=F)
  # srdb_v4$Decades <- ifelse (srdb_v4$Year <= 1999, '1990s', '2010s')
  srdb_v4$Q10_all <- coalesce(srdb_v4$Q10_0_10, srdb_v4$Q10_0_20, srdb_v4$Q10_5_15, srdb_v4$Q10_10_20, srdb_v4$Q10_other1, srdb_v4$Q10_other2)
  return(srdb_v4)
}

# load and processing srdb_v5
load_srdbv5 <- function(){
  srdb_v5 <- read.csv('../srdb/srdb-data.csv', stringsAsFactors=F)
  # srdb_v5$Decades <- ifelse (srdb_v5$Year <= 1999, '1990s', '2010s')
  return(srdb_v5)
}

# get mat from the global climate data (del)
get_del_mat <- function() {
  longterm_tm = read_file('LongTerm_tm.csv')
  file_location <- "~/Documents/PNNL/bigdata/UDel/Global2011T_XLS"
  for (i in 1:nrow(longterm_tm)) {
    start_year <- longterm_tm$StartYear[i]
    target_lat <- longterm_tm$lat[i]
    target_lon <- longterm_tm$lon[i]
    var_cols = c(which(colnames(longterm_tm) == "X1"):which(colnames(longterm_tm) == "X26")) # jth column
    
    for (j in 1:length(var_cols)){
      target_year <- start_year + var_cols[j] - var_cols[1]
      
      target_tm <- longterm_tm[i, var_cols[j]]
      
      if(is.na(target_tm) | is.na(target_lat)) {
        next
      }
      else {
        # find target years temperature data
        file_name <- paste0(file_location, "/", "air_temp.", target_year, ".txt")
        del <- read.table(file_name)
        colnames(del)[1:14] <- c("lon", "lat", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12")
        
        # calculate annual mean temperature
        del %>% 
          mutate(MAT = (M1+M2+M3+M4+M5+M6+M7+M8+M9+M10+M11+M12)/12) ->
          del
        
        # find data for the target
        ilat <- del[which.min(abs(del$lat - target_lat)), ]$lat
        ilon <- del[which.min(abs(del$lon - target_lon)), ]$lon
        
        del %>% 
          filter(lat == ilat & lon == ilon) %>% 
          dplyr::select(MAT) ->
          del_mat
        
        i_del_mat <- del_mat$MAT
        longterm_tm[i, var_cols[j]] <- i_del_mat
      }
      print(paste0("*****", i, "*****", j))
    }
  }
  return(longterm_tm)
}

#*****************************************************************************************************************
# make a drake plan 
#*****************************************************************************************************************
plan = drake_plan(
  # load data
  srdb_v4 = load_srdb(),
  srdb_v5 = load_srdbv5(),
  PT_Del = read_file('GlobalTempPrecipTimeSeries_Del.csv'),
  MGRsD = read.csv(file_in(!!file.path('~/Documents/PNNL/MGRsD/outputs/MGRsD_SRDBV5.csv'))),
  MGRsD_Q10 = mgrsd_cal_q10(),
  # longterm = read.csv(file_in(!!file.path(DATA_DIR,'LongTerm.csv'))),
  longterm = read_xlsx('LongTerm.xlsx', 1),
  longterm_tm = read_xlsx('LongTerm.xlsx', 2),
  longterm_tm_del = get_del_mat(),
  
  # COSORE data and calculate Q10
  cosore_all = get_cosore_data (),
  corsore_q10 = get_cosore_q10 (cosore_all)
)

make(plan)
