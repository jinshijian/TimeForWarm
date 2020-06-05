
pkgconfig::set_config("drake::strings_in_dots" = "literals")
library(tidyr)
library(lubridate)
library(kableExtra)
library(piecewiseSEM)
source('functions.R')
library(dplyr)   # needs to come after MASS above so select() isn't masked
library(raster)
library(drake)  # 6.1.0

OUTPUT_DIR		<- "outputs"
DATA_DIR <- 'data'
if(!file.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

#*****************************************************************************************************************
# functions 
#*****************************************************************************************************************
# Source all needed functions
source('Rcode/functions.R')
# read csv and xlsx
read_file <- function(x) read.csv(file.path(DATA_DIR, x), comment.char = "#", stringsAsFactors = FALSE)
read_xlsx <- function(x) read_excel(file.path(DATA_DIR, x), sheet = 1)
# write csv 
writ_file <- function(input, output) write.csv(input, file.path(OUT_DIR, output), row.names = FALSE)

# load and processing srdb
load_srdb <- function(){
  srdb_v4 <- read.csv('srdbv4/srdbv4.csv', stringsAsFactors=F)
  srdb_v4$Decades <- ifelse (srdb_v4$Year <= 2000, '1990s', '2010s')
  srdb_v4$Q10_all <- coalesce(srdb_v4$Q10_0_10, srdb_v4$Q10_0_20, srdb_v4$Q10_5_15, srdb_v4$Q10_10_20, srdb_v4$Q10_other1, srdb_v4$Q10_other2)
  return(srdb_v4)
}

# load and processing srdb_v5
load_srdbv5 <- function(){
  srdb_v5 <- read.csv('../srdb/srdb-data.csv', stringsAsFactors=F)
  srdb_v5$Decades <- ifelse (srdb_v5$Year <= 2000, '1990s', '2010s')
  srdb_v5$Q10_all <- coalesce(srdb_v5$Q10_0_10, srdb_v5$Q10_0_20, srdb_v5$Q10_5_15, srdb_v5$Q10_10_20, srdb_v5$Q10_other1, srdb_v5$Q10_other2)
  return(srdb_v5)
}

# load and clean MGRsD data
clean_mgrsd_obs6 <- function(){
  MGRsD_obs_6 <- read_file('MGRsDObsLargerThan6_V3.CSV')
  # some measure year from MGRsD need to be updated
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 1971 & MGRsD_obs_6$Measure_Year==1998 , 'Measure_Year'] <- 1999
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 3619 & MGRsD_obs_6$Measure_Year==2000 , 'Measure_Year'] <- 2001
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 4234 & MGRsD_obs_6$Measure_Year==1998 , 'Measure_Year'] <- 1999
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 4234 & MGRsD_obs_6$Measure_Year==2005 , 'Measure_Year'] <- 2003
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 4234 & MGRsD_obs_6$Measure_Year==2004 , 'Measure_Year'] <- 2003
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 4234 & MGRsD_obs_6$Measure_Year==2001 , 'Measure_Year'] <- 2003
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 4015 , 'Measure_Year'] <- 2002
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 4270 & MGRsD_obs_6$Measure_Year==1998 , 'Measure_Year'] <- 1999
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 4270 & MGRsD_obs_6$Measure_Year==2001 , 'Measure_Year'] <- 2002
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 4270 & MGRsD_obs_6$Measure_Year==2004 , 'Measure_Year'] <- 2003
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 4270 & MGRsD_obs_6$Measure_Year==2005 , 'Measure_Year'] <- 2003
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 4477 , 'Measure_Year'] <- 2002
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 10002 , 'Measure_Year'] <- 1999
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 10035 , 'Measure_Year'] <- 2009
  MGRsD_obs_6[MGRsD_obs_6$studynumber == 10106 , 'Measure_Year'] <- 2010
  # return data
  return(MGRsD_obs_6)
}

# calculate Q10 by site of each year
mgrsd_cal_q10 <- function(){
  MGRsD_obs_6 <- clean_mgrsd_obs6()
  MGRsD_Q10 <- Q10_cal(check_ID_year(MGRsD_obs_6))
  MGRsD_Q10 <- subset(MGRsD_Q10, MGRsD_Q10$Q10SY > 1 & MGRsD_Q10$Q10SY < 10)
  MGRsD_Q10 <- subset(MGRsD_Q10, MGRsD_Q10$MiddleClimate != "A")
  return(MGRsD_Q10)
}


#*****************************************************************************************************************
# make a drake plan 
#*****************************************************************************************************************
plan = drake_plan(
  # load data
  srdb_v4 = load_srdb(),
  srdb_v5 = load_srdbv5(),
  PT_Del = read_file('GlobalTempPrecipTimeSeries_Del.csv'),
  MGRsD = read_file('MGRsDAllData.csv'),
  MGRsD_Q10 = mgrsd_cal_q10(),
  longterm = read.csv(file_in(!!file.path(DATA_DIR,'LongTerm.csv')))
)

make(plan)
