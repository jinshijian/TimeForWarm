
movingAverage <- function(x, n=1, centered=FALSE) {
  
  if (centered) {
    before <- floor  ((n-1)/2)
    after  <- ceiling((n-1)/2)
  } else {
    before <- n-1
    after  <- 0
  }
  
  # Track the sum and count of number of non-NA items
  s     <- rep(0, length(x))
  count <- rep(0, length(x))
  
  # Add the centered data 
  new <- x
  # Add to count list wherever there isn't a 
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new
  
  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new   <- c(rep(NA, i), x[1:(length(x)-i)])
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new   <- c(x[(i+1):length(x)], rep(NA, i))
    
    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new
    
    i <- i+1
  }
  
  # return sum divided by count
  s/count
}

# filtration function
filtration <- function (sdata, Q10_type) {
  n_col <- colnames(sdata) == Q10_type
  sdata <- sdata[!is.na(sdata[, n_col]), ]
  sdata <- sdata[!is.na(sdata$Year), ]
  sdata
}

# handle study year has only few measurements
clean_mgrsd_obs6 <- function(){
  MGRsD_obs_6 <- read_file('extdata/MGRsDObsLargerThan6_V3.CSV')
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

# Calculate Q10 by site and year
check_ID_year <- function (sdata) {
  n_ID <- sort(unique(sdata$studynumber))
  for (i in 1:length(n_ID)) {
    MGRsD_ID <- sdata[sdata$studynumber == n_ID[i], ]
    n_year <- sort(unique(MGRsD_ID$Measure_Year))
    
    for (j in 1:length(n_year) ) {
      # Study_ID = 40, 1965 only has two records
      
      n_0 <- nrow( subset(MGRsD_ID, Measure_Year == n_year[j-1], c('Measure_Year')) )
      n_1 <- nrow( subset(MGRsD_ID, Measure_Year == n_year[j], c('Measure_Year')) )
      n_2 <- nrow( subset(MGRsD_ID, Measure_Year == n_year[j+1], c('Measure_Year')) )
      
      if (n_1 < 7 & n_1 < n_0) {
        sdata[sdata$studynumber == n_ID[i] & sdata$Measure_Year == n_year[j],'adj_Year'] <- 
          sdata[sdata$studynumber == n_ID[i] & sdata$Measure_Year == n_year[j],'Measure_Year'] - 1
      } 
      else if (n_1 < 7 & n_1 < n_2) { 
        sdata[sdata$studynumber == n_ID[i] & sdata$Measure_Year == n_year[j],'adj_Year'] <-
          sdata[sdata$studynumber == n_ID[i] & sdata$Measure_Year == n_year[j],'Measure_Year'] + 1
      } 
      else {  
        sdata[sdata$studynumber == n_ID[i] & sdata$Measure_Year == n_year[j],'adj_Year'] <- 
          sdata[sdata$studynumber == n_ID[i] & sdata$Measure_Year == n_year[j],'Measure_Year'] 
        
        print(paste0('*****', i, '*****', j))
      }
    }
  }
  return (sdata)
}

# Some study only have less than 3 measurements, for example, Study_ID = 40, 1965 only has two records
# Calculate Q10 for each studynumber under each year
Q10_cal <- function (sdata) {
  results <- data.frame()
  n_ID <- sort(unique(sdata$studynumber))
  
  for(i in 1 : length(n_ID) ) {
    
    MGRsD_ID <- sdata[sdata$studynumber == n_ID[i], ]
    n_year <- sort(unique(MGRsD_ID$adj_Year))
    
    for (j in 1:length(n_year) ) {
      MGRsD_Year <- MGRsD_ID[MGRsD_ID$adj_Year == n_year[j], ]
      
      # how to store this error messege in results?
      # tryCatch ({
      #   # Put the code that may generate an error here.
      #   m_nls <- nls(Rs_Norm ~ R  * exp(Q * Tm), nls.control(maxiter=5000)
      #                , data = MGRsD_Year, start = list(R = 1.4, Q = 0.035), trace = TRUE) 
      #   }, error=function(e){
      #     cat("ERROR :",conditionMessage(e), "\n")
      #     }  )
      
      m_nls <- nls(Rs_Norm ~ R  * exp(Q * Tm), nls.control(maxiter=5000)
                   , data = MGRsD_Year, start = list(R = 0.5, Q = 0.035), trace = TRUE) 
      
      sum_nls <- summary(m_nls)
      
      MiddleClimate <- MGRsD_Year$MiddleClimate[1]
      R <- sum_nls$coef[1, c(1)]
      p_R <- round(sum_nls$coef[1, c(4)], 4)
      Q <- sum_nls$coef[2, c(1)]
      p_Q <- round(sum_nls$coef[2, c(4)], 4)
      Tm <- mean(MGRsD_ID$Tm)
      Q10 <- exp(Q*10) 
      obs <- nrow(MGRsD_ID)
      
      results <- rbind(results 
                       ,data.frame(n_ID[i], n_year[j], R, p_R, Q, p_Q, Tm, Q10, obs, length(n_year), MiddleClimate) )
      
      print(paste0('*****', i,':' ,n_ID[i], '*****',j,":" ,n_year[j]))
    }
  }
  
  colnames(results) <- c( "studyNumber", "studyYear", "R", "p_R", "Q", "p_Q", "Tm_Annual", "Q10SY", "obs", "n_years", "MiddleClimate" )
  return (results)
}

# calculate Q10 by site of each year
mgrsd_cal_q10 <- function(){
  MGRsD_obs_6 <- clean_mgrsd_obs6()
  MGRsD_Q10 <- Q10_cal(check_ID_year(MGRsD_obs_6))
  MGRsD_Q10 <- subset(MGRsD_Q10, MGRsD_Q10$Q10SY > 1 & MGRsD_Q10$Q10SY < 10)
  MGRsD_Q10 <- subset(MGRsD_Q10, MGRsD_Q10$MiddleClimate != "A")
  return(MGRsD_Q10)
}
