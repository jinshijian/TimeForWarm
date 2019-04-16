# rm(list = ls())
#****************************************************************************************************
# basic  functionsd
# creat function round for bin plot
mround <- function(x,base){ 
  base*round(x/base) 
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     median = median   (xx[[col]], na.rm=na.rm),
                     #do.call("rbind", tapply(xx[[col]], measurevar, quantile, c(0.25, 0.5, 0.75)))
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#****************************************************************************************************
# functions for MGRsD

#1 Study_ID = 40, 1965 only has two records



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
      Q10 <- exp(Q*(Tm+10)) / exp(Q*(Tm))
      obs <- nrow(MGRsD_ID)
      
      results <- rbind(results 
                       ,data.frame(n_ID[i], n_year[j], R, p_R, Q, p_Q, Tm, Q10, obs, length(n_year), MiddleClimate) )
      
      print(paste0('*****', i,':' ,n_ID[i], '*****',j,":" ,n_year[j]))
    }
  }
  
  colnames(results) <- c( "studyNumber", "studyYear", "R", "p_R", "Q", "p_Q", "Tm_Annual", "Q10SY", "obs", "n_years", "MiddleClimate" )
  return (results)
}


#****************************************************************************************************
# functions for Time_for_warm analysis based on SRDB_V4

filtration <- function (sdata, Q10_type) {
  n_col <- colnames(sdata) == Q10_type
  sdata <- sdata[!is.na(sdata[, n_col]), ]
  sdata <- sdata[!is.na(sdata$Year), ]
  sdata
}


Q10_early_late <- function(sdata, Q10_type, var_title) {
  n_col <- colnames(sdata) == Q10_type
  sdata <- sdata[!is.na(sdata[, n_col]), ]
  sdata <- sdata[!is.na(sdata$Year), ]
  sdata$Q10 <- sdata[, colnames(sdata) == Q10_type]
  sdata <- sdata[sdata$Q10 > 0 & sdata$Q10 < 10, ]
  
  sdata$early_late <- ifelse(sdata$Year > 2000, 'Later', 'Early')
  n_early <- nrow(sdata[sdata$early_late == 'Early',])
  n_late <- nrow(sdata[sdata$early_late == 'Later',])
  
  print(paste0('----------', 'n_early = ', n_early, ' n_late = ', n_late, '----------'))
  print(paste0('**********ANOVA**********', Q10_type))
  anova <- aov(Q10 ~ early_late, data = sdata)
  print (summary(anova))
  
  p <- ggplot(sdata, aes(Q10, color = (Year>2000), fill = (Year>2000)) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(var_title) +
    xlab(expression(Q[10])) +
    theme(legend.position = c(0.75,0.65), legend.background = element_rect(colour = 'transparent', fill = alpha('white', 0), size = 0.75) )
}


R10_early_late <- function (sdata, R10_type, var_title) {
  print(SEPARATOR)
  sdata <- sdata[!is.na(sdata$R10), ]
  sdata <- sdata[!is.na(sdata$Year), ]
  sdata <- sdata[sdata$R10 > 0 & sdata$R10 < 10, ]
  sdata$early_late <- ifelse(sdata$Year > 2000, 'Later', 'Early')
  n_early <- nrow(sdata[sdata$early_late == 'Early',])
  n_late <- nrow(sdata[sdata$early_late == 'Later',])
  anova <- aov(R10 ~ early_late, data = sdata)
  print (summary(anova))
  
  print(paste0('----------', 'n_early = ', n_early, ' n_late = ', n_late, '----------'))
  p <- ggplot(sdata, aes(R10, color = Year>2000, fill = Year>2000) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(var_title) + xlab(expression( R[10] ))
  print(p)
}

