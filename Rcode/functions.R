# rm(list = ls())
#****************************************************************************************************
# basic  functionsd
#****************************************************************************************************

# creat function round for bin plot
mround <- function(x,base){ 
  base*round(x/base) 
}

# not in function
'%!in%' <- function(x,y)!('%in%'(x,y))

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

#****************************************************************************************************
# functions for MGRsD data processing
#****************************************************************************************************
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


#****************************************************************************************************
# time for warming and space for time
#****************************************************************************************************
# plot and show time for warm and space for time

plot_tfw_sft <- function(sdata, sdata2) {
  # Q10 vs MAT (space for time appraoch)
  sdata <- subset(sdata, sdata$Q10_all < 10)
  srdb_Q10_MAT <- ggplot(sdata, aes(TAnnual_Del, Q10_all)) +
    geom_point(alpha=0.25) + geom_smooth(method = "lm") +
    labs(x=expression(Annual~temperature~'C'~degree~C~')'), y = expression('Q '[10]))

  p_Q10_Tannual <- ggplot(sdata2, aes(Tm_Annual, Q10SY)) + 
    geom_point(alpha = 0.25) + geom_smooth(method = 'lm') +
    labs(x=expression(Annual~temperature~'('~degree~C~')'), y=expression("Q "[10])) 
  
  # time for warm (by decades)
  # Q10 by climate region
  # Q10_byRegion <- summarySE(sdata2, measurevar="Q10SY", groupvars=c("MiddleClimate"))
  # Q10_byRegion$MiddleClimate <- factor(Q10_byRegion$MiddleClimate, levels = c('E', 'Ds(w)', 'Df', 'Cw', 'Cs', 'Cf', 'B' ))
  # GQ10 <- ggplot (Q10_byRegion, aes(x = MiddleClimate, y = Q10SY))
  #  
  # GQ10 <- GQ10 + geom_bar( aes(order = desc(MiddleClimate)), stat="identity", col="black", 
  #                          fill="white",width=.75) +
  #   labs(x="Climate regions", y=element_blank()) +
  #   geom_errorbar(aes(ymin=Q10SY-se, ymax=Q10SY+se), colour="black", width=0.5)+
  #   coord_cartesian(ylim=c(0, 4.5))+
  #   
  #   annotate("text", x = 7, y = 0.5, label = paste0("D (n=", Q10_byRegion$N[1], ")" ), cex = 3.0, angle = 90, hjust = 0)+
  #   annotate("text", x = 6, y = 0.5, label = paste0("B-D (n=", Q10_byRegion$N[2], ")" ), cex = 3.0, angle = 90, hjust = 0)+
  #   annotate("text", x = 5, y = 0.5, label = paste0("C,D (n=", Q10_byRegion$N[3], ")" ), cex = 3.0, angle = 90, hjust = 0)+
  #   annotate("text", x = 4, y = 0.5, label = paste0("A-D (n=", Q10_byRegion$N[4], ")" ), cex = 3.0, angle = 90, hjust = 0)+
  #   annotate("text", x = 3, y = 0.5, label = paste0("A,B (n=", Q10_byRegion$N[5], ")" ), cex = 3.0, angle = 90, hjust = 0)+
  #   annotate("text", x = 2, y = 0.5, label = paste0("A-C (n=", Q10_byRegion$N[6], ")" ), cex = 3.0, angle = 90, hjust = 0)+
  #   annotate("text", x = 1, y = 0.5, label = paste0("A (n=", Q10_byRegion$N[7],")" ), cex = 3.0, angle = 90, hjust = 0)+
  #   theme(legend.title=element_blank())
  
  # Q10 by decades
  GQ10_dec_srdb <- ggplot(subset(sdata, !is.na(sdata$Decades)), aes(x=Decades, y=Q10_all)) + geom_violin() +
    geom_jitter(position = position_jitter(0.2), col = "gray", alpha = 0.5) +
    geom_boxplot(width = 0.2) +
    labs(y=expression("Q "[10])) 
 
  GQ10_dec_MGRsD <- ggplot(subset(sdata2, !is.na(sdata2$Decades)), aes(x=Decades, y=Q10SY)) + geom_violin() +
    geom_jitter(position = position_jitter(0.2), col = "gray", alpha = 0.5) +
    geom_boxplot(width = 0.2) +
    labs(y=expression("Q "[10])) 
  
  # Q10 relationship with Tm (mean in each climate region), not used
  # Tm.Climate <- summarySE(PT_Del, measurevar="Tm", groupvars=c("MiddleClimate"))
  # Tm.Climate <- subset(Tm.Climate, Tm.Climate$MiddleClimate != "A")
  # 
  # Q10.Tm <- cbind(Tm.Climate, Q10_byRegion)
  # Q10.Tm <- Q10.Tm[c(-8), c(-1,-2)]
  # summary(lm(Q10SY~Tm, data = Q10.Tm))
  # 
  # p_Q10_Tm <- ggplot(Q10.Tm, aes(Tm, Q10SY)) +   
  #   geom_point()+
  #   geom_smooth(method="lm")+
  #   coord_cartesian(ylim=c(0, 4.5))+
  #   labs(x=expression("Air temperature "*"("*degree~C*")")
  #        ,y=element_blank())
  
  plot_grid(srdb_Q10_MAT, GQ10_dec_srdb, p_Q10_Tannual, GQ10_dec_MGRsD
            , ncol = 2, labels = c("(a)", "(b)", "(c)", "(d)")
            , hjust = c(-3.25, -3.25, -2.5, -2.5), vjust = c(2.5, 2.5,2.5, 2.5))
}



# MLR by biome
biome_MLR_MGRsD <- function (sdata) {
  
  var_biome <- sort(unique(sdata$MiddleClimate))
  out <- data.frame()
  
  for (i in 1:length(var_biome)) {
    sub_data <- subset(sdata, sdata$MiddleClimate == var_biome[i])
    m <- lm(RsLog ~ Tm * Decades, data = sub_data)
    print(paste0(i," ********** ", var_biome[i]))
    summary(m) %>% print()
    
    inter_1990s <- summary(m)$coefficients[2,1] %>% round(3)
    n_1990s <- subset(sub_data, Decades == '1990s') %>% nrow()
    Q10_1990s <- ifelse( exp(inter_1990s*10) < 1, 1, exp(inter_1990s*10) ) %>% round(3)
    
    inter_2010s <- (summary(m)$coefficients[4,1] + summary(m)$coefficients[2,1]) %>% round(3)
    n_2010s <- subset(sub_data, Decades == '2010s') %>% nrow()
    Q10_2010s <- exp(inter_2010s*10) %>% round(3)
    p_inter <- summary(m)$coefficients[4,4] %>%  round (3)
    
    out <- rbind(out, data.frame(i, var_biome[i], inter_1990s, inter_2010s, Q10_1990s, Q10_2010s, n_1990s, n_2010s, p_inter) )
  }
  colnames(out) <- c('ID', 'Biome', 'inter_1990s', 'inter_2010s','Q10_1990s', 'Q10_2010s', 'n_1990s', 'n_2010s', 'p_inter')
  return (out)
  
}


biome_MLR_srdb <- function (sdata) {
  
  var_biome <- c("Tropical", "Subtropical", "Temperate", "Mediterranean", "Boreal", "Arctic" )
  out <- data.frame()
  
  for (i in 1:length(var_biome)) {
    sub_data <- subset(sdata, sdata$Biome == var_biome[i])
    m <- lm(log(Rs_annual) ~ TAnnual_Del * Decades, data = sub_data)
    print(paste0(i," ********** ", var_biome[i]))
    summary(m) %>% print()
    
    inter_1990s <- summary(m)$coefficients[2,1] %>% round(3)
    n_1990s <- subset(sub_data, Decades == '1990s') %>% nrow()
    inter_2010s <- (summary(m)$coefficients[4,1] + summary(m)$coefficients[2,1]) %>% round(3)
    n_2010s <- subset(sub_data, Decades == '2010s') %>% nrow()
    p_inter <- summary(m)$coefficients[4,4] %>%  round (3)
    
    out <- rbind(out, data.frame(i, var_biome[i], inter_1990s, n_1990s, inter_2010s, n_2010s, p_inter) )
  }
  colnames(out) <- c('ID', 'Biome', 'inter_1990s', 'n_1990s', 'inter_2010s', 'n_2010s', 'p_inter')
  return (out)
}


# Q10 vs PDSI regression
Q10_pdsi <- function (sdata) {
  sdata <- subset(sdata, sdata$pdsi > -999 & Q10_all < 20 )
  p <- ggplot (sdata, aes(pdsi, Q10_all))
  p <- p + geom_point(alpha = 0.2, size = 1) + geom_smooth(method = 'lm') +
    facet_wrap(~Biome, ncol = 2) + 
    labs( y = expression("Q "[10] ), x = expression( pdsi ) ) +
    theme(legend.position = 'top')
  p
}



#****************************************************************************************************
# functions for Time_for_warm analysis based on SRDB_V4

# Compare 1990s and 2010s by Biome
srdb_biome_comparisom <- function (sdata) {
  sdata$Decades <- ifelse (sdata$Year <= 2000, '1990s', '2010s')
  sdata$Biome = factor(sdata$Biome, levels=c("Tropical", "Subtropical", "Temperate", "Mediterranean", "Boreal", "Arctic" ))
  
  p <- ggplot (sdata, aes(TAnnual_Del, log(Rs_annual), colour = Decades))
  p <- p + geom_point(alpha = 0.2, size = 1) + geom_smooth(method = 'lm') +
    facet_wrap(~Biome, ncol = 2) + 
    labs( y = expression(ln(Rs)~(g~C~m^{-2}~d^{-1} ))
          ,x = expression(Air~temperature~(degree~C) ) ) +
    theme(legend.position = 'top')
  print (p)
}

# Q10 of 1990s and 2010s comparison by Biome
srdb_Q10_biome <- function (sdata) {
  sdata <- subset(sdata, sdata$Q10_all < 10 & !is.na(sdata$Decades))
  p_q10 <- ggplot(sdata, aes(Decades, Q10_all)) + geom_violin() +
    geom_jitter( position = position_jitter(0.2), col = 'gray' ) + 
    geom_boxplot(width = 0.2) +
    facet_wrap (~Biome, ncol = 2)
  
  print(p_q10)
}




#****************************************************************************************************
# functions for time for warming comparison

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



#****************************************************************************************************
# part3： functions for longterm soil respiration data: whether we can detect the Q10 change in site scale
#****************************************************************************************************
longtern_test <- function (sdata) {
  dat <- data.frame()
  for (i in 1:ncol(sdata)) {
    ID <- colnames(sdata)[i]
    y <- sdata[, i]
    n <- y %>% as.data.frame() %>% na.omit() %>% nrow()
    mk <- MannKendall(na.exclude(y))
    summary(mk)
    tau <- mk$tau
    p <- mk$sl
    dat <- rbind( dat, data.frame(ID, tau, p, n) )
  }
  
  # print (dat)
  # plot 
  x_lab1 <- paste0("n=", nrow(dat))
  plot_tau <- ggplot(dat, aes(x = x_lab1, y = tau)) + geom_violin() + 
    geom_jitter(position = position_jitter(0.2), col = "gray") + geom_boxplot(width = 0.15) +
    labs(x = "Measured years >=5", y="Mann-Kendall test tau")
  
  # n <= 5 years removed
  x_lab2 <- paste0("n=", nrow(subset(dat, n>5)))
  plot_tau_remove5 <- ggplot(subset(dat, n>5), aes(x = x_lab2, y = tau)) + geom_violin() + 
    geom_jitter(position = position_jitter(0.2), col = "gray") + geom_boxplot(width = 0.15) +
    labs(x = "Measured years >=6", y="Mann-Kendall test tau")
  
  # n <= 6 years removed
  x_lab3 <- paste0("n=", nrow(subset(dat, n>6)))
  plot_tau_remove6 <- ggplot(subset(dat, n>6), aes(x = x_lab3, y = tau)) + geom_violin() + 
    geom_jitter(position = position_jitter(0.2), col = "gray") + geom_boxplot(width = 0.15) +
    labs(x = "Measured years >=7", y="Mann-Kendall test tau")
  
  plot_n <- ggplot(dat, aes(n)) + geom_histogram (bin = 5, xlab = "years", fill = "gray", col = "black") +
    labs(x = "Measured years", y="Frequency (n)")
  
  plot_grid (plot_tau, plot_tau_remove5, plot_tau_remove6, plot_n)
}


#****************************************************************************************************
# part4： functions for cosore data
#****************************************************************************************************
# For most analyses we want to extract one or more of these pieces and combine them--for example, to get a single table of contributors. There are [various](https://cran.r-project.org/package=purrr) [packages](https://cran.r-project.org/package=rlist) for dealing with nested lists, but we can also write our own short extractor function:
  
csr_table <- function(cosore, table_name) {
  
  extract <- function(x, table_name) {
    if(is.null(x[[table_name]])) { return(NULL) }
    # Add an identifier field so we can track things as tables get combined
    x[[table_name]]$CSR_DATASET <- x$description$CSR_DATASET
    x[[table_name]]
  }
  
  dplyr::bind_rows(lapply(cosore, extract, table_name = table_name))
}


# each port / day simulate a Q10 for day and night time
# some site (e.g., i=3 in cosore), measure interval is more than 5 hour, so if in day time scale, obs < 5
DN_Q10 <- function (sdata, var_dataset) {
  outputs <- data.frame()
  var_year <- unique(sdata$Year) %>% sort ()
  
  for (i in 1:length (var_year)) {
    sub_year <- sdata[sdata$Year == var_year[i], ]
    var_month <- unique(sub_year$Month) %>% sort ()
    for (j in 1:length(var_month)) {
      sub_month <- sub_year[sub_year$Month == var_month[j], ]
      
      # only 2-3 obs if using daily timescale
      # changed to (day of 1-10, 11-20, >21)
      var_day <- unique (sub_month$Day_range) %>% sort()
      for (d in 1:length(var_day)) {
        sub_day <- subset(sub_month, Day_range == var_day[d] ) 
        var_port <- unique(sub_day$CSR_PORT) %>% sort()
        
        for (k in 1:length(var_port)) {
          sub_port <- sub_day[sub_day$CSR_PORT == var_port[k], ]
          obs_day <- subset(sub_port, sub_port$DayNight == 'Day') %>% nrow()
          obs_night <- subset(sub_port, sub_port$DayNight == 'Night') %>% nrow()
          day_TS <- mean(sub_port[sub_port$DayNight == "Day", ]$CSR_T5, na.rm = T)
          night_TS <- mean(sub_port[sub_port$DayNight == "Night", ]$CSR_T5, na.rm = T)
          
          # Calculate day time Q10 and error handle
          day_nls <- try( nls(CSR_FLUX ~ R  * exp(Q * CSR_T5), nls.control(maxiter=1000)
                              , data = subset(sub_port, sub_port$DayNight == 'Day'), start = list(R = 0.5, Q = 0.035), trace = TRUE), TRUE )
          
          if(isTRUE(class(day_nls)=="try-error" | obs_day < 6 )) { 
            R = NA
            p_R = NA
            Q = NA
            p_Q = NA
            Q10_day = NA
          } else { 
            sum_nls <- summary(day_nls)
            R = sum_nls$coef[1, c(1)]
            p_R = round(sum_nls$coef[1, c(4)], 4)
            Q = sum_nls$coef[2, c(1)]
            p_Q = round(sum_nls$coef[2, c(4)], 4)
            Q10_day = ifelse(p_Q < 0.1, exp(Q*10), NA )
            Q10_day = ifelse(Q10_day <= 10 & Q10_day >= 1, Q10_day, NA )
          } 
          
          # Calculate night time Q10 and error handle
          night_nls <- try( nls(CSR_FLUX ~ R  * exp(Q * CSR_T5), nls.control(maxiter=1000)
                                , data = subset(sub_port, sub_port$DayNight == 'Night'), start = list(R = 0.5, Q = 0.035), trace = TRUE), TRUE )
          if( isTRUE(class(night_nls)=="try-error" | obs_night < 6 ) ){
            R_night = NA
            p_R_night = NA
            Q_night <- NA
            p_Q_night <- NA
            Q10_night <- NA
          } else {
            sum_nls_night <- summary(night_nls)
            R_night <- sum_nls_night$coef[1, 1]
            p_R_night = round(sum_nls_night$coef[1, c(4)], 4)
            Q_night <- sum_nls_night$coef[2, 1]
            p_Q_night <- sum_nls_night$coef[2, 4] %>% round(4)
            Q10_night <- ifelse(p_Q_night < 0.1, exp(Q_night*10), NA)  
            Q10_night <- ifelse(Q10_night <= 10 & Q10_night >= 1, Q10_night, NA )
          }
          
          outputs <- rbind(outputs
                           ,data.frame(var_dataset, var_year[i], var_month[j], var_day[d], var_port[k], day_TS, night_TS
                                       , R, p_R, Q, p_Q, Q10_day, obs_day
                                       , Q_night, p_Q_night, Q10_night, obs_night ) )
          
          print(paste0('***** year ', i, ': ' ,var_year[i], '***** month ', j, ": " ,var_month[j], "***** day_range ", var_day[d], ": ", "***** port ", k, ": ", var_port[k] ))
        }
      }
    }
  }
  return (outputs)
}


# some dataset [e.g., i=1 in cosore] measure intervel is every half hour
DN_Q10_day <- function (sdata, var_dataset) {
  outputs <- data.frame()
  var_year <- unique(sdata$Year) %>% sort ()
  
  for (i in 1:length (var_year)) {
    sub_year <- sdata[sdata$Year == var_year[i], ]
    var_month <- unique(sub_year$Month) %>% sort ()
    for (j in 1:length(var_month)) {
      sub_month <- sub_year[sub_year$Month == var_month[j], ]
      
      # only 2-3 obs if using daily timescale
      # changed to (day of 1-10, 11-20, >21)
      var_day <- unique (sub_month$Day) %>% sort()
      for (d in 1:length(var_day)) {
        sub_day <- subset(sub_month, Day == var_day[d] ) 
        var_port <- unique(sub_day$CSR_PORT) %>% sort()
        
        for (k in 1:length(var_port)) {
          sub_port <- sub_day[sub_day$CSR_PORT == var_port[k], ]
          obs_day <- subset(sub_port, sub_port$DayNight == 'Day') %>% nrow()
          obs_night <- subset(sub_port, sub_port$DayNight == 'Night') %>% nrow()
          day_TS <- mean(sub_port[sub_port$DayNight == "Day", ]$CSR_T5, na.rm = T)
          night_TS <- mean(sub_port[sub_port$DayNight == "Night", ]$CSR_T5, na.rm = T)
          
          # Calculate day time Q10 and error handle
          day_nls <- try( nls(CSR_FLUX ~ R  * exp(Q * CSR_T5), nls.control(maxiter=1000)
                              , data = subset(sub_port, sub_port$DayNight == 'Day'), start = list(R = 0.5, Q = 0.035), trace = TRUE), TRUE )
          
          if(isTRUE(class(day_nls)=="try-error" | obs_day < 6 )) { 
            R = NA
            p_R = NA
            Q = NA
            p_Q = NA
            Q10_day = NA
          } else { 
            sum_nls <- summary(day_nls)
            R = sum_nls$coef[1, c(1)]
            p_R = round(sum_nls$coef[1, c(4)], 4)
            Q = sum_nls$coef[2, c(1)]
            p_Q = round(sum_nls$coef[2, c(4)], 4)
            Q10_day = ifelse(p_Q < 0.1, exp(Q*10), NA )
            Q10_day = ifelse(Q10_day <= 10 & Q10_day >= 1, Q10_day, NA )
          } 
          
          # Calculate night time Q10 and error handle
          night_nls <- try( nls(CSR_FLUX ~ R  * exp(Q * CSR_T5), nls.control(maxiter=1000)
                                , data = subset(sub_port, sub_port$DayNight == 'Night'), start = list(R = 0.5, Q = 0.035), trace = TRUE), TRUE )
          if( isTRUE(class(night_nls)=="try-error" | obs_night < 6 ) ){
            R_night = NA
            p_R_night = NA
            Q_night <- NA
            p_Q_night <- NA
            Q10_night <- NA
          } else {
            sum_nls_night <- summary(night_nls)
            R_night <- sum_nls_night$coef[1, 1]
            p_R_night = round(sum_nls_night$coef[1, c(4)], 4)
            Q_night <- sum_nls_night$coef[2, 1]
            p_Q_night <- sum_nls_night$coef[2, 4] %>% round(4)
            Q10_night <- ifelse(p_Q_night < 0.1, exp(Q_night*10), NA)  
            Q10_night <- ifelse(Q10_night <= 10 & Q10_night >= 1, Q10_night, NA )
          }
          
          outputs <- rbind(outputs
                           ,data.frame(var_dataset, var_year[i], var_month[j], var_day[d], var_port[k], day_TS, night_TS
                                       , R, p_R, Q, p_Q, Q10_day, obs_day
                                       , Q_night, p_Q_night, Q10_night, obs_night ) )
          
          print(paste0('***** year ', i, ': ' ,var_year[i], '***** month ', j, ": " ,var_month[j], "***** day_range ", var_day[d], ": ", "***** port ", k, ": ", var_port[k] ))
        }
      }
    }
  }
  return (outputs)
}



DN_Q10_month <- function (sdata, var_dataset) {
  outputs <- data.frame()
  var_year <- unique(sdata$Year) %>% sort ()
  
  for (i in 1:length (var_year)) {
    sub_year <- sdata[sdata$Year == var_year[i], ]
    var_month <- unique(sub_year$Month) %>% sort ()
    for (j in 1:length(var_month)) {
      sub_month <- sub_year[sub_year$Month == var_month[j], ]
      
      # only 2-3 obs if using daily timescale
      # changed to monthly time scale
     
      var_port <- unique(sub_month$CSR_PORT) %>% sort()
      
      for (k in 1:length(var_port)) {
        sub_port <- sub_month[sub_month$CSR_PORT == var_port[k], ]
        obs_day <- subset(sub_port, sub_port$DayNight == 'Day') %>% nrow()
        obs_night <- subset(sub_port, sub_port$DayNight == 'Night') %>% nrow()
        day_TS <- mean(sub_port[sub_port$DayNight == "Day", ]$CSR_T5, na.rm = T)
        night_TS <- mean(sub_port[sub_port$DayNight == "Night", ]$CSR_T5, na.rm = T)
        
        # Calculate day time Q10 and error handle
        day_nls <- try( nls(CSR_FLUX ~ R  * exp(Q * CSR_T5), nls.control(maxiter=1000)
                            , data = subset(sub_port, sub_port$DayNight == 'Day'), start = list(R = 0.5, Q = 0.035), trace = TRUE), TRUE )
        
        if(isTRUE(class(day_nls)=="try-error" | obs_day < 6 )) { 
          R = NA
          p_R = NA
          Q = NA
          p_Q = NA
          Q10_day = NA
        } else { 
          sum_nls <- summary(day_nls)
          R = sum_nls$coef[1, c(1)]
          p_R = round(sum_nls$coef[1, c(4)], 4)
          Q = sum_nls$coef[2, c(1)]
          p_Q = round(sum_nls$coef[2, c(4)], 4)
          Q10_day = ifelse(p_Q < 0.1, exp(Q*10), NA )
          Q10_day = ifelse(Q10_day <= 10 & Q10_day >= 1, Q10_day, NA )
        } 
        
        # Calculate night time Q10 and error handle
        night_nls <- try( nls(CSR_FLUX ~ R  * exp(Q * CSR_T5), nls.control(maxiter=1000)
                              , data = subset(sub_port, sub_port$DayNight == 'Night'), start = list(R = 0.5, Q = 0.035), trace = TRUE), TRUE )
        if( isTRUE(class(night_nls)=="try-error" | obs_night < 6 ) ){
          R_night = NA
          p_R_night = NA
          Q_night <- NA
          p_Q_night <- NA
          Q10_night <- NA
        } else {
          sum_nls_night <- summary(night_nls)
          R_night <- sum_nls_night$coef[1, 1]
          p_R_night = round(sum_nls_night$coef[1, c(4)], 4)
          Q_night <- sum_nls_night$coef[2, 1]
          p_Q_night <- sum_nls_night$coef[2, 4] %>% round(4)
          Q10_night <- ifelse(p_Q_night < 0.1, exp(Q_night*10), NA)  
          Q10_night <- ifelse(Q10_night <= 10 & Q10_night >= 1, Q10_night, NA )
        }
        
        outputs <- rbind(outputs
                         ,data.frame(var_dataset, var_year[i], var_month[j], sub_port$Day[1], var_port[k], day_TS, night_TS
                                     , R, p_R, Q, p_Q, Q10_day, obs_day
                                     , R_night, p_R_night, Q_night, p_Q_night, Q10_night, obs_night ) )
        
        print(paste0('***** year ', i, ': ' ,var_year[i], '***** month ', j, ": " ,var_month[j], "***** port ", k, ": ", var_port[k] ))
      }
    }
  }
  return (outputs)
}


# plot function
Q10_DN_Comp <- function (sdata) {
  
  title_day <- paste0("(a) Day time Q10, n=", nrow( subset(sdata, !is.na(sdata$Q10_day)) )
                      , ", mean=", mean(sdata$Q10_day, na.rm = TRUE) %>% round(2)
                      , ", MT5=", mean(sdata$day_TS, na.rm=TRUE) %>% round(2)) 
  
  title_night <- paste0("(b) Night time Q10, n=", nrow( subset(sdata, !is.na(sdata$Q10_night)) )
                        , ", mean=", mean(sdata$Q10_night, na.rm = TRUE) %>% round(2)
                        , ", MT5=", mean(sdata$night_TS, na.rm=TRUE) %>% round(2)) 
  
  # plot day time Q10
  plot_Q10_day <- ggplot(sdata, aes(Q10_day) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(title_day) +
    xlab(expression(Q[10])) +
    theme(legend.position = c(0.75, 0.65), legend.background = element_rect(colour = 'transparent', fill = alpha('white', 0), size = 0.75) ) +
    xlim (0, 8) +
    geom_vline(xintercept = mean(sdata$Q10_day, na.rm = TRUE), linetype="dotted", color = "blue", size=1.25)
  # geom_vline(xintercept = mean(sdata$Q10_Rs, na.rm = TRUE), col = "blue", linetype = "dotted", size = 1.25 )
  
  # plot night time Q10
  plot_Q10_night <- ggplot(sdata, aes(Q10_night) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(title_night) +
    xlab(expression(Q[10])) +
    theme(legend.position = c(0.75, 0.65), legend.background = element_rect(colour = 'transparent', fill = alpha('white', 0), size = 0.75) ) +
    xlim (0, 8) +
    geom_vline(xintercept = mean(sdata$Q10_night, na.rm = TRUE), linetype="dotted", color = "blue", size=1.25)
  # geom_vline(xintercept = mean(sdata$Q10_Rs, na.rm = TRUE), col = "blue", linetype = "dotted", size = 1.25 )
  
  plot_grid(plot_Q10_day, plot_Q10_night, ncol = 2)
}

