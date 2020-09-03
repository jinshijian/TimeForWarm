
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

# day and night time Q10 calculation
csr_rh_Q10 <- function(sdata){
  m <- try(lm(Rh_log ~ CSR_TS, data = sdata))
  if (is(m, "try-error")) {
    Q10 <- NA
  } else {
    intercept <- try(summary(m)$coefficients[2,1]) 
    Q10 <- ifelse(is(intercept, "try-error"), NA, exp(intercept*10) %>% round(3)) 
    R2 <- try(round(summary(m)$r.squared))
  }
  return(Q10)
}