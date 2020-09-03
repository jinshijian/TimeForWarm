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

# read csv and xlsx
read_file <- function(x) read.csv(file.path(DATA_DIR, x), comment.char = "#", stringsAsFactors = FALSE)
read_xlsx <- function(x, n) read_excel(file.path(DATA_DIR, x), sheet = n)
# write csv 
writ_file <- function(input, output) write.csv(input, file.path(OUT_DIR, output), row.names = FALSE)

save_agu_plot <- function(fn, plot = last_plot(), ...) {
  if(!dir.exists(plot_dir)) dir.create(plot_dir)
  suppressMessages(ggsave(file.path(plot_dir, fn), plot, ...))
}

#****************************************************************************************************
# functions for MGRsD data processing and Q10 calculation
#****************************************************************************************************
# MLR by biome
# i = 1
# sdata = MGRsD
biome_MLR_MGRsD_all <- function (sdata) {
 
  out <- data.frame()
  
  sdata %>% 
    mutate(Tm_del = round(Tm_del*8)/8) %>% 
    dplyr::select(Tm_del, RsLog, Decades) %>% 
    group_by(Decades, Tm_del) %>% 
    summarise(RsLog = mean(RsLog),
              obs = n()) ->
    sub_data
  
  m <- lm(RsLog ~ Tm_del * Decades,
          # weights = obs,
          data = sub_data)
  
  summary(m) %>% print()
  
  inter_1990s <- summary(m)$coefficients[2,1] %>% round(3)
  n_1990s <- subset(sub_data, Decades == 'Early') %>% nrow()
  Q10_1990s <- exp(inter_1990s*10) %>% round(3)
  
  inter_2010s <- (summary(m)$coefficients[4,1] + summary(m)$coefficients[2,1]) %>% round(3)
  n_2010s <- subset(sub_data, Decades == 'Later') %>% nrow()
  Q10_2010s <- exp(inter_2010s*10) %>% round(3)
  p_inter <- summary(m)$coefficients[4,4] %>%  round (3)
  p_tm <- summary(m)$coefficients[4,2] %>%  round (3)
  R2 <- summary(m)$r.squared %>% round(3)
  
  out <- rbind(out, data.frame(2001, "All_data", inter_1990s, inter_2010s, Q10_1990s, Q10_2010s, n_1990s, n_2010s, p_tm, p_inter, R2) )
  colnames(out) <- c('B_YR', 'Biome', 'inter_CT', 'inter_WM','Q10_CT', 'Q10_WM', 'n_CT', 'n_WM', 'p_tm', 'p_inter', "R2")
  return (out)
}

# MLR by biome with aggregate
# i = 1
# sdata = MGRsD
# yr = 2000
biome_MLR_MGRsD_weight <- function (sdata, yr) {
  
  sdata$Decades <- ifelse (sdata$Meas_Year <= yr, '1990s', '2010s')
  var_biome <- sort(unique(sdata$MiddleClimate))
  out <- data.frame()
  
  for (i in 1:length(var_biome)) {
    sub_data <- subset(sdata, sdata$MiddleClimate == var_biome[i])
    subset(sub_data, Decades == '1990s') -> sub_1990
    sub_data %>%
      filter(Tm_del <= max(sub_1990$Tm_del, na.rm = T) & Tm_del >= min(sub_1990$Tm_del, na.rm = T)) ->
      sub_data
    sub_data %>% 
      mutate(Tm_del = round(Tm_del*8)/8) %>% 
      dplyr::select(Tm_del, RsLog, MiddleClimate, Decades) %>% 
      group_by(MiddleClimate, Decades, Tm_del) %>% 
      summarise(RsLog = mean(RsLog),
                obs = n()) ->
      sub_data
    
    m <- try(lm(RsLog ~ Tm_del * Decades,
                # weights = obs,
                data = sub_data))
    
    # summary(m) %>% print()
    if (is(m, "try-error")) {
      inter_1990s <- NA
      n_1990s <- NA
      Q10_1990s <- NA
      
      inter_2010s <- NA
      n_2010s <- NA
      Q10_2010s <- NA
      p_tm <- NA
      p_inter <- NA
      R2 <- NA  } 
    
    else {
      inter_1990s <- summary(m)$coefficients[2,1] %>% round(3)
      n_1990s <- subset(sub_data, Decades == '1990s') %>% nrow()
      Q10_1990s <- exp(inter_1990s*10) %>% round(3)
      
      inter_2010s <- (summary(m)$coefficients[4,1] + summary(m)$coefficients[2,1]) %>% round(3)
      n_2010s <- subset(sub_data, Decades == '2010s') %>% nrow()
      Q10_2010s <- exp(inter_2010s*10) %>% round(3)
      p_tm <- summary(m)$coefficients[4,2] %>%  round (3)
      p_inter <- summary(m)$coefficients[4,4] %>%  round (3)
      R2 <- summary(m)$r.squared %>% round(3)
    }
    
    out <- rbind(out, data.frame(yr, var_biome[i], inter_1990s, inter_2010s, Q10_1990s, Q10_2010s, n_1990s, n_2010s, p_tm, p_inter, R2) )
    print(paste0(i," ********** ", var_biome[i]))
  }
  colnames(out) <- c('B_YR', 'Biome', 'inter_CT', 'inter_WM','Q10_CT', 'Q10_WM', 'n_CT', 'n_WM', 'p_tm', 'p_inter', "R2")
  return (out)
}

#****************************************************************************************************
# function for SRDB
#****************************************************************************************************
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


#****************************************************************************************************
# functions for longterm soil respiration data: whether we can detect the Q10 change in site scale
#****************************************************************************************************
# sdata = longterm
# i = 10
longtern_test <- function (sdata_rs, sdata_tm) {
  dat <- data.frame()
  for (i in 1:nrow(sdata_rs)) {
    StudyID <- sdata_rs$SRDB_study[i]
    ID <- sdata_rs$SiteID[i]
    
    # mk for annual_rs
    y <- sdata_rs[i, c(which(colnames(sdata_rs) == "X1"): which(colnames(sdata_rs) == "X26"))]
    y %>% 
      as.data.frame() %>%
      tidyr::gather() %>% 
      na.omit() ->
      y
    n <- y %>% nrow()
    mk <- MannKendall(na.exclude(y$value))
    tau <- mk$tau
    p <- mk$sl
    
    # mk for temperature
    y_tm <- sdata_tm[i, c(which(colnames(sdata_tm) == "X1"): which(colnames(sdata_tm) == "X26"))]
    y_tm %>% 
      as.data.frame() %>%
      tidyr::gather() %>% 
      na.omit() ->
      y_tm
    n_tm <- y_tm %>% nrow()
    mk_tm <- MannKendall(na.exclude(y_tm$value))
    tau_tm <- mk_tm$tau
    p_tm <- mk_tm$sl 
    
    # put all results together
    print(paste0("*****", i))
    dat <- rbind( dat, data.frame(StudyID, ID, tau, p, n, tau_tm, p_tm, n_tm) )
  }
  return(dat)
}

# linear regression
# sdata_rs = longterm
# sdata_tm = longterm_tm_del
# i = 10
longtern_lm <- function (sdata_rs, sdata_tm) {
  dat <- data.frame()
  for (i in 1:nrow(sdata_rs)) {
    StudyID <- sdata_rs$SRDB_study[i]
    ID <- sdata_rs$SiteID[i]
    
    # mk for annual_rs
    y <- sdata_rs[i, c(which(colnames(sdata_rs) == "X1"): which(colnames(sdata_rs) == "X26"))]
    y %>% 
      as.data.frame() %>%
      tidyr::gather() %>% 
      na.omit() ->
      y
    y %>% 
      mutate(yr = c(1:nrow(y))) ->
      y
    n <- y %>% nrow()
    # linear regression
    first_lm <- lm(value ~ yr, data = y)
    first_a <- summary(first_lm)$coefficients[1,1] %>% round(3)
    first_b <- summary(first_lm)$coefficients[2,1] %>% round(3)
    p_b <- summary(first_lm)$coefficients[2,4]%>% round(3)
    first_R2 <- summary(first_lm)$r.squared %>% round(3)
    
    # linear model for temperature
    y_tm <- sdata_tm[i, c(which(colnames(sdata_tm) == "X1"): which(colnames(sdata_tm) == "X26"))]
    y_tm %>% 
      as.data.frame() %>%
      tidyr::gather() %>% 
      na.omit() ->
      y_tm
    y_tm %>% 
      mutate(yr = c(1:nrow(y_tm))) ->
      y_tm
    n_tm <- y_tm %>% nrow()
    # linear regression
    first_lm_tm <- lm(value ~ yr, data = y_tm)
    first_a_tm <- summary(first_lm_tm)$coefficients[1,1] %>% round(3)
    first_b_tm <- summary(first_lm_tm)$coefficients[2,1] %>% round(3)
    p_b_tm <- summary(first_lm_tm)$coefficients[2,4]%>% round(3)
    first_R2_tm <- summary(first_lm_tm)$r.squared %>% round(3)
    
    # put all results together
    print(paste0("*****", i))
    dat <- rbind( dat, data.frame(StudyID, ID, first_a, first_b, p_b, first_R2, n,
                                  first_a_tm, first_b_tm, p_b_tm, first_R2_tm, n_tm) )
  }
  return(dat)
}

#****************************************************************************************************
# time for warming and space for time analysis and plot functions
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
  sdata$early_late <- ifelse(sdata$Year > 2001, 'Later', 'Early')
  n_early <- nrow(sdata[sdata$early_late == 'Early',])
  n_late <- nrow(sdata[sdata$early_late == 'Later',])
  anova <- aov(R10 ~ early_late, data = sdata)
  print (summary(anova))
  
  print(paste0('----------', 'n_early = ', n_early, ' n_late = ', n_late, '----------'))
  p <- ggplot(sdata, aes(R10, color = Year>2001, fill = Year>2001) ) + geom_density(stat = 'density', alpha = 0.25 ) +
    ggtitle(var_title) + xlab(expression( R[10] ))
  print(p)
}

#****************************************************************************************************
# functions for Theil-Sen trend analysis
#****************************************************************************************************
fuzz <- function(x, error) {
  x * rnorm(length(x), mean = 1, sd = error)
}

