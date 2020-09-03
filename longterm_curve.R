
#****************************************************************************************************
# prepare data
#****************************************************************************************************
# * Each site (at least 5 years continous measurement), Mann-Kendall test was used to detect the trend change
# * tau > 0 means increase trend, tau > 0.184 means significant increase trend
# * tau < 0 means decrease trend, tau < -0.184 means significant decrease trenc
# * the results show that overall, there are no trend (tau values not significantly differ from 0)
mk_results <- longtern_test(longterm, longterm_Tm)
lm_results <- longtern_lm(longterm, longterm_Tm)

lm_results$first_a %>% mean()
lm_results$first_b %>% mean()
lm_results$n %>% mean()
lm_results$first_a_tm %>% mean()
lm_results$first_b_tm %>% mean()

mk_results %>% 
  dplyr::select(tau_tm, tau) %>% 
  tidyr::gather() %>% 
  ggplot(aes(key, value, fill = key, group = key)) +
  geom_violin() 

#****************************************************************************************************
# plot
#****************************************************************************************************

par( mar=c(2, 0.2, 0.2, 0.2)
     , mai=c(0.2, 0.3, 0.0, 0.1)  # by inches, inner margin
     , omi = c(0.3, 0.4, 0.2, 0.1)  # by inches, outer margin
     , mgp = c(0.5, 0.5, 0) # set distance of axis
     , tcl = 0.4
     , cex.axis = 1.0
     , las = 1
     , mfrow=c(2,2) )

# plot temperature long term trend
tibble(x = -13:13,
       y = 0 + 0.075*x) ->
  tm_data

plot(y ~ x,
     main = "",
     xlab = "",
     ylab = "",
     xaxt = "n",
     pch = 16,
     col = "white",
     data = tm_data
)
Axis(side=1, labels=FALSE)

for(i in 1:nrow(lm_results)){
  # add SLR curve
  curve(0 + lm_results$first_b_tm[i] * x, 0-lm_results$n_tm[i]/2, 0+lm_results$n_tm[i]/2,
        col = ifelse(lm_results$p_b_tm[i] > 0.1, "gray", "black"),
        lty = ifelse(lm_results$p_b_tm[i] > 0.1, 3, 1),
        lwd = 1, add = T )
  # add average
  curve(0 + lm_results$first_b_tm %>% mean() * x, -4, 4, col = "red", lty = 1, lwd = 3, add = T )
}

mtext(side = 2, text = expression(T[Air]~anomaly~"("~degree~C~")"), line = 2.25, cex=1.0, outer = F, las = 0)
text(-10.5, 0.75, "(a)", cex = 1.5, adj = 0)

# only including sites with increase temperature trend
#****************************************************************************************************
lm_results %>% 
  filter(first_b_tm > 0) ->
  sub_lm_results

plot(y ~ x,
     main = "",
     xlab = "",
     ylab = "",
     xaxt = "n",
     # yaxt = "n",
     pch = 16,
     col = "white",
     data = tm_data
)
Axis(side=1, labels=FALSE)
# Axis(side=2, at = c(-1,-0.5,0,0.5,1), labels=c(0,0.5,1,1.5,2))

for(i in 1:nrow(sub_lm_results)){
  # add SLR curve
  curve(0 + sub_lm_results$first_b_tm[i] * x, 0-sub_lm_results$n_tm[i]/2, 0+sub_lm_results$n_tm[i]/2,
        col = ifelse(sub_lm_results$p_b_tm[i] > 0.1, "gray", "black"),
        lty = ifelse(sub_lm_results$p_b_tm[i] > 0.1, 3, 1),
        lwd = 1, add = T )
  # add average
  curve(0 + sub_lm_results$first_b_tm %>% mean() * x, -4, 4, col = "red", lty = 1, lwd = 3, add = T )
}

text(-10.5, 0.75, "(b)", cex = 1.5, adj = 0)

#****************************************************************************************************
# plot Rs long term trend
#****************************************************************************************************
tibble(x = -13:13,
       y = 0 + 35*x) ->
  rs_data

plot(y ~ x,
     main = "",
     xlab = "",
     ylab = "",
     xaxt = "n",
     pch = 16,
     col = "white",
     data = rs_data
)

Axis(side=1, at = c(-10,-5,0,5,10), labels=c(0,5,10,15,20))

for(i in 1:nrow(lm_results)){
  # add SLR curve
  curve(0 + lm_results$first_b[i] * x, 0-lm_results$n[i]/2, 0+lm_results$n[i]/2,
        col = ifelse(lm_results$p_b[i] > 0.1, "gray", "black"),
        lty = ifelse(lm_results$p_b[i] > 0.1, 3, 1),
        lwd = 1, add = T )
  # add average
  curve(0 + lm_results$first_b %>% mean() * x, -4, 4, col = "red", lty = 1, lwd = 3, add = T )
}

mtext(side = 1, text = paste0("Experiment length (year)" ), line = 1.75, cex=1, outer = F)
mtext(side = 2, text = expression(R[S]~anomaly~"("~g~C~m^{-2}~day^{-1}~")"), line = 2.25, cex=1.0, outer = F, las = 0)
text(-10.5, 375, "(c)", cex = 1.5, adj = 0)

# only including sites with increase temperature trend
#****************************************************************************************************
# plot Rs long term trend
plot(y ~ x,
     main = "",
     xlab = "",
     ylab = "",
     xaxt = "n",
     # yaxt = "n",
     pch = 16,
     col = "white",
     data = rs_data
)
Axis(side=1, at = c(-10,-5,0,5,10), labels=c(0,5,10,15,20))
# Axis(side=2, at = seq(-400,400,200), labels=c(0,200,400,600,800))

for(i in 1:nrow(sub_lm_results)){
  # add SLR curve
  curve(0 + sub_lm_results$first_b[i] * x, 0-sub_lm_results$n[i]/2, 0+sub_lm_results$n[i]/2,
        col = ifelse(sub_lm_results$p_b[i] > 0.1, "gray", "black"),
        lty = ifelse(sub_lm_results$p_b[i] > 0.1, 3, 1),
        lwd = 1, add = T )
  # add average
  curve(0 + sub_lm_results$first_b %>% mean() * x, -4, 4, col = "red", lty = 1, lwd = 3, add = T )
}

mtext(side = 1, text = paste0("Experiment length (year)" ), line = 1.75, cex=1, outer = F)
text(-10.5, 375, "(d)", cex = 1.5, adj = 0)


