
rm(list=ls())

#*******************************************************************************************************
# test code
#*******************************************************************************************************

source('functions.R')

#*******************************************************************************************************
# create some data
set.seed(993)
x <- 1:300
y <- sin(x/20) + rnorm(300,sd=.1)
y[251:255] <- NA

# Plot the unsmoothed data (gray)
plot(x, y, type="l", col=grey(.5))
# Draw gridlines
grid()

# Smoothed with lag:
# average of current sample and 19 previous samples (red)
f20 <- rep(1/20, 20)
f20
#>  [1] 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.05
#> [18] 0.05 0.05 0.05
y_lag <- filter(y, f20, sides=1)
lines(x, y_lag, col="red")

# Smoothed symmetrically:
# average of current sample, 10 future samples, and 10 past samples (blue)
f21 <- rep(1/21,21)
f21
#>  [1] 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905
#>  [8] 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905
#> [15] 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905 0.04761905
y_sym <- filter(y, f21, sides=2)
lines(x, y_sym, col="blue")




# Make same plots from before, with thicker lines
plot(x, y, type="l", col=grey(.5))
grid()
y_lag <- filter(y, rep(1/20, 20), sides=2)
lines(x, y_lag, col="red", lwd=4)         # Lagged average in red
y_sym <- filter(y, rep(1/21,21), sides=2)
lines(x, y_sym, col="blue", lwd=4)        # Symmetrical average in blue

# Calculate lagged moving average with new method and overplot with green
y_lag_na.rm <- movingAverage(y, 4)
lines(x, y_lag_na.rm, col="red", lwd=2)

# Calculate symmetrical moving average  with new method and overplot with green
y_sym_na.rm <- movingAverage(y, 4, TRUE)
lines(x, y_sym_na.rm, col="green", lwd=2)




