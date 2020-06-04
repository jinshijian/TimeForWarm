
rm(list=ls())
#*******************************************************************************************************
# 1. First test, plot 4 years rolling mean
#*******************************************************************************************************
# install.packages("zoo")
library(zoo)
library("ggplot2")

setwd(dirname(rstudioapi::getSourceEditorContext()$path))
getwd()

longterm <- read.csv("LongTerm.csv", header = T)



