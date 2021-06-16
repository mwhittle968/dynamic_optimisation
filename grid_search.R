######################## Grid search ##########################

# Set working directory
setwd("C:/Users/mw14794/OneDrive - University of Bristol/Documents/PhD/Notes/Model/R_code")

# Packages
library(dplyr)
library(reshape2)
library(utils)
library(pracma)

# Source code
source("Model_4.R")
source("interpolate.R")
load("Obs.F.data.Rdata")

# Create a dataframe for rss values
model.fit <- data.frame(model = character(), model.rss = integer())

I <- c(0.01, 0.05, 0.1, 0.5, 1)
J <- c(1, 2, 6)
K <- c(1, 2, 6)
L <- c(0.2, 1, 1.5, 3)

for (a in 1:5){
  for (b in 1:3){
    for (c in 1:3){
      model.name = paste("4", a, b, c, sep = "")
      dynamic.optimisation.4(model.name, i = I[a], j = J[b], k = K[b], l = L[c])
      load(paste("Model_4/rss.", model.name, ".RData", sep = ""))
      model.fit <- add_row(model.fit, model = model.name, model.rss = rss)
    }
  }
}
model.fit <- arrange(model.fit, model.rss)