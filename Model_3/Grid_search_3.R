######################## Grid search ##########################

# Set working directory
setwd("C:/Users/mw14794/OneDrive - University of Bristol/Documents/PhD/Notes/Model/R_code")

# Packages
library(tidyverse)
library(reshape2)
library(utils)
library(pracma)

# Source code
source("Model_3.R")
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
      model.name = paste("3", a, b, c, sep = "")
      dynamic.optimisation.3(model.name, i = I[a], j = J[b], k = K[b], l = L[c])
      load(paste("Model_3/rss.", model.name, ".RData", sep = ""))
      model.fit <- add_row(model.fit, model = model.name, model.rss = rss)
    }
  }
}

dynamic.optimisation.3("test", i = 0.01, j = 1, k = 1, l = 0.2)
