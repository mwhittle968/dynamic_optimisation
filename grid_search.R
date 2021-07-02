######################## Grid search ##########################

# Set working directory
setwd("C:/Users/mw14794/OneDrive - University of Bristol/Documents/PhD/Notes/Model/dynamic_optimisation")

# Source code
source("Model.R")

# Create a dataframe for rss values
model.fit <- data.frame(model = character(), model.rss = integer())

I <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.5, 2, 5)
J <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.5, 2, 5)
K <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 1.5, 2, 5)

for (a in 1:5){
  for (b in 1:5){
    for (c in 1:5){
    model.name = paste("", a, b, c, sep = "")
    print(model.name)
    dynamic.optimisation("8", model.name, i = I[a], j = J[b], k = K[c])
    load(paste("Model_8/rss.", model.name, ".RData", sep = ""))
    model.fit <- add_row(model.fit, model = model.name, model.rss = rss)
  }
  }
}

model.fit <- arrange(model.fit, model.rss)
