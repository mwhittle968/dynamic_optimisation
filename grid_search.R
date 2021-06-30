######################## Grid search ##########################

# Set working directory
setwd("C:/Users/mw14794/OneDrive - University of Bristol/Documents/PhD/Notes/Model/dynamic_optimisation")

# Source code
source("Model_b.R")

# Create a dataframe for rss values
model.fit <- data.frame(model = character(), model.rss = integer())

I <- c(0.01, 0.05, 0.1, 0.5, 1)
J <- c(0.5, 1, 1.5)
K <- c(0.5, 1, 1.5)

for (a in 1:5){
  for (b in 1:3){
      model.name = paste("5", a, b, sep = "")
      print(model.name)
      dynamic.optimisation("6", model.name, i = I[a], j = J[b], k = K[b])
      load(paste("Model_6/rss.", model.name, ".RData", sep = ""))
      model.fit <- add_row(model.fit, model = model.name, model.rss = rss)
    }
}

model.fit <- arrange(model.fit, model.rss)