######################## Data extraction from figures ###########################

library(metaDigitise)

data <- metaDigitise("C:/Users/mw14794/OneDrive - University of Bristol/Documents/PhD/Notes/Model/Data")

stages <- c("L1", "L2", "L3", "early.pupae", "late.pupae", "teneral.adult", "2.week.adult", "4.week.adult")
W.density <- c(data$mean[2], data$mean[4], data$mean[6], data$mean[8], data$mean[10], data$mean[12], data$mean[14], data$mean[16])
day <- c(data$mean[1], data$mean[3], data$mean[5], data$mean[7], data$mean[9], data$mean[11], data$mean[13], data$mean[15])

Obs.F.data <- data.frame(stages, day, W.density) %>%
  mutate(time.step = day/(9/4))
save(Obs.F.data, file = "Obs.F.data.Rdata")

Obs.F.plot <- ggplot(data = Obs.F.data, aes(x = time.step, y = W.density))+
  geom_line()
Obs.F.plot
