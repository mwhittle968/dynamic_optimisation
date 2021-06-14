######################## Data extraction from figures ###########################

library(metaDigitise)

data <- metaDigitise("C:/Users/mw14794/OneDrive - University of Bristol/Documents/PhD/Notes/Model/Data")
data.x <- metaDigitise("C:/Users/mw14794/OneDrive - University of Bristol/Documents/PhD/Notes/Model/Data")

stages <- c("L1", "L2", "L3", "early.pupae", "late.pupae", "teneral.adult", "2week.adult", "4week.adult")
W.density <- c(data$mean[1:8])
day <- c(data.x$mean[9], data.x$mean[11], data.x$mean[13], data.x$mean[15], data.x$mean[17], data.x$mean[19], data.x$mean[21], data.x$mean[23])

Obs.F.data <- data.frame(stages, day, W.density) %>%
  mutate(time.step = day/(9/4))
save(Obs.F.data, file = "Model_3/Obs.F.data.Rdata")

Obs.F.plot <- ggplot(data = F.data, aes(x = time.step, y = W.density))+
  geom_line()
Obs.F.plot