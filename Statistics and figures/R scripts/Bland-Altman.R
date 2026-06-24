# 20260623

setwd("~/Desktop/NN-Assisted-Image-Analysis-for-Quantifying-Intracellular-Trypanosoma-cruzi-Infection/Statistics and figures/Data")

#Libraries
library(tidyverse)

#Load data
auto <- read.csv("fov_results.csv") #Summary stats per image using the automated pipeline
manual <- read.csv("fov_manual.csv") #Summary stats for the manual quantification
total <- rbind(auto, manual)

#Difference in percent infected cells
auto$Difper <- auto$Porcentaje - manual$Porcentaje
auto$avgper <- (auto$Porcentaje + manual$Porcentaje)/2

auto$Difper %>% mean()
auto$Difper %>% sd()
 
#Using 1 SD for now.
bland_altman_percent <- ggplot(auto, aes(x = avgper, y = Difper)) + geom_point() +
  geom_hline(yintercept = 9.4, colour = "blue", linewidth = 1) +
  geom_hline(yintercept = c(27.9, -9.1), colour = "red", linetype = "dashed", linewidth = 1)
bland_altman_percent

#Is it any different if we calculate these parameters directly from the single cell results, as was modeled before?
sc_auto <- read.csv("sc_results.csv")
sc_manual <- read.csv("sc_manual.csv")

#Auto
sc_auto_grouped <- sc_auto %>% 
                  group_by(Cepa, Linea) %>%
                  summarize(
                    infected = sum(Amastigotes > 1),
                    total = n()
                  ) %>%
                  mutate(
                    percent = (infected / total)*100
                  )

# Manual             
sc_manual_grouped <- sc_manual %>% 
  group_by(Cepa, Linea) %>%
  summarize(
    infected = sum(Amastigotes > 1),
    total = n()
  ) %>%
  mutate(
    percent = (infected / total)*100
  )

#Difference
sc_auto_grouped$Dif <- sc_auto_grouped$percent - sc_manual_grouped$percent
sc_auto_grouped$avg <- (sc_auto_grouped$percent + sc_manual_grouped$percent)/2

#The differences are smaller if we evaluate the number of infected cells by group instead of by image 
#(fields with few cells have large but less meaningful differences)
sc_auto_grouped$Dif %>% mean()
sc_auto_grouped$Dif %>% sd()

#Fewer points, so fewer outliers
bland_altman_percent_grouped <- ggplot(sc_auto_grouped, aes(x = avg, y = Dif)) + geom_point() +
  geom_hline(yintercept = 1.5, colour = "blue", linewidth = 1) +
  geom_hline(yintercept = c(9.6, -6.6), colour = "red", linetype = "dashed", linewidth = 1)
bland_altman_percent_grouped

#Percent infected vs number of cells in img
Dif_vs_cells <- ggplot(auto, aes(x = Células, y = Difper)) + geom_point()
Dif_vs_cells

Percent_vs_cells <- ggplot(total, aes(x = Células, y = Porcentaje)) + geom_point() + 
  facet_grid(rows = vars(Method))
Percent_vs_cells
