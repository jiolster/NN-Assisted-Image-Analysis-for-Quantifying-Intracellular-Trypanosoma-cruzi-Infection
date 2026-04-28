rm(list = ls())

df <- read.csv("truth_results_cured.csv")

library(tidyverse)

str(df)

# Reemplazar AC16 por C2C12
df$Linea <- as.factor(df$Linea)
levels(df$Linea)[levels(df$Linea) == "AC16"] <- "C2C12"
df$Linea <- factor(df$Linea, levels = c( "BeWo", "C2C12", "Caco2",  "HeLa",   "THP1",  "Vero" ))
levels(df$Linea)

#Considerando los campos donde no se detectaron amastigotes y que no habian amastigoets como 100% de coincidencia
box <- ggplot(df, aes(x = Linea, y = Porcentaje, fill = Cepa )) +
  geom_boxplot()
box

#Eliminando esos campos
df1 <- df[which(df$Notas != "None detected"),]

box1 <- ggplot(df1, aes(x = Linea, y = Porcentaje, fill = Cepa )) +
  geom_boxplot() + scale_fill_manual(values = c('#009988', '#BBBBBB')) +
  theme(panel.spacing = unit(.2, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid= element_line(colour = "lightgray"),
        text = element_text(family = "sans",color = "black", face = "plain", size = 40),
        axis.title = element_text(face = 'plain', size = 30, hjust = 0.5, margin = margin(t=0, r=0, b=10, l=0)),
        axis.text = element_text(color = "black", face = "plain", size = 24),
        legend.title.position = "top",
        legend.title = element_text(face = 'plain', size = 30, hjust = 0.5),
        #legend.title = element_blank(),
        legend.text = element_text(face = 'plain', size = 24),
        legend.position = "top",
        legend.direction = "vertical"
  ) +
  ylab("Percent of masks colocalizing") + xlab("Cell line") +  labs(fill = "Strain")
box1

ggsave("True_Tc.svg", plot = box1, width = 10, height = 9)

009988#Considerando esos campos como 0
df2 <- df
df2$Porcentaje[which(df$Notas == "None detected")] <- 0

box2 <- ggplot(df2, aes(x = Linea, y = Porcentaje, fill = Cepa )) +
  geom_boxplot()
box2
