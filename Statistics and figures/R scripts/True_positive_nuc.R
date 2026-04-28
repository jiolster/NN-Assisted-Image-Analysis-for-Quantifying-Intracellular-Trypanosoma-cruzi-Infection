df <- read.csv("truth_results_Nuc.csv")
df1 <- read.csv("truth_results_Nuc_og.csv")

library(tidyverse)

str(df)

#Superposicion entre mascaras de nucelos y anti-Tc para el modelo entrenado
box <- ggplot(df, aes(x = Linea, y = Superposicion, fill = Cepa )) +
  geom_boxplot() +
  theme_minimal() +
  scale_fill_manual(values = c('#117733', '#CC6677')) +
  theme(panel.spacing = unit(.2, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
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
  ylab("Overlap percentage") + xlab("Cell line") +  labs(fill = "Strain")
box

ggsave("Overlap_nuc.svg", plot = box, width = 14, height = 8)


#Superposicion entre mascaras de nucelos y anti-Tc para el modelo original (cp_sam)
layer_scales(box)$y$get_limits()

box1 <- ggplot(df1, aes(x = Linea, y = Superposicion, fill = Cepa )) +
  geom_boxplot() + coord_cartesian(ylim = c(0, 85.43)) +
  scale_fill_manual(values = c('#117733', '#CC6677')) +
  theme(panel.spacing = unit(.2, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
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
  ylab("Overlap percentage") + xlab("Cell line") +  labs(fill = "Strain")
box1

ggsave("Overlap_nuc_og.svg", plot = box1, width = 14, height = 8)
