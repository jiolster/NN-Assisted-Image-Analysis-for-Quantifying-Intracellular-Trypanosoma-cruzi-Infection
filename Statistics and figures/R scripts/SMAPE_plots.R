ls()
rm(list = ls())
ls()

setwd("/home/joaquin/Desktop/20250602 - Tul y Dm wt 43a1")

#Libraries
library(tidyverse)


#Data
NN <- read.csv("sc_results.csv") #Data obtained from the automated pipeline
yo <- read.csv("sc_manual.csv") #Data obtained manually

#Joining the dataframes
NN$Method <- "Auto"
yo$Method <- "Manual"

names <- colnames(NN)
colnames(yo) <- names

df <- rbind(NN, yo)

str(df)

#Initial comparison
plot <- ggplot(df, aes(x = Linea, y = Amastigotes, fill = Cepa)) + geom_boxplot() + 
  facet_grid(cols = vars(Method))
plot

# Total amastigotes divided by total cells
avg <- aggregate(df$Amastigotes, list(df$Cepa, df$Linea, df$Method), FUN = mean)
amastot <- aggregate(df$Amastigotes, list(df$Cepa, df$Linea, df$Method), FUN = sum)
celtot <- aggregate(df$Amastigotes, list(df$Cepa, df$Linea, df$Method), FUN = max)
dev <- aggregate(df$Amastigotes, list(df$Cepa, df$Linea, df$Method), FUN = sd)

df1 <- merge(avg, dev, by = c("Group.1", "Group.2", "Group.3"))

str(avg)

#Plot
amas_por_cel <- ggplot(df1, aes(x = Group.1 , y = x.x, fill = Group.3)) + geom_col(position = 'dodge', color = 'black', size = 0.1) +
  xlab("") + ylab("Amastigotes por celula") + guides(fill=guide_legend(title="Método")) +
  facet_grid(cols = vars(Group.2)) +
  geom_errorbar(aes(x=Group.1, ymin=x.x-sqrt(x.x), ymax=x.x+sqrt(x.x)), alpha = 0.5, position = 'dodge')
amas_por_cel


#Custom theme
theme_custom <- theme(
  text = element_text(size = 12),
  panel.grid = element_blank(),
  panel.grid.major.y = element_line(colour = "#e3e1e1",
                                    linetype = 2),
  plot.title.position = 'plot',
  legend.position = 'top',
  legend.title = element_blank()
)

theme_set(theme_minimal() + theme_custom) #Sets the theme for all plots



#Total cells for every cell line
cntcelfov <- aggregate(df$Célula, list(df$Cepa, df$Linea, df$Campo), FUN = max) #Total de celulas por campo
cnttotal <- aggregate(cntcelfov$x, list(cntcelfov$Group.1, cntcelfov$Group.2), FUN = sum) #Total de celulas de cada linea

#Total amastigtes for every cell line
amas <- aggregate(df$Amastigotes, list(df$Cepa, df$Linea), FUN = sum) 

amas$rel <- amas$x / cnttotal$x #Da lo mismo que hacer el promedio como antes


#Subset amastigotes per infected cell
dfInf <- subset(df, df$Amastigotes > 0)
#Average amastigotes per infected cell
infavg <- aggregate(dfInf$Amastigotes, list(dfInf$Cepa, dfInf$Linea, dfInf$Method), FUN = mean) 

amas_por_celInf <- ggplot(infavg, aes(x = Group.1 , y = x, fill = Group.3)) + geom_col(position = 'dodge2') +
  xlab("Linea celular") + ylab("Amastigotes por celula") + guides(fill=guide_legend(title="MOI")) + 
  facet_grid(cols = vars(Group.2))+
  geom_errorbar(aes(x=Group.1, ymin=x-sqrt(x), ymax=x+sqrt(x)), alpha = 0.5, position = 'dodge')
amas_por_celInf + ggtitle('Amastigotes por celula infectada')


#Percent infected cells
dfsinf <- subset(df, df$Amastigotes == 0) #Celulas sin infeccion

cnttotalINF <- dfInf %>% count(Cepa, Linea, Method, sort = TRUE) #Celulas infectadas por linea

cnttotalsinINF <- dfsinf %>% count(Cepa, Linea, Method, sort = TRUE) #Numero de celulas no infectadas por linea

Infeccion <- merge(cnttotalINF, cnttotalsinINF, by=c("Cepa","Linea", "Method")) #df con el total de celulas infectadas (n.x) y no infectadas (n.y) para cada cepa y linea

colnames(Infeccion)[which(names(Infeccion) == "n.x")] <- "Celulas Infectadas"
colnames(Infeccion)[which(names(Infeccion) == "n.y")] <- "Celulas No infectadas"

Infeccion$`Celulas Totales` <- Infeccion$`Celulas Infectadas` + Infeccion$`Celulas No infectadas`
Infeccion$Porcentaje <- (Infeccion$`Celulas Infectadas` / Infeccion$`Celulas Totales`) * 100

bar <- ggplot(Infeccion, aes(x = Cepa, y = Porcentaje, fill = Method)) + geom_col(position = 'dodge2') +
  coord_cartesian(ylim=c(0, 100)) + geom_hline(yintercept = 60, linetype  = 'dashed') + facet_grid(cols = vars(Linea))
bar + xlab('Línea celular') + ylab('Porcentaje de células infectadas') + ggtitle('')

##Histogram
histo <- ggplot(dfInf, aes(x = Amastigotes, fill = Method)) + 
  geom_histogram(binwidth = 1, alpha = 0.7, width = 1, color = 'black', position = 'identity') + 
  coord_cartesian(ylim = c(0, 75)) + ylab("Células") + ggtitle("Número de células con x amastigotes") +
  facet_grid(rows = vars(Cepa))
histo

amas <- aggregate(df$Amastigotes, list(df$Method), FUN=sum) 


dens <- ggplot(dfInf, aes(x = x)) +
  geom_histogram( binwidth = 1, aes(x = subset(Amastigotes, Method == "Auto"), y = ..density..), fill="#69b3a2" ) +
  geom_label( aes(x=4.5, y=0.25, label="Auto"), color="#69b3a2") +
  geom_histogram( aes(x = subset(Amastigotes, Method == "Manual"), y = -..density..), fill= "#404080") +
  geom_label( aes(x=4.5, y=-0.25, label="Manual"), color="#404080") +
  xlab("value of x")
dens

subset(dfInf$Amastigotes, dfInf$Method == "Auto")


##Comparar campo por campo

fovNN <- read.csv("fov_results.csv") #Resumen de campos auto
fovNN <- fovNN[-25,] #No se conto DmTHP1 5 a mano
fovNN[is.na(fovNN)] <- 0

names1 <- colnames(fovNN)

fovMan1 <- aggregate(yo$Amastigotes, list(yo$Cepa, yo$Linea, yo$Campo), FUN = sum) #Amas totales por campo a mano
fovMan1 <- fovMan1[order(fovMan1$Group.1, fovMan1$Group.2, fovMan1$Group.3),]

fovMan2 <- aggregate(yo$Amastigotes, list(yo$Cepa, yo$Linea, yo$Campo), FUN = max) #Celulas totales por campo
fovMan2 <- fovMan2[order(fovMan2$Group.1, fovMan2$Group.2, fovMan2$Group.3),]

#Porcentaje de celulas infectadas
dfInf  <- subset(yo, yo$Amastigotes > 0) #Celulas con infeccion
dfsinf <- subset(yo, yo$Amastigotes == 0) #Celulas sin infeccion

cnttotalINF <- dfInf %>% count(Cepa, Linea, Campo, sort = TRUE) #Celulas infectadas por linea
cnttotalsinINF <- dfsinf %>% count(Cepa, Linea, Campo, sort = TRUE) #Numero de celulas no infectadas por linea
Infeccion <- merge(cnttotalINF, cnttotalsinINF, by=c("Cepa","Linea", "Campo"), all = TRUE) #df con el total de celulas infectadas (n.x) y no infectadas (n.y) para cada cepa y linea
Infeccion[is.na(Infeccion)] <- 0
colnames(Infeccion)[which(names(Infeccion) == "n.x")] <- "Celulas Infectadas"
colnames(Infeccion)[which(names(Infeccion) == "n.y")] <- "Celulas No infectadas"
Infeccion$`Celulas Totales` <- Infeccion$`Celulas Infectadas` + Infeccion$`Celulas No infectadas`
Infeccion$Porcentaje <- (Infeccion$`Celulas Infectadas` / Infeccion$`Celulas Totales`) * 100


fovMan <- Infeccion[order(Infeccion$Cepa, Infeccion$Linea, Infeccion$Campo),]
fovMan$Amastigotes <- fovMan1$x
fovMan$Células <- fovMan2$x
fovMan$Amas.por.célula <- fovMan$Amastigotes / fovMan$`Celulas Totales`
fovMan$Amas.por.célula.infectada <- fovMan$Amastigotes / fovMan$`Celulas Infectadas`

fovManual <- fovMan[,-4]
fovManual <- fovManual[,-4]
fovManual <- fovManual[,-7]

fovManual <- fovManual[,c(1, 2, 3, 6, 4, 7, 5 , 8)]
#write.csv(fovManual, file="fov_manual.csv")

names1 <- colnames(fovManual)
colnames(fovNN) <- names1

fovManual$Method <- "Manual"
fovNN$Method <- "Auto"

#Diferencia entre auto y manual para cada campo
dif <- fovManual[,c(1,2,3)]

#Symmetric mean absolute percentage error
dif$Amastigotes <- (fovNN$Amastigotes - fovManual$Amastigotes) / ((fovManual$Amastigotes + fovNN$Amastigotes)/2)
dif$`Celulas Totales` <- (fovNN$`Celulas Totales` - fovManual$`Celulas Totales`) / ((fovNN$`Celulas Totales` + fovManual$`Celulas Totales`)/2)
dif$Amas.por.célula <- (fovNN$Amas.por.célula - fovManual$Amas.por.célula) / ((fovManual$Amas.por.célula + fovNN$Amas.por.célula)/2)
dif$Porcentaje <- (fovNN$Porcentaje - fovManual$Porcentaje) / ((fovManual$Porcentaje + fovNN$Porcentaje)/2) 
dif$Amas.por.célula.infectada <- (fovNN$Amas.por.célula.infectada - fovManual$Amas.por.célula.infectada) / ((fovManual$Amas.por.célula.infectada + fovNN$Amas.por.célula.infectada)/2) 

dif[is.na(dif)] <- 0

#SMAPE plots
boxamas <- ggplot(dif, aes(x = Linea, y = Amastigotes, fill = Cepa)) + geom_boxplot() + 
  geom_hline(yintercept= c(0.1, -0.1), linetype = 3) + ylab("SAMPE") + 
  ggtitle("Difference in total amastigotes") +
  coord_cartesian(ylim = c(-2,2)) + scale_fill_manual(values = c('#009988', '#BBBBBB')) +
  theme(panel.spacing = unit(.2, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "lightgray"),
        text = element_text(family = "sans",color = "black", face = "plain", size = 24),
        axis.title = element_text(face = 'plain', size = 30, hjust = 0.5, margin = margin(t=0, r=0, b=10, l=0)),
        axis.text = element_text(color = "black", face = "plain", size = 24),
        legend.title.position = "top",
        legend.title = element_text(face = 'plain', size = 30, hjust = 0.5),
        #legend.title = element_blank(),
        legend.text = element_text(face = 'plain', size = 24),
        legend.position = "top",
        legend.direction = "horizontal"
  ) +  labs(fill = "Strain") + xlab("Cell line") 
boxamas

ggsave("SAMPE_Tc.svg", plot = boxamas, width = 14, height = 8)

boxcel <- ggplot(dif, aes(x = Linea, y = `Celulas Totales`, fill = Cepa)) + geom_boxplot() + 
  geom_hline(yintercept= c(0.1, -0.1), linetype = 3) + ylab("SAMPE") + ggtitle("Difference in total host cells") +
  coord_cartesian(ylim = c(-2,2)) + scale_fill_manual(values = c('#009988', '#BBBBBB')) +
  theme(panel.spacing = unit(.2, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "lightgray"),
        text = element_text(family = "sans",color = "black", face = "plain", size = 24),
        axis.title = element_text(face = 'plain', size = 30, hjust = 0.5, margin = margin(t=0, r=0, b=10, l=0)),
        axis.text = element_text(color = "black", face = "plain", size = 24),
        legend.title.position = "top",
        legend.title = element_text(face = 'plain', size = 30, hjust = 0.5),
        #legend.title = element_blank(),
        legend.text = element_text(face = 'plain', size = 24),
        legend.position = "top",
        legend.direction = "horizontal"
  ) +  labs(fill = "Strain") + xlab("Cell line") 
boxcel

ggsave("SAMPE_Nuc.svg", plot = boxcel, width = 14, height = 8)

#Tukey mean difference plot
scatters <-  fovManual[,c(1,2,3)]
scatters$DifAmas <- fovNN$Amastigotes - fovManual$Amastigotes
scatters$AvgAmas <- (fovNN$Amastigotes + fovManual$Amastigotes)/2

meaAmas <- mean(scatters$DifAmas)
sdAmas <- sd(scatters$DifAmas)

scatterama <- ggplot(scatters, aes(x = AvgAmas, y = DifAmas, colour = Linea, shape = Cepa)) + geom_point(size = 3)+
  geom_hline(yintercept = c(-3.44), color = 'steelblue') + geom_hline(yintercept = c(18.67, -25.55), linetype = 3) +
  coord_cartesian(ylim = c(-50,50))
scatterama

scatters$DifCel <- fovNN$`Celulas Totales` - fovManual$`Celulas Totales`
scatters$AvgCel <- (fovNN$`Celulas Totales` + fovManual$`Celulas Totales`)/2

meaCel <- mean(scatters$DifAmas)
sdCel <- sd(scatters$DifAmas)

scatterCel<- ggplot(scatters, aes(x = AvgCel, y = DifCel, colour = Linea, shape = Cepa)) + geom_point(size = 3)+
  geom_hline(yintercept = c(-0.34), color = 'steelblue') + geom_hline(yintercept = c(1.84, -2.52), linetype = 3) +
  coord_cartesian(ylim = c(-5,5))
scatterCel

#SMAPE INSpect

ins <- read.csv("fov_INsPECT.csv")
ins <- ins[-13,]
#Diferencia entre auto y manual para cada campo
difins <- fovManual[,c(1,2,3)]

#Symmetric mean absolute percentage error
difins$Amastigotes <- (ins$Amastigotes - fovManual$Amastigotes) / ((fovManual$Amastigotes + ins$Amastigotes)/2)
difins$`Celulas Totales` <- (ins$Células - fovManual$`Celulas Totales`) / ((ins$Células + fovManual$`Celulas Totales`)/2)
difins$Amas.por.célula <- (ins$Amas.por.célula - fovManual$Amas.por.célula) / ((fovManual$Amas.por.célula + ins$Amas.por.célula)/2)
difins$Porcentaje <- (ins$Porcentaje - fovManual$Porcentaje) / ((fovManual$Porcentaje + ins$Porcentaje)/2) 
difins$Amas.por.célula.infectada <- (ins$Amas.por.célula.infectada - fovManual$Amas.por.célula.infectada) / ((fovManual$Amas.por.célula.infectada + ins$Amas.por.célula.infectada)/2) 

difins[is.na(difins)] <- 0

#SMAPE plots ins
boxamasins <- ggplot(difins, aes(x = Linea, y = Amastigotes, fill = Cepa)) + geom_boxplot() + 
  geom_hline(yintercept= c(0.1, -0.1), linetype = 3) + ylab("SAMPE") + 
  ggtitle("Difference in total amastigotes (Morphological method)") +
  coord_cartesian(ylim = c(-2,2)) + scale_fill_manual(values = c('#009988', '#BBBBBB')) +
  theme(panel.spacing = unit(.2, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "lightgray"),
        text = element_text(family = "sans",color = "black", face = "plain", size = 24),
        axis.title = element_text(face = 'plain', size = 30, hjust = 0.5, margin = margin(t=0, r=0, b=10, l=0)),
        axis.text = element_text(color = "black", face = "plain", size = 24),
        legend.title.position = "top",
        legend.title = element_text(face = 'plain', size = 30, hjust = 0.5),
        #legend.title = element_blank(),
        legend.text = element_text(face = 'plain', size = 24),
        legend.position = "top",
        legend.direction = "horizontal"
  ) +  labs(fill = "Strain") + xlab("Cell line") 
boxamasins

ggsave("SAMPE_Tc_ins.svg", plot = boxamasins, width = 14, height = 8)

boxcelins <- ggplot(difins, aes(x = Linea, y = `Celulas Totales`, fill = Cepa)) + geom_boxplot() + 
  geom_hline(yintercept= c(0.1, -0.1), linetype = 3) + ylab("Error") + ylab("SAMPE") + 
  ggtitle("Relative difference in total host cells (Morphological mehtod)") +
  coord_cartesian(ylim = c(-2,2)) + scale_fill_manual(values = c('#009988', '#BBBBBB')) +
  theme(panel.spacing = unit(.2, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "lightgray"),
        text = element_text(family = "sans",color = "black", face = "plain", size = 24),
        axis.title = element_text(face = 'plain', size = 30, hjust = 0.5, margin = margin(t=0, r=0, b=10, l=0)),
        axis.text = element_text(color = "black", face = "plain", size = 24),
        legend.title.position = "top",
        legend.title = element_text(face = 'plain', size = 30, hjust = 0.5),
        #legend.title = element_blank(),
        legend.text = element_text(face = 'plain', size = 24),
        legend.position = "top",
        legend.direction = "horizontal"
  ) +  labs(fill = "Strain") + xlab("Cell line") 
boxcelins

ggsave("SAMPE_Nuc_ins.svg", plot = boxcelins, width = 14, height = 8)
