ls()
rm(list=ls())
ls()

###
# Analisis de recuento manual

setwd("/home/joaquin/Desktop/Poster 2025/20250602 - Tul y Dm wt 43a1")

yo <- read.csv("sc_manual.csv")

str(yo)

summary(yo)
#Overdispersion
mean(yo$Amastigotes)
var(yo$Amastigotes)

library(tidyverse)

#For Binomial modelling
yo$Infectado <- 0 
yo$Infectado[which(yo$Amastigotes > 0)] <-1

###Add Auto results
auto <- read.csv("sc_results.csv")
auto$Infectado <- 0 
auto$Infectado[which(auto$Amastigotes > 0)] <-1


auto$Method <- "Auto"
yo$Method <- "Manual"
df <- rbind(auto, yo)



theme_custom <- theme(
  text = element_text(size = 12),
  panel.grid = element_blank(),
  panel.grid.major.y = element_line(colour = "#e3e1e1", linetype = 2),
  plot.title.position = 'plot',
  legend.position = 'top',
  legend.title = element_blank()
)

#theme_set(theme_minimal() + theme_custom)


violinx <- ggplot(df, aes(x= Cepa, y=Amastigotes, fill = Method)) + 
  geom_violin(alpha = 1, colour = "black", linewidth = 0.2) + 
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.9), size = 0.1) +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", position = position_dodge(width = 0.9), linewidth = 0.5) +
  facet_grid(.~Linea) + 
  scale_fill_manual(values = c("#ee7733","#0077bb")) + 
  theme_minimal() +
  theme(panel.spacing = unit(.2, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        text = element_text(family = "sans",color = "black", face = "plain", size = 40),
        axis.title = element_text(face = 'plain', size = 30, hjust = 0.5, margin = margin(t=0, r=0, b=10, l=0)),
        axis.text = element_text(color = "black", face = "plain", size = 24),
        legend.title.position = "top",
        legend.title = element_text(face = 'plain', size = 30, hjust = 0.5),
        #legend.title = element_blank(),
        legend.text = element_text(face = 'plain', size = 24),
        legend.position = "none",
        legend.direction = "vertical"
        ) +
  ylab("Amastigotes per cell") + xlab("Strain") 
violinx
ggsave("violin.svg", plot = violinx, width = 14, height = 8)

#%celulas infectadas
col_prob <- ggplot(df, aes(x = Cepa, y = Infectado, fill = Method)) +
  #geom_jitter(width = 0.1, alpha = 0.3) +
  stat_summary(fun = "mean", geom = "col", aes(fill = Method), position = position_dodge2(), linewidth = 0.2, colour = 'black') +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge2(), linewidth = 0.2) +
  facet_grid(~ Linea) + 
  theme_minimal() + ylab("Proportion of infected cells") + xlab("Strain")+
  coord_cartesian(ylim = c(0, 1.2)) +
  scale_fill_manual(values = c('#EE7733', '#0077BB')) +
  theme(panel.spacing = unit(.2, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        text = element_text(family = "sans",color = "black", face = "plain", size = 40),
        axis.title = element_text(face = 'plain', size = 30, hjust = 0.5, margin = margin(t=0, r=0, b=10, l=0)),
        axis.text = element_text(color = "black", face = "plain", size = 24),
        legend.title.position = "top",
        legend.title = element_text(face = 'plain', size = 30, hjust = 0.5),
        #legend.title = element_blank(),
        legend.text = element_text(face = 'plain', size = 24),
        legend.position = "none",
        legend.direction = "vertical"
  ) +
  coord_cartesian(ylim = c(0,1))
col_prob
ggsave("Percent.svg", plot = col_prob, width = 14, height = 8)

histograma <- ggplot(df1, aes(x = Amastigotes, fill = Method)) +
  geom_histogram(binwidth = 1, colour = 'black', position = 'identity', alpha = 0.9, linewidth = 0.2) +
  scale_fill_manual(values = c('#EE7733', '#0077BB')) + theme_minimal() + 
  theme(panel.spacing = unit(.2, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        text = element_text(family = "sans",color = "black", face = "plain", size = 40),
        axis.title = element_text(face = 'plain', size = 30, hjust = 0.5, margin = margin(t=0, r=0, b=10, l=0)),
        axis.text = element_text(color = "black", face = "plain", size = 24),
        legend.title.position = "top",
        legend.title = element_text(face = 'plain', size = 30, hjust = 0.5),
        #legend.title = element_blank(),
        legend.text = element_text(face = 'plain', size = 24),
        legend.position = "none",
        legend.direction = "vertical"
  ) +
  ylab("Frequency")
histograma
ggsave("histograma.svg", plot = histograma, width = 16, height = 8)

#Modelling
library(glmmTMB)

#Poisson model with zero inflation
mxp <- glmmTMB(Amastigotes ~ Method + Cepa + Linea + (1|Campo) , family = 'poisson', data = df)
summary(mxp)

#Poisson model with zero inflation
mxp1 <- glmmTMB(Amastigotes ~ Method + Cepa + Linea + (1|Campo) , family = 'poisson', data = df, ziformula = ~1)
summary(mxp1)

mx1 <- glmmTMB(Amastigotes ~Method, family = 'nbinom2', data = df)
summary(mx1)

mx2 <- glmmTMB(Amastigotes ~Method + Cepa , family = 'nbinom2', data = df)
summary(mx2)

mx3 <- glmmTMB(Amastigotes ~ Method + Cepa + Linea, family = 'nbinom2', data = df)
summary(mx3)

mx4 <- glmmTMB(Amastigotes ~ Method + Cepa + Linea + (1|Campo), family = 'nbinom2', data = df)
summary(mx4)

mx5 <- glmmTMB(Amastigotes ~ Method + Cepa + Linea + (1|Campo) , family = 'nbinom2', data = df, ziformula = ~1)
summary(mx5)

mx6 <- glmmTMB(Amastigotes ~ Method * Cepa * Linea + (1|Campo) , family = 'nbinom2', data = df, ziformula = ~1)
summary(mx6, ddf = 'Kenward-Roger') 

#Bayesian
library(brms)
bayes_model <- brm(
  formula = Amastigotes ~ Method * Cepa * Linea + (1 | Campo),
  data = df,
  family = zero_inflated_negbinomial(),
  chains = 3,           # Number of MCMC chains
  iter = 2000,          # Total iterations per chain (you can lower it for testing)
  warmup = 1000,        # Burn-in samples
  cores = 2,            # Use multiple cores if available
  seed = 1234,          # For reproducibility
  control = list(adapt_delta = 0.95)  # Helps avoid divergent transitions
)

summary(bayes_model)
posterior_summary(bayes_model)
pp_check(bayes_model)
conditional_effects(bayes_model)
conditional_effects(bayes_model, effects = "Method:Cepa:Linea")
emmbayes <- emmeans(bayes_model, ~ Method | Cepa * Linea)
plot(emmbayes)

#Choose best model
AIC(mx1, mx2, mx3, mx4, mx5, mxp, mx6)

#Check residuals\
library(DHARMa)
res <- simulateResiduals(mx6)
plot(res) #m1 and m6 have bad qq plot and heterogeneous variance

#Comparisons. Simple effects for method
library(emmeans)
emm <- emmeans(mx6, ~ Method | Cepa * Linea) #Marginal means
emm

plot(pairs(emm, adjust = "fdr"), contrasts = TRUE) #comparisons
contrast_df <- as.data.frame(pairs(emm, adjust = "fdr"))

#95% Confidence intervals (CL = estimate +- tcritical * SE)
alpha <- 0.05
contrast_df <- contrast_df %>%
  mutate(
    t_crit = qt(1 - alpha/2, df),
    lower.CL = estimate - t_crit * SE,
    upper.CL = estimate + t_crit * SE
  )

library(ggthemes)
contrast_plot <- ggplot(contrast_df, aes(x = interaction(Cepa, Linea), y = estimate, color = interaction(Cepa, Linea))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),width = 0.3, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Difference in estimated means (Automatic - Manual)",
       x = "Strain and Cell Line",
       title = "Method diffrences between infections") +
  theme_minimal() + coord_cartesian(ylim = c(-2, 2))+ scale_color_pander()
contrast_plot 
ggsave("amaspercel_contrasts.svg", plot = contrast_plot, width = 10, height = 6)

#Column plots
avg <- aggregate(df$Amastigotes, list(df$Cepa, df$Linea, df$Method), FUN = mean)
column_amas <- ggplot(avg, aes(x = Group.1, y = x, fill = Group.3)) + geom_col(position = "dodge") + 
  facet_grid(~Group.2)
column_amas + scale_fill_pander()

# Plotting the predictions from the model
#emm_interactions <- emmeans(mx6, ~ Method * Cepa * Linea, type = "response")
emm_df1 <- as.data.frame(emm)

#Model means with comparison CI
model_amaspercel <- ggplot(emm_df1, aes(x = Cepa, y = emmean, fill = Method)) +
  geom_col(position = position_dodge2(width = 0.9), colour = 'black') +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width = 0.9)) +
  facet_grid( ~ Linea) +
  labs(y = "Mean amastigotes per cell", x = "Parasite Strain",title = "Model-estimated means with 95% CI")+
  theme_minimal() + scale_fill_stata()
model_amaspercel
ggsave("amas_por_cel_model.svg", plot = model_amaspercel, width = 10, height = 6)

#Raw means with CI
raw_amaspercel <- ggplot(df, aes(x = Cepa, y = Amastigotes, fill = Method)) +
  #geom_jitter(width = 0.1, alpha = 0.3) +
  stat_summary(fun = "mean", geom = "col", aes(fill = Method), position = position_dodge2(), colour = 'black') +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge2()) +
  facet_grid(~ Linea) + theme_minimal() + coord_cartesian(ylim = c(0, 20)) + ylab("Amastigotes per cell") + xlab("Strain")
raw_amaspercel
ggsave("amas_por_cel_raw.svg", plot = raw_amaspercel, width = 10, height = 6)


####
#Just infected cells
df1 <- subset(df, df$Amastigotes > 0)

median(df1$Amastigotes)#6.570147
var(df1$Amastigotes)#40.15875
min(df1$Amastigotes)#1
max(df1$Amastigotes)#47

violinfected <- ggplot(df1, aes(x = Cepa, y = Amastigotes, fill = Method)) + geom_violin() + facet_grid(~Linea)
violinfected

violinfected <- ggplot(df1, aes(x= Cepa, y=Amastigotes, fill = Method)) + 
  geom_violin(alpha = 1, colour = "black", linewidth = 0.2) + 
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.9), size = 0.1) +
  stat_summary(fun.data = mean_cl_normal, geom = "pointrange", position = position_dodge(width = 0.9), linewidth = 0.5) +
  facet_grid(.~Linea) + 
  scale_fill_manual(values = c("#ee7733","#0077bb")) + 
  theme_minimal() +
  theme(panel.spacing = unit(.2, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        text = element_text(family = "sans",color = "black", face = "plain", size = 40),
        axis.title = element_text(face = 'plain', size = 30, hjust = 0.5, margin = margin(t=0, r=0, b=10, l=0)),
        axis.text = element_text(color = "black", face = "plain", size = 24),
        legend.title.position = "top",
        legend.title = element_text(face = 'plain', size = 30, hjust = 0.5),
        #legend.title = element_blank(),
        legend.text = element_text(face = 'plain', size = 24),
        legend.position = "none",
        legend.direction = "vertical"
  ) +
  ylab("Amastigotes per infected cell") + xlab("Strain") 
violinfected
ggsave("violinfected.svg", plot = violinfected, width = 14, height = 8)

raw_amasperinfcel <- ggplot(df1, aes(x = Cepa, y = Amastigotes, fill = Method)) +
  #geom_jitter(width = 0.1, alpha = 0.3) +
  stat_summary(fun = "mean", geom = 'col', position = position_dodge2(width = 0.9), colour = 'black') +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge2(width = 0.9)) +
  facet_grid(~ Linea) + 
  theme_minimal() + 
  coord_cartesian(ylim = c(0, 20)) + ylab("Amastigotes per infected cell") + xlab("Strain")
raw_amasperinfcel
ggsave("amas_por_celinf_raw.svg", plot = raw_amasperinfcel, width = 10, height = 6)


annotations <- data.frame(
  x = c(round(min(df1$Amastigotes), 2), round(median(df1$Amastigotes), 2), round(max(df1$Amastigotes), 2)),
  x_pos = c(round(min(df1$Amastigotes), 2), 6, round(max(df1$Amastigotes), 2)),
  y = c(135, 80, 5),
  label = c("Min:", "Median:", "Max:"),
  Method = c("Both", "Both","Both")
) 
#df1$Method <- factor(df1$Method, levels = c("Manual", "Auto"))
histinf <- ggplot(df1, aes(x = Amastigotes, fill = Method)) +
  geom_histogram(binwidth = 1, colour = 'black', position = 'identity', alpha = 0.8) + #facet_grid(rows  = vars(Method))+
  xlab("Number of intracellular amastigotes") + ylab("Number of cells") +
  #geom_text(data = annotations, aes(x = x_pos, y = y, label = paste(label, x)), size = 5, fontface = "bold") +
  theme_minimal() +
  theme(element_text(face = 'bold'))
histinf
ggsave("hist_inf_overlap.svg", plot = histinf, width = 14, height = 6)

manual <- subset(df1, df1$Method == "Manual")
auto <- subset(df1, df1$Method == "Auto")

densinf_top <- ggplot(auto, aes(x = Amastigotes, color = Method, fill = Method)) +
  geom_density(colour = "#F8766D", fill = "#F8766D", alpha = 0.4) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
densinf_top

densinf_bottom <- ggplot(manual, aes(x = Amastigotes, color = Method, fill = Method)) +
  geom_density(colour = "#00BFC4", fill = "#00BFC4",alpha = 0.4) +
  theme_minimal() +
  scale_y_reverse()
densinf_bottom

(densinf_top + densinf_bottom) + plot_layout(ncol = 1, heights = c(1, 1))

(densinf_top /plot_spacer()/ densinf_bottom) + plot_layout(heights = c(1, 0, 1))

denshistinf <- ggplot(df1, aes(x = Amastigotes, color = Method)) +
  geom_histogram(aes(y = ..density..), binwidth = 1, fill = "gray", color = "black") +
  geom_density(color = "red", size = 1.2) +
  labs(title = "Histogram + Density Curve") 
denshistinf #+  facet_grid(Cepa~Linea)

modes <- density(df1$Amastigotes)
plot(modes)

#Models
#Poisson model
mip <- glmmTMB(Amastigotes ~ Method + Cepa + Linea + (1|Campo) , family = 'poisson', data = df1)
summary(mip)

#Poisson model with zero inflation
mip1 <- glmmTMB(Amastigotes ~ Method + Cepa + Linea + (1|Campo) , family = 'poisson', data = df1, ziformula = ~1)
summary(mip1)

mip2 <- glmmTMB(Amastigotes ~ Method + Cepa + Linea + (1|Campo) + (1|Campo:Célula) , family = truncated_poisson(), data = df1, ziformula = ~1)
summary(mip2)

mi1 <- glmmTMB(Amastigotes ~ Method, family = 'nbinom2', data = df1)
summary(mi1)

mi2 <- glmmTMB(Amastigotes ~Method + Cepa , family = 'nbinom2', data = df1)
summary(mi2)

mi3 <- glmmTMB(Amastigotes ~ Method + Cepa + Linea, family = 'nbinom2', data = df1)
summary(mi3)

mi4 <- glmmTMB(Amastigotes ~ Method + Cepa + Linea + (1|Campo), family = 'nbinom2', data = df1)
summary(mi4)

mi5 <- glmmTMB(Amastigotes ~ Method * Cepa * Linea + (1|Campo) , family = 'nbinom2', data = df1)
summary(mi5)

mi6 <- glmmTMB(Amastigotes ~ Method + Cepa + Linea + (1|Campo) , family = 'nbinom2', data = df1, ziformula = ~1)
summary(mi6)

mi7 <- glmmTMB(Amastigotes ~ Method * Cepa * Linea + (1|Campo) , family = 'nbinom2', data = df1, ziformula = ~1)
summary(mi7)

#micom <- glmmTMB(Amastigotes ~ Method * Cepa * Linea + (1|Campo) , family = 'compois', data = df1)
#summary(micom)

miq <- glm(Amastigotes ~ Method * Cepa * Linea + (1|Campo) , family = 'quasipoisson', data = df1)
summary(miq)

mi8 <- glmmTMB(Amastigotes ~ Method * Cepa * Linea + (1|Campo), family = truncated_nbinom2(), data = df1)
summary(mi8)

mi9 <- glmmTMB(Amastigotes ~ Method * Cepa * Linea + (1|Campo) + (1|Célula), family = truncated_nbinom2(), data = df1)
summary(mi9)

mi10 <- glmmTMB(Amastigotes ~ Method * Cepa * Linea + (1|Campo) + (1|Campo:Célula), family = truncated_nbinom2(), data = df1)
summary(mi10)

mi11 <-  glmmTMB(Amastigotes ~ Method * Cepa * Linea + (1|Campo) + (1|Campo:Célula), family = nbinom2(), data = df1)
summary(mi11)

#micom1 <- glmmTMB(Amastigotes ~ Method * Cepa * Linea + (1|Campo), family = truncated_compois(), data = df1)
#summary(micom1)

library(quantreg)

model_qr <- rq(Amastigotes ~ Method * Cepa * Linea, tau = 0.5, data = df1)
summary(model_qr)


AIC(mip, mip1, mip2, mi1, mi2, mi3, mi4, mi5, mi6, mi7,mi8, mi9, mi10, miq, model_qr, model_gam) #mi5

#Check residuals
res <- simulateResiduals(mi9)
plotQQunif(res) #m1 and m6 have bad qq plot and heterogeneous variance
plotResiduals(res)

plot(
  model_qr$fitted.values,
  residuals(model_qr, type = "pearson"),
  xlab = "Fitted Values",
  ylab = "Pearson Residuals",
  main = "Residuals vs. Fitted Values (Quasi-Poisson)"
)
abline(h = 0, col = "red")

#Check observed vs predicted values

emminf <- emmeans(mi10, ~ Method * Cepa * Linea, type = "response")
fitted_df1<- as.data.frame(emminf)
obs_means <- df1 %>%
  group_by(Method, Cepa, Linea) %>%
  summarise(observed_mean = mean(Amastigotes, na.rm = TRUE),
            .groups = "drop")
plot_df1 <- left_join(fitted_df1, obs_means, by = c("Method", "Cepa", "Linea"))
ggplot(plot_df1, aes(x = response, y = observed_mean)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = asymp.LCL, xmax = asymp.UCL), height = 0.2, color = "gray40") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Fitted mean (from model)",
       y = "Observed mean (raw data)",
       title = "Fitted vs Observed Means per Group") +
  theme_minimal() + coord_cartesian(ylim = c(0,20), xlim = c(0,20))


#No consigo mejor que mi10

#Modelos no parametricos
library(mgcv)

model_gam <- gam(Amastigotes ~ s(Method) + s(Cepa) + s(Linea) + s(Campo, bs = "re"),
                 family = nb(),
                 data = df1,
                 method = "REML")
summary(model_gam)


#Efectos simples por metodo
emm_inf <- emmeans(mi10, ~ Method | Cepa * Linea)
emm_inf

plot(pairs(emm_inf, adjust = "fdr"), contrasts = TRUE) #comparisons
pairs(emm_inf, adjutst = 'fdr')
contrast_df <- as.data.frame(pairs(emm_inf, adjust = "fdr"))

#95% Confidence intervals (CL = estimate +- tcritical * SE)
alpha <- 0.05
contrast_df <- contrast_df %>%
  mutate(
    t_crit = qt(1 - alpha/2, df),
    lower.CL = estimate - t_crit * SE,
    upper.CL = estimate + t_crit * SE
  )


contrast_plot <- ggplot(contrast_df, aes(x = interaction(Cepa, Linea), y = estimate, color = interaction(Cepa, Linea))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),width = 0.3, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Difference in estimated means (Automatic - Manual)",
       x = "Strain and Cell Line",
       title = "Method diffrences between infections") +
  theme_minimal() + coord_cartesian(ylim = c(-2, 2)) + scale_color_pander()
contrast_plot 
ggsave("amaspercelinf_contrasts.svg", plot = contrast_plot, width = 10, height = 6)


#Model means with simple effects CI
emm_inf_df <- as.data.frame(emm_inf)
model_amaspercelinf <- ggplot(emm_inf_df, aes(x = Cepa, y = emmean, fill = Method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),width = 0.2, position = position_dodge(width = 0.9)) +
  facet_grid( ~ Linea) +
  labs(y = "Mean amastigotes per cell", x = "Parasite Strain",title = "Model-estimated means with 95% CI")+
  theme_minimal()
model_amaspercelinf
ggsave("amas_por_celinf_model.svg", plot = model_amaspercelinf, width = 10, height = 6)


### Modelos binomiales (%cel inf)

mb1 <- glmmTMB(Infectado ~ Method, family = 'binomial', data = df)
summary(mb1)

mb2 <- glmmTMB(Infectado ~ Method + Cepa + Linea, family = 'binomial', data = df)
summary(mb2)

mb3 <- glmmTMB(Infectado ~ Method * Cepa * Linea, family = 'binomial', data = df)
summary(mb3)

mb4 <- glmmTMB(Infectado ~ Method + Cepa + Linea + (1|Campo), family = 'binomial', data = df)
summary(mb4)

mb5 <- glmmTMB(Infectado ~ Method * Cepa * Linea + (1|Campo), family = 'binomial', data = df)
summary(mb5)


#Best model mb4
AIC(mb1, mb2, mb3, mb4, mb5)

#Check fit
res <- simulateResiduals(mb4)
plot(res)

#Contrasts
emm_prob_mean <- emmeans(mb5, ~ Method, type = "response")
pairs(emm_prob_mean)
emm_prob <-emmeans(mb5, ~ Method | Cepa * Linea, type = 'response') #Percent infected
pairs(emm_prob, adjust = 'fdr')

porbs_df <- predict(mb5, df, type = 'response' )

contrast_prob_df <- as.data.frame(pairs(emm_prob, adjust = "fdr"))

#95% Confidence intervals (CL = estimate +- tcritical * SE)
alpha <- 0.05
contrast_prob_df <- contrast_prob_df %>%
  mutate(
    t_crit = qt(1 - alpha/2, df),
    lower.CL = odds.ratio - t_crit * SE,
    upper.CL = odds.ratio + t_crit * SE
  )


contrast_plot <- ggplot(contrast_prob_df, aes(x = interaction(Cepa, Linea), y = odds.ratio, color = interaction(Cepa, Linea))) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),width = 0.3, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(y = "Odds ratio (Automatic / Manual)",
       x = "Strain and Cell Line",
       title = "Method diffrences between infections") +
  scale_color_pander()+ 
  theme_minimal()# + coord_cartesian(ylim = c(-2, 2)) 
contrast_plot 
ggsave("percent_contrasts.svg", plot = contrast_plot, width = 10, height = 6)

#Raw values
col_prob <- ggplot(df, aes(x = Cepa, y = Infectado, fill = Method)) +
  #geom_jitter(width = 0.1, alpha = 0.3) +
  stat_summary(fun = "mean", geom = "col", aes(fill = Method), position = position_dodge2(), , colour = 'black') +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge2()) +
  facet_grid(~ Linea) + 
  theme_minimal() + ylab("Percent infected cells") + coord_cartesian(ylim = c(0, 1.2)) +
  scale_fill_manual(values = c('#EE7733', '#0077BB')) 
col_prob
ggsave("Percent.svg", plot = col_prob, width = 10, height = 6)

scatter <- ggplot(df, aes(y = Infectado, x = Method, fill = Cepa)) + 
  geom_jitter(height = 0.05, alpha = 0.5, size = 3, shape = 21) + 
  facet_grid(~Linea)
scatter

###Resumen por campo

#Amas por campo
fov <- aggregate(df$Amastigotes, list(df$Cepa, df$Linea, df$Campo, df$Method), FUN = sum)

#Renombrar columnas
colnames(fov) <- c("Cepa","Linea", "Campo", "Método", "Amastigotes")

#Numero de Celulas
cel <- aggregate(df$Amastigotes, list(df$Cepa, df$Linea, df$Campo, df$Method), FUN = length)

fov$Células <- cel$x #Agrego numero de celulas a base de datos

#Amas por celula
amas <- aggregate(df$Amastigotes, list(df$Cepa, df$Linea, df$Campo, df$Method), FUN = mean) 

#Da lo mismo dividir amas totales por celulas totales que hacer el promedio de amas
fov$Amas_Célula <- amas$x #Agrego numero de celulas a base de datos

#Base de datos de solo celulas infectadas
inf <- subset(df, df$Amastigotes > 0)

#Numero de celulas infectadas
celinf <- aggregate(inf$Amastigotes, list(inf$Cepa, inf$Linea, inf$Campo, inf$Method), FUN = length)

#Union de resumenes de numero de celulas totales (x.x) e infectadas(x.y)
cels <- merge(cel, celinf, by= c("Group.1", "Group.2", "Group.3", "Group.4"), all.x = TRUE)
cels$x.y <- replace_na(cels$x.y, 0) #Campos con ninguna celula infectada

#Ordenar para que los datos queden en las filas correctas
cels <- cels[order(cels$Group.1, cels$Group.2, cels$Group.3, cels$Group.4),]
fov <- fov[order(fov$Cepa, fov$Linea, fov$Campo, fov$'Método'),]

#Copiar datos
fov$'Células infectadas' <- cels$x.y

#Amas por celula infectada
fov$'Amastigotes por célula infectada' <- fov$Amastigotes / fov$`Células infectadas`
fov$`Amastigotes por célula infectada` <- replace_na(fov$`Amastigotes por célula infectada`, 0)

#Porcentaje de celulas infectadas
fov$'Porcentaje de células infectadas' <- (fov$`Células infectadas` / fov$Células) *100

#Indice de infeccion
fov$index <- fov$`Porcentaje de células infectadas` * fov$`Amastigotes por célula infectada`



fov_wide <- fov %>%
  pivot_wider(
    id_cols = c(Cepa, Linea, Campo),
    names_from = Método,     # Creates new columns from the 'Método' values
    values_from = Amastigotes  # Fills those new columns with 'Amastigotes' values
  )

#Amas esperados por campo vs amas obtenidos
scatter <- ggplot(fov_wide, aes(x = Manual, y = Auto, colour = Cepa)) +
           geom_point(size = 3) + geom_abline(slope = 1)
scatter

indice <- ggplot(fov, aes(x = Cepa, y = `Índice de infección`, fill = Método)) + 
  stat_summary(fun = "mean", geom = "col", position = position_dodge2(), colour = 'black', linewidth = 0.2) + 
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge2(), linewidth = 0.2) + 
  facet_grid(~Linea) + coord_cartesian(ylim = c(0, NA)) + 
  scale_fill_manual(values = c("#ee7733","#0077bb")) + 
  theme_minimal() +
  theme(panel.spacing = unit(.2, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        text = element_text(family = "sans",color = "black", face = "plain", size = 40),
        axis.title = element_text(face = 'plain', size = 30, hjust = 0.5, margin = margin(t=0, r=0, b=10, l=0)),
        axis.text = element_text(color = "black", face = "plain", size = 24),
        legend.title.position = "top",
        legend.title = element_text(face = 'plain', size = 30, hjust = 0.5),
        #legend.title = element_blank(),
        legend.text = element_text(face = 'plain', size = 24),
        legend.position = "none",
        legend.direction = "vertical"
  ) +
  ylab("Infection index") + xlab("Strain") 
indice
ggsave("Index.svg", plot = indice, width = 14, height = 8)

#Modelo indice
str(fov)
fov$Cepa <- as.factor(fov$Cepa)
fov$Linea <- as.factor(fov$Linea)
fov$Método <- as.factor(fov$Método)

mix1 <- lm(`Índice de infección` ~ Método, data = fov)
summary(mix1)
mix2 <- lm(`Índice de infección` ~ Método + Cepa + Linea, data = fov)
summary(mix2)
mix3 <- lm(`Índice de infección` ~ Método * Cepa * Linea, data = fov)
summary(mix3)
library(nlme)
mix4 <- gls(index~ Método * Cepa * Linea, weights = varIdent(form~1|Linea), data = fov)
summary(mix4)
mix5 <- gls(index~ Método * Cepa * Linea, weights = varIdent(form~1|Método * Cepa * Linea), data = fov)
summary(mix5)
AIC(mix1, mix2, mix3, mix4, mix5)

library(np)
bw <- npregbw(formula = `Índice de infección` ~ Método + Cepa + Linea, regtype = 'll', bmmethod="cv.aic", data = fov)
mix6 <- mpreg <- (bws = bw)
residuals <- resid(mix6)
fitted <- fitted(mix6)
# Residual plot
plot(fitted, residuals, pch = 16, col = "darkgray")
abline(h = 0, col = "red")
hist(residuals, main = "Residuals", xlab = "Residual", breaks = 20)

plot(mix6, plot.errors.method = "bootstrap", plot.errors.boot.num = 25)

res <- simulateResiduals(mix6)
plot(res)

shapiro.test(res)

res <- resid(mix6)
plot_data <- data.frame(residuals = res, category = fov$Linea)
boxplot(residuals ~ category, data = plot_data)
leveneTest(index ~ Linea, data = fov)

#Supuestos
e<-resid(mix6) # residuos de pearson
pre<-predict(mix6) #predichos
par(mfrow = c(1, 2))
plot(pre, e, xlab="Predichos", ylab="Residuos de pearson",main="Gr?fico de dispersi?n de RE vs PRED",cex.main=.8 )
abline(0,0)
qqnorm(e, cex.main=.8)
qqline(e)
par(mfrow = c(1, 1))
plot(mix4) #RP vs PRED
shapiro.test(e)
plot_model(mix4, type="diag")[[2]]

##Graficar todas las medidas juntas
library(hrbrthemes)
library(GGally)
library(viridis)


# Plot
coord_plot <- ggparcoord(fov,
           columns = c(7,9,10), groupColumn = 4, order = "skewness",
           showPoints = TRUE, 
           scale = "uniminmax",
           title = "Parallel Coordinate Plot for Infection Parameters",
           alphaLines = 0.5) + 
  scale_color_manual(values = c('#EE7733', '#0077BB')) +
  #theme_ipsum()+
  facet_grid(Linea ~ Cepa)+
  theme(plot.title = element_text(size=10)) 
coord_plot


###Incluyendo metodo morfologico
ins <- read.csv("sc_INsPECT.csv")
ins$Infectado <- 0 
ins$Infectado[which(ins$Amastigotes > 0)] <-1

ins$Method <- "INsPECT"
df <- rbind(df, ins)

df$Method <- factor(df$Method, levels = c("Auto", "Manual", "INsPECT"))

violins <- ggplot(df, aes(x = Cepa, y = Amastigotes, fill = Method)) + geom_violin(alpha = 0.5)+ facet_grid(.~Linea)
violins

mb1 <- glmmTMB(Amastigotes ~ Method, family = 'nbinom2', data = df)
summary(mb1)

mb2 <- glmmTMB(Amastigotes ~ Cepa + Method + (1|Campo), family = "nbinom2", data = df)
summary(mb2)

mb3 <- glmmTMB(Amastigotes ~ Cepa * Method * Linea + (1|Campo), family = "nbinom2", data = df)
summary(mb3)

comps <- emmeans(mb3, pairwise ~ Cepa + Method + Linea, type = "response")
comps
comps2 <-emmeans(mb3, ~ Method, type = "response")
comps2
estad2<-as.data.frame(comps$emmeans)
# Las barras de error representan IC
ggplot(estad2, aes(x=Method, y=response, color = Cepa)) +
  geom_errorbar(aes(ymin=asymp.LCL, ymax=asymp.UCL), width=.1) +# facet_grid(.~Method) + 
  geom_point(size = 3)+ ylab("log Amastigotes") +   theme_grey(base_size = 16) #+  annotate("text", x = c(1,2,3), y = c(1, 5.3, 6.4), label = c("A", "B", "B"))


emm <- emmeans(mb3, ~ Method | Cepa * Linea) #Marginal means
emm

plot(pairs(emm, adjust = "fdr"), contrasts = TRUE) #comparisons
pairs(emm, adjust = "fdr")
contrast_df <- as.data.frame(pairs(emm, adjust = "fdr"))

#95% Confidence intervals (CL = estimate +- tcritical * SE)
alpha <- 0.05
contrast_df <- contrast_df %>%
  mutate(
    t_crit = qt(1 - alpha/2, df),
    lower.CL = estimate - t_crit * SE,
    upper.CL = estimate + t_crit * SE
  )

contrast_plot <- ggplot(contrast_df, aes(x = interaction(Cepa, Linea), y = estimate, color = contrast)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),width = 0.3, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Difference in estimated means (Automatic - Manual)",
       x = "Strain and Cell Line",
       title = "Method diffrences between infections") +
  theme_minimal() + coord_cartesian(ylim = c(-2, 2))+ scale_color_pander()
contrast_plot 
ggsave("amaspercel_contrasts.svg", plot = contrast_plot, width = 10, height = 6)

#Column plots
avg <- aggregate(df$Amastigotes, list(df$Cepa, df$Linea, df$Method), FUN = mean)
column_amas <- ggplot(avg, aes(x = Group.1, y = x, fill = Group.3)) + geom_col(position = "dodge") + 
  facet_grid(~Group.2)
column_amas + scale_fill_pander()

# Plotting the predictions from the model
#emm_interactions <- emmeans(mx6, ~ Method * Cepa * Linea, type = "response")
emm_df1 <- as.data.frame(emm)

#Model means with comparison CI
model_amaspercel <- ggplot(emm_df1, aes(x = Cepa, y = emmean, fill = Method)) +
  geom_col(position = position_dodge2(width = 0.9), colour = 'black') +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), position = position_dodge2(width = 0.9)) +
  facet_grid( ~ Linea) +
  labs(y = "Mean amastigotes per cell", x = "Parasite Strain",title = "Model-estimated means with 95% CI")+
  theme_minimal() + scale_fill_stata()
model_amaspercel
ggsave("amas_por_cel_model.svg", plot = model_amaspercel, width = 10, height = 6)

#Raw means with CI
raw_amaspercel <- ggplot(df, aes(x = Cepa, y = Amastigotes, fill = Method)) +
  #geom_jitter(width = 0.1, alpha = 0.3) +
  stat_summary(fun = "mean", geom = "col", aes(fill = Method), position = position_dodge2(), colour = 'black') +
  stat_summary(fun.data = mean_cl_normal, geom = "errorbar", position = position_dodge2()) +
  facet_grid(~ Linea) + theme_minimal() + coord_cartesian(ylim = c(0, 20)) + ylab("Amastigotes per cell") + xlab("Strain")
raw_amaspercel
ggsave("amas_por_cel_raw.svg", plot = raw_amaspercel, width = 10, height = 6)


##Celulas inefctadas
inf2 <- subset(df, df$Infectado == 1)

violins <- ggplot(inf2, aes(x = Cepa, y = Amastigotes, fill = Method)) + geom_violin(alpha = 0.5)+ facet_grid(.~Linea)
violins

mb1 <- glmmTMB(Amastigotes ~ Method, family = 'nbinom2', data = df)
summary(mb1)

mb2 <- glmmTMB(Amastigotes ~ Cepa + Method + (1|Campo), family = "nbinom2", data = df)
summary(mb2)

mbi3 <- glmmTMB(Amastigotes ~ Cepa * Method * Linea + (1|Campo), family = "nbinom2", data = inf2)
summary(mbi3)

emm <- emmeans(mbi3, ~ Method | Cepa * Linea) #Marginal means
emm

plot(pairs(emm, adjust = "fdr"), contrasts = TRUE) #comparisons
pairs(emm, adjust = "fdr")
contrast_df <- as.data.frame(pairs(emm, adjust = "fdr"))

#95% Confidence intervals (CL = estimate +- tcritical * SE)
alpha <- 0.05
contrast_df <- contrast_df %>%
  mutate(
    t_crit = qt(1 - alpha/2, df),
    lower.CL = estimate - t_crit * SE,
    upper.CL = estimate + t_crit * SE
  )

contrast_plot <- ggplot(contrast_df, aes(x = interaction(Cepa, Linea), y = estimate, color = contrast)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) + 
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),width = 0.3, position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "Difference in estimated means (Automatic - Manual)",
       x = "Strain and Cell Line",
       title = "Method diffrences between infections") +
  theme_minimal() + coord_cartesian(ylim = c(-2, 2))+ scale_color_pander()
contrast_plot 

###### Comparacion de operadores
ops <- read.csv("sc_results_all.csv")

ops$Cepa <- as.factor(ops$Cepa)
ops$Operador <- as.factor(ops$Operador)

ops$Operador <- factor(comp$Operador, levels = c("Auto", "Joaquin", "Salo", "Inspect"))

library(ggthemes)
plot <- ggplot(comp, aes(x = Cepa, y = Parasitos, fill = Operador)) + geom_violin()
plot

#Statistical models
mops <- glmmTMB(Parasitos ~ Cepa + Operador + (1|Replica) + (1|Replica:Campo), family = "nbinom2", data = ops)
summary(mops)

emm1 <- emmeans(mops, "Parasitos")

compsops <- emmeans(mops,  ~ Cepa + Operador, type = "response")
compsops
estadops <- as.data.frame(compsops$emmeans)

ggplot(estadops, aes(x=Operador, y=response, color = Cepa)) +
  geom_errorbar(aes(ymin=asymp.LCL, ymax=asymp.UCL), width=.1) +# facet_grid(.~Method) + 
  geom_point(size = 3)+ ylab("log Amastigotes") +   theme_grey(base_size = 16) #+  annotate("text", x = c(1,2,3), y = c(1, 5.3, 6.4), label = c("A", "B", "B"))
contrast(compsops, method = "pairwise")
