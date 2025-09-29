#PLANT HEIGHT DEVELOPMENT DURING SIX WEEKS
#Import Height_data excel sheet

Data_1 <-subset(Height_data, Substrate!="Only maize" & Substrate!="Maize + min." & Substrate!="Hornum + min." & ID!="89" & ID!="98")

library(dplyr)
# Rename 
Data_1 <- Data_1 %>%
  mutate(Degree = case_when(
    Degree == 1 ~ 1,
    Degree == 3 ~ 3,
    Degree == 5 ~ 5,
    Degree == 7 ~ 7,
    Degree == 10 ~ 9,
    TRUE ~ Degree  # If none of the conditions match, keep the original value
  ))

Data_2 <-subset(Data_1, Day=="39")
Data_3 <-subset(Data_1, Day!="6")


# Rename 
Data_1 <- Data_1 %>%
  mutate(Degree = case_when(
    Degree == 1 ~ 1,
    Degree == 3 ~ 3,
    Degree == 5 ~ 5,
    Degree == 7 ~ 7,
    Degree == 10 ~ 9,
    TRUE ~ Degree  # If none of the conditions match, keep the original value
  ))


#mixed linear model testing the effects on plant height (cm)

library(lme4)

Height_model <- lmer(Height ~ Substrate * Cultivar * Degree * Day + (1 | ID), data = Data_1)

AIC(Height_model)
BIC(Height_model)

residuals <- resid(Height_model)

#produce residual vs. fitted plot
plot(fitted(Height_model), residuals)
#add a horizontal line at 0 
abline(0,0)

qqnorm(residuals)
abline(0,sd(residuals))

# Assuming 'Height_model' is your linear mixed-effects model object
summary_output <- summary(Height_model)

library(lmerTest)

summary(Height_model)

(aov <- anova(Height_model))


library(emmeans)

joint_tests(Height_model)

#Tukey test for each cultivar at each fertilizer at day 39

Height_model_2 <- lmer(Height ~ Substrate * Cultivar * as.factor(Degree) * as.factor(Day) + (1 | ID), data = Data_1)

day_emm <- emmeans(Height_model_2, ~ Degree | Day*Substrate*Cultivar, type = "response")

pwpm(day_emm)

library(multcomp)
# Compute the CLDs using the cld() function
cld_results <- cld(day_emm, Letters = letters, reversed = TRUE)
print(cld_results)






#Height for each cultivar at day 39

Height_model_day39 <- lm(Height ~ Substrate * Cultivar * Degree, data=Data_2)

summary(Height_model_day39)

residuals <- resid(Height_model_day39)

#produce residual vs. fitted plot
plot(fitted(Height_model_day39), residuals)
#add a horizontal line at 0 
abline(0,0)

qqnorm(residuals)
abline(0,sd(residuals))

joint_tests(Height_model_day39)

#Tukey test for each cultivar at each fertilizer at day 39

day39_emm <- emmeans(Height_model_day39, ~ Dilution | Cultivar*Substrate, type = "response")
pwpm(day39_emm)

library(multcomp)
# Compute the CLDs using the cld() function
cld_results <- cld(day39_emm, Letters = letters, reversed = TRUE)
print(cld_results)



#Combined model testing the effects of Cultivar*Fertilizer*Microbiome on shoot_dry_weight

Data_1 <-subset(Dw_data, Sample_ID!="89" & Sample_ID!="98")

LM_1000 <- lm(log(Shoot_corr_mg) ~ Cultivar * Fertilizer * new_microbiome, data = Data_1)

residuals <- resid(LM_1000)

#produce residual vs. fitted plot
plot(fitted(LM_1000), residuals)
#add a horizontal line at 0 
abline(0,0)

qqnorm(residuals)
abline(0,sd(residuals))

summary(LM_1000)

joint_tests(LM_1000)

#seperate linear regressions for each cultivar at each fertilizer treatment

# For Babushka and Chitin
LM_bab_chi <- lm(Shoot_corr_mg ~ new_microbiome, data = Data_1, subset = Cultivar == "Babushka" & Fertilizer == "Chitin")

summary(LM_bab_chi)

# For Babushka and mineral
LM_bab_min <- lm(Shoot_corr_mg ~ new_microbiome, data = Data_1, subset = Cultivar == "Babushka" & Fertilizer == "Mineral")

summary(LM_bab_min)


# For RGT and chitin
LM_RGT_chi <- lm(Shoot_corr_mg ~ new_microbiome, data = Data_1, subset = Cultivar == "RGT" & Fertilizer == "Chitin")

summary(LM_RGT_chi)

# For RGT and mineral
LM_RGT_min <- lm(Shoot_corr_mg ~ new_microbiome, data = Data_1, subset = Cultivar == "RGT" & Fertilizer == "Mineral")

summary(LM_RGT_min)


#PLOT WITH GEOM SMOOTH (LM)
library(ggplot2)

# Define axis labels
x.axis.label <- "Microbiome dilution"
y.axis.label <- "Shoot dry weight (mg)"

# Create the plot using facet grid
plot <- ggplot(Data_1, aes(new_microbiome, Shoot_corr_mg, shape = Cultivar)) +
  geom_smooth(data = subset(Data_1, Cultivar == "Babushka"),color = "#ff791f", method = "lm", se = TRUE, size = 2) +
  geom_smooth(data = subset(Data_1, Cultivar == "RGT"), color = "#000000", method = "lm", se = TRUE, size = 2) +
  geom_smooth(data = subset(Data_1, Cultivar == "Babushka"),color = "#000000", method = "lm", se = FALSE, size = 1) +
  geom_smooth(data = subset(Data_1, Cultivar == "RGT"), color = "#FFFFFF", method = "lm", se = FALSE, size = 1) +
  
  geom_point(data = subset(Data_1, Cultivar == "Babushka"), shape = 16, size = 4, color = "#ff791f") +
  geom_point(data = subset(Data_1, Cultivar == "RGT"), shape = 17, size = 4, color = "#000000") +
  geom_point(data = subset(Data_1, Cultivar == "Babushka"), shape = 16, size = 3, color = "#000000") +
  geom_point(data = subset(Data_1, Cultivar == "RGT"), shape = 17, size = 3, color = "#FFFFFF") +
  labs(x = x.axis.label, y = y.axis.label) +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  facet_grid(~ Fertilizer) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.justification = "right",
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 16, color = "black"),
        strip.text.y = element_text(size = 16, color = "black"),
        panel.spacing = unit(1.5, "lines")) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9),
                     labels = c(expression(10^-1),
                                expression(10^-3),
                                expression(10^-5),
                                expression(10^-7),
                                "No")) +
  scale_shape_manual(values = c("Babushka" = 16, "RGT" = 17))

# Display the plot
plot




#Total shoot N

MODEL_N_lm <-lm(log(Total_N_Shoot) ~ Cultivar*Fertilizer*Microbiome, data = Data_1)

summary(MODEL_N_lm)

residuals <- resid(MODEL_N_lm)

#produce residual vs. fitted plot
plot(fitted(MODEL_N_lm), residuals)
#add a horizontal line at 0 
abline(0,0)

qqnorm(residuals)
abline(0,sd(residuals))

joint_tests(MODEL_N_lm)







Data_Chi <- subset(Data_1, Fertilizer=="Chitin")
Data_Min <- subset(Data_1, Fertilizer=="Mineral")


#Linear model testing the effects of microbiome dilution on combined plants in treatment chitin

lm_chi_combined <- lm(Total_N_Shoot ~ new_microbiome, data=Data_Chi)

summary(lm_chi_combined)



#Graphs Total N Shoot

#PLOT WITH GEOM loess and LM

# Define axis labels
x.axis.label <- "Microbiome dilution"
y.axis.label <- "Total shoot nitrogen (mg)"

Total_N_loess <-ggplot(Data_1, aes(x = new_microbiome, y = Total_N_Shoot)) +
  geom_point(data = subset(Data_1, Cultivar == "Babushka"), shape = 16, size = 4, color = "#ff791f") +
  geom_point(data = subset(Data_1, Cultivar == "RGT"), shape = 17, size = 4, color = "#000000") +
  geom_point(data = subset(Data_1, Cultivar == "Babushka"), shape = 16, size = 3, color = "#000000") +
  geom_point(data = subset(Data_1, Cultivar == "RGT"), shape = 17, size = 3, color = "#FFFFFF") +
  
  geom_smooth(data = subset(Data_1, Cultivar == "Babushka"),color = "#ff791f", method = "lm", se = FALSE, linewidth = 1) +
  geom_smooth(data = subset(Data_1, Cultivar == "Babushka"),color = "#000000", method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.5) +
  geom_smooth(data = subset(Data_1, Cultivar == "RGT"),color = "#000000", method = "lm", se = FALSE, linewidth = 1) +
  geom_smooth(data = subset(Data_1, Cultivar == "RGT"),color = "#FFFFFF", method = "lm", se = FALSE, linetype = "dashed", linewidth = 0.5) +
  facet_grid(~ Fertilizer) +
    
    
  geom_smooth(data = subset(Data_Chi, Cultivar == "Babushka"),color = "#ff791f", method = "loess", se = FALSE, linewidth = 2) +
  geom_smooth(data = subset(Data_Chi, Cultivar == "Babushka"),color = "#000000", method = "loess", se = FALSE, linewidth = 1) +
  geom_smooth(data = subset(Data_Chi, Cultivar == "RGT"),color = "#000000", method = "loess", se = FALSE, linewidth = 2) +
  geom_smooth(data = subset(Data_Chi, Cultivar == "RGT"),color = "#FFFFFF", method = "loess", se = FALSE, linewidth = 1) +
  
  labs(x = x.axis.label, y = y.axis.label) +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  facet_grid(~ Fertilizer) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.justification = "right",
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 16, color = "black"),
        strip.text.y = element_text(size = 16, color = "black"),
        panel.spacing = unit(1.5, "lines")) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9),
                     labels = c(expression(10^-1),
                                expression(10^-3),
                                expression(10^-5),
                                expression(10^-7),
                                "No")) +
  scale_shape_manual(values = c("Babushka" = 16, "RGT" = 17))+
  coord_cartesian(ylim = c(0, max(6)))  # Adjust multiplier as needed

# Display the plot
Total_N_loess

ggsave('Total_N_chi_loess_lm.pdf',Total_N_loess,width=8,height=6)
ggsave('Total_N_chi_loess_lm.png',Total_N_loess,width=8,height=6,type="cairo")


#Regression for RGT in chitin

lm_chi_RGT <- lm(Total_N_Shoot ~ new_microbiome, subset(Data_Chi, Cultivar == "RGT"))

summary(lm_min_RGT)

#Regression for Babushka in chitin

lm_chi_Bab <- lm(Total_N_Shoot ~ new_microbiome, subset(Data_Chi, Cultivar == "Babushka"))

summary(lm_min_Bab)

#Regression for RGT in min

lm_min_RGT <- lm(Total_N_Shoot ~ new_microbiome, subset(Data_Min, Cultivar == "RGT"))

summary(lm_min_RGT)

#Regression for Babushka in min

lm_min_Bab <- lm(Total_N_Shoot ~ new_microbiome, subset(Data_Min, Cultivar == "Babushka"))

summary(lm_min_Bab)



#Regression for both cultivars in chitin

lm_chi_both <- lm(Total_N_Shoot ~ new_microbiome, data = Data_Chi)

summary(lm_chi_both)


###NO3NH4###

NO3NH4_bulk_sub <-subset(NO3_NH4_N_bulk, Degree!="2" & Degree!="4" & Degree!="6")

library(dplyr)

# Rename Degree to Microbiome 
NO3NH4_bulk_sub <- NO3NH4_bulk_sub %>%
  rename(Microbiome = Dilution)

# Rename 
NO3NH4_bulk_sub <- NO3NH4_bulk_sub %>%
  mutate(Microbiome = case_when(
    Microbiome == "a" ~ "a",
    Microbiome == "c" ~ "b",
    Microbiome == "e" ~ "c",
    Microbiome == "g" ~ "d",
    Microbiome == "h" ~ "f",
    TRUE ~ Microbiome  # If none of the conditions match, keep the original value
  ))


#Linear model testing the effect of microbiome on NH4 in treatments without plants.

NH4_bulk <- lm(Percent_of_total_N_NO3NH4 ~ Degree, data=NO3NH4_bulk_sub)

summary(NH4_bulk)


residuals <- resid(NH4_bulk)

#produce residual vs. fitted plot
plot(fitted(NH4_bulk), residuals)
#add a horizontal line at 0 
abline(0,0)

qqnorm(residuals)
abline(0,sd(residuals))

joint_tests(NH4_bulk)


###GRAPHS###

library(ggplot2)

# Define axis labels
x.axis.label <- "Microbiome dilution"
y.axis.label <- expression(bold("% of chitin-N mineralized"))


NO3NH4_per <- ggplot(data = NO3NH4_bulk_sub, aes(x = Degree, y = Percent_of_total_N_NO3NH4)) +
  geom_point(data = NO3NH4_bulk_sub, shape = 16, size = 3, color = "#000000") +
  geom_point(data = NO3NH4_bulk_sub, shape = 16, size = 2, color = "#FFFFFF") +
  
  geom_smooth(data = NO3NH4_bulk_sub,color = "#000000", method = "loess", se = TRUE, linewidth = 2) +
  geom_smooth(data = NO3NH4_bulk_sub, color = "#FFFFFF", method = "loess", se = FALSE, linewidth = 1) +
  
  geom_smooth(data = NO3NH4_bulk_sub, color = "#ff0000", method = "lm", se = FALSE, linewidth = 1, linetype="dashed") +
  
  labs(x = x.axis.label, y = y.axis.label) +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.justification = "right",
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 16, color = "black"),
        strip.text.y = element_text(size = 16, color = "black"),
        panel.spacing = unit(1.5, "lines")) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9),
                     labels = c(expression(10^-1),
                                expression(10^-3),
                                expression(10^-5),
                                expression(10^-7),
                                "No")) 
  
# Display the plot
  NO3NH4_per

ggsave('NO3NH4_percent_of_chitin_N.pdf',NO3NH4_per,width=6,height=5)
ggsave('NO3NH4_percent_of_chitin_N.png',NO3NH4_per,width=6,height=5,type="cairo")


#PLOT Height day 39 WITH GEOM SMOOTH (LM)

# Rename 
Data_2 <- Data_2 %>%
  mutate(Degree = case_when(
    Degree == 1 ~ 1,
    Degree == 3 ~ 3,
    Degree == 5 ~ 5,
    Degree == 7 ~ 7,
    Degree == 10 ~ 9,
    TRUE ~ Degree  # If none of the conditions match, keep the original value
  ))

Data_2 <- Data_2 %>%
  mutate(Substrate = case_when(
    Substrate == "Chitin + min." ~ "Chitin",
    Substrate == "Full min." ~ "Mineral",
   TRUE ~ Substrate  # If none of the conditions match, keep the original value
  ))


library(ggplot2)

# Define axis labels
x.axis.label <- "Microbiome dilution"
y.axis.label <- "Shoot height (cm)"

# Create the plot using facet grid
(plot <- ggplot(Data_2, aes(Degree, Height, shape = Cultivar)) +
  geom_smooth(data = subset(Data_2, Cultivar == "Babushka"),color = "#ff791f", method = "lm", se = TRUE, linewidth = 2) +
  geom_smooth(data = subset(Data_2, Cultivar == "RGT"), color = "#000000", method = "lm", se = TRUE, linewidth = 2) +
  geom_smooth(data = subset(Data_2, Cultivar == "Babushka"),color = "#000000", method = "lm", se = FALSE, linewidth = 1) +
  geom_smooth(data = subset(Data_2, Cultivar == "RGT"), color = "#FFFFFF", method = "lm", se = FALSE, linewidth = 1) +
  
  geom_point(data = subset(Data_2, Cultivar == "Babushka"), shape = 16, size = 4, color = "#ff791f") +
  geom_point(data = subset(Data_2, Cultivar == "RGT"), shape = 17, size = 4, color = "#000000") +
  geom_point(data = subset(Data_2, Cultivar == "Babushka"), shape = 16, size = 3, color = "#000000") +
  geom_point(data = subset(Data_2, Cultivar == "RGT"), shape = 17, size = 3, color = "#FFFFFF") +
  labs(x = x.axis.label, y = y.axis.label) +
  theme(axis.title.x = element_text(face = "bold", size = 12),
        axis.title.y = element_text(face = "bold", size = 12)) +
  facet_grid(~ Substrate) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.justification = "right",
        axis.text.x = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 16, color = "black", face = "bold"),
        axis.title.y = element_text(size = 16, color = "black", face = "bold"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 16, color = "black"),
        strip.text.y = element_text(size = 16, color = "black"),
        panel.spacing = unit(1.5, "lines")) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9),
                     labels = c(expression(10^-1),
                                expression(10^-3),
                                expression(10^-5),
                                expression(10^-7),
                                "No")) +
  scale_shape_manual(values = c("Babushka" = 16, "RGT" = 17)))

# Display the plot
plot


#Try making bar chart of the height at day 39 still using facet grid = Substrate

# Rename 
Data_1 <- Data_1 %>%
  mutate(Degree = case_when(
    Degree == 1 ~ 1,
    Degree == 3 ~ 3,
    Degree == 5 ~ 5,
    Degree == 7 ~ 7,
    Degree == 10 ~ 9,
    TRUE ~ Degree  # If none of the conditions match, keep the original value
  ))



Data_1 <- Data_1 %>%
  mutate(Substrate = case_when(
    Substrate == "Chitin + min." ~ "Chitin",
    Substrate == "Full min." ~ "Mineral",
    TRUE ~ Substrate  # If none of the conditions match, keep the original value
  ))


Data_1_na <- na.omit(Data_1)

#Calculte average and SE

sum_Data_1 <- Data_1_na %>%
  group_by(Substrate, Dilution, Day, Cultivar) %>%
  summarise( 
    n=n(),
    mean=mean(Height),
    sd=sd(Height)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

sum_Data_1 <- sum_Data_1 %>%
  mutate(Day = as.factor(Day))

########BAR PLOT Dw########## 
(BC_height <- ggplot(data=sum_Data_1, aes(x=Day, y=mean, fill=Dilution)) +
  ggtitle(element_blank()) +
  geom_bar(position = position_dodge(0.85), width=0.85, stat="identity", color="black") +
  ylab(expression(bold(Shoot~height~(cm)))) +
  scale_fill_manual(
    values=c("#DDAA33", "#4CDFD2", "#004488", "#FF6600", "#AA4499"),
    labels=c(a = expression(10^-1), b = expression(10^-3), c = expression(10^-5), d = expression(10^-7), f = "No")
  ) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width = 0.85), stat="identity", width=0.4, colour="black", alpha=1, linewidth=0.5) +
  facet_grid(Cultivar~Substrate) +
  theme_minimal() +
  
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.title.x = element_text(size = 20, color = "black", face = "bold"),
        axis.title.y = element_text(size = 20, color = "black", face = "bold"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_blank(),
        strip.text.x = element_text(size = 18, color = "black"),
        strip.text.y = element_text(size = 18, color = "black")) +
  
  guides(fill = guide_legend(nrow = 1)))


ggsave('BC_height.pdf',BC_height,width=12,height=8)
ggsave('BC_height.png',BC_height,width=12,height=8,type="cairo")


# Create the plot with logarithmic trendlines, non-overlapping points, and facets for Cultivar and Fertilizer
x.axis.label <- "Day" # Define the label for the x-axis
y.axis.label <- "Growth rate (cm)" # Define the label for the y-axis

# Create custom labeller function to change facet labels
custom_labeller <- labeller(
  Substrate = c("Chitin + min." = "Chitin", "Full min." = "Mineral")
)

# Convert Degree to a factor
Data_3$Degree <- as.factor(Data_3$Degree)

Rate_plot <- ggplot(Data_3, aes(x = Day, y = Rate)) +
  geom_smooth(aes(colour = Degree, fill = Degree), method = "loess", formula = y ~ x, se = FALSE, fill = "#DDDDDD", linewidth=1) +
  geom_point(aes(colour = Degree, fill = Degree), size = 3, position = "jitter") +  # Add points with the same colors and jitter them
  labs(x = x.axis.label, y = y.axis.label) +
  facet_grid(Cultivar ~ Substrate, labeller = custom_labeller) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.title.x = element_text(size = 20, color = "black", face = "bold"),
        axis.title.y = element_text(size = 20, color = "black", face = "bold"),
        legend.text = element_text(size = 12, color = "black"),
        legend.title = element_text(size = 14, color = "black", face = "bold"),
        strip.text.x = element_text(size = 18, color = "black"),
        strip.text.y = element_text(size = 18, color = "black"))

# Customize the line colors and legend labels in the legend
Rate_plot <- Rate_plot +
  scale_color_manual(
    values = c("#DDAA33", "#4CDFD2", "#004488", "#FF6600", "#AA4499"),
    labels = expression(paste("10"^-1, " 10"^-3, " 10"^-5, " 10"^-7, " 10"^-10)),
    name = "Microbiome dilution",
    guide = guide_legend(override.aes = list(fill = c("#DDAA33", "#4CDFD2", "#004488", "#FF6600", "#AA4499")))
  )

# Display the plot
print(Height_plot)

ggsave('Height.pdf',Height_plot,width=8,height=6)
ggsave('Height.png',Height_plot,width=8,height=6,type="cairo")





