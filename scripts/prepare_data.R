###-----------------------------------------------------------------------###
###   Script to find minimum degree m of Bernstein Polynomials            ###
###   Author: Juliana Freitas de Mello e Silva                            ###
###                                                                       ###
###   Reading R database "growth". This data base contains longitudinal   ###
###   measurements for boys and girls. Since we are interested in         ###
###   working with growth, we calculate the adjacent differences which,   ###
###   in some cases, resulted in negative values.                         ###
###-----------------------------------------------------------------------###

###-------------------------------
###   Removing objects
###-------------------------------

rm(list = ls())

###-------------------------------
###   Setting seed
###-------------------------------

set.seed(123456789)

###-------------------------------
###   Requiring packages
###-------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)
library(sitar)

###-------------------------------
###   Working directories
###-------------------------------

directory <- paste0(getwd(), "\\")
data_directory <- paste0(directory, "data\\")
figures_directory <- paste0(directory, "figures\\")
tables_directory <- paste0(directory, "tables\\")

###-------------------------------###-------------------------------###


###------------------------------###------------------------------###

### Graphical definitions
colour_plot <- c("#999999", "#F0E442", "#E69F00", "#D55E00", "#56B4E9", "#0072B2", "#009E73", "#CC79A7")
size_plot <- 12
vjust_plot <- 1.8
theme_plot <- theme_bw() +
  theme(axis.text.x = element_text(size = size_plot),
        axis.text.y = element_text(size = size_plot),
        axis.title.x = element_text(size = size_plot, vjust = -vjust_plot), 
        axis.title.y = element_text(size = size_plot, vjust = -vjust_plot), 
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size = size_plot),
        legend.title = element_text(size = size_plot),
        legend.position = "bottom",
        plot.margin = unit(x = c(1, 1, 1, 0.5), units = "cm"),
        plot.title = element_text(hjust = 0.5)
  )

###------------------------------###------------------------------###

###-------------------------------
###   Important settings for this analysis
###-------------------------------



###------------------------------###------------------------------###

###-------------------------------
###   Calling functions
###-------------------------------



###------------------------------###------------------------------###

###-------------------------------
###   Reading data
###-------------------------------

data("berkeley")
head(berkeley)

### Removing NAs
data <- berkeley[which(!is.na(berkeley$height)),]

### Creating growth variable
data <- data %>%
  group_by(id) %>%
  mutate(growth = c(0, diff(height))) %>%
  ungroup()

### Creating data frame with both observed mean and sd variables
df_summary <- data %>%
  group_by(age, sex) %>%
  summarise(obs_mean_growth = mean(growth, na.rm = TRUE),
            obs_sd_growth = sd(growth, na.rm = TRUE))

###------------------------------###------------------------------###

###------------------------------
###   Observed trajectories
###------------------------------

ggarrange(
  
  ### Boys
  ggplot(data = subset(data, data$sex %in% 1)) +
    geom_line(mapping = aes(x = age, y = growth, group = id), colour = colour_plot[1], linewidth = 0.3) +
    geom_line(mapping = aes(x = age, y = obs_mean_growth), colour = "black", linetype = 2,
              linewidth = 1,
              data = subset(df_summary, df_summary$sex %in% 1)) +
    ggtitle("Boys") +
    scale_x_continuous(name = "Age (in years)", breaks = seq(from = 0, to = 21, by = 2.5)) +
    scale_y_continuous(name = "Growth (in cm)\n", breaks = seq(from = -2, to = 14, by = 2),
                       limits = c(-0.5, 14.0)) +
    theme_plot,
  
  ### Girls
  ggplot(data = subset(data, data$sex %in% 2)) +
    geom_line(mapping = aes(x = age, y = growth, group = id), colour = colour_plot[1], linewidth = 0.3) +
    geom_line(mapping = aes(x = age, y = obs_mean_growth), colour = "black", linetype = 2,
              linewidth = 1,
              data = subset(df_summary, df_summary$sex %in% 2)) +
    ggtitle("Girls") +
    scale_x_continuous(name = "Age (in years)", breaks = seq(from = 0, to = 21, by = 2.5)) +
    scale_y_continuous(name = "Growth (in cm)\n", breaks = seq(from = -2, to = 14, by = 2),
                       limits = c(-0.5, 14.0)) +
    theme_plot,
  
  ncol = 2, nrow = 1)


###------------------------------###------------------------------###

###------------------------------
###   Standard deviation along time
###------------------------------

ggplot(data = df_summary) +
  geom_line(mapping = aes(x = age, y = obs_sd_growth, colour = sex), linewidth = 1) +
  scale_colour_manual(name = "", values = colour_plot[7:8], labels = c("males", "females")) +
  scale_x_continuous(name = "Age (in years)", breaks = seq(from = 0, to = 21, by = 2.5)) +
  scale_y_continuous(name = "Standard deviation (in cm)\n", breaks = seq(from = 0, to = 4, by = 0.5),
                     limits = c(0, 3.6)) +
  theme_plot


###------------------------------###------------------------------###

###------------------------------
###   Exporting data to export for modeling in JAGS
###------------------------------

data_export <- data.frame(id = as.numeric(data$id),
                          age = data$age,
                          gender = data$sex, 
                          height = data$height
)

write.table(x = data_export, file = paste0(data_directory, "data_growth_long.txt"))

