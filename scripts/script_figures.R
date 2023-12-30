###-------------------------------------------------------------###
###   File #5: figures of application to Berkeley Growth Study  ###
###   Longitudinal data - Bernstein Polynomials                 ###
###   Author: Juliana Freitas de Mello e Silva                  ###
###-------------------------------------------------------------###

###--------------------------------
###   Removing objects
###---------------- ----------------

rm(list = ls())

###--------------------------------
###   Setting seed
###--------------------------------

set.seed(123456789)

###--------------------------------
###   Requiring packages
###--------------------------------

library(coda)
library(dplyr)
library(ggplot2)
library(ggpubr)

###--------------------------------
###   Working directories
###--------------------------------

directory <- paste0(getwd(), "\\")
data_directory <- paste0(directory, "data\\data_application\\")
figures_directory <- paste0(directory, "figures\\application\\")
results_directory <- paste0(directory, "results\\results_application\\")
tables_directory <- paste0(directory, "tables\\")

###--------------------------------###--------------------------------###

###--------------------------------
###   Graphical definitions
###--------------------------------

colour_plot <- c("#999999", "#F0E442", "#E69F00", "#D55E00", "#56B4E9", "#0072B2", "#009E73", "#CC79A7")
size_plot <- 18
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
        plot.title = element_text(hjust = 0.5, size = size_plot)
  )

###--------------------------------###--------------------------------###

###-------------------------------
###   Important settings for this analysis
###-------------------------------



###--------------------------------###--------------------------------###

###--------------------------------
###   Calling functions
###--------------------------------

source(paste0(directory, "utils\\", "functions.R"))

###------------------------------###------------------------------###

###-------------------------------
###   Reading data
###-------------------------------

data_long <- read.table(file = paste0(data_directory, "growth_long.txt"))
T_max <- max(data_long$age)


head(data_long)

###------------------------------###------------------------------###

###--------------------------------
###   Important specifications
###--------------------------------

m_boys <- 7
m_girls <- 8

t_plot <- seq(from = min(data_long$age), to = T_max, length.out = 100)

###--------------------------------
###   Observed mean curves
###--------------------------------

df_summary <- data_long %>%
  group_by(age, gender) %>%
  summarise(mu = mean(height, na.rm = TRUE)) %>%
  ungroup()
df_summary$group <- "observed mean"

###--------------------------------
###   Estimated mean curve for boys
###--------------------------------

### BP basis
basis_plot_boys <- sapply(1:m_boys,
                          function(l)(dbeta(x = t_plot / T_max, shape1 = l,
                                            shape2 = m_boys - l + 1) / m_boys)) 

### Posterior sample for mu, boys
posterior_mu_xi_boys <- as.matrix(read.table(file = paste0(results_directory, "m",
                                                           m_boys,
                                                           "/posterior_mu_xi_boys.txt")))

names_mu_xi_boys <-  paste0("mu_xi_boys", 1:m_boys)
colnames(posterior_mu_xi_boys) <- names_mu_xi_boys

### Posterior sample of estimated mean curve, boys
posterior_mean_curve_boys <- basis_plot_boys %*% t(posterior_mu_xi_boys)



###--------------------------------
###   Estimated mean curve for girls
###--------------------------------

### BP basis
basis_plot_girls <- sapply(1:m_girls,
                           function(l)(dbeta(x = t_plot / T_max, shape1 = l,
                                             shape2 = m_girls - l + 1) / m_girls))

### Posterior sample for mu, girls
posterior_mu_xi_girls <- as.matrix(read.table(file = paste0(results_directory, "m",
                                                            m_girls,
                                                            "/posterior_mu_xi_girls.txt")))
names_mu_xi_girls <-  paste0("mu_xi_girls", 1:m_girls)
colnames(posterior_mu_xi_girls) <- names_mu_xi_girls

### Posterior sample of estimated mean curve, girls
posterior_mean_curve_girls <- basis_plot_girls %*% t(posterior_mu_xi_girls)


###--------------------------------
###   Adding estimated mean curves in df_summary
###-------------------------------

df_summary <- df_summary %>%
  add_row(data.frame(age = rep(x = t_plot, times = 2),
                     gender = c(rep(x = 0, times = length(t_plot)),
                                rep(x = 1, times = length(t_plot))),
                     mu = c(apply(posterior_mean_curve_boys, 1, median),
                                    apply(posterior_mean_curve_girls, 1, median)),
                     group = rep(x = "estimated mean", times = 2*length(t_plot)))
  )
  
df_summary <- df_summary %>%
  add_row(data.frame(age = rep(x = t_plot, times = 2),
                     gender = 0,
                     mu = c(as.numeric(HPDinterval(obj = as.mcmc(t(posterior_mean_curve_boys)),
                                                   prob = 0.95)[,"lower"]),
                            as.numeric(HPDinterval(obj = as.mcmc(t(posterior_mean_curve_boys)),
                                                   prob = 0.95)[,"upper"])),
                     group = c(rep(x = "HPD_lower", times = length(t_plot)),
                               rep(x = "HPD_upper", times = length(t_plot)))
  ))

df_summary <- df_summary %>%
  add_row(data.frame(age = rep(x = t_plot, times = 2),
                     gender = 1,
                     mu = c(as.numeric(HPDinterval(obj = as.mcmc(t(posterior_mean_curve_girls)),
                                                   prob = 0.95)[,"lower"]),
                            as.numeric(HPDinterval(obj = as.mcmc(t(posterior_mean_curve_girls)),
                                                   prob = 0.95)[,"upper"])),
                     group = c(rep(x = "HPD_lower", times = length(t_plot)),
                               rep(x = "HPD_upper", times = length(t_plot)))
  ))

###--------------------------------
###   Plots
###--------------------------------

### Observed trajectories, observed mean curve, estimated mean curve

pdf(file = paste0(figures_directory, "observed_trajectories_estimated_observed_mean_curve_boys_m",
                  m_boys, "_girls_m", m_girls, ".pdf"),
    width = 10*2, height = 6)
ggarrange(
  
  ### Boys
  ggplot(data = subset(df_summary, df_summary$gender %in% 0)) +
    geom_line(mapping = aes(x = age, y = height, group = id),
              colour = colour_plot[1], linewidth = 0.6,
              data = data_long[which(data_long$gender %in% 0),]) +
    geom_line(mapping = aes(x = age, y = mu, colour = group, linetype = group), linewidth = 1.2,
              data = subset(df_summary, df_summary$gender %in% 0)) +
    ggtitle(paste0("Boys (m = ", m_boys, ")")) +
    scale_colour_manual(name = "", values = colour_plot[c(4, 4, 4, 6)],
                        labels = c("estimated mean", "HPD_lower", "HPD_upper", "observed mean")) +
    scale_linetype_manual(name = "", values = c(1, 2, 2, 1),
                          labels = c("estimated mean", "HPD_lower", "HPD_upper", "observed mean")) +
    scale_x_continuous(name = "Age (in years)",
                       breaks = seq(from = 0, to = T_max, by = 2.5)) +
    scale_y_continuous(name = "Growth (in cm / year)\n",
                       breaks = seq(from = 0, to = 20, by = 2.5), limits = c(-0.3, 14.0)) +
    theme_plot,
  
  ### Girls
  
  ggplot(data = subset(df_summary, df_summary$gender %in% 1)) +
    geom_line(mapping = aes(x = age, y = height, group = id),
              colour = colour_plot[1], linewidth = 0.6,
              data = data_long[which(data_long$gender %in% 1),]) +
    geom_line(mapping = aes(x = age, y = mu, colour = group, linetype = group), linewidth = 1.2,
              data = subset(df_summary, df_summary$gender %in% 1)) +
    ggtitle(paste0("Girls (m = ", m_girls, ")")) +
    scale_colour_manual(name = "", values = colour_plot[c(4, 4, 4, 6)],
                        labels = c("estimated mean", "HPD_lower", "HPD_upper", "observed mean")) +
    scale_linetype_manual(name = "", values = c(1, 2, 2, 1),
                          labels = c("estimated mean", "HPD_lower", "HPD_upper", "observed mean")) +
    scale_x_continuous(name = "Age (in years)",
                       breaks = seq(from = 0, to = T_max, by = 2.5)) +
    scale_y_continuous(name = "Growth (in cm / year)\n",
                       breaks = seq(from = 0, to = 20, by = 2.5),  limits = c(-0.3, 14.0)) +
    theme_plot,
  
  ncol = 2, nrow = 1, legend = "bottom", common.legend = TRUE
  
  
)

dev.off()


### Observed mean curve, estimated mean curve

pdf(file = paste0(figures_directory, "estimated_observed_mean_curve_boys_m",
                  m_boys, "_girls_m", m_girls, ".pdf"),
    width = 10*2, height = 8)
ggarrange(
  
  ### Boys
  ggplot(data = subset(df_summary, df_summary$gender %in% 0)) +
    geom_line(mapping = aes(x = age, y = mu, colour = group, linetype = group),
              linewidth = 1) +
    scale_colour_manual(name = "", values = colour_plot[c(4, 4, 4, 6)],
                        labels = c("estimated mean", "HPD", "HPD", "observed mean")) +
    scale_linetype_manual(name = "", values = c(1, 2, 2, 1),
                          labels = c("estimated mean", "HPD", "HPD", "observed mean")) +
    scale_x_continuous(name = "Age (in years)",
                       breaks = seq(from = 0, to = T_max, by = 1.5)) +
    scale_y_continuous(name = "Growth (in cm / year)\n",
                       breaks = seq(from = 0, to = 20, by = 2.5), limits = c(-0.3, 8.6)) +
    theme_plot,
  
  ### Girls
  ggplot(data = subset(df_summary, df_summary$gender %in% 1)) +
    geom_line(mapping = aes(x = age, y = mu, colour = group, linetype = group), linewidth = 1) +
    scale_colour_manual(name = "", values = colour_plot[c(4, 4, 4, 6)],
                        labels = c("estimated mean", "HPD_lower", "HPD_upper", "observed mean")) +
    scale_linetype_manual(name = "", values = c(1, 2, 2, 1),
                          labels = c("estimated mean", "HPD_lower", "HPD_upper", "observed mean")) +
    scale_x_continuous(name = "Age (in years)",
                       breaks = seq(from = 0, to = T_max, by = 1.5)) +
    scale_y_continuous(name = "Growth (in cm / year)\n",
                       breaks = seq(from = 0, to = 20, by = 2.5), limits = c(-0.3, 8.6)) +
    theme_plot,
  
  ncol = 2, nrow = 1, legend = "bottom", common.legend = TRUE
  
  
)
dev.off()

### Now you can go to file #4b or to file #6







