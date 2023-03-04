###-------------------------------------------------###
###   Degree selection - figure os basis functions  ###
###   Author: Juliana Freitasd de Mello e Silva     ###
###-------------------------------------------------###

### Removing objects
rm(list = ls())

### Setting seed
set.seed(123456789)

### Loading packages
library(ggplot2)
library(scales)

### Directories
directory <- "C:\\Users\\julia\\Desktop\\Projetos\\Doutorado\\artigo_degreeselection\\"
results_directory <- paste0(directory, "figures\\")


### Graphical definitions
colour_plot <- c("#999999", "#F0E442", "#E69F00", "#D55E00", "#56B4E9", "#0072B2", "#009E73", "#CC79A7")
size_plot <- 25
vjust_plot <- 1.8
theme_plot <- theme_bw() + # white backgroud
  theme(axis.text.x = element_text(size = size_plot), # optional settings for x-axis
        axis.text.y = element_text(size = size_plot), # optional settings for y-axis
        axis.title.x = element_text(size = size_plot, vjust = -vjust_plot), # optional settings for x-axis
        axis.title.y = element_text(size = size_plot, vjust = vjust_plot), # optional settings for y-axis
        axis.title.y.right = element_text(size = size_plot, vjust = vjust_plot), # optional settings for y-axis
        legend.key.size = unit(1, "line"),
        legend.text = element_text(size = size_plot),
        legend.title = element_text(size = size_plot),
        plot.margin = unit(x = c(1, 1, 1, 0.5), units = "cm"), # optional settings for margins
        plot.title = element_text(hjust = 0.5)
        # text = element_text(family = "NimbusMon", size = size_plot)
  )

###------------------------------###------------------------------###


m <- 4 # therefore, BP degree is m-1

t <- seq(from = 0, to = 1, length.out = 100)
y <- c(sapply(1:m, function(k)(dbinom(x = k-1, size = m-1, prob = t))))
basis <- factor(x = c(sapply(1:m, function(k)(rep(x = paste0("k = ", k), times = length(t))))),
                levels = paste0("k = ", 1:m), 
                labels = paste0("k=", 1:m))

df_plot <- data.frame(t = t,
                      y = y,
                      basis = basis)

pdf(file = paste0(results_directory, "bernstein_basis_m", m, "_thick_points.pdf"), width = 10, height = 6)
ggplot(data = df_plot) +
  geom_line(mapping = aes(x = t, y = y, group = basis, colour = basis), size = 1.2) +
  geom_vline(xintercept = sapply(1:m, function(k)((k-1) / (m-1))), colour = grey(0.5), linetype = 2, size = 0.8) +
  scale_x_continuous(breaks = ((1:m) - 1) / (m-1), labels = round(x = ((1:m) - 1) / (m-1), digits = 2)) +
  xlab("t") + ylab(expression("b"["k,m-1" ](t))) +
  theme_plot
dev.off()







