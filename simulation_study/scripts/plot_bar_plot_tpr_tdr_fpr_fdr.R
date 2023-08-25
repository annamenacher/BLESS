##############
### Script ###
##############

# Script to plot the barplots for evaluating the inference results, such as TPR, TDR, FPR, and FDR, for 
# simulation study setting that vary in sample size (N=500, N=1,000, N=5,000) and base rate intensities
# (lambda=1, lambda=2, lambda=3).

#################
### Libraries ###
#################

# Libraries
library(ggplot2)
library(grid)
library(latex2exp)

#################
### Constants ###
#################

# Define output path.
path = ''

# Set directory to output path.
setwd(path)

# TPR
png("TPR.png", width = 10, height = 10, units = 'in', res = 300)

value = c(0.6635, 0.9080, 0.9770, 
          0.9458, 0.9983, 1.0000,
          1.0000, 1.0000, 1.0000,
          0.9972, 1.0000, 1.0000, 
          1.0000, 1.0000, 1.0000,
          1.0000, 1.0000, 1.0000,
          0.6600, 0.9629, 0.9963,
          0.9738, 0.9999, 1.0000,
          1.0000, 1.0000, 1.0000)
lambda = rep(c('lambda=1', 'lambda=2', 'lambda=3'), 9*3)
N = rep(c(c('N1', 'N1', 'N1'), c('N2', 'N2', 'N2'), c('N3', 'N3', 'N3')), 3)
model = c(rep('BLESS', 9), rep('BSGLMM', 9), rep('Firth', 9))

data = data.frame(value = value, lambda = lambda, N = N, model = model)

g1 = ggplot(data = data, aes(x = interaction(lambda, N), y = value, fill = factor(model))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_cartesian(ylim = c(0, 1)) +
  annotate("text", x = 1:9, y = -0.1,
           label = rep(c(TeX(r"($\lambda$=1)"), TeX(r"($\lambda$=2)"), TeX(r"($\lambda$=3)")), 3), size = 10) +
  annotate("text", c(2, 5, 8), y = -0.18, label = c("N=500", "N=1,000", "N=5,000"), size = 10) +
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 6.5, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", colour = "white"),
        legend.position = "none",
        plot.title = element_text(size = 40, face = "bold")) +
  ylab('') + 
  guides(fill = guide_legend(title = "Model"))

g2 <- ggplot_gtable(ggplot_build(g1))
g2$layout$clip[g2$layout$name == "panel"] <- "off"
grid.draw(g2)
dev.off()

# TDR
#pdf("TDR.pdf", width=10, height=10)
png("TDR.png", width = 10, height = 10, units = 'in', res = 300)

value = c(0.9520, 0.9808, 0.9893,
          0.9713, 0.9882, 0.9936,
          0.9873, 0.9973, 0.9968,
          0.8829, 0.8941, 0.8949 , 
          0.8906,  0.8993, 0.9092 ,
          0.9516,  0.9548, 0.9515 ,
          0.8911, 0.9368, 0.9489,
          0.9406, 0.9607, 0.9644,
          0.9702, 0.9743, 0.9716)
lambda = rep(c('lambda=1', 'lambda=2', 'lambda=3'), 9*3)
N = rep(c(c('N1', 'N1', 'N1'), c('N2', 'N2', 'N2'), c('N3', 'N3', 'N3')), 3)
model = c(rep('BLESS', 9), rep('BSGLMM', 9), rep('Firth', 9))

data = data.frame(value = value, lambda = lambda, N = N, model = model)

g1 = ggplot(data = data, aes(x = interaction(lambda, N), y = value, fill = factor(model))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_cartesian(ylim = c(0, 1)) +
  annotate("text", x = 1:9, y = -0.1,
           label = rep(c(TeX(r"($\lambda$=1)"), TeX(r"($\lambda$=2)"), TeX(r"($\lambda$=3)")), 3), size = 10) +
  annotate("text", c(2, 5, 8), y = -0.18, label = c("N=500", "N=1,000", "N=5,000"), size = 10) +
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 6.5, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20),
        axis.text.y = element_text(size = 30), 
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", colour = "white"),
        legend.position = "none",
        plot.title = element_text(size = 40, face = "bold")) +
  ylab('') + 
  guides(fill = guide_legend(title = "Model"))

g2 <- ggplot_gtable(ggplot_build(g1))
g2$layout$clip[g2$layout$name == "panel"] <- "off"
grid.draw(g2)
dev.off()

# FPR
#pdf("FPR.pdf", width=10, height=10)
png("FPR.png", width = 10, height = 10, units = 'in', res = 300)

value = c(0.0334, 0.0178, 0.0106, 
          0.0282, 0.0120, 0.0065,
          0.0129, 0.0027, 0.0032,
          0.1340, 0.1196, 0.1188, 
          0.1243, 0.1132, 0.1009,
          0.1019, 0.0815, 0.0790,
          0.0802, 0.0652, 0.0540,
          0.0619, 0.0413, 0.0372,
          0.0309, 0.0265, 0.0293)
lambda = rep(c('lambda=1', 'lambda=2', 'lambda=3'), 9*3)
N = rep(c(c('N1', 'N1', 'N1'), c('N2', 'N2', 'N2'), c('N3', 'N3', 'N3')), 3)
model = c(rep('BLESS', 9), rep('BSGLMM', 9), rep('Firth', 9))

data = data.frame(value = value, lambda = lambda, N = N, model = model)

g1 = ggplot(data = data, aes(x = interaction(lambda, N), y = value, fill = factor(model))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_cartesian(ylim = c(0, 0.2), expand = T) +
  annotate("text", x = 1:9, y =  -0.02,
           label = rep(c(TeX(r"($\lambda$=1)"), TeX(r"($\lambda$=2)"), TeX(r"($\lambda$=3)")), 3), size = 10) +
  annotate("text", c(2, 5, 8), y = -0.035, label = c("N=500", "N=1,000", "N=5,000"), size = 10) +
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 6.5, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size = 20),
        axis.text.y = element_text(size = 30),
        plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill = "white", colour = "white"),
        legend.position = "none",
        plot.title = element_text(size = 40, face = "bold")) +
  ylab('') + 
  guides(fill = guide_legend(title = "Model"))

g2 <- ggplot_gtable(ggplot_build(g1))
g2$layout$clip[g2$layout$name == "panel"] <- "off"
grid.draw(g2)
dev.off()

# FDR
png("FDR.png", width = 10, height = 10, units = 'in', res = 300)

value = c(0.0480, 0.0192, 0.0107, 
          0.0287, 0.0118, 0.0064,
          0.0127, 0.0027, 0.0032,
          0.1171, 0.1059, 0.1051, 
          0.1094, 0.1007, 0.0908 ,
          0.0916, 0.0748, 0.0728,
          0.1089, 0.0632, 0.0511,
          0.0594, 0.0393, 0.0356,
          0.0298, 0.0257, 0.0284)
lambda = rep(c('lambda=1', 'lambda=2', 'lambda=3'), 9*3)
N = rep(c(c('N1', 'N1', 'N1'), c('N2', 'N2', 'N2'), c('N3', 'N3', 'N3')), 3)
model = c(rep('BLESS', 9), rep('BSGLMM', 9), rep('Firth', 9))

data = data.frame(value = value, lambda = lambda, N = N, model = model)

g1 = ggplot(data = data, aes(x = interaction(lambda, N), y = value, fill = factor(model))) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_cartesian(ylim = c(0, 0.2), expand = T) +
  annotate("text", x = 1:9, y =  -0.02,
           label = rep(c(TeX(r"($\lambda$=1)"), TeX(r"($\lambda$=2)"), TeX(r"($\lambda$=3)")), 3), size = 10) +
  annotate("text", c(2, 5, 8), y = -0.035, label = c("N=500", "N=1,000", "N=5,000"), size = 10) +
  theme_classic() +
  theme(plot.margin = unit(c(1, 1, 6.5, 1), "lines"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 30),
        text = element_text(size = 20),
        legend.position=c(0.85, 0.8),
        legend.text = element_text(size = 30),
        legend.title = element_text(size = 30),
        legend.key.height = unit(1.5, 'cm'),
        legend.key.width = unit(1.5, 'cm'),
        plot.title = element_text(size = 40, face = "bold")) +
  ylab('') + 
  guides(fill=guide_legend(title="Model"))

g2 <- ggplot_gtable(ggplot_build(g1))
g2$layout$clip[g2$layout$name == "panel"] <- "off"
grid.draw(g2)
dev.off()




