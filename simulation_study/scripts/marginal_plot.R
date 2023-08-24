##############
### Script ###
##############

# Script to get marginal posterior of gamma under prior v0=0 plots across range of spike variances in backwards DPE procedure.
# (Code for Simulation Study)

#################
### Libraries ###
#################

library(ggplot2)
library(latex2exp)

############
### Plot ###
############

# path: path for file containing marginal posterior of gamma under prior v0=0 values
path_marginal = ''
# Range of spike variances DPE was evaluated over.
v0 = exp(seq(-20,-1, length.out=15))
# Marginal posterior of gamma values.
marg_back = as.numeric(unlist(read.csv(path_marginal)[,2]))

# Transform to spike variances to log-range.
x = log(v0)
# Combine log spike variance values and marginal values in a dataframe.
df = data.frame(cbind(x, marg_back))
colnames(df) = c('log_v0', 'marginal')
df$log_v0 <- as.numeric(df$log_v0)
df$marginal <- as.numeric(df$marginal)

# Plot marginal posterior of gamma under prior v0=0.
ggplot(df, aes(x=log_v0, y=marginal)) +
  geom_line(size=1.5) + 
  geom_point(size=3) + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        text=element_text(size=40) ) +
  scale_color_manual(values=c('black'))+
  xlab(TeX('$log(\\nu_0)$')) + 
  ylab(TeX('$log(p(\\gamma | y)$')) 

