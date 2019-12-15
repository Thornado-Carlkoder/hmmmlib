rm(list = ls())
library(tidyverse)
library(ggpubr)
library(svglite)

figurewidth = 8
setwd("~/bioinformatics/hmm/git_tc_hmmmlib/test_framework/plots")

# This figure shows that the matrix implementation is faster than the conventional implementation
# All fwbw based algorithms
# scaled by theoretical running time
# full density

state = "1"

# Consider varying statespace and varying alphabet
# Yet i don't know how to show that the iterations in baum-welch are linearly scaled

# The data is produced by running_time_figure_A.py


## Inputsize ##

data = read_csv(paste0("../running_times/running_time_A", state, ".csv"), col_names = T, comment = '#')
#names(data) = c('test', 'observations', 'time', 'algorithm', 'variant', 'iterations', 'statespace')

data$algorithm = data %>% pull(algorithm) %>% recode(forward = "Forward",
                                                     backward_time = "Backward",
                                                     baumWelch = "Baum-Welch (1 iteration)",
                                                     posteriorDecoding = "Posterior Decoding",
                                                     viterbi = "Viterbi")

# Control the order of variables in geom_grid
data$algorithm = factor(data$algorithm, levels=c("Forward", "Backward", "Posterior Decoding", "Baum-Welch (1 iteration)", "Viterbi"))


data_input = data %>% filter(test == "inputsize")


data_input_grouped = data_input %>%
    #mutate_at(vars(iterations), factor) %>%
    group_by(observations, algorithm, variant, statespace) %>%
    summarise(mean = mean(time), sd = sd(time))



data_input_grouped %>% 
    filter(variant != "Conventional sparse") %>%
    ggplot(aes(observations, mean / observations, color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = (mean - sd)/observations, ymax = (mean + sd)/observations), size = 0.3, alpha = .65) +
    facet_grid(statespace ~ algorithm, scales = "free") +
    labs(y = "mean time [s] scaled"
         #caption = "error bars: standard deviation of 3 replicates\nalphabet size: 4"
         #title = "Running time (scaled) for increasing input size",
         #subtitle = "Linear algebra based implementations are faster."
         ) +
    geom_hline(yintercept = 0, alpha = 0) + 
    theme_light()
ggsave(paste0("pdf/figure_A", state, ".pdf"), width = figurewidth)
ggsave(paste0("svg/figure_A", state, ".svg"), width = figurewidth)













