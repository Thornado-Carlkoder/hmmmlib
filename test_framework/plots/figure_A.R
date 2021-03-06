library(tidyverse)
library(ggpubr)
library(svglite)

setwd("~/bioinformatics/hmm/git_tc_hmmmlib/test_framework/plots")

# This figure shows that the matrix implementation is faster than the conventional implementation
# All fwbw based algorithms
# scaled by theoretical running time
# full density


# Consider varying statespace and varying alphabet
# Yet i don't know how to show that the iterations in baum-welch are linearly scaled

# The data is produced by running_time_figure_A.py


## Inputsize ##

data = read_csv("../data/A.csv", col_names = F)
names(data) = c('test',
                'observations',
                'time',
                'algorithm',
                'variant',
                'iterations',
                'statespace')

data$algorithm = data %>% pull(algorithm) %>% recode(forward = "Forward",
                                                     backward_time = "Backward",
                                                     baumWelch = "Baum-Welch (1 iteration)",
                                                     posteriorDecoding = "Posterior Decoding")

# Control the order of variables in geom_grid
data$algorithm = factor(data$algorithm, levels=c("Forward", "Backward", "Posterior Decoding", "Baum-Welch (1 iteration)"))


data_input = data %>% filter(test == "inputsize")


data_input_grouped = data_input %>%
    mutate_at(vars(iterations), factor) %>%
    group_by(observations, algorithm, variant, iterations, statespace) %>%
    summarise(mean = mean(time), sd = sd(time))



data_input_grouped %>% filter(iterations == "1" |
                                  is.na(iterations)) %>%
    filter(variant != "Conventional sparse") %>%
    ggplot(aes(observations, mean / observations, color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = (mean - sd)/observations, ymax = (mean + sd)/observations), size = 0.3, alpha = .65) +
    facet_grid(statespace ~ algorithm, scales = "free") +
    labs(y = "mean time [s] scaled",
         caption = "error bars: standard deviation of 5 replicates\nalphabet size: 4",
         title = "Running time (scaled) for increasing input size",
         subtitle = "Linear algebra based implementations are faster.") +
    geom_hline(yintercept = 0, alpha = 0)
ggsave("figure_A.pdf", height = 6, width = 10)











