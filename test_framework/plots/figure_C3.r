library(tidyverse)
library(ggpubr)
library(svglite)

setwd("~/bioinformatics/hmm/git_tc_hmmmlib/test_framework/plots")

data = read_csv("../running_times/running_time_C3.csv", col_names = T)
#names(data) = c('test', 'observations', 'time', 'algorithm', 'variant', 'iterations')
#data_statespace = data %>% filter(test == "statespace")
data$algorithm = data %>% pull(algorithm) %>% recode(forward = "Forward",
                                                     backward_time = "Backward",
                                                     baumWelch = "Baum-Welch (1 iteration)",
                                                     posteriorDecoding = "Posterior Decoding",
                                                     viterbi = "Viterbi")

caption = "error bars: standard deviation of 3 replicates\ninput size: 100000, alphabet size: 4"

data_statespace_grouped = data %>%
    group_by(observations, algorithm, variant) %>%
    summarise(mean = mean(time), sd = sd(time))

data_statespace_grouped %>% ggplot(aes(observations, mean, color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), size = 0.3, alpha = .65) +
    facet_wrap(. ~ algorithm, scales = "free") +
    labs(
        x = "hidden states",
        y = "mean time [s]",
        caption = caption,
        title = "Running time of increasing hidden state space"
    )
ggsave("figure_C3_raw.pdf", height = 5, width = 9)

data_statespace_grouped %>% ggplot(aes(observations, mean / (observations ^ 2), color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(
        aes(
            ymin = (mean - sd) / observations ^ 2,
            ymax = (mean + sd) / observations ^ 2
        ),
        size = 0.3,
        alpha = .65
    ) +
    facet_wrap(. ~ algorithm, scales = "free") +
    labs(
        x = "hidden states",
        y = "mean time [s] scaled",
        caption = caption,
        title = "Running time of increasing hidden state space"
    )
ggsave("figure_C3_scaled.pdf",
       height = 5,
       width = 9)

