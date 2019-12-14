library(tidyverse)
library(ggpubr)
library(svglite)

setwd("~/bioinformatics/hmm/git_tc_hmmmlib/test_framework/plots")

data = read_csv("../running_times/running_time_C4.csv", col_names = T) %>% rename(iterations = observations)
#names(data) = c('test', 'iterations', 'time', 'algorithm', 'variant', 'iterations')


caption = "error bars: standard deviation of 3 replicates\nhidden states: 10, alphabet size: 4,input size: 100000"

data_grouped = data %>%
    group_by(iterations, algorithm, variant) %>%
    summarise(mean = mean(time), sd = sd(time))

data_grouped %>% ggplot(aes(iterations, mean, color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), size = 0.3, alpha = .65) +
    #facet_wrap(. ~ algorithm, scales = "free") +
    labs(
        x = "iterations",
        y = "mean time [s]",
        caption = caption,
        title = "Baum-Welch"
    )
ggsave("figure_C4_raw.pdf", height = 5, width = 9)

data_grouped %>% ggplot(aes(iterations, mean / iterations, color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(
        aes(
            ymin = (mean - sd) / iterations,
            ymax = (mean + sd) / iterations
        ),
        size = 0.3,
        alpha = .65
    ) +
    #facet_wrap(. ~ algorithm, scales = "free") +
    labs(
        x = "iterations",
        y = "mean time [s] scaled",
        caption = caption,
        title = "Baum-Welch"
    )
ggsave("figure_C4_scaled.pdf",
       height = 5,
       width = 9)
