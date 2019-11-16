library(tidyverse)
library(ggpubr)
library(svglite)

setwd("~/bioinformatics/hmm/git_tc_hmmmlib/test_framework/plots")
#data <- read_csv("../time2.csv",)
#data <- read_csv("../short.csv")
#data <- read_csv("../medium.csv")
#data <- read_csv("../oct-31.csv")
#data <- read_csv("../oct-31-2.csv")

#data = read_csv("../../hmmmlib/build/withblas.csv")
data = read_csv("../nov_16.csv")

# running time tests on varying input size
data_input = data %>% filter(test == "inputsize")

data_input_grouped = data_input %>%
    mutate_at(vars(iterations), factor) %>% 
    group_by(observations, algorithm, variant, iterations) %>% 
    summarise(mean = mean(time), sd = sd(time)) 





p_mean = data_input_grouped %>% filter(iterations == "1" | is.na(iterations)) %>% ggplot(aes(observations, mean, color = variant)) +
    geom_point() + 
    geom_line() + 
    geom_errorbar(aes(ymin=(mean-sd), ymax=(mean+sd) ), size = 0.3, alpha = .65) +
    facet_wrap(.~algorithm, scales = "free") +
    labs(y = "mean time [s]", caption = "error bars: standard deviation", title = "Running time of varying input size") + 
    geom_hline(yintercept = 0, alpha = 0)
p_mean

p_mean_normalized = data_input_grouped %>% filter(iterations == "1" | is.na(iterations)) %>% ggplot(aes(observations, mean/observations, color = variant)) +
    geom_point() + 
    geom_line() + 
    geom_errorbar(aes(ymin=(mean-sd)/observations, ymax=(mean+sd)/observations ), width = 10000, size = 0.3, alpha = .65) +
    facet_wrap(.~algorithm, scales = "free") +
    labs(y = "normalized mean time [s]", caption = "error bars: standard deviation", title = "Running time of varying input size") +
    geom_hline(yintercept = 0, alpha = 0)
p_mean_normalized



ggsave("inputsize_raw.pdf", height = 5, width = 9, plot = p_mean)
ggsave("inputsize_normalized.pdf", height = 5, width = 9, plot = p_mean_normalized)


















