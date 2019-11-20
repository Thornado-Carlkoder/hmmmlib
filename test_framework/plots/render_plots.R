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
#data = read_csv("../nov_16.csv")
data = read_csv("../sparse_data.csv")





## Running time tests on varying input size ##
data_input = data %>% filter(test == "inputsize")

data_input_grouped = data_input %>%
    mutate_at(vars(iterations), factor) %>% 
    group_by(observations, algorithm, variant, iterations) %>% 
    summarise(mean = mean(time), sd = sd(time)) 



data_input_grouped %>% filter(iterations == "1" | is.na(iterations)) %>% ggplot(aes(observations, mean, color = variant)) +
    geom_point() + 
    geom_line() + 
    geom_errorbar(aes(ymin=(mean-sd), ymax=(mean+sd) ), size = 0.3, alpha = .65) +
    facet_wrap(.~algorithm, scales = "free") +
    labs(y = "mean time [s]", caption = "error bars: standard deviation of 5 replicates", title = "Running time of increasing input size") + 
    geom_hline(yintercept = 0, alpha = 0)
ggsave("inputsize_raw.pdf", height = 5, width = 9)


data_input_grouped %>% filter(iterations == "1" | is.na(iterations)) %>% ggplot(aes(observations, mean/observations, color = variant)) +
    geom_point() + 
    geom_line() + 
    geom_errorbar(aes(ymin=(mean-sd)/observations, ymax=(mean+sd)/observations ), size = 0.3, alpha = .65) +
    facet_wrap(.~algorithm, scales = "free") +
    labs(y = "normalized mean time [s]", caption = "error bars: standard deviation of 5 replicates", title = "Normalized running time of increasing input size") +
    geom_hline(yintercept = 0, alpha = 0)
ggsave("inputsize_normalized.pdf", height = 5, width = 9)





## Running time of varying state space
data_statespace = data %>% filter(test == "statespace")



data_statespace_grouped = data_statespace %>%
    group_by(observations, algorithm, variant, iterations) %>% 
    summarise(mean = mean(time), sd = sd(time)) 

data_statespace_grouped %>% ggplot(aes(observations, mean, color = variant)) + 
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=(mean-sd), ymax=(mean+sd) ), size = 0.3, alpha = .65) +
    facet_wrap(.~algorithm, scales = "free") + 
    labs(x = "hidden states", y = "mean time [s]", caption = "error bars: standard deviation of 5 replicates", title = "Running time of increasing hidden state space")
ggsave("state_space_raw.pdf", height = 5, width = 9)

data_statespace_grouped %>% ggplot(aes(observations, mean/(observations^2), color = variant)) + 
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=(mean-sd)/observations^2, ymax=(mean+sd)/observations^2 ), size = 0.3, alpha = .65) +
    facet_wrap(.~algorithm, scales = "free") + 
    labs(x = "hidden states", y = "mean time [s]", caption = "error bars: standard deviation of 5 replicates", title = "Running time of increasing hidden state space")
ggsave("state_space_normalized.pdf", height = 5, width = 9)


# Sparseness
data_sparseness = data %>% filter(test == "sparseness") %>% filter(is.na(iterations) | iterations == 1)

data_sparseness_grouped = data_sparseness %>%
    group_by(observations, algorithm, variant, iterations) %>% 
    summarise(mean = mean(time), sd = sd(time)) 

data_sparseness_grouped %>% ggplot(aes(1-observations, mean, color = variant)) + 
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=(mean-sd), ymax=(mean+sd) ), size = 0.3, alpha = .65) +
    facet_wrap(.~algorithm, scales = "free") + 
    labs(x = "hidden states", y = "mean time [s]", caption = "error bars: standard deviation of 5 replicates", title = "Running time of increasing denseness")
ggsave("denseness_raw.pdf", height = 5, width = 9)

data_sparseness_grouped %>% ggplot(aes(observations, mean/(observations^2), color = variant)) + 
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin=(mean-sd)/observations^2, ymax=(mean+sd)/observations^2 ), size = 0.3, alpha = .65) +
    facet_wrap(.~algorithm, scales = "free") + 
    labs(x = "hidden states", y = "mean time [s]", caption = "error bars: standard deviation of 5 replicates", title = "Running time of increasing denseness")
ggsave("denseness_normalized.pdf", height = 5, width = 9)











