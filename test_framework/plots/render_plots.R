library(tidyverse)
library(ggpubr)
library(svglite)

setwd("~/bioinformatics/hmm/git_tc_hmmmlib/test_framework/plots")



## Inputsize ##

data = read_csv("../newdata/input.csv", col_names = F)
names(data) = c('test',
                'observations',
                'time',
                'algorithm',
                'variant',
                'iterations')
data_input = data %>% filter(test == "inputsize")

data_input_grouped = data_input %>%
    mutate_at(vars(iterations), factor) %>%
    group_by(observations, algorithm, variant, iterations) %>%
    summarise(mean = mean(time), sd = sd(time))



data_input_grouped %>% filter(iterations == "1" |
                                  is.na(iterations)) %>% ggplot(aes(observations, mean, color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), size = 0.3, alpha = .65) +
    facet_wrap(. ~ algorithm, scales = "free") +
    labs(y = "mean time [s]", caption = "error bars: standard deviation of 5 replicates", title = "Running time of increasing input size") +
    geom_hline(yintercept = 0, alpha = 0)
ggsave("inputsize_raw.pdf", height = 6, width = 10)


data_input_grouped %>% filter(iterations == "1" |
                                  is.na(iterations)) %>% ggplot(aes(observations, mean / observations, color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(
        ymin = (mean - sd) / observations,
        ymax = (mean + sd) / observations
    ),
    size = 0.3,
    alpha = .65) +
    facet_wrap(. ~ algorithm, scales = "free") +
    labs(y = "normalized mean time [s]", caption = "error bars: standard deviation of 5 replicates", title = "Normalized running time of increasing input size") +
    geom_hline(yintercept = 0, alpha = 0)
ggsave("inputsize_normalized.pdf",
       height = 6,
       width = 10)





## State space ##

data = read_csv("../newdata/statespace.csv", col_names = F)
names(data) = c('test',
                'observations',
                'time',
                'algorithm',
                'variant',
                'iterations')
data_statespace = data %>% filter(test == "statespace")


data_statespace_grouped = data_statespace %>%
    group_by(observations, algorithm, variant, iterations) %>%
    summarise(mean = mean(time), sd = sd(time))

data_statespace_grouped %>% ggplot(aes(observations, mean, color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), size = 0.3, alpha = .65) +
    facet_wrap(. ~ algorithm, scales = "free") +
    labs(
        x = "hidden states",
        y = "mean time [s]",
        caption = "error bars: standard deviation of 5 replicates\ninput size: 120000, alphabet size: 4",
        title = "Running time of increasing hidden state space"
    )
ggsave("state_space_raw.pdf", height = 5, width = 9)

data_statespace_grouped %>% ggplot(aes(observations, mean / (observations ^
                                                                 2), color = variant)) +
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
        y = "mean time [s]",
        caption = "error bars: standard deviation of 5 replicates\ninput size: 120000, alphabet size: 4",
        title = "Running time of increasing hidden state space"
    )
ggsave("state_space_normalized.pdf",
       height = 5,
       width = 9)



## Inputsize ##

data = read_csv("../newdata/sparse.csv", col_names = F)
names(data) = c('test',
                'observations',
                'time',
                'algorithm',
                'variant',
                'iterations')

data_sparseness = data %>% filter(test == "sparseness_1500") %>% filter(is.na(iterations) |
                                                                       iterations == 1)

data_sparseness_grouped = data_sparseness %>%
    group_by(observations, algorithm, variant, iterations) %>%
    summarise(mean = mean(time), sd = sd(time))

data_sparseness_grouped %>% ggplot(aes(1 - observations, mean, color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), size = 0.3, alpha = .65) +
    facet_wrap(. ~ algorithm, scales = "free") +
    labs(
        x = "denseness",
        y = "mean time [s]",
        caption = "error bars: standard deviation of 5 replicates\ninput size: 90000, alphabet size: 4, hidden states: 16",
        title = "Running time of increasing denseness"
    )
ggsave("denseness_raw.pdf", height = 5, width = 9)







## Alphabet ##
data = read_csv("../newdata/alphabet.csv", col_names = F)
names(data) = c('test',
                'observations',
                'time',
                'algorithm',
                'variant',
                'iterations')
data_alphabet = data %>% filter(test == "alphabetsize") %>% filter(is.na(iterations) |
                                                                       iterations == 1)

data_alphabet_grouped = data_alphabet %>%
    group_by(observations, algorithm, variant, iterations) %>%
    summarise(mean = mean(time), sd = sd(time))

data_alphabet_grouped %>% ggplot(aes(observations, mean, color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), size = 0.3, alpha = .65) +
    facet_wrap(. ~ algorithm, scales = "free") +
    labs(
        x = "alphabet size",
        y = "mean time [s]",
        caption = "error bars: standard deviation of 5 replicates\ninput size = 90000, hidden states: 8",
        title = "Running time of increasing alphabet size"
    )
ggsave("alphabet_raw.pdf", height = 5, width = 9)






## Statesparse

data = read_csv("../newdata/statesparse.csv", col_names = F)
names(data) = c('test',
                'observations',
                'time',
                'algorithm',
                'variant',
                'statespace')
data$algorithm = data %>% pull(algorithm) %>% recode(forward = "Algorithm: Forward",
                                                     backward_time = "Algorithm: Backward")


data_grouped = data %>%
    group_by(observations, algorithm, variant, statespace) %>%
    summarise(mean = mean(time), sd = sd(time))

data_grouped %>% ggplot(aes(1-observations, mean, color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), width = 0.05, alpha = .65) +
    facet_grid(statespace ~ algorithm, scales = "free") +
    
    labs(
        x = "denseness",
        y = "mean time [s]",
        caption = "error bars: standard deviation of 5 replicates \ninputsize: 90000 characters, alphabet size: 4",
        title = "Denseness for different state spaces"
    )
#geom_hline(yintercept = 0, alpha = 0)
ggsave("denseness_vs_statespace.pdf", height = 7, width = 10)

