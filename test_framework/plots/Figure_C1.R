library(tidyverse)
library(ggpubr)
library(svglite)

## Alphabet ##
data = read_csv("../running_times/running_time_C1.csv", col_names = T)
#names(data) = c('test', 'observations', 'time', 'algorithm', 'variant', 'iterations')
#data_alphabet = data %>% filter(test == "alphabetsize") %>% filter(is.na(iterations) |
                                                                       #iterations == 1)

data_alphabet_grouped = data %>%
    group_by(observations, algorithm, variant) %>%
    summarise(mean = mean(time), sd = sd(time))

data_alphabet_grouped %>% ggplot(aes(observations, mean, color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), size = 0.3, alpha = .65) +
    facet_wrap(. ~ algorithm, scales = "free") +
    labs(
        x = "alphabet size",
        y = "mean time [s]",
        caption = "error bars: standard deviation of 3 replicates\ninput size = 100000, hidden states: 8"
        #title = "Running time of increasing alphabet size"
    )
ggsave("pdf/figure_C1.pdf", height = 5, width = 9)




