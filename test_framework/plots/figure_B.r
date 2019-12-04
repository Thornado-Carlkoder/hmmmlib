library(tidyverse)
library(ggpubr)
library(svglite)


# This is for plotting figure B
# For the algorithms Forward and Backward, we want to show that the sparse implementations are faster than the dense implementations.




setwd("~/bioinformatics/hmm/git_tc_hmmmlib/test_framework/plots")

n_edges = function(n, d) {
    # n: number of hidden states. 
    # d: density
    n+(n*n-n)*d
}



## Statesparse

#data = read_csv("../newdata/statesparse.csv", col_names = F)
#data = read_csv("../bigvs.csv", col_names = F)
data = read_csv("../data/B.csv")
names(data) = c('test',
                'observations',
                'time',
                'algorithm',
                'variant',
                'statespace')
data$algorithm = data %>% pull(algorithm) %>% recode(forward = "Forward",
                                                     backward_time = "Backward")

caption = "error bars: standard deviation of 1 replicates \ninputsize: 90000 characters, alphabet size: 4"
data_grouped = data %>%
    #filter(variant != "CSR") %>% 
    mutate(observations = 1-observations) %>% # turn sparseness into density
    group_by(observations, algorithm, variant, statespace) %>%
    summarise(mean = mean(time), sd = sd(time))

# raw
data_grouped %>% filter(variant %in% c("RSB", "BLAS")) %>% 
    ggplot(aes(observations, mean, mean, color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), width = 0.05, alpha = .65) +
    facet_grid(statespace ~ algorithm, scales = "free") +
    
    labs(
        x = "density",
        y = "mean time [s]",
        caption = caption,
        title = "Running time for density versus state space",
        subtitle = "Higher state space prefers sparse implementation"
    )
#geom_hline(yintercept = 0, alpha = 0)
ggsave("figure_B.pdf", height = 7, width = 10)



# scaled by the number of edges
data_grouped %>% filter(observations >= 0.1) %>%  ggplot(aes(observations, mean/n_edges(statespace, observations), mean, color = variant)) +
    geom_point() +
    geom_line() +
    #geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), width = 0.05, alpha = .65) +
    facet_grid(statespace ~ algorithm, scales = "free") +
    labs(
        x = "density",
        y = "mean time [s] scaled to number of edges",
        caption = caption,
        title = "Running time for density versus state space"
    )
#geom_hline(yintercept = 0, alpha = 0)
ggsave("figure_B_scaled.pdf", height = 7, width = 10)



