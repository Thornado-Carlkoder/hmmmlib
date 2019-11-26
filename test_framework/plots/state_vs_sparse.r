library(tidyverse)

setwd("~/bioinformatics/hmm/git_tc_hmmmlib/test_framework/plots")

<<<<<<< HEAD
data = read_csv("../statevssparse.csv", col_names = F)
=======
data = read_csv("../fw_bw_statesparse.csv", col_names = F)
>>>>>>> master
names(data) = c('test',
                'observations',
                'time',
                'algorithm',
                'variant',
                'statespace')

## Varying state space
# data %>% View 


data_grouped = data %>%
    group_by(observations, algorithm, variant, statespace) %>%
    summarise(mean = mean(time), sd = sd(time))

data_grouped %>% ggplot(aes(1-observations, mean, color = variant)) +
    geom_point() +
    geom_line() +
    geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), size = 0.3, alpha = .65) +
<<<<<<< HEAD
    facet_grid(statespace ~ algorithm, scales = "free") +
=======
    facet_grid(statespace ~ algorithm, scales = "free_y") +
>>>>>>> master
    labs(
        x = "denseness",
        y = "mean time [s]",
        caption = "error bars: standard deviation of 5 replicates",
        title = "Denseness for different state spaces"
    )
ggsave("denseness_vs_statespace_problem.pdf", height = 5, width = 9)
