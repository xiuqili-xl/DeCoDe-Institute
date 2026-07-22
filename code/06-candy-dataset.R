# Goal ----
# 2017 survey data for likely trick or treaters
# Katie has wrangled the data so that each we know the % of people that 
# rated each candy JOY, MEH, DESPAIR, Not.Answered 


# Load libaries ----
library(tidyverse)
library(here)


# Import data ----
candy_data <- read_csv(here("data_raw", "candy_wrangled_data_forR.csv"),
                       name_repair = "universal") %>%
  rename(candy = column_name)


# Inspect data ----
candy_data

candy_data %>%
  mutate(total = DESPAIR + JOY + MEH + Not.Answered) %>%
  filter(round(total, digits = 1) != 1)
## so the four columns do add up to 100% 


# Graph ranking ----
ggplot(data = candy_data,
       mapping = aes(x = reorder(candy, rev(JOY)), y = JOY)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3))


# Color code by Chocolate ----
ggplot(data = candy_data,
       mapping = aes(x = reorder(candy, rev(JOY)), y = JOY, fill = Chocolate)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3))


# Correlation between Joy and Dispair ----
ggplot(data = candy_data,
       mapping = aes(x = JOY, y = DESPAIR, color = Chocolate, fill = Chocolate)) +
  geom_point(shape = 21, alpha = 0.5) +
  theme_bw()


ggplot(data = candy_data,
       mapping = aes(x = JOY, y = DESPAIR, color = Chocolate, fill = Chocolate)) +
  geom_point(shape = 21, alpha = 0.5) +
  xlim(0, 1) +
  ylim(1, 0) +                  # this is so favored candies are on top
  coord_fixed(ratio = 1) +
  theme_bw()


# Stacked bar graph? ----
candy_data %>%
  pivot_longer(cols = -c("candy", "Chocolate"), names_to = "Feelings", values_to = "Percentage") %>%
  mutate(Feelings = factor(Feelings, levels = c("JOY", "MEH", "DESPAIR", "Not.Answered"))) %>%
  left_join(arrange(candy_data, JOY) %>% mutate(order = row_number()) %>% select(candy, order), by = "candy") %>%
  ggplot(mapping = aes(x = reorder(candy, rev(order)), y = Percentage, fill = fct_rev(Feelings))) +
  geom_col() +
  ggokabeito::scale_fill_okabe_ito(order = c(3, 2, 4, 6)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 3))





  



