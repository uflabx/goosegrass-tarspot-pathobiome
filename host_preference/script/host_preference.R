library(tidyverse)
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(car)
library(ggplot2)
library(patchwork)

select <- dplyr::select

leaf <- read_csv(
  "host_preference/data/leaf_sheath_germination.csv",
  show_col_types = FALSE
)

glass <- read_csv(
  "host_preference/data/glass_germination.csv",
  show_col_types = FALSE
)

leaf <- leaf %>%
  mutate(
    Experiment = factor(Experiment),
    Inoculum_source = factor(
      Inoculum_source,
      levels = c("Corn", "Goosegrass")
    ),
    Leaf_sheath_host = factor(
      Leaf_sheath_host,
      levels = c("Corn", "Goosegrass")
    ),
    Time_hai = factor(
      Time_hai,
      levels = c(24, 48),
      labels = c("24 hai", "48 hai")
    )
  )

glass <- glass %>%
  mutate(
    Experiment = factor(Experiment),
    Inoculum_source = factor(
      Inoculum_source,
      levels = c("Corn", "Goosegrass")
    ),
    Time_hai = factor(Time_hai)
  )

# QA

glimpse(leaf)
glimpse(glass)

dim(leaf)
dim(glass)

summary(leaf)
summary(glass)

anyNA(leaf)
anyNA(glass)

leaf_duplicates <- leaf %>%
  count(
    Experiment,
    Inoculum_source,
    Leaf_sheath_host,
    Rep,
    Time_hai
  ) %>%
  filter(n > 1)

glass_duplicates <- glass %>%
  count(
    Experiment,
    Inoculum_source,
    Rep,
    Time_hai
  ) %>%
  filter(n > 1)

leaf_duplicates
glass_duplicates

stopifnot(
  nrow(leaf_duplicates) == 0,
  nrow(glass_duplicates) == 0,
  all(
    leaf$Total ==
      leaf$Germinated +
      leaf$Non_germinated
  ),
  all(
    glass$Total ==
      glass$Germinated +
      glass$Non_germinated
  ),
  all(leaf$Germinated >= 0),
  all(leaf$Non_germinated >= 0),
  all(glass$Germinated >= 0),
  all(glass$Non_germinated >= 0),
  all(leaf$Germinated <= leaf$Total),
  all(glass$Germinated <= glass$Total),
  all(leaf$Total > 0),
  all(glass$Total > 0)
)

leaf %>%
  count(
    Experiment,
    Leaf_sheath_host,
    Inoculum_source,
    Time_hai
  ) %>%
  print(n = Inf)

glass %>%
  count(
    Experiment,
    Inoculum_source,
    Time_hai
  ) %>%
  print(n = Inf)

leaf_observed_summary <- leaf %>%
  group_by(
    Leaf_sheath_host,
    Inoculum_source,
    Time_hai
  ) %>%
  summarise(
    n = n(),
    Germinated = sum(Germinated),
    Total = sum(Total),
    observed_germination = Germinated / Total,
    .groups = "drop"
  )

leaf_observed_summary

glass_observed_summary <- glass %>%
  group_by(
    Experiment,
    Inoculum_source
  ) %>%
  summarise(
    n = n(),
    Germinated = sum(Germinated),
    Total = sum(Total),
    observed_germination = Germinated / Total,
    .groups = "drop"
  )

glass_observed_summary

# Full leaf-sheath model

model_germ_leaf_time <- glmmTMB(
  cbind(
    Germinated,
    Non_germinated
  ) ~
    Leaf_sheath_host *
    Inoculum_source *
    Time_hai +
    (1 | Experiment),
  family = betabinomial(link = "logit"),
  data = leaf
)

summary(model_germ_leaf_time)

anova_germ_leaf_time <- car::Anova(
  model_germ_leaf_time,
  type = 3
)

anova_germ_leaf_time

set.seed(123)

dharma_leaf_time <- simulateResiduals(
  fittedModel = model_germ_leaf_time,
  n = 1000
)

plot(dharma_leaf_time)

testUniformity(dharma_leaf_time)
testDispersion(dharma_leaf_time)
testOutliers(dharma_leaf_time)

plotResiduals(
  dharma_leaf_time,
  leaf$Leaf_sheath_host
)

plotResiduals(
  dharma_leaf_time,
  leaf$Inoculum_source
)

plotResiduals(
  dharma_leaf_time,
  leaf$Time_hai
)

emm_leaf_time <- emmeans(
  model_germ_leaf_time,
  ~ Leaf_sheath_host |
    Inoculum_source * Time_hai,
  type = "response"
)

emm_leaf_time

contrast_leaf_time <- pairs(
  emm_leaf_time,
  by = c(
    "Inoculum_source",
    "Time_hai"
  ),
  adjust = "sidak"
)

contrast_leaf_time

# Corn-derived inoculum model

model_time_corn_inoculum <- glmmTMB(
  cbind(
    Germinated,
    Non_germinated
  ) ~
    Leaf_sheath_host *
    Time_hai +
    (1 | Experiment),
  family = betabinomial(link = "logit"),
  data = leaf %>%
    filter(Inoculum_source == "Corn")
)

summary(model_time_corn_inoculum)

anova_time_corn <- car::Anova(
  model_time_corn_inoculum,
  type = 3
)

anova_time_corn

set.seed(123)

dharma_time_corn <- simulateResiduals(
  fittedModel = model_time_corn_inoculum,
  n = 1000
)

plot(dharma_time_corn)

testUniformity(dharma_time_corn)
testDispersion(dharma_time_corn)
testOutliers(dharma_time_corn)

plotResiduals(
  dharma_time_corn,
  leaf %>%
    filter(Inoculum_source == "Corn") %>%
    pull(Leaf_sheath_host)
)

plotResiduals(
  dharma_time_corn,
  leaf %>%
    filter(Inoculum_source == "Corn") %>%
    pull(Time_hai)
)

emm_time_corn_inoculum <- emmeans(
  model_time_corn_inoculum,
  ~ Leaf_sheath_host * Time_hai,
  type = "response"
)

emm_time_corn_inoculum

# Goosegrass-derived inoculum model

model_time_goosegrass_inoculum <- glmmTMB(
  cbind(
    Germinated,
    Non_germinated
  ) ~
    Leaf_sheath_host *
    Time_hai +
    (1 | Experiment),
  family = betabinomial(link = "logit"),
  data = leaf %>%
    filter(Inoculum_source == "Goosegrass")
)

summary(model_time_goosegrass_inoculum)

anova_time_goosegrass <- car::Anova(
  model_time_goosegrass_inoculum,
  type = 3
)

anova_time_goosegrass

set.seed(123)

dharma_time_goosegrass <- simulateResiduals(
  fittedModel = model_time_goosegrass_inoculum,
  n = 1000
)

plot(dharma_time_goosegrass)

testUniformity(dharma_time_goosegrass)
testDispersion(dharma_time_goosegrass)
testOutliers(dharma_time_goosegrass)

plotResiduals(
  dharma_time_goosegrass,
  leaf %>%
    filter(Inoculum_source == "Goosegrass") %>%
    pull(Leaf_sheath_host)
)

plotResiduals(
  dharma_time_goosegrass,
  leaf %>%
    filter(Inoculum_source == "Goosegrass") %>%
    pull(Time_hai)
)

emm_time_goosegrass_inoculum <- emmeans(
  model_time_goosegrass_inoculum,
  ~ Leaf_sheath_host * Time_hai,
  type = "response"
)

emm_time_goosegrass_inoculum

# Glass reference model

model_glass_reference <- glmmTMB(
  cbind(
    Germinated,
    Non_germinated
  ) ~
    Inoculum_source +
    (1 | Experiment),
  family = betabinomial(link = "logit"),
  data = glass
)

summary(model_glass_reference)

anova_glass_reference <- car::Anova(
  model_glass_reference,
  type = 3
)

anova_glass_reference

set.seed(123)

dharma_glass_reference <- simulateResiduals(
  fittedModel = model_glass_reference,
  n = 1000
)

plot(dharma_glass_reference)

testUniformity(dharma_glass_reference)
testDispersion(dharma_glass_reference)
testOutliers(dharma_glass_reference)

plotResiduals(
  dharma_glass_reference,
  glass$Inoculum_source
)

emm_glass_reference <- emmeans(
  model_glass_reference,
  ~ Inoculum_source,
  type = "response"
)

emm_glass_reference

# Outputs

fig_6_leaf_data <- bind_rows(
  as.data.frame(
    emm_time_corn_inoculum
  ) %>%
    mutate(
      Inoculum_source = "Corn"
    ),
  as.data.frame(
    emm_time_goosegrass_inoculum
  ) %>%
    mutate(
      Inoculum_source = "Goosegrass"
    )
) %>%
  mutate(
    Inoculum_label = factor(
      Inoculum_source,
      levels = c("Corn", "Goosegrass"),
      labels = c(
        "Corn-derived inoculum",
        "Goosegrass-derived inoculum"
      )
    ),
    Leaf_sheath_host = factor(
      Leaf_sheath_host,
      levels = c("Corn", "Goosegrass")
    ),
    Display_group = as.character(
      Leaf_sheath_host
    ),
    Time_numeric = case_when(
      Time_hai == "24 hai" ~ 1,
      Time_hai == "48 hai" ~ 2
    ),
    x_plot = case_when(
      Leaf_sheath_host == "Corn" ~
        Time_numeric - 0.10,
      Leaf_sheath_host == "Goosegrass" ~
        Time_numeric + 0.10
    )
  )

fig_6_glass_data <- as.data.frame(
  emm_glass_reference
) %>%
  mutate(
    Inoculum_label = factor(
      Inoculum_source,
      levels = c("Corn", "Goosegrass"),
      labels = c(
        "Corn-derived inoculum",
        "Goosegrass-derived inoculum"
      )
    ),
    Display_group = "Glass",
    x_plot = 0.8
  )

host_time_pvalues <- tibble(
  Inoculum_source = c(
    "Corn",
    "Goosegrass"
  ),
  p_value = c(
    anova_time_corn[
      "Leaf_sheath_host:Time_hai",
      "Pr(>Chisq)"
    ],
    anova_time_goosegrass[
      "Leaf_sheath_host:Time_hai",
      "Pr(>Chisq)"
    ]
  )
) %>%
  mutate(
    Inoculum_label = factor(
      Inoculum_source,
      levels = c(
        "Corn",
        "Goosegrass"
      ),
      labels = c(
        "Corn-derived inoculum",
        "Goosegrass-derived inoculum"
      )
    ),
    p_label = paste0(
      "'Host'~'\u00d7'~'time'~~italic(p)==",
      sprintf(
        "%.3f",
        p_value
      )
    ),
    x = 1.5,
    y = 1.05
  )

anova_germ_leaf_time_df <- as.data.frame(
  anova_germ_leaf_time
) %>%
  rownames_to_column(
    "Term"
  )

anova_time_corn_df <- as.data.frame(
  anova_time_corn
) %>%
  rownames_to_column(
    "Term"
  )

anova_time_goosegrass_df <- as.data.frame(
  anova_time_goosegrass
) %>%
  rownames_to_column(
    "Term"
  )

anova_glass_reference_df <- as.data.frame(
  anova_glass_reference
) %>%
  rownames_to_column(
    "Term"
  )

# Figures

host_colours_host <- c(
  "Corn" = "#648FFF",
  "Goosegrass" = "#FE6100"
)

condition_colours <- c(
  "Corn" = "#648FFF",
  "Goosegrass" = "#FE6100",
  "Glass" = "grey50"
)


fig_6_host_p <- as.data.frame(
  contrast_leaf_time
) %>%
  filter(
    Time_hai == "24 hai"
  ) %>%
  mutate(
    Inoculum_label = factor(
      Inoculum_source,
      levels = c(
        "Corn",
        "Goosegrass"
      ),
      labels = c(
        "Corn-derived inoculum",
        "Goosegrass-derived inoculum"
      )
    ),
    p_label = case_when(
      p.value < 0.0001 ~ "italic(P) < 0.0001",
      TRUE ~ paste0(
        "italic(p) == ",
        formatC(
          p.value,
          format = "f",
          digits = 3
        )
      )
    ),
    x1 = 0.90,
    x2 = 1.10,
    bracket_y = case_when(
      Inoculum_source == "Corn" ~ 0.94,
      Inoculum_source == "Goosegrass" ~ 0.72
    ),
    bracket_tip_y = bracket_y - 0.025,
    p_y = bracket_y + 0.035
  )

fig_6 <- ggplot(
  fig_6_leaf_data,
  aes(
    x = x_plot,
    y = prob,
    color = Display_group,
    group = Leaf_sheath_host
  )
) +
  geom_line(
    linewidth = 0.8
  ) +
  geom_point(
    size = 3.4
  ) +
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    width = 0.05,
    linewidth = 0.6
  ) +
  geom_segment(
    data = fig_6_glass_data,
    aes(
      x = x_plot,
      xend = x_plot,
      y = asymp.LCL,
      yend = asymp.UCL,
      color = Display_group
    ),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_segment(
    data = fig_6_glass_data,
    aes(
      x = x_plot - 0.025,
      xend = x_plot + 0.025,
      y = asymp.LCL,
      yend = asymp.LCL,
      color = Display_group
    ),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_segment(
    data = fig_6_glass_data,
    aes(
      x = x_plot - 0.025,
      xend = x_plot + 0.025,
      y = asymp.UCL,
      yend = asymp.UCL,
      color = Display_group
    ),
    inherit.aes = FALSE,
    linewidth = 0.6
  ) +
  geom_point(
    data = fig_6_glass_data,
    aes(
      x = x_plot,
      y = prob,
      color = Display_group
    ),
    inherit.aes = FALSE,
    size = 3.4
  ) +
  geom_text(
    data = host_time_pvalues,
    aes(
      x = x,
      y = y,
      label = p_label
    ),
    inherit.aes = FALSE,
    parse = TRUE,
    color = "black",
    size = 3.3
  ) +
  geom_segment(
    data = fig_6_host_p,
    aes(
      x = x1,
      xend = x2,
      y = bracket_y,
      yend = bracket_y
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  geom_segment(
    data = fig_6_host_p,
    aes(
      x = x1,
      xend = x1,
      y = bracket_y,
      yend = bracket_tip_y
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  geom_segment(
    data = fig_6_host_p,
    aes(
      x = x2,
      xend = x2,
      y = bracket_y,
      yend = bracket_y
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  geom_segment(
    data = fig_6_host_p,
    aes(
      x = x2,
      xend = x2,
      y = bracket_y,
      yend = bracket_tip_y
    ),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.5
  ) +
  geom_text(
    data = fig_6_host_p,
    aes(
      x = 1,
      y = p_y,
      label = p_label
    ),
    inherit.aes = FALSE,
    parse = TRUE,
    color = "black",
    size = 3.4
  ) +
  facet_wrap(
    ~ Inoculum_label
  ) +
  scale_color_manual(
    values = condition_colours,
    breaks = c(
      "Glass",
      "Corn",
      "Goosegrass"
    ),
    labels = c(
      "Glass slide",
      "Corn leaf sheath",
      "Goosegrass leaf sheath"
    )
  ) +
  scale_x_continuous(
    breaks = c(
      1,
      2
    ),
    labels = c(
      "24 h",
      "48 h"
    ),
    limits = c(
      0.65,
      2.3
    )
  ) +
  scale_y_continuous(
    limits = c(
      0,
      1.05
    ),
    breaks = seq(
      0,
      1,
      0.2
    ),
    labels = scales::percent_format(
      accuracy = 1
    )
  ) +
  labs(
    x = "Time after inoculation",
    y = "Ascospore germination",
    color = "Germination environment"
  ) +
  theme_classic(
    base_family = "Arial",
    base_size = 12
  ) +
  theme(
    legend.position = "right",
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold"
    )
  )

fig_6

ggsave(
  "Figure6.png",
  plot = fig_6,
  width = 7.2,
  height = 6,
  units = "in"
)

sessionInfo()
