# County-level epidemiological analysis of preventive healthcare access and hypertension

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(scales)
})

url <- "https://data.cdc.gov/api/views/cwsq-ngmh/rows.csv?accessType=DOWNLOAD"

raw <- readr::read_csv(url, show_col_types = FALSE)

# Keep 2022 county-level records and pivot indicators into modeling columns
model_df <- raw %>%
  janitor::clean_names() %>%
  filter(
    year == 2022,
    locationdesc == "County"
  ) %>%
  select(
    stateabbr,
    locationname,
    topic,
    question,
    data_value_unit,
    data_value_type,
    data_value
  ) %>%
  mutate(indicator = case_when(
    str_detect(question, regex("Hypertension", ignore_case = TRUE)) ~ "hypertension_prev",
    str_detect(question, regex("Routine checkup|Visited doctor|Preventive", ignore_case = TRUE)) ~ "preventive_care_access",
    str_detect(question, regex("health insurance|uninsured", ignore_case = TRUE)) ~ "insurance_status",
    str_detect(question, regex("obesity", ignore_case = TRUE)) ~ "obesity_prev",
    str_detect(question, regex("smoking|smoker", ignore_case = TRUE)) ~ "smoking_prev",
    str_detect(question, regex("physical inactivity|no leisure-time physical activity", ignore_case = TRUE)) ~ "physical_inactivity_prev",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(indicator)) %>%
  group_by(stateabbr, locationname, indicator) %>%
  summarize(value = mean(data_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = indicator,
    values_from = value
  ) %>%
  drop_na(
    hypertension_prev,
    preventive_care_access,
    insurance_status,
    obesity_prev,
    smoking_prev,
    physical_inactivity_prev
  )

# EDA plots
p1 <- ggplot(model_df, aes(hypertension_prev)) +
  geom_histogram(bins = 30, fill = "#2c7fb8", color = "white") +
  labs(
    title = "Distribution of county-level hypertension prevalence",
    x = "Hypertension prevalence (%)",
    y = "County count"
  ) +
  theme_minimal()

p2 <- ggplot(model_df, aes(preventive_care_access, hypertension_prev)) +
  geom_point(alpha = 0.35, color = "#7fcdbb") +
  geom_smooth(method = "lm", se = TRUE, color = "#253494") +
  labs(
    title = "Preventive care access vs hypertension prevalence",
    x = "Preventive care access (%)",
    y = "Hypertension prevalence (%)"
  ) +
  theme_minimal()

p3 <- ggplot(model_df, aes(insurance_status, hypertension_prev)) +
  geom_point(alpha = 0.35, color = "#f03b20") +
  geom_smooth(method = "lm", se = TRUE, color = "#7f0000") +
  labs(
    title = "Insurance status vs hypertension prevalence",
    x = "Insurance coverage (%)",
    y = "Hypertension prevalence (%)"
  ) +
  theme_minimal()

ggsave("eda_hypertension_distribution.png", p1, width = 8, height = 5, dpi = 300)
ggsave("eda_preventive_vs_hypertension.png", p2, width = 8, height = 5, dpi = 300)
ggsave("eda_insurance_vs_hypertension.png", p3, width = 8, height = 5, dpi = 300)

# Multivariable model with interaction for effect modification
fit <- lm(
  hypertension_prev ~ preventive_care_access * insurance_status +
    obesity_prev + smoking_prev + physical_inactivity_prev,
  data = model_df
)

model_results <- tidy(fit, conf.int = TRUE) %>%
  mutate(across(where(is.numeric), ~ round(.x, 4)))

write_csv(model_results, "hypertension_model_results.csv")

# Interaction visualization
pred_grid <- tidyr::expand_grid(
  preventive_care_access = seq(min(model_df$preventive_care_access), max(model_df$preventive_care_access), length.out = 100),
  insurance_status = quantile(model_df$insurance_status, probs = c(0.25, 0.5, 0.75)),
  obesity_prev = mean(model_df$obesity_prev),
  smoking_prev = mean(model_df$smoking_prev),
  physical_inactivity_prev = mean(model_df$physical_inactivity_prev)
) %>%
  mutate(
    insurance_group = factor(
      insurance_status,
      labels = c("Low insurance (25th pct)", "Median insurance (50th pct)", "High insurance (75th pct)")
    ),
    pred_hypertension = predict(fit, newdata = .)
  )

p_int <- ggplot(pred_grid, aes(preventive_care_access, pred_hypertension, color = insurance_group)) +
  geom_line(linewidth = 1.2) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Interaction effect: preventive care access × insurance status",
    x = "Preventive care access (%)",
    y = "Predicted hypertension prevalence (%)",
    color = "Insurance status"
  ) +
  theme_minimal()

ggsave("interaction_effect_plot.png", p_int, width = 8, height = 5, dpi = 300)

cat("Rows in analytical dataset:", nrow(model_df), "\n")
cat("\nModel summary:\n")
print(summary(fit))
cat("\nTidy coefficient table saved to hypertension_model_results.csv\n")
