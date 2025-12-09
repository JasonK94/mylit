# Script to add SET and SET_DETAIL columns to config
library(dplyr)
library(readr)

config_path <- "config/config_complete2.csv"
config <- read_csv(config_path, show_col_types = FALSE)

# Define mapping logic
config <- config %>%
    mutate(
        SET = case_when(
            gem_name %in% c("GEM1", "GEM2", "GEM3", "GEM4") ~ "SET1",
            gem_name %in% c("GEM5", "GEM6", "GEM7", "GEM8") ~ "SET2",
            gem_name %in% c("GEM9", "GEM10", "GEM11", "GEM12") ~ "SET3",
            TRUE ~ "Other"
        ),
        SET_DETAIL = case_when(
            gem_name %in% c("GEM1", "GEM2", "GEM3", "GEM4") ~ "SET1",
            gem_name %in% c("GEM5", "GEM6", "GEM7", "GEM8") ~ "SET2",
            gem_name == "GEM9" ~ "SET3_1",
            gem_name == "GEM10" ~ "SET3_2",
            gem_name %in% c("GEM11", "GEM12") ~ "SET3_3",
            TRUE ~ "Other"
        )
    )

# Write back
write_csv(config, config_path)
print(paste("Updated", config_path, "with SET and SET_DETAIL columns."))
