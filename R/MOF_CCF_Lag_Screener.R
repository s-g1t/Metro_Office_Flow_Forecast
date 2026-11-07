# *************************************************************************************************************************
# Metro Market Office Flow Forecast - CCF Lag Screener
# *************************************************************************************************************************

# This script screens cross-correlations (on first differences) between each office fundamental and the forecast drivers 
# over a pre-defined lag grid. The strongest correlations on the grid define the lag search space for the tuning agent.
#
# Stephen Gardner, smg263@georgetown.edu

# *************************************************************************************************************************
# CONFIGURATION
# *************************************************************************************************************************

# clear the environment prior to starting
rm(list = ls())

# set the formatting for decimals
options(scipen = 999)

# load libraries
library(tidyverse)
library(lubridate)
library(glue)
library(yaml)

# read the project config - all data paths and driver lists come from here
if (!file.exists("config.yml"))
  stop("config.yml not found - open the .Rproj so the working directory is the project root.")
config <- yaml::read_yaml("config.yml")

# *************************************************************************************************************************
# DATA IMPORT
# *************************************************************************************************************************

offc.fundamentals <- read_csv(config$data$offc.fundamentals)
offc.metro_fcast_vars <- read_csv(config$data$offc.metro_fcast_vars)

# not every fundamental needs econometric tuning
target_fundamentals <- c("total_inventory", "total_vacancy", "total_absorption", "wgtavg_rent")

# create a lag list
lag_list <- c(0:6, 8, 12, 16, 20)
corr_list <- list()

# *************************************************************************************************************************
# AGGREGATE TO TOTAL MARKET
# *************************************************************************************************************************
# tune once per fundamental on the aggregate series, since fable is running the same model per submarket/class desc

offc.funds_agg <- offc.fundamentals %>%
                  group_by(date) %>%
                  summarise(total_inventory  = sum(total_inventory, na.rm = TRUE),
                            total_absorption = sum(total_absorption, na.rm = TRUE),
                            total_vacancy    = sum(total_vacancy, na.rm = TRUE),
                            wgtavg_rent      = sum(wgtavg_rent * total_inventory, na.rm = TRUE) / sum(total_inventory, na.rm = TRUE),
                           .groups = "drop") %>%
                  arrange(date) %>%
                  pivot_longer(cols = all_of(target_fundamentals),
                  names_to = "offc_fundamental", values_to = "fundamental_val")

# split into one data frame per fundamental and store them in a list, so the correlation loop can work through them one at a time
fund_split <- split(offc.funds_agg, offc.funds_agg$offc_fundamental)
fund_split <- keep(fund_split, ~ nrow(.x) > 0)

# *************************************************************************************************************************
# RUN CORRELATION LOOP
# *************************************************************************************************************************
# ccf on first differences

for (i in names(fund_split)) 
  {var_corr <- list()
   current_fund <- i
# fundamental drivers are in the config file
   fund_vars <- names(config$lags[[current_fund]])
  
   for (x in fund_vars) 
     {if (!x %in% names(offc.metro_fcast_vars)) next
# join the fundamental to the driver on date and drop incomplete rows
      comparison_data <- fund_split[[i]] %>%
                         select(date, fundamental_val) %>%
                         inner_join(offc.metro_fcast_vars %>% select(date, all_of(x)), by = "date") %>%
                         drop_na(fundamental_val, all_of(x)) %>%
                         arrange(date)
    
      if (nrow(comparison_data) < 30) next
      ccf_matrix <- ccf(diff(comparison_data$fundamental_val),
                        diff(comparison_data[[x]]),
                        lag.max = 20,
                        plot = FALSE,
                        na.action = na.omit)
    
      var_corr[[x]] <- tibble(variable = x,
                              lag = as.vector(ccf_matrix$lag),
                              correlation = as.vector(ccf_matrix$acf)) %>%
                       mutate(lag = glue("Lag {lag}")) %>%
                       filter(lag %in% glue("Lag {lag_list}")) %>%
                       pivot_wider(names_from = lag, values_from = correlation)
      } # end x in fund_vars
   
   if (length(var_corr) > 0) 
     {corr_list[[i]] <- bind_rows(var_corr) %>% mutate(offc_fundamental = current_fund)}
} # end i in fund_split

# *************************************************************************************************************************
# BIND AND EXPORT
# *************************************************************************************************************************

# bind the per-fundamental results into one data frame
corr.offc_fcast_vars <- bind_rows(corr_list) %>%
                        select(offc_fundamental, variable, everything())

write_csv(corr.offc_fcast_vars, glue("{config$paths$output_dir}/corr_offc_fcast_vars.csv"))
print(glue("Correlation table written to {config$paths$output_dir}/corr_offc_fcast_vars.csv"))

# pull the top 3 |correlation| lags per fundamental x driver - these peaks are
# the tuner's search space, and this table is the same one a person would
# use to hand-pick lags, so the process on disk matches the manual process
lag_peaks <- corr.offc_fcast_vars %>%
             pivot_longer(cols = starts_with("Lag"), 
                          names_to = "lag", values_to = "correlation") %>%
             mutate(lag = as.integer(str_remove(lag, "Lag ")),
                    abs_corr = abs(correlation)) %>%
             group_by(offc_fundamental, variable) %>%
             arrange(desc(abs_corr), .by_group = TRUE) %>%
             slice_head(n = 3) %>%
             mutate(peak_rank = row_number()) %>%
             ungroup() %>%
             select(offc_fundamental, variable, peak_rank, lag, correlation)

write_csv(lag_peaks, glue("{config$paths$output_dir}/corr_offc_fcast_vars_peaks.csv"))
print(glue("Peak lags written to {config$paths$output_dir}/corr_offc_fcast_vars_peaks.csv"))

# store the peaks in a list, one entry per fundamental x driver, so they can be pulled straight into other code blocks
peak_list <- list()
for (i in seq_len(nrow(lag_peaks))) 
  {
  peak_list[[lag_peaks$offc_fundamental[i]]][[lag_peaks$variable[i]]] <- c(peak_list[[lag_peaks$offc_fundamental[i]]][[lag_peaks$variable[i]]], lag_peaks$lag[i])
}
