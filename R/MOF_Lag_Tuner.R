# *************************************************************************************************************************
# Metro Market Office Flow - Lag Tuner
# *************************************************************************************************************************
# This script tunes the lag lengths of one fundamental's forecast drivers by
# coordinate descent, scored on retro-vintage backtests - the same
# specifications as the agent's tune_lags tool, in R, so a person can run the
# whole process (screen, tune, forecast, assess) from RStudio without the
# agent. Candidates come from the CCF screener's peak lags plus the incumbent
# config values; the winner must also beat the initial lags on a long
# out-of-sample window from a 2021 Q1 origin (the validation gate). The
# agent's tune_lags calls this script through Rscript and reads the results
# back from tune_history.csv.
# Stephen Gardner, smg263@georgetown.edu
# One regex dialect throughout: the portable core every engine shares
# (python, R with or without perl, grep) - explicit ranges like [0-9] and
# [A-Za-z0-9_] instead of the \d and \w escapes, [ \t] for whitespace.
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
# trailing arguments from the command line - how Rscript (and so the agent)
# names the fundamental to tune; empty when the script runs in RStudio
args <- commandArgs(trailingOnly = TRUE)
# search constants - the same values the agent's tuner uses
max_evals <- 40L
lag_min <- 0L
lag_max <- 12L
search_steps <- 4L
gate_origin_date <- "2021-01-01"
gate_steps <- 20L
# read the project config - paths, lags and the vintage setup all come from here
if (!file.exists("config.yml"))
  stop("config.yml not found - open the .Rproj so the working directory is the project root.")
config <- yaml::read_yaml("config.yml")
# *************************************************************************************************************************
# FORMAT LAGS
# *************************************************************************************************************************
# one lag set as a printable string - loop the drivers, collect one
# "name = value" piece per driver, then join the pieces into a single line
format_lags <- function(lags) {
  parts <- c()
  # walk the drivers in the set, one printable piece each
  for (x in names(lags)) parts <- c(parts, glue("{x} = {lags[[x]]}"))
  glue_collapse(parts, sep = ", ")
}
# a WAPE value as printable text, n/a when the percentage is undefined
wape_text <- function(wape) {
  if (is.na(wape)) "n/a" else sprintf("%.2f%%", wape)
}
# *************************************************************************************************************************
# SET LAGS
# *************************************************************************************************************************
# swap one fundamental's lag values in the config text - locate the
# fundamental's block by its header line, then one str_replace per driver swaps
# its number. everything else in the file (comments included) is untouched
set_lags <- function(lines, fundamental, lags) {
  # the block runs from the fundamental's header to the next value-less header line
  start <- str_which(lines, glue("^[ \t]*{fundamental}:"))[1]
  headers <- str_which(lines, "^[ \t]*[A-Za-z0-9_.]+:[ \t]*$")
  end <- min(headers[headers > start], length(lines) + 1)
  block <- lines[start:(end - 1)]
  # one swap per driver - keep everything before the number, replace the number
  for (nm in names(lags)) {
    block <- str_replace(block, glue("^([ \t]*{nm}:[ \t]*)[0-9]+"), glue("\\1{as.integer(lags[[nm]])}"))
  }
  lines[start:(end - 1)] <- block
  lines
}
# the two vintage swaps for a retro run - keep the indent and the key, swap the quoted value
set_vintage <- function(lines, retro_date) {
  lines <- str_replace(lines, "^([ \t]*desc:[ \t]*)\".*\"", "\\1\"retro\"")
  lines <- str_replace(lines, "^([ \t]*retro_date:[ \t]*)\".*\"", glue("\\1\"{retro_date}\""))
  lines
}
# *************************************************************************************************************************
# RETRO LABEL
# *************************************************************************************************************************
# the vintage label a retro_date produces, in the forecast script's own format (2025-10-01 -> 2025Q4)
retro_label <- function(date_str) {
  d <- ymd(date_str)
  glue("{year(d)}Q{quarter(d)}")
}
# *************************************************************************************************************************
# VINTAGE LABELS
# *************************************************************************************************************************
# every vintage label on disk for the finalcast files, newest first - sorted
# by quarter, not write time, so a retro run never reorders the comparison
vintage_labels <- function() {
  # start with the finalcast files and pipe to the end - the label is the piece after the last _ and before the .csv
  labels <- list.files(config$paths$forecast_dir, pattern = "^finalcast\\.offc_fund_.*\\.csv$") %>%
    str_extract("[^_]+(?=\\.csv$)")
  # sort by the quarter each label names - yq turns a label back into its quarter's first day
  labels[order(yq(labels), decreasing = TRUE)]
}
# *************************************************************************************************************************
# SCORE RUN
# *************************************************************************************************************************
# score one retro run: its forecasts against the newest vintage's actuals,
# for one fundamental only, up to `steps` quarters ahead - RMSE plus WAPE
# (the volume weighted absolute percentage error, the human readable one).
# actuals live in the newest report file - the finalcast files are
# forecast-only by construction
score_run <- function(retro_lab, fundamental, steps) {
  empty <- list(rmse = NA_real_, wape = NA_real_)
  labels <- vintage_labels()
  if (length(labels) < 2) return(empty)
  actuals_path <- glue("{config$paths$report_dir}/offc_fund_forecast_{labels[1]}.csv")
  fcast_path <- glue("{config$paths$forecast_dir}/finalcast.offc_fund_{retro_lab}.csv")
  if (!file.exists(actuals_path) || !file.exists(fcast_path)) return(empty)
  # the vintage quarter as an ordinal - a forecast's horizon is its own quarter minus this
  origin_ord <- year(yq(retro_lab)) * 4 + quarter(yq(retro_lab))
  # start with the actuals df and pipe to the end df - this fundamental's actual rows from the newest report file
  actuals.df <- read_csv(actuals_path, show_col_types = FALSE) %>%
    filter(cast_desc == "actual", offc_fundamental_desc == fundamental) %>%
    select(msa_desc, class_desc, submkt_desc, offc_fundamental_desc, date, actual_val = offc_fundamental_val)
  # start with the forecast df and pipe to the end df - this fundamental's forecast rows inside the horizon
  fcast.df <- read_csv(fcast_path, show_col_types = FALSE) %>%
    filter(cast_desc == "forecast", offc_fundamental_desc == fundamental) %>%
    mutate(q_ord = year(date) * 4 + quarter(date)) %>%
    filter(q_ord - origin_ord <= steps) %>%
    select(msa_desc, class_desc, submkt_desc, offc_fundamental_desc, date, forecast_val = offc_fundamental_val)
  # join the forecasts to the actuals on the segment keys and the quarter, then summarise the errors
  scored.df <- fcast.df %>%
    inner_join(actuals.df, by = c("msa_desc", "class_desc", "submkt_desc", "offc_fundamental_desc", "date")) %>%
    mutate(err = actual_val - forecast_val)
  if (nrow(scored.df) == 0) return(empty)
  denom <- sum(abs(scored.df$actual_val))
  list(rmse = sqrt(mean(scored.df$err^2)),
       wape = if (denom > 0) 100 * sum(abs(scored.df$err)) / denom else NA_real_)
}
# *************************************************************************************************************************
# CANDIDATES FOR
# *************************************************************************************************************************
# the candidate lags for one driver: {peak - 1, peak, peak + 1} around each
# screener peak, plus the incumbent - a hand edited config lag is a manual
# override that gets scored head to head, never silently dropped.
# no peaks for the driver -> the +-2 fallback around the incumbent
candidates_for <- function(peaks.df, fundamental, regressor, current) {
  # this driver's peak lags from the screener table, if it has any
  peaks <- peaks.df %>%
    filter(offc_fundamental == fundamental, variable == regressor) %>%
    pull(lag) %>%
    unique()
  cands <- current
  if (length(peaks) == 0) {
    # fallback - walk two either side of the incumbent, inside the allowed range
    for (d in c(-2, -1, 1, 2)) if (current + d >= lag_min && current + d <= lag_max) cands <- c(cands, current + d)
  } else {
    # walk the neighbourhood of every peak, inside the allowed range
    for (p in peaks) for (d in c(-1, 0, 1)) if (p + d >= lag_min && p + d <= lag_max) cands <- c(cands, p + d)
  }
  sort(unique(cands))
}
# *************************************************************************************************************************
# TUNE FUNDAMENTAL
# *************************************************************************************************************************
# tune one fundamental's lags: back up the config, score the baseline, run
# the coordinate descent, put the winner through the 2021 gate, adopt what
# survives, and append every evaluation to the tune history
tune_fundamental <- function(fundamental) {
  if (!fundamental %in% names(config$lags)) {
    stop(glue("Unknown fundamental {fundamental}. Choices: {glue_collapse(names(config$lags), sep = ', ')}"))
  }
  regressors <- names(config$lags[[fundamental]])
  current <- unlist(config$lags[[fundamental]])
  standing_retro_date <- as.character(config$vintage$retro_date)
  # snapshot the config text and back the file up before touching anything
  original_lines <- readLines("config.yml", warn = FALSE)
  backup_name <- glue("config.yml.bak-{format(Sys.time(), '%Y%m%d-%H%M%S')}")
  file.copy("config.yml", backup_name)
  completed <- FALSE
  # if the tune ends abnormally, put the exact pre-tune config back
  on.exit(if (!completed) {
    writeLines(original_lines, "config.yml")
    cat("  [tune] tune ended abnormally - config.yml restored to its exact pre-tune state\n")
  }, add = TRUE)
  # the screener's peaks, if it has been run - the candidate source. an empty
  # table makes the search a fallback search, and that is said out loud
  peaks_path <- glue("{config$paths$output_dir}/corr_offc_fcast_vars_peaks.csv")
  peaks.df <- if (file.exists(peaks_path)) read_csv(peaks_path, show_col_types = FALSE) else tibble(offc_fundamental = character(), variable = character(), peak_rank = integer(), lag = integer(), correlation = numeric())
  if (nrow(peaks.df) == 0) {
    spec_note <- "WARNING: no CCF screener spec loaded - candidates fall back to +-2 around the current lags"
    cat(glue("  [tune] {spec_note}"), "
", sep = "")
  } else {
    spec_note <- "candidates: CCF-guided (screener peaks + incumbent)"
  }
  history_list <- list()
  log_lines <- spec_note
  eval_count <- 0
  # one evaluation - write the lag set into the config (retro mode, at the
  # given origin date), run the forecast script, score the run, log the row
  run_eval <- function(lags, phase, steps, retro_date) {
    eval_count <<- eval_count + 1
    # start from the pre-tune text every time - swap this evaluation's lags
    # and the retro vintage in, write it, and the next evaluation starts clean
    eval_lines <- set_lags(original_lines, fundamental, lags)
    eval_lines <- set_vintage(eval_lines, retro_date)
    writeLines(eval_lines, "config.yml")
    # this run's vintage label and its finalcast file - note when the run
    # starts, so a failed run is never scored off a stale file afterwards
    lab <- retro_label(retro_date)
    fcast_path <- glue("{config$paths$forecast_dir}/finalcast.offc_fund_{lab}.csv")
    run_started <- Sys.time()
    status <- system2("Rscript", "MOF_Forecast_Fundamentals.R", stdout = FALSE, stderr = FALSE)
    fresh <- file.exists(fcast_path) && file.mtime(fcast_path) > run_started
    scored <- if (status == 0 && fresh) score_run(lab, fundamental, steps) else list(rmse = NA_real_, wape = NA_real_)
    # keep every evaluation in a list - the rows are bound into one df at the end
    history_list[[length(history_list) + 1]] <<- tibble(
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
      fundamental = fundamental,
      eval = eval_count,
      phase = phase,
      lags = format_lags(lags),
      rmse = scored$rmse,
      wape_pct = if (is.na(scored$wape)) NA_real_ else round(scored$wape, 4),
      winner = FALSE)
    score_txt <- if (is.na(scored$rmse)) "FAIL" else format(scored$rmse, scientific = FALSE, trim = TRUE)
    cat(glue("  [tune] {phase} eval {eval_count}: {format_lags(lags)} -> RMSE {score_txt} (WAPE {wape_text(scored$wape)})"), "
", sep = "")
    scored
  }
  # baseline - the incumbent set, scored before any candidate is tried
  best <- current
  baseline <- run_eval(current, "search", search_steps, standing_retro_date)
  if (is.na(baseline$rmse)) stop("Baseline run failed or had nothing to score - check the R pipeline.")
  best_score <- baseline$rmse
  best_wape <- baseline$wape
  baseline_score <- best_score
  baseline_wape <- best_wape
  # the evaluation numbers of the baseline and of the running best - how the
  # adopted set is matched back to its history row at the end, no string keys
  baseline_eval <- eval_count
  best_eval <- baseline_eval
  log_lines <- c(log_lines, glue("baseline RMSE={format(baseline_score, digits = 4)} (WAPE {wape_text(baseline_wape)}) @ {format_lags(current)}"))
  # coordinate descent - one pass over the drivers; a candidate that scores
  # better than the running best is adopted on the spot, so later drivers
  # are tried against the improved set
  for (reg in regressors) {
    if (eval_count >= max_evals) break
    # walk this driver's candidates against the running best
    for (cand in candidates_for(peaks.df, fundamental, reg, best[[reg]])) {
      if (cand == best[[reg]] || eval_count >= max_evals) next
      trial <- best
      trial[[reg]] <- cand
      scored <- run_eval(trial, "search", search_steps, standing_retro_date)
      score_txt <- if (is.na(scored$rmse)) "FAIL" else format(scored$rmse, scientific = FALSE, trim = TRUE)
      log_lines <- c(log_lines, glue("  {reg}={cand}: RMSE={score_txt} (WAPE {wape_text(scored$wape)})"))
      if (!is.na(scored$rmse) && scored$rmse < best_score) {
        best <- trial
        best_score <- scored$rmse
        best_wape <- scored$wape
        best_eval <- eval_count
      }
    }
    log_lines <- c(log_lines, glue("  -> {reg}: kept {best[[reg]]} (RMSE {format(best_score, digits = 4)})"))
  }
  # the 2021 validation gate - the winner must also beat the initial lags on
  # the long window from the 2021 Q1 origin; both sets are run at that origin
  adopted <- best
  if (identical(as.numeric(best), as.numeric(current))) {
    gate_note <- "gate skipped - the search kept the initial lags, nothing to validate"
  } else {
    gate_winner <- run_eval(best, "gate_2021", gate_steps, gate_origin_date)
    gate_initial <- run_eval(current, "gate_2021", gate_steps, gate_origin_date)
    if (is.na(gate_winner$rmse) || is.na(gate_initial$rmse)) {
      gate_note <- "gate could not be scored - winner adopted on the search score alone"
    } else if (gate_winner$rmse < gate_initial$rmse) {
      gate_note <- glue("gate PASSED - winner {format(gate_winner$rmse, digits = 4)} (WAPE {wape_text(gate_winner$wape)}) beats initial {format(gate_initial$rmse, digits = 4)} (WAPE {wape_text(gate_initial$wape)}) on the 2021 window")
    } else {
      adopted <- current
      gate_note <- glue("gate REJECTED the search winner - winner {format(gate_winner$rmse, digits = 4)} (WAPE {wape_text(gate_winner$wape)}) vs initial {format(gate_initial$rmse, digits = 4)} (WAPE {wape_text(gate_initial$wape)}) on the 2021 window; initial lags kept")
    }
  }
  log_lines <- c(log_lines, glue("  GATE: {gate_note}"))
  # adopt - write the surviving lag set into the config; the original text
  # with only the lags changed, so the vintage setup ends exactly as it began
  writeLines(set_lags(original_lines, fundamental, adopted), "config.yml")
  # the adopted set's evaluation - the search winner's, or the baseline's
  # when the gate rejected the winner and the initial lags were kept
  adopted_eval <- if (identical(as.numeric(adopted), as.numeric(best))) best_eval else baseline_eval
  # bind the history list into one df - flag the adopted set's evaluation -
  # and append it to the tune history
  history.df <- bind_rows(history_list)
  history.df$winner[history.df$eval == adopted_eval] <- TRUE
  hist_path <- glue("{config$paths$output_dir}/tune_history.csv")
  write_csv(history.df, hist_path, append = file.exists(hist_path), na = "")
  search_runs <- sum(history.df$phase == "search")
  headline <- glue("TUNE DONE for {fundamental} ({search_runs} search runs). Baseline RMSE={format(baseline_score, digits = 4)} (WAPE {wape_text(baseline_wape)}) -> search winner: {format_lags(best)} with RMSE={format(best_score, digits = 4)} (WAPE {wape_text(best_wape)}).")
  if (!identical(as.numeric(adopted), as.numeric(best))) headline <- glue("{headline} Adopted: {format_lags(adopted)} (gate rejection - see log).")
  completed <- TRUE
  cat(glue("{headline}"), "
", sep = "")
  cat(glue("Backup of the pre-tune config: {backup_name}"), "
", sep = "")
  cat(glue("Tune history appended to tune_history.csv ({nrow(history.df)} evaluations)."), "
", sep = "")
  cat(glue_collapse(log_lines, sep = "\n"), "\n", sep = "")
  invisible(history.df)
}
# *************************************************************************************************************************
# RUN
# *************************************************************************************************************************
# with a trailing argument, tune just that fundamental; with none (RStudio),
# tune every modelled fundamental in turn - each tune picks the config up
# exactly where the last one left it
fundamentals_to_run <- if (length(args) >= 1) args else c("total_inventory", "total_vacancy", "total_absorption", "wgtavg_rent")
# loop the fundamentals, one full tune each
for (f in fundamentals_to_run) tune_fundamental(f)
