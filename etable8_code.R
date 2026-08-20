# eTable 8. Incremental prognostic performance of dentition
# Full reproduction: Uno C, time-dependent AUC, Brier score,
# continuous NRI, Rubin-pooled Wald tests, and B=200 bootstrap
# optimism correction for delta Uno C.

library(readxl)
library(dplyr)
library(survival)
library(tibble)

DATA_PATH <- "KFACS_analytic_dataset_MICE_v2.xlsx"
SHEET_NAME <- "analytic_dataset_MI"

OUTPUT_DIR <- file.path(dirname(DATA_PATH), "KFACS_incremental_prognostic")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

B <- 200
MASTER_SEED <- 20260819
RUN_BOOTSTRAP <- TRUE

TIME_POINTS <- c(5, 8)
NRI_TIME <- 5
G_MIN <- 1e-6

required_columns <- c(
  "_imputation_id",
  "cohort_start_year",
  "followup_years",
  "death_event",
  "age",
  "sex",
  "edu",
  "income_level",
  "smoking_history_ge100",
  "ever_alcohol_use",
  "systemic_disease_count",
  "bmi",
  "mna_scr_gr",
  "frailty_3cat",
  "mmse_kc_score",
  "systemic_reserve_adverse_count",
  "functional_dentition",
  "natural_teeth"
)

dat <- read_excel(DATA_PATH, sheet = SHEET_NAME)

missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns) > 0) {
  stop(
    "Required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

dat <- dat %>%
  filter(!is.na(cohort_start_year)) %>%
  mutate(
    death_event = as.integer(death_event),
    sex = factor(sex),
    income_level = factor(income_level),
    smoking_history_ge100 = factor(smoking_history_ge100),
    ever_alcohol_use = factor(ever_alcohol_use),
    mna_scr_gr = factor(mna_scr_gr),
    frailty_3cat = factor(frailty_3cat),
    fd_lt10 = as.integer(functional_dentition < 10)
  )

mi_list <- split(dat, dat$`_imputation_id`)
mi_ids <- names(mi_list)

if (length(mi_list) != 20) {
  warning("Expected 20 imputed datasets; found ", length(mi_list), ".")
}

cat("Number of imputations =", length(mi_list), "\n")
cat("N in imputation 1 =", nrow(mi_list[[1]]), "\n")
cat(
  "Deaths in imputation 1 =",
  sum(mi_list[[1]]$death_event == 1, na.rm = TRUE),
  "\n"
)

# -------------------------------------------------------------------------
# Model definitions
# -------------------------------------------------------------------------

f_base <- Surv(followup_years, death_event) ~
  age + sex + edu + income_level +
  smoking_history_ge100 + ever_alcohol_use +
  systemic_disease_count + bmi + mna_scr_gr +
  frailty_3cat + mmse_kc_score +
  systemic_reserve_adverse_count

extended_models <- list(
  "Continuous functional dentition" =
    update(f_base, . ~ . + functional_dentition),
  "Continuous natural teeth" =
    update(f_base, . ~ . + natural_teeth),
  "Functional dentition <10" =
    update(f_base, . ~ . + fd_lt10)
)

added_terms <- c(
  "Continuous functional dentition" = "functional_dentition",
  "Continuous natural teeth" = "natural_teeth",
  "Functional dentition <10" = "fd_lt10"
)

analysis_order <- names(extended_models)

# -------------------------------------------------------------------------
# Core functions
# -------------------------------------------------------------------------

fit_cox <- function(formula, data) {
  coxph(
    formula,
    data = data,
    ties = "efron",
    x = TRUE,
    y = TRUE,
    model = TRUE
  )
}

get_uno_c <- function(fit, newdata = NULL) {
  cc <- if (is.null(newdata)) {
    survival::concordance(fit, timewt = "n/G2")
  } else {
    survival::concordance(
      fit,
      newdata = newdata,
      timewt = "n/G2"
    )
  }
  as.numeric(cc$concordance)
}

get_baseline_hazard_at <- function(fit, t) {
  bh <- survival::basehaz(fit, centered = FALSE)
  eligible <- bh$time <= t

  if (!any(eligible)) {
    return(0)
  }

  max(bh$hazard[eligible], na.rm = TRUE)
}

predict_event_risk <- function(fit, newdata, t) {
  h0 <- get_baseline_hazard_at(fit, t)

  lp <- as.numeric(
    predict(
      fit,
      newdata = newdata,
      type = "lp",
      reference = "zero"
    )
  )

  risk <- 1 - exp(-h0 * exp(lp))
  pmin(pmax(risk, 0), 1)
}

make_censoring_survival <- function(time, event) {
  fit_g <- survfit(Surv(time, 1L - event) ~ 1)

  km_time <- fit_g$time
  km_surv <- fit_g$surv

  g_at <- function(x) {
    vapply(
      x,
      function(z) {
        idx <- which(km_time <= z)
        g <- if (length(idx) == 0) 1 else km_surv[max(idx)]
        max(g, G_MIN)
      },
      numeric(1)
    )
  }

  g_left <- function(x) {
    vapply(
      x,
      function(z) {
        idx <- which(km_time < z)
        g <- if (length(idx) == 0) 1 else km_surv[max(idx)]
        max(g, G_MIN)
      },
      numeric(1)
    )
  }

  list(at = g_at, left = g_left)
}

ipcw_auc <- function(time, event, risk, t, gfun = NULL) {
  if (is.null(gfun)) {
    gfun <- make_censoring_survival(time, event)
  }

  cases <- which(event == 1 & time <= t)
  controls <- which(time > t)

  if (length(cases) == 0 || length(controls) == 0) {
    return(NA_real_)
  }

  w_case <- 1 / gfun$left(time[cases])
  w_control <- rep(1 / gfun$at(t), length(controls))

  comparison <- outer(
    risk[cases],
    risk[controls],
    FUN = "-"
  )

  concordant <- (comparison > 0) + 0.5 * (comparison == 0)
  pair_weights <- outer(w_case, w_control, FUN = "*")

  denom <- sum(pair_weights)
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }

  sum(concordant * pair_weights) / denom
}

ipcw_brier <- function(time, event, risk, t, gfun = NULL) {
  if (is.null(gfun)) {
    gfun <- make_censoring_survival(time, event)
  }

  n <- length(time)
  contribution <- numeric(n)

  cases <- event == 1 & time <= t
  controls <- time > t

  if (any(cases)) {
    contribution[cases] <-
      ((1 - risk[cases])^2) / gfun$left(time[cases])
  }

  if (any(controls)) {
    contribution[controls] <-
      (risk[controls]^2) / gfun$at(t)
  }

  mean(contribution)
}

ipcw_continuous_nri <- function(
  time,
  event,
  risk_base,
  risk_extended,
  t,
  gfun = NULL
) {
  if (is.null(gfun)) {
    gfun <- make_censoring_survival(time, event)
  }

  cases <- which(event == 1 & time <= t)
  controls <- which(time > t)

  if (length(cases) == 0 || length(controls) == 0) {
    return(
      c(
        NRI = NA_real_,
        Event_NRI = NA_real_,
        Nonevent_NRI = NA_real_
      )
    )
  }

  delta <- risk_extended - risk_base

  w_case <- 1 / gfun$left(time[cases])
  w_control <- rep(1 / gfun$at(t), length(controls))

  event_direction <- sign(delta[cases])
  nonevent_direction <- -sign(delta[controls])

  event_nri <- sum(w_case * event_direction) / sum(w_case)
  nonevent_nri <- sum(w_control * nonevent_direction) / sum(w_control)

  c(
    NRI = event_nri + nonevent_nri,
    Event_NRI = event_nri,
    Nonevent_NRI = nonevent_nri
  )
}

extract_added_term <- function(fit, term) {
  beta <- coef(fit)
  vc <- vcov(fit)

  if (!term %in% names(beta)) {
    stop("Added model term not found: ", term)
  }

  list(
    beta = as.numeric(beta[term]),
    variance = as.numeric(vc[term, term])
  )
}

rubin_pool_scalar <- function(q, u) {
  keep <- is.finite(q) & is.finite(u) & u >= 0
  q <- q[keep]
  u <- u[keep]

  m <- length(q)
  if (m == 0) {
    return(
      tibble(
        beta = NA_real_,
        se = NA_real_,
        df = NA_real_,
        p_value = NA_real_
      )
    )
  }

  q_bar <- mean(q)
  u_bar <- mean(u)
  b <- if (m > 1) var(q) else 0
  total_var <- u_bar + (1 + 1 / m) * b
  se <- sqrt(total_var)

  if (!is.finite(se) || se <= 0) {
    p <- NA_real_
    df <- NA_real_
  } else if (b <= .Machine$double.eps) {
    df <- Inf
    p <- 2 * pnorm(-abs(q_bar / se))
  } else {
    r <- ((1 + 1 / m) * b) / u_bar

    if (!is.finite(r) || r <= 0) {
      df <- Inf
      p <- 2 * pnorm(-abs(q_bar / se))
    } else {
      df <- (m - 1) * (1 + 1 / r)^2
      p <- 2 * pt(-abs(q_bar / se), df = df)
    }
  }

  tibble(
    beta = q_bar,
    se = se,
    df = df,
    p_value = p
  )
}

format_p <- function(p) {
  if (!is.finite(p)) return(NA_character_)
  if (p < 0.001) return("<.001")
  sub("^0", "", sprintf("%.3f", p))
}

# -------------------------------------------------------------------------
# Apparent performance in each imputed dataset
# -------------------------------------------------------------------------

performance_by_imp <- list()
wald_by_imp <- list()

for (i in seq_along(mi_list)) {
  d <- mi_list[[i]]

  time <- as.numeric(d$followup_years)
  event <- as.integer(d$death_event)

  fit_base <- fit_cox(f_base, d)
  base_uno <- get_uno_c(fit_base)

  gfun <- make_censoring_survival(time, event)

  base_risk_5 <- predict_event_risk(fit_base, d, 5)
  base_risk_8 <- predict_event_risk(fit_base, d, 8)

  base_auc_5 <- ipcw_auc(
    time, event, base_risk_5, 5, gfun
  )
  base_auc_8 <- ipcw_auc(
    time, event, base_risk_8, 8, gfun
  )

  base_brier_5 <- ipcw_brier(
    time, event, base_risk_5, 5, gfun
  )
  base_brier_8 <- ipcw_brier(
    time, event, base_risk_8, 8, gfun
  )

  for (analysis in analysis_order) {
    fit_ext <- fit_cox(extended_models[[analysis]], d)
    ext_uno <- get_uno_c(fit_ext)

    ext_risk_5 <- predict_event_risk(fit_ext, d, 5)
    ext_risk_8 <- predict_event_risk(fit_ext, d, 8)

    ext_auc_5 <- ipcw_auc(
      time, event, ext_risk_5, 5, gfun
    )
    ext_auc_8 <- ipcw_auc(
      time, event, ext_risk_8, 8, gfun
    )

    ext_brier_5 <- ipcw_brier(
      time, event, ext_risk_5, 5, gfun
    )
    ext_brier_8 <- ipcw_brier(
      time, event, ext_risk_8, 8, gfun
    )

    nri_5 <- ipcw_continuous_nri(
      time = time,
      event = event,
      risk_base = base_risk_5,
      risk_extended = ext_risk_5,
      t = NRI_TIME,
      gfun = gfun
    )

    performance_by_imp[[length(performance_by_imp) + 1]] <-
      tibble(
        Imputation = as.integer(mi_ids[i]),
        Analysis = analysis,
        Reference_UnoC = base_uno,
        Extended_UnoC = ext_uno,
        Delta_UnoC = ext_uno - base_uno,
        Reference_AUC_5y = base_auc_5,
        Extended_AUC_5y = ext_auc_5,
        Delta_AUC_5y = ext_auc_5 - base_auc_5,
        Reference_AUC_8y = base_auc_8,
        Extended_AUC_8y = ext_auc_8,
        Delta_AUC_8y = ext_auc_8 - base_auc_8,
        Reference_Brier_5y = base_brier_5,
        Extended_Brier_5y = ext_brier_5,
        Delta_Brier_5y = ext_brier_5 - base_brier_5,
        Reference_Brier_8y = base_brier_8,
        Extended_Brier_8y = ext_brier_8,
        Delta_Brier_8y = ext_brier_8 - base_brier_8,
        Continuous_NRI_5y = unname(nri_5["NRI"]),
        Event_NRI_5y = unname(nri_5["Event_NRI"]),
        Nonevent_NRI_5y = unname(nri_5["Nonevent_NRI"])
      )

    term_result <- extract_added_term(
      fit_ext,
      added_terms[[analysis]]
    )

    wald_by_imp[[length(wald_by_imp) + 1]] <-
      tibble(
        Imputation = as.integer(mi_ids[i]),
        Analysis = analysis,
        beta = term_result$beta,
        variance = term_result$variance
      )
  }

  cat(
    "Apparent performance completed for imputation",
    i,
    "of",
    length(mi_list),
    "\n"
  )
}

performance_by_imp <- bind_rows(performance_by_imp)
wald_by_imp <- bind_rows(wald_by_imp)

performance_summary <- performance_by_imp %>%
  group_by(Analysis) %>%
  summarise(
    Reference_UnoC = mean(Reference_UnoC, na.rm = TRUE),
    Extended_UnoC = mean(Extended_UnoC, na.rm = TRUE),
    Delta_UnoC = mean(Delta_UnoC, na.rm = TRUE),

    Reference_AUC_5y = mean(Reference_AUC_5y, na.rm = TRUE),
    Extended_AUC_5y = mean(Extended_AUC_5y, na.rm = TRUE),
    Delta_AUC_5y = mean(Delta_AUC_5y, na.rm = TRUE),

    Reference_AUC_8y = mean(Reference_AUC_8y, na.rm = TRUE),
    Extended_AUC_8y = mean(Extended_AUC_8y, na.rm = TRUE),
    Delta_AUC_8y = mean(Delta_AUC_8y, na.rm = TRUE),

    Reference_Brier_5y = mean(Reference_Brier_5y, na.rm = TRUE),
    Extended_Brier_5y = mean(Extended_Brier_5y, na.rm = TRUE),
    Delta_Brier_5y = mean(Delta_Brier_5y, na.rm = TRUE),

    Reference_Brier_8y = mean(Reference_Brier_8y, na.rm = TRUE),
    Extended_Brier_8y = mean(Extended_Brier_8y, na.rm = TRUE),
    Delta_Brier_8y = mean(Delta_Brier_8y, na.rm = TRUE),

    Continuous_NRI_5y = mean(Continuous_NRI_5y, na.rm = TRUE),
    Event_NRI_5y = mean(Event_NRI_5y, na.rm = TRUE),
    Nonevent_NRI_5y = mean(Nonevent_NRI_5y, na.rm = TRUE),
    .groups = "drop"
  )

wald_pooled <- wald_by_imp %>%
  group_by(Analysis) %>%
  group_modify(
    ~ rubin_pool_scalar(.x$beta, .x$variance)
  ) %>%
  ungroup()

# -------------------------------------------------------------------------
# B=200 bootstrap internal validation for delta Uno C
# -------------------------------------------------------------------------

if (RUN_BOOTSTRAP) {
  base_boot_all <- list()
  extended_boot_all <- list()

  for (i in seq_along(mi_list)) {
    d <- mi_list[[i]]
    n <- nrow(d)

    set.seed(MASTER_SEED + i)

    cat(
      "\nBOOTSTRAP: IMPUTATION",
      i,
      "OF",
      length(mi_list),
      "\n"
    )

    for (b in seq_len(B)) {
      if (b == 1 || b %% 20 == 0) {
        cat("Bootstrap", b, "of", B, "\n")
      }

      idx <- sample(
        seq_len(n),
        size = n,
        replace = TRUE
      )

      boot <- d[idx, , drop = FALSE]

      base_result <- tryCatch(
        {
          fit_base_boot <- fit_cox(f_base, boot)

          c_base_boot <- get_uno_c(fit_base_boot)
          c_base_test <- get_uno_c(
            fit_base_boot,
            newdata = d
          )

          if (
            !is.finite(c_base_boot) ||
            !is.finite(c_base_test)
          ) {
            stop("Non-finite reference Uno C")
          }

          tibble(
            Imputation = as.integer(mi_ids[i]),
            Bootstrap = b,
            C_Base_Boot = c_base_boot,
            C_Base_Test = c_base_test,
            Optimism_Base = c_base_boot - c_base_test
          )
        },
        error = function(e) NULL
      )

      if (is.null(base_result)) {
        next
      }

      base_boot_all[[length(base_boot_all) + 1]] <-
        base_result

      for (analysis in analysis_order) {
        ext_result <- tryCatch(
          {
            fit_ext_boot <- fit_cox(
              extended_models[[analysis]],
              boot
            )

            c_ext_boot <- get_uno_c(fit_ext_boot)
            c_ext_test <- get_uno_c(
              fit_ext_boot,
              newdata = d
            )

            if (
              !is.finite(c_ext_boot) ||
              !is.finite(c_ext_test)
            ) {
              stop("Non-finite extended Uno C")
            }

            delta_boot <-
              c_ext_boot - base_result$C_Base_Boot

            delta_test <-
              c_ext_test - base_result$C_Base_Test

            tibble(
              Imputation = as.integer(mi_ids[i]),
              Bootstrap = b,
              Analysis = analysis,
              C_Base_Boot = base_result$C_Base_Boot,
              C_Base_Test = base_result$C_Base_Test,
              C_Extended_Boot = c_ext_boot,
              C_Extended_Test = c_ext_test,
              Optimism_Extended =
                c_ext_boot - c_ext_test,
              Delta_Boot = delta_boot,
              Delta_Test = delta_test,
              Optimism_Delta =
                delta_boot - delta_test
            )
          },
          error = function(e) NULL
        )

        if (!is.null(ext_result)) {
          extended_boot_all[[
            length(extended_boot_all) + 1
          ]] <- ext_result
        }
      }
    }

    saveRDS(
      list(
        base = base_boot_all,
        extended = extended_boot_all
      ),
      file.path(
        OUTPUT_DIR,
        paste0(
          "etable8_B200_checkpoint_after_imputation_",
          i,
          ".rds"
        )
      )
    )
  }

  base_boot_results <- bind_rows(base_boot_all)
  extended_boot_results <- bind_rows(extended_boot_all)

  base_optimism_by_imp <- base_boot_results %>%
    group_by(Imputation) %>%
    summarise(
      Mean_Optimism_Base =
        mean(Optimism_Base, na.rm = TRUE),
      Successful_Base_Bootstraps = n(),
      .groups = "drop"
    )

  extended_optimism_by_imp <- extended_boot_results %>%
    group_by(Imputation, Analysis) %>%
    summarise(
      Mean_Optimism_Extended =
        mean(Optimism_Extended, na.rm = TRUE),
      Mean_Optimism_Delta =
        mean(Optimism_Delta, na.rm = TRUE),
      Successful_Bootstraps = n(),
      .groups = "drop"
    )

  validation_by_imp <- performance_by_imp %>%
    select(
      Imputation,
      Analysis,
      Reference_UnoC,
      Extended_UnoC,
      Delta_UnoC
    ) %>%
    left_join(
      base_optimism_by_imp,
      by = "Imputation"
    ) %>%
    left_join(
      extended_optimism_by_imp,
      by = c("Imputation", "Analysis")
    ) %>%
    mutate(
      Corrected_Base_UnoC =
        Reference_UnoC - Mean_Optimism_Base,
      Corrected_Extended_UnoC =
        Extended_UnoC - Mean_Optimism_Extended,
      Corrected_Delta_UnoC =
        Delta_UnoC - Mean_Optimism_Delta
    )

  bootstrap_summary <- validation_by_imp %>%
    group_by(Analysis) %>%
    summarise(
      Corrected_Base_UnoC =
        mean(Corrected_Base_UnoC, na.rm = TRUE),
      Corrected_Extended_UnoC =
        mean(Corrected_Extended_UnoC, na.rm = TRUE),
      Mean_Optimism_Delta =
        mean(Mean_Optimism_Delta, na.rm = TRUE),
      Corrected_Delta_UnoC =
        mean(Corrected_Delta_UnoC, na.rm = TRUE),
      Mean_Successful_Bootstraps =
        mean(Successful_Bootstraps, na.rm = TRUE),
      .groups = "drop"
    )
} else {
  base_boot_results <- tibble()
  extended_boot_results <- tibble()
  validation_by_imp <- tibble()

  bootstrap_summary <- tibble(
    Analysis = analysis_order,
    Corrected_Base_UnoC = NA_real_,
    Corrected_Extended_UnoC = NA_real_,
    Mean_Optimism_Delta = NA_real_,
    Corrected_Delta_UnoC = NA_real_,
    Mean_Successful_Bootstraps = NA_real_
  )
}

# -------------------------------------------------------------------------
# Final eTable 8
# -------------------------------------------------------------------------

etable8_numeric <- performance_summary %>%
  left_join(
    wald_pooled %>%
      select(
        Analysis,
        beta,
        se,
        df,
        p_value
      ),
    by = "Analysis"
  ) %>%
  left_join(
    bootstrap_summary,
    by = "Analysis"
  ) %>%
  mutate(
    row_order = match(Analysis, analysis_order)
  ) %>%
  arrange(row_order) %>%
  select(-row_order)

etable8_publication <- etable8_numeric %>%
  transmute(
    `Added oral measure` = Analysis,
    `Reference Uno C` = Reference_UnoC,
    `Extended Uno C` = Extended_UnoC,
    `Delta Uno C` = Delta_UnoC,
    `5-y AUC` = Extended_AUC_5y,
    `Delta AUC (5 y)` = Delta_AUC_5y,
    `8-y AUC` = Extended_AUC_8y,
    `Delta AUC (8 y)` = Delta_AUC_8y,
    `5-y Brier score` = Extended_Brier_5y,
    `8-y Brier score` = Extended_Brier_8y,
    `5-y continuous NRI` = Continuous_NRI_5y,
    `Event NRI` = Event_NRI_5y,
    `Nonevent NRI` = Nonevent_NRI_5y,
    `P value` = p_value,
    `Optimism-corrected Delta Uno C` =
      Corrected_Delta_UnoC
  )

etable8_display <- etable8_publication %>%
  mutate(
    `Reference Uno C` =
      sprintf("%.4f", `Reference Uno C`),
    `Extended Uno C` =
      sprintf("%.4f", `Extended Uno C`),
    `Delta Uno C` =
      sprintf("%+.4f", `Delta Uno C`),
    `5-y AUC` =
      sprintf("%.3f", `5-y AUC`),
    `Delta AUC (5 y)` =
      sprintf("%+.4f", `Delta AUC (5 y)`),
    `8-y AUC` =
      sprintf("%.3f", `8-y AUC`),
    `Delta AUC (8 y)` =
      sprintf("%+.4f", `Delta AUC (8 y)`),
    `5-y Brier score` =
      sprintf("%.4f", `5-y Brier score`),
    `8-y Brier score` =
      sprintf("%.4f", `8-y Brier score`),
    `5-y continuous NRI` =
      sprintf("%+.4f", `5-y continuous NRI`),
    `Event NRI` =
      sprintf("%+.4f", `Event NRI`),
    `Nonevent NRI` =
      sprintf("%+.4f", `Nonevent NRI`),
    `P value` =
      vapply(`P value`, format_p, character(1)),
    `Optimism-corrected Delta Uno C` =
      sprintf(
        "%+.4f",
        `Optimism-corrected Delta Uno C`
      )
  )

# -------------------------------------------------------------------------
# Output
# -------------------------------------------------------------------------

write.csv(
  performance_by_imp,
  file.path(
    OUTPUT_DIR,
    "etable8_performance_by_imputation.csv"
  ),
  row.names = FALSE
)

write.csv(
  performance_summary,
  file.path(
    OUTPUT_DIR,
    "etable8_performance_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  wald_by_imp,
  file.path(
    OUTPUT_DIR,
    "etable8_wald_by_imputation.csv"
  ),
  row.names = FALSE
)

write.csv(
  wald_pooled,
  file.path(
    OUTPUT_DIR,
    "etable8_wald_pooled.csv"
  ),
  row.names = FALSE
)

if (RUN_BOOTSTRAP) {
  write.csv(
    base_boot_results,
    file.path(
      OUTPUT_DIR,
      "etable8_B200_base_bootstrap_all.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    extended_boot_results,
    file.path(
      OUTPUT_DIR,
      "etable8_B200_extended_bootstrap_all.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    validation_by_imp,
    file.path(
      OUTPUT_DIR,
      "etable8_B200_validation_by_imputation.csv"
    ),
    row.names = FALSE
  )

  write.csv(
    bootstrap_summary,
    file.path(
      OUTPUT_DIR,
      "etable8_B200_optimism_corrected_summary.csv"
    ),
    row.names = FALSE
  )
}

write.csv(
  etable8_numeric,
  file.path(
    OUTPUT_DIR,
    "etable8_full_numeric.csv"
  ),
  row.names = FALSE
)

write.csv(
  etable8_publication,
  file.path(
    OUTPUT_DIR,
    "etable8_publication_numeric.csv"
  ),
  row.names = FALSE
)

write.csv(
  etable8_display,
  file.path(
    OUTPUT_DIR,
    "etable8_publication_display.csv"
  ),
  row.names = FALSE
)

if (requireNamespace("openxlsx", quietly = TRUE)) {
  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "eTable 8")
  openxlsx::writeData(
    wb,
    "eTable 8",
    etable8_display
  )

  openxlsx::addWorksheet(wb, "Numeric")
  openxlsx::writeData(
    wb,
    "Numeric",
    etable8_numeric
  )

  openxlsx::addWorksheet(wb, "By imputation")
  openxlsx::writeData(
    wb,
    "By imputation",
    performance_by_imp
  )

  openxlsx::addWorksheet(wb, "Wald pooling")
  openxlsx::writeData(
    wb,
    "Wald pooling",
    wald_pooled
  )

  if (RUN_BOOTSTRAP) {
    openxlsx::addWorksheet(wb, "Bootstrap summary")
    openxlsx::writeData(
      wb,
      "Bootstrap summary",
      bootstrap_summary
    )
  }

  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "center",
    valign = "center",
    wrapText = TRUE
  )

  openxlsx::addStyle(
    wb,
    "eTable 8",
    header_style,
    rows = 1,
    cols = seq_len(ncol(etable8_display)),
    gridExpand = TRUE
  )

  openxlsx::setColWidths(
    wb,
    "eTable 8",
    cols = seq_len(ncol(etable8_display)),
    widths = "auto"
  )

  openxlsx::saveWorkbook(
    wb,
    file.path(
      OUTPUT_DIR,
      "etable8_full_reproduction.xlsx"
    ),
    overwrite = TRUE
  )
}

cat("\nFINAL eTABLE 8\n")
print(as_tibble(etable8_display), n = Inf, width = Inf)

cat("\nAnalysis complete.\n")
cat("Results saved in:\n", OUTPUT_DIR, "\n")
