# B=200 bootstrap internal validation for Uno C
# 20 multiply imputed datasets

library(readxl)
library(dplyr)
library(survival)
library(tibble)

DATA_PATH <- "KFACS_analytic_dataset_MICE_v2.xlsx"
SHEET_NAME <- "analytic_dataset_MI"
B <- 200
MASTER_SEED <- 20260819

OUTPUT_DIR <- file.path(dirname(DATA_PATH), "KFACS_incremental_prognostic")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

dat <- read_excel(DATA_PATH, sheet = SHEET_NAME) %>%
  filter(!is.na(cohort_start_year)) %>%
  mutate(
    sex = factor(sex),
    income_level = factor(income_level),
    smoking_history_ge100 = factor(smoking_history_ge100),
    ever_alcohol_use = factor(ever_alcohol_use),
    mna_scr_gr = factor(mna_scr_gr),
    frailty_3cat = factor(frailty_3cat),
    fd_lt10 = as.integer(functional_dentition < 10)
  )

mi_list <- split(dat, dat$`_imputation_id`)

cat("Number of imputations =", length(mi_list), "\n")
cat("N in imputation 1 =", nrow(mi_list[[1]]), "\n")
cat("Deaths in imputation 1 =", sum(mi_list[[1]]$death_event == 1, na.rm = TRUE), "\n")

f_base <- Surv(followup_years, death_event) ~
  age + sex + edu + income_level +
  smoking_history_ge100 + ever_alcohol_use +
  systemic_disease_count + bmi + mna_scr_gr +
  frailty_3cat + mmse_kc_score +
  systemic_reserve_adverse_count

f_fd   <- update(f_base, . ~ . + functional_dentition)
f_nt   <- update(f_base, . ~ . + natural_teeth)
f_fd10 <- update(f_base, . ~ . + fd_lt10)

extended_models <- list(
  "Continuous functional dentition" = f_fd,
  "Continuous natural teeth" = f_nt,
  "Functional dentition <10" = f_fd10
)

get_uno_c <- function(fit, newdata = NULL) {
  cc <- if (is.null(newdata)) {
    survival::concordance(fit, timewt = "n/G2")
  } else {
    survival::concordance(fit, newdata = newdata, timewt = "n/G2")
  }
  as.numeric(cc$concordance)
}

# Apparent performance
apparent_results <- list()

for (i in seq_along(mi_list)) {
  d <- mi_list[[i]]

  fit_base <- coxph(f_base, data = d, x = TRUE, y = TRUE, model = TRUE)
  C_base <- get_uno_c(fit_base)

  for (j in seq_along(extended_models)) {
    model_name <- names(extended_models)[j]
    fit_ext <- coxph(extended_models[[j]], data = d, x = TRUE, y = TRUE, model = TRUE)
    C_ext <- get_uno_c(fit_ext)

    apparent_results[[length(apparent_results) + 1]] <- data.frame(
      Imputation = i,
      Analysis = model_name,
      Apparent_Base_UnoC = C_base,
      Apparent_Extended_UnoC = C_ext,
      Apparent_Delta_UnoC = C_ext - C_base
    )
  }
}

apparent_results <- bind_rows(apparent_results)

apparent_summary <- apparent_results %>%
  group_by(Analysis) %>%
  summarise(
    Apparent_Base_UnoC = mean(Apparent_Base_UnoC, na.rm = TRUE),
    Apparent_Extended_UnoC = mean(Apparent_Extended_UnoC, na.rm = TRUE),
    Apparent_Delta_UnoC = mean(Apparent_Delta_UnoC, na.rm = TRUE),
    .groups = "drop"
  )

print(apparent_summary, n = Inf, width = Inf)

# Bootstrap optimism correction
base_boot_all <- list()
extended_boot_all <- list()

for (i in seq_along(mi_list)) {
  d <- mi_list[[i]]
  n <- nrow(d)
  set.seed(MASTER_SEED + i)

  cat("\nIMPUTATION", i, "OF", length(mi_list), "\n")

  for (b in seq_len(B)) {
    if (b == 1 || b %% 20 == 0) cat("Bootstrap", b, "of", B, "\n")

    idx <- sample(seq_len(n), size = n, replace = TRUE)
    boot <- d[idx, , drop = FALSE]

    base_result <- tryCatch({
      fit_base_boot <- coxph(f_base, data = boot, x = TRUE, y = TRUE, model = TRUE)
      C_base_boot <- get_uno_c(fit_base_boot)
      C_base_test <- get_uno_c(fit_base_boot, newdata = d)

      if (!is.finite(C_base_boot) || !is.finite(C_base_test)) stop("Non-finite base Uno C")

      data.frame(
        Imputation = i,
        Bootstrap = b,
        C_Base_Boot = C_base_boot,
        C_Base_Test = C_base_test,
        Optimism_Base = C_base_boot - C_base_test
      )
    }, error = function(e) NULL)

    if (is.null(base_result)) next

    base_boot_all[[length(base_boot_all) + 1]] <- base_result

    for (j in seq_along(extended_models)) {
      model_name <- names(extended_models)[j]

      ext_result <- tryCatch({
        fit_ext_boot <- coxph(extended_models[[j]], data = boot, x = TRUE, y = TRUE, model = TRUE)
        C_ext_boot <- get_uno_c(fit_ext_boot)
        C_ext_test <- get_uno_c(fit_ext_boot, newdata = d)

        if (!is.finite(C_ext_boot) || !is.finite(C_ext_test)) stop("Non-finite extended Uno C")

        Delta_boot <- C_ext_boot - base_result$C_Base_Boot
        Delta_test <- C_ext_test - base_result$C_Base_Test

        data.frame(
          Imputation = i,
          Bootstrap = b,
          Analysis = model_name,
          C_Base_Boot = base_result$C_Base_Boot,
          C_Base_Test = base_result$C_Base_Test,
          C_Extended_Boot = C_ext_boot,
          C_Extended_Test = C_ext_test,
          Optimism_Extended = C_ext_boot - C_ext_test,
          Delta_Boot = Delta_boot,
          Delta_Test = Delta_test,
          Optimism_Delta = Delta_boot - Delta_test
        )
      }, error = function(e) NULL)

      if (!is.null(ext_result)) {
        extended_boot_all[[length(extended_boot_all) + 1]] <- ext_result
      }
    }
  }

  saveRDS(
    list(base = base_boot_all, extended = extended_boot_all),
    file.path(OUTPUT_DIR, paste0("FINAL_B200_checkpoint_after_imputation_", i, ".rds"))
  )
}

base_boot_results <- bind_rows(base_boot_all)
extended_boot_results <- bind_rows(extended_boot_all)

base_optimism_by_imp <- base_boot_results %>%
  group_by(Imputation) %>%
  summarise(
    Mean_Optimism_Base = mean(Optimism_Base, na.rm = TRUE),
    Successful_Base_Bootstraps = n(),
    .groups = "drop"
  )

extended_optimism_by_imp <- extended_boot_results %>%
  group_by(Imputation, Analysis) %>%
  summarise(
    Mean_Optimism_Extended = mean(Optimism_Extended, na.rm = TRUE),
    Mean_Optimism_Delta = mean(Optimism_Delta, na.rm = TRUE),
    Successful_Bootstraps = n(),
    .groups = "drop"
  )

validation_by_imp <- apparent_results %>%
  left_join(base_optimism_by_imp, by = "Imputation") %>%
  left_join(extended_optimism_by_imp, by = c("Imputation", "Analysis")) %>%
  mutate(
    Corrected_Base_UnoC = Apparent_Base_UnoC - Mean_Optimism_Base,
    Corrected_Extended_UnoC = Apparent_Extended_UnoC - Mean_Optimism_Extended,
    Corrected_Delta_UnoC = Apparent_Delta_UnoC - Mean_Optimism_Delta
  )

bootstrap_summary <- validation_by_imp %>%
  group_by(Analysis) %>%
  summarise(
    Apparent_Base_UnoC = mean(Apparent_Base_UnoC, na.rm = TRUE),
    Corrected_Base_UnoC = mean(Corrected_Base_UnoC, na.rm = TRUE),
    Apparent_Extended_UnoC = mean(Apparent_Extended_UnoC, na.rm = TRUE),
    Corrected_Extended_UnoC = mean(Corrected_Extended_UnoC, na.rm = TRUE),
    Apparent_Delta_UnoC = mean(Apparent_Delta_UnoC, na.rm = TRUE),
    Mean_Optimism_Delta = mean(Mean_Optimism_Delta, na.rm = TRUE),
    Corrected_Delta_UnoC = mean(Corrected_Delta_UnoC, na.rm = TRUE),
    Mean_Successful_Bootstraps = mean(Successful_Bootstraps, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    row_order = case_when(
      Analysis == "Continuous functional dentition" ~ 1,
      Analysis == "Continuous natural teeth" ~ 2,
      Analysis == "Functional dentition <10" ~ 3,
      TRUE ~ 99
    )
  ) %>%
  arrange(row_order) %>%
  select(-row_order)

options(pillar.sigfig = 7)

cat("\nFINAL B = 200 BOOTSTRAP OPTIMISM-CORRECTED UNO C\n")
print(as_tibble(bootstrap_summary), n = Inf, width = Inf)

bootstrap_display <- bootstrap_summary %>%
  mutate(
    Apparent_Base_UnoC = sprintf("%.4f", Apparent_Base_UnoC),
    Corrected_Base_UnoC = sprintf("%.4f", Corrected_Base_UnoC),
    Apparent_Extended_UnoC = sprintf("%.4f", Apparent_Extended_UnoC),
    Corrected_Extended_UnoC = sprintf("%.4f", Corrected_Extended_UnoC),
    Apparent_Delta_UnoC = sprintf("%+.5f", Apparent_Delta_UnoC),
    Mean_Optimism_Delta = sprintf("%+.5f", Mean_Optimism_Delta),
    Corrected_Delta_UnoC = sprintf("%+.5f", Corrected_Delta_UnoC),
    Mean_Successful_Bootstraps = sprintf("%.1f", Mean_Successful_Bootstraps)
  )

print(as_tibble(bootstrap_display), n = Inf, width = Inf)

write.csv(apparent_results,
          file.path(OUTPUT_DIR, "FINAL_B200_apparent_UnoC_by_imputation.csv"),
          row.names = FALSE)

write.csv(apparent_summary,
          file.path(OUTPUT_DIR, "FINAL_B200_apparent_UnoC_summary.csv"),
          row.names = FALSE)

write.csv(base_boot_results,
          file.path(OUTPUT_DIR, "FINAL_B200_base_bootstrap_all.csv"),
          row.names = FALSE)

write.csv(extended_boot_results,
          file.path(OUTPUT_DIR, "FINAL_B200_extended_bootstrap_all.csv"),
          row.names = FALSE)

write.csv(validation_by_imp,
          file.path(OUTPUT_DIR, "FINAL_B200_validation_by_imputation.csv"),
          row.names = FALSE)

write.csv(bootstrap_summary,
          file.path(OUTPUT_DIR, "FINAL_B200_optimism_corrected_UnoC_summary.csv"),
          row.names = FALSE)

write.csv(bootstrap_display,
          file.path(OUTPUT_DIR, "FINAL_B200_optimism_corrected_UnoC_display.csv"),
          row.names = FALSE)

cat("\nANALYSIS COMPLETE\n")
cat("Results saved in:\n", OUTPUT_DIR, "\n")
