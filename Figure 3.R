# ============================================================
# eFigure 2
# Deficit-based dose-response association below 20 units
# Two-panel figure: Natural teeth / Functional dentition
# ============================================================


# ============================================================
# 0. 필요한 패키지 설치 및 로딩
# ============================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "survival",
  "ggplot2",
  "patchwork"
)

missing_packages <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )
}

library(readxl)
library(dplyr)
library(survival)
library(ggplot2)
library(patchwork)


# ============================================================
# 1. Excel 데이터 불러오기
# ============================================================

file_path <- "KFACS_analytic_dataset_MICE_v2.xlsx"

mi_long <- read_excel(
  path = file_path,
  sheet = "analytic_dataset_MI"
)

cat("전체 행 수:", nrow(mi_long), "\n")
cat(
  "다중대치 자료 수:",
  length(unique(mi_long$`_imputation_id`)),
  "\n"
)


# ============================================================
# 2. 분석 변수 정리
# ============================================================

mi_analysis <- mi_long %>%
  
  # cohort year가 확인되지 않는 3개 기록 제외
  filter(!is.na(cohort_start_year)) %>%
  
  mutate(
    death_event = as.integer(death_event),
    followup_years = as.numeric(followup_years),
    
    # Deficit below 20
    # Count가 20 이상이면 deficit = 0
    # Count가 20 미만이면 deficit = 20 - count
    natural_teeth_deficit = pmax(
      20 - natural_teeth,
      0
    ),
    
    functional_dentition_deficit = pmax(
      20 - functional_dentition,
      0
    ),
    
    # Cox model 범주형 변수
    sex_factor = factor(sex),
    
    center_factor = factor(center),
    
    cohort_year_factor = factor(
      cohort_start_year
    ),
    
    income_level_factor = factor(
      income_level,
      levels = c(
        "High",
        "Middle",
        "Low/no income"
      )
    )
  )


# 각 대치본에 2,715명이 있는지 확인
imputation_counts <- mi_analysis %>%
  count(`_imputation_id`)

print(imputation_counts)

cat(
  "각 대치본의 최소 인원:",
  min(imputation_counts$n),
  "\n"
)

cat(
  "각 대치본의 최대 인원:",
  max(imputation_counts$n),
  "\n"
)


# ============================================================
# 3. Cox Model 2에 필요한 변수의 결측 확인
# ============================================================

model_variables <- c(
  "followup_years",
  "death_event",
  "natural_teeth_deficit",
  "functional_dentition_deficit",
  "age",
  "sex_factor",
  "center_factor",
  "cohort_year_factor",
  "edu",
  "income_level_factor",
  "smoking_history_ge100",
  "ever_alcohol_use",
  "systemic_disease_count",
  "systemic_reserve_adverse_count"
)

missing_check <- sapply(
  mi_analysis[model_variables],
  function(x) sum(is.na(x))
)

print(missing_check)


# ============================================================
# 4. 각 대치본에서 Cox model을 적합하고
#    Rubin's rules로 통합하는 함수
# ============================================================

pool_cox_exposure <- function(
    data,
    exposure_name
) {
  
  imputation_ids <- sort(
    unique(data$`_imputation_id`)
  )
  
  result_list <- lapply(
    imputation_ids,
    function(imp_id) {
      
      one_imp <- data %>%
        filter(`_imputation_id` == imp_id)
      
      model_formula <- as.formula(
        paste0(
          "Surv(followup_years, death_event) ~ ",
          exposure_name,
          " + age",
          " + sex_factor",
          " + center_factor",
          " + cohort_year_factor",
          " + edu",
          " + income_level_factor",
          " + smoking_history_ge100",
          " + ever_alcohol_use",
          " + systemic_disease_count",
          " + systemic_reserve_adverse_count"
        )
      )
      
      fit <- coxph(
        formula = model_formula,
        data = one_imp,
        ties = "efron",
        x = TRUE,
        y = TRUE
      )
      
      beta_value <- unname(
        coef(fit)[exposure_name]
      )
      
      variance_value <- unname(
        vcov(fit)[
          exposure_name,
          exposure_name
        ]
      )
      
      data.frame(
        imputation = imp_id,
        beta = beta_value,
        variance = variance_value,
        n = fit$n,
        events = fit$nevent
      )
    }
  )
  
  individual_results <- bind_rows(
    result_list
  )
  
  m <- nrow(individual_results)
  
  # Rubin's rules
  pooled_beta <- mean(
    individual_results$beta
  )
  
  within_variance <- mean(
    individual_results$variance
  )
  
  between_variance <- var(
    individual_results$beta
  )
  
  total_variance <- (
    within_variance +
      (1 + 1 / m) * between_variance
  )
  
  pooled_se <- sqrt(
    total_variance
  )
  
  # Rubin's rules에 따른 자유도
  if (
    is.na(between_variance) ||
    between_variance < .Machine$double.eps
  ) {
    
    degrees_freedom <- Inf
    
  } else {
    
    relative_increase <- (
      (1 + 1 / m) *
        between_variance /
        within_variance
    )
    
    degrees_freedom <- (
      (m - 1) *
        (1 + 1 / relative_increase)^2
    )
  }
  
  critical_value <- ifelse(
    is.finite(degrees_freedom),
    qt(
      0.975,
      df = degrees_freedom
    ),
    qnorm(0.975)
  )
  
  test_statistic <- pooled_beta / pooled_se
  
  p_value <- if (
    is.finite(degrees_freedom)
  ) {
    
    2 * pt(
      abs(test_statistic),
      df = degrees_freedom,
      lower.tail = FALSE
    )
    
  } else {
    
    2 * pnorm(
      abs(test_statistic),
      lower.tail = FALSE
    )
  }
  
  # 1-unit deficit HR
  hr_1 <- exp(
    pooled_beta
  )
  
  lower_1 <- exp(
    pooled_beta -
      critical_value * pooled_se
  )
  
  upper_1 <- exp(
    pooled_beta +
      critical_value * pooled_se
  )
  
  # 5-unit deficit HR
  hr_5 <- exp(
    5 * pooled_beta
  )
  
  lower_5 <- exp(
    5 * pooled_beta -
      critical_value * 5 * pooled_se
  )
  
  upper_5 <- exp(
    5 * pooled_beta +
      critical_value * 5 * pooled_se
  )
  
  list(
    beta = pooled_beta,
    se = pooled_se,
    df = degrees_freedom,
    critical_value = critical_value,
    p = p_value,
    
    hr_1 = hr_1,
    lower_1 = lower_1,
    upper_1 = upper_1,
    
    hr_5 = hr_5,
    lower_5 = lower_5,
    upper_5 = upper_5,
    
    within_variance = within_variance,
    between_variance = between_variance,
    total_variance = total_variance,
    
    individual_results = individual_results
  )
}


# ============================================================
# 5. Natural teeth deficit model
# ============================================================

natural_result <- pool_cox_exposure(
  data = mi_analysis,
  exposure_name = "natural_teeth_deficit"
)


# ============================================================
# 6. Functional dentition deficit model
# ============================================================

functional_result <- pool_cox_exposure(
  data = mi_analysis,
  exposure_name =
    "functional_dentition_deficit"
)


# ============================================================
# 7. Pooled Cox 결과 확인
# ============================================================

pooled_results <- data.frame(
  Exposure = c(
    "Natural teeth deficit below 20",
    "Functional dentition deficit below 20"
  ),
  
  HR_per_1_unit = c(
    natural_result$hr_1,
    functional_result$hr_1
  ),
  
  Lower_95CI_per_1_unit = c(
    natural_result$lower_1,
    functional_result$lower_1
  ),
  
  Upper_95CI_per_1_unit = c(
    natural_result$upper_1,
    functional_result$upper_1
  ),
  
  P_value = c(
    natural_result$p,
    functional_result$p
  ),
  
  HR_per_5_units = c(
    natural_result$hr_5,
    functional_result$hr_5
  ),
  
  Lower_95CI_per_5_units = c(
    natural_result$lower_5,
    functional_result$lower_5
  ),
  
  Upper_95CI_per_5_units = c(
    natural_result$upper_5,
    functional_result$upper_5
  )
)

print(
  pooled_results %>%
    mutate(
      across(
        where(is.numeric),
        ~ round(.x, 4)
      )
    )
)


# ============================================================
# 8. 결과 테이블 CSV 저장
# ============================================================

output_folder <- "."

write.csv(
  pooled_results,
  file = file.path(
    output_folder,
    "eFigure2_pooled_Cox_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 9. 그래프용 predicted HR 데이터 생성
# ============================================================

make_prediction_data <- function(
    pooled_result
) {
  
  count_values <- seq(
    from = 0,
    to = 32,
    by = 0.05
  )
  
  deficit_values <- pmax(
    20 - count_values,
    0
  )
  
  beta <- pooled_result$beta
  se <- pooled_result$se
  critical_value <- pooled_result$critical_value
  
  data.frame(
    count = count_values,
    deficit = deficit_values,
    
    hr = exp(
      beta * deficit_values
    ),
    
    lower = exp(
      beta * deficit_values -
        critical_value *
        se *
        deficit_values
    ),
    
    upper = exp(
      beta * deficit_values +
        critical_value *
        se *
        deficit_values
    )
  )
}

natural_prediction <- make_prediction_data(
  natural_result
)

functional_prediction <- make_prediction_data(
  functional_result
)


# ============================================================
# 10. P-value 표시 함수
# ============================================================

format_p_value <- function(p) {
  
  if (p < 0.001) {
    return("P < .001")
  }
  
  paste0(
    "P = ",
    sub(
      "^0",
      "",
      sprintf("%.3f", p)
    )
  )
}


# ============================================================
# 11. 패널 생성 함수
# ============================================================

make_deficit_panel <- function(
    prediction_data,
    pooled_result,
    panel_title,
    x_axis_title,
    line_color,
    ribbon_color
) {
  
  model_annotation <- paste0(
    "HR per 1-unit deficit below 20:\n",
    sprintf(
      "%.2f",
      pooled_result$hr_1
    ),
    " (",
    sprintf(
      "%.2f",
      pooled_result$lower_1
    ),
    ", ",
    sprintf(
      "%.2f",
      pooled_result$upper_1
    ),
    "), ",
    format_p_value(
      pooled_result$p
    )
  )
  
  selected_points <- data.frame(
    count = c(
      0,
      10,
      19
    )
  ) %>%
    mutate(
      deficit = pmax(
        20 - count,
        0
      ),
      
      hr = exp(
        pooled_result$beta *
          deficit
      ),
      
      label = paste0(
        count,
        ": HR ",
        sprintf(
          "%.2f",
          hr
        )
      )
    )
  
  ggplot(
    prediction_data,
    aes(
      x = count,
      y = hr
    )
  ) +
    
    geom_ribbon(
      aes(
        ymin = lower,
        ymax = upper
      ),
      fill = ribbon_color,
      alpha = 0.42
    ) +
    
    geom_line(
      color = line_color,
      linewidth = 1.15
    ) +
    
    # HR = 1
    geom_hline(
      yintercept = 1,
      linetype = "dashed",
      linewidth = 0.55,
      color = "grey25"
    ) +
    
    # Severe loss threshold
    geom_vline(
      xintercept = 10,
      linetype = "dotted",
      linewidth = 0.55,
      color = "grey25"
    ) +
    
    # Functional threshold
    geom_vline(
      xintercept = 20,
      linetype = "dashed",
      linewidth = 0.65,
      color = "grey20"
    ) +
    
    geom_point(
      data = selected_points,
      aes(
        x = count,
        y = hr
      ),
      inherit.aes = FALSE,
      size = 2.4,
      color = "grey15"
    ) +
    
    geom_text(
      data = selected_points,
      aes(
        x = count,
        y = hr,
        label = label
      ),
      inherit.aes = FALSE,
      hjust = c(
        -0.05,
        -0.05,
        1.05
      ),
      vjust = -0.75,
      size = 3.2,
      color = "grey20"
    ) +
    
    annotate(
      geom = "text",
      x = 0.7,
      y = 2.27,
      label = model_annotation,
      hjust = 0,
      vjust = 1,
      size = 3.4,
      color = "grey20"
    ) +
    
    annotate(
      geom = "text",
      x = 20.15,
      y = 0.965,
      label = "≥20: reference",
      hjust = 0,
      vjust = 1,
      size = 3.2,
      color = "grey20"
    ) +
    
    scale_x_continuous(
      limits = c(
        0,
        32
      ),
      breaks = c(
        0,
        10,
        20,
        28,
        32
      ),
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    
    scale_y_continuous(
      limits = c(
        0.65,
        2.40
      ),
      breaks = seq(
        0.8,
        2.4,
        by = 0.2
      ),
      expand = expansion(
        mult = c(
          0,
          0
        )
      )
    ) +
    
    labs(
      title = panel_title,
      x = x_axis_title,
      y = paste0(
        "Adjusted hazard ratio for ",
        "all-cause mortality"
      )
    ) +
    
    theme_bw(
      base_size = 11
    ) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        size = 13,
        hjust = 0.5,
        margin = margin(
          b = 5
        )
      ),
      
      axis.title = element_text(
        size = 10.5
      ),
      
      axis.text = element_text(
        size = 9.5,
        color = "black"
      ),
      
      panel.grid.major = element_line(
        color = "grey88",
        linewidth = 0.35
      ),
      
      panel.grid.minor = element_blank(),
      
      panel.border = element_rect(
        color = "grey20",
        linewidth = 0.7
      ),
      
      plot.margin = margin(
        t = 8,
        r = 10,
        b = 8,
        l = 8
      )
    )
}


# ============================================================
# 12. Panel A: Natural teeth
# ============================================================

panel_A <- make_deficit_panel(
  prediction_data = natural_prediction,
  pooled_result = natural_result,
  
  panel_title = "Natural teeth",
  x_axis_title = "Natural teeth count",
  
  line_color = "#00527A",
  ribbon_color = "#6F9BB5"
)


# ============================================================
# 13. Panel B: Functional dentition
# ============================================================

panel_B <- make_deficit_panel(
  prediction_data = functional_prediction,
  pooled_result = functional_result,
  
  panel_title = "Functional dentition",
  x_axis_title = "Functional dentition count",
  
  line_color = "#85492F",
  ribbon_color = "#B89C8D"
)


# ============================================================
# 14. A-B 패널 결합
# 전체 figure title은 넣지 않음
# ============================================================

eFigure_2 <- (
  panel_A | panel_B
) +
  plot_annotation(
    tag_levels = "A"
  ) &
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 17
    ),
    
    plot.tag.position = c(
      0.01,
      0.99
    )
  )

print(eFigure_2)


# ============================================================
# 15. Figure 저장
# ============================================================

ggsave(
  filename = file.path(
    output_folder,
    "eFigure2_deficit_two_panels.png"
  ),
  plot = eFigure_2,
  width = 14,
  height = 6.2,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_folder,
    "eFigure2_deficit_two_panels.tiff"
  ),
  plot = eFigure_2,
  width = 14,
  height = 6.2,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

ggsave(
  filename = file.path(
    output_folder,
    "eFigure2_deficit_two_panels.pdf"
  ),
  plot = eFigure_2,
  width = 14,
  height = 6.2,
  units = "in",
  bg = "white"
)

cat(
  "\n완료: Figure와 pooled Cox 결과가 다음 폴더에 저장되었습니다.\n",
  output_folder,
  "\n"
)