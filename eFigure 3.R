# ============================================================
# eFigure 3
# Sequential pathway attenuation
#
# - No figure title
# - Blue dots: actual HRs from 20 imputed datasets
# - Violin: distribution of the 20 actual HRs
# - Black dots and line: Rubin-pooled HRs
# - Large text and white halo
# ============================================================


# ============================================================
# 0. 필요한 패키지 설치 및 로딩
# ============================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "survival",
  "ggplot2"
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

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(survival)
  library(ggplot2)
})


# ============================================================
# 1. 파일 경로 설정
# ============================================================

file_path <- "KFACS_analytic_dataset_MICE_v2.xlsx"
output_folder <- "."

dir.create(
  output_folder,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Excel 데이터 불러오기
# ============================================================

mi_long <- read_excel(
  path = file_path,
  sheet = "analytic_dataset_MI"
)

cat("Excel 데이터 불러오기 완료\n")
cat("전체 행 수:", nrow(mi_long), "\n")
cat(
  "다중대치 데이터셋 수:",
  length(unique(mi_long$`_imputation_id`)),
  "\n\n"
)


# ============================================================
# 3. 필요한 변수 존재 여부 확인
# ============================================================

required_variables <- c(
  "_imputation_id",
  "functional_dentition_lt10",
  "cohort_start_year",
  "center",
  "sex",
  "age",
  "edu",
  "income_level",
  "smoking_history_ge100",
  "ever_alcohol_use",
  "chs_wtloss",
  "bmi",
  "mna_scr_gr",
  "mmse_kc_score",
  "systemic_disease_count",
  "chs_total",
  "systemic_reserve_adverse_count",
  "death_event",
  "followup_years"
)

absent_variables <- setdiff(
  required_variables,
  names(mi_long)
)

if (length(absent_variables) > 0) {
  stop(
    paste0(
      "Excel에 다음 변수가 없습니다: ",
      paste(absent_variables, collapse = ", ")
    )
  )
}


# ============================================================
# 4. 분석용 변수 생성
# ============================================================

mi_analysis <- mi_long %>%
  
  # cohort year를 결정할 수 없는 3개 기록 제외
  filter(!is.na(cohort_start_year)) %>%
  
  mutate(
    death_event = as.integer(death_event),
    followup_years = as.numeric(followup_years),
    
    # Primary exposure:
    # functional dentition <10 vs ≥10
    fd_lt10 = as.integer(
      functional_dentition_lt10
    ),
    
    # CHS weight-loss coding
    # 원자료: 1 = present, 2 = absent
    weight_loss_flag = if_else(
      chs_wtloss == 1,
      1,
      0
    ),
    
    # MNA screening risk
    # 1 = normal
    # 2 또는 3 = nutritional risk
    mna_risk_flag = if_else(
      mna_scr_gr >= 2,
      1,
      0
    ),
    
    # Weight loss를 제외한 CHS frailty burden
    frailty_no_weight_loss = pmax(
      chs_total - weight_loss_flag,
      0
    ),
    
    # 범주형 보정변수
    sex_factor = factor(sex),
    
    # Center가 결측인 3명은 Unknown level로 유지
    # 이를 통해 각 대치본 N = 2,715 유지
    center_factor = factor(
      ifelse(
        is.na(center),
        "Unknown",
        as.character(center)
      )
    ),
    
    cohort_year_factor = factor(
      cohort_start_year
    ),
    
    income_factor = factor(
      income_level,
      levels = c(
        "High",
        "Middle",
        "Low/no income"
      )
    )
  )


# ============================================================
# 5. 각 대치본의 분석대상 수 확인
# ============================================================

imputation_counts <- mi_analysis %>%
  count(`_imputation_id`)

print(imputation_counts)

cat(
  "\n각 대치본의 최소 N:",
  min(imputation_counts$n),
  "\n"
)

cat(
  "각 대치본의 최대 N:",
  max(imputation_counts$n),
  "\n\n"
)

if (
  min(imputation_counts$n) != 2715 ||
  max(imputation_counts$n) != 2715
) {
  warning(
    "각 대치본의 분석대상 수가 모두 2,715인지 확인하세요."
  )
}


# ============================================================
# 6. 분석변수 결측 확인
# ============================================================

analysis_variables <- c(
  "followup_years",
  "death_event",
  "fd_lt10",
  "age",
  "sex_factor",
  "center_factor",
  "cohort_year_factor",
  "edu",
  "income_factor",
  "smoking_history_ge100",
  "ever_alcohol_use",
  "weight_loss_flag",
  "bmi",
  "mna_risk_flag",
  "mmse_kc_score",
  "systemic_disease_count",
  "frailty_no_weight_loss",
  "systemic_reserve_adverse_count"
)

missing_check <- sapply(
  mi_analysis[analysis_variables],
  function(x) sum(is.na(x))
)

cat("분석변수별 결측 수\n")
print(missing_check)

if (any(missing_check > 0)) {
  stop(
    paste0(
      "다음 변수에 결측이 남아 있습니다: ",
      paste(
        names(missing_check)[missing_check > 0],
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 7. Sequential pathway model 정의
# ============================================================

model_terms <- list(
  
  "1. Demographic" = c(
    "age",
    "sex_factor",
    "center_factor",
    "cohort_year_factor"
  ),
  
  "2. +Socioeconomic/\nbehavioral" = c(
    "age",
    "sex_factor",
    "center_factor",
    "cohort_year_factor",
    "edu",
    "income_factor",
    "smoking_history_ge100",
    "ever_alcohol_use"
  ),
  
  "3. +Weight loss/\nnutrition" = c(
    "age",
    "sex_factor",
    "center_factor",
    "cohort_year_factor",
    "edu",
    "income_factor",
    "smoking_history_ge100",
    "ever_alcohol_use",
    "weight_loss_flag",
    "bmi",
    "mna_risk_flag"
  ),
  
  "4. +Cognition" = c(
    "age",
    "sex_factor",
    "center_factor",
    "cohort_year_factor",
    "edu",
    "income_factor",
    "smoking_history_ge100",
    "ever_alcohol_use",
    "weight_loss_flag",
    "bmi",
    "mna_risk_flag",
    "mmse_kc_score"
  ),
  
  "5. +Systemic\ndisease" = c(
    "age",
    "sex_factor",
    "center_factor",
    "cohort_year_factor",
    "edu",
    "income_factor",
    "smoking_history_ge100",
    "ever_alcohol_use",
    "weight_loss_flag",
    "bmi",
    "mna_risk_flag",
    "mmse_kc_score",
    "systemic_disease_count"
  ),
  
  "6A. +Frailty &\nreserve (full)" = c(
    "age",
    "sex_factor",
    "center_factor",
    "cohort_year_factor",
    "edu",
    "income_factor",
    "smoking_history_ge100",
    "ever_alcohol_use",
    "weight_loss_flag",
    "bmi",
    "mna_risk_flag",
    "mmse_kc_score",
    "systemic_disease_count",
    "frailty_no_weight_loss",
    "systemic_reserve_adverse_count"
  )
)

model_order <- names(model_terms)


# ============================================================
# 8. 각 대치본에서 Cox model 적합
# ============================================================

fit_imputation_models <- function(
    data,
    model_terms
) {
  
  imputation_ids <- sort(
    unique(data$`_imputation_id`)
  )
  
  all_results <- list()
  result_index <- 1
  
  for (model_name in names(model_terms)) {
    
    adjustment_terms <- model_terms[[model_name]]
    
    model_formula <- as.formula(
      paste0(
        "Surv(followup_years, death_event) ~ ",
        "fd_lt10 + ",
        paste(
          adjustment_terms,
          collapse = " + "
        )
      )
    )
    
    for (imp_id in imputation_ids) {
      
      one_imp <- data %>%
        filter(
          `_imputation_id` == imp_id
        )
      
      fit <- coxph(
        formula = model_formula,
        data = one_imp,
        ties = "efron",
        na.action = na.fail
      )
      
      beta_value <- unname(
        coef(fit)["fd_lt10"]
      )
      
      variance_value <- unname(
        vcov(fit)[
          "fd_lt10",
          "fd_lt10"
        ]
      )
      
      all_results[[result_index]] <- data.frame(
        Model = model_name,
        Imputation = imp_id,
        Beta = beta_value,
        Variance = variance_value,
        HR = exp(beta_value),
        N = fit$n,
        Events = fit$nevent
      )
      
      result_index <- result_index + 1
    }
  }
  
  bind_rows(all_results)
}


imputation_results <- fit_imputation_models(
  data = mi_analysis,
  model_terms = model_terms
)

imputation_results$Model <- factor(
  imputation_results$Model,
  levels = model_order
)


# 각 모델에 실제 점이 20개인지 확인
dot_check <- imputation_results %>%
  count(Model)

cat("\n각 모델별 실제 HR 점 개수\n")
print(dot_check)


# ============================================================
# 9. Rubin's rules pooling 함수
# ============================================================

pool_one_model <- function(
    beta_values,
    variance_values
) {
  
  m <- length(beta_values)
  
  pooled_beta <- mean(
    beta_values
  )
  
  within_variance <- mean(
    variance_values
  )
  
  between_variance <- var(
    beta_values
  )
  
  total_variance <- (
    within_variance +
      (1 + 1 / m) *
      between_variance
  )
  
  pooled_se <- sqrt(
    total_variance
  )
  
  # Rubin 자유도
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
  
  critical_value <- if (
    is.finite(degrees_freedom)
  ) {
    
    qt(
      0.975,
      df = degrees_freedom
    )
    
  } else {
    
    qnorm(0.975)
  }
  
  test_statistic <- (
    pooled_beta /
      pooled_se
  )
  
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
  
  data.frame(
    Beta = pooled_beta,
    SE = pooled_se,
    DF = degrees_freedom,
    
    HR = exp(
      pooled_beta
    ),
    
    Lower95 = exp(
      pooled_beta -
        critical_value *
        pooled_se
    ),
    
    Upper95 = exp(
      pooled_beta +
        critical_value *
        pooled_se
    ),
    
    P_value = p_value
  )
}


# ============================================================
# 10. 각 model의 pooled 결과 계산
# ============================================================

pooled_list <- lapply(
  model_order,
  function(model_name) {
    
    model_data <- imputation_results %>%
      filter(
        Model == model_name
      )
    
    pooled_one <- pool_one_model(
      beta_values = model_data$Beta,
      variance_values = model_data$Variance
    )
    
    pooled_one$Model <- model_name
    
    pooled_one
  }
)

pooled_results <- bind_rows(
  pooled_list
)

pooled_results$Model <- factor(
  pooled_results$Model,
  levels = model_order
)


# ============================================================
# 11. Attenuation 계산
# ============================================================

reference_beta <- pooled_results$Beta[
  pooled_results$Model == model_order[1]
]

pooled_results <- pooled_results %>%
  mutate(
    Attenuation_percent = (
      reference_beta - Beta
    ) / reference_beta * 100
  )


# ============================================================
# 12. 벡터화된 P-value 표시 함수
# ============================================================

format_p_value <- function(p) {
  
  ifelse(
    p < 0.001,
    "P<.001",
    paste0(
      "P=",
      sub(
        "^0",
        "",
        sprintf("%.3f", p)
      )
    )
  )
}


# ============================================================
# 13. HR 및 P-value annotation 생성
# ============================================================

pooled_results <- pooled_results %>%
  mutate(
    P_text = format_p_value(
      P_value
    ),
    
    Annotation = paste0(
      "HR=",
      sprintf("%.2f", HR),
      "\n",
      P_text
    ),
    
    Annotation_y = HR + 0.095
  )


# ============================================================
# 14. 결과 확인
# ============================================================

cat("\nRubin-pooled pathway 결과\n")

print(
  pooled_results %>%
    select(
      Model,
      HR,
      Lower95,
      Upper95,
      P_value,
      Attenuation_percent
    ) %>%
    mutate(
      across(
        where(is.numeric),
        ~ round(.x, 4)
      )
    )
)


# ============================================================
# 15. 결과 CSV 저장
# ============================================================

write.csv(
  imputation_results,
  file = file.path(
    output_folder,
    "eFigure3_actual_imputation_HRs.csv"
  ),
  row.names = FALSE
)

write.csv(
  pooled_results,
  file = file.path(
    output_folder,
    "eFigure3_pooled_pathway_results.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 16. Y축 상한 설정
# ============================================================

y_upper <- max(
  c(
    imputation_results$HR,
    pooled_results$Annotation_y
  ),
  na.rm = TRUE
) + 0.22

y_upper <- max(
  y_upper,
  1.90
)


# ============================================================
# 17. eFigure 3 생성
# ============================================================

set.seed(20260801)

eFigure_3 <- ggplot() +
  
  # 실제 20개 HR의 분포
  geom_violin(
    data = imputation_results,
    aes(
      x = Model,
      y = HR
    ),
    fill = "#B9C9E3",
    color = NA,
    alpha = 0.66,
    width = 0.75,
    trim = TRUE,
    scale = "width",
    adjust = 1.1
  ) +
  
  # 20개 대치본의 실제 HR
  geom_point(
    data = imputation_results,
    aes(
      x = Model,
      y = HR
    ),
    position = position_jitter(
      width = 0.085,
      height = 0,
      seed = 20260801
    ),
    shape = 21,
    size = 3.4,
    stroke = 0.65,
    fill = "#527DB8",
    color = "#315D96",
    alpha = 0.85
  ) +
  
  # HR = 1 기준선
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.80,
    color = "grey35"
  ) +
  
  # Rubin-pooled HR 연결선
  geom_line(
    data = pooled_results,
    aes(
      x = Model,
      y = HR,
      group = 1
    ),
    linewidth = 1.35,
    color = "black"
  ) +
  
  # Rubin-pooled HR 점
  geom_point(
    data = pooled_results,
    aes(
      x = Model,
      y = HR
    ),
    size = 5.2,
    shape = 21,
    fill = "black",
    color = "black"
  ) +
  
  # HR 및 P-value 흰색 halo
  geom_text(
    data = pooled_results,
    aes(
      x = Model,
      y = Annotation_y,
      label = Annotation
    ),
    size = 7.0,
    fontface = "bold",
    lineheight = 0.92,
    color = "white"
  ) +
  
  # HR 및 P-value 검은색 글자
  geom_text(
    data = pooled_results,
    aes(
      x = Model,
      y = Annotation_y,
      label = Annotation
    ),
    size = 5.6,
    fontface = "bold",
    lineheight = 0.92,
    color = "black"
  ) +
  
  scale_x_discrete(
    drop = FALSE
  ) +
  
  scale_y_continuous(
    breaks = seq(
      1.0,
      ceiling(y_upper * 10) / 10,
      by = 0.2
    ),
    labels = function(x) {
      sprintf("%.1f", x)
    },
    expand = expansion(
      mult = c(
        0.02,
        0.06
      )
    )
  ) +
  
  coord_cartesian(
    ylim = c(
      0.95,
      y_upper
    ),
    clip = "off"
  ) +
  
  labs(
    title = NULL,
    subtitle = NULL,
    x = NULL,
    y = "Hazard ratio"
  ) +
  
  theme_bw(
    base_size = 18
  ) +
  
  theme(
    # Title 완전 제거
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    
    # Y축 제목
    axis.title.y = element_text(
      size = 21,
      color = "black",
      margin = margin(
        r = 14
      )
    ),
    
    # Y축 숫자
    axis.text.y = element_text(
      size = 18,
      color = "black"
    ),
    
    # X축 model 이름
    axis.text.x = element_text(
      size = 17,
      color = "black",
      lineheight = 0.95,
      margin = margin(
        t = 11
      )
    ),
    
    axis.ticks = element_line(
      linewidth = 0.75,
      color = "black"
    ),
    
    axis.ticks.length = grid::unit(
      0.22,
      "cm"
    ),
    
    panel.grid.major.x = element_blank(),
    
    panel.grid.major.y = element_line(
      color = "grey88",
      linewidth = 0.55
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "black",
      linewidth = 1.1
    ),
    
    plot.margin = margin(
      t = 24,
      r = 24,
      b = 22,
      l = 22
    )
  )


# Figure 화면 출력
print(eFigure_3)


# ============================================================
# 18. 고해상도 Figure 저장
# ============================================================

# PNG 900 dpi
ggsave(
  filename = file.path(
    output_folder,
    "eFigure3_actual_distribution_large_text.png"
  ),
  plot = eFigure_3,
  width = 15,
  height = 8,
  units = "in",
  dpi = 900,
  bg = "white"
)

# TIFF 900 dpi
ggsave(
  filename = file.path(
    output_folder,
    "eFigure3_actual_distribution_large_text.tiff"
  ),
  plot = eFigure_3,
  width = 15,
  height = 8,
  units = "in",
  dpi = 900,
  compression = "lzw",
  bg = "white"
)

# 벡터 PDF
ggsave(
  filename = file.path(
    output_folder,
    "eFigure3_actual_distribution_large_text.pdf"
  ),
  plot = eFigure_3,
  width = 15,
  height = 8,
  units = "in",
  bg = "white"
)


cat(
  "\neFigure 3 생성 완료\n",
  "파란 점: 실제 20개 imputation-specific HR\n",
  "검은 점: Rubin-pooled HR\n",
  "저장 폴더:\n",
  output_folder,
  "\n"
)