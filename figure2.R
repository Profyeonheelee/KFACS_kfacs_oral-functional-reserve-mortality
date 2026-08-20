# ============================================================
# Figure 2
# Kaplan-Meier Survival Estimates According to
# Natural Teeth and Functional Dentition
#
# - Representative imputation: 1 of 20
# - Analytic cohort: N = 2,715
# - survminer와 ggpubr는 사용하지 않음
# ============================================================


# ============================================================
# 0. 필요한 패키지 설치 및 불러오기
# ============================================================

required_packages <- c(
  "readxl",
  "dplyr",
  "survival",
  "ggplot2",
  "patchwork",
  "scales"
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
  library(patchwork)
  library(scales)
})

cat("All required packages were loaded successfully.\n")


# ============================================================
# 1. 파일 경로 및 저장 폴더
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


# Figure 2에는 대표 다중대치본 1번 사용
dat <- mi_long %>%
  filter(`_imputation_id` == 1) %>%
  filter(!is.na(cohort_start_year)) %>%
  mutate(
    death_event = as.integer(death_event),
    followup_years = as.numeric(followup_years)
  )


# 기본 정보 확인
cat("Analytic N =", nrow(dat), "\n")
cat(
  "Deaths =",
  sum(dat$death_event, na.rm = TRUE),
  "\n"
)

cat(
  "Mean follow-up =",
  round(
    mean(dat$followup_years, na.rm = TRUE),
    2
  ),
  "years\n"
)

if (nrow(dat) != 2715) {
  warning(
    paste0(
      "Analytic N이 2,715가 아닙니다. 현재 N = ",
      nrow(dat)
    )
  )
}


# ============================================================
# 3. Figure 2 분석 그룹 생성
# ============================================================

dat <- dat %>%
  mutate(
    
    # --------------------------------------------------------
    # Panel A: Natural teeth group
    # --------------------------------------------------------
    
    natural_teeth_group = case_when(
      natural_teeth >= 20 ~ "≥20",
      natural_teeth >= 10 ~ "10–19",
      natural_teeth >= 1  ~ "1–9",
      natural_teeth < 1   ~ "0",
      TRUE ~ NA_character_
    ),
    
    natural_teeth_group = factor(
      natural_teeth_group,
      levels = c(
        "≥20",
        "10–19",
        "1–9",
        "0"
      )
    ),
    
    
    # --------------------------------------------------------
    # Panel B: Functional dentition group
    # --------------------------------------------------------
    
    functional_dentition_group = case_when(
      functional_dentition >= 20 ~ "≥20",
      functional_dentition >= 10 ~ "10–19",
      functional_dentition >= 1  ~ "1–9",
      functional_dentition < 1   ~ "0",
      TRUE ~ NA_character_
    ),
    
    functional_dentition_group = factor(
      functional_dentition_group,
      levels = c(
        "≥20",
        "10–19",
        "1–9",
        "0"
      )
    ),
    
    
    # --------------------------------------------------------
    # Panel C: Severe natural tooth loss
    # --------------------------------------------------------
    
    severe_tooth_loss = case_when(
      natural_teeth < 10  ~ "<10 natural teeth",
      natural_teeth >= 10 ~ "≥10 natural teeth",
      TRUE ~ NA_character_
    ),
    
    severe_tooth_loss = factor(
      severe_tooth_loss,
      levels = c(
        "≥10 natural teeth",
        "<10 natural teeth"
      )
    ),
    
    
    # --------------------------------------------------------
    # Panel D: Functional reclassification
    # --------------------------------------------------------
    
    reclassification_group = case_when(
      
      natural_teeth >= 20 ~
        "Natural teeth ≥20",
      
      natural_teeth < 20 &
        functional_dentition >= 20 ~
        "Natural <20, functional ≥20",
      
      natural_teeth < 20 &
        functional_dentition < 20 ~
        "Natural <20, functional <20",
      
      TRUE ~ NA_character_
    ),
    
    reclassification_group = factor(
      reclassification_group,
      levels = c(
        "Natural teeth ≥20",
        "Natural <20, functional ≥20",
        "Natural <20, functional <20"
      )
    )
  )


# ============================================================
# 4. 그룹별 표본 수 확인
# ============================================================

cat("\nPanel A: Natural teeth\n")
print(
  table(
    dat$natural_teeth_group,
    useNA = "ifany"
  )
)

cat("\nPanel B: Functional dentition\n")
print(
  table(
    dat$functional_dentition_group,
    useNA = "ifany"
  )
)

cat("\nPanel C: Severe natural tooth loss\n")
print(
  table(
    dat$severe_tooth_loss,
    useNA = "ifany"
  )
)

cat("\nPanel D: Functional reclassification\n")
print(
  table(
    dat$reclassification_group,
    useNA = "ifany"
  )
)


# ============================================================
# 5. survfit 결과를 ggplot용 데이터로 변환하는 함수
# ============================================================

survfit_to_df <- function(fit) {
  
  s <- summary(fit)
  
  km_df <- data.frame(
    time = s$time,
    surv = s$surv,
    lower = s$lower,
    upper = s$upper,
    strata = as.character(s$strata)
  )
  
  # 각 그룹에 time = 0, survival = 1 추가
  strata_names <- names(fit$strata)
  
  baseline_df <- data.frame(
    time = 0,
    surv = 1,
    lower = 1,
    upper = 1,
    strata = strata_names
  )
  
  km_df <- bind_rows(
    baseline_df,
    km_df
  )
  
  # "변수명=그룹명"에서 그룹명만 남김
  km_df$strata <- sub(
    pattern = "^[^=]+=",
    replacement = "",
    x = km_df$strata
  )
  
  km_df
}


# ============================================================
# 6. Log-rank P value 함수
# ============================================================

get_logrank_p <- function(
    formula,
    data
) {
  
  test <- survdiff(
    formula,
    data = data
  )
  
  degrees_freedom <- length(test$n) - 1
  
  p_value <- pchisq(
    test$chisq,
    df = degrees_freedom,
    lower.tail = FALSE
  )
  
  if (p_value < 0.001) {
    return("Log-rank P < .001")
  }
  
  paste0(
    "Log-rank P = ",
    sub(
      "^0",
      "",
      sprintf("%.3f", p_value)
    )
  )
}


# ============================================================
# 7. 공통 Kaplan-Meier 그래프 함수
# ============================================================

make_km_plot <- function(
    data,
    group_var,
    title_text,
    color_values,
    legend_position = c(0.23, 0.18)
) {
  
  formula_km <- as.formula(
    paste0(
      "Surv(followup_years, death_event) ~ ",
      group_var
    )
  )
  
  fit <- survfit(
    formula = formula_km,
    data = data,
    conf.type = "log-log"
  )
  
  km_df <- survfit_to_df(fit)
  
  p_text <- get_logrank_p(
    formula = formula_km,
    data = data
  )
  
  ggplot(
    km_df,
    aes(
      x = time,
      y = surv,
      color = strata,
      fill = strata,
      group = strata
    )
  ) +
    
    # --------------------------------------------------------
  # 95% confidence interval
  # --------------------------------------------------------
  
  geom_ribbon(
    aes(
      ymin = lower,
      ymax = upper
    ),
    alpha = 0.20,
    color = NA,
    show.legend = FALSE
  ) +
    
    # --------------------------------------------------------
  # Kaplan-Meier step curves
  # --------------------------------------------------------
  
  geom_step(
    linewidth = 0.9,
    direction = "hv"
  ) +
    
    # --------------------------------------------------------
  # Log-rank P value
  # --------------------------------------------------------
  
  annotate(
    geom = "text",
    x = 8.4,
    y = 0.72,
    label = p_text,
    hjust = 1,
    size = 3.3,
    color = "black"
  ) +
    
    scale_color_manual(
      values = color_values
    ) +
    
    scale_fill_manual(
      values = color_values
    ) +
    
    scale_x_continuous(
      limits = c(0, 9),
      breaks = 0:9,
      expand = expansion(
        mult = c(0, 0)
      )
    ) +
    
    scale_y_continuous(
      limits = c(0.70, 1.005),
      breaks = seq(
        0.70,
        1.00,
        by = 0.05
      ),
      labels = function(x) {
        sprintf("%.2f", x)
      },
      expand = expansion(
        mult = c(0, 0.01)
      )
    ) +
    
    labs(
      title = title_text,
      x = "Follow-up time, years",
      y = "Survival probability",
      color = NULL,
      fill = NULL
    ) +
    
    theme_bw(
      base_size = 11
    ) +
    
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 12
      ),
      
      axis.title = element_text(
        size = 10
      ),
      
      axis.text = element_text(
        size = 9,
        color = "black"
      ),
      
      axis.ticks = element_line(
        color = "black",
        linewidth = 0.5
      ),
      
      legend.position = legend_position,
      
      legend.background = element_rect(
        fill = scales::alpha(
          "white",
          0.78
        ),
        color = NA
      ),
      
      legend.text = element_text(
        size = 8
      ),
      
      legend.key.width = grid::unit(
        0.65,
        "cm"
      ),
      
      legend.key.height = grid::unit(
        0.38,
        "cm"
      ),
      
      panel.grid.major = element_line(
        color = "grey88",
        linewidth = 0.3
      ),
      
      panel.grid.minor = element_blank(),
      
      panel.border = element_rect(
        color = "black",
        linewidth = 0.7
      ),
      
      plot.margin = margin(
        t = 9,
        r = 9,
        b = 9,
        l = 9
      )
    )
}


# ============================================================
# 8. Panel A
# ============================================================

plot_A <- make_km_plot(
  data = dat,
  group_var = "natural_teeth_group",
  title_text = "Survival by natural teeth group",
  
  color_values = c(
    "≥20"  = "#0072B2",
    "10–19" = "#009E73",
    "1–9"   = "#E69F00",
    "0"     = "#6A3D9A"
  ),
  
  legend_position = c(
    0.20,
    0.18
  )
)


# ============================================================
# 9. Panel B
# ============================================================

plot_B <- make_km_plot(
  data = dat,
  group_var = "functional_dentition_group",
  title_text = "Survival by functional dentition group",
  
  color_values = c(
    "≥20"  = "#0072B2",
    "10–19" = "#009E73",
    "1–9"   = "#E69F00",
    "0"     = "#6A3D9A"
  ),
  
  legend_position = c(
    0.22,
    0.18
  )
)


# ============================================================
# 10. Panel C
# ============================================================

plot_C <- make_km_plot(
  data = dat,
  group_var = "severe_tooth_loss",
  title_text = "Survival by severe natural tooth loss",
  
  color_values = c(
    "≥10 natural teeth" = "#0072B2",
    "<10 natural teeth" = "#D55E5E"
  ),
  
  legend_position = c(
    0.23,
    0.16
  )
)


# ============================================================
# 11. Panel D
# ============================================================

plot_D <- make_km_plot(
  data = dat,
  group_var = "reclassification_group",
  title_text = "Survival by functional reclassification",
  
  color_values = c(
    "Natural teeth ≥20" = "#0072B2",
    "Natural <20, functional ≥20" = "#2E8B57",
    "Natural <20, functional <20" = "#D55E5E"
  ),
  
  legend_position = c(
    0.30,
    0.18
  )
)


# ============================================================
# 12. 2 × 2 Figure 구성
# ============================================================

figure_2 <- (
  plot_A | plot_B
) / (
  plot_C | plot_D
) +
  
  plot_annotation(
    tag_levels = "A"
  ) &
  
  theme(
    plot.tag = element_text(
      face = "bold",
      size = 20,
      color = "black"
    ),
    
    plot.tag.position = c(
      0.01,
      0.99
    )
  )


# Figure 화면 출력
print(figure_2)


# ============================================================
# 13. Figure 2 저장
# ============================================================

# ------------------------------------------------------------
# PNG: 고해상도 600 dpi
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_folder,
    "Figure2_Kaplan_Meier_survival_updated.png"
  ),
  plot = figure_2,
  width = 13,
  height = 9.5,
  units = "in",
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)


# ------------------------------------------------------------
# TIFF: 고해상도 600 dpi, LZW 압축
# 저널 제출용
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_folder,
    "Figure2_Kaplan_Meier_survival_updated.tiff"
  ),
  plot = figure_2,
  width = 13,
  height = 9.5,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white",
  limitsize = FALSE
)


# ------------------------------------------------------------
# PDF: 벡터 파일
# 확대해도 글자와 선이 깨지지 않음
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_folder,
    "Figure2_Kaplan_Meier_survival_updated.pdf"
  ),
  plot = figure_2,
  width = 13,
  height = 9.5,
  units = "in",
  bg = "white",
  limitsize = FALSE
)


cat(
  "\nFigure 2 생성 및 저장 완료\n",
  "저장 폴더:\n",
  output_folder,
  "\n\n",
  "생성 파일:\n",
  "Figure2_Kaplan_Meier_survival_updated.png\n",
  "Figure2_Kaplan_Meier_survival_updated.tiff\n",
  "Figure2_Kaplan_Meier_survival_updated.pdf\n"
)
