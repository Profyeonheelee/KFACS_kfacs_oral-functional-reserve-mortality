"""
etable2_code.py -- Variable definitions and derivation rules for eTable 2.

This table is a prespecified data dictionary and does not require the analytic
datasets to be loaded. Running the script writes both CSV and Excel versions
of eTable 2.
"""

from pathlib import Path
import pandas as pd


COLUMNS = [
    "Domain",
    "Final variable",
    "Concept",
    "Source/raw variable(s)",
    "Source/time point",
    "Definition/derivation rule",
    "Primary role in analysis",
]


ROWS = [
    {
        "Domain": "Oral exposure",
        "Final variable": "natural_teeth",
        "Concept": "Remaining natural teeth",
        "Source/raw variable(s)": "Baseline dental/panoramic tooth-count variable",
        "Source/time point": "Baseline dental/panoramic assessment, 2016–2017",
        "Definition/derivation rule": (
            "Continuous count of remaining natural teeth; values outside the "
            "biologically plausible range of 0–32 were treated as missing."
        ),
        "Primary role in analysis": "Primary oral exposure; Table 1, Table 2, survival curves",
    },
    {
        "Domain": "Oral exposure",
        "Final variable": "dental_implant_count",
        "Concept": "Implant-supported restorations",
        "Source/raw variable(s)": "Baseline dental/panoramic implant-supported restoration count",
        "Source/time point": "Baseline dental/panoramic assessment, 2016–2017",
        "Definition/derivation rule": (
            "Continuous count of implant-supported restorations; implausible values "
            "were treated as missing."
        ),
        "Primary role in analysis": "Restored oral function marker; Table 1, Table 2",
    },
    {
        "Domain": "Oral exposure",
        "Final variable": "any_dental_implant",
        "Concept": "Any implant-supported restoration",
        "Source/raw variable(s)": "dental_implant_count",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Binary indicator: 1 if implant-supported restoration count ≥1; 0 if count = 0."
        ),
        "Primary role in analysis": "Secondary implant exposure; Table 1, Table 2",
    },
    {
        "Domain": "Oral exposure",
        "Final variable": "functional_dentition",
        "Concept": "Functional dentition / oral functional reserve",
        "Source/raw variable(s)": "natural_teeth + dental_implant_count",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Sum of remaining natural teeth and implant-supported restorations."
        ),
        "Primary role in analysis": (
            "Primary constructed oral functional reserve variable; Table 1, Table 2"
        ),
    },
    {
        "Domain": "Oral exposure",
        "Final variable": "functional_dentition_lt20",
        "Concept": "Reduced functional dentition",
        "Source/raw variable(s)": "functional_dentition",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Binary indicator: 1 if functional dentition <20; 0 if functional dentition ≥20."
        ),
        "Primary role in analysis": "Baseline Table 1 grouping variable",
    },
    {
        "Domain": "Oral exposure",
        "Final variable": "natural_teeth_lt10",
        "Concept": "Severe natural tooth loss",
        "Source/raw variable(s)": "natural_teeth",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Binary indicator: 1 if natural teeth <10; 0 if natural teeth ≥10."
        ),
        "Primary role in analysis": "Clinically interpretable oral exposure; Table 2",
    },
    {
        "Domain": "Oral exposure",
        "Final variable": "functional_dentition_lt10",
        "Concept": "Low functional dentition",
        "Source/raw variable(s)": "functional_dentition",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Binary indicator: 1 if functional dentition <10; 0 if functional dentition ≥10."
        ),
        "Primary role in analysis": "Main pathway exposure; Table 2, Table 3",
    },
    {
        "Domain": "Mortality outcome",
        "Final variable": "death_event",
        "Concept": "All-cause mortality during follow-up",
        "Source/raw variable(s)": (
            "2016–2024: center_pal_5, home_pal_5, call_pal_5; "
            "2025 update: center_pal_5_1, home_pal_5_1, call_pal_5_1"
        ),
        "Source/time point": "Serial KFACS follow-up records linked by participant ID",
        "Definition/derivation rule": (
            "Binary indicator: 1 if death was recorded at least once in center, "
            "home-visit, or telephone follow-up proxy-response current-status "
            "variables; otherwise 0."
        ),
        "Primary role in analysis": (
            "Primary outcome; Table 1, Table 2, Table 3, survival analyses"
        ),
    },
    {
        "Domain": "Mortality outcome",
        "Final variable": "death_date",
        "Concept": "Date of death for Cox analysis",
        "Source/raw variable(s)": (
            "2016–2024: center_pal_6, home_pal_6, call_pal_6; "
            "2025 update: center_pal_6_1, home_pal_6_1, call_pal_6_1"
        ),
        "Source/time point": "Serial KFACS follow-up records",
        "Definition/derivation rule": (
            "Death date was derived from year-month death-date variables. When death "
            "was reported in multiple sources, the earliest usable death date was "
            "retained. Incomplete or unknown death dates were flagged in "
            "death_date_quality and approximated according to predefined rules."
        ),
        "Primary role in analysis": "Cox time-to-event calculation",
    },
    {
        "Domain": "Mortality outcome",
        "Final variable": "death_source",
        "Concept": "Death ascertainment source",
        "Source/raw variable(s)": "Center, home-visit, or telephone death-status source",
        "Source/time point": "Serial KFACS follow-up records",
        "Definition/derivation rule": (
            "Source corresponding to the earliest usable death record."
        ),
        "Primary role in analysis": "Death ascertainment audit",
    },
    {
        "Domain": "Mortality outcome",
        "Final variable": "followup_years",
        "Concept": "Follow-up time",
        "Source/raw variable(s)": "Baseline dental assessment date, death_date, censor date",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Time from baseline dental/panoramic assessment to death or censoring. "
            "Survivors were censored at the last available mortality follow-up period, "
            "including the 2025 mortality update when available."
        ),
        "Primary role in analysis": "Cox model time variable",
    },
    {
        "Domain": "Core covariate",
        "Final variable": "age_baseline",
        "Concept": "Baseline age",
        "Source/raw variable(s)": "KFACS baseline age variable",
        "Source/time point": "Baseline examination",
        "Definition/derivation rule": "Continuous age in years.",
        "Primary role in analysis": "Model 1 and Model 2 covariate; Table 1",
    },
    {
        "Domain": "Core covariate",
        "Final variable": "female",
        "Concept": "Female sex",
        "Source/raw variable(s)": "KFACS sex variable",
        "Source/time point": "Baseline examination",
        "Definition/derivation rule": "Binary indicator for female sex.",
        "Primary role in analysis": "Model 1 and Model 2 covariate; Table 1",
    },
    {
        "Domain": "Core covariate",
        "Final variable": "center",
        "Concept": "Study center",
        "Source/raw variable(s)": "KFACS center/location variable",
        "Source/time point": "Baseline examination",
        "Definition/derivation rule": "Categorical study-center variable.",
        "Primary role in analysis": "Model 1 and Model 2 covariate",
    },
    {
        "Domain": "Core covariate",
        "Final variable": "cohort_start_year",
        "Concept": "Cohort start year",
        "Source/raw variable(s)": "Participant ID prefix / baseline year",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "2016 for kf16 participants and 2017 for kf17 participants."
        ),
        "Primary role in analysis": "Model 1 and Model 2 covariate",
    },
    {
        "Domain": "Socioeconomic covariate",
        "Final variable": "education_years",
        "Concept": "Education",
        "Source/raw variable(s)": "KFACS baseline education variable",
        "Source/time point": "Baseline questionnaire",
        "Definition/derivation rule": "Continuous years of education.",
        "Primary role in analysis": "Table 1; Model 2 covariate",
    },
    {
        "Domain": "Socioeconomic covariate",
        "Final variable": "income_level",
        "Concept": "Income level",
        "Source/raw variable(s)": "KFACS baseline income variable",
        "Source/time point": "Baseline questionnaire",
        "Definition/derivation rule": (
            "Categorical variable with 3 levels: high, middle, and low/no income. "
            "Missing values were multiply imputed as part of the multiple-imputation "
            "procedure (m = 20 datasets) and were not retained as a separate analytic category."
        ),
        "Primary role in analysis": "Table 1; Model 2 covariate",
    },
    {
        "Domain": "Behavioral covariate",
        "Final variable": "smoking_history_ge100",
        "Concept": "Smoking history",
        "Source/raw variable(s)": "KFACS baseline smoking variable",
        "Source/time point": "Baseline questionnaire",
        "Definition/derivation rule": (
            "Binary indicator for lifetime smoking history ≥100 cigarettes."
        ),
        "Primary role in analysis": "Table 1; Model 2 covariate",
    },
    {
        "Domain": "Behavioral covariate",
        "Final variable": "ever_alcohol_use",
        "Concept": "Ever alcohol use",
        "Source/raw variable(s)": "KFACS baseline alcohol-use variable",
        "Source/time point": "Baseline questionnaire",
        "Definition/derivation rule": "Binary indicator for ever alcohol use.",
        "Primary role in analysis": "Table 1; Model 2 covariate",
    },
    {
        "Domain": "Systemic disease covariate",
        "Final variable": "systemic_disease_count",
        "Concept": "Baseline systemic disease burden",
        "Source/raw variable(s)": "Baseline physician-reported systemic disease category indicators",
        "Source/time point": "Baseline health-status questionnaire",
        "Definition/derivation rule": (
            "Count of baseline physician-reported systemic disease categories coded "
            "as present. Baseline dementia was not included in the main Table 1 or "
            "main adjustment models."
        ),
        "Primary role in analysis": (
            "Table 1; Table 2 Model 2 covariate; Table 3 final systemic disease/reserve domain"
        ),
    },
    {
        "Domain": "Weight loss/nutrition",
        "Final variable": "weight_loss_chs",
        "Concept": "CHS weight-loss component",
        "Source/raw variable(s)": "KFACS CHS weight-loss item",
        "Source/time point": "Baseline frailty questionnaire",
        "Definition/derivation rule": (
            "Binary CHS weight-loss component according to the KFACS coding rule."
        ),
        "Primary role in analysis": "Table 1; Table 3 weight loss/nutrition domain",
    },
    {
        "Domain": "Weight loss/nutrition",
        "Final variable": "bmi",
        "Concept": "Body mass index",
        "Source/raw variable(s)": "Baseline height and weight / BMI variable",
        "Source/time point": "Baseline examination",
        "Definition/derivation rule": "Continuous BMI in kg/m².",
        "Primary role in analysis": "Table 1; Table 3 weight loss/nutrition domain",
    },
    {
        "Domain": "Weight loss/nutrition",
        "Final variable": "mna_screen_risk",
        "Concept": "MNA screening risk",
        "Source/raw variable(s)": "MNA screening score/category",
        "Source/time point": "Baseline nutrition assessment",
        "Definition/derivation rule": (
            "Binary indicator for nutritional risk based on MNA screening classification."
        ),
        "Primary role in analysis": "Table 1; Table 3 weight loss/nutrition domain",
    },
    {
        "Domain": "Weight loss/nutrition / systemic reserve descriptor",
        "Final variable": "albumin_g_dl",
        "Concept": "Serum albumin",
        "Source/raw variable(s)": "Baseline serum albumin variable",
        "Source/time point": "Baseline laboratory examination",
        "Definition/derivation rule": "Continuous albumin in g/dL.",
        "Primary role in analysis": (
            "Table 1 descriptor; component of adverse systemic reserve marker count"
        ),
    },
    {
        "Domain": "CHS frailty",
        "Final variable": "chs_frailty_status_kfacs",
        "Concept": "KFACS-based CHS frailty status",
        "Source/raw variable(s)": (
            "KFACS CHS components: weight loss, exhaustion, weakness, slowness, "
            "low physical activity"
        ),
        "Source/time point": (
            "Derived from baseline frailty questionnaire and performance measures"
        ),
        "Definition/derivation rule": (
            "Five CHS components were scored according to KFACS coding rules and summed. "
            "Frailty status was classified as robust for score 0, prefrail for score 1–2, "
            "and frail for score 3–5."
        ),
        "Primary role in analysis": "Main Table 1 frailty variable",
    },
    {
        "Domain": "CHS frailty component",
        "Final variable": "chs_weakness_kfacs",
        "Concept": "Weakness",
        "Source/raw variable(s)": "Grip strength, sex, BMI",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Weakness was defined using KFACS sex- and BMI-specific grip-strength thresholds."
        ),
        "Primary role in analysis": "Component of chs_frailty_status_kfacs",
    },
    {
        "Domain": "CHS frailty component",
        "Final variable": "chs_slowness_kfacs",
        "Concept": "Slowness",
        "Source/raw variable(s)": "4-m gait speed, sex, height",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Slowness was defined using KFACS sex- and height-specific 4-m gait-speed thresholds."
        ),
        "Primary role in analysis": "Component of chs_frailty_status_kfacs",
    },
    {
        "Domain": "CHS frailty component",
        "Final variable": "chs_low_activity_kfacs",
        "Concept": "Low physical activity",
        "Source/raw variable(s)": "IPAQ-based physical activity variables",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Low physical activity was defined as <494.65 kcal/week in men and "
            "<283.50 kcal/week in women according to KFACS coding rules."
        ),
        "Primary role in analysis": "Component of chs_frailty_status_kfacs",
    },
    {
        "Domain": "Physical frailty pathway",
        "Final variable": "physical_frailty_burden_no_weight_loss_kfacs",
        "Concept": "Physical frailty burden excluding weight loss",
        "Source/raw variable(s)": (
            "KFACS CHS components excluding weight loss: exhaustion, weakness, "
            "slowness, low physical activity"
        ),
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Sum of non-weight-loss CHS components. Weight loss was excluded to avoid "
            "overlap with the preceding weight loss/nutrition domain in the pathway model."
        ),
        "Primary role in analysis": "Table 3 physical frailty burden domain",
    },
    {
        "Domain": "Physical function descriptor",
        "Final variable": "grip_strength_max",
        "Concept": "Grip strength",
        "Source/raw variable(s)": "Right/left grip-strength measurements",
        "Source/time point": "Baseline physical performance assessment",
        "Definition/derivation rule": (
            "Maximum valid grip strength across available trials, kg."
        ),
        "Primary role in analysis": "Table 1 descriptor",
    },
    {
        "Domain": "Physical function descriptor",
        "Final variable": "gait_speed_m_per_s",
        "Concept": "Gait speed",
        "Source/raw variable(s)": "4-m walking-time variables",
        "Source/time point": "Baseline physical performance assessment",
        "Definition/derivation rule": (
            "Gait speed calculated in m/s from valid walking-time measures."
        ),
        "Primary role in analysis": "Table 1 descriptor",
    },
    {
        "Domain": "Cognition",
        "Final variable": "mmse_kc_score",
        "Concept": "Cognitive function",
        "Source/raw variable(s)": "MMSE-KC score",
        "Source/time point": "Baseline cognitive assessment",
        "Definition/derivation rule": "Continuous MMSE-KC score.",
        "Primary role in analysis": "Table 1; Table 3 cognition domain",
    },
    {
        "Domain": "Systemic reserve",
        "Final variable": "systemic_reserve_adverse_count",
        "Concept": "Adverse systemic reserve marker count",
        "Source/raw variable(s)": "Albumin, hemoglobin, hs-CRP, vitamin D, eGFR",
        "Source/time point": "Baseline laboratory examination",
        "Definition/derivation rule": (
            "Count of four laboratory-based adverse indicators: low albumin <3.8 g/dL, "
            "anemia/low hemoglobin <12.0 g/dL in women and <13.0 g/dL in men, elevated "
            "hs-CRP >3.0 mg/L, and vitamin D deficiency <20 ng/mL. The count was calculated "
            "when at least 3 of the 4 components were non-missing. Reduced eGFR <60 "
            "mL/min/1.73 m² is shown for reference but was not included in the count because "
            "creatinine was measured at only three of the nine centers. Missing laboratory "
            "values were not assumed to be normal and were not coded as 0."
        ),
        "Primary role in analysis": (
            "Table 1; Table 2 Model 2 covariate; Table 3 final systemic disease/reserve domain"
        ),
    },
    {
        "Domain": "Systemic reserve component",
        "Final variable": "low_albumin",
        "Concept": "Low albumin",
        "Source/raw variable(s)": "Serum albumin",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Binary indicator: 1 if albumin <3.8 g/dL; 0 if ≥3.8 g/dL; "
            "missing if unavailable or invalid."
        ),
        "Primary role in analysis": "Component of adverse systemic reserve marker count",
    },
    {
        "Domain": "Systemic reserve component",
        "Final variable": "anemia_low_hemoglobin",
        "Concept": "Anemia / low hemoglobin",
        "Source/raw variable(s)": "Hemoglobin, sex",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Binary indicator: 1 if hemoglobin <12.0 g/dL in women or <13.0 g/dL "
            "in men; 0 otherwise; missing if hemoglobin or sex unavailable."
        ),
        "Primary role in analysis": "Component of adverse systemic reserve marker count",
    },
    {
        "Domain": "Systemic reserve component",
        "Final variable": "elevated_hscrp",
        "Concept": "Elevated systemic inflammation marker",
        "Source/raw variable(s)": "hs-CRP",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Binary indicator: 1 if hs-CRP >3.0 mg/L; 0 if ≤3.0 mg/L; "
            "missing if unavailable or invalid."
        ),
        "Primary role in analysis": "Component of adverse systemic reserve marker count",
    },
    {
        "Domain": "Systemic reserve component",
        "Final variable": "vitamin_d_deficiency",
        "Concept": "Vitamin D deficiency",
        "Source/raw variable(s)": "25-hydroxyvitamin D",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Binary indicator: 1 if vitamin D <20 ng/mL; 0 if ≥20 ng/mL; "
            "missing if unavailable or invalid."
        ),
        "Primary role in analysis": "Component of adverse systemic reserve marker count",
    },
    {
        "Domain": "Systemic reserve component",
        "Final variable": "reduced_egfr",
        "Concept": "Reduced kidney function",
        "Source/raw variable(s)": "eGFR",
        "Source/time point": "Derived variable",
        "Definition/derivation rule": (
            "Binary indicator: 1 if eGFR <60 mL/min/1.73 m²; 0 if ≥60 mL/min/1.73 m²; "
            "missing if unavailable or invalid."
        ),
        "Primary role in analysis": "Component of adverse systemic reserve marker count",
    },
]


def build_table():
    return pd.DataFrame(ROWS, columns=COLUMNS)


def main():
    output_dir = Path(".")
    table = build_table()

    csv_path = output_dir / "etable2_output.csv"
    xlsx_path = output_dir / "etable2_output.xlsx"

    table.to_csv(csv_path, index=False, encoding="utf-8-sig")

    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        table.to_excel(writer, sheet_name="eTable 2", index=False)
        ws = writer.book["eTable 2"]
        ws.freeze_panes = "A2"
        ws.auto_filter.ref = ws.dimensions

        widths = {
            "A": 32,
            "B": 44,
            "C": 38,
            "D": 62,
            "E": 48,
            "F": 95,
            "G": 58,
        }
        for col, width in widths.items():
            ws.column_dimensions[col].width = width

        from openpyxl.styles import Alignment, Font

        for cell in ws[1]:
            cell.font = Font(bold=True)
            cell.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)

        for row in ws.iter_rows(min_row=2):
            for cell in row:
                cell.alignment = Alignment(vertical="top", wrap_text=True)

    print(f"eTable 2 rows: {len(table)}")
    print(f"Wrote: {csv_path}")
    print(f"Wrote: {xlsx_path}")


if __name__ == "__main__":
    main()
