# 9545 Capstone: Dyadic Confirmatory Factor Analysis

# Part 1. Introduction and Research Questions
## 1.1 Background
## 1.2 Research Questions
## 1.3 Expectations and Analytic Approach
## 1.4 Dataset
## 1.5 Conceptual Dyadic Model

# Part 2. Data Preparation and Validation
## 2.1 Load Packages and Import Data
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(psych)
library(gt)
library(lavaan)
library(ggplot2)
library(semPlot)

options(dplyr.summarise.inform = FALSE)

piq_raw <- readr::read_csv(
  "piq_scored_all.csv",
  show_col_types = FALSE
)


## 2.2 Define Variables and Analysis Helpers
#| label: define-variables
qi_items_9 <- paste0("QI", 1:9)
qi_items_7 <- paste0("QI", 1:7)
id_columns <- c(
  "StudySessionID",
  "ConversationNum",
  "PID",
  "ptrID"
)

required_columns <- c(id_columns, qi_items_9)

missing_columns <- setdiff(
  required_columns,
  names(piq_raw)
)

if (length(missing_columns) > 0L) {
  stop(
    "Required columns are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

fit_indices <- c(
  "chisq", "df", "pvalue",
  "cfi", "tli", "rmsea",
  "rmsea.ci.lower", "rmsea.ci.upper",
  "srmr"
)

fit_row <- function(fit, label) {
  values <- lavaan::fitMeasures(fit, fit_indices)
  
  tibble::as_tibble_row(as.list(values)) |>
    mutate(
      Model = label,
      N = as.integer(lavaan::lavInspect(fit, "nobs")),
      .before = 1
    )
}

show_fit <- function(table) {
  table |>
    gt::gt() |>
    gt::cols_label(
      N = "Records",
      chisq = "Chi-square",
      df = "df",
      pvalue = "p",
      cfi = "CFI",
      tli = "TLI",
      rmsea = "RMSEA",
      rmsea.ci.lower = "90% CI lower",
      rmsea.ci.upper = "90% CI upper",
      srmr = "SRMR"
    ) |>
    gt::fmt_integer(columns = c(N, df)) |>
    gt::fmt_number(
      columns = chisq,
      decimals = 2
    ) |>
    gt::fmt_number(
      columns = c(
        pvalue, cfi, tli, rmsea,
        rmsea.ci.lower, rmsea.ci.upper, srmr
      ),
      decimals = 3
    ) |>
    gt::sub_missing(missing_text = "—")
}

partner_correlation <- function(fit) {
  unname(
    lavaan::lavInspect(fit, "cor.lv")[
      "quality_A", "quality_B"
    ]
  )
}

fmt3 <- function(x) {
  sprintf("%.3f", x)
}

fmt2 <- function(x) {
  sprintf("%.2f", x)
}

fmt_p <- function(x) {
  if (is.na(x)) {
    return("p unavailable")
  }
  
  if (x < .001) {
    return("p < .001")
  }
  
  paste0("p = ", sprintf("%.3f", x))
}

fit_text <- function(fit) {
  x <- lavaan::fitMeasures(fit, fit_indices)
  
  paste0(
    "χ²(", as.integer(x["df"]), ") = ",
    fmt2(x["chisq"]), ", ",
    fmt_p(x["pvalue"]),
    ", CFI = ", fmt3(x["cfi"]),
    ", TLI = ", fmt3(x["tli"]),
    ", RMSEA = ", fmt3(x["rmsea"]),
    ", 90% CI [",
    fmt3(x["rmsea.ci.lower"]), ", ",
    fmt3(x["rmsea.ci.upper"]), "]",
    ", SRMR = ", fmt3(x["srmr"])
  )
}

check_fit <- function(fit) {
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) {
    stop("A model did not converge. Inspect it before interpreting results.")
  }
  
  invisible(fit)
}


## 2.3 Select IDs and Item Responses
#| label: select-data
piq_dat <- piq_raw |>
  select(all_of(required_columns))

numeric_items <- vapply(
  piq_dat[qi_items_9],
  is.numeric,
  logical(1)
)

if (!all(numeric_items)) {
  stop("All QI items must be numeric before analysis.")
}

glimpse(piq_dat)


## 2.4 Check ID Completeness and Self-Ratings
#| label: check-ids
id_check <- piq_dat |>
  summarise(
    rating_records = n(),
    missing_session = sum(is.na(StudySessionID)),
    missing_round = sum(is.na(ConversationNum)),
    missing_PID = sum(is.na(PID)),
    missing_partner = sum(is.na(ptrID)),
    self_ratings = sum(PID == ptrID, na.rm = TRUE)
  )

id_check |>
  gt::gt() |>
  gt::tab_header(
    title = "ID and Rating-Record Checks"
  )

if (anyNA(piq_dat[id_columns])) {
  stop("Missing IDs must be resolved before matching partners.")
}

if (any(piq_dat$PID == piq_dat$ptrID)) {
  stop("Self-ratings were found. Resolve them before matching partners.")
}


## 2.5 Check Duplicate Directed Ratings
#| label: check-duplicates
#| include: false
duplicate_ratings <- piq_dat |>
  count(
    StudySessionID,
    ConversationNum,
    PID,
    ptrID,
    name = "n_records"
  ) |>
  filter(n_records > 1L)

tibble(
  Check = "Duplicate directed-rating keys",
  Count = nrow(duplicate_ratings)
) |>
  gt::gt()

if (nrow(duplicate_ratings) > 0L) {
  print(duplicate_ratings)
  
  stop(
    "Duplicate directed ratings were found. ",
    "Resolve them before reshaping."
  )
}


## 2.6 Check Item Ranges and Missing Responses
#| label: check-item-values
#| include: false
item_check <- piq_dat |>
  select(all_of(qi_items_9)) |>
  pivot_longer(
    cols = everything(),
    names_to = "Item",
    values_to = "Response"
  ) |>
  group_by(Item) |>
  summarise(
    Nonmissing = sum(!is.na(Response)),
    Missing = sum(is.na(Response)),
    Missing_percent = mean(is.na(Response)) * 100,
    Outside_1_to_7 = sum(
      !is.na(Response) &
        (Response < 1 | Response > 7)
    ),
    Noninteger = sum(
      !is.na(Response) &
        abs(Response - round(Response)) > 1e-8
    ),
    .groups = "drop"
  )

item_check |>
  gt::gt() |>
  gt::tab_header(
    title = "Item Range and Missingness Checks"
  ) |>
  gt::fmt_number(
    columns = Missing_percent,
    decimals = 2
  )

if (
  any(item_check$Outside_1_to_7 > 0L) ||
  any(item_check$Noninteger > 0L)
) {
  stop("Unexpected QI response values were found.")
}


## 2.7 Summarize Item Responses
#| label: item-descriptives

item_descriptives <- piq_dat |>
  select(all_of(qi_items_9)) |>
  psych::describe() |>
  as.data.frame() |>
  tibble::rownames_to_column("Item") |>
  select(
    Item, n, mean, sd, median,
    min, max, skew, kurtosis
  )

item_descriptives |>
  gt::gt() |>
  gt::tab_header(
    title = "Item Descriptive Statistics",
    subtitle = "Perceived interaction quality"
  ) |>
  gt::cols_label(
    n = "N",
    mean = "M",
    sd = "SD",
    median = "Median",
    min = "Min",
    max = "Max",
    skew = "Skewness",
    kurtosis = "Kurtosis"
  ) |>
  gt::fmt_number(
    columns = c(mean, sd, median, skew, kurtosis),
    decimals = 2
  ) |>
  gt::tab_source_note(
    source_note = paste(
      "N represents nonmissing rating records,",
      "not independent participants."
    )
  )


# Part 3. Preliminary One-Factor CFA
## 3.1 Purpose and Estimation

## 3.2 Test 1: Initial Nine-Item Model
#| label: fit-individual-nine

qi_model1 <- "
  perceived_interaction_quality =~
    QI1 + QI2 + QI3 + QI4 + QI5 +
    QI6 + QI7 + QI8 + QI9
"

qi_fit1 <- lavaan::cfa(
  model = qi_model1,
  data = piq_dat,
  estimator = "ML",
  missing = "ML",
  std.lv = TRUE
)

check_fit(qi_fit1)

summary(
  qi_fit1,
  standardized = TRUE,
  fit.measures = TRUE,
  rsquare = TRUE
)


## 3.3 Inspect Potential Local Misfit
#| label: individual-nine-mi

lavaan::modindices(qi_fit1) |>
  filter(op == "~~", lhs != rhs) |>
  arrange(desc(mi)) |>
  select(lhs, op, rhs, mi, epc, sepc.all) |>
  slice_head(n = 10) |>
  gt::gt() |>
  gt::tab_header(
    title = "Largest Residual-Covariance Modification Indices",
    subtitle = "Initial nine-item CFA"
  ) |>
  gt::fmt_number(
    columns = c(mi, epc, sepc.all),
    decimals = 3
  )


## 3.4 Test 2: Nine-Item Model With QI8–QI9 Residual Covariance
#| label: fit-individual-nine-residual

qi_model2 <- "
  perceived_interaction_quality =~
    QI1 + QI2 + QI3 + QI4 + QI5 +
    QI6 + QI7 + QI8 + QI9

  QI8 ~~ QI9
"

qi_fit2 <- lavaan::cfa(
  model = qi_model2,
  data = piq_dat,
  estimator = "ML",
  missing = "ML",
  std.lv = TRUE
)

check_fit(qi_fit2)

summary(
  qi_fit2,
  standardized = TRUE,
  fit.measures = TRUE,
  rsquare = TRUE
)


## 3.5 Test 3: Seven-Item Sensitivity Model
#| label: fit-individual-seven

qi_model3 <- "
  perceived_interaction_quality =~
    QI1 + QI2 + QI3 + QI4 + QI5 +
    QI6 + QI7
"

qi_fit3 <- lavaan::cfa(
  model = qi_model3,
  data = piq_dat,
  estimator = "ML",
  missing = "ML",
  std.lv = TRUE
)

check_fit(qi_fit3)

summary(
  qi_fit3,
  standardized = TRUE,
  fit.measures = TRUE,
  rsquare = TRUE
)


## 3.6 Summarize Preliminary Model Fit
#| label: individual-fit-table

individual_fit_table <- bind_rows(
  fit_row(
    qi_fit1,
    "Test 1: Nine-item one-factor CFA"
  ),
  fit_row(
    qi_fit2,
    "Test 2: Nine-item CFA with QI8–QI9 residual covariance"
  ),
  fit_row(
    qi_fit3,
    "Test 3: Seven-item one-factor CFA"
  )
)

individual_fit_table |>
  show_fit() |>
  gt::tab_header(
    title = "Preliminary One-Factor CFA Results"
  ) |>
  gt::tab_source_note(
    source_note = paste(
      "Seven- and nine-item models contain different indicators.",
      "Their fit indices are presented descriptively;",
      "they are not treated as a conventional nested-model comparison.",
      "Rows represent individual rating records."
    )
  )

## 3.7 Inspect Standardized Factor Loadings
#| label: individual-loadings

individual_models <- list(
  "Test 1: Nine items" = qi_fit1,
  "Test 2: Nine items + residual covariance" = qi_fit2,
  "Test 3: Seven items" = qi_fit3
)

individual_loadings <- bind_rows(
  lapply(names(individual_models), function(label) {
    lavaan::standardizedSolution(
      individual_models[[label]]
    ) |>
      filter(op == "=~") |>
      transmute(
        Model = label,
        Item = rhs,
        Standardized_loading = est.std
      )
  })
)

individual_loadings |>
  gt::gt(groupname_col = "Model") |>
  gt::tab_header(
    title = "Standardized Factor Loadings"
  ) |>
  gt::fmt_number(
    columns = Standardized_loading,
    decimals = 3
  )


# Part 4. Construct the Dyadic Dataset
## 4.1 Assign Partner Labels
#| label: assign-partner-labels

piq_pair_long <- piq_dat |>
  mutate(
    A_ID = pmin(PID, ptrID),
    B_ID = pmax(PID, ptrID),
    member = if_else(PID == A_ID, "A", "B")
  )


## 4.2 Check Reciprocal Pair Completeness
#| label: check-reciprocal-pairs

pair_check <- piq_pair_long |>
  group_by(
    StudySessionID,
    ConversationNum,
    A_ID,
    B_ID
  ) |>
  summarise(
    n_records = n(),
    n_A = sum(member == "A"),
    n_B = sum(member == "B"),
    .groups = "drop"
  )

pair_check |>
  count(n_records, n_A, n_B, name = "Conversations") |>
  gt::gt() |>
  gt::tab_header(
    title = "Reciprocal Pair Completeness"
  )

incomplete_pairs <- pair_check |>
  filter(
    n_records != 2L |
      n_A != 1L |
      n_B != 1L
  )

if (nrow(incomplete_pairs) > 0L) {
  print(incomplete_pairs)

  stop(
    "Some conversations do not have exactly one rating per partner. ",
    "Review these records before continuing."
  )
}

## 4.3 Reshape to One Row per Conversation
#| label: reshape-dyads

piq_dyad <- piq_pair_long |>
  select(
    StudySessionID,
    ConversationNum,
    A_ID,
    B_ID,
    member,
    all_of(qi_items_9)
  ) |>
  pivot_wider(
    id_cols = c(
      StudySessionID,
      ConversationNum,
      A_ID,
      B_ID
    ),
    names_from = member,
    values_from = all_of(qi_items_9),
    names_glue = "{.value}_{member}"
  ) |>
  arrange(
    StudySessionID,
    ConversationNum,
    A_ID,
    B_ID
  )

stopifnot(
  nrow(piq_dyad) == nrow(pair_check),
  nrow(piq_dat) == 2L * nrow(piq_dyad)
)

glimpse(piq_dyad)


## 4.4 Summarize the Analytic Sample
#| label: sample-summary

n_rating_records <- nrow(piq_dat)

n_participants <- piq_dat |>
  distinct(StudySessionID, PID) |>
  nrow()

n_sessions <- n_distinct(piq_dat$StudySessionID)

n_conversations <- nrow(piq_dyad)

sample_summary <- tibble(
  Unit = c(
    "Individual rating records",
    "Participants (session–PID combinations)",
    "Sessions",
    "Reciprocally matched conversations"
  ),
  N = c(
    n_rating_records,
    n_participants,
    n_sessions,
    n_conversations
  )
)

sample_summary |>
  gt::gt() |>
  gt::tab_header(
    title = "Analytic Sample"
  ) |>
  gt::fmt_integer(columns = N)


# Part 5. Fit the Dyadic CFA Models
## 5.1 Dyadic Measurement Structure
## 5.2 Model D1: Initial Nine-Item Dyadic CFA
#| label: fit-dyadic-nine

qi_dyad_model1 <- "
  quality_A =~
    QI1_A + QI2_A + QI3_A + QI4_A + QI5_A +
    QI6_A + QI7_A + QI8_A + QI9_A

  quality_B =~
    QI1_B + QI2_B + QI3_B + QI4_B + QI5_B +
    QI6_B + QI7_B + QI8_B + QI9_B

  quality_A ~~ quality_B
"

qi_dyad_fit1 <- lavaan::cfa(
  model = qi_dyad_model1,
  data = piq_dyad,
  estimator = "ML",
  missing = "ML",
  std.lv = TRUE
)

check_fit(qi_dyad_fit1)

summary(
  qi_dyad_fit1,
  standardized = TRUE,
  fit.measures = TRUE,
  rsquare = TRUE
)


## 5.3 Inspect QI8–QI9 Correlation Residuals
#| label: dyadic-correlation-residuals

dyadic_residuals <- lavaan::lavResiduals(
  qi_dyad_fit1,
  type = "cor"
)

qi89_correlation_residuals <- tibble(
  Partner = c("A", "B"),
  Item_pair = c(
    "QI8_A and QI9_A",
    "QI8_B and QI9_B"
  ),
  Correlation_residual = c(
    dyadic_residuals$cov["QI8_A", "QI9_A"],
    dyadic_residuals$cov["QI8_B", "QI9_B"]
  )
)

qi89_correlation_residuals |>
  gt::gt() |>
  gt::tab_header(
    title = "QI8–QI9 Correlation Residuals",
    subtitle = "Initial nine-item dyadic CFA"
  ) |>
  gt::fmt_number(
    columns = Correlation_residual,
    decimals = 3
  ) |>
  gt::tab_source_note(
    source_note = paste(
      "Observed correlation minus model-implied correlation.",
      "These values are not estimated error correlations."
    )
  )


## 5.4 Model D2: Nine-Item Dyadic CFA With Residual Covariances
#| label: fit-dyadic-nine-residual

qi_dyad_model2 <- "
  quality_A =~
    QI1_A + QI2_A + QI3_A + QI4_A + QI5_A +
    QI6_A + QI7_A + QI8_A + QI9_A

  quality_B =~
    QI1_B + QI2_B + QI3_B + QI4_B + QI5_B +
    QI6_B + QI7_B + QI8_B + QI9_B

  quality_A ~~ quality_B

  QI8_A ~~ QI9_A
  QI8_B ~~ QI9_B
"

qi_dyad_fit2 <- lavaan::cfa(
  model = qi_dyad_model2,
  data = piq_dyad,
  estimator = "ML",
  missing = "ML",
  std.lv = TRUE
)

check_fit(qi_dyad_fit2)

summary(
  qi_dyad_fit2,
  standardized = TRUE,
  fit.measures = TRUE,
  rsquare = TRUE
)


## 5.5 Model D3: Seven-Item Dyadic CFA
#| label: fit-dyadic-seven

qi_dyad_model_7 <- "
  quality_A =~
    QI1_A + QI2_A + QI3_A + QI4_A + QI5_A +
    QI6_A + QI7_A

  quality_B =~
    QI1_B + QI2_B + QI3_B + QI4_B + QI5_B +
    QI6_B + QI7_B

  quality_A ~~ quality_B
"

qi_dyad_fit_7 <- lavaan::cfa(
  model = qi_dyad_model_7,
  data = piq_dyad,
  estimator = "ML",
  missing = "ML",
  std.lv = TRUE
)

check_fit(qi_dyad_fit_7)

summary(
  qi_dyad_fit_7,
  standardized = TRUE,
  fit.measures = TRUE,
  rsquare = TRUE
)


## 5.6 Compare Dyadic Model Fit
#| label: dyadic-fit-table

dyadic_fit_table <- bind_rows(
  fit_row(
    qi_dyad_fit1,
    "D1: Nine-item dyadic CFA"
  ),
  fit_row(
    qi_dyad_fit2,
    "D2: Nine-item dyadic CFA with residual covariances"
  ),
  fit_row(
    qi_dyad_fit_7,
    "D3: Seven-item dyadic CFA"
  )
)

dyadic_fit_table |>
  show_fit() |>
  gt::tab_header(
    title = "Dyadic CFA Model Fit"
  ) |>
  gt::tab_source_note(
    source_note = paste(
      "Records represent conversations.",
      "Seven- and nine-item specifications use different indicators;",
      "comparisons across these item sets are descriptive.",
      "Dependence across conversations has not been modeled."
    )
  )


## 5.7 Compare Partner Latent Correlations
#| label: partner-correlations

partner_r_table <- tibble(
  Model = c(
    "D1: Nine-item dyadic CFA",
    "D2: Nine-item dyadic CFA with residual covariances",
    "D3: Seven-item dyadic CFA"
  ),
  Partner_latent_r = c(
    partner_correlation(qi_dyad_fit1),
    partner_correlation(qi_dyad_fit2),
    partner_correlation(qi_dyad_fit_7)
  )
)

partner_r_table |>
  gt::gt() |>
  gt::tab_header(
    title = "Partner Latent Correlations"
  ) |>
  gt::cols_label(
    Partner_latent_r = "Partner latent r"
  ) |>
  gt::fmt_number(
    columns = Partner_latent_r,
    decimals = 3
  ) |>
  gt::tab_source_note(
    source_note = paste(
      "Similar point estimates do not establish statistical equivalence.",
      "Correlations are interpreted conditional on each measurement model."
    )
  )


## 5.8 Inspect Remaining Local Misfit
#| label: revised-dyadic-mi

lavaan::modindices(qi_dyad_fit2) |>
  filter(op == "~~", lhs != rhs) |>
  arrange(desc(mi)) |>
  select(lhs, op, rhs, mi, epc, sepc.all) |>
  slice_head(n = 10) |>
  gt::gt() |>
  gt::tab_header(
    title = "Remaining Potential Local Misfit",
    subtitle = "Revised nine-item dyadic CFA"
  ) |>
  gt::fmt_number(
    columns = c(mi, epc, sepc.all),
    decimals = 3
  )


# Part 6. Visualize the Findings
## 6.1 Observed Partner Mean Scores
#| label: fig-partner-scores
#| fig-cap: "Observed nine-item mean scores for both partners. Each point represents one conversation; the dashed line marks equal observed means."
#| fig-width: 6
#| fig-height: 5.5

partner_scores <- piq_dyad |>
  mutate(
    mean_A = rowMeans(
      as.matrix(
        pick(all_of(paste0(qi_items_9, "_A")))
      ),
      na.rm = TRUE
    ),
    mean_B = rowMeans(
      as.matrix(
        pick(all_of(paste0(qi_items_9, "_B")))
      ),
      na.rm = TRUE
    )
  ) |>
  filter(
    is.finite(mean_A),
    is.finite(mean_B)
  )

ggplot(
  partner_scores,
  aes(x = mean_A, y = mean_B)
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "grey50"
  ) +
  geom_point(
    alpha = .40,
    colour = "#216869",
    size = 2
  ) +
  scale_x_continuous(
    breaks = 1:7
  ) +
  scale_y_continuous(
    breaks = 1:7
  ) +
  coord_fixed(
    xlim = c(1, 7),
    ylim = c(1, 7)
  ) +
  labs(
    x = "Partner A observed mean QI score",
    y = "Partner B observed mean QI score"
  ) +
  theme_minimal(base_size = 12)


## 6.2 Plot the Revised Nine-Item Dyadic CFA
#| label: fig-dyadic-nine
#| fig-cap: "Standardized measurement paths from the revised nine-item dyadic CFA. The fitted QI8–QI9 residual covariances are omitted from this simplified diagram."
#| fig-width: 14
#| fig-height: 20

semPlot::semPaths(
  object = qi_dyad_fit2,
  what = "paths",
  whatLabels = "std",
  layout = "tree",
  style = "ram",
  rotation = 2,
  nCharNodes = 0,
  intercepts = FALSE,
  residuals = FALSE,
  exoCov = TRUE,
  sizeMan = 5,
  sizeLat = 8,
  edge.label.cex = .8,
  nDigits = 2,
  fade = FALSE,
  mar = c(5, 5, 5, 5)
)


## 6.3 Plot the Seven-Item Dyadic CFA
#| label: fig-dyadic-seven
#| fig-cap: "Standardized estimates from the seven-item dyadic CFA. Residual variances are omitted for readability."
#| fig-width: 14
#| fig-height: 20

semPlot::semPaths(
  object = qi_dyad_fit_7,
  what = "paths",
  whatLabels = "std",
  layout = "tree",
  style = "ram",
  rotation = 2,
  nCharNodes = 0,
  intercepts = FALSE,
  residuals = FALSE,
  exoCov = TRUE,
  sizeMan = 5,
  sizeLat = 8,
  edge.label.cex = .8,
  nDigits = 2,
  fade = FALSE,
  mar = c(5, 5, 5, 5)
)


#| label: model-selection-summary

chisq_d1 <- unname(lavaan::fitMeasures(qi_dyad_fit1, "chisq"))
chisq_d2 <- unname(lavaan::fitMeasures(qi_dyad_fit2, "chisq"))
df_d1 <- unname(lavaan::fitMeasures(qi_dyad_fit1, "df"))
df_d2 <- unname(lavaan::fitMeasures(qi_dyad_fit2, "df"))

chisq_diff_d1_d2 <- chisq_d1 - chisq_d2
df_diff_d1_d2 <- df_d1 - df_d2
p_diff_d1_d2 <- pchisq(
  chisq_diff_d1_d2,
  df = df_diff_d1_d2,
  lower.tail = FALSE
)

selection_table <- dyadic_fit_table |>
  select(Model, chisq, df, cfi, tli, rmsea, srmr) |>
  mutate(Partner_latent_r = partner_r_table$Partner_latent_r)

selection_table |>
  gt::gt() |>
  gt::tab_header(title = "Dyadic Model Selection Summary") |>
  gt::fmt_number(columns = chisq, decimals = 2) |>
  gt::fmt_number(
    columns = c(cfi, tli, rmsea, srmr, Partner_latent_r),
    decimals = 3
  )

# Part 7. Results, Interpretation, and Limitations
## 7.1 Selection of the Final Model
## 7.2 Preliminary Measurement Results
## 7.3 Dyadic Measurement Results
## 7.4 Partner Latent Association
## 7.5 Interpretation of QI8 and QI9
## 7.6 Limitations
