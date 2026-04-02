########################################################################
# Long COVID Symptom Network Analysis — Unified Pipeline
#
# Sections:
#   0.  Setup & output directories
#   1.  User configuration
#   2.  Data loading & preprocessing
#   3.  Helper functions
#   4A. Baseline network (unweighted co-occurrence)
#   4B. IPW computation & Table 1
#   4C. IPW-weighted network
#   4D. Time-bin network analysis
#   4E. Bootstrap & permutation inference
########################################################################

########################################################################
# 0. SETUP
########################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readxl)
  library(igraph)
  library(ggplot2)
  library(ggrepel)
  library(ineq)       # Gini coefficient
  library(tableone)
  library(survey)
})

out_dir <- "C:/Users/young/UConn Health Center/Zeyu/P6_LongCOVID_Network/Output"
dir.create(file.path(out_dir, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "data"),    recursive = TRUE, showWarnings = FALSE)

########################################################################
# 1. USER CONFIGURATION
########################################################################

data_path   <- "C:/Users/young/UConn Health Center/Zeyu/P4_ LongCOVID/data/combined_longcovid.xlsx"

id_col      <- "SUBJ_ID"
infect_col  <- "infected"      # 0 = uninfected, 1 = infected

# Time column for the time-bin section (set to NULL to skip that section)
# Use "months_since_latest" if already in months; use "days_since_infection" if in days
time_col    <- "months_since_latest"

# Symptom columns (positions 29–58 in the spreadsheet; verify after loading)
# These will be overwritten by the auto-detect below; hard-code here if needed
symptom_cols_hardcoded <- c(
  "ALOPECIA","CONSTIPATION","DIARRHEA","DIZZ","HEARING_DIS","LOW_APPE",
  "LOW_CONCENTRATE","NAUSEA","OTH_GI","PAIN","SLE_DIS","ageusia",
  "anosmia","arthritis","chest_pain","chill","cognition_disorder",
  "cough","dyspnea","fatigue","fever","headache","myalgia",
  "palpitation","rash","rhinorrhea_throat","sore_throat","sputum",
  "weakness"
)

# Propensity-score covariates (used for IPW)
ps_vars <- c(
  "Age_group", "SEX", "CCI_group", "BMI_group",
  "smoking", "alcohol",
  "Accumulate_vac_group", "marital_status",
  "househole_member_group", "education_status", "JOB_PREINF"
)
cat_vars <- ps_vars   # all treated as categorical for Table 1

# Bootstrap / permutation iterations
n_boot <- 500
n_perm <- 1000
set_seed_boot <- 123
set_seed_perm <- 456

########################################################################
# 2. DATA LOADING & PREPROCESSING
########################################################################

df_raw <- read_excel(data_path)
message("Loaded: ", nrow(df_raw), " rows, ", ncol(df_raw), " cols")

# -- Derived group variables -------------------------------------------
df <- df_raw %>%
  mutate(
    Accumulate_vac_group = case_when(
      Accumulate_vac == 0              ~ "0",
      Accumulate_vac %in% 1:2          ~ "1-2",
      Accumulate_vac == 3              ~ "3",
      Accumulate_vac == 4              ~ "4",
      Accumulate_vac >= 5              ~ "5+",
      TRUE ~ NA_character_
    ),
    househole_member_group = case_when(
      househole_member == 1 ~ "1",
      househole_member == 2 ~ "2",
      househole_member == 3 ~ "3",
      househole_member == 4 ~ "4",
      househole_member >= 5 ~ "5+",
      TRUE ~ NA_character_
    ),
    BMI_group = case_when(
      BMI < 18.5            ~ "Underweight (<18.5)",
      BMI >= 18.5 & BMI < 23 ~ "Normal (18.5-22.9)",
      BMI >= 23             ~ "Overweight (>=23)",
      TRUE ~ NA_character_
    ),
    CCI_group = case_when(
      CCI_TOTAL == 0              ~ "None (0)",
      CCI_TOTAL >= 1 & CCI_TOTAL <= 2 ~ "Moderate (1-2)",
      CCI_TOTAL >= 3              ~ "Severe/High (>=3)",
      TRUE ~ NA_character_
    ),
    Age_group = case_when(
      AGE >= 20 & AGE <= 39 ~ "20-39",
      AGE >= 40 & AGE <= 59 ~ "40-59",
      AGE >= 60             ~ "60+",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    BMI_group            = factor(BMI_group,
                                  levels = c("Underweight (<18.5)", "Normal (18.5-22.9)", "Overweight (>=23)")),
    CCI_group            = factor(CCI_group,
                                  levels = c("None (0)", "Moderate (1-2)", "Severe/High (>=3)")),
    Age_group            = factor(Age_group, levels = c("20-39", "40-59", "60+")),
    Accumulate_vac_group = factor(Accumulate_vac_group, levels = c("0","1-2","3","4","5+")),
    househole_member_group = factor(househole_member_group, levels = c("1","2","3","4","5+"))
  )

# -- Auto-detect or use hard-coded symptom columns --------------------
symptom_cols <- if (all(symptom_cols_hardcoded %in% names(df))) {
  symptom_cols_hardcoded
} else {
  names(df)[29:58]
}
message("Symptom columns: ", paste(symptom_cols, collapse = ", "))

########################################################################
# 3. HELPER FUNCTIONS
########################################################################

# -- Binary coercion --------------------------------------------------
to01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(!is.na(x) & x != 0))
  x <- tolower(trimws(as.character(x)))
  x[x %in% c("", "na", "n/a", "null")] <- NA_character_
  out <- rep(NA_integer_, length(x))
  out[x %in% c("1","y","yes","true","t","positive","pos","present")] <- 1L
  out[x %in% c("0","n","no","false","f","negative","neg","absent")]  <- 0L
  out[is.na(out)] <- 0L
  out
}

# -- Symptom display names --------------------------------------------
clean_symptom_names <- function(x) {
  dplyr::recode(x,
    "LOW_CONCENTRATE"  = "Difficulty concentrating",
    "SLE_DIS"          = "Sleep disturbance",
    "LOW_APPE"         = "Loss of appetite",
    "HEARING_DIS"      = "Hearing difficulty",
    "OTH_GI"           = "Other GI",
    "PAIN"             = "Pain",
    "DIZZ"             = "Dizziness",
    "ALOPECIA"         = "Hair loss",
    "ageusia"          = "Loss of taste",
    "anosmia"          = "Loss of smell",
    "cognition_disorder" = "Cognitive impairment",
    "chest_pain"       = "Chest pain",
    "rhinorrhea_throat" = "Runny nose/throat",
    "sore_throat"      = "Sore throat",
    "dyspnea"          = "Shortness of breath",
    "palpitation"      = "Palpitations",
    "myalgia"          = "Muscle pain",
    "fatigue"          = "Fatigue",
    "weakness"         = "Weakness",
    "headache"         = "Headache",
    "fever"            = "Fever",
    "cough"            = "Cough",
    "sputum"           = "Sputum",
    "rash"             = "Rash",
    "chill"            = "Chills",
    "arthritis"        = "Joint pain",
    "CONSTIPATION"     = "Constipation",
    "DIARRHEA"         = "Diarrhea",
    "NAUSEA"           = "Nausea",
    .default = str_replace_all(x, "_", " ")
  )
}

# -- Symptom category -------------------------------------------------
symptom_category <- function(x) {
  case_when(
    x %in% c("cough","sputum","dyspnea","rhinorrhea_throat","sore_throat") ~ "Respiratory",
    x %in% c("ageusia","anosmia","LOW_CONCENTRATE","cognition_disorder",
              "DIZZ","headache","HEARING_DIS","SLE_DIS") ~ "Neurologic/Sensory",
    x %in% c("NAUSEA","DIARRHEA","CONSTIPATION","OTH_GI","LOW_APPE") ~ "Gastrointestinal",
    x %in% c("fatigue","weakness","myalgia","fever","chill") ~ "Systemic",
    x %in% c("chest_pain","palpitation") ~ "Cardiopulmonary",
    TRUE ~ "Other"
  )
}

# -- Weighted network density (from adjacency matrix) -----------------
weighted_density <- function(A) {
  n <- nrow(A)
  if (n < 2) return(NA_real_)
  sum(A[upper.tri(A)], na.rm = TRUE) / choose(n, 2)
}

# -- Weighted strength assortativity ----------------------------------
weighted_assortativity_strength <- function(g) {
  edf <- as_data_frame(g, what = "edges")
  if (nrow(edf) == 0) return(NA_real_)
  s  <- strength(g, weights = E(g)$weight)
  su <- s[edf$from];  sv <- s[edf$to];  w <- edf$weight
  mu_u <- weighted.mean(su, w, na.rm = TRUE)
  mu_v <- weighted.mean(sv, w, na.rm = TRUE)
  cov_w <- sum(w * (su - mu_u) * (sv - mu_v), na.rm = TRUE) / sum(w, na.rm = TRUE)
  sd_u  <- sqrt(sum(w * (su - mu_u)^2, na.rm = TRUE) / sum(w, na.rm = TRUE))
  sd_v  <- sqrt(sum(w * (sv - mu_v)^2, na.rm = TRUE) / sum(w, na.rm = TRUE))
  if (sd_u == 0 || sd_v == 0) return(NA_real_)
  cov_w / (sd_u * sd_v)
}

# -- Global network metrics (igraph object) ---------------------------
global_network_metrics <- function(g) {
  if (is.null(g) || vcount(g) < 2 || ecount(g) == 0) {
    return(tibble(nodes=NA, edges=NA, weighted_density=NA,
                  mean_strength=NA, max_strength=NA, strength_gini=NA,
                  mean_pagerank=NA, max_pagerank=NA, pagerank_gini=NA,
                  n_communities=NA, modularity=NA, assortativity_strength=NA))
  }
  st   <- strength(g, weights = E(g)$weight)
  pr   <- page_rank(g, weights = E(g)$weight)$vector
  comm <- cluster_louvain(g, weights = E(g)$weight)
  n    <- vcount(g)
  tibble(
    nodes                = n,
    edges                = ecount(g),
    weighted_density     = sum(E(g)$weight) / (n * (n - 1) / 2),
    mean_strength        = mean(st),
    max_strength         = max(st),
    strength_gini        = ineq(st,  type = "Gini"),
    mean_pagerank        = mean(pr),
    max_pagerank         = max(pr),
    pagerank_gini        = ineq(pr,  type = "Gini"),
    n_communities        = length(unique(membership(comm))),
    modularity           = modularity(comm),
    assortativity_strength = weighted_assortativity_strength(g)
  )
}

# -- Node-level metrics -----------------------------------------------
calc_node_metrics <- function(g, group_label) {
  if (is.null(g) || vcount(g) == 0) return(tibble())
  if (ecount(g) == 0) {
    return(tibble(
      symptom   = V(g)$name,
      strength  = 0,
      pagerank  = 1 / vcount(g),
      community = NA_integer_,
      group     = group_label
    ))
  }
  comm <- cluster_louvain(g, weights = E(g)$weight)
  tibble(
    symptom   = V(g)$name,
    strength  = strength(g,   weights = E(g)$weight),
    pagerank  = page_rank(g,  weights = E(g)$weight)$vector,
    community = as.integer(membership(comm)),
    group     = group_label
  )
}

########################################################################
# 4A. BASELINE NETWORK (unweighted co-occurrence)
########################################################################

# -- Build adjacency matrix (raw co-occurrence counts) ----------------
# weight_col = NULL → each person contributes weight 1
build_adj <- function(data, symptom_cols, weight_col = NULL,
                      normalize = c("none", "jaccard", "cosine", "expected")) {
  normalize <- match.arg(normalize)
  X <- data %>% select(all_of(symptom_cols)) %>%
       mutate(across(everything(), to01)) %>% as.matrix()
  w <- if (is.null(weight_col)) rep(1, nrow(data)) else data[[weight_col]]
  w <- ifelse(is.na(w), 0, w)
  Xw  <- X * sqrt(w)
  O   <- t(Xw) %*% Xw;  diag(O) <- 0
  prev <- colSums(X * w);  Nw <- sum(w)
  A <- O
  if (normalize == "jaccard") {
    J <- matrix(0, ncol(X), ncol(X), dimnames = list(symptom_cols, symptom_cols))
    for (i in seq_len(ncol(X))) for (j in seq_len(ncol(X))) if (i != j) {
      u <- prev[i] + prev[j] - O[i, j]
      J[i, j] <- ifelse(u > 0, O[i, j] / u, 0)
    }
    A <- J
  }
  if (normalize == "cosine") {
    C <- matrix(0, ncol(X), ncol(X), dimnames = list(symptom_cols, symptom_cols))
    for (i in seq_len(ncol(X))) for (j in seq_len(ncol(X))) if (i != j) {
      d <- sqrt(prev[i] * prev[j])
      C[i, j] <- ifelse(d > 0, O[i, j] / d, 0)
    }
    A <- C
  }
  if (normalize == "expected") {
    p <- prev / Nw;  E_mat <- outer(p, p) * Nw;  diag(E_mat) <- 0
    A <- ifelse(E_mat > 0, O / E_mat, 0);  diag(A) <- 0
  }
  dimnames(A) <- list(symptom_cols, symptom_cols)
  A
}

# -- Compute network metrics from adjacency matrix --------------------
compute_network_metrics_adj <- function(A) {
  g <- graph_from_adjacency_matrix(A, mode = "undirected", weighted = TRUE, diag = FALSE)
  if (vcount(g) < 2 || ecount(g) == 0)
    return(list(graph = g, pagerank = NA, strength = NA, community = NA,
                modularity = NA_real_, density = NA_real_, assortativity = NA_real_))
  pr   <- page_rank(g, directed = FALSE, weights = E(g)$weight)$vector
  st   <- strength(g, weights = E(g)$weight)
  cl   <- cluster_louvain(g, weights = E(g)$weight)
  list(graph = g, pagerank = pr, strength = st,
       community = membership(cl), modularity = modularity(cl),
       density = weighted_density(A),
       assortativity = weighted_assortativity_strength(g))
}

# -- Fit two-group networks from adjacency matrices -------------------
fit_two_group_networks <- function(data, symptom_cols, group_col = "infected",
                                   weight_col = NULL, normalize = "none") {
  dat0 <- data %>% filter(.data[[group_col]] == 0)
  dat1 <- data %>% filter(.data[[group_col]] == 1)
  A0   <- build_adj(dat0, symptom_cols, weight_col, normalize)
  A1   <- build_adj(dat1, symptom_cols, weight_col, normalize)
  m0   <- compute_network_metrics_adj(A0)
  m1   <- compute_network_metrics_adj(A1)
  node_df <- tibble(
    symptom        = symptom_cols,
    strength_0     = unname(m0$strength[symptom_cols]),
    strength_1     = unname(m1$strength[symptom_cols]),
    pagerank_0     = unname(m0$pagerank[symptom_cols]),
    pagerank_1     = unname(m1$pagerank[symptom_cols]),
    community_0    = unname(m0$community[symptom_cols]),
    community_1    = unname(m1$community[symptom_cols]),
    delta_strength = strength_1 - strength_0,
    delta_pagerank = pagerank_1 - pagerank_0
  )%>%
    mutate(symptom_clean = clean_symptom_names(symptom),
           category      = symptom_category(symptom))
  global_df <- tibble(
    density_0       = m0$density,       density_1       = m1$density,
    delta_density   = m1$density   - m0$density,
    modularity_0    = m0$modularity,    modularity_1    = m1$modularity,
    delta_modularity = m1$modularity - m0$modularity,
    assortativity_0 = m0$assortativity, assortativity_1 = m1$assortativity,
    delta_assortativity = m1$assortativity - m0$assortativity
  )
  list(node = node_df, global = global_df, A0 = A0, A1 = A1, m0 = m0, m1 = m1)
}

# ==========  Run baseline analysis  ==================================

df_clean <- df %>%
  select(all_of(c(id_col, infect_col, symptom_cols))) %>%
  mutate(across(all_of(symptom_cols), to01)) %>%
  group_by(.data[[id_col]]) %>%
  summarise(
    !!infect_col := dplyr::first(.data[[infect_col]]),
    across(all_of(symptom_cols), ~ max(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(!!infect_col := as.integer(as.character(.data[[infect_col]])))

# --- Diagnostic: verify group split ----------------------------------
message("\n--- Baseline group counts ---")
print(table(df_clean[[infect_col]], useNA = "ifany"))
stopifnot("No uninfected (0) rows found!" = sum(df_clean[[infect_col]] == 0, na.rm = TRUE) > 0)
stopifnot("No infected (1) rows found!"   = sum(df_clean[[infect_col]] == 1, na.rm = TRUE) > 0)

baseline_fit <- fit_two_group_networks(
  data         = df_clean,
  symptom_cols = symptom_cols,
  group_col    = infect_col,
  weight_col   = NULL,
  normalize    = "none"
)

# Node-level baseline table
baseline_node <- baseline_fit$node %>%
  mutate(z_strength = as.numeric(scale(delta_strength)),
         z_pagerank = as.numeric(scale(delta_pagerank)),
         change_score = abs(z_strength) + abs(z_pagerank)) %>%
  arrange(desc(delta_pagerank))

# Global baseline metrics
baseline_global <- baseline_fit$global

message("\n=== BASELINE GLOBAL METRICS ===")
print(as.data.frame(t(baseline_global)))

# Build igraph objects for plotting
g_base_inf  <- baseline_fit$m1$graph
g_base_ctrl <- baseline_fit$m0$graph

# Comprehensive global metrics table (using igraph helper)
global_metrics_table <- bind_rows(
  global_network_metrics(g_base_ctrl) %>% mutate(group = "Uninfected"),
  global_network_metrics(g_base_inf)  %>% mutate(group = "Infected")
) %>% relocate(group)



print(global_metrics_table)
write.csv(global_metrics_table,
          file.path(out_dir, "tables", "baseline_global_metrics.csv"),
          row.names = FALSE)

# Node comparison table
# Replace the current node_compare_baseline block with this

node_compare_baseline <- baseline_fit$node %>%
  dplyr::select(
    symptom,
    strength_0,
    strength_1,
    pagerank_0,
    pagerank_1,
    community_0,
    community_1,
    delta_strength,
    delta_pagerank,
    symptom_clean,
    category
  ) %>%
  dplyr::rename(
    strength_Uninfected   = strength_0,
    strength_Infected     = strength_1,
    pagerank_Uninfected   = pagerank_0,
    pagerank_Infected     = pagerank_1,
    community_Uninfected  = community_0,
    community_Infected    = community_1
  ) %>%
  arrange(desc(delta_pagerank))

write.csv(
  node_compare_baseline,
  file.path(out_dir, "tables", "baseline_node_comparison.csv"),
  row.names = FALSE
)

print(head(node_compare_baseline, 15))



summary(baseline_fit$m0$strength)
summary(baseline_fit$m1$strength)
table(baseline_fit$m0$strength == 0, useNA = "ifany")
table(baseline_fit$m1$strength == 0, useNA = "ifany")

# Scatter: delta strength vs delta pagerank
p_scatter_base <- ggplot(baseline_node,
       aes(x = z_strength, y = z_pagerank,
           color = category, label = symptom_clean)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(alpha = 0.85, size = 2.5) +
  geom_text_repel(size = 3.5, max.overlaps = 100) +
  labs(x = "Delta Strength (z-score)", y = "Delta PageRank (z-score)",
       title = "Baseline: Differential symptom network rewiring (Infected - Uninfected)",
       color = "Category") +
  theme_classic(base_size = 12)
ggsave(file.path(out_dir, "figures", "baseline_delta_scatter.png"),
       p_scatter_base, width = 10, height = 7, dpi = 300)
print(p_scatter_base)

########################################################################
# 4B. IPW COMPUTATION & TABLE 1
########################################################################

# -- Build propensity-score dataset -----------------------------------
df_ps <- df %>%
  select(all_of(c(id_col, infect_col, ps_vars))) %>%
  filter(complete.cases(.)) %>%
  mutate(
    infected             = as.integer(as.character(.data[[infect_col]])),
    Age_group            = factor(Age_group),
    SEX                  = factor(SEX),
    CCI_group            = factor(CCI_group),
    BMI_group            = factor(BMI_group),
    smoking              = factor(smoking),
    alcohol              = factor(alcohol),
    Accumulate_vac_group = factor(Accumulate_vac_group),
    marital_status       = factor(marital_status),
    househole_member_group = factor(househole_member_group),
    education_status     = factor(education_status),
    JOB_PREINF           = factor(JOB_PREINF)
  )

message("PS dataset: ", nrow(df_ps), " complete rows")
table(df_ps$infected)

# -- Propensity score model -------------------------------------------
ps_formula <- as.formula(paste(
  infect_col, "~",
  paste(ps_vars, collapse = " + ")
))

ps_model <- glm(ps_formula, data = df_ps, family = binomial())
summary(ps_model)

df_ps$pscore <- predict(ps_model, type = "response")
summary(df_ps$pscore)

# -- Stabilized IPW (truncated at 99th percentile) --------------------
p_infect <- mean(df_ps$infected == 1)
df_ps <- df_ps %>%
  mutate(
    ipw = ifelse(
      infected == 1,
      p_infect / pscore,
      (1 - p_infect) / (1 - pscore)
    ),
    ipw = pmin(ipw, quantile(ipw, 0.99, na.rm = TRUE))
  )

summary(df_ps$ipw)

# -- Table 1 (unweighted) --------------------------------------------
df_table1 <- df_ps %>%
  group_by(.data[[id_col]]) %>%
  slice(1) %>% ungroup() %>%
  mutate(infected = factor(infected, levels = c(0, 1),
                           labels = c("Uninfected", "Infected")))

tab_unw <- CreateTableOne(
  vars      = ps_vars,
  strata    = infect_col,
  data      = df_table1,
  factorVars = cat_vars,
  test      = FALSE
)
print(tab_unw, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)

# -- Table 1 (IPW-weighted) ------------------------------------------
des_w  <- svydesign(ids = ~1, weights = ~ipw, data = df_table1)
tab_w  <- svyCreateTableOne(
  vars      = ps_vars,
  strata    = infect_col,
  data      = des_w,
  factorVars = cat_vars,
  test      = FALSE
)
print(tab_w, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)

# -- Save Table 1 to CSV ----------------------------------------------
tab_unw_df <- as.data.frame(print(tab_unw, showAllLevels=TRUE, quote=FALSE, noSpaces=TRUE, test=FALSE))
tab_w_df   <- as.data.frame(print(tab_w,   showAllLevels=TRUE, quote=FALSE, noSpaces=TRUE, test=FALSE))

write.csv(tab_unw_df, file.path(out_dir, "tables", "Table1_unweighted.csv"), row.names = TRUE)
write.csv(tab_w_df,   file.path(out_dir, "tables", "Table1_weighted.csv"),   row.names = TRUE)

########################################################################
# 4C. IPW-WEIGHTED NETWORK
########################################################################

# -- IPW network builder (edge weight = sum of IPW for co-reporters) --
build_symptom_network_ipw <- function(df_sub, id_col, symptom_cols,
                                      weight_col = "ipw",
                                      na_weight_action = c("drop","one","zero")) {
  na_weight_action <- match.arg(na_weight_action)
  required <- c(id_col, symptom_cols, weight_col)
  missing  <- setdiff(required, names(df_sub))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse=", "))

  # Person-level weights (one row per person)
  w_by_id <- df_sub %>%
    select(all_of(c(id_col, weight_col))) %>%
    group_by(across(all_of(id_col))) %>%
    summarise(w = { ww <- .data[[weight_col]][!is.na(.data[[weight_col]])];
                    if (length(ww) == 0) NA_real_ else ww[1] }, .groups = "drop")

  if (na_weight_action == "drop")  w_by_id <- w_by_id %>% filter(!is.na(w))
  if (na_weight_action == "one")   w_by_id <- w_by_id %>% mutate(w = if_else(is.na(w), 1, w))
  if (na_weight_action == "zero")  w_by_id <- w_by_id %>% mutate(w = if_else(is.na(w), 0, w))

  # Long symptom table (positive only)
  long_sym <- df_sub %>%
    pivot_longer(cols = all_of(symptom_cols), names_to = "symptom", values_to = "value_raw") %>%
    mutate(value = to01(value_raw)) %>%
    filter(value == 1) %>%
    select(all_of(id_col), symptom) %>%
    distinct() %>%
    left_join(w_by_id %>% select(all_of(id_col), w), by = id_col)
  if (na_weight_action == "drop") long_sym <- long_sym %>% filter(!is.na(w))

  if (nrow(long_sym) == 0) {
    g <- make_empty_graph(n = length(symptom_cols), directed = FALSE)
    V(g)$name <- symptom_cols
    return(list(graph = g, edges = tibble(), long_sym = long_sym))
  }

  # Per-person symptom sets + IPW weight
  person_syms <- long_sym %>%
    group_by(across(all_of(id_col))) %>%
    summarise(symptoms = list(sort(unique(symptom))), w = first(w),
              k = length(unique(symptom)), .groups = "drop")

  # Weighted edge list
  pairs_raw <- person_syms %>%
    filter(k >= 2) %>%
    mutate(pairs = map(symptoms, ~{
      cmb <- combn(.x, 2)
      tibble(source = cmb[1,], target = cmb[2,])
    })) %>%
    select(all_of(id_col), w, pairs) %>%
    unnest(pairs)

  edge_raw <- if (nrow(pairs_raw) == 0) {
    tibble(source=character(), target=character(), strength_n=numeric(), strength_w=numeric())
  } else {
    pairs_raw %>%
      group_by(source, target) %>%
      summarise(strength_n = n_distinct(.data[[id_col]]),
                strength_w = sum(w, na.rm = TRUE), .groups = "drop")
  }

  # Build igraph (include all symptom nodes even if isolated)
  vertices <- tibble(name = symptom_cols)
  g <- graph_from_data_frame(edge_raw, directed = FALSE, vertices = vertices)
  if (nrow(edge_raw) > 0) {
    E(g)$weight <- edge_raw$strength_w
    E(g)$n      <- edge_raw$strength_n
  }
  list(graph = g, edges = edge_raw, long_sym = long_sym)
}

# ==========  Merge IPW weights into main data  =======================

df_ipw <- df %>%
  left_join(df_ps %>% select(all_of(c(id_col, "ipw"))), by = id_col) %>%
  select(all_of(c(id_col, infect_col, symptom_cols, "ipw"))) %>%
  mutate(across(all_of(symptom_cols), to01)) %>%
  group_by(.data[[id_col]]) %>%
  summarise(
    !!infect_col := dplyr::first(.data[[infect_col]]),
    ipw = {
      tmp <- ipw[!is.na(ipw)]
      if (length(tmp) == 0) NA_real_ else tmp[1]
    },
    across(all_of(symptom_cols), ~ max(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    !!infect_col := as.integer(as.character(.data[[infect_col]])),
    ipw = ifelse(is.na(ipw), 1, ipw)
  )

# --- Diagnostic: verify IPW group split ------------------------------
message("\n--- IPW group counts ---")
print(table(df_ipw[[infect_col]], useNA = "ifany"))
message("IPW non-missing: ", sum(!is.na(df_ipw$ipw)),
        " / total: ", nrow(df_ipw))
message("IPW non-missing by group:")
print(df_ipw %>% group_by(.data[[infect_col]]) %>%
        summarise(n = n(), n_ipw = sum(!is.na(ipw)), .groups = "drop"))

# -- Build IPW-weighted networks for each group -----------------------
# na_weight_action = "one": people missing IPW (dropped by complete.cases
# in the PS model) still contribute with weight = 1.
# Change to "drop" if you want to exclude them entirely.
net_inf_ipw  <- build_symptom_network_ipw(
  df_ipw %>% filter(.data[[infect_col]] == 1),
  id_col, symptom_cols, weight_col = "ipw", na_weight_action = "one")

net_ctrl_ipw <- build_symptom_network_ipw(
  df_ipw %>% filter(.data[[infect_col]] == 0),
  id_col, symptom_cols, weight_col = "ipw", na_weight_action = "one")

message("Infected IPW network:   ", ecount(net_inf_ipw$graph),  " edges, ",
        vcount(net_inf_ipw$graph),  " nodes")
message("Uninfected IPW network: ", ecount(net_ctrl_ipw$graph), " edges, ",
        vcount(net_ctrl_ipw$graph), " nodes")

g_inf_ipw  <- net_inf_ipw$graph
g_ctrl_ipw <- net_ctrl_ipw$graph

# -- Node-level IPW metrics -------------------------------------------
# -- Node-level IPW metrics -------------------------------------------
node_inf_ipw  <- calc_node_metrics(g_inf_ipw,  "Infected_IPW")
node_ctrl_ipw <- calc_node_metrics(g_ctrl_ipw, "Uninfected_IPW")

message("Summary of infected IPW node strengths:")
print(summary(node_inf_ipw$strength))

message("Summary of uninfected IPW node strengths:")
print(summary(node_ctrl_ipw$strength))

message("Number of edges in infected IPW graph: ", igraph::ecount(g_inf_ipw))
message("Number of edges in uninfected IPW graph: ", igraph::ecount(g_ctrl_ipw))

driver_table_ipw <- bind_rows(node_inf_ipw, node_ctrl_ipw) %>%
  pivot_wider(
    names_from = group,
    values_from = c(strength, pagerank, community)
  ) %>%
  mutate(
    delta_strength_ipw = strength_Infected_IPW - strength_Uninfected_IPW,
    delta_pagerank_ipw = pagerank_Infected_IPW - pagerank_Uninfected_IPW
  ) %>%
  arrange(desc(delta_pagerank_ipw))

print(head(driver_table_ipw, 15))
summary(driver_table_ipw$strength_Uninfected_IPW)
table(driver_table_ipw$strength_Uninfected_IPW == 0, useNA = "ifany")

message("\n=== TOP IPW DRIVER SYMPTOMS ===")
print(head(driver_table_ipw, 15))

write.csv(
  driver_table_ipw,
  file.path(out_dir, "tables", "ipw_node_driver_table.csv"),
  row.names = FALSE
)




# -- Global IPW metrics -----------------------------------------------
ipw_global <- bind_rows(
  global_network_metrics(g_ctrl_ipw) %>% mutate(group = "Uninfected_IPW"),
  global_network_metrics(g_inf_ipw)  %>% mutate(group = "Infected_IPW")
) %>% relocate(group)

print(ipw_global)
write.csv(ipw_global,
          file.path(out_dir, "tables", "ipw_global_metrics.csv"),
          row.names = FALSE)

# -- Node-level IPW metrics -------------------------------------------

########################################################################
# 4D. TIME-BIN NETWORK ANALYSIS
########################################################################

if (!is.null(time_col) && time_col %in% names(df)) {

  is_days <- grepl("day", time_col, ignore.case = TRUE)

  df2 <- df %>%
    mutate(
      infected01  = to01(.data[[infect_col]]),
      time_raw    = suppressWarnings(as.numeric(.data[[time_col]])),
      months_since = if (is_days) time_raw / 30.44 else time_raw,
      timebin = cut(
        months_since,
        breaks = c(0, 3, 6, 12, 24, 50),
        include.lowest = TRUE, right = TRUE,
        labels = c("0-3","3-6","6-12","12-24","24+")
      )
    )

  message("Time-bin distribution:")
  print(df2 %>% count(timebin, infected01))

  # -- Time-sliced network builder ------------------------------------
  build_symptom_network_timeslice <- function(df_sub, id_col, symptom_cols,
                                              timebin_col = "timebin") {
    df_sub <- df_sub %>%
      filter(!is.na(.data[[timebin_col]])) %>%
      mutate(unit_id = paste0(.data[[id_col]], "_", .data[[timebin_col]]))
    long_sym <- df_sub %>%
      pivot_longer(cols = all_of(symptom_cols), names_to = "symptom", values_to = "value_raw") %>%
      mutate(value = to01(value_raw)) %>%
      filter(value == 1) %>%
      select(unit_id, symptom) %>% distinct()
    if (nrow(long_sym) == 0) {
      return(list(graph = make_empty_graph(), edges = tibble(), long_sym = long_sym, comm = NULL))
    }
    pairs <- long_sym %>%
      inner_join(long_sym, by = "unit_id", suffix = c("_1","_2")) %>%
      filter(symptom_1 < symptom_2) %>%
      transmute(source = symptom_1, target = symptom_2, unit_id)
    edge_raw <- pairs %>%
      group_by(source, target) %>%
      summarise(weight = n_distinct(unit_id), .groups = "drop")
    g <- graph_from_data_frame(edge_raw, directed = FALSE)
    E(g)$weight <- edge_raw$weight
    comm <- cluster_louvain(g, weights = E(g)$weight)
    V(g)$community <- membership(comm)
    V(g)$strength  <- strength(g, weights = E(g)$weight)
    list(graph = g, edges = edge_raw, long_sym = long_sym, comm = comm)
  }

  pagerank_gini_safe <- function(g) {
    if (is.null(g) || ecount(g) == 0) return(NA_real_)
    pr <- page_rank(g, weights = E(g)$weight)$vector
    ineq(pr, type = "Gini")
  }

  modularity_safe <- function(comm_obj) {
    if (is.null(comm_obj)) return(NA_real_)
    modularity(comm_obj)
  }

  # -- Bootstrap CIs for time-bin network metrics -----------------------
  bootstrap_timebin_metrics <- function(df2, id_col, symptom_cols,
                                        infect_col = "infected01",
                                        timebin_col = "timebin",
                                        B = 300,
                                        seed = 123) {
    set.seed(seed)
    
    tb_levels <- c("0-3", "3-6", "6-12", "12-24")
    out <- list()
    k <- 1
    
    for (grp in c(0, 1)) {
      for (tb in tb_levels) {
        
        df_sub <- df2 %>%
          filter(.data[[infect_col]] == grp,
                 as.character(.data[[timebin_col]]) == tb)
        
        n_sub <- nrow(df_sub)
        if (n_sub < 2) next
        
        boot_list <- vector("list", B)
        
        for (b in seq_len(B)) {
          samp <- df_sub[sample(seq_len(n_sub), size = n_sub, replace = TRUE), , drop = FALSE]
          
          net <- build_symptom_network_timeslice(samp, id_col, symptom_cols, timebin_col)
          g   <- net$graph
          
          boot_list[[b]] <- tibble(
            infected = grp,
            timebin  = tb,
            weighted_density = if (ecount(g) > 0) {
              sum(E(g)$weight) / (vcount(g) * (vcount(g) - 1) / 2)
            } else NA_real_,
            modularity = modularity_safe(net$comm),
            pagerank_gini = pagerank_gini_safe(g),
            strength_assortativity = if (ecount(g) > 0) weighted_assortativity_strength(g) else NA_real_
          )
        }
        
        boot_df <- bind_rows(boot_list)
        
        out[[k]] <- boot_df %>%
          summarise(
            infected = grp,
            timebin  = tb,
            
            wd_lcl = quantile(weighted_density, 0.025, na.rm = TRUE),
            wd_ucl = quantile(weighted_density, 0.975, na.rm = TRUE),
            
            mod_lcl = quantile(modularity, 0.025, na.rm = TRUE),
            mod_ucl = quantile(modularity, 0.975, na.rm = TRUE),
            
            pg_lcl = quantile(pagerank_gini, 0.025, na.rm = TRUE),
            pg_ucl = quantile(pagerank_gini, 0.975, na.rm = TRUE),
            
            assort_lcl = quantile(strength_assortativity, 0.025, na.rm = TRUE),
            assort_ucl = quantile(strength_assortativity, 0.975, na.rm = TRUE)
          )
        
        k <- k + 1
      }
    }
    
    bind_rows(out)
  }
  
  
  
  
  
  
  # -- Loop over groups × timebins ------------------------------------
  timebins <- c("0-3","3-6","6-12","12-24")
  results_metrics <- list()
  results_nodes   <- list()

  for (grp in c(0, 1)) {
    for (tb in timebins) {
      df_sub <- df2 %>% filter(infected01 == grp, as.character(timebin) == tb)
      net    <- build_symptom_network_timeslice(df_sub, id_col, symptom_cols, "timebin")
      g      <- net$graph

      met <- tibble(
        infected               = grp,
        timebin                = tb,
        n_units                = n_distinct(paste0(df_sub[[id_col]], "_", tb)),
        nodes                  = vcount(g),
        edges                  = ecount(g),
        weighted_density       = if (ecount(g) > 0) sum(E(g)$weight) / (vcount(g) * (vcount(g)-1) / 2) else NA_real_,
        modularity             = modularity_safe(net$comm),
        pagerank_gini          = pagerank_gini_safe(g),
        strength_assortativity = if (ecount(g) > 0) weighted_assortativity_strength(g) else NA_real_
      )
      results_metrics[[paste0("g", grp, "_", tb)]] <- met

      if (ecount(g) > 0) {
        pr <- page_rank(g, weights = E(g)$weight)$vector
        st <- strength(g, weights = E(g)$weight)
        results_nodes[[paste0("nodes_g", grp, "_", tb)]] <- tibble(
          infected  = grp, timebin = tb,
          symptom   = names(pr),
          pagerank  = as.numeric(pr),
          strength  = as.numeric(st),
          community = as.integer(V(g)$community)
        )
      }
    }
  }

  metrics_df <- bind_rows(results_metrics) %>%
    mutate(group  = ifelse(infected == 1, "Infected", "Not infected"),
           timebin = factor(timebin, levels = timebins))
  

  
  
  # -- Bootstrap CIs for time-bin metrics -------------------------------
  message("Running bootstrap for time-bin metrics ...")
  metrics_ci <- bootstrap_timebin_metrics(
    df2 = df2,
    id_col = id_col,
    symptom_cols = symptom_cols,
    infect_col = "infected01",
    timebin_col = "timebin",
    B = 300,   # use 100 first for testing, then 300-500 for final
    seed = 123
  )
  
  metrics_df <- metrics_df %>%
    left_join(metrics_ci, by = c("infected", "timebin"))
  
  metrics_df_formatted <- metrics_df %>%
    mutate(
      weighted_density_ci = sprintf("%.2f (%.2f, %.2f)", weighted_density, wd_lcl, wd_ucl),
      modularity_ci = sprintf("%.4f (%.4f, %.4f)", modularity, mod_lcl, mod_ucl),
      pagerank_gini_ci = sprintf("%.3f (%.3f, %.3f)", pagerank_gini, pg_lcl, pg_ucl),
      strength_assortativity_ci = sprintf("%.4f (%.4f, %.4f)", strength_assortativity, assort_lcl, assort_ucl)
    ) %>%
    select(
      infected, timebin, n_units, nodes, edges, group,
      weighted_density_ci,
      modularity_ci,
      pagerank_gini_ci,
      strength_assortativity_ci
    )
  
  print(metrics_df_formatted)
  
  write.csv(
    metrics_df_formatted,
    file.path(out_dir, "tables", "metrics_timebin_by_group_with_bootstrap_ci.csv"),
    row.names = FALSE
  )

  nodes_df <- bind_rows(results_nodes) %>%
    mutate(group  = ifelse(infected == 1, "Infected", "Not infected"),
           timebin = factor(as.character(timebin), levels = timebins),
           community = factor(community))

  # Save
  saveRDS(metrics_df, file.path(out_dir, "data", "metrics_timebin_by_group.rds"))
  saveRDS(nodes_df,   file.path(out_dir, "data", "nodes_timebin_by_group.rds"))
  write.csv(metrics_df, file.path(out_dir, "tables", "metrics_timebin_by_group.csv"), row.names = FALSE)

  message("\n=== TIME-BIN METRICS ===")
  print(metrics_df)

  # -- Temporal trend plots -------------------------------------------
  p_density <- ggplot(metrics_df, aes(x = timebin, y = weighted_density,
                                       group = group, color = group)) +
    geom_line() + geom_point(size = 2) +
    theme_minimal(base_size = 13) +
    labs(x = "Months since infection", y = "Weighted density",
         title = "Temporal change in network coupling", color = NULL)

  p_mod <- ggplot(metrics_df, aes(x = timebin, y = modularity,
                                   group = group, color = group)) +
    geom_line() + geom_point(size = 2) +
    theme_minimal(base_size = 13) +
    labs(x = "Months since infection", y = "Modularity (Louvain)",
         title = "Temporal change in community separation", color = NULL)

  p_gini <- ggplot(metrics_df, aes(x = timebin, y = pagerank_gini,
                                    group = group, color = group)) +
    geom_line() + geom_point(size = 2) +
    theme_minimal(base_size = 13) +
    labs(x = "Months since infection", y = "Gini(PageRank)",
         title = "Temporal change in hub dominance", color = NULL)

  ggsave(file.path(out_dir, "figures", "trend_density.png"),    p_density, width=8, height=4.5, dpi=300)
  ggsave(file.path(out_dir, "figures", "trend_modularity.png"), p_mod,     width=8, height=4.5, dpi=300)
  ggsave(file.path(out_dir, "figures", "trend_gini.png"),       p_gini,    width=8, height=4.5, dpi=300)

  print(p_density); print(p_mod); print(p_gini)

  # -- Lollipop plots per group × timebin ----------------------------
  plot_lollipop_timebin <- function(nodes_df, grp, tb, K = 15, max_pr = NULL) {
    dat <- nodes_df %>%
      filter(group == grp, timebin == tb) %>%
      arrange(desc(pagerank)) %>% slice(1:K) %>%
      mutate(symptom = factor(symptom, levels = rev(symptom)))
    if (nrow(dat) == 0) return(NULL)
    p <- ggplot(dat, aes(x = pagerank, y = symptom, color = community)) +
      geom_segment(aes(x = 0, xend = pagerank, yend = symptom), linewidth = 1.1, alpha = 0.85) +
      geom_point(aes(size = strength), alpha = 0.95) +
      theme_minimal(base_size = 13) +
      labs(title  = paste0(grp, " | ", tb, " months: Top symptoms by PageRank"),
           subtitle = "Line = PageRank (global influence); point size = strength (local)",
           x = "Weighted PageRank", y = NULL, color = "Community", size = "Strength")
    if (!is.null(max_pr)) p <- p + coord_cartesian(xlim = c(0, max_pr))
    p
  }

  max_pr <- max(nodes_df$pagerank, na.rm = TRUE)
  lollipop_dir <- file.path(out_dir, "figures", "lollipop_timebin")
  dir.create(lollipop_dir, recursive = TRUE, showWarnings = FALSE)

  for (grp in c("Infected", "Not infected")) {
    for (tb in timebins) {
      p <- plot_lollipop_timebin(nodes_df, grp, tb, K = 15, max_pr = max_pr)
      if (!is.null(p)) {
        fname <- paste0("lollipop_", gsub(" ", "", tolower(grp)), "_", tb, ".png")
        ggsave(file.path(lollipop_dir, fname), p, width = 7.5, height = 6, dpi = 300)
      }
    }
  }

  # -- Heatmap: node metrics over time --------------------------------
  zscore <- function(x) as.numeric(scale(x))

  make_metric_heatmap <- function(dat, metric_col, grp_label) {
    dat_g <- dat %>%
      filter(group == grp_label, !is.na(timebin)) %>%
      group_by(symptom, timebin) %>%
      summarise(value = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop") %>%
      mutate(value_z = zscore(value))
    sym_order <- dat_g %>%
      group_by(symptom) %>%
      summarise(overall = mean(value, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(overall)) %>% pull(symptom)
    dat_g <- dat_g %>% mutate(symptom = factor(symptom, levels = rev(sym_order)))
    ggplot(dat_g, aes(x = timebin, y = symptom, fill = value_z)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(low="#f0f0f0", mid="#ffffff", high="#2166ac",
                           midpoint=0, name="z-score") +
      theme_minimal(base_size = 12) +
      labs(x = "Months since infection", y = NULL,
           title = paste0(grp_label, ": ",
                          ifelse(metric_col=="pagerank","PageRank","Strength"), " over time"),
           subtitle = "z-scored within group; darker = higher")
  }

  for (grp in c("Infected", "Not infected")) {
    for (met in c("pagerank", "strength")) {
      p_heat <- make_metric_heatmap(nodes_df, met, grp)
      fname <- paste0("heatmap_", gsub(" ","",tolower(grp)), "_", met, "_timebin.png")
      ggsave(file.path(out_dir, "figures", fname), p_heat, width=7.5, height=8, dpi=300)
    }
  }

  # -- Per-symptom PageRank trajectory (infected group) ---------------
  pr_dir <- file.path(out_dir, "figures", "pagerank_over_time")
  dir.create(pr_dir, recursive = TRUE, showWarnings = FALSE)

  nodes_inf_traj <- nodes_df %>%
    filter(group == "Infected") %>%
    group_by(symptom, timebin) %>%
    summarise(pagerank = mean(pagerank, na.rm = TRUE), .groups = "drop")

  ymax_traj <- max(nodes_inf_traj$pagerank, na.rm = TRUE)

  for (sym in sort(unique(nodes_inf_traj$symptom))) {
    dat  <- nodes_inf_traj %>% filter(symptom == sym)
    p    <- ggplot(dat, aes(x = timebin, y = pagerank, group = 1)) +
      geom_line() + geom_point(size = 2) +
      coord_cartesian(ylim = c(0, ymax_traj)) +
      theme_minimal(base_size = 13) +
      labs(title = paste0(sym, ": PageRank over time (infected)"),
           x = "Months since infection", y = "Weighted PageRank")
    fname <- paste0("pr_infected_", gsub("[^A-Za-z0-9_-]+", "_", sym), ".png")
    ggsave(file.path(pr_dir, fname), p, width = 7, height = 4.5, dpi = 300)
  }
  message("Saved ", n_distinct(nodes_inf_traj$symptom), " per-symptom PageRank plots")

  # Faceted overview
  p_facet <- ggplot(nodes_inf_traj, aes(x = timebin, y = pagerank, group = 1)) +
    geom_line() + geom_point(size = 1.2) +
    facet_wrap(~ symptom, scales = "free_y", ncol = 5) +
    theme_minimal(base_size = 10) +
    labs(x = "Time bin", y = "Weighted PageRank",
         title = "PageRank trajectories (infected group)")
  ggsave(file.path(out_dir, "figures", "pagerank_all_symptoms_facet.png"),
         p_facet, width = 14, height = 10, dpi = 300)

} else {
  message("Skipping time-bin section: time_col not found in data.")
}

########################################################################
# 4E. BOOTSTRAP & PERMUTATION INFERENCE
########################################################################
# Uses build_adj / fit_two_group_networks (Section 4A framework).
# To run with IPW weights, pass weight_col = "ipw" and use df_ipw.
# Here we run on the unweighted baseline data; change `use_weight_col`
# and `data_for_inference` to switch to IPW-weighted inference.
########################################################################

# Toggle: set use_weight_col <- "ipw" and data_for_inference <- df_ipw
#         to run bootstrap on IPW-weighted networks
use_weight_col    <- NULL      # NULL = unweighted baseline
data_for_inference <- df_clean  # or df_ipw for IPW-weighted

# -- Bootstrap --------------------------------------------------------
bootstrap_network_diff <- function(data, symptom_cols, group_col = "infected",
                                   weight_col = NULL, normalize = "none",
                                   B = 500, seed = 123) {
  set.seed(seed)
  dat0 <- data %>% filter(.data[[group_col]] == 0)
  dat1 <- data %>% filter(.data[[group_col]] == 1)
  n0   <- nrow(dat0);  n1 <- nrow(dat1)
  res  <- vector("list", B)
  for (b in seq_len(B)) {
    samp0 <- dat0 %>% slice_sample(n = n0, replace = TRUE)
    samp1 <- dat1 %>% slice_sample(n = n1, replace = TRUE)
    fitb  <- tryCatch(
      fit_two_group_networks(bind_rows(samp0, samp1), symptom_cols,
                             group_col, weight_col, normalize),
      error = function(e) NULL
    )
    if (!is.null(fitb)) {
      res[[b]] <- fitb$node %>%
        select(symptom, delta_strength, delta_pagerank) %>%
        mutate(b = b,
               delta_density       = fitb$global$delta_density,
               delta_modularity    = fitb$global$delta_modularity,
               delta_assortativity = fitb$global$delta_assortativity)
    }
  }
  bind_rows(res)
}

summarise_bootstrap <- function(boot_res) {
  node_summary <- boot_res %>%
    group_by(symptom) %>%
    summarise(
      ds_mean     = mean(delta_strength,  na.rm = TRUE),
      ds_lcl      = quantile(delta_strength,  0.025, na.rm = TRUE),
      ds_ucl      = quantile(delta_strength,  0.975, na.rm = TRUE),
      dpr_mean    = mean(delta_pagerank,  na.rm = TRUE),
      dpr_lcl     = quantile(delta_pagerank,  0.025, na.rm = TRUE),
      dpr_ucl     = quantile(delta_pagerank,  0.975, na.rm = TRUE),
      p_sign_ds   = 2 * pmin(mean(delta_strength <= 0, na.rm=TRUE),
                              mean(delta_strength >= 0, na.rm=TRUE)),
      p_sign_dpr  = 2 * pmin(mean(delta_pagerank <= 0, na.rm=TRUE),
                              mean(delta_pagerank >= 0, na.rm=TRUE)),
      prop_pos_ds  = mean(delta_strength > 0, na.rm = TRUE),
      prop_pos_dpr = mean(delta_pagerank > 0, na.rm = TRUE),
      .groups = "drop"
    )
  global_summary <- boot_res %>%
    summarise(
      ddens_mean   = mean(delta_density,       na.rm = TRUE),
      ddens_lcl    = quantile(delta_density,    0.025, na.rm = TRUE),
      ddens_ucl    = quantile(delta_density,    0.975, na.rm = TRUE),
      dmod_mean    = mean(delta_modularity,     na.rm = TRUE),
      dmod_lcl     = quantile(delta_modularity, 0.025, na.rm = TRUE),
      dmod_ucl     = quantile(delta_modularity, 0.975, na.rm = TRUE),
      dassort_mean = mean(delta_assortativity,  na.rm = TRUE),
      dassort_lcl  = quantile(delta_assortativity, 0.025, na.rm = TRUE),
      dassort_ucl  = quantile(delta_assortativity, 0.975, na.rm = TRUE)
    )
  list(node = node_summary, global = global_summary)
}

# -- Permutation test -------------------------------------------------
permute_network_diff <- function(data, symptom_cols, group_col = "infected",
                                 weight_col = NULL, normalize = "none",
                                 B = 1000, seed = 456) {
  set.seed(seed)
  res <- vector("list", B)
  for (b in seq_len(B)) {
    datb <- data
    datb[[group_col]] <- sample(datb[[group_col]])
    fitb <- tryCatch(
      fit_two_group_networks(datb, symptom_cols, group_col, weight_col, normalize),
      error = function(e) NULL
    )
    if (!is.null(fitb)) {
      res[[b]] <- fitb$node %>%
        select(symptom, delta_strength, delta_pagerank) %>%
        mutate(b = b,
               delta_density       = fitb$global$delta_density,
               delta_modularity    = fitb$global$delta_modularity,
               delta_assortativity = fitb$global$delta_assortativity)
    }
  }
  bind_rows(res)
}

summarise_permutation <- function(obs_fit, perm_res) {
  obs_node <- obs_fit$node
  node_p <- perm_res %>%
    group_by(symptom) %>%
    summarise(null_ds  = list(delta_strength),
              null_dpr = list(delta_pagerank), .groups = "drop") %>%
    left_join(obs_node %>% select(symptom, delta_strength, delta_pagerank), by = "symptom") %>%
    rowwise() %>%
    mutate(
      p_perm_ds  = mean(abs(unlist(null_ds))  >= abs(delta_strength)),
      p_perm_dpr = mean(abs(unlist(null_dpr)) >= abs(delta_pagerank))
    ) %>%
    ungroup() %>%
    mutate(q_perm_ds  = p.adjust(p_perm_ds,  method = "BH"),
           q_perm_dpr = p.adjust(p_perm_dpr, method = "BH")) %>%
    select(symptom, p_perm_ds, p_perm_dpr, q_perm_ds, q_perm_dpr)
  global_p <- tibble(
    p_delta_density =
      mean(abs(perm_res$delta_density)       >= abs(obs_fit$global$delta_density)),
    p_delta_modularity =
      mean(abs(perm_res$delta_modularity)    >= abs(obs_fit$global$delta_modularity)),
    p_delta_assortativity =
      mean(abs(perm_res$delta_assortativity) >= abs(obs_fit$global$delta_assortativity))
  )
  list(node = node_p, global = global_p)
}

# ==========  Run bootstrap & permutation  ============================
message("\nRunning bootstrap (B = ", n_boot, ") ...")
boot_res <- bootstrap_network_diff(
  data         = data_for_inference,
  symptom_cols = symptom_cols,
  group_col    = infect_col,
  weight_col   = use_weight_col,
  normalize    = "none",
  B            = n_boot,
  seed         = set_seed_boot
)
boot_sum <- summarise_bootstrap(boot_res)

message("Running permutation (B = ", n_perm, ") ...")
perm_res <- permute_network_diff(
  data         = data_for_inference,
  symptom_cols = symptom_cols,
  group_col    = infect_col,
  weight_col   = use_weight_col,
  normalize    = "none",
  B            = n_perm,
  seed         = set_seed_perm
)
obs_fit_for_perm <- fit_two_group_networks(
  data         = data_for_inference,
  symptom_cols = symptom_cols,
  group_col    = infect_col,
  weight_col   = use_weight_col,
  normalize    = "none"
)
perm_sum <- summarise_permutation(obs_fit_for_perm, perm_res)

# -- Final combined results table -------------------------------------
final_node_results <- obs_fit_for_perm$node %>%
  left_join(boot_sum$node, by = "symptom") %>%
  left_join(perm_sum$node, by = "symptom") %>%
  mutate(
    z_strength   = as.numeric(scale(delta_strength)),
    z_pagerank   = as.numeric(scale(delta_pagerank)),
    change_score = abs(z_strength) + abs(z_pagerank)
  ) %>%
  arrange(desc(delta_pagerank))

message("\n=== FINAL NODE RESULTS (top 15 by delta PageRank) ===")
print(head(final_node_results %>%
           select(symptom_clean, delta_strength, delta_pagerank,
                  ds_lcl, ds_ucl, dpr_lcl, dpr_ucl,
                  p_sign_dpr, q_perm_dpr), 15))

final_global_results <- bind_cols(
  obs_fit_for_perm$global, boot_sum$global, perm_sum$global
)
message("\n=== FINAL GLOBAL RESULTS ===")
print(as.data.frame(t(final_global_results)))

# Save results
write.csv(final_node_results,
          file.path(out_dir, "tables", "bootstrap_node_results.csv"),
          row.names = FALSE)
write.csv(final_global_results,
          file.path(out_dir, "tables", "bootstrap_global_results.csv"),
          row.names = FALSE)
saveRDS(boot_res, file.path(out_dir, "data", "bootstrap_raw.rds"))
saveRDS(perm_res, file.path(out_dir, "data", "permutation_raw.rds"))

# -- Forest plot: bootstrap CIs for delta PageRank --------------------
p_forest <- final_node_results %>%
  mutate(symptom_clean = factor(symptom_clean, levels = rev(symptom_clean)),
         sig = ifelse(q_perm_dpr < 0.05, "FDR < 0.05", "ns")) %>%
  ggplot(aes(x = dpr_mean, y = symptom_clean, color = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = dpr_lcl, xmax = dpr_ucl), height = 0.3, alpha = 0.7) +
  geom_point(size = 2) +
  scale_color_manual(values = c("FDR < 0.05" = "#d73027", "ns" = "grey50")) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  labs(x = "Delta PageRank (95% bootstrap CI)",
       y = NULL, color = NULL,
       title = "Bootstrap: Infection-associated change in symptom PageRank")


ggsave(file.path(out_dir, "figures", "bootstrap_forest_pagerank.png"),
       p_forest, width = 9, height = 8, dpi = 300)
print(p_forest)

message("\nDone. All outputs saved to: ", normalizePath(out_dir))


###

# -- Network plot helper ----------------------------------------------
plot_two_group_network <- function(g_left, g_right,
                                   left_title = "A. Uninfected",
                                   right_title = "B. Infected",
                                   file = NULL,
                                   width = 14, height = 6.5, res = 300,
                                   seed = 123) {
  set.seed(seed)
  
  # make sure both graphs contain all symptom nodes
  all_nodes <- union(V(g_left)$name, V(g_right)$name)
  
  add_missing_vertices <- function(g, all_nodes) {
    missing_nodes <- setdiff(all_nodes, V(g)$name)
    if (length(missing_nodes) > 0) {
      g <- add_vertices(g, nv = length(missing_nodes), name = missing_nodes)
    }
    g
  }
  
  g_left  <- add_missing_vertices(g_left,  all_nodes)
  g_right <- add_missing_vertices(g_right, all_nodes)
  
  # reorder vertices identically
  g_left  <- permute(g_left,  match(all_nodes, V(g_left)$name))
  g_right <- permute(g_right, match(all_nodes, V(g_right)$name))
  
  # use a common layout based on combined graph for comparability
  g_comb <- igraph::union(g_left, g_right, byname = TRUE)
  lay <- layout_with_fr(g_comb, weights = E(g_comb)$weight)
  
  # common community colors
  get_membership <- function(g) {
    if (ecount(g) == 0) return(rep(1L, vcount(g)))
    membership(cluster_louvain(g, weights = E(g)$weight))
  }
  
  mem_left  <- get_membership(g_left)
  mem_right <- get_membership(g_right)
  
  # consistent palette across panels
  all_comms <- sort(unique(c(mem_left, mem_right)))
  pal <- setNames(
    grDevices::hcl.colors(length(all_comms), palette = "Dark 3"),
    all_comms
  )
  
  # node sizes based on strength
  get_vertex_size <- function(g, min_size = 4, max_size = 18) {
    if (ecount(g) == 0) return(rep(min_size, vcount(g)))
    s <- strength(g, weights = E(g)$weight)
    if (max(s, na.rm = TRUE) == min(s, na.rm = TRUE)) {
      return(rep((min_size + max_size) / 2, length(s)))
    }
    scales::rescale(s, to = c(min_size, max_size))
  }
  
  # edge widths based on weight
  get_edge_width <- function(g, min_w = 0.4, max_w = 4) {
    if (ecount(g) == 0) return(numeric(0))
    w <- E(g)$weight
    if (max(w, na.rm = TRUE) == min(w, na.rm = TRUE)) {
      return(rep((min_w + max_w) / 2, length(w)))
    }
    scales::rescale(w, to = c(min_w, max_w))
  }
  
  vsize_left  <- get_vertex_size(g_left)
  vsize_right <- get_vertex_size(g_right)
  ewidth_left  <- get_edge_width(g_left)
  ewidth_right <- get_edge_width(g_right)
  
  vcol_left  <- pal[as.character(mem_left)]
  vcol_right <- pal[as.character(mem_right)]
  
  # save to file if requested
  if (!is.null(file)) {
    png(file, width = width, height = height, units = "in", res = res)
    on.exit(dev.off(), add = TRUE)
  }
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  par(mfrow = c(1, 2), mar = c(1, 1, 3, 1), xpd = TRUE)
  
  # left panel
  plot(
    g_left,
    layout = lay,
    vertex.label = V(g_left)$name,
    vertex.label.cex = 0.55,
    vertex.label.family = "sans",
    vertex.label.color = "black",
    vertex.size = vsize_left,
    vertex.color = vcol_left,
    vertex.frame.color = NA,
    edge.width = ewidth_left,
    edge.color = grDevices::adjustcolor("grey40", alpha.f = 0.35),
    main = left_title
  )
  
  legend("topright",
         legend = c("strength", "", "", ""),
         bty = "n", cex = 0.9)
  
  if (ecount(g_left) > 0) {
    s_left <- strength(g_left, weights = E(g_left)$weight)
    s_breaks <- round(pretty(range(s_left), n = 4))
    s_breaks <- unique(s_breaks[s_breaks > 0])
    if (length(s_breaks) > 0) {
      legend("right",
             legend = s_breaks,
             pch = 16,
             pt.cex = scales::rescale(s_breaks, to = c(0.8, 2.5)),
             bty = "n", title = "strength", inset = c(-0.02, 0))
    }
    w_left <- E(g_left)$weight
    w_breaks <- round(pretty(range(w_left), n = 4))
    w_breaks <- unique(w_breaks[w_breaks > 0])
    if (length(w_breaks) > 0) {
      legend("bottomright",
             legend = w_breaks,
             lwd = scales::rescale(w_breaks, to = c(0.4, 4)),
             col = "grey50",
             bty = "n", title = "weight", inset = c(-0.02, 0))
    }
  }
  
  legend("bottomleft",
         legend = paste("community", sort(unique(mem_left))),
         pch = 16,
         col = pal[as.character(sort(unique(mem_left)))],
         bty = "n",
         cex = 0.8)
  
  # right panel
  plot(
    g_right,
    layout = lay,
    vertex.label = V(g_right)$name,
    vertex.label.cex = 0.55,
    vertex.label.family = "sans",
    vertex.label.color = "black",
    vertex.size = vsize_right,
    vertex.color = vcol_right,
    vertex.frame.color = NA,
    edge.width = ewidth_right,
    edge.color = grDevices::adjustcolor("grey40", alpha.f = 0.35),
    main = right_title
  )
  
  if (ecount(g_right) > 0) {
    s_right <- strength(g_right, weights = E(g_right)$weight)
    s_breaks <- round(pretty(range(s_right), n = 4))
    s_breaks <- unique(s_breaks[s_breaks > 0])
    if (length(s_breaks) > 0) {
      legend("right",
             legend = s_breaks,
             pch = 16,
             pt.cex = scales::rescale(s_breaks, to = c(0.8, 2.5)),
             bty = "n", title = "strength", inset = c(-0.02, 0))
    }
    w_right <- E(g_right)$weight
    w_breaks <- round(pretty(range(w_right), n = 4))
    w_breaks <- unique(w_breaks[w_breaks > 0])
    if (length(w_breaks) > 0) {
      legend("bottomright",
             legend = w_breaks,
             lwd = scales::rescale(w_breaks, to = c(0.4, 4)),
             col = "grey50",
             bty = "n", title = "weight", inset = c(-0.02, 0))
    }
  }
  
  legend("bottomleft",
         legend = paste("community", sort(unique(mem_right))),
         pch = 16,
         col = pal[as.character(sort(unique(mem_right)))],
         bty = "n",
         cex = 0.8)
}

library(scales)
g_base_inf  <- baseline_fit$m1$graph
g_base_ctrl <- baseline_fit$m0$graph

# -- Baseline two-panel network plot ----------------------------------
plot_two_group_network(
  g_left  = g_base_ctrl,
  g_right = g_base_inf,
  left_title  = "A. Uninfected",
  right_title = "B. Infected",
  file = file.path(out_dir, "figures", "baseline_network_two_panel.png"),
  width = 14,
  height = 6.5,
  res = 300,
  seed = 123
)

g_inf_ipw  <- net_inf_ipw$graph
g_ctrl_ipw <- net_ctrl_ipw$graph

# -- IPW two-panel network plot ---------------------------------------
plot_two_group_network(
  g_left  = g_ctrl_ipw,
  g_right = g_inf_ipw,
  left_title  = "A. Uninfected (IPW)",
  right_title = "B. Infected (IPW)",
  file = file.path(out_dir, "figures", "ipw_network_two_panel.png"),
  width = 14,
  height = 6.5,
  res = 300,
  seed = 123
)


########################


plot_timebin_symptom_influence <- function(nodes_df,
                                           top_n = 12,
                                           file = NULL,
                                           width = 16,
                                           height = 9,
                                           dpi = 300) {
  
  nodes_plot <- nodes_df %>%
    mutate(
      symptom = clean_symptom_names(symptom),
      group = factor(group, levels = c("Infected", "Not infected")),
      timebin = factor(as.character(timebin), levels = c("0-3", "3-6", "6-12", "12-24"))
    ) %>%
    group_by(group, timebin) %>%
    arrange(desc(pagerank), .by_group = TRUE) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(
      panel = paste0(group, ": ", timebin, " months"),
      panel = factor(
        panel,
        levels = c(
          "Infected: 0-3 months", "Infected: 3-6 months",
          "Infected: 6-12 months", "Infected: 12-24 months",
          "Not infected: 0-3 months", "Not infected: 3-6 months",
          "Not infected: 6-12 months", "Not infected: 12-24 months"
        )
      )
    )
  
  # reorder within panel
  nodes_plot <- nodes_plot %>%
    group_by(panel) %>%
    arrange(pagerank, .by_group = TRUE) %>%
    mutate(symptom_panel = factor(
      paste(panel, symptom, sep = "___"),
      levels = paste(panel, symptom, sep = "___")
    )) %>%
    ungroup()
  
  # stable community colors
  comm_vals <- sort(unique(as.character(nodes_plot$community)))
  comm_cols <- setNames(
    grDevices::hcl.colors(length(comm_vals), palette = "Dark 3"),
    comm_vals
  )
  
  p <- ggplot(nodes_plot, aes(x = pagerank, y = symptom_panel, color = community)) +
    geom_segment(aes(x = 0, xend = pagerank, yend = symptom_panel), linewidth = 0.8) +
    geom_point(aes(size = strength)) +
    facet_wrap(~ panel, ncol = 4, scales = "free_y") +
    scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
    scale_color_manual(values = comm_cols) +
    labs(
      title = "Time-bin-specific symptom influence by infection status",
      x = "Weighted PageRank (global influence)",
      y = NULL,
      color = "Community",
      size = "Strength"
    ) +
    theme_bw()
  
  if (!is.null(file)) {
    ggsave(file, p, width = width, height = height, dpi = dpi)
  }
  
  return(p)
}

p <- plot_timebin_symptom_influence(nodes_df)

print(p)

plot_timebin_symptom_influence(
  nodes_df,
  file = file.path(out_dir, "figures", "timebin_plot.png")
)

########

# -- Forest-style time-bin symptom influence plot ---------------------
plot_timebin_forest_full <- function(nodes_df,
                                     file = NULL,
                                     width = 16,
                                     height = 10,
                                     dpi = 300) {
  
  nodes_plot <- nodes_df %>%
    mutate(
      symptom = clean_symptom_names(symptom),
      group = factor(group, levels = c("Infected", "Not infected")),
      timebin = factor(as.character(timebin), levels = c("0-3", "3-6", "6-12", "12-24")),
      panel = paste0(group, ": ", timebin, " months"),
      panel = factor(
        panel,
        levels = c(
          "Infected: 0-3 months", "Infected: 3-6 months",
          "Infected: 6-12 months", "Infected: 12-24 months",
          "Not infected: 0-3 months", "Not infected: 3-6 months",
          "Not infected: 6-12 months", "Not infected: 12-24 months"
        )
      )
    ) %>%
    group_by(panel) %>%
    arrange(pagerank, .by_group = TRUE) %>%
    mutate(
      symptom_panel = factor(
        paste(panel, symptom, sep = "___"),
        levels = paste(panel, symptom, sep = "___")
      )
    ) %>%
    ungroup()
  
  comm_vals <- sort(unique(as.character(nodes_plot$community)))
  comm_cols <- setNames(
    grDevices::hcl.colors(length(comm_vals), palette = "Dark 3"),
    comm_vals
  )
  
  p <- ggplot(nodes_plot, aes(x = pagerank, y = symptom_panel)) +
    geom_segment(
      aes(x = 0, xend = pagerank, yend = symptom_panel, color = community),
      linewidth = 0.7, alpha = 0.8
    ) +
    geom_point(aes(size = strength, color = community), alpha = 0.95) +
    facet_wrap(~ panel, ncol = 4, scales = "free_y") +
    scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
    scale_color_manual(values = comm_cols) +
    labs(
      title = "Supplementary Figure. Time-bin-specific symptom influence by infection status",
      subtitle = "Top row: Infected; Bottom row: Not infected. Line length indicates weighted PageRank; point size indicates strength.",
      x = "Weighted PageRank (global influence)",
      y = NULL,
      color = "Community",
      size = "Strength"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    )
  
  if (!is.null(file)) {
    ggsave(file, p, width = width, height = height, dpi = dpi)
  }
  
  return(p)
}
print(metrics_df)

# -- Forest-style full symptom plot across time bins ------------------
p_timebin_forest_full <- plot_timebin_forest_full(
  nodes_df = nodes_df,
  file = file.path(out_dir, "figures", "timebin_symptom_influence_forest_full.png"),
  width = 16,
  height = 10,
  dpi = 300
)

print(p_timebin_forest_full)

########


# -- Overall two-panel lollipop plot by group -------------------------
plot_overall_group_lollipop <- function(g_ctrl, g_inf,
                                        top_n = NULL,
                                        file = NULL,
                                        width = 12,
                                        height = 7,
                                        dpi = 300) {
  
  ctrl_df <- calc_node_metrics(g_ctrl, "Control") %>%
    mutate(group = "Panel A. Control")
  
  inf_df <- calc_node_metrics(g_inf, "Infected") %>%
    mutate(group = "Panel B. Infected")
  
  plot_df <- bind_rows(ctrl_df, inf_df) %>%
    mutate(symptom = clean_symptom_names(symptom))
  
  if (!is.null(top_n)) {
    plot_df <- plot_df %>%
      group_by(group) %>%
      arrange(desc(pagerank), .by_group = TRUE) %>%
      slice_head(n = top_n) %>%
      ungroup()
  }
  
  plot_df <- plot_df %>%
    group_by(group) %>%
    arrange(pagerank, .by_group = TRUE) %>%
    mutate(symptom_panel = factor(
      paste(group, symptom, sep = "___"),
      levels = paste(group, symptom, sep = "___")
    )) %>%
    ungroup()
  
  comm_vals <- sort(unique(as.character(plot_df$community)))
  comm_cols <- setNames(
    grDevices::hcl.colors(length(comm_vals), palette = "Dark 3"),
    comm_vals
  )
  
  p <- ggplot(plot_df, aes(x = pagerank, y = symptom_panel, color = factor(community))) +
    geom_segment(
      aes(x = 0, xend = pagerank, yend = symptom_panel),
      linewidth = 0.8, alpha = 0.9
    ) +
    geom_point(aes(size = strength), alpha = 0.95) +
    facet_wrap(~ group, ncol = 2, scales = "free_y") +
    scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
    scale_color_manual(values = comm_cols, name = "Community") +
    labs(
      x = "Weighted PageRank (global influence)",
      y = NULL,
      size = "Strength"
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      axis.text.y = element_text(size = 9)
    )
  
  if (!is.null(file)) {
    ggsave(file, p, width = width, height = height, dpi = dpi)
  }
  
  return(p)
}

g_base_inf  <- baseline_fit$m1$graph
g_base_ctrl <- baseline_fit$m0$graph

# -- Overall group-level lollipop plot --------------------------------
p_overall_group_lollipop <- plot_overall_group_lollipop(
  g_ctrl = g_base_ctrl,
  g_inf  = g_base_inf,
  top_n  = NULL,  # use NULL for full symptom list; use 20 if you want top 20 only
  file   = file.path(out_dir, "figures", "overall_group_lollipop.png"),
  width  = 12,
  height = 7,
  dpi    = 300
)

print(p_overall_group_lollipop)

ctrl_df <- calc_node_metrics(g_ctrl, "Control") %>%
  mutate(group = "Panel A. Control")

ctrl_df <- calc_node_metrics(g_ctrl, "Control") %>%
  mutate(group = "Panel A. Control")

###

# -- Domain mapping for plot labels -----------------------------------
symptom_domain_plot <- function(x) {
  dplyr::case_when(
    x %in% c("rhinorrhea_throat", "sore_throat", "anosmia", "ageusia", "HEARING_DIS") ~ "Sensory/ENT",
    x %in% c("sputum", "cough", "dyspnea") ~ "Respiratory",
    x %in% c("SLE_DIS", "DIZZ", "LOW_CONCENTRATE", "cognition_disorder", "OTH_NEURO_SX") ~ "Cognitive/Neuro",
    x %in% c("fatigue", "arthritis", "headache", "weakness", "myalgia", "PAIN", "chill", "fever") ~ "Systemic/MSK",
    x %in% c("palpitation", "LOW_APPE", "CONSTIPATION", "ALOPECIA", "DIARRHEA", "chest_pain", "rash", "OTH_GI", "NAUSEA") ~ "Other/GI",
    TRUE ~ "Other"
  )
}

# -- Plot symptom trajectories by domain ------------------------------
plot_pagerank_trajectories_by_domain <- function(nodes_df,
                                                 out_dir_fig,
                                                 metric_col = "pagerank",
                                                 ncol = 4,
                                                 width_per_col = 3.2,
                                                 height_per_row = 2.3,
                                                 base_size = 11) {
  
  df_plot <- nodes_df %>%
    dplyr::mutate(
      symptom_raw = symptom,
      symptom = clean_symptom_names(symptom),
      domain = symptom_domain_plot(symptom_raw),
      timebin = factor(as.character(timebin), levels = c("0-3", "3-6", "6-12", "12-24")),
      group = factor(group, levels = c("Infected", "Not infected"))
    ) %>%
    dplyr::group_by(domain, symptom, group, timebin) %>%
    dplyr::summarise(value = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop")
  
  # order symptoms within each domain by average infected value
  symptom_order_df <- df_plot %>%
    dplyr::filter(group == "Infected") %>%
    dplyr::group_by(domain, symptom) %>%
    dplyr::summarise(avg_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(domain, dplyr::desc(avg_value))
  
  domains <- unique(symptom_order_df$domain)
  
  for (dm in domains) {
    dm_syms <- symptom_order_df %>%
      dplyr::filter(domain == dm) %>%
      dplyr::pull(symptom)
    
    df_dm <- df_plot %>%
      dplyr::filter(domain == dm) %>%
      dplyr::mutate(symptom = factor(symptom, levels = dm_syms))
    
    n_panels <- length(unique(df_dm$symptom))
    n_rows <- ceiling(n_panels / ncol)
    
    p <- ggplot(df_dm, aes(x = timebin, y = value, group = group, color = group, linetype = group)) +
      geom_line(linewidth = 0.6) +
      geom_point(size = 1.6) +
      facet_wrap(~ symptom, ncol = ncol, scales = "fixed") +
      labs(
        title = paste0(dm, ": temporal ", ifelse(metric_col == "pagerank", "PageRank", metric_col), " trajectories"),
        x = "Months since infection (time bin)",
        y = ifelse(metric_col == "pagerank", "Weighted PageRank", metric_col),
        color = NULL,
        linetype = NULL
      ) +
      theme_bw(base_size = base_size) +
      theme(
        strip.text = element_text(face = "bold", size = base_size - 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(face = "bold")
      )
    
    outfile <- file.path(
      out_dir_fig,
      paste0("traj_", gsub("[^A-Za-z0-9]+", "_", tolower(dm)), "_", metric_col, ".png")
    )
    
    ggsave(
      filename = outfile,
      plot = p,
      width = width_per_col * ncol,
      height = height_per_row * n_rows,
      dpi = 300
    )
    
    print(p)
  }
}
write.csv(metrics_df, file.path(out_dir, "tables", "metrics_timebin_by_group.csv"), row.names = FALSE)

# -- Domain-grouped temporal trajectory plots -------------------------
traj_dir <- file.path(out_dir, "figures", "trajectories_by_domain")
dir.create(traj_dir, recursive = TRUE, showWarnings = FALSE)

plot_pagerank_trajectories_by_domain(
  nodes_df = nodes_df,
  out_dir_fig = traj_dir,
  metric_col = "pagerank",
  ncol = 4,
  width_per_col = 3.2,
  height_per_row = 2.3,
  base_size = 11
)

plot_pagerank_trajectories_by_domain(
  nodes_df = nodes_df,
  out_dir_fig = traj_dir,
  metric_col = "strength",
  ncol = 4,
  width_per_col = 3.2,
  height_per_row = 2.3,
  base_size = 11
)

plot_pagerank_trajectories_by_domain(
  nodes_df = nodes_df,
  out_dir_fig = traj_dir,
  metric_col = "strength",
  ncol = 4,
  width_per_col = 3.2,
  height_per_row = 2.3,
  base_size = 11
)

# -- Combined figure across all domains -------------------------------
p_traj_all_domains <- plot_pagerank_trajectories_all_domains(
  nodes_df = nodes_df,
  metric_col = "pagerank",
  file = file.path(out_dir, "figures", "pagerank_trajectories_all_domains.png"),
  base_size = 10,
  width = 18,
  height = 11,
  free_y = FALSE
)

print(p_traj_all_domains)

p_strength_all_domains <- plot_pagerank_trajectories_all_domains(
  nodes_df = nodes_df,
  metric_col = "strength",
  file = file.path(out_dir, "figures", "strength_trajectories_all_domains.png"),
  base_size = 10,
  width = 18,
  height = 11,
  free_y = TRUE
)

print(p_strength_all_domains)

###
plot_pagerank_trajectories_domain_matrix <- function(nodes_df,
                                                     metric_col = "pagerank",
                                                     max_rows_per_domain = 9,
                                                     domain_levels = c("Sensory/ENT",
                                                                       "Respiratory",
                                                                       "Cognitive/Neuro",
                                                                       "Systemic/MSK",
                                                                       "Other/GI"),
                                                     file = NULL,
                                                     width = 18,
                                                     height = 14,
                                                     base_size = 10) {
  
  # Prepare data
  df_plot <- nodes_df %>%
    dplyr::mutate(
      symptom_raw = symptom,
      symptom = clean_symptom_names(symptom),
      domain = symptom_domain_plot(symptom_raw),
      timebin = factor(as.character(timebin), levels = c("0-3", "3-6", "6-12", "12-24")),
      group = factor(group, levels = c("Infected", "Not infected"))
    ) %>%
    dplyr::group_by(domain, symptom, group, timebin) %>%
    dplyr::summarise(value = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(domain = factor(domain, levels = domain_levels))
  
  # Order symptoms within each domain by average infected value
  symptom_order <- df_plot %>%
    dplyr::filter(group == "Infected") %>%
    dplyr::group_by(domain, symptom) %>%
    dplyr::summarise(avg_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(domain, dplyr::desc(avg_value)) %>%
    dplyr::group_by(domain) %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # Keep only up to max_rows_per_domain symptoms per domain
  symptom_order <- symptom_order %>%
    dplyr::filter(row_id <= max_rows_per_domain)
  
  # Attach row positions
  df_plot2 <- df_plot %>%
    dplyr::inner_join(
      symptom_order %>% dplyr::select(domain, symptom, row_id),
      by = c("domain", "symptom")
    )
  
  # Build empty slots so each domain column has the same number of rows
  blank_slots <- expand.grid(
    domain = factor(domain_levels, levels = domain_levels),
    row_id = seq_len(max_rows_per_domain)
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::left_join(
      symptom_order %>% dplyr::select(domain, row_id, symptom),
      by = c("domain", "row_id")
    ) %>%
    dplyr::mutate(
      panel_label = dplyr::if_else(is.na(symptom), "", symptom),
      panel_id = paste(domain, row_id, sep = "___")
    )
  
  # Plotting data
  df_plot2 <- df_plot2 %>%
    dplyr::mutate(
      panel_id = paste(domain, row_id, sep = "___")
    )
  
  blank_slots <- blank_slots %>%
    dplyr::mutate(
      domain = factor(domain, levels = domain_levels),
      row_id = factor(row_id, levels = seq_len(max_rows_per_domain)),
      panel_id = factor(panel_id, levels = blank_slots$panel_id)
    )
  
  df_plot2 <- df_plot2 %>%
    dplyr::mutate(
      domain = factor(domain, levels = domain_levels),
      row_id = factor(row_id, levels = seq_len(max_rows_per_domain)),
      panel_id = factor(panel_id, levels = levels(blank_slots$panel_id))
    )
  
  # Dummy layer so blank facets are retained
  dummy_df <- blank_slots %>%
    dplyr::mutate(
      timebin = factor("0-3", levels = c("0-3", "3-6", "6-12", "12-24")),
      value = NA_real_,
      group = factor("Infected", levels = c("Infected", "Not infected"))
    )
  
  p <- ggplot() +
    geom_blank(
      data = dummy_df,
      aes(x = timebin, y = value)
    ) +
    geom_line(
      data = df_plot2,
      aes(x = timebin, y = value, group = group, color = group, linetype = group),
      linewidth = 0.55
    ) +
    geom_point(
      data = df_plot2,
      aes(x = timebin, y = value, color = group),
      size = 1.4
    ) +
    facet_grid(
      rows = vars(row_id),
      cols = vars(domain),
      scales = "fixed",
      switch = "x"
    ) +
    labs(
      title = paste0("Temporal ", ifelse(metric_col == "pagerank", "PageRank", metric_col),
                     " trajectories grouped by domain"),
      x = "Months since infection (time bin)",
      y = ifelse(metric_col == "pagerank", "Weighted PageRank", metric_col),
      color = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = base_size) +
    theme(
      strip.text.x = element_text(face = "bold", size = base_size),
      strip.text.y = element_blank(),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(0.5, "lines"),
      plot.title = element_text(face = "bold")
    )
  
  # Add symptom labels manually inside each panel strip area via labeller workaround
  # Simpler approach: use geom_text with fixed x/y for panel titles
  label_df <- blank_slots %>%
    dplyr::filter(panel_label != "") %>%
    dplyr::mutate(
      x = factor("0-3", levels = c("0-3", "3-6", "6-12", "12-24")),
      y = Inf
    )
  
  p <- p +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = panel_label),
      inherit.aes = FALSE,
      vjust = 1.4,
      hjust = 0,
      size = 2.8,
      fontface = "bold"
    )
  
  if (!is.null(file)) {
    ggsave(file, p, width = width, height = height, dpi = 300)
  }
  
  p
}

plot_pagerank_trajectories_domain_matrix <- function(nodes_df,
                                                     metric_col = "pagerank",
                                                     max_rows_per_domain = 9,
                                                     domain_levels = c("Sensory/ENT",
                                                                       "Respiratory",
                                                                       "Cognitive/Neuro",
                                                                       "Systemic/MSK",
                                                                       "Other/GI"),
                                                     file = NULL,
                                                     width = 18,
                                                     height = 14,
                                                     base_size = 10) {
  
  # Prepare data
  df_plot <- nodes_df %>%
    dplyr::mutate(
      symptom_raw = symptom,
      symptom = clean_symptom_names(symptom),
      domain = symptom_domain_plot(symptom_raw),
      timebin = factor(as.character(timebin), levels = c("0-3", "3-6", "6-12", "12-24")),
      group = factor(group, levels = c("Infected", "Not infected"))
    ) %>%
    dplyr::group_by(domain, symptom, group, timebin) %>%
    dplyr::summarise(value = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(domain = factor(domain, levels = domain_levels))
  
  # Order symptoms within each domain by average infected value
  symptom_order <- df_plot %>%
    dplyr::filter(group == "Infected") %>%
    dplyr::group_by(domain, symptom) %>%
    dplyr::summarise(avg_value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(domain, dplyr::desc(avg_value)) %>%
    dplyr::group_by(domain) %>%
    dplyr::mutate(row_id = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # Keep only up to max_rows_per_domain symptoms per domain
  symptom_order <- symptom_order %>%
    dplyr::filter(row_id <= max_rows_per_domain)
  
  # Attach row positions
  df_plot2 <- df_plot %>%
    dplyr::inner_join(
      symptom_order %>% dplyr::select(domain, symptom, row_id),
      by = c("domain", "symptom")
    )
  
  # Build empty slots so each domain column has the same number of rows
  blank_slots <- expand.grid(
    domain = factor(domain_levels, levels = domain_levels),
    row_id = seq_len(max_rows_per_domain)
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::left_join(
      symptom_order %>% dplyr::select(domain, row_id, symptom),
      by = c("domain", "row_id")
    ) %>%
    dplyr::mutate(
      panel_label = dplyr::if_else(is.na(symptom), "", symptom),
      panel_id = paste(domain, row_id, sep = "___")
    )
  
  # Plotting data
  df_plot2 <- df_plot2 %>%
    dplyr::mutate(
      panel_id = paste(domain, row_id, sep = "___")
    )
  
  blank_slots <- blank_slots %>%
    dplyr::mutate(
      domain = factor(domain, levels = domain_levels),
      row_id = factor(row_id, levels = seq_len(max_rows_per_domain)),
      panel_id = factor(panel_id, levels = blank_slots$panel_id)
    )
  
  df_plot2 <- df_plot2 %>%
    dplyr::mutate(
      domain = factor(domain, levels = domain_levels),
      row_id = factor(row_id, levels = seq_len(max_rows_per_domain)),
      panel_id = factor(panel_id, levels = levels(blank_slots$panel_id))
    )
  
  # Dummy layer so blank facets are retained
  dummy_df <- blank_slots %>%
    dplyr::mutate(
      timebin = factor("0-3", levels = c("0-3", "3-6", "6-12", "12-24")),
      value = NA_real_,
      group = factor("Infected", levels = c("Infected", "Not infected"))
    )
  
  p <- ggplot() +
    geom_blank(
      data = dummy_df,
      aes(x = timebin, y = value)
    ) +
    geom_line(
      data = df_plot2,
      aes(x = timebin, y = value, group = group, color = group, linetype = group),
      linewidth = 0.55
    ) +
    geom_point(
      data = df_plot2,
      aes(x = timebin, y = value, color = group),
      size = 1.4
    ) +
    facet_grid(
      rows = vars(row_id),
      cols = vars(domain),
      scales = "fixed",
      switch = "x"
    ) +
    labs(
      title = paste0("Temporal ", ifelse(metric_col == "pagerank", "PageRank", metric_col),
                     " trajectories grouped by domain"),
      x = "Months since infection (time bin)",
      y = ifelse(metric_col == "pagerank", "Weighted PageRank", metric_col),
      color = NULL,
      linetype = NULL
    ) +
    theme_bw(base_size = base_size) +
    theme(
      strip.text.x = element_text(face = "bold", size = base_size),
      strip.text.y = element_blank(),
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      panel.spacing = grid::unit(0.5, "lines"),
      plot.title = element_text(face = "bold")
    )
  
  # Add symptom labels manually inside each panel strip area via labeller workaround
  # Simpler approach: use geom_text with fixed x/y for panel titles
  label_df <- blank_slots %>%
    dplyr::filter(panel_label != "") %>%
    dplyr::mutate(
      x = factor("0-3", levels = c("0-3", "3-6", "6-12", "12-24")),
      y = Inf
    )
  
  p <- p +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = panel_label),
      inherit.aes = FALSE,
      vjust = 1.4,
      hjust = 0,
      size = 2.8,
      fontface = "bold"
    )
  
  if (!is.null(file)) {
    ggsave(file, p, width = width, height = height, dpi = 300)
  }
  
  p
}

p_domain_matrix <- plot_pagerank_trajectories_domain_matrix(
  nodes_df = nodes_df,
  metric_col = "pagerank",
  max_rows_per_domain = 9,
  file = file.path(out_dir, "figures", "pagerank_trajectories_domain_matrix.png"),
  width = 18,
  height = 14,
  base_size = 10
)

print(p_domain_matrix)

p_domain_matrix_strength <- plot_pagerank_trajectories_domain_matrix(
  nodes_df = nodes_df,
  metric_col = "strength",
  max_rows_per_domain = 9,
  file = file.path(out_dir, "figures", "strength_trajectories_domain_matrix.png"),
  width = 18,
  height = 14,
  base_size = 10
)

print(p_domain_matrix_strength)
###
# -- GLOBAL metrics (overall, not time-binned) ------------------------
global_metrics_clean <- global_metrics_table %>%
  mutate(
    infected = ifelse(group == "Infected", 1, 0),
    group_label = ifelse(group == "Infected", "Infected_all", "Uninfected_all")
  ) %>%
  select(
    infected,
    group_label,
    nodes,
    edges,
    weighted_density,
    modularity,
    assortativity_strength
  )

metrics_long <- metrics_df %>%
  mutate(
    group_label = ifelse(infected == 1, "Infected_all", "Uninfected_all")
  ) %>%
  select(
    infected,
    group_label,
    timebin,
    n_units,
    nodes,
    edges,
    weighted_density,
    modularity,
    strength_assortativity
  )


# --- N unit row ---
tab_n_global <- df_ipw %>%
  mutate(group_label = ifelse(.data[[infect_col]] == 1, "Infected_all", "Uninfected_all")) %>%
  group_by(group_label) %>%
  summarise(Global = n(), .groups = "drop")

tab_n <- metrics_long %>%
  select(group_label, timebin, n_units) %>%
  pivot_wider(names_from = timebin, values_from = n_units) %>%
  left_join(tab_n_global, by = "group_label") %>%
  mutate(metric = "N_unit")

# --- Nodes ---
tab_nodes <- metrics_long %>%
  select(group_label, timebin, nodes) %>%
  pivot_wider(names_from = timebin, values_from = nodes) %>%
  left_join(
    global_metrics_clean %>% select(group_label, nodes) %>% rename(Global = nodes),
    by = "group_label"
  ) %>%
  mutate(metric = "Nodes")

# --- Edges ---
tab_edges <- metrics_long %>%
  select(group_label, timebin, edges) %>%
  pivot_wider(names_from = timebin, values_from = edges) %>%
  left_join(
    global_metrics_clean %>% select(group_label, edges) %>% rename(Global = edges),
    by = "group_label"
  ) %>%
  mutate(metric = "Edges")

# --- Density ---
tab_density <- metrics_long %>%
  select(group_label, timebin, weighted_density) %>%
  pivot_wider(names_from = timebin, values_from = weighted_density) %>%
  left_join(
    global_metrics_clean %>% select(group_label, weighted_density) %>%
      rename(Global = weighted_density),
    by = "group_label"
  ) %>%
  mutate(metric = "Weighted density")

# --- Modularity ---
tab_mod <- metrics_long %>%
  select(group_label, timebin, modularity) %>%
  pivot_wider(names_from = timebin, values_from = modularity) %>%
  left_join(
    global_metrics_clean %>% select(group_label, modularity) %>%
      rename(Global = modularity),
    by = "group_label"
  ) %>%
  mutate(metric = "Modularity")

# --- Assortativity ---
tab_assort <- metrics_long %>%
  select(group_label, timebin, strength_assortativity) %>%
  pivot_wider(names_from = timebin, values_from = strength_assortativity) %>%
  left_join(
    global_metrics_clean %>% select(group_label, assortativity_strength) %>%
      rename(Global = assortativity_strength),
    by = "group_label"
  ) %>%
  mutate(metric = "Strength assortativity")

final_table <- bind_rows(
  tab_n,
  tab_nodes,
  tab_edges,
  tab_density,
  tab_mod,
  tab_assort
) %>%
  select(metric, group_label, Global, `0-3`, `3-6`, `6-12`, `12-24`) %>%
  arrange(metric, group_label)

print(final_table)

write.csv(
  final_table,
  file.path(out_dir, "tables", "global_timebin_summary_table.csv"),
  row.names = FALSE
)

