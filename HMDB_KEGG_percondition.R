# =============================================================================
# Metabolic Pathway Enrichment Bar Plots — KEGG & HMDB (separately)
# =============================================================================
# Reads an Excel metabolomics file, averages intensities per condition,
# maps compound names to KEGG and HMDB pathways via live API queries,
# then computes per-pathway mean fold enrichment (no control comparison —
# absolute enrichment within each condition) and plots the top-N pathways
# as bar graphs, coloured by condition, for both databases separately.
#
# CONDITIONS (edit CONDITION_MAP to match your experiment):
#   NC  = Cancer
#   NN, C = Control
#   TC  = Torpor + Cancer
#   TN  = Torpor
#
# REQUIRED PACKAGES:
#   CRAN:         readxl, dplyr, tidyr, ggplot2, stringr, purrr,
#                 scales, httr, jsonlite
#   Bioconductor: KEGGREST, metaboliteIDmapping   (for HMDB)
#
# Install:
#   install.packages(c("readxl","dplyr","tidyr","ggplot2","stringr",
#                      "purrr","scales","httr","jsonlite"))
#   if (!require("BiocManager")) install.packages("BiocManager")
#   BiocManager::install(c("KEGGREST", "metaboliteIDmapping"))
# =============================================================================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(purrr)
  library(scales)
  library(httr)
  library(jsonlite)
})

install.packages(c("readxl","dplyr","tidyr","ggplot2","stringr",
                   "purrr","scales","httr","jsonlite"))
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("KEGGREST", "metaboliteIDmapping"))

# =============================================================================
# ── USER CONFIGURATION ────────────────────────────────────────────────────────
# =============================================================================

INPUT_FILE    <- "RP_Plasma_Pos_Neg_Curated_alldata.xlsx"
SHEET_NAME    <- "Plasma RP PosNeg"
NAME_COL      <- "Name"
SAMPLE_PREFIX <- "Plasma_"

CONDITION_MAP <- list(
  "^(NN|C)$" = "Control",
  "^NC$"      = "Cancer",
  "^TC$"      = "Torpor + Cancer",
  "^TN$"      = "Torpor"
)

TOP_N_PATHWAYS <- 15    # top N pathways per condition per bar chart
MIN_COMPOUNDS  <- 2     # minimum compounds per pathway to include
API_PAUSE_SEC  <- 0.35  # pause between API calls to respect rate limits

# Output files
OUT_KEGG <- "barplot_KEGG_pathways.png"
OUT_HMDB <- "barplot_HMDB_pathways.png"
PLOT_WIDTH  <- 18   # inches
PLOT_HEIGHT <- 11
PLOT_DPI    <- 180

# =============================================================================
# ── HELPERS ───────────────────────────────────────────────────────────────────
# =============================================================================

parse_condition <- function(col_name, prefix, cond_map) {
  base <- sub(paste0("^", prefix), "", col_name)
  code <- str_extract(base, "^[A-Za-z]+")
  if (is.na(code)) return(NA_character_)
  for (pat in names(cond_map)) {
    if (grepl(pat, toupper(code), perl = TRUE)) return(cond_map[[pat]])
  }
  NA_character_
}

safe_message <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

# =============================================================================
# ── LOAD & AVERAGE INTENSITIES ────────────────────────────────────────────────
# =============================================================================

safe_message("Reading: ", INPUT_FILE)
raw <- read_excel(INPUT_FILE, sheet = SHEET_NAME)

sample_cols <- names(raw)[str_starts(names(raw), SAMPLE_PREFIX)]
col_conds   <- setNames(
  sapply(sample_cols, parse_condition, prefix = SAMPLE_PREFIX, cond_map = CONDITION_MAP),
  sample_cols
)

# Keep only named (string) compounds
df <- raw %>%
  filter(sapply(.data[[NAME_COL]], is.character)) %>%
  mutate(across(all_of(sample_cols),
                ~ suppressWarnings(as.numeric(.)) %>% na_if(0)))

# Average intensities per condition
conditions <- unique(na.omit(col_conds))
for (cond in conditions) {
  cc <- names(col_conds[!is.na(col_conds) & col_conds == cond])
  df[[paste0("mean_", make.names(cond))]] <- rowMeans(df[, cc, drop = FALSE], na.rm = TRUE)
}

mean_cols <- paste0("mean_", make.names(conditions))
df_means  <- df %>%
  select(all_of(c(NAME_COL, "Molecular Formula")), all_of(mean_cols)) %>%
  rename(compound_name = all_of(NAME_COL),
         mol_formula   = "Molecular Formula") %>%
  filter(rowSums(!is.na(across(all_of(mean_cols)))) > 0)

safe_message("Named compounds with intensities: ", nrow(df_means))

# =============================================================================
# ── SECTION 1: KEGG MAPPING ──────────────────────────────────────────────────
# =============================================================================
# Uses KEGGREST to search compound names → KEGG compound IDs → pathway IDs

safe_message("\n=== KEGG PATHWAY MAPPING ===")

if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  stop("KEGGREST not installed. Run: BiocManager::install('KEGGREST')")
}
library(KEGGREST)

# Query KEGG for each compound name
# keggFind searches the KEGG COMPOUND database by name
query_kegg_compound <- function(name) {
  tryCatch({
    Sys.sleep(API_PAUSE_SEC)
    hits <- keggFind("compound", name)
    if (length(hits) == 0) return(character(0))
    names(hits)  # returns KEGG compound IDs like "cpd:C00041"
  }, error = function(e) character(0))
}

# Get pathways for a KEGG compound ID
get_kegg_pathways <- function(cpd_id) {
  tryCatch({
    Sys.sleep(API_PAUSE_SEC)
    info <- keggGet(cpd_id)[[1]]
    pws  <- info$PATHWAY
    if (is.null(pws)) return(character(0))
    names(pws)  # pathway IDs like "hsa00250"
  }, error = function(e) character(0))
}

# Get human-readable pathway name from KEGG pathway ID
get_kegg_pathway_name <- function(pw_id) {
  tryCatch({
    Sys.sleep(API_PAUSE_SEC)
    info <- keggGet(pw_id)[[1]]
    info$NAME %||% pw_id
  }, error = function(e) pw_id)
}

safe_message("Querying KEGG for ", nrow(df_means), " compounds (this may take several minutes)...")

# Map each compound → KEGG compound IDs → pathway IDs
# Cache results to avoid repeat API calls for the same name
kegg_cache <- list()

kegg_map <- df_means %>%
  mutate(
    kegg_cpd = map(compound_name, function(nm) {
      if (!is.null(kegg_cache[[nm]])) return(kegg_cache[[nm]])
      # Try exact name first, then first word of name
      ids <- query_kegg_compound(nm)
      if (length(ids) == 0) {
        short_nm <- word(nm, 1)
        ids <- query_kegg_compound(short_nm)
      }
      kegg_cache[[nm]] <<- ids
      ids
    })
  ) %>%
  filter(lengths(kegg_cpd) > 0) %>%
  mutate(kegg_cpd_first = map_chr(kegg_cpd, 1))  # use first hit

safe_message("Compounds matched to KEGG: ", nrow(kegg_map), " / ", nrow(df_means))

# Fetch pathways for each matched compound
safe_message("Fetching KEGG pathways for matched compounds...")

pathway_cache <- list()
kegg_map <- kegg_map %>%
  mutate(
    kegg_pathway_ids = map(kegg_cpd_first, function(cpd) {
      if (!is.null(pathway_cache[[cpd]])) return(pathway_cache[[cpd]])
      pws <- get_kegg_pathways(cpd)
      pathway_cache[[cpd]] <<- pws
      pws
    })
  ) %>%
  filter(lengths(kegg_pathway_ids) > 0)

# Unnest: one row per compound × pathway
kegg_long <- kegg_map %>%
  select(compound_name, all_of(mean_cols), kegg_pathway_ids) %>%
  unnest(kegg_pathway_ids) %>%
  rename(kegg_pathway_id = kegg_pathway_ids)

# Fetch human pathway names (cache them)
pw_ids_unique   <- unique(kegg_long$kegg_pathway_id)
pw_name_cache   <- list()

safe_message("Fetching pathway names for ", length(pw_ids_unique), " unique KEGG pathways...")
pw_names <- setNames(
  sapply(pw_ids_unique, function(pid) {
    if (!is.null(pw_name_cache[[pid]])) return(pw_name_cache[[pid]])
    nm <- get_kegg_pathway_name(pid)
    pw_name_cache[[pid]] <<- nm
    nm
  }),
  pw_ids_unique
)

kegg_long <- kegg_long %>%
  mutate(pathway_name = pw_names[kegg_pathway_id]) %>%
  # Clean pathway names: remove trailing "[PATH:hsa...]" style suffixes
  mutate(pathway_name = str_remove(pathway_name, "\\s*\\[.*\\]$") %>% trimws())

# =============================================================================
# ── SECTION 2: HMDB MAPPING ──────────────────────────────────────────────────
# =============================================================================
# Uses metaboliteIDmapping (Bioconductor) for name → HMDB ID mapping,
# then queries the HMDB REST API for pathway information.

safe_message("\n=== HMDB PATHWAY MAPPING ===")

if (!requireNamespace("metaboliteIDmapping", quietly = TRUE)) {
  stop("metaboliteIDmapping not installed. Run: BiocManager::install('metaboliteIDmapping')")
}
library(metaboliteIDmapping)

# metaboliteIDmapping provides mappingTable() which maps compound names → HMDB IDs
safe_message("Loading metaboliteIDmapping table...")
id_map <- mappingTable()  # large data frame with name, HMDB, InChIKey, etc.

# Match our compound names to the mapping table (case-insensitive)
our_names_lower <- tolower(df_means$compound_name)
id_map_lower    <- id_map %>% mutate(name_lower = tolower(Name))

matched_hmdb <- df_means %>%
  mutate(name_lower = tolower(compound_name)) %>%
  left_join(
    id_map_lower %>% select(name_lower, HMDB) %>% distinct(),
    by = "name_lower"
  ) %>%
  filter(!is.na(HMDB), HMDB != "")

safe_message("Compounds matched to HMDB ID: ", nrow(matched_hmdb), " / ", nrow(df_means))

# Query HMDB REST API for pathway information per HMDB ID
# HMDB API endpoint: https://hmdb.ca/metabolites/{HMDB_ID}.json
query_hmdb_pathways <- function(hmdb_id) {
  tryCatch({
    Sys.sleep(API_PAUSE_SEC)
    url  <- paste0("https://hmdb.ca/metabolites/", hmdb_id, ".json")
    resp <- GET(url, timeout(15))
    if (status_code(resp) != 200) return(character(0))
    data <- fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = TRUE)
    # Pathways are in data$biological_properties$pathways
    pws <- data$biological_properties$pathways
    if (is.null(pws) || nrow(pws) == 0) return(character(0))
    # Return pathway names (combine KEGG and SMPDB/MetaCyc pathway names)
    unique(pws$name)
  }, error = function(e) character(0))
}

safe_message("Querying HMDB API for ", nrow(matched_hmdb), " matched compounds...")

hmdb_pw_cache <- list()
matched_hmdb <- matched_hmdb %>%
  mutate(
    hmdb_pathways = map(HMDB, function(hid) {
      if (!is.null(hmdb_pw_cache[[hid]])) return(hmdb_pw_cache[[hid]])
      pws <- query_hmdb_pathways(hid)
      hmdb_pw_cache[[hid]] <<- pws
      pws
    })
  ) %>%
  filter(lengths(hmdb_pathways) > 0)

safe_message("Compounds with HMDB pathway data: ", nrow(matched_hmdb))

# Unnest: one row per compound × pathway
hmdb_long <- matched_hmdb %>%
  select(compound_name, all_of(mean_cols), hmdb_pathways) %>%
  unnest(hmdb_pathways) %>%
  rename(pathway_name = hmdb_pathways) %>%
  filter(!is.na(pathway_name), nchar(pathway_name) > 0)

# =============================================================================
# ── SECTION 3: COMPUTE PATHWAY ENRICHMENT (WITHIN-CONDITION, NO CONTROL) ─────
# =============================================================================
# For each condition, for each pathway:
#   enrichment_score = sum of mean intensities of all compounds in that pathway
#   (represents total pathway "activity" / abundance in that condition)
# We then normalise to the total signal in that condition so pathways are
# comparable across conditions with different overall intensity levels.
# Finally we rank by enrichment score and take TOP_N_PATHWAYS.

compute_pathway_enrichment <- function(long_df, mean_cols, conditions) {
  long_df %>%
    pivot_longer(cols = all_of(mean_cols),
                 names_to  = "condition_raw",
                 values_to = "intensity") %>%
    mutate(
      condition = str_remove(condition_raw, "^mean\\.") %>%
        str_replace_all("\\.", " ") %>% trimws()
    ) %>%
    filter(!is.na(intensity), intensity > 0) %>%
    group_by(condition, pathway_name) %>%
    summarise(
      n_compounds       = n_distinct(compound_name),
      sum_intensity     = sum(intensity, na.rm = TRUE),
      mean_intensity    = mean(intensity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_compounds >= MIN_COMPOUNDS) %>%
    group_by(condition) %>%
    mutate(
      # Fractional enrichment: pathway's share of total condition signal
      total_signal        = sum(sum_intensity),
      enrichment_fraction = sum_intensity / total_signal * 100,
      # Log-normalised score for display
      log_mean_intensity  = log10(mean_intensity + 1)
    ) %>%
    ungroup()
}

safe_message("\nComputing KEGG pathway enrichment scores...")
kegg_enrichment <- compute_pathway_enrichment(kegg_long, mean_cols, conditions)

safe_message("Computing HMDB pathway enrichment scores...")
hmdb_enrichment <- compute_pathway_enrichment(hmdb_long, mean_cols, conditions)

# =============================================================================
# ── SECTION 4: BAR PLOTS ──────────────────────────────────────────────────────
# =============================================================================

COND_COLORS <- c(
  "Control"        = "#2C7BB6",
  "Cancer"         = "#C0392B",
  "Torpor + Cancer"= "#8E44AD",
  "Torpor"         = "#117A65"
)

make_barplot <- function(enrichment_df, db_label, top_n = TOP_N_PATHWAYS) {

  # For each condition, get top-N pathways by enrichment_fraction
  top_pws <- enrichment_df %>%
    group_by(condition) %>%
    slice_max(enrichment_fraction, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    pull(pathway_name) %>%
    unique()

  plot_df <- enrichment_df %>%
    filter(pathway_name %in% top_pws) %>%
    # Order pathways by mean enrichment across all conditions
    group_by(pathway_name) %>%
    mutate(mean_enrich_all = mean(enrichment_fraction, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(pathway_name = fct_reorder(pathway_name, mean_enrich_all))

  # Fill in zeros for missing condition × pathway combos
  plot_df <- plot_df %>%
    complete(condition, pathway_name,
             fill = list(enrichment_fraction = 0, n_compounds = 0))

  cond_colors_used <- COND_COLORS[unique(plot_df$condition)]

  ggplot(plot_df, aes(x = enrichment_fraction, y = pathway_name, fill = condition)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.75),
             width = 0.68, alpha = 0.88) +
    geom_text(
      data = filter(plot_df, enrichment_fraction > 0),
      aes(label = paste0("n=", n_compounds)),
      position = position_dodge(width = 0.75),
      hjust = -0.12, size = 2.6, color = "#444444"
    ) +
    scale_fill_manual(values = cond_colors_used, name = "Condition") +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.15)),
      labels = function(x) paste0(round(x, 2), "%")
    ) +
    labs(
      title    = paste0(db_label, " — Pathway Enrichment by Condition"),
      subtitle = paste0(
        "Enrichment = pathway's % share of total condition signal  ·  ",
        "Top ", top_n, " pathways per condition  ·  n = compounds detected"
      ),
      x = "Enrichment Score (% of total condition signal)",
      y = NULL
    ) +
    theme_minimal(base_size = 10.5) +
    theme(
      plot.title       = element_text(face = "bold", size = 13, margin = margin(b = 4)),
      plot.subtitle    = element_text(size = 8.5, color = "#555555", margin = margin(b = 10)),
      axis.text.y      = element_text(size = 9),
      axis.text.x      = element_text(size = 9),
      axis.title.x     = element_text(size = 10, margin = margin(t = 6)),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = 9),
      legend.text      = element_text(size = 9),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "#EEEEEE", linewidth = 0.5),
      plot.margin        = margin(12, 20, 8, 8)
    )
}

safe_message("\nBuilding KEGG bar plot...")
p_kegg <- make_barplot(kegg_enrichment, "KEGG")

safe_message("Building HMDB bar plot...")
p_hmdb <- make_barplot(hmdb_enrichment, "HMDB")

# =============================================================================
# ── SAVE ─────────────────────────────────────────────────────────────────────
# =============================================================================

safe_message("Saving KEGG plot → ", OUT_KEGG)
ggsave(OUT_KEGG, plot = p_kegg,
       width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI, bg = "white")

safe_message("Saving HMDB plot → ", OUT_HMDB)
ggsave(OUT_HMDB, plot = p_hmdb,
       width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI, bg = "white")

safe_message("\nAll done.")
print(p_kegg)
print(p_hmdb)

# Print summary tables
safe_message("\n=== Top KEGG Pathways per Condition ===")
kegg_enrichment %>%
  group_by(condition) %>%
  slice_max(enrichment_fraction, n = 10) %>%
  select(condition, pathway_name, n_compounds, enrichment_fraction) %>%
  mutate(enrichment_fraction = round(enrichment_fraction, 3)) %>%
  print(n = 40)

safe_message("\n=== Top HMDB Pathways per Condition ===")
hmdb_enrichment %>%
  group_by(condition) %>%
  slice_max(enrichment_fraction, n = 10) %>%
  select(condition, pathway_name, n_compounds, enrichment_fraction) %>%
  mutate(enrichment_fraction = round(enrichment_fraction, 3)) %>%
  print(n = 40)
