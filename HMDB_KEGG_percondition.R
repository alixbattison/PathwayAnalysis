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

  # =============================================================================
# EXPORT ADD-ON MODULE — Excel + SVG output
# =============================================================================
# Drop this block at the END of either:
#   • pathway_enrichment_dotplot.R   (differential analysis vs Control)
#   • pathway_enrichment_barplot.R   (within-condition analysis)
#
# It auto-detects which script it is running inside and exports:
#   1. A formatted Excel workbook with GraphPad-ready pivot sheets
#   2. SVG versions of every plot already produced by the parent script
#
# ADDITIONAL PACKAGES NEEDED (beyond each script's own requirements):
#   install.packages("openxlsx")
#
# HOW TO USE:
#   Option A — append inline:
#     Open your script and paste this entire file at the bottom, after the
#     existing save/print block.
#
#   Option B — source at the end:
#     Add this one line to the bottom of your script:
#       source("export_addon.R")
#
# OUTPUT FILES (written to the same directory as your input file):
#   Differential script  →  <stem>_differential_enrichment.xlsx
#                            <stem>_dotplot_<condition>.svg  (one per panel)
#   Within-condition script → <stem>_within_condition_enrichment.xlsx
#                              <stem>_barplot_KEGG.svg
#                              <stem>_barplot_HMDB.svg
# =============================================================================

# ── Package check ─────────────────────────────────────────────────────────────
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  stop(
    "openxlsx is required for Excel export.\n",
    "Install it with: install.packages('openxlsx')"
  )
}
suppressPackageStartupMessages(library(openxlsx))

# ── Detect which script this module is running inside ────────────────────────
# The differential script produces `results` (a named list of data frames,
# one per non-control condition). The within-condition script produces
# `kegg_enrichment` / `hmdb_enrichment`.

.addon_mode <- if (exists("results") && is.list(results) &&
                   all(sapply(results, is.data.frame))) {
  "differential"
} else if (exists("kegg_enrichment") || exists("hmdb_enrichment")) {
  "within"
} else {
  stop(
    "export_addon.R could not detect which script it is running inside.\n",
    "Make sure you source() or paste this file AFTER the main analysis code."
  )
}

message("\n[export_addon] Detected mode: ", .addon_mode)

# ── Output stem (based on INPUT_FILE if it exists, otherwise "output") ────────
.stem <- if (exists("INPUT_FILE") && nchar(INPUT_FILE) > 0) {
  tools::file_path_sans_ext(basename(INPUT_FILE))
} else {
  "metabolomics_output"
}

# ── Shared Excel styling helpers ──────────────────────────────────────────────

.header_style <- function(wb, bg_hex, font_color = "white") {
  createStyle(
    fontName     = "Arial",
    fontSize     = 10,
    fontColour   = font_color,
    fgFill       = bg_hex,
    halign       = "CENTER",
    valign       = "CENTER",
    textDecoration = "bold",
    wrapText     = TRUE,
    border       = "Bottom",
    borderColour = "#CCCCCC"
  )
}

.even_row_style <- createStyle(fgFill = "#F5F5F5", fontName = "Arial", fontSize = 9)
.odd_row_style  <- createStyle(fgFill = "#FFFFFF", fontName = "Arial", fontSize = 9)
.num_style      <- createStyle(numFmt = "0.0000",  fontName = "Arial", fontSize = 9)
.pval_style     <- createStyle(numFmt = "0.000000", fontName = "Arial", fontSize = 9)

.write_sheet <- function(wb, sheet_name, df, header_color,
                          col_widths = NULL, note = NULL) {
  # Truncate sheet name to Excel's 31-char limit
  sheet_name <- substr(sheet_name, 1, 31)

  addWorksheet(wb, sheet_name, gridLines = TRUE)

  start_row <- 1L
  if (!is.null(note) && nchar(note) > 0) {
    writeData(wb, sheet_name, note, startRow = 1, startCol = 1)
    addStyle(wb, sheet_name,
             createStyle(fontColour = "#666666", fontSize = 9,
                         fontName = "Arial", italic = TRUE),
             rows = 1, cols = 1)
    mergeCells(wb, sheet_name, cols = 1:ncol(df), rows = 1)
    start_row <- 2L
  }

  writeDataTable(wb, sheet_name, df,
                 startRow = start_row, startCol = 1,
                 tableStyle = "none", withFilter = TRUE,
                 headerStyle = .header_style(wb, header_color))

  # Alternate row shading
  n_data <- nrow(df)
  if (n_data > 0) {
    even_rows <- start_row + 1 + seq(1, n_data, 2)   # 1-indexed data rows
    odd_rows  <- start_row + 1 + seq(0, n_data, 2)
    even_rows <- even_rows[even_rows <= start_row + n_data]
    odd_rows  <- odd_rows[odd_rows  <= start_row + n_data]
    if (length(even_rows) > 0)
      addStyle(wb, sheet_name, .even_row_style, rows = even_rows,
               cols = 1:ncol(df), gridExpand = TRUE, stack = TRUE)
    if (length(odd_rows) > 0)
      addStyle(wb, sheet_name, .odd_row_style, rows = odd_rows,
               cols = 1:ncol(df), gridExpand = TRUE, stack = TRUE)
  }

  # Numeric column formatting
  num_cols <- which(sapply(df, is.numeric))
  for (ci in num_cols) {
    col_name <- names(df)[ci]
    style    <- if (grepl("pvalue|p_value|P_Value", col_name, ignore.case = TRUE)) {
      .pval_style
    } else {
      .num_style
    }
    addStyle(wb, sheet_name, style,
             rows = (start_row + 1):(start_row + n_data),
             cols = ci, gridExpand = TRUE, stack = TRUE)
  }

  # Column widths
  if (!is.null(col_widths)) {
    setColWidths(wb, sheet_name, cols = seq_along(col_widths), widths = col_widths)
  } else {
    # Auto-size: cap at 50 characters
    auto_widths <- pmin(
      sapply(seq_along(df), function(i) {
        max(nchar(as.character(c(names(df)[i], df[[i]]))), na.rm = TRUE) * 1.1
      }),
      50
    )
    setColWidths(wb, sheet_name, cols = seq_along(df), widths = auto_widths)
  }

  freezePane(wb, sheet_name, firstActiveRow = start_row + 1)
  invisible(wb)
}

# =============================================================================
# ── MODE: DIFFERENTIAL (pathway_enrichment_dotplot.R) ─────────────────────────
# =============================================================================

if (.addon_mode == "differential") {

  message("[export_addon] Building differential enrichment Excel workbook...")

  wb <- createWorkbook()

  # Condition colours (reuse from parent script if available)
  .diff_colors <- list(
    "Cancer vs Control"         = "C0392B",
    "Torpor vs Control"         = "117A65",
    "Torpor + Cancer vs Control"= "8E44AD",
    "Torpor+Cancer vs Control"  = "8E44AD"
  )
  .fallback_color <- "2C3E50"

  # ── Sheet per comparison ──────────────────────────────────────────────────
  all_results_combined <- bind_rows(
    lapply(names(results), function(cond) {
      r <- results[[cond]]
      r$comparison <- paste0(cond, " vs Control")
      r
    })
  )

  for (cond in names(results)) {
    r <- results[[cond]] %>%
      arrange(pvalue) %>%
      mutate(
        significant_p05 = ifelse(pvalue < 0.05, "Yes", "No"),
        direction       = ifelse(log2_fold > 0, "Enriched", "Depleted")
      ) %>%
      rename(
        Pathway             = pathway,
        N_Compounds         = n_compounds,
        Fold_Enrichment     = fold_enrichment,
        Log2_Fold_Change    = log2_fold,
        P_Value             = pvalue,
        Neg_Log10_P         = neg_log10_p,
        Significant_p0.05   = significant_p05,
        Direction           = direction
      ) %>%
      select(Pathway, N_Compounds, Fold_Enrichment, Log2_Fold_Change,
             P_Value, Neg_Log10_P, Significant_p0.05, Direction)

    comp_label <- paste0(cond, " vs Control")
    bg <- .diff_colors[[comp_label]] %||% .fallback_color

    .write_sheet(
      wb, comp_label, r, bg,
      col_widths = c(42, 14, 16, 18, 14, 14, 18, 12),
      note = paste0(
        comp_label,
        " | Mann-Whitney U (two-sided) | Log2_Fold_Change: +ve = enriched in ",
        cond, ", -ve = depleted | keyword-based HMDB/KEGG classification"
      )
    )
  }

  # ── Combined long-format sheet ────────────────────────────────────────────
  combined_long <- all_results_combined %>%
    arrange(comparison, pvalue) %>%
    mutate(
      significant_p05 = ifelse(pvalue < 0.05, "Yes", "No"),
      direction       = ifelse(log2_fold > 0, "Enriched", "Depleted")
    ) %>%
    rename(
      Comparison          = comparison,
      Pathway             = pathway,
      N_Compounds         = n_compounds,
      Fold_Enrichment     = fold_enrichment,
      Log2_Fold_Change    = log2_fold,
      P_Value             = pvalue,
      Neg_Log10_P         = neg_log10_p,
      Significant_p0.05   = significant_p05,
      Direction           = direction
    )

  .write_sheet(
    wb, "All Comparisons", combined_long, "34495E",
    note = "All comparisons combined. Filter by Comparison column in GraphPad or Excel."
  )

  # ── GraphPad-ready pivots ─────────────────────────────────────────────────
  # Log2 fold change pivot
  pivot_log2 <- combined_long %>%
    select(Comparison, Pathway, Log2_Fold_Change) %>%
    pivot_wider(names_from = Comparison, values_from = Log2_Fold_Change,
                values_fill = 0) %>%
    mutate(abs_max = apply(across(-Pathway), 1, function(x) max(abs(x), na.rm = TRUE))) %>%
    arrange(desc(abs_max)) %>%
    select(-abs_max)

  .write_sheet(
    wb, "Log2FC — GraphPad Ready", pivot_log2, "E67E22",
    note = paste0(
      "GraphPad-ready pivot: paste into a Grouped table. ",
      "Rows = pathways, Columns = comparisons. ",
      "Values = log2(fold change) vs Control. Positive = enriched, negative = depleted."
    )
  )

  # -log10(p) pivot
  pivot_p <- combined_long %>%
    select(Comparison, Pathway, Neg_Log10_P) %>%
    pivot_wider(names_from = Comparison, values_from = Neg_Log10_P,
                values_fill = 0) %>%
    mutate(max_p = apply(across(-Pathway), 1, max, na.rm = TRUE)) %>%
    arrange(desc(max_p)) %>%
    select(-max_p)

  .write_sheet(
    wb, "NegLog10P — GraphPad Ready", pivot_p, "C0392B",
    note = paste0(
      "GraphPad-ready pivot: paste into a Grouped table. ",
      "Rows = pathways, Columns = comparisons. ",
      "Values = -log10(p-value). Higher = more statistically significant."
    )
  )

  # Fold enrichment pivot
  pivot_fold <- combined_long %>%
    select(Comparison, Pathway, Fold_Enrichment) %>%
    pivot_wider(names_from = Comparison, values_from = Fold_Enrichment,
                values_fill = 1) %>%
    mutate(max_fe = apply(across(-Pathway), 1, max, na.rm = TRUE)) %>%
    arrange(desc(max_fe)) %>%
    select(-max_fe)

  .write_sheet(
    wb, "FoldEnrichment — GraphPad", pivot_fold, "117A65",
    note = paste0(
      "GraphPad-ready pivot. ",
      "Values = median(condition) / median(Control). Values > 1 = enriched."
    )
  )

  # N compounds pivot
  pivot_n <- combined_long %>%
    select(Comparison, Pathway, N_Compounds) %>%
    pivot_wider(names_from = Comparison, values_from = N_Compounds,
                values_fill = 0)

  .write_sheet(
    wb, "N Compounds", pivot_n, "5D9B3F",
    note = "Number of classified compounds per pathway per comparison (for reference)."
  )

  # ── Save workbook ─────────────────────────────────────────────────────────
  .excel_out <- paste0(.stem, "_differential_enrichment.xlsx")
  saveWorkbook(wb, .excel_out, overwrite = TRUE)
  message("[export_addon] Saved: ", .excel_out)

  # ── SVG export — one file per comparison panel ────────────────────────────
  message("[export_addon] Saving SVG dot plots...")

  # Rebuild panels using the same make_panel() function from the parent script
  p_min <- min(sapply(results, function(r) min(r$neg_log10_p, na.rm = TRUE)))
  p_max <- max(sapply(results, function(r) max(r$neg_log10_p, na.rm = TRUE)))

  for (i in seq_along(non_control)) {
    cond <- non_control[i]
    cols <- if (!is.null(CONDITION_COLORS[[cond]])) {
      CONDITION_COLORS[[cond]]
    } else {
      FALLBACK_COLORS[[min(i, length(FALLBACK_COLORS))]]
    }
    p <- make_panel(results[[cond]], cond, cols, p_min, p_max)
    svg_file <- paste0(.stem, "_dotplot_",
                       gsub("[^A-Za-z0-9]", "_", cond), ".svg")
    ggsave(svg_file, plot = p,
           width  = OUTPUT_WIDTH / length(non_control),
           height = OUTPUT_HEIGHT,
           device = "svg", bg = "white")
    message("[export_addon]   Saved: ", svg_file)
  }

  # Combined SVG (all panels side by side)
  combined_svg <- combine_plots(
    lapply(seq_along(non_control), function(i) {
      cond <- non_control[i]
      cols <- if (!is.null(CONDITION_COLORS[[cond]])) {
        CONDITION_COLORS[[cond]]
      } else {
        FALLBACK_COLORS[[min(i, length(FALLBACK_COLORS))]]
      }
      make_panel(results[[cond]], cond, cols, p_min, p_max)
    }),
    title = "Metabolic Pathway Enrichment"
  )
  svg_combined <- paste0(.stem, "_dotplot_all_panels.svg")
  ggsave(svg_combined, plot = combined_svg,
         width = OUTPUT_WIDTH, height = OUTPUT_HEIGHT,
         device = "svg", bg = "white")
  message("[export_addon]   Saved: ", svg_combined)
}

# =============================================================================
# ── MODE: WITHIN-CONDITION (pathway_enrichment_barplot.R) ─────────────────────
# =============================================================================

if (.addon_mode == "within") {

  message("[export_addon] Building within-condition enrichment Excel workbook...")

  # Detect which enrichment tables exist
  .has_kegg <- exists("kegg_enrichment") && is.data.frame(kegg_enrichment) && nrow(kegg_enrichment) > 0
  .has_hmdb <- exists("hmdb_enrichment") && is.data.frame(hmdb_enrichment) && nrow(hmdb_enrichment) > 0

  # Standardise column names — handle either the API-based or keyword-based
  # versions of the enrichment data frames
  .normalise_enrichment <- function(df) {
    # Rename whichever column holds the pathway name
    if ("pathway_name" %in% names(df) && !"Pathway" %in% names(df))
      df <- rename(df, Pathway = pathway_name)
    if ("pathway" %in% names(df) && !"Pathway" %in% names(df))
      df <- rename(df, Pathway = pathway)
    # Rename condition column
    if ("Condition" %in% names(df) && !"condition" %in% names(df))
      df <- rename(df, condition = Condition)
    # Rename score column
    if ("Enrichment_Score_log10_mean" %in% names(df))
      df <- rename(df, enrich_score = Enrichment_Score_log10_mean)
    if ("enrichment_fraction" %in% names(df) && !"enrich_score" %in% names(df))
      df <- rename(df, enrich_score = enrichment_fraction)
    if ("log_mean_intensity" %in% names(df) && !"enrich_score" %in% names(df))
      df <- rename(df, enrich_score = log_mean_intensity)
    # Rename n_compounds
    if ("N_Compounds" %in% names(df) && !"n_compounds" %in% names(df))
      df <- rename(df, n_compounds = N_Compounds)
    df
  }

  # Detect condition levels (use CONDITION_MAP order if available)
  .conditions_ordered <- if (exists("conditions") && length(conditions) > 0) {
    as.character(conditions)
  } else if (.has_kegg) {
    unique(kegg_enrichment$condition %||% kegg_enrichment$Condition)
  } else {
    unique(hmdb_enrichment$condition %||% hmdb_enrichment$Condition)
  }

  wb <- createWorkbook()

  .db_colors <- list(KEGG = "2C7BB6", HMDB = "8E44AD")
  .db_n_colors <- list(KEGG = "5D9B3F", HMDB = "117A65")

  for (.db in c("KEGG", "HMDB")) {
    .enrich_raw <- if (.db == "KEGG" && .has_kegg) kegg_enrichment else
                   if (.db == "HMDB" && .has_hmdb) hmdb_enrichment else next

    .enrich <- .normalise_enrichment(.enrich_raw)

    # Full long-format sheet
    full_export <- .enrich %>%
      select(condition, Pathway, n_compounds, enrich_score) %>%
      arrange(condition, desc(enrich_score)) %>%
      rename(
        Condition         = condition,
        N_Compounds       = n_compounds,
        Enrichment_Score  = enrich_score
      )

    .write_sheet(
      wb, paste0(.db, " — Full Data"), full_export,
      .db_colors[[.db]],
      note = paste0(
        .db, "-style keyword classification. ",
        "Enrichment_Score = mean log10(intensity+1) per compound (n-normalised). ",
        "Sorted by condition then score descending."
      )
    )

    # GraphPad-ready enrichment score pivot (pathways × conditions)
    pivot_score <- .enrich %>%
      select(condition, Pathway, enrich_score) %>%
      pivot_wider(names_from = condition, values_from = enrich_score,
                  values_fill = 0) %>%
      # Reorder to desired condition order where possible
      {
        cond_cols <- intersect(.conditions_ordered, names(.))
        other_cols <- setdiff(names(.), c("Pathway", cond_cols))
        select(., Pathway, all_of(cond_cols), all_of(other_cols))
      } %>%
      mutate(mean_score = rowMeans(across(-Pathway), na.rm = TRUE)) %>%
      arrange(desc(mean_score)) %>%
      select(-mean_score)

    .write_sheet(
      wb, paste0(.db, " — Score GraphPad"), pivot_score,
      .db_colors[[.db]],
      note = paste0(
        "GraphPad-ready: paste into a Grouped table. ",
        "Rows = pathways, Columns = conditions. ",
        "Values = mean log10(intensity+1) per compound (n-normalised enrichment score)."
      )
    )

    # GraphPad-ready N compounds pivot
    pivot_n <- .enrich %>%
      select(condition, Pathway, n_compounds) %>%
      pivot_wider(names_from = condition, values_from = n_compounds,
                  values_fill = 0L) %>%
      {
        cond_cols <- intersect(.conditions_ordered, names(.))
        other_cols <- setdiff(names(.), c("Pathway", cond_cols))
        select(., Pathway, all_of(cond_cols), all_of(other_cols))
      }

    .write_sheet(
      wb, paste0(.db, " — N Compounds"), pivot_n,
      .db_n_colors[[.db]],
      note = "Number of detected compounds per pathway per condition (for reference / error bars in GraphPad)."
    )
  }

  # ── Save workbook ─────────────────────────────────────────────────────────
  .excel_out <- paste0(.stem, "_within_condition_enrichment.xlsx")
  saveWorkbook(wb, .excel_out, overwrite = TRUE)
  message("[export_addon] Saved: ", .excel_out)

  # ── SVG bar plots ─────────────────────────────────────────────────────────
  message("[export_addon] Saving SVG bar plots...")

  # The parent script may use p_kegg / p_hmdb (API version) or
  # the keyword version that doesn't store plot objects — rebuild if needed.
  .save_barplot_svg <- function(plot_obj, plot_name, width = 18, height = 11) {
    svg_file <- paste0(.stem, "_barplot_", plot_name, ".svg")
    ggsave(svg_file, plot = plot_obj,
           width = width, height = height,
           device = "svg", bg = "white")
    message("[export_addon]   Saved: ", svg_file)
  }

  if (exists("p_kegg") && inherits(p_kegg, "gg")) {
    .save_barplot_svg(p_kegg, "KEGG",
                      width  = ifelse(exists("PLOT_WIDTH"),  PLOT_WIDTH,  18),
                      height = ifelse(exists("PLOT_HEIGHT"), PLOT_HEIGHT, 11))
  } else if (.has_kegg && exists("make_barplot")) {
    .save_barplot_svg(
      make_barplot(kegg_enrichment,
                   if (exists("TOP_N_PATHWAYS")) TOP_N_PATHWAYS else 15),
      "KEGG"
    )
  }

  if (exists("p_hmdb") && inherits(p_hmdb, "gg")) {
    .save_barplot_svg(p_hmdb, "HMDB",
                      width  = ifelse(exists("PLOT_WIDTH"),  PLOT_WIDTH,  18),
                      height = ifelse(exists("PLOT_HEIGHT"), PLOT_HEIGHT, 11))
  } else if (.has_hmdb && exists("make_barplot")) {
    .save_barplot_svg(
      make_barplot(hmdb_enrichment,
                   if (exists("TOP_N_PATHWAYS")) TOP_N_PATHWAYS else 15),
      "HMDB"
    )
  }
}

message("\n[export_addon] Export complete.")
message("[export_addon] Files written:")
message("  Excel : ", .excel_out)
if (.addon_mode == "differential") {
  for (cond in non_control)
    message("  SVG   : ", paste0(.stem, "_dotplot_",
                                  gsub("[^A-Za-z0-9]", "_", cond), ".svg"))
  message("  SVG   : ", paste0(.stem, "_dotplot_all_panels.svg"))
} else {
  if (.has_kegg) message("  SVG   : ", paste0(.stem, "_barplot_KEGG.svg"))
  if (.has_hmdb) message("  SVG   : ", paste0(.stem, "_barplot_HMDB.svg"))
}
