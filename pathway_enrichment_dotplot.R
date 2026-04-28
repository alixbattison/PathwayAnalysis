# =============================================================================
# Metabolic Pathway Enrichment Dot Plot
# =============================================================================
# Reads a metabolomics Excel file, classifies named compounds into HMDB/KEGG-
# style pathways, computes Mann-Whitney U enrichment vs a control group, and
# produces a multi-panel dot plot (one panel per non-control condition).
#
# INPUT FILE FORMAT EXPECTED:
#   - One sheet containing all data (set SHEET_NAME below)
#   - A column called "Name" with compound names (strings or numeric IDs)
#   - Sample intensity columns prefixed with a common string (e.g. "Plasma_")
#   - Sample names encode condition as a letter prefix before the replicate
#     number, e.g. Plasma_NC7, Plasma_NN1, Plasma_TC10, Plasma_TN5, Plasma_C1
#
# CONDITION MAPPING (edit to match your experiment):
#   The CONDITION_MAP list maps regex patterns (matched against the prefix
#   extracted from each sample column) to a condition label.
#   One condition must be named "Control" — all others are compared to it.

#  Each condition (Cancer, Torpor, Torpor+Cancer) is compared to the Control group (NN + C samples)
#  For each metabolic pathway, it takes all compound intensities in that pathway across all replicates in the condition vs. all replicates in Control, and runs a Mann-Whitney U test
#  The fold enrichment is the ratio of medians (condition / control)
#
# REQUIRED R PACKAGES:
#   readxl, dplyr, tidyr, ggplot2, stringr, purrr, scales
#   Install with: install.packages(c("readxl","dplyr","tidyr","ggplot2",
#                                    "stringr","purrr","scales"))
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(purrr)
library(scales)

# =============================================================================
# ── USER CONFIGURATION ────────────────────────────────────────────────────────
# =============================================================================

INPUT_FILE  <- "RP_Plasma_Pos_Neg_Curated_alldata.xlsx"   # path to your Excel file
SHEET_NAME  <- "Plasma RP PosNeg"                          # sheet to read
NAME_COL    <- "Name"                                       # column with compound names
SAMPLE_PREFIX <- "Plasma_"                                  # prefix identifying sample columns

# Map condition-code regex (matched against prefix before the replicate number)
# to a human-readable label. Exactly one entry must map to "Control".
CONDITION_MAP <- list(
  "^(NN|C)$"  = "Control",         # NN and C samples = control
  "^NC$"       = "Cancer",          # NC = cancer
  "^TC$"       = "Torpor + Cancer", # TC = torpor + cancer
  "^TN$"       = "Torpor"           # TN = torpor
)

# Minimum number of non-NA intensity values required in each group
MIN_N <- 3

# Minimum number of compounds per pathway to include in the plot
MIN_COMPOUNDS <- 2

# Output PNG path (set to NULL to skip saving, plot will still be shown)
OUTPUT_FILE <- "pathway_enrichment_dotplot.png"
OUTPUT_WIDTH  <- 20   # inches
OUTPUT_HEIGHT <- 10   # inches
OUTPUT_DPI    <- 180

# =============================================================================
# ── PATHWAY CLASSIFIER ────────────────────────────────────────────────────────
# =============================================================================
# Keyword-based classifier mapping compound names to HMDB/KEGG-style pathways.
# Each entry: list(pattern = regex_pattern, pathway = "Pathway Name")
# Rules are tested IN ORDER — first match wins.

PATHWAY_RULES <- list(
  list(pattern = "carnitine",
       pathway = "Fatty Acid β-Oxidation\n(Acylcarnitines)"),
  list(pattern = "cholic|chenodeoxychol|deoxycholic|lithocholic|ursodeoxychol|glycocholic|taurocholic|glycocheno|taurochenodeoxy|glycodeoxycholic|tauroursodeoxy|cholanoic|cholanate|cholate|glycolithocholic|sulfolithocholic|obeticholic|hyodeoxycholic",
       pathway = "Bile Acid Metabolism"),
  list(pattern = "cholesterol|cholest|oxysterol|epoxycholest|hydroxycholest|ketocholest|norlanost|lanost|7alpha-hydroxy|7beta-hydroxy|7-dehydro|24-hydroxycholest|25-hydroxy|3beta-hydroxy|ergosterol|desmosterol|lanosterol|hopane",
       pathway = "Steroid / Cholesterol\nMetabolism"),
  list(pattern = "cortisol|cortisone|corticosterone|aldosterone|testosterone|dehydroepiandrosterone|dhea|estradiol|estriol|estrone|progesterone|pregnanediol|pregnenolone|androstan|androst|vitamin d|calciferol|cholecalciferol|ergocalciferol",
       pathway = "Steroid Hormone Metabolism"),
  list(pattern = "sphingo|ceramide|sphinganine|sphingosine|ganglioside|cerebroside|sulfatide|phytosphingosine|sphingomyelin|glucosylceramide|galactosylceramide|lactosylceramide",
       pathway = "Sphingolipid Metabolism"),
  list(pattern = "phosphatidylcholine|lysophosphatidylcholine|phosphatidylethanolamine|lysophosphatidylethanolamine|phosphatidylinositol|lysophosphatidylinositol|phosphatidylserine|phosphatidylglycerol|lysophosphatidic|phosphatidic|glycerophosphocholine|glycerophosphoethanolamine|1-palmitoyl|1-stearoyl|1-oleoyl|1-linoleoyl|sn-glycero-3-phospho",
       pathway = "Glycerophospholipid\nMetabolism"),
  list(pattern = "^lyso|lysophospho",
       pathway = "Glycerophospholipid\nMetabolism"),
  list(pattern = "prostaglandin|thromboxane|leukotriene|lipoxin|\\bhete\\b|\\bhpete\\b|\\beete\\b|hydroxyoctadeca|hydroxyeicosa|eicosanoid|resolvin|protectin|maresin|isoprostane|\\bpgd\\b|\\bpge\\b|\\bpgf\\b|prostacyclin|5-hete|12-hete|15-hete|20-hete|\\bhode\\b|\\bkode\\b",
       pathway = "Eicosanoid / Oxylipin\nMetabolism"),
  list(pattern = "triacylglycerol|diacylglycerol|monoacylglycerol|triglyceride|diglyceride|monoglyceride",
       pathway = "Glycerolipid Metabolism"),
  list(pattern = "hexadecanoic|octadecanoic|\\bdecanoic|dodecanoic|tetradecanoic|palmitic|stearic|\\boleic|linoleic|linolenic|myristic|lauric|capric|caprylic|butyric|valeric|caproic|heptanoic|pelargonic|undecanoic|tridecanoic|pentadecanoic|nonadecanoic|eicosanoic|behenic|nervonic|gadoleic|gondoic|erucic|adipic|suberic|azelaic|sebacic|dodecanedioic|hexadecanedioic|octadecanedioic|traumatic|arachidonic|docosahexaenoic|eicosapentaenoic|\\bdha\\b|\\bepa\\b|fatty acid",
       pathway = "Fatty Acid Metabolism"),
  list(pattern = "glutathione|dipeptide|tripeptide|carnosine|anserine|balenine|homocarnosine",
       pathway = "Peptide / Glutathione\nMetabolism"),
  list(pattern = "taurine|hypotaurine|cysteinesulfinate|homocysteine|cystathionine|s-adenosyl|\\bmethionine",
       pathway = "Sulfur / Methionine\nMetabolism"),
  list(pattern = "\\balanine|\\barginine|asparagine|aspartate|aspartic|\\bcysteine|glutamine|\\bglutamate|glutamic|\\bglycine|histidine|isoleucine|\\bleucine|\\blysine|phenylalanine|\\bproline|\\bserine|threonine|tryptophan|\\btyrosine|\\bvaline|ornithine|citrulline|\\bcreatine|creatinine|sarcosine|\\bbetaine|gamma-aminobutyric|\\bgaba\\b|kynurenine|kynurenic|4-hydroxyphenyl|phenylpyruvate|phenylacetate|homogentisate|vanillylmandelic|homovanillic|dopamine|norepinephrine|epinephrine|\\bserotonin|melatonin|5-hydroxyindole|\\bindole|indican|indoxyl|tryptamine|anthranilate|quinolinate|picolinic|hippurate|3-hydroxyhippurate",
       pathway = "Amino Acid Metabolism"),
  list(pattern = "adenosine|guanosine|\\binosine|xanthosine|hypoxanthine|\\bxanthine|uric acid|\\badenine|\\bguanine|\\bpurine|\\bcaffeine|theobromine|theophylline|paraxanthine|nicotinamide|\\bnad\\b|\\bnadh\\b|\\bnadp\\b|\\bfad\\b|riboflavin",
       pathway = "Purine Metabolism"),
  list(pattern = "cytidine|uridine|thymidine|\\buracil|\\bcytosine|\\bthymine|orotic|pyrimidine|dihydrouracil|beta-ureidopropionate|n-carbamoyl|n-carbamyl",
       pathway = "Pyrimidine Metabolism"),
  list(pattern = "citric acid|\\bcitrate|isocitrate|aconitate|aconitic|oxaloacetate|\\bsuccinate|succinic|\\bfumarate|fumaric|\\bmalate|\\bmalic|\\bmalonate|alpha-ketoglutarate|2-oxoglutarate|\\bpyruvate|pyruvic|acetoacetate|3-hydroxybutyrate|\\boxalic|glyoxylate|methylmalonate|methylsuccinate|methylcitrate|2-hydroxyglutarate",
       pathway = "TCA Cycle / Organic Acids"),
  list(pattern = "\\bglucose|\\bfructose|galactose|\\bmannose|\\bribose|\\bxylose|gluconate|glucuronate|glucuronide|\\blactate|\\blactic|glucosamine|n-acetylglucosamine|n-acetylgalactosamine|sialic|neuraminic|sorbitol|mannitol|xylitol|arabitol|erythritol|inositol|\\bmaltose|cellobiose|trehalose|\\bsucrose|glucarate",
       pathway = "Glycolysis / Carbohydrate\nMetabolism"),
  list(pattern = "\\bvitamin|\\bretinol|retinoic|\\bretinal|tocopherol|tocotrienol|ascorbic|thiamine|pyridoxine|pyridoxal|pyridoxamine|cobalamin|\\bfolate|\\bfolic|\\bbiotin|\\bniacin|pantothenic|pantothenate|coenzyme a|coenzyme q|ubiquinol|ubiquinone|menaquinone|phylloquinone|lipoic",
       pathway = "Vitamins & Cofactors"),
  list(pattern = "flavon|flavonoid|quercetin|kaempferol|luteolin|apigenin|naringenin|hesperetin|catechin|epicatechin|gallocatechin|anthocyanidin|anthocyanin|resveratrol|gallic|ellagic|chlorogenic|caffeic|ferulic|sinapic|coumaric|coumarin|stilbene|\\blignin|lignan|phenylpropanoid|protocatechuate|protocatechuic|vanillic|syringic|p-cresol|hippurate|benzoate",
       pathway = "Polyphenol / Flavonoid\nMetabolism"),
  list(pattern = "terpen|monoterpene|sesquiterpen|diterpen|triterpen|carotenoid|carotene|xanthophyll|lutein|zeaxanthin|lycopene|geraniol|limonene|menthol|camphor|farnesyl|geranylgeranyl|squalene|mevalonate|isoprene|phytol",
       pathway = "Terpenoid / Isoprenoid\nMetabolism")
)

classify_compound <- function(name) {
  if (!is.character(name) || is.na(name) || nchar(trimws(name)) == 0) return(NA_character_)
  n <- tolower(name)
  for (rule in PATHWAY_RULES) {
    if (grepl(rule$pattern, n, perl = TRUE)) return(rule$pathway)
  }
  return(NA_character_)
}

# =============================================================================
# ── HELPER: PARSE CONDITION FROM SAMPLE COLUMN NAME ─────────────────────────
# =============================================================================

parse_condition <- function(col_name, sample_prefix, condition_map) {
  base <- sub(paste0("^", sample_prefix), "", col_name)
  # Extract leading letters (condition code) before any digit or underscore
  code <- str_extract(base, "^[A-Za-z]+")
  if (is.na(code)) return(NA_character_)
  code_upper <- toupper(code)
  for (pattern in names(condition_map)) {
    if (grepl(pattern, code_upper, perl = TRUE)) {
      return(condition_map[[pattern]])
    }
  }
  return(NA_character_)
}

# =============================================================================
# ── LOAD DATA ─────────────────────────────────────────────────────────────────
# =============================================================================

message("Reading: ", INPUT_FILE, " [sheet: ", SHEET_NAME, "]")
raw <- read_excel(INPUT_FILE, sheet = SHEET_NAME)

# Identify sample columns and their conditions
sample_cols <- names(raw)[str_starts(names(raw), SAMPLE_PREFIX)]
if (length(sample_cols) == 0) stop("No sample columns found with prefix: ", SAMPLE_PREFIX)

col_conditions <- setNames(
  sapply(sample_cols, parse_condition,
         sample_prefix = SAMPLE_PREFIX,
         condition_map = CONDITION_MAP),
  sample_cols
)

# Report
message("\nSample → Condition mapping:")
for (s in names(col_conditions)) {
  message("  ", s, " → ", col_conditions[s])
}

condition_labels <- unique(na.omit(col_conditions))
if (!"Control" %in% condition_labels) stop("No samples mapped to 'Control'. Check CONDITION_MAP.")

non_control <- setdiff(condition_labels, "Control")
message("\nConditions to compare vs Control: ", paste(non_control, collapse=", "))

control_cols <- names(col_conditions[col_conditions == "Control" & !is.na(col_conditions)])

# =============================================================================
# ── FILTER & CLASSIFY ─────────────────────────────────────────────────────────
# =============================================================================

df <- raw %>%
  filter(is.character(.data[[NAME_COL]]) | !is.na(.data[[NAME_COL]])) %>%
  filter(sapply(.data[[NAME_COL]], is.character)) %>%
  mutate(pathway = sapply(.data[[NAME_COL]], classify_compound)) %>%
  filter(!is.na(pathway))

# Convert sample columns to numeric, replace 0 with NA
df <- df %>%
  mutate(across(all_of(sample_cols), ~ suppressWarnings(as.numeric(.)) %>%
                  na_if(0)))

message("\nCompounds with pathway annotation: ", nrow(df), " / ", nrow(raw))
message("\nPathway counts:")
print(sort(table(df$pathway), decreasing = TRUE))

# =============================================================================
# ── ENRICHMENT ANALYSIS ───────────────────────────────────────────────────────
# =============================================================================

compute_enrichment <- function(df, group_cols, control_cols, min_n = 3) {
  df %>%
    group_by(pathway) %>%
    summarise(
      n_compounds = n(),
      ctrl_vals   = list(na.omit(unlist(across(all_of(control_cols))))),
      grp_vals    = list(na.omit(unlist(across(all_of(group_cols))))),
      .groups = "drop"
    ) %>%
    filter(
      n_compounds >= MIN_COMPOUNDS,
      lengths(ctrl_vals) >= min_n,
      lengths(grp_vals)  >= min_n
    ) %>%
    mutate(
      pvalue = map2_dbl(grp_vals, ctrl_vals, ~ {
        tryCatch(wilcox.test(.x, .y, exact = FALSE)$p.value, error = function(e) 1)
      }),
      fold_enrichment = map2_dbl(grp_vals, ctrl_vals, ~ median(.x) / pmax(median(.y), 1e-9)),
      log2_fold       = log2(pmax(fold_enrichment, 1e-6)),
      neg_log10_p     = -log10(pvalue + 1e-300)
    ) %>%
    select(-ctrl_vals, -grp_vals) %>%
    arrange(pvalue)
}

# Compute results for every non-control condition
results <- lapply(non_control, function(cond) {
  gcols <- names(col_conditions[col_conditions == cond & !is.na(col_conditions)])
  message("\nComputing enrichment: ", cond, " (n=", length(gcols), " samples) vs Control (n=", length(control_cols), ")")
  res <- compute_enrichment(df, gcols, control_cols)
  res$condition <- cond
  res
})
names(results) <- non_control

all_results <- bind_rows(results)

# =============================================================================
# ── PLOT ──────────────────────────────────────────────────────────────────────
# =============================================================================

# Global p-value range for consistent bubble sizing across panels
p_min <- min(all_results$neg_log10_p, na.rm = TRUE)
p_max <- max(all_results$neg_log10_p, na.rm = TRUE)

# Colour scheme: one up-colour + one down-colour per condition
CONDITION_COLORS <- list(
  "Cancer"          = c(up = "#C0392B", dn = "#2874A6"),
  "Torpor + Cancer" = c(up = "#8E44AD", dn = "#1A5276"),
  "Torpor"          = c(up = "#117A65", dn = "#784212")
)
# Fallback palette for extra conditions
FALLBACK_COLORS <- list(
  c(up = "#D35400", dn = "#1F618D"),
  c(up = "#1E8449", dn = "#6C3483"),
  c(up = "#B7950B", dn = "#2C3E50")
)

make_panel <- function(res, condition_label, colors, p_min, p_max) {
  res <- res %>%
    filter(n_compounds >= MIN_COMPOUNDS) %>%
    arrange(log2_fold) %>%
    mutate(
      pw_label    = pathway,
      pw_label    = factor(pw_label, levels = pw_label),  # preserve order
      direction   = ifelse(log2_fold > 0, "up", "down"),
      dot_color   = ifelse(log2_fold > 0, colors["up"], colors["dn"]),
      dot_alpha   = ifelse(pvalue < 0.05, 1.0, 0.35),
      dot_size    = 2 + (neg_log10_p - p_min) / (p_max - p_min + 1e-9) * 10,
      significant = pvalue < 0.05
    )

  x_abs_max <- max(abs(res$log2_fold), na.rm = TRUE) + 0.5

  ggplot(res, aes(x = log2_fold, y = pw_label)) +
    # Alternating row bands
    geom_rect(aes(xmin = -Inf, xmax = Inf,
                  ymin = as.numeric(pw_label) - 0.5,
                  ymax = as.numeric(pw_label) + 0.5,
                  fill = as.numeric(pw_label) %% 2 == 0),
              inherit.aes = FALSE) +
    scale_fill_manual(values = c("TRUE" = "#F5F5F5", "FALSE" = "white"), guide = "none") +
    # Reference line
    geom_vline(xintercept = 0, linetype = "dashed", color = "#888888", linewidth = 0.6) +
    # Bubbles
    geom_point(aes(size = dot_size, color = dot_color, alpha = dot_alpha)) +
    # Black ring for significant
    geom_point(data = filter(res, significant),
               aes(size = dot_size),
               shape = 21, fill = NA, color = "black", stroke = 0.8) +
    # n= annotation
    geom_text(aes(x = x_abs_max * 0.96, label = paste0("n=", n_compounds)),
              hjust = 1, size = 2.8, color = "#666666", fontface = "italic") +
    scale_color_identity() +
    scale_alpha_identity() +
    scale_size_identity() +
    xlim(-x_abs_max, x_abs_max) +
    labs(
      title = paste0(condition_label, "  vs  Control"),
      x     = "log₂ Fold Enrichment",
      y     = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title        = element_text(face = "bold", size = 11, hjust = 0.5, margin = margin(b = 8)),
      axis.text.y       = element_text(size = 9, color = "#222222"),
      axis.text.x       = element_text(size = 9),
      axis.title.x      = element_text(size = 9.5, margin = margin(t = 6)),
      panel.grid.major.x = element_line(color = "#DDDDDD", linewidth = 0.4),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.border       = element_blank(),
      axis.line.x        = element_line(color = "#AAAAAA", linewidth = 0.4),
      plot.margin        = margin(8, 16, 8, 8)
    )
}

# Build panels
panels <- lapply(seq_along(non_control), function(i) {
  cond <- non_control[i]
  cols <- if (!is.null(CONDITION_COLORS[[cond]])) {
    CONDITION_COLORS[[cond]]
  } else {
    FALLBACK_COLORS[[min(i, length(FALLBACK_COLORS))]]
  }
  make_panel(results[[cond]], cond, cols, p_min, p_max)
})

# Combine panels with patchwork (if available) or gridExtra
combine_plots <- function(panels, title) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    combined <- Reduce(`|`, panels) +
      plot_annotation(
        title    = title,
        subtitle = paste0(
          "Bubble size = \u2212log\u2081\u2080(p)  \u00b7  ",
          "Black ring = p < 0.05  \u00b7  ",
          "Faded = non-significant  \u00b7  ",
          "Mann\u2013Whitney U test  \u00b7  HMDB/KEGG-style classification"
        ),
        theme = theme(
          plot.title    = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 4)),
          plot.subtitle = element_text(size = 9,  hjust = 0.5, color = "#555555", margin = margin(b = 10))
        )
      )
    return(combined)
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    library(gridExtra)
    return(gridExtra::arrangeGrob(grobs = panels, nrow = 1,
           top = grid::textGrob(title, gp = grid::gpar(fontface = "bold", fontsize = 13))))
  } else {
    message("Install 'patchwork' or 'gridExtra' to combine panels. Printing last panel only.")
    return(panels[[length(panels)]])
  }
}

combined <- combine_plots(
  panels,
  title = "Metabolic Pathway Enrichment — Plasma Metabolomics"
)

# =============================================================================
# ── SAVE & DISPLAY ────────────────────────────────────────────────────────────
# =============================================================================

if (!is.null(OUTPUT_FILE)) {
  message("\nSaving plot to: ", OUTPUT_FILE)
  ggsave(OUTPUT_FILE, plot = combined,
         width = OUTPUT_WIDTH, height = OUTPUT_HEIGHT,
         dpi = OUTPUT_DPI, bg = "white")
  message("Done.")
}

print(combined)

# Also print enrichment tables
message("\n=== Enrichment Results ===")
for (cond in non_control) {
  message("\n--- ", cond, " vs Control ---")
  print(results[[cond]] %>%
    select(pathway, n_compounds, fold_enrichment, log2_fold, pvalue) %>%
    mutate(across(c(fold_enrichment, log2_fold), round, 3),
           pvalue = signif(pvalue, 3)) %>%
    arrange(pvalue),
    n = Inf)
}

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
