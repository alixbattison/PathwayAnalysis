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
