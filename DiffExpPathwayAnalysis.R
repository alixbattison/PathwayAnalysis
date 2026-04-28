# =============================================================================
# Metabolic Pathway Differential Enrichment — KEGG API & HMDB API
# =============================================================================
# Reads a metabolomics Excel file, maps compound names to KEGG and HMDB
# pathways via live API calls (run separately for each database), computes
# Mann-Whitney U differential enrichment of each pathway vs a Control group,
# and produces:
#   • Multi-panel dot plots (SVG + PNG), one panel per non-control condition
#   • A formatted Excel workbook with GraphPad-ready pivot sheets
#
# DATABASES (run independently):
#   KEGG  — uses KEGGREST (Bioconductor): name → compound ID → pathway
#   HMDB  — uses metaboliteIDmapping (Bioconductor) + HMDB REST API:
#            name → HMDB ID → pathway via https://hmdb.ca/metabolites/{ID}.json
#
# CONDITIONS (edit CONDITION_MAP to match your experiment):
#   NC  = Cancer          (no torpor, cancer)
#   NN, C = Control       (no torpor, no cancer)
#   TC  = Torpor + Cancer
#   TN  = Torpor
#
# REQUIRED PACKAGES:
#   CRAN:         readxl, dplyr, tidyr, ggplot2, stringr, purrr,
#                 scales, openxlsx, httr, jsonlite
#   Bioconductor: KEGGREST, metaboliteIDmapping
#
# Install:
#   install.packages(c("readxl","dplyr","tidyr","ggplot2","stringr",
#                      "purrr","scales","openxlsx","httr","jsonlite"))
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
  library(openxlsx)
  library(httr)
  library(jsonlite)
})

# =============================================================================
# ── USER CONFIGURATION ────────────────────────────────────────────────────────
# =============================================================================

INPUT_FILE    <- "RP_Plasma_Pos_Neg_Curated_alldata.xlsx"
SHEET_NAME    <- "Plasma RP PosNeg"
NAME_COL      <- "Name"
SAMPLE_PREFIX <- "Plasma_"

# Condition mapping: regex matched against leading letters of each sample column.
# Exactly ONE entry must map to "Control".
CONDITION_MAP <- list(
  "^(NN|C)$" = "Control",
  "^NC$"      = "Cancer",
  "^TC$"      = "Torpor + Cancer",
  "^TN$"      = "Torpor"
)

# Minimum non-NA values per group to run the test
MIN_N         <- 3
# Minimum compounds per pathway to include in results
MIN_COMPOUNDS <- 2
# Pause between API calls (seconds) — KEGG allows ~10 req/s via KEGGREST
API_PAUSE_SEC <- 0.35

# Which database to run: "KEGG", "HMDB", or "BOTH"
RUN_DATABASE  <- "BOTH"

# Top N pathways to show per panel in the dot plot
TOP_N_PLOT    <- 20

# Output files (NULL = skip saving that format)
OUT_PNG_KEGG  <- "differential_dotplot_KEGG.png"
OUT_PNG_HMDB  <- "differential_dotplot_HMDB.png"
OUT_SVG_KEGG  <- "differential_dotplot_KEGG.svg"
OUT_SVG_HMDB  <- "differential_dotplot_HMDB.svg"
OUT_EXCEL     <- "differential_enrichment_API.xlsx"
PLOT_WIDTH    <- 20   # inches (total, shared across panels)
PLOT_HEIGHT   <- 11
PLOT_DPI      <- 180

# =============================================================================
# ── HELPERS ───────────────────────────────────────────────────────────────────
# =============================================================================

safe_message <- function(...) message("[", format(Sys.time(), "%H:%M:%S"), "] ", ...)

parse_condition <- function(col_name, prefix, cond_map) {
  base <- sub(paste0("^", prefix), "", col_name)
  code <- str_extract(base, "^[A-Za-z]+")
  if (is.na(code)) return(NA_character_)
  for (pat in names(cond_map)) {
    if (grepl(pat, toupper(code), perl = TRUE)) return(cond_map[[pat]])
  }
  NA_character_
}

# =============================================================================
# ── LOAD & PREPARE DATA ───────────────────────────────────────────────────────
# =============================================================================

safe_message("Reading: ", INPUT_FILE, " [", SHEET_NAME, "]")
raw <- read_excel(INPUT_FILE, sheet = SHEET_NAME)

sample_cols <- names(raw)[str_starts(names(raw), SAMPLE_PREFIX)]
if (length(sample_cols) == 0)
  stop("No sample columns found with prefix: ", SAMPLE_PREFIX)

col_conds <- setNames(
  sapply(sample_cols, parse_condition,
         prefix = SAMPLE_PREFIX, cond_map = CONDITION_MAP),
  sample_cols
)

safe_message("Sample → Condition mapping:")
for (s in names(col_conds)) message("  ", s, " → ", col_conds[[s]])

if (!"Control" %in% col_conds)
  stop("No samples mapped to 'Control'. Check CONDITION_MAP.")

# Filter to string-named compounds; zero → NA
df <- raw %>%
  filter(sapply(.data[[NAME_COL]], is.character)) %>%
  mutate(across(all_of(sample_cols),
                ~ suppressWarnings(as.numeric(.)) %>% na_if(0)))

safe_message("Named compounds: ", nrow(df))

all_conditions  <- unique(na.omit(col_conds))
non_control     <- setdiff(all_conditions, "Control")
control_cols    <- names(col_conds[!is.na(col_conds) & col_conds == "Control"])

safe_message("Non-control conditions: ", paste(non_control, collapse = ", "))

# =============================================================================
# ── KEGG PATHWAY MAPPING ──────────────────────────────────────────────────────
# =============================================================================

map_kegg <- function(df, name_col, api_pause = API_PAUSE_SEC) {

  if (!requireNamespace("KEGGREST", quietly = TRUE))
    stop("KEGGREST not installed. Run: BiocManager::install('KEGGREST')")
  library(KEGGREST)

  compound_names <- df[[name_col]]
  safe_message("KEGG: querying ", length(compound_names), " compound names...")

  # ── Step 1: name → KEGG compound ID ──────────────────────────────────────
  cpd_cache <- list()

  query_kegg_cpd <- function(nm) {
    if (!is.null(cpd_cache[[nm]])) return(cpd_cache[[nm]])
    Sys.sleep(api_pause)
    ids <- tryCatch({
      hits <- keggFind("compound", nm)
      if (length(hits) > 0) names(hits) else {
        # Retry with first word of name
        short <- word(nm, 1)
        hits2 <- keggFind("compound", short)
        if (length(hits2) > 0) names(hits2) else character(0)
      }
    }, error = function(e) character(0))
    cpd_cache[[nm]] <<- ids
    ids
  }

  df_cpd <- df %>%
    mutate(
      kegg_cpd_ids = map(.data[[name_col]], query_kegg_cpd),
      kegg_cpd     = map_chr(kegg_cpd_ids,
                             ~ if (length(.x) > 0) .x[[1]] else NA_character_)
    ) %>%
    filter(!is.na(kegg_cpd))

  safe_message("KEGG: matched ", nrow(df_cpd), " / ", nrow(df), " compounds to KEGG IDs")

  # ── Step 2: compound ID → pathway IDs ────────────────────────────────────
  pw_id_cache <- list()

  get_cpd_pathways <- function(cpd_id) {
    if (!is.null(pw_id_cache[[cpd_id]])) return(pw_id_cache[[cpd_id]])
    Sys.sleep(api_pause)
    pws <- tryCatch({
      info <- keggGet(cpd_id)[[1]]
      if (!is.null(info$PATHWAY)) names(info$PATHWAY) else character(0)
    }, error = function(e) character(0))
    pw_id_cache[[cpd_id]] <<- pws
    pws
  }

  safe_message("KEGG: fetching pathways for ", nrow(df_cpd), " matched compounds...")

  df_pw <- df_cpd %>%
    mutate(pathway_ids = map(kegg_cpd, get_cpd_pathways)) %>%
    filter(lengths(pathway_ids) > 0) %>%
    unnest(pathway_ids) %>%
    rename(pathway_id = pathway_ids)

  # ── Step 3: pathway ID → human-readable name ─────────────────────────────
  pw_ids_unique <- unique(df_pw$pathway_id)
  safe_message("KEGG: resolving names for ", length(pw_ids_unique), " unique pathways...")

  pw_name_cache <- list()
  pw_name_map <- setNames(
    sapply(pw_ids_unique, function(pid) {
      if (!is.null(pw_name_cache[[pid]])) return(pw_name_cache[[pid]])
      Sys.sleep(api_pause)
      nm <- tryCatch({
        info <- keggGet(pid)[[1]]
        nm   <- info$NAME %||% pid
        # Strip trailing " - Homo sapiens (human)" style suffixes
        str_remove(nm, "\\s*-\\s*[A-Z][a-z].*$") %>% trimws()
      }, error = function(e) pid)
      pw_name_cache[[pid]] <<- nm
      nm
    }),
    pw_ids_unique
  )

  df_pw <- df_pw %>%
    mutate(pathway_name = pw_name_map[pathway_id])

  safe_message("KEGG: final long table: ", nrow(df_pw), " compound × pathway rows")
  df_pw
}

# =============================================================================
# ── HMDB PATHWAY MAPPING ──────────────────────────────────────────────────────
# =============================================================================

map_hmdb <- function(df, name_col, api_pause = API_PAUSE_SEC) {

  if (!requireNamespace("metaboliteIDmapping", quietly = TRUE))
    stop("metaboliteIDmapping not installed. Run: BiocManager::install('metaboliteIDmapping')")
  library(metaboliteIDmapping)

  safe_message("HMDB: loading metaboliteIDmapping table...")
  id_tbl <- mappingTable()

  # ── Step 1: compound name → HMDB ID ──────────────────────────────────────
  name_lower <- tolower(df[[name_col]])
  id_lower   <- id_tbl %>%
    mutate(.name_lower = tolower(Name)) %>%
    select(.name_lower, HMDB) %>%
    distinct() %>%
    filter(!is.na(HMDB), HMDB != "")

  df_hmdb <- df %>%
    mutate(.name_lower = tolower(.data[[name_col]])) %>%
    left_join(id_lower, by = ".name_lower") %>%
    select(-.name_lower) %>%
    filter(!is.na(HMDB))

  safe_message("HMDB: matched ", nrow(df_hmdb), " / ", nrow(df), " compounds to HMDB IDs")

  # ── Step 2: HMDB ID → pathway names via REST API ─────────────────────────
  hmdb_ids_unique <- unique(df_hmdb$HMDB)
  safe_message("HMDB: querying REST API for ", length(hmdb_ids_unique), " unique HMDB IDs...")

  pw_cache <- list()

  query_hmdb_pws <- function(hmdb_id) {
    if (!is.null(pw_cache[[hmdb_id]])) return(pw_cache[[hmdb_id]])
    Sys.sleep(api_pause)
    pws <- tryCatch({
      url  <- paste0("https://hmdb.ca/metabolites/", hmdb_id, ".json")
      resp <- GET(url, timeout(20))
      if (status_code(resp) != 200) return(character(0))
      data <- fromJSON(content(resp, "text", encoding = "UTF-8"), flatten = TRUE)
      bio  <- data$biological_properties$pathways
      if (is.null(bio) || nrow(bio) == 0) return(character(0))
      unique(bio$name[!is.na(bio$name)])
    }, error = function(e) character(0))
    pw_cache[[hmdb_id]] <<- pws
    pws
  }

  df_pw <- df_hmdb %>%
    mutate(pathway_names = map(HMDB, query_hmdb_pws)) %>%
    filter(lengths(pathway_names) > 0) %>%
    unnest(pathway_names) %>%
    rename(pathway_name = pathway_names) %>%
    filter(!is.na(pathway_name), nchar(pathway_name) > 0)

  safe_message("HMDB: final long table: ", nrow(df_pw), " compound × pathway rows")
  df_pw
}

# =============================================================================
# ── DIFFERENTIAL ENRICHMENT TEST ─────────────────────────────────────────────
# =============================================================================
# For each condition vs Control, for each pathway:
#   - Pool all compound intensities within the pathway across replicates
#   - Mann-Whitney U test (two-sided)
#   - Fold enrichment = median(condition) / median(control)

compute_differential <- function(df_long, sample_cols, col_conds,
                                  non_control, control_cols,
                                  min_n = MIN_N, min_cpds = MIN_COMPOUNDS) {

  name_col_in_long <- if ("compound_name" %in% names(df_long)) "compound_name" else NAME_COL

  results <- list()

  for (cond in non_control) {
    cond_cols <- names(col_conds[!is.na(col_conds) & col_conds == cond])
    safe_message("  Testing: ", cond, " vs Control (",
                 length(cond_cols), " vs ", length(control_cols), " samples)")

    rows <- df_long %>%
      group_by(pathway_name) %>%
      summarise(
        n_compounds  = n_distinct(.data[[name_col_in_long]]),
        ctrl_vals    = list(na.omit(unlist(across(all_of(control_cols))))),
        cond_vals    = list(na.omit(unlist(across(all_of(cond_cols))))),
        .groups      = "drop"
      ) %>%
      filter(
        n_compounds >= min_cpds,
        lengths(ctrl_vals) >= min_n,
        lengths(cond_vals) >= min_n
      ) %>%
      mutate(
        pvalue = map2_dbl(cond_vals, ctrl_vals, ~ {
          tryCatch(
            wilcox.test(.x, .y, exact = FALSE)$p.value,
            error = function(e) 1
          )
        }),
        fold_enrichment = map2_dbl(cond_vals, ctrl_vals, ~ {
          med_c <- median(.y, na.rm = TRUE)
          median(.x, na.rm = TRUE) / max(med_c, 1e-9)
        }),
        log2_fold   = log2(pmax(fold_enrichment, 1e-6)),
        neg_log10_p = -log10(pvalue + 1e-300),
        direction   = ifelse(log2_fold > 0, "Enriched", "Depleted"),
        significant = pvalue < 0.05
      ) %>%
      select(-ctrl_vals, -cond_vals) %>%
      arrange(pvalue)

    results[[cond]] <- rows
  }
  results
}

# =============================================================================
# ── DOT PLOT ──────────────────────────────────────────────────────────────────
# =============================================================================

CONDITION_COLORS <- list(
  "Cancer"          = c(up = "#C0392B", dn = "#2874A6"),
  "Torpor + Cancer" = c(up = "#8E44AD", dn = "#1A5276"),
  "Torpor"          = c(up = "#117A65", dn = "#784212")
)
FALLBACK_COLORS <- list(
  c(up = "#D35400", dn = "#1F618D"),
  c(up = "#1E8449", dn = "#6C3483"),
  c(up = "#B7950B", dn = "#2C3E50")
)

make_dotplot <- function(results, db_label, top_n = TOP_N_PLOT,
                          out_png = NULL, out_svg = NULL,
                          width = PLOT_WIDTH, height = PLOT_HEIGHT,
                          dpi = PLOT_DPI) {

  all_nlp <- unlist(lapply(results, `[[`, "neg_log10_p"))
  p_min   <- min(all_nlp, na.rm = TRUE)
  p_max   <- max(all_nlp, na.rm = TRUE)

  make_panel <- function(res, cond_label, colors) {
    res <- res %>%
      filter(n_compounds >= MIN_COMPOUNDS) %>%
      slice_max(neg_log10_p, n = top_n, with_ties = FALSE) %>%
      arrange(log2_fold) %>%
      mutate(
        pw_label  = factor(pathway_name, levels = pathway_name),
        dot_color = ifelse(log2_fold > 0, colors["up"], colors["dn"]),
        dot_alpha = ifelse(significant, 1.0, 0.35),
        dot_size  = 2 + (neg_log10_p - p_min) / (p_max - p_min + 1e-9) * 10
      )

    x_lim <- max(abs(res$log2_fold), na.rm = TRUE) + 0.5

    ggplot(res, aes(x = log2_fold, y = pw_label)) +
      geom_rect(
        aes(xmin = -Inf, xmax = Inf,
            ymin = as.numeric(pw_label) - 0.5,
            ymax = as.numeric(pw_label) + 0.5,
            fill = as.numeric(pw_label) %% 2 == 0),
        inherit.aes = FALSE
      ) +
      scale_fill_manual(values = c("TRUE" = "#F5F5F5", "FALSE" = "white"),
                        guide = "none") +
      geom_vline(xintercept = 0, linetype = "dashed",
                 color = "#888888", linewidth = 0.6) +
      geom_point(aes(size = dot_size, color = dot_color, alpha = dot_alpha)) +
      geom_point(
        data = filter(res, significant),
        aes(size = dot_size),
        shape = 21, fill = NA, color = "black", stroke = 0.9
      ) +
      geom_text(
        aes(x = x_lim * 0.96,
            label = paste0("n=", n_compounds)),
        hjust = 1, size = 2.8, color = "#666666", fontface = "italic"
      ) +
      scale_color_identity() +
      scale_alpha_identity() +
      scale_size_identity() +
      xlim(-x_lim, x_lim) +
      labs(
        title = paste0(cond_label, "  vs  Control"),
        x     = "log\u2082 Fold Enrichment",
        y     = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(
        plot.title         = element_text(face = "bold", size = 11,
                                          hjust = 0.5, margin = margin(b = 8)),
        axis.text.y        = element_text(size = 8.5),
        axis.text.x        = element_text(size = 9),
        axis.title.x       = element_text(size = 9.5, margin = margin(t = 6)),
        panel.grid.major.x = element_line(color = "#DDDDDD", linewidth = 0.4),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.border       = element_blank(),
        axis.line.x        = element_line(color = "#AAAAAA", linewidth = 0.4),
        plot.margin        = margin(8, 16, 8, 8)
      )
  }

  panels <- lapply(seq_along(results), function(i) {
    cond  <- names(results)[i]
    cols  <- CONDITION_COLORS[[cond]] %||% FALLBACK_COLORS[[min(i, 3)]]
    make_panel(results[[cond]], cond, cols)
  })

  # Combine panels
  combined <- if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    Reduce(`|`, panels) +
      plot_annotation(
        title    = paste0(db_label, " — Metabolic Pathway Differential Enrichment"),
        subtitle = paste0(
          "Bubble size = \u2212log\u2081\u2080(p)  \u00b7  ",
          "Black ring = p < 0.05  \u00b7  Faded = non-significant  \u00b7  ",
          "Mann\u2013Whitney U  \u00b7  Top ", top_n, " pathways by significance per panel"
        ),
        theme = theme(
          plot.title    = element_text(face = "bold", size = 13,
                                       hjust = 0.5, margin = margin(b = 4)),
          plot.subtitle = element_text(size = 9, hjust = 0.5,
                                       color = "#555555", margin = margin(b = 10))
        )
      )
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    library(gridExtra)
    gridExtra::arrangeGrob(grobs = panels, nrow = 1)
  } else {
    panels[[length(panels)]]
  }

  if (!is.null(out_png)) {
    ggsave(out_png, plot = combined, width = width, height = height,
           dpi = dpi, bg = "white")
    safe_message("Saved PNG: ", out_png)
  }
  if (!is.null(out_svg)) {
    ggsave(out_svg, plot = combined, width = width, height = height,
           device = "svg", bg = "white")
    safe_message("Saved SVG: ", out_svg)
  }

  print(combined)
  invisible(combined)
}

# =============================================================================
# ── EXCEL EXPORT ──────────────────────────────────────────────────────────────
# =============================================================================

export_excel <- function(results_list, out_file) {

  wb <- createWorkbook()

  .hdr <- function(bg, fc = "white")
    createStyle(fontName = "Arial", fontSize = 10, fontColour = fc,
                fgFill = bg, halign = "CENTER", valign = "CENTER",
                textDecoration = "bold", wrapText = TRUE,
                border = "Bottom", borderColour = "#CCCCCC")

  .write <- function(wb, sname, df, bg, note = NULL) {
    sname <- substr(sname, 1, 31)
    addWorksheet(wb, sname, gridLines = TRUE)
    r0 <- 1L
    if (!is.null(note)) {
      writeData(wb, sname, note, startRow = 1, startCol = 1)
      addStyle(wb, sname,
               createStyle(fontColour = "#666666", fontSize = 9,
                           fontName = "Arial", italic = TRUE),
               rows = 1, cols = 1)
      mergeCells(wb, sname, cols = 1:ncol(df), rows = 1)
      r0 <- 2L
    }
    writeDataTable(wb, sname, df, startRow = r0, startCol = 1,
                   tableStyle = "none", withFilter = TRUE,
                   headerStyle = .hdr(bg))
    n <- nrow(df)
    if (n > 0) {
      for (r in seq_len(n)) {
        fg <- if (r %% 2 == 0) "#F5F5F5" else "#FFFFFF"
        addStyle(wb, sname,
                 createStyle(fgFill = fg, fontName = "Arial", fontSize = 9),
                 rows = r0 + r, cols = 1:ncol(df),
                 gridExpand = TRUE, stack = TRUE)
      }
    }
    # Numeric formatting
    for (ci in which(sapply(df, is.numeric))) {
      fmt <- if (grepl("pvalue|p_value|P_Value", names(df)[ci], ignore.case = TRUE))
               "0.000000" else "0.0000"
      addStyle(wb, sname,
               createStyle(numFmt = fmt, fontName = "Arial", fontSize = 9),
               rows = (r0 + 1):(r0 + n), cols = ci,
               gridExpand = TRUE, stack = TRUE)
    }
    # Auto column widths
    widths <- pmin(
      sapply(seq_along(df), function(i)
        max(nchar(as.character(c(names(df)[i], df[[i]]))), na.rm = TRUE) * 1.15),
      52
    )
    setColWidths(wb, sname, cols = seq_along(df), widths = widths)
    freezePane(wb, sname, firstActiveRow = r0 + 1)
  }

  .comp_colors <- c(
    "Cancer vs Control"          = "C0392B",
    "Torpor vs Control"          = "117A65",
    "Torpor + Cancer vs Control" = "8E44AD"
  )

  for (db in names(results_list)) {
    results <- results_list[[db]]

    all_combined <- bind_rows(lapply(names(results), function(cond) {
      results[[cond]] %>%
        mutate(Comparison = paste0(cond, " vs Control"),
               Pathway    = pathway_name) %>%
        rename(N_Compounds      = n_compounds,
               Fold_Enrichment  = fold_enrichment,
               Log2_Fold_Change = log2_fold,
               P_Value          = pvalue,
               Neg_Log10_P      = neg_log10_p,
               Direction        = direction,
               Significant_p05  = significant) %>%
        select(Comparison, Pathway, N_Compounds, Fold_Enrichment,
               Log2_Fold_Change, P_Value, Neg_Log10_P,
               Significant_p05, Direction)
    }))

    # Per-comparison sheets
    for (cond in names(results)) {
      comp  <- paste0(cond, " vs Control")
      color <- .comp_colors[comp] %||% "34495E"
      sub   <- filter(all_combined, Comparison == comp) %>% select(-Comparison)
      .write(wb, paste0(db, " \u2014 ", comp), sub, color,
             note = paste0(
               db, " | ", comp,
               " | Mann-Whitney U (two-sided) | Log2_Fold_Change: ",
               "+ve = enriched in condition, -ve = depleted | ",
               db, " API-based pathway mapping"
             ))
    }

    # All comparisons combined
    .write(wb, paste0(db, " \u2014 All"), all_combined, "34495E",
           note = paste0(db, " | All comparisons. Filter by Comparison column."))

    # Log2FC GraphPad pivot
    pivot_l2 <- all_combined %>%
      select(Comparison, Pathway, Log2_Fold_Change) %>%
      pivot_wider(names_from = Comparison, values_from = Log2_Fold_Change,
                  values_fill = 0) %>%
      mutate(.abs = apply(across(-Pathway), 1, function(x) max(abs(x), na.rm = TRUE))) %>%
      arrange(desc(.abs)) %>% select(-.abs)
    .write(wb, paste0(db, " \u2014 Log2FC GraphPad"), pivot_l2, "E67E22",
           note = paste0(
             "GraphPad Grouped table. ",
             "Rows = pathways, Columns = comparisons. ",
             "Values = log2(fold change) vs Control. Sorted by max absolute change."
           ))

    # -log10(p) GraphPad pivot
    pivot_p <- all_combined %>%
      select(Comparison, Pathway, Neg_Log10_P) %>%
      pivot_wider(names_from = Comparison, values_from = Neg_Log10_P,
                  values_fill = 0) %>%
      mutate(.mx = apply(across(-Pathway), 1, max, na.rm = TRUE)) %>%
      arrange(desc(.mx)) %>% select(-.mx)
    .write(wb, paste0(db, " \u2014 NegLog10P GraphPad"), pivot_p, "C0392B",
           note = paste0(
             "GraphPad Grouped table. ",
             "Rows = pathways, Columns = comparisons. ",
             "Values = -log10(p-value). Higher = more significant."
           ))

    # Fold enrichment pivot
    pivot_fe <- all_combined %>%
      select(Comparison, Pathway, Fold_Enrichment) %>%
      pivot_wider(names_from = Comparison, values_from = Fold_Enrichment,
                  values_fill = 1) %>%
      mutate(.mx = apply(across(-Pathway), 1, max, na.rm = TRUE)) %>%
      arrange(desc(.mx)) %>% select(-.mx)
    .write(wb, paste0(db, " \u2014 FoldEnrich GraphPad"), pivot_fe, "2E86C1",
           note = paste0(
             "GraphPad Grouped table. ",
             "Values = median(condition) / median(Control). >1 = enriched."
           ))

    # N compounds pivot
    pivot_n <- all_combined %>%
      select(Comparison, Pathway, N_Compounds) %>%
      pivot_wider(names_from = Comparison, values_from = N_Compounds,
                  values_fill = 0L)
    .write(wb, paste0(db, " \u2014 N Compounds"), pivot_n, "5D9B3F",
           note = "Compound counts per pathway per comparison (reference / error bar data).")
  }

  saveWorkbook(wb, out_file, overwrite = TRUE)
  safe_message("Saved Excel: ", out_file)
}

# =============================================================================
# ── RUN PIPELINE ──────────────────────────────────────────────────────────────
# =============================================================================

all_results <- list()   # will hold list(KEGG = ..., HMDB = ...)

# ── KEGG ──────────────────────────────────────────────────────────────────────
if (RUN_DATABASE %in% c("KEGG", "BOTH")) {

  safe_message("\n===== KEGG PIPELINE =====")
  kegg_long <- map_kegg(df, NAME_COL, API_PAUSE_SEC)

  safe_message("Running KEGG differential enrichment tests...")
  kegg_results <- compute_differential(
    kegg_long, sample_cols, col_conds, non_control, control_cols
  )
  all_results[["KEGG"]] <- kegg_results

  safe_message("Plotting KEGG dot plot...")
  make_dotplot(
    kegg_results, "KEGG",
    top_n   = TOP_N_PLOT,
    out_png = OUT_PNG_KEGG,
    out_svg = OUT_SVG_KEGG
  )

  safe_message("KEGG enrichment results:")
  for (cond in names(kegg_results)) {
    cat("\n---", cond, "vs Control ---\n")
    print(kegg_results[[cond]] %>%
      select(pathway_name, n_compounds, fold_enrichment, log2_fold, pvalue) %>%
      mutate(across(c(fold_enrichment, log2_fold), round, 3),
             pvalue = signif(pvalue, 3)),
      n = 20)
  }
}

# ── HMDB ──────────────────────────────────────────────────────────────────────
if (RUN_DATABASE %in% c("HMDB", "BOTH")) {

  safe_message("\n===== HMDB PIPELINE =====")
  hmdb_long <- map_hmdb(df, NAME_COL, API_PAUSE_SEC)

  safe_message("Running HMDB differential enrichment tests...")
  hmdb_results <- compute_differential(
    hmdb_long, sample_cols, col_conds, non_control, control_cols
  )
  all_results[["HMDB"]] <- hmdb_results

  safe_message("Plotting HMDB dot plot...")
  make_dotplot(
    hmdb_results, "HMDB",
    top_n   = TOP_N_PLOT,
    out_png = OUT_PNG_HMDB,
    out_svg = OUT_SVG_HMDB
  )

  safe_message("HMDB enrichment results:")
  for (cond in names(hmdb_results)) {
    cat("\n---", cond, "vs Control ---\n")
    print(hmdb_results[[cond]] %>%
      select(pathway_name, n_compounds, fold_enrichment, log2_fold, pvalue) %>%
      mutate(across(c(fold_enrichment, log2_fold), round, 3),
             pvalue = signif(pvalue, 3)),
      n = 20)
  }
}

# ── Excel export ───────────────────────────────────────────────────────────────
if (!is.null(OUT_EXCEL) && length(all_results) > 0) {
  safe_message("\n===== EXCEL EXPORT =====")
  export_excel(all_results, OUT_EXCEL)
}

safe_message("\nAll done. Files written:")
if (RUN_DATABASE %in% c("KEGG","BOTH")) {
  if (!is.null(OUT_PNG_KEGG)) message("  PNG : ", OUT_PNG_KEGG)
  if (!is.null(OUT_SVG_KEGG)) message("  SVG : ", OUT_SVG_KEGG)
}
if (RUN_DATABASE %in% c("HMDB","BOTH")) {
  if (!is.null(OUT_PNG_HMDB)) message("  PNG : ", OUT_PNG_HMDB)
  if (!is.null(OUT_SVG_HMDB)) message("  SVG : ", OUT_SVG_HMDB)
}
if (!is.null(OUT_EXCEL)) message("  XLSX: ", OUT_EXCEL)
