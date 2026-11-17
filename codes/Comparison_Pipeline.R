# ==============================================================================
# Spatial Proteomics
# User-selected GF vs SPF SINGLE Joint-MBR Comparison Pipeline
#
# INPUT
#   ONE DIA-NN report.parquet generated from ALL retained RAW files analyzed
#   together with MBR ON:
#
#     GF1-2, GF1-3, GF1-4, GF1-6
#     SPF1-2, SPF1-3, SPF1-4, SPF1-6
#
# ANALYSIS
#   - q-value filtering
#   - cRAP / mixed-cRAP audit and removal before downstream matrices
#   - automatic assignment of every ROI to one of the selected sets
#   - detected-protein and completeness summaries
#   - GF/SPF ROI QC, PCA, correlation and clustering
#   - main GF-vs-SPF protein statistics using SET-LEVEL summaries
#
# STATISTICAL UNIT
#   ROI-level tests are descriptive QC only. Main protein-level inference uses
#   one summary value per selected set: GF set values vs SPF set values.
#
# cRAP HANDLING
#   PURE_CRAP  -> removed
#   MIXED_CRAP -> ambiguous non-cRAP + cRAP group; removed from main matrix and
#                 saved separately for audit
#   CLEAN      -> retained
# ==============================================================================


# ============================================================================== 
# 1. Shiny/standalone entry point
# ============================================================================== 

COMPARISON_REQUIRED_PACKAGES <- c(
  "DBI", "duckdb", "dplyr", "tidyr", "tibble", "ggplot2", "pheatmap", "patchwork"
)

# Source this file, then call run_comparison(). Nothing is executed on source.
run_comparison <- function(report_path, output_dir, gf_set_patterns,
                           spf_set_patterns, db_memory_limit = "1GB",
                           png_dpi = 600L, excluded_rois = NULL) {
  missing_packages <- COMPARISON_REQUIRED_PACKAGES[
    !vapply(COMPARISON_REQUIRED_PACKAGES, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Install comparison packages first: install.packages(c(",
         paste(sprintf('"%s"', missing_packages), collapse = ", "), "))",
         call. = FALSE)
  }
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(tibble)
    library(ggplot2)
    library(pheatmap)
    library(patchwork)
  })

  normalize_patterns <- function(value, condition) {
    original_names <- names(value)
    value <- as.character(value)
    if (!length(value) || is.null(original_names) ||
        length(original_names) != length(value) || anyNA(value) ||
        anyNA(original_names) || any(!nzchar(trimws(value))) ||
        any(!nzchar(trimws(original_names)))) {
      stop(condition, " set patterns must be a non-empty named character vector.",
           call. = FALSE)
    }
    names(value) <- original_names
    names(value) <- toupper(gsub("_", "-", trimws(names(value)), fixed = TRUE))
    value <- toupper(gsub("_", "-", trimws(value), fixed = TRUE))
    expected <- paste0("^", condition, "[0-9]+(?:-[0-9]+)?$")
    if (any(!grepl(expected, names(value), perl = TRUE)) ||
        any(!grepl(expected, value, perl = TRUE))) {
      stop(condition, " set names must look like ", condition, "1-2.", call. = FALSE)
    }
    if (anyDuplicated(names(value)) || anyDuplicated(value)) {
      stop(condition, " set patterns must be unique.", call. = FALSE)
    }
    value
  }

  GF_SET_PATTERNS <- normalize_patterns(gf_set_patterns, "GF")
  SPF_SET_PATTERNS <- normalize_patterns(spf_set_patterns, "SPF")
  if (length(GF_SET_PATTERNS) < 3L || length(SPF_SET_PATTERNS) < 3L) {
    stop("Select at least three GF sets and three SPF sets for set-level inference.",
         call. = FALSE)
  }

  empty_exclusions <- data.frame(
    Dataset_ID = character(), ROI = character(), ROI_Number = integer(),
    Final_QC_Action = character(), Exclusion_Source = character(),
    Original_Run = character(), stringsAsFactors = FALSE
  )
  normalize_exclusions <- function(value) {
    if (is.null(value) || (is.data.frame(value) && !nrow(value))) {
      return(empty_exclusions)
    }
    if (!is.data.frame(value) || !"Dataset_ID" %in% names(value)) {
      stop("excluded_rois must be a data.frame containing Dataset_ID.",
           call. = FALSE)
    }
    value$Dataset_ID <- toupper(gsub("_", "-", trimws(
      as.character(value$Dataset_ID)), fixed = TRUE))
    if (!"ROI_Number" %in% names(value)) {
      if (!"ROI" %in% names(value)) {
        stop("excluded_rois must contain ROI_Number or ROI.", call. = FALSE)
      }
      value$ROI_Number <- suppressWarnings(as.integer(sub(
        "^ROI0*([0-9]+)$", "\\1", toupper(trimws(as.character(value$ROI))),
        perl = TRUE
      )))
    } else {
      value$ROI_Number <- suppressWarnings(as.integer(value$ROI_Number))
    }
    if (anyNA(value$Dataset_ID) || any(!nzchar(value$Dataset_ID)) ||
        anyNA(value$ROI_Number) || any(value$ROI_Number < 1L)) {
      stop("excluded_rois contains an invalid Dataset_ID or ROI number.",
           call. = FALSE)
    }
    value$ROI <- paste0("ROI", value$ROI_Number)
    for (column in c("Final_QC_Action", "Exclusion_Source", "Original_Run")) {
      if (!column %in% names(value)) value[[column]] <- NA_character_
      value[[column]] <- as.character(value[[column]])
    }
    selected_sets <- c(names(GF_SET_PATTERNS), names(SPF_SET_PATTERNS))
    value <- value[value$Dataset_ID %in% selected_sets, c(
      "Dataset_ID", "ROI", "ROI_Number", "Final_QC_Action",
      "Exclusion_Source", "Original_Run"
    ), drop = FALSE]
    empty_run <- is.na(value$Original_Run) |
      !nzchar(trimws(ifelse(is.na(value$Original_Run), "", value$Original_Run)))
    value$Original_Run[empty_run] <- NA_character_
    unique(value)
  }
  EXCLUDED_ROIS <- normalize_exclusions(excluded_rois)

  REPORT_FILE <- normalizePath(report_path, winslash = "/", mustWork = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  OUTPUT_DIR <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  WORK_DIR <- dirname(REPORT_FILE)
  old_wd <- getwd()
  setwd(WORK_DIR)
  on.exit(setwd(old_wd), add = TRUE)

  TOTAL_SET_COUNT <- length(GF_SET_PATTERNS) + length(SPF_SET_PATTERNS)
  COMPARISON_LABEL <- paste0(
    length(GF_SET_PATTERNS), " GF sets vs ", length(SPF_SET_PATTERNS), " SPF sets"
  )


# ============================================================================== 
# 2. ANALYSIS SETTINGS
# ============================================================================== 


# ------------------------------------------------------------------------------
# Protein completeness
#
# "EACH_SET":
#   protein must pass >=70% completeness in EACH selected set.
#   This is the strictest and is the default.
#
# "EACH_CONDITION":
#   protein must pass >=70% completeness across all GF ROIs and separately
#   across all SPF ROIs.
# ------------------------------------------------------------------------------

COMPLETENESS_CUTOFF <- 0.70

CORE_COMPLETENESS_MODE <- "EACH_CONDITION"


# ------------------------------------------------------------------------------
# DIA-NN confidence filtering
# ------------------------------------------------------------------------------

APPLY_Q_FILTER <- TRUE
Q_VALUE_CUTOFF <- 0.01

Q_COLUMN_PRIORITY <- c(
  "Global.PG.Q.Value",
  "PG.Q.Value",
  "Q.Value"
)


# ------------------------------------------------------------------------------
# DIA-NN abundance column
#
# PG.MaxLFQ is preferred.
# Both GF and SPF reports MUST resolve to the same abundance column.
# ------------------------------------------------------------------------------

ABUNDANCE_COLUMN_PRIORITY <- c(
  "PG.MaxLFQ",
  "PG.Quantity",
  "Protein.Group.Quantity",
  "Precursor.Quantity"
)


# ------------------------------------------------------------------------------
# cRAP filtering
# ------------------------------------------------------------------------------

REMOVE_CRAP <- TRUE

# Mixed mouse/non-cRAP + cRAP groups are ambiguous.
# TRUE = exclude from main quantitative analysis but save separately.
REMOVE_MIXED_CRAP_GROUPS <- TRUE

CRAP_PREFIXES <- c(
  "CON__",
  "CON_",
  "CONTAM_",
  "CONTAMINANT_",
  "CRAP_",
  "CRAP|",
  "CONTAM|"
)

USE_BUILTIN_GPM_CRAP <- TRUE

# Additional exact accessions if needed.
CRAP_EXACT_ACCESSIONS <- character(0)


# ------------------------------------------------------------------------------
# GF-vs-SPF protein comparison
#
# Statistics are performed on the selected GF set summaries vs SPF set summaries,
# NOT on all individual ROIs as independent biological replicates.
# ------------------------------------------------------------------------------

MIN_SETS_PER_CONDITION_FOR_TEST <- 3

DIFF_FDR_CUTOFF <- 0.05

# abs(log2FC) >= 1 = >=2-fold difference
DIFF_MIN_ABS_LOG2FC <- 1.0


# ------------------------------------------------------------------------------
# PCA / heatmap
# ------------------------------------------------------------------------------

PCA_SCALE <- TRUE

TOP_N_VARIABLE_PROTEINS <- 500

TOP_N_DIFFERENTIAL_PROTEINS <- 50


# ------------------------------------------------------------------------------
# Figure settings
# ------------------------------------------------------------------------------

if (!is.numeric(png_dpi) || length(png_dpi) != 1L ||
    !is.finite(png_dpi) || png_dpi < 72 || png_dpi > 1200) {
  stop("png_dpi must be between 72 and 1200.", call. = FALSE)
}
PNG_DPI <- as.integer(png_dpi)

PLOT_WIDTH <- 10
PLOT_HEIGHT <- 7

HEATMAP_WIDTH <- 12
HEATMAP_HEIGHT <- 10

ROI_COR_SHOW_NAMES <- FALSE

VARIABLE_HEATMAP_SHOW_ROW_NAMES <- FALSE
VARIABLE_HEATMAP_SHOW_COL_NAMES <- FALSE

DIFF_HEATMAP_SHOW_ROW_NAMES <- TRUE

ALL_PLOTS_PDF_WIDTH <- 12
ALL_PLOTS_PDF_HEIGHT <- 9


# ==============================================================================
# 3. Helper functions
# ==============================================================================

output_path <- function(file_name) {
  file.path(
    OUTPUT_DIR,
    file_name
  )
}


first_existing_column <- function(
    priority_vector,
    available_columns
) {

  hit <- priority_vector[
    priority_vector %in%
      available_columns
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}

set_pattern_regex <- function(pattern) {
  canonical <- toupper(gsub("_", "-", trimws(as.character(pattern)), fixed = TRUE))
  token <- gsub("-", "[-_]", canonical, fixed = TRUE)
  paste0("(^|[^A-Z0-9])", token, "($|[^A-Z0-9-])")
}

set_pattern_matches <- function(pattern, values) {
  grepl(set_pattern_regex(pattern), toupper(as.character(values)), perl = TRUE)
}


safe_median <- function(x) {

  x <- suppressWarnings(
    as.numeric(x)
  )

  x <- x[
    is.finite(x) &
      x > 0
  ]

  if (length(x) == 0) {
    return(NA_real_)
  }

  median(
    x,
    na.rm = TRUE
  )
}


safe_sum <- function(x) {

  x <- suppressWarnings(
    as.numeric(x)
  )

  x <- x[
    is.finite(x) &
      x > 0
  ]

  if (length(x) == 0) {
    return(NA_real_)
  }

  sum(
    x,
    na.rm = TRUE
  )
}


median_impute_by_protein <- function(mat) {

  imputed <- t(
    apply(
      mat,
      1,
      function(x) {

        finite_x <- x[
          is.finite(x)
        ]

        if (length(finite_x) == 0) {
          return(x)
        }

        replacement <- median(
          finite_x,
          na.rm = TRUE
        )

        x[
          !is.finite(x)
        ] <- replacement

        x
      }
    )
  )

  rownames(imputed) <- rownames(mat)
  colnames(imputed) <- colnames(mat)

  imputed
}


save_ggplot_png <- function(
    plot_object,
    file_name,
    width = PLOT_WIDTH,
    height = PLOT_HEIGHT
) {

  ggplot2::ggsave(
    filename = file_name,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = PNG_DPI,
    bg = "white"
  )
}


save_pheatmap <- function(
    ph,
    file_stem,
    width = HEATMAP_WIDTH,
    height = HEATMAP_HEIGHT
) {

  if (is.null(ph)) {
    return(
      invisible(FALSE)
    )
  }

  grDevices::pdf(
    paste0(
      file_stem,
      ".pdf"
    ),
    width = width,
    height = height
  )

  grid::grid.newpage()
  grid::grid.draw(
    ph$gtable
  )

  grDevices::dev.off()


  grDevices::png(
    paste0(
      file_stem,
      ".png"
    ),
    width = width,
    height = height,
    units = "in",
    res = PNG_DPI
  )

  grid::grid.newpage()
  grid::grid.draw(
    ph$gtable
  )

  grDevices::dev.off()

  invisible(TRUE)
}


format_p <- function(p) {

  if (!is.finite(p)) {
    return("NA")
  }

  format.pval(
    p,
    digits = 3,
    eps = 0.001
  )
}


make_metric_plot <- function(
    data,
    metric_column,
    y_label,
    title
) {

  test_data <- data[
    is.finite(
      data[[metric_column]]
    ) &
      !is.na(
        data$Condition
      ),
    ,
    drop = FALSE
  ]

  kw_p <- NA_real_

  if (
    nrow(test_data) >= 4 &&
      length(
        unique(
          test_data$Condition
        )
      ) >= 2
  ) {

    kw_p <- tryCatch(
      kruskal.test(
        test_data[[metric_column]] ~
          test_data$Condition
      )$p.value,
      error = function(e) NA_real_
    )
  }


  p <- ggplot(
    test_data,
    aes(
      x = Condition,
      y = .data[[metric_column]],
      color = Condition
    )
  ) +
    geom_boxplot(
      outlier.shape = NA,
      color = "black"
    ) +
    geom_jitter(
      width = 0.15,
      alpha = 0.70,
      size = 1.4
    ) +
    theme_bw() +
    theme(
      legend.position = "none"
    ) +
    labs(
      title = title,
      subtitle = paste0(
        "ROI-level descriptive comparison; Kruskal-Wallis p = ",
        format_p(
          kw_p
        )
      ),
      x = NULL,
      y = y_label
    )


  list(
    plot = p,
    p_value = kw_p
  )
}


# ==============================================================================
# 4. Built-in GPM cRAP accession reference
#
# Same reference logic as the GF-only pipeline.
# The filtering logic below distinguishes pure cRAP from mixed ambiguous groups.
# ==============================================================================

GPM_CRAP_ACCESSIONS <- c(

  # Laboratory proteins / reagents
  "P02769", "P00766", "P00767", "P00711", "Q7M135",
  "P00792", "P00791", "Q10735", "P30879", "P04188",
  "P00760", "Q29463", "P00761", "P02662", "P02663",
  "P02666", "P02668", "P42212", "P81054",

  # Dust / contact proteins
  "P04745", "P13645", "P35527", "P04264", "P35908",
  "Q15323", "Q14532", "O76011", "Q92764", "O76013",
  "O76014", "O76015", "O76009", "Q14525", "Q14533",
  "Q9NSB4", "P78385", "Q9NSB2", "P78386", "O43790",
  "O77727", "P02534", "P25690", "P02539", "P15241",
  "P25691", "P02444", "P02445", "P02443", "P02441",
  "Q02958", "P02438", "P02439", "P02440", "P08131",
  "P26372", "O82803", "P15252",

  # Molecular-weight / MS standards
  "P00004", "P00921", "P00330", "P00883", "P00698",
  "P68082", "P01012", "P00722", "P00366",

  # Sigma-Aldrich Universal Protein Standard proteins
  "P02768", "P01008", "P08758", "P61769", "P55957",
  "P00915", "P00918", "P04040", "P07339", "P08311",
  "P01031", "P02741", "P00167", "P99999", "P01133",
  "P05413", "P06396", "P08263", "P09211", "P69905",
  "P68871", "P01344", "P10145", "P06732", "P00709",
  "P41159", "P61626", "P02144", "Q15843", "P15559",
  "P16083", "P01127", "P62937", "Q06830", "P01112",
  "P02753", "P00441", "P63165", "P12081", "P10636",
  "P10599", "P01375", "P02787", "P02788", "P51965",
  "O00762", "P63279", "P62979",

  # Common viral contaminant
  "P32503"
)


crap_accessions_upper <- unique(
  toupper(
    c(
      if (USE_BUILTIN_GPM_CRAP) {
        GPM_CRAP_ACCESSIONS
      } else {
        character(0)
      },
      CRAP_EXACT_ACCESSIONS
    )
  )
)

crap_prefixes_upper <- unique(
  toupper(
    CRAP_PREFIXES
  )
)


# ------------------------------------------------------------------------------
# Determine whether ONE protein member is cRAP.
#
# Examples handled:
#   P00698
#   sp|P00698|LYSC_CHICK
#   CON__P00698
# ------------------------------------------------------------------------------

member_is_crap <- function(member) {

  if (
    is.na(member) ||
      !nzchar(
        trimws(
          as.character(member)
        )
      )
  ) {
    return(FALSE)
  }

  member <- toupper(
    trimws(
      as.character(member)
    )
  )


  prefix_hit <- any(
    vapply(
      crap_prefixes_upper,
      function(prefix) {
        startsWith(
          member,
          prefix
        )
      },
      logical(1)
    )
  )


  pipe_tokens <- unlist(
    strsplit(
      member,
      split = "\\|",
      perl = TRUE
    ),
    use.names = FALSE
  )

  pipe_tokens <- trimws(
    pipe_tokens
  )

  exact_hit <- any(
    pipe_tokens %in%
      crap_accessions_upper
  )


  prefix_hit || exact_hit
}


# ------------------------------------------------------------------------------
# Classify one DIA-NN protein group as:
#   CLEAN
#   PURE_CRAP
#   MIXED_CRAP
#
# P00687 + P00698 is MIXED_CRAP, not "P00687 is cRAP".
# ------------------------------------------------------------------------------

classify_crap_group <- function(
    protein_group,
    protein_ids = NA_character_
) {

  strings_to_use <- c(
    protein_group,
    protein_ids
  )

  strings_to_use <- strings_to_use[
    !is.na(
      strings_to_use
    ) &
      nzchar(
        trimws(
          strings_to_use
        )
      )
  ]


  if (length(strings_to_use) == 0) {
    return("CLEAN")
  }


  members <- unique(
    unlist(
      strsplit(
        strings_to_use,
        split = "[;,]",
        perl = TRUE
      ),
      use.names = FALSE
    )
  )

  members <- trimws(
    members
  )

  members <- members[
    nzchar(
      members
    )
  ]


  if (length(members) == 0) {
    return("CLEAN")
  }


  flags <- vapply(
    members,
    member_is_crap,
    logical(1)
  )


  if (
    any(flags) &&
      all(flags)
  ) {
    return("PURE_CRAP")
  }


  if (
    any(flags) &&
      any(!flags)
  ) {
    return("MIXED_CRAP")
  }


  "CLEAN"
}


# ==============================================================================
# 5. Process one condition
# ==============================================================================

process_condition <- function(
    report_all,
    condition_name,
    set_patterns
) {

  cat("\n=====================================================\n")
  cat("PROCESSING CONDITION:", condition_name, "\n")
  cat("=====================================================\n")


  # Work on an independent copy of the ONE joint-MBR report.
  report <- as.data.frame(
    report_all,
    check.names = FALSE
  )


  RUN_COL <- first_existing_column(
    c(
      "Run",
      "File.Name",
      "File"
    ),
    names(report)
  )

  PROTEIN_COL <- first_existing_column(
    c(
      "Protein.Group",
      "Protein.Ids"
    ),
    names(report)
  )

  PROTEIN_IDS_COL <- first_existing_column(
    c(
      "Protein.Ids",
      "Protein.Group"
    ),
    names(report)
  )

  ABUNDANCE_COL <- first_existing_column(
    ABUNDANCE_COLUMN_PRIORITY,
    names(report)
  )

  Q_COL <- first_existing_column(
    Q_COLUMN_PRIORITY,
    names(report)
  )


  if (is.na(RUN_COL)) {
    stop(
      paste0(
        condition_name,
        ": DIA-NN Run column not found."
      )
    )
  }

  if (is.na(PROTEIN_COL)) {
    stop(
      paste0(
        condition_name,
        ": Protein.Group / Protein.Ids not found."
      )
    )
  }

  if (is.na(ABUNDANCE_COL)) {
    stop(
      paste0(
        condition_name,
        ": supported abundance column not found."
      )
    )
  }

  condition_keep <- Reduce(
    `|`,
    lapply(unname(set_patterns), function(pattern) {
      set_pattern_matches(pattern, report[[RUN_COL]])
    })
  )
  condition_keep[is.na(condition_keep)] <- FALSE
  report <- report[condition_keep, , drop = FALSE]
  if (!nrow(report)) {
    stop(condition_name, ": no rows matched the selected set patterns.", call. = FALSE)
  }


  cat("Run column      :", RUN_COL, "\n")
  cat("Protein column  :", PROTEIN_COL, "\n")
  cat("Protein IDs col :", PROTEIN_IDS_COL, "\n")
  cat("Abundance column:", ABUNDANCE_COL, "\n")
  cat(
    "Q-value column  :",
    ifelse(
      is.na(Q_COL),
      "not found",
      Q_COL
    ),
    "\n"
  )


  # ---------------------------------------------------------------------------
  # Q-value filtering
  # ---------------------------------------------------------------------------

  n_before_q <- nrow(
    report
  )


  if (
    APPLY_Q_FILTER &&
      !is.na(Q_COL)
  ) {

    q_values <- suppressWarnings(
      as.numeric(
        report[[Q_COL]]
      )
    )

    keep_q <-
      is.finite(q_values) &
      q_values <= Q_VALUE_CUTOFF

    report <- report[
      keep_q,
      ,
      drop = FALSE
    ]
  }


  n_after_q <- nrow(
    report
  )


  # ---------------------------------------------------------------------------
  # cRAP classification using UNIQUE Protein.Group + Protein.Ids pairs
  # ---------------------------------------------------------------------------

  pg_values <- as.character(
    report[[PROTEIN_COL]]
  )

  pid_values <- as.character(
    report[[PROTEIN_IDS_COL]]
  )


  crap_key_df <- unique(
    data.frame(
      Protein.Group.Internal = pg_values,
      Protein.Ids.Internal = pid_values,
      stringsAsFactors = FALSE
    )
  )


  crap_key_df$cRAP_Class <- mapply(
    classify_crap_group,
    protein_group =
      crap_key_df$Protein.Group.Internal,
    protein_ids =
      crap_key_df$Protein.Ids.Internal,
    USE.NAMES = FALSE
  )


  row_key <- paste(
    pg_values,
    pid_values,
    sep = "\t"
  )

  audit_key <- paste(
    crap_key_df$Protein.Group.Internal,
    crap_key_df$Protein.Ids.Internal,
    sep = "\t"
  )


  report$cRAP.Class.Internal <-
    crap_key_df$cRAP_Class[
      match(
        row_key,
        audit_key
      )
    ]


  pure_crap_audit <- crap_key_df %>%
    filter(
      cRAP_Class == "PURE_CRAP"
    )

  mixed_crap_audit <- crap_key_df %>%
    filter(
      cRAP_Class == "MIXED_CRAP"
    )


  write.csv(
    pure_crap_audit,
    output_path(
      paste0(
        "00_",
        condition_name,
        "_Pure_cRAP_Removed.csv"
      )
    ),
    row.names = FALSE
  )


  write.csv(
    mixed_crap_audit,
    output_path(
      paste0(
        "00_",
        condition_name,
        "_Mixed_cRAP_Ambiguous_Groups.csv"
      )
    ),
    row.names = FALSE
  )


  remove_flag <- rep(
    FALSE,
    nrow(report)
  )


  if (REMOVE_CRAP) {

    remove_flag <-
      report$cRAP.Class.Internal ==
      "PURE_CRAP"

    if (REMOVE_MIXED_CRAP_GROUPS) {

      remove_flag <-
        remove_flag |
        report$cRAP.Class.Internal ==
        "MIXED_CRAP"
    }
  }


  n_before_crap <- nrow(
    report
  )

  n_pure_crap_rows <- sum(
    report$cRAP.Class.Internal ==
      "PURE_CRAP",
    na.rm = TRUE
  )

  n_mixed_rows <- sum(
    report$cRAP.Class.Internal ==
      "MIXED_CRAP",
    na.rm = TRUE
  )


  report <- report[
    !remove_flag,
    ,
    drop = FALSE
  ]


  n_after_crap <- nrow(
    report
  )


  # Hard audit
  if (
    REMOVE_CRAP &&
      any(
        report$cRAP.Class.Internal ==
          "PURE_CRAP",
        na.rm = TRUE
      )
  ) {
    stop(
      paste0(
        condition_name,
        ": pure cRAP remained after filtering."
      )
    )
  }


  if (
    REMOVE_CRAP &&
      REMOVE_MIXED_CRAP_GROUPS &&
      any(
        report$cRAP.Class.Internal ==
          "MIXED_CRAP",
        na.rm = TRUE
      )
  ) {
    stop(
      paste0(
        condition_name,
        ": mixed cRAP groups remained after filtering."
      )
    )
  }


  # ---------------------------------------------------------------------------
  # Assign runs to the selected sets
  # ---------------------------------------------------------------------------

  all_runs <- unique(
    as.character(
      report[[RUN_COL]]
    )
  )


  match_matrix <- sapply(
    set_patterns,
    function(pattern) {
      set_pattern_matches(pattern, all_runs)
    }
  )


  if (is.null(dim(match_matrix))) {

    match_matrix <- matrix(
      match_matrix,
      ncol = 1
    )

    colnames(
      match_matrix
    ) <- names(
      set_patterns
    )
  }


  n_matches <- rowSums(
    match_matrix
  )


  if (any(n_matches > 1)) {

    ambiguous_runs <- all_runs[
      n_matches > 1
    ]

    stop(
      paste0(
        condition_name,
        ": some Run names matched multiple set patterns:\n",
        paste(
          ambiguous_runs,
          collapse = "\n"
        )
      )
    )
  }


  assigned_set <- rep(
    NA_character_,
    length(
      all_runs
    )
  )


  for (i in seq_along(all_runs)) {

    if (n_matches[i] == 1) {

      assigned_set[i] <- names(
        set_patterns
      )[
        which(
          match_matrix[i, ]
        )
      ]
    }
  }


  run_map <- data.frame(
    Run.Original = all_runs,
    Set = assigned_set,
    stringsAsFactors = FALSE
  )


  run_map$Run.Short <- basename(
    gsub(
      "\\\\",
      "/",
      run_map$Run.Original
    )
  )


  run_map$Run.Short <- make.unique(
    run_map$Run.Short
  )


  run_map <- run_map %>%
    filter(
      !is.na(Set)
    )


  expected_sets <- names(
    set_patterns
  )

  observed_sets <- unique(
    run_map$Set
  )

  missing_sets <- setdiff(
    expected_sets,
    observed_sets
  )


  if (length(missing_sets) > 0) {

    stop(
      paste0(
        condition_name,
        ": no DIA-NN runs found for set(s): ",
        paste(
          missing_sets,
          collapse = ", "
        ),
        "\nCheck ",
        condition_name,
        "_SET_PATTERNS."
      )
    )
  }


  run_map$Condition <- condition_name

  run_map$Sample <- paste0(
    condition_name,
    "__",
    run_map$Run.Short
  )


  write.csv(
    run_map,
    output_path(
      paste0(
        "00_",
        condition_name,
        "_Run_Set_Assignment.csv"
      )
    ),
    row.names = FALSE
  )


  # ---------------------------------------------------------------------------
  # Restrict to target sets
  # ---------------------------------------------------------------------------

  target_report <- report[
    as.character(
      report[[RUN_COL]]
    ) %in%
      run_map$Run.Original,
    ,
    drop = FALSE
  ]


  target_report$Sample.Internal <-
    run_map$Sample[
      match(
        as.character(
          target_report[[RUN_COL]]
        ),
        run_map$Run.Original
      )
    ]


  target_report$Set.Internal <-
    run_map$Set[
      match(
        as.character(
          target_report[[RUN_COL]]
        ),
        run_map$Run.Original
      )
    ]


  target_report$Condition.Internal <-
    condition_name


  target_report$Protein.Internal <-
    as.character(
      target_report[[PROTEIN_COL]]
    )


  target_report$Abundance.Internal <-
    suppressWarnings(
      as.numeric(
        target_report[[ABUNDANCE_COL]]
      )
    )


  target_report <- target_report %>%
    filter(
      !is.na(
        Sample.Internal
      ),
      !is.na(
        Set.Internal
      ),
      !is.na(
        Protein.Internal
      ),
      Protein.Internal != "",
      is.finite(
        Abundance.Internal
      ),
      Abundance.Internal > 0
    )


  # ---------------------------------------------------------------------------
  # Protein x run aggregation
  # ---------------------------------------------------------------------------

  if (ABUNDANCE_COL == "Precursor.Quantity") {

    warning(
      paste0(
        condition_name,
        ": PG.MaxLFQ / PG.Quantity not found. ",
        "Precursor.Quantity will be summed to protein level."
      )
    )


    protein_long <- target_report %>%
      group_by(
        Protein.Group = Protein.Internal,
        Sample = Sample.Internal,
        Set = Set.Internal,
        Condition = Condition.Internal
      ) %>%
      summarise(
        Abundance = safe_sum(
          Abundance.Internal
        ),
        .groups = "drop"
      )

  } else {

    protein_long <- target_report %>%
      group_by(
        Protein.Group = Protein.Internal,
        Sample = Sample.Internal,
        Set = Set.Internal,
        Condition = Condition.Internal
      ) %>%
      summarise(
        Abundance = safe_median(
          Abundance.Internal
        ),
        .groups = "drop"
      )
  }


  sample_info <- run_map %>%
    transmute(
      Sample,
      Condition,
      Set,
      Run.Original,
      Run.Short
    )


  filter_summary <- data.frame(
    Condition = condition_name,
    Input_Rows = n_before_q,
    Rows_After_Q_Filter = n_after_q,
    Pure_cRAP_Rows_Detected = n_pure_crap_rows,
    Mixed_cRAP_Rows_Detected = n_mixed_rows,
    Rows_After_cRAP_Filter = n_after_crap,
    Number_Target_Runs = nrow(
      run_map
    ),
    Abundance_Column = ABUNDANCE_COL,
    Q_Column = ifelse(
      is.na(Q_COL),
      "not found",
      Q_COL
    ),
    stringsAsFactors = FALSE
  )


  list(
    protein_long = protein_long,
    sample_info = sample_info,
    filter_summary = filter_summary,
    abundance_col = ABUNDANCE_COL,
    q_col = Q_COL
  )
}


# ==============================================================================
# 6. Read ONE joint-MBR report, then process GF and SPF subsets
# ==============================================================================

if (!file.exists(REPORT_FILE)) {
  stop(
    paste0(
      "Joint-MBR report.parquet not found:\n",
      REPORT_FILE
    )
  )
}

selected_sets_manifest <- data.frame(
  Condition = c(rep("GF", length(GF_SET_PATTERNS)),
                rep("SPF", length(SPF_SET_PATTERNS))),
  Set = c(names(GF_SET_PATTERNS), names(SPF_SET_PATTERNS)),
  Run_Pattern = c(unname(GF_SET_PATTERNS), unname(SPF_SET_PATTERNS)),
  stringsAsFactors = FALSE
)
utils::write.csv(
  selected_sets_manifest,
  output_path("00_Selected_GF_SPF_Sets.csv"),
  row.names = FALSE
)
utils::write.csv(
  EXCLUDED_ROIS,
  output_path("00_QC_Excluded_ROIs_Applied.csv"),
  row.names = FALSE,
  na = ""
)

cat("\nReading joint-MBR DIA-NN report for selected GF/SPF sets:\n")
cat(REPORT_FILE, "\n")
cat("QC-excluded ROIs applied:", nrow(EXCLUDED_ROIS), "\n")

sql_work <- tempfile("comparison_sql_", tmpdir = OUTPUT_DIR)
dir.create(sql_work)
sql_spill <- file.path(sql_work, "spill")
dir.create(sql_spill)
sql_con <- DBI::dbConnect(
  duckdb::duckdb(), dbdir = file.path(sql_work, "comparison.duckdb"),
  bigint = "numeric",
  config = list(
    memory_limit = db_memory_limit, threads = "2",
    preserve_insertion_order = "false",
    temp_directory = normalizePath(sql_spill, winslash = "/", mustWork = TRUE)
  )
)
on.exit({
  if (!is.null(sql_con)) try(DBI::dbDisconnect(sql_con, shutdown = TRUE), silent = TRUE)
  unlink(sql_work, recursive = TRUE)
}, add = TRUE)
sql_path <- as.character(DBI::dbQuoteString(sql_con, REPORT_FILE))
schema <- DBI::dbGetQuery(
  sql_con, paste0("DESCRIBE SELECT * FROM read_parquet(", sql_path, ")")
)$column_name
LOAD_RUN_COL <- first_existing_column(c("Run", "File.Name", "File"), schema)
LOAD_PROTEIN_COL <- first_existing_column(c("Protein.Group", "Protein.Ids"), schema)
LOAD_PROTEIN_IDS_COL <- first_existing_column(c("Protein.Ids", "Protein.Group"), schema)
LOAD_ABUNDANCE_COL <- first_existing_column(ABUNDANCE_COLUMN_PRIORITY, schema)
LOAD_Q_COL <- first_existing_column(Q_COLUMN_PRIORITY, schema)
required_load_columns <- c(LOAD_RUN_COL, LOAD_PROTEIN_COL, LOAD_ABUNDANCE_COL)
if (anyNA(required_load_columns)) {
  stop("The report is missing a supported Run, protein, or abundance column.",
       call. = FALSE)
}
load_columns <- unique(na.omit(c(
  LOAD_RUN_COL, LOAD_PROTEIN_COL, LOAD_PROTEIN_IDS_COL,
  LOAD_ABUNDANCE_COL, LOAD_Q_COL
)))
selected_patterns <- unique(c(unname(GF_SET_PATTERNS), unname(SPF_SET_PATTERNS)))
run_identifier <- as.character(DBI::dbQuoteIdentifier(sql_con, LOAD_RUN_COL))
pattern_sql <- vapply(selected_patterns, function(pattern) {
  paste0(
    "regexp_matches(upper(CAST(", run_identifier, " AS VARCHAR)), ",
    as.character(DBI::dbQuoteString(sql_con, set_pattern_regex(pattern))), ")"
  )
}, character(1))
exclusion_predicates <- character()
exact_runs <- unique(EXCLUDED_ROIS$Original_Run[!is.na(EXCLUDED_ROIS$Original_Run)])
if (length(exact_runs)) {
  exclusion_predicates <- c(exclusion_predicates, vapply(exact_runs, function(run_name) {
    paste0(
      "CAST(", run_identifier, " AS VARCHAR) = ",
      as.character(DBI::dbQuoteString(sql_con, run_name))
    )
  }, character(1)))
}
fallback_exclusions <- EXCLUDED_ROIS[is.na(EXCLUDED_ROIS$Original_Run), , drop = FALSE]
if (nrow(fallback_exclusions)) {
  fallback_sql <- vapply(seq_len(nrow(fallback_exclusions)), function(i) {
    dataset_regex <- set_pattern_regex(fallback_exclusions$Dataset_ID[[i]])
    roi_regex <- paste0(
      "(^|[^A-Z0-9])ROI0*", fallback_exclusions$ROI_Number[[i]],
      "($|[^0-9])"
    )
    paste0(
      "(regexp_matches(upper(CAST(", run_identifier, " AS VARCHAR)), ",
      as.character(DBI::dbQuoteString(sql_con, dataset_regex)), ") AND ",
      "regexp_matches(upper(CAST(", run_identifier, " AS VARCHAR)), ",
      as.character(DBI::dbQuoteString(sql_con, roi_regex)), "))"
    )
  }, character(1))
  exclusion_predicates <- c(exclusion_predicates, fallback_sql)
}
exclusion_sql <- if (length(exclusion_predicates)) {
  paste0(" AND NOT (", paste(exclusion_predicates, collapse = " OR "), ")")
} else {
  ""
}
select_sql <- paste(
  as.character(DBI::dbQuoteIdentifier(sql_con, load_columns)), collapse = ", "
)
report_all <- DBI::dbGetQuery(
  sql_con,
  paste0(
    "SELECT ", select_sql, " FROM read_parquet(", sql_path, ") WHERE (",
    paste(pattern_sql, collapse = " OR "), ")", exclusion_sql
  )
)
DBI::dbDisconnect(sql_con, shutdown = TRUE)
sql_con <- NULL
unlink(sql_work, recursive = TRUE)
if (!nrow(report_all)) {
  stop("No report rows matched the selected GF/SPF sets.", call. = FALSE)
}
report_all <- as.data.frame(report_all, check.names = FALSE)

cat(
  "Rows:", nrow(report_all),
  "| Columns:", ncol(report_all),
  "\n"
)


gf <- process_condition(
  report_all = report_all,
  condition_name = "GF",
  set_patterns = GF_SET_PATTERNS
)

spf <- process_condition(
  report_all = report_all,
  condition_name = "SPF",
  set_patterns = SPF_SET_PATTERNS
)


# Both subsets originate from the same DIA-NN analysis and therefore should
# resolve to the same quantitative field.
if (
  gf$abundance_col !=
    spf$abundance_col
) {

  stop(
    paste0(
      "GF and SPF unexpectedly resolved to different abundance columns:\n",
      "GF  = ",
      gf$abundance_col,
      "\nSPF = ",
      spf$abundance_col
    )
  )
}


filter_summary <- bind_rows(
  gf$filter_summary,
  spf$filter_summary
)

write.csv(
  filter_summary,
  output_path(
    "00_GF_SPF_Filter_Summary.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 7. Combine GF + SPF protein-level data
# ==============================================================================

protein_long <- bind_rows(
  gf$protein_long,
  spf$protein_long
)


sample_info <- bind_rows(
  gf$sample_info,
  spf$sample_info
)


sample_info$Condition <- factor(
  sample_info$Condition,
  levels = c(
    "GF",
    "SPF"
  )
)


all_set_names <- c(
  names(
    GF_SET_PATTERNS
  ),
  names(
    SPF_SET_PATTERNS
  )
)

sample_info$Set <- factor(
  sample_info$Set,
  levels = all_set_names
)


rownames(
  sample_info
) <- sample_info$Sample


write.csv(
  sample_info,
  output_path(
    "01_All_Sample_Set_Information.csv"
  ),
  row.names = FALSE
)


write.csv(
  protein_long,
  output_path(
    "02_GF_SPF_Protein_Run_Long_Abundance.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 8. Build combined Protein x ROI matrix
# ==============================================================================

matrix_df <- protein_long %>%
  select(
    Protein.Group,
    Sample,
    Abundance
  ) %>%
  pivot_wider(
    names_from = Sample,
    values_from = Abundance
  )


protein_ids <- matrix_df$Protein.Group


linear_matrix <- as.matrix(
  matrix_df[
    ,
    setdiff(
      names(
        matrix_df
      ),
      "Protein.Group"
    ),
    drop = FALSE
  ]
)


rownames(
  linear_matrix
) <- protein_ids


storage.mode(
  linear_matrix
) <- "numeric"


missing_samples <- setdiff(
  sample_info$Sample,
  colnames(
    linear_matrix
  )
)


if (length(missing_samples) > 0) {

  stop(
    paste0(
      "Some assigned samples have no usable protein abundance data:\n",
      paste(
        missing_samples,
        collapse = "\n"
      )
    )
  )
}


linear_matrix <- linear_matrix[
  ,
  sample_info$Sample,
  drop = FALSE
]


linear_matrix[
  !is.finite(
    linear_matrix
  ) |
    linear_matrix <= 0
] <- NA_real_


log2_matrix <- log2(
  linear_matrix
)


write.csv(
  cbind(
    Protein.Group = rownames(
      linear_matrix
    ),
    as.data.frame(
      linear_matrix,
      check.names = FALSE
    )
  ),
  output_path(
    "03_GF_SPF_Protein_ROI_Linear_Matrix.csv"
  ),
  row.names = FALSE
)


write.csv(
  cbind(
    Protein.Group = rownames(
      log2_matrix
    ),
    as.data.frame(
      log2_matrix,
      check.names = FALSE
    )
  ),
  output_path(
    "04_GF_SPF_Protein_ROI_Log2_Matrix.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 9. Column lists for conditions and sets
# ==============================================================================

condition_columns <- list(
  GF = sample_info$Sample[
    sample_info$Condition == "GF"
  ],
  SPF = sample_info$Sample[
    sample_info$Condition == "SPF"
  ]
)


set_columns <- lapply(
  all_set_names,
  function(s) {

    sample_info$Sample[
      sample_info$Set == s
    ]
  }
)

names(
  set_columns
) <- all_set_names


# ==============================================================================
# 10. Detected protein and completeness counts by set
# ==============================================================================

set_summary_list <- list()

detected_by_set <- list()

completeness_by_set <- list()


for (s in all_set_names) {

  cols <- set_columns[[s]]

  m <- log2_matrix[
    ,
    cols,
    drop = FALSE
  ]


  detected <- rownames(m)[
    rowSums(
      is.finite(m)
    ) > 0
  ]


  completeness <- rowMeans(
    is.finite(m)
  )


  detected_by_set[[s]] <- detected

  completeness_by_set[[s]] <- completeness


  condition_this_set <- as.character(
    unique(
      sample_info$Condition[
        sample_info$Set == s
      ]
    )
  )


  set_summary_list[[s]] <- data.frame(
    Condition = condition_this_set,
    Set = s,
    N_ROIs = length(
      cols
    ),
    Detected_Proteins = length(
      detected
    ),
    Proteins_70pct_Completeness = sum(
      completeness >=
        COMPLETENESS_CUTOFF,
      na.rm = TRUE
    ),
    Median_Proteins_per_ROI = median(
      colSums(
        is.finite(m)
      ),
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}


set_summary <- bind_rows(
  set_summary_list
)


set_summary$Condition <- factor(
  set_summary$Condition,
  levels = c(
    "GF",
    "SPF"
  )
)


set_summary$Set <- factor(
  set_summary$Set,
  levels = all_set_names
)


write.csv(
  set_summary,
  output_path(
    "05_GF_SPF_Set_Summary.csv"
  ),
  row.names = FALSE
)


p_detected <- ggplot(
  set_summary,
  aes(
    x = Set,
    y = Detected_Proteins,
    fill = Condition
  )
) +
  geom_col() +
  geom_text(
    aes(
      label = Detected_Proteins
    ),
    vjust = -0.25,
    size = 3
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = "Detected Protein Groups in GF and SPF Sets",
    subtitle = "Final joint-MBR data within each condition",
    x = NULL,
    y = "Detected protein groups"
  )


p_complete70 <- ggplot(
  set_summary,
  aes(
    x = Set,
    y = Proteins_70pct_Completeness,
    fill = Condition
  )
) +
  geom_col() +
  geom_text(
    aes(
      label =
        Proteins_70pct_Completeness
    ),
    vjust = -0.25,
    size = 3
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = paste0(
      "Protein Groups with >= ",
      COMPLETENESS_CUTOFF * 100,
      "% Completeness by Set"
    ),
    x = NULL,
    y = "Protein groups"
  )


# ==============================================================================
# 11. GF vs SPF detection overlap
# ==============================================================================

gf_detected <- rownames(
  log2_matrix
)[
  rowSums(
    is.finite(
      log2_matrix[
        ,
        condition_columns$GF,
        drop = FALSE
      ]
    )
  ) > 0
]


spf_detected <- rownames(
  log2_matrix
)[
  rowSums(
    is.finite(
      log2_matrix[
        ,
        condition_columns$SPF,
        drop = FALSE
      ]
    )
  ) > 0
]


detected_common <- intersect(
  gf_detected,
  spf_detected
)

detected_gf_only <- setdiff(
  gf_detected,
  spf_detected
)

detected_spf_only <- setdiff(
  spf_detected,
  gf_detected
)


condition_overlap_summary <- data.frame(
  Category = c(
    "GF only",
    "GF + SPF",
    "SPF only"
  ),
  N_Proteins = c(
    length(
      detected_gf_only
    ),
    length(
      detected_common
    ),
    length(
      detected_spf_only
    )
  ),
  stringsAsFactors = FALSE
)


write.csv(
  condition_overlap_summary,
  output_path(
    "06_GF_SPF_Protein_Detection_Overlap_Summary.csv"
  ),
  row.names = FALSE
)


write.csv(
  data.frame(
    Protein.Group = detected_common
  ),
  output_path(
    "07_Proteins_Detected_in_Both_GF_and_SPF.csv"
  ),
  row.names = FALSE
)


p_overlap <- ggplot(
  condition_overlap_summary,
  aes(
    x = Category,
    y = N_Proteins
  )
) +
  geom_col() +
  geom_text(
    aes(
      label = N_Proteins
    ),
    vjust = -0.25
  ) +
  theme_bw() +
  labs(
    title = "GF vs SPF Protein Detection Overlap",
    subtitle = "Detected at least once within the corresponding condition",
    x = NULL,
    y = "Protein groups"
  )


# ==============================================================================
# 12. Core protein completeness filter
# ==============================================================================

if (
  !CORE_COMPLETENESS_MODE %in%
    c(
      "EACH_SET",
      "EACH_CONDITION"
    )
) {

  stop(
    "CORE_COMPLETENESS_MODE must be EACH_SET or EACH_CONDITION."
  )
}


if (
  CORE_COMPLETENESS_MODE ==
    "EACH_SET"
) {

  completeness_table <- data.frame(
    Protein.Group = rownames(
      log2_matrix
    ),
    stringsAsFactors = FALSE
  )


  for (s in all_set_names) {

    completeness_table[[s]] <-
      rowMeans(
        is.finite(
          log2_matrix[
            ,
            set_columns[[s]],
            drop = FALSE
          ]
        )
      )
  }


  completeness_table$Min_Completeness <-
    apply(
      completeness_table[
        ,
        all_set_names,
        drop = FALSE
      ],
      1,
      min,
      na.rm = TRUE
    )


  completeness_table$Keep_Core <-
    apply(
      completeness_table[
        ,
        all_set_names,
        drop = FALSE
      ],
      1,
      function(x) {

        all(
          is.finite(x) &
            x >=
            COMPLETENESS_CUTOFF
        )
      }
    )

} else {

  completeness_table <- data.frame(
    Protein.Group = rownames(
      log2_matrix
    ),
    GF = rowMeans(
      is.finite(
        log2_matrix[
          ,
          condition_columns$GF,
          drop = FALSE
        ]
      )
    ),
    SPF = rowMeans(
      is.finite(
        log2_matrix[
          ,
          condition_columns$SPF,
          drop = FALSE
        ]
      )
    ),
    stringsAsFactors = FALSE
  )


  completeness_table$Min_Completeness <-
    pmin(
      completeness_table$GF,
      completeness_table$SPF
    )


  completeness_table$Keep_Core <-
    completeness_table$GF >=
    COMPLETENESS_CUTOFF &
    completeness_table$SPF >=
    COMPLETENESS_CUTOFF
}


write.csv(
  completeness_table,
  output_path(
    "08_GF_SPF_Protein_Completeness_Table.csv"
  ),
  row.names = FALSE
)


core_proteins <- completeness_table$Protein.Group[
  completeness_table$Keep_Core
]


core_matrix <- log2_matrix[
  core_proteins,
  ,
  drop = FALSE
]


if (nrow(core_matrix) < 10) {

  stop(
    paste0(
      "Fewer than 10 proteins passed the core completeness filter. ",
      "Current mode: ",
      CORE_COMPLETENESS_MODE
    )
  )
}


cat("\n=====================================================\n")
cat("CORE PROTEIN FILTER\n")
cat("=====================================================\n")
cat(
  "Mode:",
  CORE_COMPLETENESS_MODE,
  "\n"
)
cat(
  "Cutoff:",
  COMPLETENESS_CUTOFF * 100,
  "%\n"
)
cat(
  "Core proteins:",
  nrow(
    core_matrix
  ),
  "\n"
)
cat("=====================================================\n")


# ==============================================================================
# 13. ROI-level QC metrics
# ==============================================================================

common_matrix <- log2_matrix[
  detected_common,
  ,
  drop = FALSE
]


roi_metrics <- data.frame(
  Sample = colnames(
    log2_matrix
  ),
  Condition = sample_info[
    colnames(
      log2_matrix
    ),
    "Condition",
    drop = TRUE
  ],
  Set = sample_info[
    colnames(
      log2_matrix
    ),
    "Set",
    drop = TRUE
  ],
  Protein_ID_Count = colSums(
    is.finite(
      log2_matrix
    )
  ),
  Missing_Pct_Common_GF_SPF = colMeans(
    !is.finite(
      common_matrix
    )
  ) * 100,
  Median_Log2_Abundance = apply(
    log2_matrix,
    2,
    function(x) {

      finite_x <- x[
        is.finite(x)
      ]

      if (length(finite_x) == 0) {
        return(NA_real_)
      }

      median(
        finite_x,
        na.rm = TRUE
      )
    }
  ),
  Total_Log2_Abundance = log2(
    colSums(
      linear_matrix,
      na.rm = TRUE
    )
  ),
  stringsAsFactors = FALSE
)


roi_metrics$Condition <- factor(
  roi_metrics$Condition,
  levels = c(
    "GF",
    "SPF"
  )
)


roi_metrics$Set <- factor(
  roi_metrics$Set,
  levels = all_set_names
)


write.csv(
  roi_metrics,
  output_path(
    "09_ROI_Level_GF_SPF_QC_Metrics.csv"
  ),
  row.names = FALSE
)


metric1 <- make_metric_plot(
  roi_metrics,
  "Protein_ID_Count",
  "Detected protein groups",
  "Protein IDs per ROI"
)

metric2 <- make_metric_plot(
  roi_metrics,
  "Missing_Pct_Common_GF_SPF",
  "Missing (%)",
  "ROI Missingness"
)

metric3 <- make_metric_plot(
  roi_metrics,
  "Median_Log2_Abundance",
  "Median log2 abundance",
  "Median Protein Abundance"
)

metric4 <- make_metric_plot(
  roi_metrics,
  "Total_Log2_Abundance",
  "log2 total abundance",
  "Total Protein Abundance"
)


roi_metric_test_summary <- data.frame(
  Metric = c(
    "Protein_ID_Count",
    "Missing_Pct_Common_GF_SPF",
    "Median_Log2_Abundance",
    "Total_Log2_Abundance"
  ),
  ROI_Level_Kruskal_P = c(
    metric1$p_value,
    metric2$p_value,
    metric3$p_value,
    metric4$p_value
  ),
  stringsAsFactors = FALSE
)


roi_metric_test_summary$FDR_BH <- p.adjust(
  roi_metric_test_summary$ROI_Level_Kruskal_P,
  method = "BH"
)


write.csv(
  roi_metric_test_summary,
  output_path(
    "10_ROI_Level_GF_SPF_QC_Metric_Tests.csv"
  ),
  row.names = FALSE
)


p_roi_qc <- (
  metric1$plot |
    metric2$plot
) / (
  metric3$plot |
    metric4$plot
)


# ==============================================================================
# 14. Set-level mean protein abundance and 8 x 8 correlation
# ==============================================================================

set_mean_matrix <- sapply(
  all_set_names,
  function(s) {

    rowMeans(
      core_matrix[
        ,
        set_columns[[s]],
        drop = FALSE
      ],
      na.rm = TRUE
    )
  }
)


set_mean_matrix[
  !is.finite(
    set_mean_matrix
  )
] <- NA_real_


rownames(
  set_mean_matrix
) <- rownames(
  core_matrix
)

colnames(
  set_mean_matrix
) <- all_set_names


write.csv(
  cbind(
    Protein.Group = rownames(
      set_mean_matrix
    ),
    as.data.frame(
      set_mean_matrix,
      check.names = FALSE
    )
  ),
  output_path(
    "11_Core_Protein_Set_Mean_Log2_Abundance.csv"
  ),
  row.names = FALSE
)


set_mean_cor <- cor(
  set_mean_matrix,
  use = "pairwise.complete.obs",
  method = "pearson"
)


write.csv(
  set_mean_cor,
  output_path(
    "12_Selected_Set_Mean_Correlation_Matrix.csv"
  ),
  row.names = TRUE
)


set_annotation <- data.frame(
  Condition = c(
    rep(
      "GF",
      length(
        GF_SET_PATTERNS
      )
    ),
    rep(
      "SPF",
      length(
        SPF_SET_PATTERNS
      )
    )
  )
)

rownames(
  set_annotation
) <- all_set_names


ph_set_cor <- pheatmap::pheatmap(
  set_mean_cor,
  annotation_col = set_annotation,
  annotation_row = set_annotation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.3f",
  border_color = NA,
  main = paste0(
    paste0(TOTAL_SET_COUNT, "-Set Mean Protein-Abundance Correlation\n"),
    "Core mode: ",
    CORE_COMPLETENESS_MODE
  ),
  silent = TRUE
)


# ==============================================================================
# 15. GF mean vs SPF mean protein correlation
# ==============================================================================

gf_mean <- rowMeans(
  set_mean_matrix[
    ,
    names(
      GF_SET_PATTERNS
    ),
    drop = FALSE
  ],
  na.rm = TRUE
)


spf_mean <- rowMeans(
  set_mean_matrix[
    ,
    names(
      SPF_SET_PATTERNS
    ),
    drop = FALSE
  ],
  na.rm = TRUE
)


gf_spf_mean_df <- data.frame(
  Protein.Group = rownames(
    set_mean_matrix
  ),
  GF_Mean_Log2 = gf_mean,
  SPF_Mean_Log2 = spf_mean,
  stringsAsFactors = FALSE
)


keep_mean_cor <-
  is.finite(
    gf_spf_mean_df$GF_Mean_Log2
  ) &
  is.finite(
    gf_spf_mean_df$SPF_Mean_Log2
  )


mean_cor_test <- cor.test(
  gf_spf_mean_df$GF_Mean_Log2[
    keep_mean_cor
  ],
  gf_spf_mean_df$SPF_Mean_Log2[
    keep_mean_cor
  ],
  method = "pearson"
)


write.csv(
  gf_spf_mean_df,
  output_path(
    "13_GF_vs_SPF_Mean_Protein_Abundance.csv"
  ),
  row.names = FALSE
)


p_mean_cor <- ggplot(
  gf_spf_mean_df[
    keep_mean_cor,
    ,
    drop = FALSE
  ],
  aes(
    x = GF_Mean_Log2,
    y = SPF_Mean_Log2
  )
) +
  geom_point(
    alpha = 0.55,
    size = 1.5
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  theme_bw() +
  labs(
    title = "GF vs SPF Mean Protein-Abundance Correlation",
    subtitle = paste0(
      "Pearson r = ",
      round(
        unname(
          mean_cor_test$estimate
        ),
        3
      ),
      "; p = ",
      format_p(
        mean_cor_test$p.value
      )
    ),
    x = "GF mean log2 abundance",
    y = "SPF mean log2 abundance"
  )


# ==============================================================================
# 16. PCA of all retained GF + SPF ROIs
#
# Median imputation ONLY for PCA / clustering.
# ==============================================================================

core_imputed <- median_impute_by_protein(
  core_matrix
)


protein_sd <- apply(
  core_imputed,
  1,
  sd,
  na.rm = TRUE
)


valid_pca_rows <-
  is.finite(
    protein_sd
  ) &
  protein_sd > 0 &
  rowSums(
    !is.finite(
      core_imputed
    )
  ) == 0


pca_matrix <- core_imputed[
  valid_pca_rows,
  ,
  drop = FALSE
]


if (nrow(pca_matrix) < 3) {

  stop(
    "Too few proteins remain for PCA."
  )
}


pca_res <- prcomp(
  t(
    pca_matrix
  ),
  center = TRUE,
  scale. = PCA_SCALE
)


variance_explained <-
  pca_res$sdev^2 /
  sum(
    pca_res$sdev^2
  ) *
  100


pca_scores <- as.data.frame(
  pca_res$x
) %>%
  rownames_to_column(
    "Sample"
  ) %>%
  left_join(
    sample_info %>%
      select(
        Sample,
        Condition,
        Set
      ),
    by = "Sample"
  )


write.csv(
  pca_scores,
  output_path(
    "14_GF_SPF_PCA_Scores.csv"
  ),
  row.names = FALSE
)


p_pca <- ggplot(
  pca_scores,
  aes(
    x = PC1,
    y = PC2,
    color = Condition
  )
) +
  geom_point(
    size = 2.5,
    alpha = 0.80
  ) +
  theme_bw() +
  labs(
    title = "PCA of GF and SPF Final-MBR ROIs",
    subtitle = paste0(
      "Core mode: ",
      CORE_COMPLETENESS_MODE,
      "; median imputation only for PCA"
    ),
    x = paste0(
      "PC1 (",
      round(
        variance_explained[1],
        1
      ),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(
        variance_explained[2],
        1
      ),
      "%)"
    )
  )


condition_counts <- table(
  pca_scores$Condition
)


eligible_conditions <- names(
  condition_counts[
    condition_counts >= 3
  ]
)


if (length(eligible_conditions) > 0) {

  p_pca <- p_pca +
    stat_ellipse(
      data = pca_scores %>%
        filter(
          Condition %in%
            eligible_conditions
        ),
      aes(
        group = Condition,
        color = Condition
      ),
      level = 0.95,
      linewidth = 0.7,
      show.legend = FALSE
    )
}


# Set-level PCA centroids
pca_set_centroids <- pca_scores %>%
  group_by(
    Condition,
    Set
  ) %>%
  summarise(
    PC1 = mean(
      PC1,
      na.rm = TRUE
    ),
    PC2 = mean(
      PC2,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


write.csv(
  pca_set_centroids,
  output_path(
    "15_PCA_Set_Centroids.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 17. ROI-to-ROI correlation heatmap
# ==============================================================================

roi_cor <- suppressWarnings(
  cor(
    core_matrix,
    use = "pairwise.complete.obs",
    method = "pearson"
  )
)


roi_annotation <- sample_info[
  colnames(
    roi_cor
  ),
  c(
    "Condition",
    "Set"
  ),
  drop = FALSE
]


ph_roi_cor <- pheatmap::pheatmap(
  roi_cor,
  annotation_col = roi_annotation,
  annotation_row = roi_annotation,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_row_names = ROI_COR_SHOW_NAMES,
  show_col_names = ROI_COR_SHOW_NAMES,
  border_color = NA,
  main = paste0(
    "ROI-to-ROI Pearson Correlation: GF vs SPF\n",
    "Non-imputed core protein matrix"
  ),
  silent = TRUE
)


# ==============================================================================
# 18. Top-variable-protein heatmap
# ==============================================================================

protein_variance <- apply(
  pca_matrix,
  1,
  var,
  na.rm = TRUE
)


protein_variance <- protein_variance[
  is.finite(
    protein_variance
  )
]


n_top_variable <- min(
  TOP_N_VARIABLE_PROTEINS,
  length(
    protein_variance
  )
)


top_variable_proteins <- names(
  sort(
    protein_variance,
    decreasing = TRUE
  )
)[
  seq_len(
    n_top_variable
  )
]


variable_heatmap_matrix <- pca_matrix[
  top_variable_proteins,
  ,
  drop = FALSE
]


variable_heatmap_z <- t(
  scale(
    t(
      variable_heatmap_matrix
    )
  )
)


valid_variable_rows <-
  rowSums(
    !is.finite(
      variable_heatmap_z
    )
  ) == 0


variable_heatmap_z <- variable_heatmap_z[
  valid_variable_rows,
  ,
  drop = FALSE
]


ph_variable <- pheatmap::pheatmap(
  variable_heatmap_z,
  annotation_col = roi_annotation,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_row_names =
    VARIABLE_HEATMAP_SHOW_ROW_NAMES,
  show_col_names =
    VARIABLE_HEATMAP_SHOW_COL_NAMES,
  border_color = NA,
  main = paste0(
    "Top ",
    nrow(
      variable_heatmap_z
    ),
    " Variable Proteins: GF vs SPF\n",
    "Protein-wise Z-score"
  ),
  silent = TRUE
)


# ==============================================================================
# 19. SET-LEVEL protein matrix
#
# Each protein receives ONE summary value per original set.
#
# This produces:
#   one observation per selected GF/SPF set and protein.
# ==============================================================================

set_median_matrix <- sapply(
  all_set_names,
  function(s) {

    apply(
      core_matrix[
        ,
        set_columns[[s]],
        drop = FALSE
      ],
      1,
      function(x) {

        x <- x[
          is.finite(x)
        ]

        if (length(x) == 0) {
          return(NA_real_)
        }

        median(
          x,
          na.rm = TRUE
        )
      }
    )
  }
)


rownames(
  set_median_matrix
) <- rownames(
  core_matrix
)

colnames(
  set_median_matrix
) <- all_set_names


write.csv(
  cbind(
    Protein.Group = rownames(
      set_median_matrix
    ),
    as.data.frame(
      set_median_matrix,
      check.names = FALSE
    )
  ),
  output_path(
    "16_Set_Level_Median_Log2_Abundance_Matrix.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 20. Selected GF-set vs SPF-set protein comparison (MAIN INFERENCE)
#
# Primary:
#   Welch t-test on SET-level summaries.
#
# Sensitivity:
#   Wilcoxon rank-sum test on SET-level summaries.
#
# Effect:
#   log2FC = mean(SPF set values) - mean(GF set values)
# ==============================================================================

gf_sets <- names(
  GF_SET_PATTERNS
)

spf_sets <- names(
  SPF_SET_PATTERNS
)


diff_result_list <- vector(
  "list",
  nrow(
    set_median_matrix
  )
)


for (i in seq_len(
  nrow(
    set_median_matrix
  )
)) {

  protein_id <- rownames(
    set_median_matrix
  )[i]


  gf_values <- as.numeric(
    set_median_matrix[
      i,
      gf_sets,
      drop = TRUE
    ]
  )


  spf_values <- as.numeric(
    set_median_matrix[
      i,
      spf_sets,
      drop = TRUE
    ]
  )


  gf_values <- gf_values[
    is.finite(
      gf_values
    )
  ]

  spf_values <- spf_values[
    is.finite(
      spf_values
    )
  ]


  n_gf <- length(
    gf_values
  )

  n_spf <- length(
    spf_values
  )


  mean_gf <- ifelse(
    n_gf > 0,
    mean(
      gf_values,
      na.rm = TRUE
    ),
    NA_real_
  )


  mean_spf <- ifelse(
    n_spf > 0,
    mean(
      spf_values,
      na.rm = TRUE
    ),
    NA_real_
  )


  median_gf <- ifelse(
    n_gf > 0,
    median(
      gf_values,
      na.rm = TRUE
    ),
    NA_real_
  )


  median_spf <- ifelse(
    n_spf > 0,
    median(
      spf_values,
      na.rm = TRUE
    ),
    NA_real_
  )


  log2fc <- mean_spf -
    mean_gf


  welch_p <- NA_real_
  wilcox_p <- NA_real_


  if (
    n_gf >=
      MIN_SETS_PER_CONDITION_FOR_TEST &&
      n_spf >=
      MIN_SETS_PER_CONDITION_FOR_TEST
  ) {

    welch_p <- tryCatch(
      t.test(
        spf_values,
        gf_values,
        var.equal = FALSE
      )$p.value,
      error = function(e) NA_real_
    )


    wilcox_p <- tryCatch(
      suppressWarnings(
        wilcox.test(
          spf_values,
          gf_values,
          exact = FALSE
        )$p.value
      ),
      error = function(e) NA_real_
    )
  }


  diff_result_list[[i]] <- data.frame(
    Protein.Group = protein_id,
    N_GF_Sets = n_gf,
    N_SPF_Sets = n_spf,
    Mean_Log2_GF = mean_gf,
    Mean_Log2_SPF = mean_spf,
    Median_Log2_GF = median_gf,
    Median_Log2_SPF = median_spf,
    Log2FC_SPF_vs_GF = log2fc,
    FoldChange_SPF_vs_GF = ifelse(
      is.finite(
        log2fc
      ),
      2^log2fc,
      NA_real_
    ),
    Welch_P = welch_p,
    Wilcoxon_P = wilcox_p,
    stringsAsFactors = FALSE
  )
}


diff_results <- bind_rows(
  diff_result_list
)


diff_results$Welch_FDR_BH <- p.adjust(
  diff_results$Welch_P,
  method = "BH"
)


diff_results$Wilcoxon_FDR_BH <- p.adjust(
  diff_results$Wilcoxon_P,
  method = "BH"
)


diff_results$Strong_GF_SPF_Difference <-
  is.finite(
    diff_results$Welch_FDR_BH
  ) &
  diff_results$Welch_FDR_BH <
  DIFF_FDR_CUTOFF &
  is.finite(
    diff_results$Log2FC_SPF_vs_GF
  ) &
  abs(
    diff_results$Log2FC_SPF_vs_GF
  ) >=
  DIFF_MIN_ABS_LOG2FC


diff_results <- diff_results %>%
  arrange(
    Welch_FDR_BH,
    desc(
      abs(
        Log2FC_SPF_vs_GF
      )
    )
  )


write.csv(
  diff_results,
  output_path(
    "17_GF_vs_SPF_Set_Level_Protein_Comparison.csv"
  ),
  row.names = FALSE
)


strong_diff <- diff_results %>%
  filter(
    Strong_GF_SPF_Difference
  )


write.csv(
  strong_diff,
  output_path(
    "18_Strong_GF_vs_SPF_Differential_Proteins.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 21. Volcano-style plot
# ==============================================================================

volcano_df <- diff_results %>%
  mutate(
    MinusLog10FDR =
      -log10(
        pmax(
          Welch_FDR_BH,
          .Machine$double.xmin
        )
      )
  )


p_volcano <- ggplot(
  volcano_df,
  aes(
    x = Log2FC_SPF_vs_GF,
    y = MinusLog10FDR,
    shape = Strong_GF_SPF_Difference
  )
) +
  geom_point(
    alpha = 0.65,
    size = 1.7
  ) +
  geom_vline(
    xintercept = c(
      -DIFF_MIN_ABS_LOG2FC,
      DIFF_MIN_ABS_LOG2FC
    ),
    linetype = 2
  ) +
  geom_hline(
    yintercept =
      -log10(
        DIFF_FDR_CUTOFF
      ),
    linetype = 2
  ) +
  theme_bw() +
  labs(
    title = paste0("GF vs SPF Protein Comparison (", COMPARISON_LABEL, ")"),
    subtitle = "Set-level Welch test; positive log2FC = higher in SPF",
    x = "log2FC (SPF - GF)",
    y = "-log10(BH-FDR)",
    shape = "Strong difference"
  )


# ==============================================================================
# 22. Differential-protein set-level heatmap
# ==============================================================================

ph_diff <- NULL


diff_heatmap_candidates <-
  strong_diff$Protein.Group


# If no protein passes both FDR and fold cutoff, show the top proteins by FDR
# as an exploratory heatmap, but do NOT relabel them as significant.
if (length(
  diff_heatmap_candidates
) < 2) {

  diff_heatmap_candidates <-
    diff_results$Protein.Group[
      is.finite(
        diff_results$Welch_FDR_BH
      )
    ]

  diff_heatmap_candidates <- head(
    diff_heatmap_candidates,
    TOP_N_DIFFERENTIAL_PROTEINS
  )

  diff_heatmap_title <-
    "Top GF-vs-SPF Proteins by Welch FDR\nExploratory; not necessarily FDR-significant"

} else {

  diff_heatmap_candidates <- head(
    diff_heatmap_candidates,
    TOP_N_DIFFERENTIAL_PROTEINS
  )

  diff_heatmap_title <-
    paste0(
      "Strong GF-vs-SPF Differential Proteins\n",
      "Welch FDR < ",
      DIFF_FDR_CUTOFF,
      "; |log2FC| >= ",
      DIFF_MIN_ABS_LOG2FC
    )
}


if (
  length(
    diff_heatmap_candidates
  ) >= 2
) {

  diff_heatmap_matrix <-
    set_median_matrix[
      diff_heatmap_candidates,
      ,
      drop = FALSE
    ]


  diff_heatmap_z <- t(
    scale(
      t(
        diff_heatmap_matrix
      )
    )
  )


  valid_diff_rows <-
    rowSums(
      !is.finite(
        diff_heatmap_z
      )
    ) == 0


  diff_heatmap_z <- diff_heatmap_z[
    valid_diff_rows,
    ,
    drop = FALSE
  ]


  if (
    nrow(
      diff_heatmap_z
    ) >= 2
  ) {

    ph_diff <- pheatmap::pheatmap(
      diff_heatmap_z,
      annotation_col = set_annotation,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      show_row_names =
        DIFF_HEATMAP_SHOW_ROW_NAMES,
      border_color = NA,
      main = diff_heatmap_title,
      silent = TRUE
    )
  }
}


# ==============================================================================
# 23. Save figures
# ==============================================================================

save_ggplot_png(
  p_detected,
  output_path(
    "Figure_01_Detected_Protein_Counts_Selected_Sets.png"
  ),
  width = 12,
  height = 7
)


save_ggplot_png(
  p_complete70,
  output_path(
    "Figure_02_70pct_Completeness_Protein_Counts_Selected_Sets.png"
  ),
  width = 12,
  height = 7
)


save_ggplot_png(
  p_overlap,
  output_path(
    "Figure_03_GF_SPF_Protein_Detection_Overlap.png"
  )
)


save_ggplot_png(
  p_roi_qc,
  output_path(
    "Figure_04_GF_SPF_ROI_QC_Metrics.png"
  ),
  width = 13,
  height = 10
)


save_pheatmap(
  ph_set_cor,
  output_path(
    "Figure_05_Selected_Set_Mean_Correlation"
  ),
  width = 9,
  height = 8
)


save_ggplot_png(
  p_mean_cor,
  output_path(
    "Figure_06_GF_vs_SPF_Mean_Protein_Correlation.png"
  )
)


save_ggplot_png(
  p_pca,
  output_path(
    "Figure_07_GF_SPF_PCA.png"
  ),
  width = 10,
  height = 8
)


save_pheatmap(
  ph_roi_cor,
  output_path(
    "Figure_08_GF_SPF_ROI_Pearson_Correlation"
  )
)


save_pheatmap(
  ph_variable,
  output_path(
    "Figure_09_GF_SPF_Top_Variable_Protein_Heatmap"
  )
)


save_ggplot_png(
  p_volcano,
  output_path(
    "Figure_10_GF_vs_SPF_Set_Level_Volcano.png"
  )
)


if (!is.null(
  ph_diff
)) {

  save_pheatmap(
    ph_diff,
    output_path(
      "Figure_11_GF_vs_SPF_Differential_Protein_Heatmap"
    ),
    width = 9,
    height = 10
  )
}


# ==============================================================================
# 24. Analysis source manifest
# ==============================================================================

manifest <- data.frame(
  Analysis = c(
    "Detected protein counts",
    "70% completeness counts",
    "GF-SPF detection overlap",
    "ROI QC metrics",
    "Selected-set correlation",
    "GF-SPF mean correlation",
    "PCA",
    "ROI correlation",
    "Top-variable heatmap",
    "GF-vs-SPF protein comparison",
    "Differential-protein heatmap"
  ),
  Statistical_Unit = c(
    "Set descriptive",
    "Set descriptive",
    "Condition descriptive",
    "ROI descriptive",
    "Set mean",
    "Condition mean",
    "ROI descriptive",
    "ROI descriptive",
    "ROI descriptive",
    COMPARISON_LABEL,
    "Set median"
  ),
  Imputation = c(
    "None",
    "None",
    "None",
    "None",
    "None",
    "None",
    "Protein-median imputation",
    "None",
    "Protein-median imputation",
    "None",
    "None"
  ),
  Notes = c(
    "Protein detected at least once in each set",
    "Calculated independently within each set",
    "Detected at least once within GF/SPF",
    "Do not interpret ROI p-values as biological replication",
    "Pearson correlation of set-mean protein abundance",
    "Pearson correlation of GF mean vs SPF mean abundance",
    "Combined GF + SPF core protein matrix",
    "Pairwise-complete Pearson correlation",
    "Protein-wise Z-score",
    "Welch + Wilcoxon on set-level summaries; log2FC = SPF - GF",
    paste0("Protein-wise Z-score across ", TOTAL_SET_COUNT, " set summaries")
  ),
  stringsAsFactors = FALSE
)


write.csv(
  manifest,
  output_path(
    "19_Analysis_Source_Manifest.csv"
  ),
  row.names = FALSE
)


# ==============================================================================
# 25. Save all major figures into one PDF
# ==============================================================================

all_plots_pdf <- output_path(
  "GF_vs_SPF_Selected_Sets_All_Comparison_Plots.pdf"
)


grDevices::pdf(
  all_plots_pdf,
  width = ALL_PLOTS_PDF_WIDTH,
  height = ALL_PLOTS_PDF_HEIGHT,
  onefile = TRUE
)


print(
  p_detected
)

print(
  p_complete70
)

print(
  p_overlap
)

print(
  p_roi_qc
)


grid::grid.newpage()
grid::grid.draw(
  ph_set_cor$gtable
)


print(
  p_mean_cor
)

print(
  p_pca
)


grid::grid.newpage()
grid::grid.draw(
  ph_roi_cor$gtable
)


grid::grid.newpage()
grid::grid.draw(
  ph_variable$gtable
)


print(
  p_volcano
)


if (!is.null(
  ph_diff
)) {

  grid::grid.newpage()
  grid::grid.draw(
    ph_diff$gtable
  )
}


grDevices::dev.off()


# ==============================================================================
# 26. Final console summary
# ==============================================================================

cat("\n=====================================================\n")
cat("FINAL GF vs SPF SELECTED-SET JOINT-MBR SUMMARY\n")
cat("=====================================================\n")

cat(
  "Single joint-MBR report:",
  REPORT_FILE,
  "\n"
)

cat(
  "Abundance column:",
  gf$abundance_col,
  "\n"
)

cat(
  "Core completeness mode:",
  CORE_COMPLETENESS_MODE,
  "\n"
)

cat(
  "Completeness cutoff:",
  COMPLETENESS_CUTOFF * 100,
  "%\n"
)


cat("\nGF sets / ROI counts:\n")

for (s in names(
  GF_SET_PATTERNS
)) {

  cat(
    s,
    ":",
    length(
      set_columns[[s]]
    ),
    "\n"
  )
}


cat("\nSPF sets / ROI counts:\n")

for (s in names(
  SPF_SET_PATTERNS
)) {

  cat(
    s,
    ":",
    length(
      set_columns[[s]]
    ),
    "\n"
  )
}


cat(
  "\nProteins detected in GF:",
  length(
    gf_detected
  ),
  "\n"
)

cat(
  "Proteins detected in SPF:",
  length(
    spf_detected
  ),
  "\n"
)

cat(
  "Proteins detected in both:",
  length(
    detected_common
  ),
  "\n"
)

cat(
  "Core proteins used for comparison:",
  nrow(
    core_matrix
  ),
  "\n"
)

cat(
  "Proteins used for PCA:",
  nrow(
    pca_matrix
  ),
  "\n"
)

cat(
  "Strong GF-vs-SPF proteins:",
  nrow(
    strong_diff
  ),
  "\n"
)

cat(
  "Output directory:",
  OUTPUT_DIR,
  "\n"
)

cat(
  "All-plots PDF:",
  all_plots_pdf,
  "\n"
)


cat("\nINTERPRETATION:\n")
cat(
  paste0("- All ", TOTAL_SET_COUNT,
         " selected sets originate from the same DIA-NN joint-MBR analysis.\n")
)
cat(
  "- ROI-level p-values are descriptive QC only.\n"
)
cat(
  paste0("- Main protein statistics use ", length(GF_SET_PATTERNS),
         " GF set summaries vs ", length(SPF_SET_PATTERNS),
         " SPF set summaries.\n")
)
cat(
  "- Positive Log2FC_SPF_vs_GF indicates higher abundance in SPF.\n"
)
cat(
  "- Negative Log2FC_SPF_vs_GF indicates higher abundance in GF.\n"
)

cat("=====================================================\n")

comparison_summary <- data.frame(
  GF_Sets = length(GF_SET_PATTERNS),
  SPF_Sets = length(SPF_SET_PATTERNS),
  GF_ROIs = sum(vapply(set_columns[names(GF_SET_PATTERNS)], length, integer(1))),
  SPF_ROIs = sum(vapply(set_columns[names(SPF_SET_PATTERNS)], length, integer(1))),
  Detected_GF = length(gf_detected),
  Detected_SPF = length(spf_detected),
  Detected_Both = length(detected_common),
  Core_Proteins = nrow(core_matrix),
  PCA_Proteins = nrow(pca_matrix),
  Strong_Differential_Proteins = nrow(strong_diff),
  GF_Selected = paste(names(GF_SET_PATTERNS), collapse = ";"),
  SPF_Selected = paste(names(SPF_SET_PATTERNS), collapse = ";"),
  stringsAsFactors = FALSE
)
utils::write.csv(
  comparison_summary,
  output_path("20_Comparison_Run_Summary.csv"),
  row.names = FALSE
)

comparison_files <- list.files(
  OUTPUT_DIR, recursive = TRUE, full.names = TRUE, include.dirs = FALSE
)
invisible(list(summary = comparison_summary, files = comparison_files))
}
