# ==============================================================================
# DIA-NN Spatial Proteomics Dual QC — automatic GF/SPF dataset IDs
# Paired report-first-pass.parquet + report.parquet; Shiny entry point: run_qc().
#
# Dataset IDs are read from Run names (GF1-4, SPF1-2, GF12-10, GF1_4, ...).
# Different datasets are QC'ed independently; ROI1 from different IDs is not pooled.
# Source this file to define functions; it does not start analysis or change folders.
#
# Standalone use:
# source("QC_Code.R")
# run_qc("report-first-pass.parquet", "report.parquet", "QC_results")
#
# QC cutoffs, 70% completeness, no imputation and plot settings match the supplied
# pipeline. The original manual exclusions are retained ONLY for GF1-4.
# To use no manual exclusions: run_qc(..., manual_exclude_by_dataset = list()).
# ==============================================================================

QC_REQUIRED_PACKAGES <- c("DBI", "duckdb", "dplyr", "tidyr", "tibble",
                          "ggplot2", "pheatmap", "patchwork")
QC_REQUIRED_COLUMNS <- c("Run", "Protein.Group", "Q.Value", "PG.Q.Value", "PG.MaxLFQ")
QC_DATASET_PATTERN <- "(?<![[:alnum:]])(?:SPF|GF)[0-9]+(?:[-_][0-9]+)?(?![[:alnum:]-])"
QC_ROI_PATTERN <- "(?<![[:alnum:]])ROI0*([0-9]+)"
QC_MANUAL_EXCLUSIONS <- list(
  "GF1-4" = c("ROI5", "ROI63", "ROI79", "ROI34", "ROI38", "ROI40", "ROI45",
              "ROI76", "ROI29", "ROI30", "ROI31", "ROI39", "ROI41", "ROI74")
)

qc_check_packages <- function() {
  missing <- QC_REQUIRED_PACKAGES[!vapply(
    QC_REQUIRED_PACKAGES, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop(
    "Install QC packages first: install.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "), "))", call. = FALSE)
}

parse_qc_runs <- function(run_names, report_label, roi_min = 1L, roi_max = 96L) {
  run_names <- unique(as.character(run_names[!is.na(run_names)]))
  if (!length(run_names)) stop("[", report_label, "] No Run names.", call. = FALSE)
  # Use the file basename, not folder names, to avoid attributing a parent ID
  # to a differently named Run. Both Windows and Unix path separators are accepted.
  base <- basename(gsub("\\", "/", run_names, fixed = TRUE))
  upper <- toupper(base)
  tokens <- regmatches(upper, gregexpr(QC_DATASET_PATTERN, upper, perl = TRUE))
  ids <- vapply(tokens, function(x) {
    x <- unique(gsub("_", "-", x, fixed = TRUE))
    if (length(x) == 1L) x else NA_character_
  }, character(1))
  ambiguous <- vapply(tokens, function(x) {
    length(unique(gsub("_", "-", x, fixed = TRUE))) > 1L
  }, logical(1))
  if (any(ambiguous)) stop("[", report_label, "] Ambiguous dataset IDs in Run: ",
                            paste(head(run_names[ambiguous], 5L), collapse = "; "), call. = FALSE)
  roi_matches <- regmatches(upper, regexec(QC_ROI_PATTERN, upper, perl = TRUE))
  roi <- vapply(roi_matches, function(x) {
    if (length(x) >= 2L) suppressWarnings(as.integer(x[[2]])) else NA_integer_
  }, integer(1))
  in_range <- !is.na(roi) & roi >= roi_min & roi <= roi_max
  if (any(in_range & is.na(ids))) stop(
    "[", report_label, "] Cannot determine GF/SPF dataset ID from ROI Run: ",
    paste(head(run_names[in_range & is.na(ids)], 5L), collapse = "; "),
    ". Include an ID such as GF1-4 or SPF1-2 in each Run name.", call. = FALSE)
  reason <- ifelse(is.na(roi), "No_ROI_token",
                    ifelse(!in_range, "Outside_ROI_range", "Included"))
  out <- data.frame(Report = report_label, Original_Run = run_names,
                    Dataset_ID = ids, ROI_Number = roi, Parsing_Status = reason,
                    stringsAsFactors = FALSE)
  included <- out[out$Parsing_Status == "Included", , drop = FALSE]
  if (!nrow(included)) stop("[", report_label, "] No valid GF/SPF ROI Runs.", call. = FALSE)
  key <- paste(included$Dataset_ID, included$ROI_Number, sep = "/")
  if (anyDuplicated(key)) stop(
    "[", report_label, "] Multiple original Runs map to the same dataset/ROI: ",
    paste(unique(key[duplicated(key) | duplicated(key, fromLast = TRUE)]), collapse = ", "),
    ". Separate technical reruns before QC.", call. = FALSE)
  out
}

qc_register_report <- function(con, path, view) {
  sql_view <- as.character(DBI::dbQuoteIdentifier(con, view))
  sql_path <- as.character(DBI::dbQuoteString(con, normalizePath(
    path, winslash = "/", mustWork = TRUE)))
  DBI::dbExecute(con, paste("CREATE OR REPLACE VIEW", sql_view,
                           "AS SELECT * FROM read_parquet(", sql_path, ")"))
  columns <- DBI::dbGetQuery(con, paste("DESCRIBE", sql_view))$column_name
  missing <- setdiff(QC_REQUIRED_COLUMNS, columns)
  if (length(missing)) stop("[", view, "] Missing columns: ",
                            paste(missing, collapse = ", "), call. = FALSE)
  DBI::dbGetQuery(con, paste(
    'SELECT DISTINCT CAST("Run" AS VARCHAR) AS run_name FROM', sql_view))$run_name
}

qc_report_aggregate_sql <- function(view) {
  # view is an internal fixed name, not user-supplied SQL.
  if (!view %in% c("qc_first_input", "qc_final_input")) stop("Invalid QC source view.")
  paste0(
    "WITH typed AS (SELECT m.Original_Run, m.ROI_Number, ",
    'CAST(r."Protein.Group" AS VARCHAR) AS protein_group, ',
    'TRY_CAST(r."Q.Value" AS DOUBLE) AS q_value, ',
    'TRY_CAST(r."PG.Q.Value" AS DOUBLE) AS pg_q_value, ',
    'TRY_CAST(r."PG.MaxLFQ" AS DOUBLE) AS intensity ',
    "FROM ", view, " r JOIN qc_roi_map m ",
    'ON CAST(r."Run" AS VARCHAR) = m.Original_Run), ',
    "passing AS (SELECT * FROM typed WHERE q_value <= 0.01 AND pg_q_value <= 0.01 ",
    "AND protein_group IS NOT NULL ",
    "AND NOT regexp_matches(protein_group, 'cRAP|CON__|CONTAMINANT', 'i')) ",
    "SELECT MIN(Original_Run) AS Original_Run, ROI_Number, ",
    '\'ROI\' || CAST(ROI_Number AS VARCHAR) AS "Run", ',
    'protein_group AS "Protein.Group", MIN(q_value) AS "Q.Value", ',
    'MIN(pg_q_value) AS "PG.Q.Value", ',
    'MAX(CASE WHEN isnan(intensity) THEN NULL ELSE intensity END) AS "PG.MaxLFQ" ',
    "FROM passing GROUP BY ROI_Number, protein_group ORDER BY ROI_Number, protein_group"
  )
}

qc_read_report <- function(con, view, run_map) {
  # The same FDR and contaminant rules as make_detection_matrix / make_protein_matrix.
  # Aggregate max(PG.MaxLFQ) per ROI/Protein.Group on disk, before collecting in R.
  # The complete Run map separately retains ROIs with zero FDR-passing proteins.
  DBI::dbWriteTable(con, "qc_roi_map",
                    run_map[, c("Original_Run", "ROI_Number"), drop = FALSE],
                    temporary = TRUE, overwrite = TRUE, row.names = FALSE)
  DBI::dbGetQuery(con, qc_report_aggregate_sql(view))
}

run_qc <- function(first_pass_path, report_path, output_dir,
                   dataset_ids = NULL, manual_exclude_by_dataset = QC_MANUAL_EXCLUSIONS,
                   db_memory_limit = "1GB", png_dpi = 600L) {
  qc_check_packages()
  # Attach the same plotting/data-manipulation packages used by the supplied code.
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr); library(tibble)
    library(ggplot2); library(pheatmap); library(patchwork)
  })
  first_pass_path <- normalizePath(first_pass_path, winslash = "/", mustWork = TRUE)
  report_path <- normalizePath(report_path, winslash = "/", mustWork = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  if (!is.list(manual_exclude_by_dataset) ||
      (length(manual_exclude_by_dataset) && is.null(names(manual_exclude_by_dataset)))) {
    stop("manual_exclude_by_dataset must be a named list keyed by dataset ID.", call. = FALSE)
  }
  if (!is.numeric(png_dpi) || length(png_dpi) != 1L ||
      !is.finite(png_dpi) || png_dpi < 72 || png_dpi > 1200) {
    stop("png_dpi must be between 72 and 1200.", call. = FALSE)
  }
  temp_work <- tempfile("qc_sql_", tmpdir = output_dir)
  dir.create(temp_work)
  spill <- file.path(temp_work, "spill"); dir.create(spill)
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = file.path(temp_work, "qc.duckdb"),
    bigint = "numeric", config = list(
      memory_limit = db_memory_limit, threads = "2", preserve_insertion_order = "false",
      temp_directory = normalizePath(spill, winslash = "/", mustWork = TRUE)))
  on.exit({
    try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
    unlink(temp_work, recursive = TRUE)
  }, add = TRUE)
  first_runs <- qc_register_report(con, first_pass_path, "qc_first_input")
  final_runs <- qc_register_report(con, report_path, "qc_final_input")
  first_map <- parse_qc_runs(first_runs, "First-pass")
  final_map <- parse_qc_runs(final_runs, "Final-MBR")
  run_map <- rbind(first_map, final_map)
  utils::write.csv(run_map, file.path(output_dir, "QC_Run_Name_Map.csv"), row.names = FALSE)
  first_ids <- sort(unique(first_map$Dataset_ID[first_map$Parsing_Status == "Included"]))
  final_ids <- sort(unique(final_map$Dataset_ID[final_map$Parsing_Status == "Included"]))
  if (!setequal(first_ids, final_ids)) stop(
    "Dataset IDs differ between the paired reports. First-pass only: ",
    paste(setdiff(first_ids, final_ids), collapse = ", "),
    "; Final-MBR only: ", paste(setdiff(final_ids, first_ids), collapse = ", "), call. = FALSE)
  ids <- first_ids
  if (!is.null(dataset_ids)) {
    dataset_ids <- unique(gsub("_", "-", toupper(as.character(dataset_ids)), fixed = TRUE))
    if (!all(dataset_ids %in% ids)) stop("Requested dataset ID was not present in both reports.")
    ids <- ids[ids %in% dataset_ids]
  }
  if (!length(ids)) stop("No datasets selected for QC.", call. = FALSE)
  cat("Auto-detected datasets:", paste(ids, collapse = ", "), "\n")
  cat("Each dataset is analysed independently; ROI numbers are not pooled across IDs.\n")
  cat("Ignored non-ROI/out-of-range Runs:",
      sum(run_map$Parsing_Status != "Included"), "\n")
  summaries <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    id <- ids[[i]]
    first_sub <- first_map[first_map$Dataset_ID %in% id &
                             first_map$Parsing_Status == "Included", , drop = FALSE]
    final_sub <- final_map[final_map$Dataset_ID %in% id &
                             final_map$Parsing_Status == "Included", , drop = FALSE]
    manual <- manual_exclude_by_dataset[[id]]
    if (is.null(manual)) manual <- character()
    manual <- unique(as.character(manual))
    if (anyNA(manual) || any(!grepl("^ROI[0-9]+$", manual))) {
      stop("Invalid manual ROI names for dataset ", id, call. = FALSE)
    }
    cat("\nDataset", i, "of", length(ids), ":", id, "\n")
    cat("Manual exclusions:", if (length(manual)) paste(manual, collapse = ", ") else "None", "\n")
    summaries[[i]] <- run_qc_one_dataset(
      con, id, first_pass_path, report_path, first_sub, final_sub,
      file.path(output_dir, id), manual, as.integer(png_dpi))
    invisible(gc(verbose = FALSE))
  }
  summary <- do.call(rbind, summaries)
  utils::write.csv(summary, file.path(output_dir, "QC_Dataset_Summary.csv"), row.names = FALSE)
  files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  files <- files[!startsWith(normalizePath(files, winslash = "/", mustWork = TRUE),
                            paste0(normalizePath(temp_work, winslash = "/", mustWork = TRUE), "/"))]
  list(summary = summary, files = files)
}



run_qc_one_dataset <- function(QC_CON, dataset_id, first_pass_path, report_path,
                                FIRST_RUN_MAP, FINAL_RUN_MAP, dataset_output_dir,
                                manual_exclude, png_dpi = 600L) {
  dir.create(dataset_output_dir, recursive = TRUE, showWarnings = FALSE)
  dataset_output_dir <- normalizePath(dataset_output_dir, winslash = "/", mustWork = TRUE)
  old_wd <- getwd()
  old_devices <- as.integer(grDevices::dev.list())
  setwd(dataset_output_dir)
  # Suppress intermediate interactive displays; the supplied export code writes
  # the actual PNG/PDF results. Close only devices opened by this dataset.
  grDevices::pdf(file = NULL)
  on.exit({
    opened <- setdiff(as.integer(grDevices::dev.list()), old_devices)
    for (device in rev(opened)) try(grDevices::dev.off(device), silent = TRUE)
    setwd(old_wd)
  }, add = TRUE)
  FIRST_PASS_FILE <- first_pass_path
  FINAL_MBR_FILE <- report_path
  OUTPUT_PREFIX <- dataset_id
  MANUAL_EXCLUDE <- manual_exclude

  # ROI range
  ROI_MIN <- 1
  ROI_MAX <- 96

  # Primary QC thresholds
  FAIL_NMAD    <- 3.0
  WARNING_NMAD <- 2.5

  # Secondary QC warning threshold
  SECONDARY_NMAD <- 3.0

  # Protein completeness for downstream zonation
  PROTEIN_COMPLETENESS <- 0.70

  # Optional manual exclusions after reviewing REVIEW samples.
  # Example:
  # MANUAL_EXCLUDE <- c("ROI12", "ROI44") or character(0)



  # Correlation heatmap output settings ------------------------------------------
  # 두 heatmap의 ROI 순서는 아래 reference dataset의 hierarchical clustering
  # 순서를 공통으로 사용합니다.
  #
  # "Final-MBR" 권장:
  #   Final-MBR의 completeness가 더 높기 때문에 일반적으로 안정적인
  #   ROI ordering reference로 사용하기 좋습니다.
  HEATMAP_ORDER_REFERENCE <- "Final-MBR"

  # Figure size / resolution
  HEATMAP_SINGLE_WIDTH  <- 10
  HEATMAP_SINGLE_HEIGHT <- 9

  HEATMAP_COMBINED_WIDTH  <- 20
  HEATMAP_COMBINED_HEIGHT <- 9

  HEATMAP_PNG_DPI <- png_dpi

  # General QC plot output settings ----------------------------------------------
  # PNG는 600 dpi로 저장합니다.
  # PDF는 vector format이므로 dpi 설정이 필요하지 않습니다.
  PLOT_PNG_DPI <- png_dpi

  PLOT_SINGLE_WIDTH  <- 10
  PLOT_SINGLE_HEIGHT <- 7

  # ROI가 많은 bar plot은 세로 길이를 조금 더 크게 저장
  PLOT_TALL_WIDTH  <- 10
  PLOT_TALL_HEIGHT <- 12

  # 모든 plot을 하나의 multipage PDF에 저장할 때 사용할 크기
  ALL_PLOTS_PDF_WIDTH  <- 12
  ALL_PLOTS_PDF_HEIGHT <- 9


  safe_nmad <- function(x) {

    result <- mad(
      x,
      center = median(x, na.rm = TRUE),
      constant = 1.4826,
      na.rm = TRUE
    )

    if (is.na(result) || result == 0) {
      result <- 1e-6
    }

    return(result)
  }


  check_required_columns <- function(df, file_label) {

    required_cols <- c(
      "Run",
      "Protein.Group",
      "Q.Value",
      "PG.Q.Value",
      "PG.MaxLFQ"
    )

    missing_cols <- setdiff(required_cols, colnames(df))

    if (length(missing_cols) > 0) {
      stop(
        paste0(
          "[", file_label, "] Missing required columns: ",
          paste(missing_cols, collapse = ", ")
        )
      )
    }
  }


  prepare_diann_data <- function(df) {
    df %>%
      filter(!is.na(ROI_Number), ROI_Number >= ROI_MIN, ROI_Number <= ROI_MAX) %>%
      mutate(Run = paste0("ROI", ROI_Number))
  }

  # ==============================================================================
  # Protein detection matrix for QC
  #
  # First-pass report may contain PG.MaxLFQ as all NA.
  # Therefore:
  #   - Identified protein count
  #   - Missing percentage
  # are calculated from FDR-passing Protein.Group PRESENCE, not PG.MaxLFQ.
  #
  # 1 = detected
  # 0 = not detected
  # ==============================================================================

  make_detection_matrix <- function(df, all_runs) {

    detection_summary <- df %>%

      filter(
        Q.Value <= 0.01,
        PG.Q.Value <= 0.01,
        !is.na(Protein.Group),
        !grepl(
          "cRAP|CON__|CONTAMINANT",
          Protein.Group,
          ignore.case = TRUE
        )
      ) %>%

      distinct(
        Run,
        Protein.Group
      ) %>%

      mutate(
        Detected = 1L
      )


    if (nrow(detection_summary) == 0) {

      return(
        matrix(
          0L,
          nrow = 0,
          ncol = length(all_runs),
          dimnames = list(
            character(0),
            all_runs
          )
        )
      )
    }


    det_mat <- detection_summary %>%

      pivot_wider(
        names_from = Run,
        values_from = Detected,
        values_fill = 0L
      ) %>%

      column_to_rownames(
        "Protein.Group"
      ) %>%

      as.matrix()


    storage.mode(det_mat) <- "integer"


    # Add ROIs with no FDR-passing proteins
    missing_runs <- setdiff(
      all_runs,
      colnames(det_mat)
    )


    if (length(missing_runs) > 0) {

      add_mat <- matrix(
        0L,
        nrow = nrow(det_mat),
        ncol = length(missing_runs),
        dimnames = list(
          rownames(det_mat),
          missing_runs
        )
      )

      det_mat <- cbind(
        det_mat,
        add_mat
      )
    }


    det_mat <- det_mat[
      ,
      all_runs,
      drop = FALSE
    ]


    return(det_mat)
  }

  make_protein_matrix <- function(df, all_runs) {

    pg_summary <- df %>%

      filter(
        Q.Value <= 0.01,
        PG.Q.Value <= 0.01,
        !is.na(Protein.Group),
        !grepl(
          "cRAP|CON__|CONTAMINANT",
          Protein.Group,
          ignore.case = TRUE
        )
      ) %>%

      group_by(
        Run,
        Protein.Group
      ) %>%

      summarise(
        Intensity = {
          x <- PG.MaxLFQ
          if (all(is.na(x))) {
            NA_real_
          } else {
            max(x, na.rm = TRUE)
          }
        },
        .groups = "drop"
      ) %>%

      mutate(
        Intensity = ifelse(
          !is.finite(Intensity) | Intensity <= 0,
          NA_real_,
          Intensity
        )
      )


    expr_mat <- pg_summary %>%

      pivot_wider(
        names_from = Run,
        values_from = Intensity
      ) %>%

      column_to_rownames("Protein.Group") %>%

      as.matrix()


    storage.mode(expr_mat) <- "numeric"


    # Numerical ROI order
    roi_order <- order(
      as.integer(
        sub(
          "ROI",
          "",
          colnames(expr_mat)
        )
      )
    )

    expr_mat <- expr_mat[
      ,
      roi_order,
      drop = FALSE
    ]


    # Add ROIs missing completely from the PG.MaxLFQ matrix as NA columns
    missing_runs <- setdiff(
      all_runs,
      colnames(expr_mat)
    )

    if (length(missing_runs) > 0) {

      add_mat <- matrix(
        NA_real_,
        nrow = nrow(expr_mat),
        ncol = length(missing_runs),
        dimnames = list(
          rownames(expr_mat),
          missing_runs
        )
      )

      expr_mat <- cbind(
        expr_mat,
        add_mat
      )
    }


    expr_mat <- expr_mat[
      ,
      all_runs,
      drop = FALSE
    ]


    return(expr_mat)
  }


  run_roi_qc <- function(file_path, label) {

    cat("\n=====================================================\n")
    cat("QC DATASET:", label, "\n")
    cat("FILE:", file_path, "\n")
    cat("=====================================================\n")


    run_map <- if (identical(label, "First-pass")) FIRST_RUN_MAP else FINAL_RUN_MAP
    source_view <- if (identical(label, "First-pass")) "qc_first_input" else "qc_final_input"
    df <- qc_read_report(QC_CON, source_view, run_map)
    check_required_columns(df, label)
    df <- prepare_diann_data(df)
    # Use every mapped ROI, including ROIs with no FDR-passing proteins.
    detected_rois <- sort(unique(run_map$ROI_Number))
    if (nrow(df) == 0L) stop("[", dataset_id, " / ", label,
      "] No protein groups passed the original FDR/contaminant filters.", call. = FALSE)
    cat("FDR-passing ROI/Protein.Group aggregate rows:", nrow(df), "\n")
    cat("ROIs detected:", length(detected_rois), "\n")

    if (length(detected_rois) > 0) {
      cat(
        "ROI range:",
        min(detected_rois),
        "-",
        max(detected_rois),
        "\n"
      )
    }


    # Use the same explicit ROI order for both matrices
    all_runs <- paste0("ROI", detected_rois)


    # Presence/absence matrix:
    # used for protein count and missingness in BOTH First-pass and Final-MBR
    detection_mat <- make_detection_matrix(
      df,
      all_runs
    )


    # PG.MaxLFQ quantitative matrix:
    # primarily used for Final-MBR abundance/correlation.
    # First-pass PG.MaxLFQ may be entirely NA.
    expr_mat <- make_protein_matrix(
      df,
      all_runs
    )


    if (length(all_runs) < 2) {
      stop(
        paste0(
          "[", label,
          "] At least 2 ROIs are required for QC."
        )
      )
    }


    log2_expr <- log2(
      expr_mat
    )


    total_detected_proteins <- nrow(
      detection_mat
    )


    # ---------------------------------------------------------------------------
    # ROI QC metrics
    # ---------------------------------------------------------------------------

    # Identified protein count:
    # based on FDR-passing Protein.Group presence, NOT PG.MaxLFQ.
    protein_count <- colSums(
      detection_mat > 0
    )


    # Identification missingness:
    # comparable between First-pass and Final-MBR.
    if (nrow(detection_mat) > 0) {

      missing_pct <- colMeans(
        detection_mat == 0
      ) * 100

    } else {

      missing_pct <- rep(
        100,
        length(all_runs)
      )

      names(missing_pct) <- all_runs
    }


    # Total PG.MaxLFQ intensity.
    # If First-pass contains no PG.MaxLFQ values, leave this as NA.
    total_intensity <- vapply(
      seq_len(
        ncol(expr_mat)
      ),
      function(j) {

        x <- expr_mat[, j]

        if (all(is.na(x))) {
          return(NA_real_)
        }

        sum(
          x,
          na.rm = TRUE
        )
      },
      numeric(1)
    )

    names(total_intensity) <- colnames(
      expr_mat
    )


    # ---------------------------------------------------------------------------
    # Correlation matrix
    #
    # IMPORTANT:
    #   cor_mat_plot = heatmap/clustering용. 대각선은 1로 유지.
    #   cor_mat_qc   = mean correlation 계산용. 자기상관(대각선)만 NA 처리.
    #
    # 이렇게 분리하지 않고 대각선까지 NA로 만든 matrix를 pheatmap에 넘기면,
    # off-diagonal correlation도 모두 NA인 경우 pheatmap이
    # "from must be a finite number" 오류를 발생시킬 수 있음.
    # ---------------------------------------------------------------------------

    finite_quant_values <- sum(
      is.finite(
        log2_expr
      ),
      na.rm = TRUE
    )


    if (finite_quant_values > 0) {

      cor_mat_plot <- suppressWarnings(
        cor(
          log2_expr,
          use = "pairwise.complete.obs",
          method = "pearson"
        )
      )

    } else {

      # First-pass may have no PG.MaxLFQ values at all.
      cor_mat_plot <- matrix(
        NA_real_,
        nrow = length(all_runs),
        ncol = length(all_runs),
        dimnames = list(
          all_runs,
          all_runs
        )
      )
    }


    cor_mat_qc <- cor_mat_plot

    if (nrow(cor_mat_qc) > 0) {
      diag(cor_mat_qc) <- NA_real_
    }


    mean_correlation <- rowMeans(
      cor_mat_qc,
      na.rm = TRUE
    )


    mean_correlation[
      !is.finite(mean_correlation)
    ] <- NA_real_


    qc_stats <- data.frame(

      Sample = colnames(expr_mat),

      ROI_Number = as.integer(
        sub(
          "ROI",
          "",
          colnames(expr_mat)
        )
      ),

      Identified_Proteins = protein_count,

      Protein_Coverage_Pct =
        protein_count /
        total_detected_proteins *
        100,

      Missing_Percentage = missing_pct,

      Total_Intensity = total_intensity,

      Mean_Correlation = mean_correlation,

      stringsAsFactors = FALSE
    )


    # ---------------------------------------------------------------------------
    # Primary QC thresholds
    # ---------------------------------------------------------------------------

    median_protein <- median(
      qc_stats$Identified_Proteins,
      na.rm = TRUE
    )

    nmad_protein <- safe_nmad(
      qc_stats$Identified_Proteins
    )

    fail_lower_cutoff <-
      median_protein -
      FAIL_NMAD *
      nmad_protein

    warning_lower_cutoff <-
      median_protein -
      WARNING_NMAD *
      nmad_protein

    high_warning_cutoff <-
      median_protein +
      FAIL_NMAD *
      nmad_protein


    qc_stats <- qc_stats %>%

      mutate(

        Protein_NMAD_Score =
          (
            median_protein -
            Identified_Proteins
          ) /
          nmad_protein,

        Primary_QC = case_when(

          Identified_Proteins < fail_lower_cutoff
          ~ "FAIL_LowDepth",

          Identified_Proteins < warning_lower_cutoff
          ~ "WARNING_LowDepth",

          Identified_Proteins > high_warning_cutoff
          ~ "WARNING_HighDepth",

          TRUE
          ~ "PASS"
        )
      )


    # ---------------------------------------------------------------------------
    # Secondary QC
    # ---------------------------------------------------------------------------

    intensity_available <- any(
      is.finite(
        qc_stats$Total_Intensity
      )
    )


    correlation_available <- any(
      is.finite(
        qc_stats$Mean_Correlation
      )
    )


    if (intensity_available) {

      median_intensity <- median(
        qc_stats$Total_Intensity,
        na.rm = TRUE
      )

      nmad_intensity <- safe_nmad(
        qc_stats$Total_Intensity
      )

    } else {

      median_intensity <- NA_real_
      nmad_intensity <- NA_real_
    }


    median_missing <- median(
      qc_stats$Missing_Percentage,
      na.rm = TRUE
    )

    nmad_missing <- safe_nmad(
      qc_stats$Missing_Percentage
    )


    if (correlation_available) {

      median_cor <- median(
        qc_stats$Mean_Correlation,
        na.rm = TRUE
      )

      nmad_cor <- safe_nmad(
        qc_stats$Mean_Correlation
      )

    } else {

      median_cor <- NA_real_
      nmad_cor <- NA_real_
    }


    qc_stats <- qc_stats %>%

      mutate(

        # IMPORTANT:
        # intensity_available / correlation_available are scalar TRUE/FALSE values.
        # Do NOT use scalar ifelse(test, vector, FALSE) here, because R returns
        # length(test) and can recycle the first ROI result across all ROIs.
        #
        # Scalar logical AND vector logical gives the intended ROI-wise result.

        Secondary_LowIntensity =
          intensity_available &
          is.finite(Total_Intensity) &
          Total_Intensity <
          (
            median_intensity -
            SECONDARY_NMAD *
            nmad_intensity
          ),

        Secondary_HighMissing =
          Missing_Percentage >
          (
            median_missing +
            SECONDARY_NMAD *
            nmad_missing
          ),

        Secondary_LowCorrelation =
          correlation_available &
          is.finite(Mean_Correlation) &
          Mean_Correlation <
          (
            median_cor -
            SECONDARY_NMAD *
            nmad_cor
          ),

        Secondary_Warning =
          Secondary_LowIntensity |
          Secondary_HighMissing |
          Secondary_LowCorrelation
      )


    qc_stats$Secondary_Warning_Reason <- apply(

      qc_stats,

      1,

      function(x) {

        reason <- character(0)

        if (isTRUE(as.logical(x["Secondary_LowIntensity"]))) {
          reason <- c(reason, "Low intensity")
        }

        if (isTRUE(as.logical(x["Secondary_HighMissing"]))) {
          reason <- c(reason, "High missingness")
        }

        if (isTRUE(as.logical(x["Secondary_LowCorrelation"]))) {
          reason <- c(reason, "Low correlation")
        }

        if (length(reason) == 0) {
          return("None")
        }

        paste(
          reason,
          collapse = "; "
        )
      }
    )


    # Dataset-level status only.
    # Final cross-report decision is calculated later.
    qc_stats <- qc_stats %>%

      mutate(

        Dataset_Auto_QC =
          ifelse(
            Primary_QC == "FAIL_LowDepth",
            "FAIL",
            "RETAIN"
          )
      )


    # ---------------------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------------------

    finite_pgmaxlfq_cells <- sum(
      is.finite(
        expr_mat
      ),
      na.rm = TRUE
    )

    cat("\nPrimary QC summary:", label, "\n")

    cat(
      "Finite PG.MaxLFQ matrix cells:",
      finite_pgmaxlfq_cells,
      "/",
      length(expr_mat),
      "\n"
    )

    if (finite_pgmaxlfq_cells == 0) {

      cat(
        "NOTE:",
        label,
        "contains no usable PG.MaxLFQ values. Protein count and missingness are calculated from FDR-passing Protein.Group detections instead.\n"
      )
    }
    cat("Median proteins         :", round(median_protein, 1), "\n")
    cat("NMAD                    :", round(nmad_protein, 1), "\n")
    cat("FAIL cutoff (-3 NMAD)   :", round(fail_lower_cutoff, 1), "\n")
    cat("WARNING cutoff (-2.5)   :", round(warning_lower_cutoff, 1), "\n")
    cat("High warning (+3 NMAD)  :", round(high_warning_cutoff, 1), "\n")

    cat("\nSecondary QC thresholds:", label, "\n")

    if (intensity_available) {
      cat(
        "Low intensity cutoff     :",
        signif(
          median_intensity -
            SECONDARY_NMAD *
            nmad_intensity,
          6
        ),
        "\n"
      )
    } else {
      cat("Low intensity cutoff     : NA (metric unavailable)\n")
    }

    cat(
      "High missingness cutoff :",
      round(
        median_missing +
          SECONDARY_NMAD *
          nmad_missing,
        3
      ),
      "%\n"
    )

    if (correlation_available) {
      cat(
        "Low correlation cutoff :",
        round(
          median_cor -
            SECONDARY_NMAD *
            nmad_cor,
          6
        ),
        "\n"
      )
    } else {
      cat("Low correlation cutoff : NA (metric unavailable)\n")
    }

    print(
      table(
        qc_stats$Primary_QC,
        useNA = "ifany"
      )
    )


    return(
      list(

        label = label,

        detection_mat = detection_mat,

        expr_mat = expr_mat,

        cor_mat = cor_mat_plot,
        cor_mat_qc = cor_mat_qc,

        qc_stats = qc_stats,

        cutoffs = list(
          median_protein = median_protein,
          nmad_protein = nmad_protein,
          fail_lower_cutoff = fail_lower_cutoff,
          warning_lower_cutoff = warning_lower_cutoff,
          high_warning_cutoff = high_warning_cutoff
        )
      )
    )
  }


  prefix_qc_columns <- function(qc_df, prefix) {

    qc_df %>%

      rename_with(
        ~ paste0(prefix, .x),
        -c(Sample, ROI_Number)
      )
  }


  # ------------------------------------------------------------------------------
  # Correlation heatmap helper functions
  #
  # Goals:
  #   1) First-pass / Final-MBR heatmap에 동일한 ROI 순서를 사용
  #   2) Final-MBR hierarchical clustering 순서를 기본 reference로 사용
  #   3) Final-MBR clustering이 불가능하면 First-pass로 fallback
  #   4) 둘 다 불가능하면 ROI number 순서 사용
  #   5) First-pass / Final-MBR 각각 PDF + PNG 자동 저장
  #   6) 두 heatmap을 나란히 배치한 comparison PDF + PNG 자동 저장
  #   7) 실제 correlation matrix의 NA는 0으로 변환하지 않음
  #      (clustering distance 계산용 copy에서만 NA -> 0)
  # ------------------------------------------------------------------------------


  # ------------------------------------------------------------------------------
  # A. Off-diagonal finite correlation 개수
  # ------------------------------------------------------------------------------

  count_finite_offdiag <- function(cor_mat) {

    if (is.null(cor_mat) || length(cor_mat) == 0) {
      return(0L)
    }

    x <- as.matrix(cor_mat)

    if (nrow(x) < 2 || ncol(x) < 2) {
      return(0L)
    }

    diag(x) <- NA_real_

    sum(
      is.finite(x),
      na.rm = TRUE
    )
  }


  # ------------------------------------------------------------------------------
  # B. ROI 이름을 숫자 순서로 정렬
  # ------------------------------------------------------------------------------

  sort_roi_names_numeric <- function(x) {

    x[
      order(
        suppressWarnings(
          as.integer(
            sub(
              "^ROI",
              "",
              x
            )
          )
        )
      )
    ]
  }


  # ------------------------------------------------------------------------------
  # C. First-pass와 Final-MBR에 공통으로 존재하는 ROI의 공통 순서 결정
  #
  # Default:
  #   Final-MBR correlation matrix에서 hierarchical clustering
  #
  # Fallback:
  #   First-pass -> ROI number
  # ------------------------------------------------------------------------------

  get_common_roi_order <- function(
      first_cor,
      final_cor,
      reference = "Final-MBR"
  ) {

    first_cor <- as.matrix(first_cor)
    final_cor <- as.matrix(final_cor)

    common_samples <- intersect(
      rownames(first_cor),
      rownames(final_cor)
    )

    if (length(common_samples) == 0) {

      stop(
        paste0(
          "No common ROI names were found between First-pass and Final-MBR ",
          "correlation matrices."
        )
      )
    }

    common_samples <- sort_roi_names_numeric(
      common_samples
    )


    if (length(common_samples) == 1) {

      warning(
        "Only one common ROI was found. Hierarchical clustering is not possible."
      )

      return(
        list(
          order = common_samples,
          reference_used = "ROI-number",
          n_common = 1L
        )
      )
    }


    first_sub <- first_cor[
      common_samples,
      common_samples,
      drop = FALSE
    ]

    final_sub <- final_cor[
      common_samples,
      common_samples,
      drop = FALSE
    ]


    try_cluster_order <- function(mat) {

      mat <- as.matrix(mat)

      # ROI 간 실제 finite correlation이 없으면 clustering 불가
      if (count_finite_offdiag(mat) == 0) {
        return(NULL)
      }

      # Clustering용 copy
      # 실제 plotting matrix는 수정하지 않음
      cluster_mat <- mat

      diag(cluster_mat) <- 1

      cluster_mat[
        !is.finite(cluster_mat)
      ] <- 0

      d <- as.dist(
        1 - cluster_mat
      )

      # 모든 distance가 finite인지 확인
      if (!all(is.finite(d))) {
        return(NULL)
      }

      hc <- hclust(
        d,
        method = "complete"
      )

      hc$labels[
        hc$order
      ]
    }


    # Requested reference 우선
    if (identical(reference, "Final-MBR")) {

      roi_order <- try_cluster_order(
        final_sub
      )

      if (!is.null(roi_order)) {

        return(
          list(
            order = roi_order,
            reference_used = "Final-MBR",
            n_common = length(common_samples)
          )
        )
      }


      warning(
        paste0(
          "Final-MBR correlation matrix could not provide a valid common ",
          "hierarchical-clustering order. Falling back to First-pass."
        )
      )


      roi_order <- try_cluster_order(
        first_sub
      )

      if (!is.null(roi_order)) {

        return(
          list(
            order = roi_order,
            reference_used = "First-pass",
            n_common = length(common_samples)
          )
        )
      }

    } else {

      roi_order <- try_cluster_order(
        first_sub
      )

      if (!is.null(roi_order)) {

        return(
          list(
            order = roi_order,
            reference_used = "First-pass",
            n_common = length(common_samples)
          )
        )
      }


      warning(
        paste0(
          "First-pass correlation matrix could not provide a valid common ",
          "hierarchical-clustering order. Falling back to Final-MBR."
        )
      )


      roi_order <- try_cluster_order(
        final_sub
      )

      if (!is.null(roi_order)) {

        return(
          list(
            order = roi_order,
            reference_used = "Final-MBR",
            n_common = length(common_samples)
          )
        )
      }
    }


    warning(
      paste0(
        "Neither First-pass nor Final-MBR produced enough finite ROI-to-ROI ",
        "correlations for hierarchical clustering. ROI-number order will be used."
      )
    )


    list(
      order = common_samples,
      reference_used = "ROI-number",
      n_common = length(common_samples)
    )
  }


  # ------------------------------------------------------------------------------
  # D. 같은 ROI 순서를 강제로 사용하는 pheatmap object 생성
  #
  # cluster_rows = FALSE
  # cluster_cols = FALSE
  #
  # 이유:
  #   ROI order를 이미 공통 reference에서 한 번 결정했으므로,
  #   각각의 heatmap에서 다시 clustering하면 순서가 달라질 수 있습니다.
  # ------------------------------------------------------------------------------

  make_fixed_order_correlation_heatmap <- function(
      cor_mat,
      roi_order,
      annotation_col = NULL,
      title = "ROI Pearson Correlation"
  ) {

    if (is.null(cor_mat) || length(cor_mat) == 0) {

      warning(
        paste0(
          title,
          ": correlation matrix is empty. Heatmap skipped."
        )
      )

      return(NULL)
    }


    cor_mat <- as.matrix(
      cor_mat
    )


    available_order <- roi_order[
      roi_order %in%
        rownames(cor_mat) &
        roi_order %in%
        colnames(cor_mat)
    ]


    if (length(available_order) < 2) {

      warning(
        paste0(
          title,
          ": fewer than 2 ROIs are available after common-order filtering. ",
          "Heatmap skipped."
        )
      )

      return(NULL)
    }


    plot_mat <- cor_mat[
      available_order,
      available_order,
      drop = FALSE
    ]


    n_finite_offdiag <- count_finite_offdiag(
      plot_mat
    )


    cat(
      "\n",
      title,
      "| ROIs:",
      length(available_order),
      "| finite off-diagonal correlations:",
      n_finite_offdiag,
      "\n"
    )


    if (n_finite_offdiag == 0) {

      warning(
        paste0(
          title,
          ": no finite ROI-to-ROI correlations were calculated. ",
          "Heatmap skipped."
        )
      )

      return(NULL)
    }


    # 자기상관은 1로 유지
    diag(plot_mat) <- 1


    if (!is.null(annotation_col)) {

      annotation_col <- annotation_col[
        available_order,
        ,
        drop = FALSE
      ]
    }


    pheatmap::pheatmap(
      plot_mat,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      annotation_col = annotation_col,
      na_col = "grey90",
      border_color = NA,
      fontsize_row = 5,   # y-axis ROI label
      fontsize_col = 5,   # x-axis ROI label
      main = title,
      silent = TRUE
    )
  }


  # ------------------------------------------------------------------------------
  # E. pheatmap object를 PDF / PNG로 저장
  # ------------------------------------------------------------------------------

  save_pheatmap_pdf_png <- function(
      ph,
      file_stem,
      width,
      height,
      dpi = 300
  ) {

    if (is.null(ph)) {

      warning(
        paste0(
          file_stem,
          ": heatmap object is NULL. File saving skipped."
        )
      )

      return(invisible(FALSE))
    }


    pdf_file <- paste0(
      file_stem,
      ".pdf"
    )

    png_file <- paste0(
      file_stem,
      ".png"
    )


    # PDF -----------------------------------------------------------------------
    grDevices::pdf(
      pdf_file,
      width = width,
      height = height,
      onefile = FALSE
    )

    grid::grid.newpage()
    grid::grid.draw(
      ph$gtable
    )

    grDevices::dev.off()


    # PNG -----------------------------------------------------------------------
    grDevices::png(
      png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi
    )

    grid::grid.newpage()
    grid::grid.draw(
      ph$gtable
    )

    grDevices::dev.off()


    cat(
      "Saved:",
      pdf_file,
      "\n"
    )

    cat(
      "Saved:",
      png_file,
      "\n"
    )


    invisible(TRUE)
  }


  # ------------------------------------------------------------------------------
  # F. 두 pheatmap object를 같은 페이지에 나란히 저장
  # ------------------------------------------------------------------------------

  draw_two_pheatmaps_side_by_side <- function(
      ph_left,
      ph_right
  ) {

    grid::grid.newpage()

    grid::pushViewport(
      grid::viewport(
        layout = grid::grid.layout(
          nrow = 1,
          ncol = 2,
          widths = grid::unit(
            c(1, 1),
            "null"
          )
        )
      )
    )


    # Left: First-pass ----------------------------------------------------------
    grid::pushViewport(
      grid::viewport(
        layout.pos.row = 1,
        layout.pos.col = 1
      )
    )

    grid::grid.draw(
      ph_left$gtable
    )

    grid::upViewport()


    # Right: Final-MBR ----------------------------------------------------------
    grid::pushViewport(
      grid::viewport(
        layout.pos.row = 1,
        layout.pos.col = 2
      )
    )

    grid::grid.draw(
      ph_right$gtable
    )

    grid::upViewport()


    grid::upViewport()
  }


  save_side_by_side_pdf_png <- function(
      ph_first,
      ph_final,
      file_stem,
      width,
      height,
      dpi = 300
  ) {

    if (
      is.null(ph_first) ||
      is.null(ph_final)
    ) {

      warning(
        paste0(
          file_stem,
          ": one or both heatmap objects are NULL. ",
          "Side-by-side comparison file saving skipped."
        )
      )

      return(invisible(FALSE))
    }


    pdf_file <- paste0(
      file_stem,
      ".pdf"
    )

    png_file <- paste0(
      file_stem,
      ".png"
    )


    # PDF -----------------------------------------------------------------------
    grDevices::pdf(
      pdf_file,
      width = width,
      height = height,
      onefile = FALSE
    )

    draw_two_pheatmaps_side_by_side(
      ph_first,
      ph_final
    )

    grDevices::dev.off()


    # PNG -----------------------------------------------------------------------
    grDevices::png(
      png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi
    )

    draw_two_pheatmaps_side_by_side(
      ph_first,
      ph_final
    )

    grDevices::dev.off()


    cat(
      "Saved:",
      pdf_file,
      "\n"
    )

    cat(
      "Saved:",
      png_file,
      "\n"
    )


    invisible(TRUE)
  }


  # ==============================================================================
  # 4. Run First-pass QC and Final-MBR QC independently
  # ==============================================================================

  first_qc <- run_roi_qc(
    FIRST_PASS_FILE,
    "First-pass"
  )

  final_qc <- run_roi_qc(
    FINAL_MBR_FILE,
    "Final-MBR"
  )


  # ==============================================================================
  # 4-1. Correlation diagnostics
  # ==============================================================================

  cat("\n=====================================================\n")
  cat("CORRELATION DIAGNOSTICS\n")
  cat("=====================================================\n")

  cat(
    "First-pass finite off-diagonal correlations:",
    sum(
      is.finite(first_qc$cor_mat_qc),
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Final-MBR finite off-diagonal correlations:",
    sum(
      is.finite(final_qc$cor_mat_qc),
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "First-pass proteins per ROI range:",
    min(colSums(!is.na(first_qc$expr_mat))),
    "-",
    max(colSums(!is.na(first_qc$expr_mat))),
    "\n"
  )

  cat(
    "Final-MBR proteins per ROI range:",
    min(colSums(!is.na(final_qc$expr_mat))),
    "-",
    max(colSums(!is.na(final_qc$expr_mat))),
    "\n"
  )

  cat("=====================================================\n")


  # ==============================================================================
  # 5. Merge First-pass and Final-MBR QC
  # ==============================================================================

  first_compare <- prefix_qc_columns(
    first_qc$qc_stats,
    "First_"
  )

  final_compare <- prefix_qc_columns(
    final_qc$qc_stats,
    "MBR_"
  )


  qc_compare <- full_join(
    first_compare,
    final_compare,
    by = c(
      "Sample",
      "ROI_Number"
    )
  ) %>%

    arrange(
      ROI_Number
    )


  # ==============================================================================
  # 6. Cross-report consensus / review classification
  # ==============================================================================

  qc_compare <- qc_compare %>%

    mutate(

      # Did each dataset mark this ROI as a primary low-depth failure?
      First_Is_FAIL =
        First_Primary_QC == "FAIL_LowDepth",

      MBR_Is_FAIL =
        MBR_Primary_QC == "FAIL_LowDepth",


      Cross_Report_Status = case_when(

        # Missing from either report: never auto-exclude
        is.na(First_Primary_QC) |
          is.na(MBR_Primary_QC)
        ~ "REVIEW_MISSING_REPORT",

        # Strongest evidence: fails independently in both analyses
        First_Is_FAIL &
          MBR_Is_FAIL
        ~ "FAIL_CONSENSUS",

        # First-pass low depth, but MBR restores it above FAIL cutoff
        First_Is_FAIL &
          !MBR_Is_FAIL
        ~ "REVIEW_MBR_RESCUED",

        # Rare discordance in opposite direction
        !First_Is_FAIL &
          MBR_Is_FAIL
        ~ "REVIEW_FINAL_FAIL",

        # Primary or secondary warnings without a FAIL
        First_Primary_QC != "PASS" |
          MBR_Primary_QC != "PASS" |
          First_Secondary_Warning |
          MBR_Secondary_Warning
        ~ "REVIEW_WARNING",

        TRUE
        ~ "RETAIN"
      ),


      # Automatic action:
      # Only consensus failures are automatically excluded.
      QC_Action = case_when(

        Cross_Report_Status == "FAIL_CONSENSUS"
        ~ "EXCLUDE",

        grepl(
          "^REVIEW",
          Cross_Report_Status
        )
        ~ "REVIEW",

        TRUE
        ~ "RETAIN"
      ),


      # MBR effect on identified protein count
      Protein_Count_Delta_MBR =
        MBR_Identified_Proteins -
        First_Identified_Proteins,

      Protein_Count_Change_Pct =
        ifelse(
          !is.na(First_Identified_Proteins) &
            First_Identified_Proteins > 0,

          (
            MBR_Identified_Proteins -
            First_Identified_Proteins
          ) /
            First_Identified_Proteins *
            100,

          NA_real_
        ),

      Missing_Delta_MBR =
        MBR_Missing_Percentage -
        First_Missing_Percentage
    )


  # ==============================================================================
  # 7. Apply optional manual exclusion after REVIEW
  # ==============================================================================

  qc_compare <- qc_compare %>%

    mutate(

      Manual_Exclude =
        Sample %in% MANUAL_EXCLUDE,

      Final_QC_Action = case_when(

        Manual_Exclude
        ~ "EXCLUDE_MANUAL",

        QC_Action == "EXCLUDE"
        ~ "EXCLUDE",

        QC_Action == "REVIEW"
        ~ "REVIEW",

        TRUE
        ~ "RETAIN"
      )
    )


  # ==============================================================================
  # 8. Cross-report summary
  # ==============================================================================

  cat("\n=====================================================\n")
  cat("FIRST-PASS vs FINAL-MBR QC SUMMARY\n")
  cat("=====================================================\n")

  print(
    table(
      qc_compare$Cross_Report_Status,
      useNA = "ifany"
    )
  )

  cat("\nFinal QC action:\n")

  print(
    table(
      qc_compare$Final_QC_Action,
      useNA = "ifany"
    )
  )

  cat("\nROI transitions:\n")

  print(
    table(
      First = qc_compare$First_Primary_QC,
      Final_MBR = qc_compare$MBR_Primary_QC,
      useNA = "ifany"
    )
  )

  cat("=====================================================\n")


  # Show ROIs requiring attention
  review_table <- qc_compare %>%

    filter(
      Final_QC_Action != "RETAIN"
    ) %>%

    select(
      Sample,
      ROI_Number,
      First_Identified_Proteins,
      MBR_Identified_Proteins,
      Protein_Count_Delta_MBR,
      Protein_Count_Change_Pct,
      First_Primary_QC,
      MBR_Primary_QC,
      First_Secondary_Warning_Reason,
      MBR_Secondary_Warning_Reason,
      Cross_Report_Status,
      Final_QC_Action
    )

  print(
    as_tibble(review_table),
    n = Inf
  )


  # ==============================================================================
  # 9. First-pass vs Final-MBR comparison plots
  # ==============================================================================

  # --------------------------------------------------------------------------
  # 9-1. Identified protein count: First-pass vs Final-MBR
  # --------------------------------------------------------------------------

  p_protein_compare <- ggplot(
    qc_compare,
    aes(
      x = First_Identified_Proteins,
      y = MBR_Identified_Proteins,
      shape = Cross_Report_Status
    )
  ) +

    geom_point(
      size = 3
    ) +

    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed"
    ) +

    theme_bw() +

    labs(
      title = "First-pass vs Final-MBR: Identified Protein Groups",
      subtitle = "Dashed line = y = x",
      x = "First-pass identified proteins",
      y = "Final-MBR identified proteins",
      shape = "Cross-report QC"
    )


  print(
    p_protein_compare
  )


  # --------------------------------------------------------------------------
  # 9-2. MBR change in identified proteins per ROI
  # --------------------------------------------------------------------------

  p_mbr_delta <- ggplot(
    qc_compare,
    aes(
      x = reorder(
        Sample,
        Protein_Count_Delta_MBR
      ),
      y = Protein_Count_Delta_MBR,
      fill = Cross_Report_Status
    )
  ) +

    geom_col() +

    coord_flip() +

    theme_bw() +
  
    theme(
      axis.text.y = element_text(size = 5)
    ) +

    labs(
      title = "MBR Effect on Protein Identification",
      subtitle = "Final-MBR minus First-pass identified proteins",
      x = "ROI",
      y = "Delta identified proteins",
      fill = "Cross-report QC"
    )


  print(
    p_mbr_delta
  )


  # --------------------------------------------------------------------------
  # 9-3. First-pass and Final-MBR missingness
  # --------------------------------------------------------------------------

  missing_long <- qc_compare %>%

    select(
      Sample,
      ROI_Number,
      First_Missing_Percentage,
      MBR_Missing_Percentage,
      Cross_Report_Status
    ) %>%

    pivot_longer(
      cols = c(
        First_Missing_Percentage,
        MBR_Missing_Percentage
      ),
      names_to = "Analysis",
      values_to = "Missing_Percentage"
    ) %>%

    mutate(
      Analysis = recode(
        Analysis,
        First_Missing_Percentage = "First-pass",
        MBR_Missing_Percentage = "Final-MBR"
      )
    )


  p_missing_compare <- ggplot(
    missing_long,
    aes(
      x = ROI_Number,
      y = Missing_Percentage,
      shape = Analysis,
      color = Analysis
    )
  ) +

    geom_point(
      size = 2.8,
      alpha = 0.9
    ) +

    # First-pass와 Final-MBR을 색상과 모양 모두로 구분
    scale_shape_manual(
      values = c(
        "First-pass" = 16,
        "Final-MBR" = 17
      )
    ) +

    scale_color_manual(
      values = c(
        "First-pass" = "#0072B2",
        "Final-MBR" = "#D55E00"
      )
    ) +

    theme_bw() +

    labs(
      title = "First-pass vs Final-MBR Missingness",
      x = "ROI",
      y = "Missing (%)",
      shape = "Analysis",
      color = "Analysis"
    ) +

    guides(
      shape = guide_legend(
        order = 1
      ),
      color = guide_legend(
        order = 1
      )
    )


  print(
    p_missing_compare
  )


  # --------------------------------------------------------------------------
  # 9-4. QC status counts
  # --------------------------------------------------------------------------

  status_count <- qc_compare %>%

    count(
      Cross_Report_Status,
      name = "N_ROIs"
    ) %>%

    arrange(
      N_ROIs
    ) %>%

    mutate(

      # 범례에 각 QC classification의 ROI 수를 함께 표시
      # 예: FAIL_CONSENSUS (n=3)
      Legend_Label = paste0(
        Cross_Report_Status,
        " (n=",
        N_ROIs,
        ")"
      ),

      # Bar와 legend의 표시 순서를 동일하게 유지
      Cross_Report_Status = factor(
        Cross_Report_Status,
        levels = Cross_Report_Status
      ),

      Legend_Label = factor(
        Legend_Label,
        levels = Legend_Label
      )
    )


  p_status <- ggplot(
    status_count,
    aes(
      x = Cross_Report_Status,
      y = N_ROIs,
      fill = Legend_Label
    )
  ) +

    geom_col() +

    coord_flip() +

    theme_bw() +

    labs(
      title = "Cross-report ROI QC Classification",
      x = "QC classification",
      y = "Number of ROIs",
      fill = "QC classification (Number of ROIs)"
    )


  print(
    p_status
  )


  # ==============================================================================
  # 10. Correlation heatmaps
  #
  # IMPORTANT:
  #   First-pass와 Final-MBR heatmap의 ROI 순서를 동일하게 맞춥니다.
  #
  # Default ordering:
  #   Final-MBR correlation matrix hierarchical clustering
  #
  # Fallback:
  #   First-pass clustering -> numerical ROI order
  #
  # Outputs:
  #   1) First-pass PDF
  #   2) First-pass PNG
  #   3) Final-MBR PDF
  #   4) Final-MBR PNG
  #   5) Side-by-side comparison PDF
  #   6) Side-by-side comparison PNG
  #   7) Common ROI order CSV
  # ==============================================================================


  # --------------------------------------------------------------------------
  # 10-1. Heatmap annotations
  # --------------------------------------------------------------------------

  first_annotation <- qc_compare %>%

    select(
      Sample,
      First_Primary_QC,
      Cross_Report_Status
    ) %>%

    column_to_rownames(
      "Sample"
    )


  mbr_annotation <- qc_compare %>%

    select(
      Sample,
      MBR_Primary_QC,
      Cross_Report_Status
    ) %>%

    column_to_rownames(
      "Sample"
    )


  # --------------------------------------------------------------------------
  # 10-2. Determine ONE common ROI order
  # --------------------------------------------------------------------------

  heatmap_order_info <- get_common_roi_order(

    first_cor = first_qc$cor_mat,

    final_cor = final_qc$cor_mat,

    reference = HEATMAP_ORDER_REFERENCE
  )


  common_roi_order <- heatmap_order_info$order


  cat("\n=====================================================\n")
  cat("COMMON PEARSON CORRELATION HEATMAP ORDER\n")
  cat("=====================================================\n")

  cat(
    "Common ROIs:",
    heatmap_order_info$n_common,
    "\n"
  )

  cat(
    "Order reference used:",
    heatmap_order_info$reference_used,
    "\n"
  )

  cat(
    "ROI order:\n",
    paste(
      common_roi_order,
      collapse = ", "
    ),
    "\n"
  )

  cat("=====================================================\n")


  # Save ROI order used in both heatmaps
  heatmap_order_df <- data.frame(

    Heatmap_Order =
      seq_along(
        common_roi_order
      ),

    Sample =
      common_roi_order,

    Order_Reference =
      heatmap_order_info$reference_used,

    stringsAsFactors = FALSE
  )


  write.csv(

    heatmap_order_df,

    paste0(
      OUTPUT_PREFIX,
      "_Pearson_Correlation_Common_ROI_Order.csv"
    ),

    row.names = FALSE
  )


  # --------------------------------------------------------------------------
  # 10-3. Build First-pass heatmap using the COMMON order
  # --------------------------------------------------------------------------

  ph_first <- make_fixed_order_correlation_heatmap(

    cor_mat =
      first_qc$cor_mat,

    roi_order =
      common_roi_order,

    annotation_col =
      first_annotation,

    title =
      paste0(
        "First-pass ROI Pearson Correlation\n",
        "Common ROI order: ",
        heatmap_order_info$reference_used,
        " reference"
      )
  )


  # --------------------------------------------------------------------------
  # 10-4. Build Final-MBR heatmap using EXACTLY the same order
  # --------------------------------------------------------------------------

  ph_final <- make_fixed_order_correlation_heatmap(

    cor_mat =
      final_qc$cor_mat,

    roi_order =
      common_roi_order,

    annotation_col =
      mbr_annotation,

    title =
      paste0(
        "Final-MBR ROI Pearson Correlation\n",
        "Common ROI order: ",
        heatmap_order_info$reference_used,
        " reference"
      )
  )


  # --------------------------------------------------------------------------
  # 10-5. Display both heatmaps in RStudio plot device
  # --------------------------------------------------------------------------

  if (!is.null(ph_first)) {

    grid::grid.newpage()

    grid::grid.draw(
      ph_first$gtable
    )
  }


  if (!is.null(ph_final)) {

    grid::grid.newpage()

    grid::grid.draw(
      ph_final$gtable
    )
  }


  # --------------------------------------------------------------------------
  # 10-6. Save individual First-pass heatmap: PDF + PNG
  # --------------------------------------------------------------------------

  save_pheatmap_pdf_png(

    ph =
      ph_first,

    file_stem =
      paste0(
        OUTPUT_PREFIX,
        "_FirstPass_Pearson_Correlation_SameROIOrder"
      ),

    width =
      HEATMAP_SINGLE_WIDTH,

    height =
      HEATMAP_SINGLE_HEIGHT,

    dpi =
      HEATMAP_PNG_DPI
  )


  # --------------------------------------------------------------------------
  # 10-7. Save individual Final-MBR heatmap: PDF + PNG
  # --------------------------------------------------------------------------

  save_pheatmap_pdf_png(

    ph =
      ph_final,

    file_stem =
      paste0(
        OUTPUT_PREFIX,
        "_FinalMBR_Pearson_Correlation_SameROIOrder"
      ),

    width =
      HEATMAP_SINGLE_WIDTH,

    height =
      HEATMAP_SINGLE_HEIGHT,

    dpi =
      HEATMAP_PNG_DPI
  )


  # --------------------------------------------------------------------------
  # 10-8. Save side-by-side comparison: PDF + PNG
  # --------------------------------------------------------------------------

  save_side_by_side_pdf_png(

    ph_first =
      ph_first,

    ph_final =
      ph_final,

    file_stem =
      paste0(
        OUTPUT_PREFIX,
        "_FirstPass_vs_FinalMBR_Pearson_Correlation_SideBySide"
      ),

    width =
      HEATMAP_COMBINED_WIDTH,

    height =
      HEATMAP_COMBINED_HEIGHT,

    dpi =
      HEATMAP_PNG_DPI
  )


  # ==============================================================================
  # 11. Automatic exclusion
  #
  # EXCLUDE:
  #   - FAIL_CONSENSUS
  #   - MANUAL_EXCLUDE
  #
  # REVIEW:
  #   retained automatically until manually reviewed
  # ==============================================================================

  exclude_samples <- qc_compare %>%

    filter(
      Final_QC_Action %in% c(
        "EXCLUDE",
        "EXCLUDE_MANUAL"
      )
    ) %>%

    pull(
      Sample
    )


  retained_samples <- setdiff(
    colnames(final_qc$expr_mat),
    exclude_samples
  )


  clean_expr_mat <- final_qc$expr_mat[
    ,
    retained_samples,
    drop = FALSE
  ]


  cat("\n=====================================================\n")
  cat("FINAL ROI FILTERING\n")
  cat("=====================================================\n")

  cat(
    "Final-MBR ROIs :", ncol(final_qc$expr_mat), "\n"
  )

  cat(
    "Auto/manual excluded:",
    length(exclude_samples),
    "\n"
  )

  cat(
    "Retained ROIs  :",
    ncol(clean_expr_mat),
    "\n"
  )

  if (length(exclude_samples) > 0) {
    cat(
      "Excluded samples:",
      paste(exclude_samples, collapse = ", "),
      "\n"
    )
  }

  cat("=====================================================\n")


  if (ncol(clean_expr_mat) == 0L) {
    stop("[", dataset_id, "] All ROIs were excluded. Review automatic/manual exclusions.", call. = FALSE)
  }

  # ==============================================================================
  # 12. Protein completeness after ROI QC
  # ==============================================================================

  protein_detection_count <- rowSums(
    !is.na(
      clean_expr_mat
    )
  )


  protein_completeness <- rowMeans(
    !is.na(
      clean_expr_mat
    )
  )


  protein_qc <- data.frame(

    Protein.Group =
      rownames(
        clean_expr_mat
      ),

    Detection_Count =
      protein_detection_count,

    Total_ROIs =
      ncol(
        clean_expr_mat
      ),

    Completeness =
      protein_completeness,

    Completeness_Pct =
      protein_completeness *
      100,

    stringsAsFactors = FALSE
  )


  # ==============================================================================
  # 13. >=70% completeness filtering
  # ==============================================================================

  keep_70 <-
    protein_completeness >=
    PROTEIN_COMPLETENESS


  zonation_expr_mat <- clean_expr_mat[
    keep_70,
    ,
    drop = FALSE
  ]


  cat("\n=====================================================\n")
  cat("PROTEIN COMPLETENESS FILTERING\n")
  cat("=====================================================\n")

  cat(
    "Proteins before filtering:",
    nrow(clean_expr_mat),
    "\n"
  )

  cat(
    "Completeness threshold:",
    PROTEIN_COMPLETENESS * 100,
    "%\n"
  )

  cat(
    "Proteins after filtering:",
    nrow(zonation_expr_mat),
    "\n"
  )

  cat(
    "Proteins removed:",
    sum(!keep_70),
    "\n"
  )

  cat("=====================================================\n")


  # ==============================================================================
  # 14. Protein completeness distribution
  # ==============================================================================

  p_completeness <- ggplot(
    protein_qc,
    aes(
      x = Completeness_Pct
    )
  ) +

    geom_histogram(
      bins = 40
    ) +

    geom_vline(
      xintercept =
        PROTEIN_COMPLETENESS *
        100,
      linetype = "dashed",
      linewidth = 1
    ) +

    theme_bw() +

    labs(
      title = "Protein Data Completeness after Dual ROI QC",
      subtitle = "Dashed line = 70% completeness cutoff",
      x = "Detected across QC-retained ROIs (%)",
      y = "Number of Protein Groups"
    )


  print(
    p_completeness
  )


  # ==============================================================================
  # 14-1. Export ALL QC plots
  #
  # Outputs:
  #
  # Individual 600-dpi PNG:
  #   1) First-pass vs Final-MBR protein count
  #   2) MBR protein-count delta
  #   3) First-pass vs Final-MBR missingness
  #   4) Cross-report ROI QC classification
  #   5) Protein completeness distribution
  #   6) First-pass Pearson correlation heatmap          [already generated above]
  #   7) Final-MBR Pearson correlation heatmap           [already generated above]
  #   8) First-pass vs Final-MBR Pearson side-by-side    [already generated above]
  #
  # One multipage PDF:
  #   *_All_QC_Plots.pdf
  #
  # NOTE:
  #   PDF is vector format; dpi applies only to PNG files.
  # ==============================================================================


  # ------------------------------------------------------------------------------
  # A. Helper: save ggplot as 600-dpi PNG
  # ------------------------------------------------------------------------------

  save_ggplot_png <- function(
      plot_object,
      file_name,
      width = PLOT_SINGLE_WIDTH,
      height = PLOT_SINGLE_HEIGHT,
      dpi = PLOT_PNG_DPI
  ) {

    ggplot2::ggsave(
      filename = file_name,
      plot = plot_object,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      bg = "white"
    )

    cat(
      "Saved 600-dpi PNG:",
      file_name,
      "\n"
    )
  }


  # ------------------------------------------------------------------------------
  # B. Save all ggplot objects individually at 600 dpi
  # ------------------------------------------------------------------------------

  save_ggplot_png(

    plot_object =
      p_protein_compare,

    file_name =
      paste0(
        OUTPUT_PREFIX,
        "_01_FirstPass_vs_FinalMBR_Identified_Proteins.png"
      ),

    width =
      9,

    height =
      8
  )


  save_ggplot_png(

    plot_object =
      p_mbr_delta,

    file_name =
      paste0(
        OUTPUT_PREFIX,
        "_02_MBR_Protein_Count_Delta.png"
      ),

    width =
      PLOT_TALL_WIDTH,

    height =
      PLOT_TALL_HEIGHT
  )


  save_ggplot_png(

    plot_object =
      p_missing_compare,

    file_name =
      paste0(
        OUTPUT_PREFIX,
        "_03_FirstPass_vs_FinalMBR_Missingness.png"
      ),

    width =
      12,

    height =
      7
  )


  save_ggplot_png(

    plot_object =
      p_status,

    file_name =
      paste0(
        OUTPUT_PREFIX,
        "_04_CrossReport_ROI_QC_Classification.png"
      ),

    width =
      11,

    height =
      7
  )


  save_ggplot_png(

    plot_object =
      p_completeness,

    file_name =
      paste0(
        OUTPUT_PREFIX,
        "_08_Protein_Completeness_Distribution.png"
      ),

    width =
      9,

    height =
      7
  )


  # ------------------------------------------------------------------------------
  # C. Create ONE multipage PDF containing every generated plot
  #
  # Page order:
  #   1. Identified protein comparison
  #   2. MBR protein-count delta
  #   3. Missingness comparison
  #   4. Cross-report QC classification
  #   5. First-pass Pearson heatmap
  #   6. Final-MBR Pearson heatmap
  #   7. Pearson heatmaps side-by-side
  #   8. Protein completeness
  # ------------------------------------------------------------------------------

  all_plots_pdf <- paste0(
    OUTPUT_PREFIX,
    "_All_QC_Plots.pdf"
  )


  grDevices::pdf(
    file = all_plots_pdf,
    width = ALL_PLOTS_PDF_WIDTH,
    height = ALL_PLOTS_PDF_HEIGHT,
    onefile = TRUE
  )


  # Page 1 -----------------------------------------------------------------------
  print(
    p_protein_compare
  )


  # Page 2 -----------------------------------------------------------------------
  print(
    p_mbr_delta
  )


  # Page 3 -----------------------------------------------------------------------
  print(
    p_missing_compare
  )


  # Page 4 -----------------------------------------------------------------------
  print(
    p_status
  )


  # Page 5 -----------------------------------------------------------------------
  if (!is.null(ph_first)) {

    grid::grid.newpage()

    grid::grid.draw(
      ph_first$gtable
    )

  } else {

    grid::grid.newpage()

    grid::grid.text(
      "First-pass Pearson correlation heatmap was not available."
    )
  }


  # Page 6 -----------------------------------------------------------------------
  if (!is.null(ph_final)) {

    grid::grid.newpage()

    grid::grid.draw(
      ph_final$gtable
    )

  } else {

    grid::grid.newpage()

    grid::grid.text(
      "Final-MBR Pearson correlation heatmap was not available."
    )
  }


  # Page 7 -----------------------------------------------------------------------
  if (
    !is.null(ph_first) &&
    !is.null(ph_final)
  ) {

    draw_two_pheatmaps_side_by_side(
      ph_first,
      ph_final
    )

  } else {

    grid::grid.newpage()

    grid::grid.text(
      "Side-by-side Pearson correlation heatmap was not available."
    )
  }


  # Page 8 -----------------------------------------------------------------------
  print(
    p_completeness
  )


  grDevices::dev.off()


  cat(
    "\n=====================================================\n"
  )

  cat(
    "ALL QC PLOTS EXPORTED\n"
  )

  cat(
    "=====================================================\n"
  )

  cat(
    "Multipage PDF:",
    all_plots_pdf,
    "\n"
  )

  cat(
    "Individual PNG resolution:",
    PLOT_PNG_DPI,
    "dpi\n"
  )

  cat(
    "Pearson heatmap PNG resolution:",
    HEATMAP_PNG_DPI,
    "dpi\n"
  )

  cat(
    "=====================================================\n"
  )


  # ==============================================================================
  # 15. Output files
  # ==============================================================================

  # First-pass QC
  write.csv(
    first_qc$qc_stats,
    paste0(
      OUTPUT_PREFIX,
      "_FirstPass_ROI_QC.csv"
    ),
    row.names = FALSE
  )


  # Final-MBR QC
  write.csv(
    final_qc$qc_stats,
    paste0(
      OUTPUT_PREFIX,
      "_FinalMBR_ROI_QC.csv"
    ),
    row.names = FALSE
  )


  # Direct First-pass vs Final-MBR comparison
  write.csv(
    qc_compare,
    paste0(
      OUTPUT_PREFIX,
      "_Dual_QC_Comparison.csv"
    ),
    row.names = FALSE
  )


  # ROIs requiring exclusion/review
  write.csv(
    review_table,
    paste0(
      OUTPUT_PREFIX,
      "_QC_Review_ROIs.csv"
    ),
    row.names = FALSE
  )


  # Retained Final-MBR expression matrix
  clean_expr_df <- as.data.frame(
    clean_expr_mat
  ) %>%

    rownames_to_column(
      "Protein.Group"
    )


  write.csv(
    clean_expr_df,
    paste0(
      OUTPUT_PREFIX,
      "_FinalMBR_QC_Retained_Matrix.csv"
    ),
    row.names = FALSE
  )


  # Protein completeness
  write.csv(
    protein_qc,
    paste0(
      OUTPUT_PREFIX,
      "_Protein_Completeness.csv"
    ),
    row.names = FALSE
  )


  # Final zonation matrix: >=70% completeness, no imputation
  zonation_expr_df <- as.data.frame(
    zonation_expr_mat
  ) %>%

    rownames_to_column(
      "Protein.Group"
    )


  write.csv(
    zonation_expr_df,
    paste0(
      OUTPUT_PREFIX,
      "_Zonation_70pct_NoImputation.csv"
    ),
    row.names = FALSE
  )


  # ==============================================================================
  # 16. Final summary
  # ==============================================================================

  cat("\n=====================================================\n")
  cat("FINAL DUAL-QC SUMMARY\n")
  cat("=====================================================\n")

  cat(
    "First-pass ROIs         :",
    ncol(first_qc$expr_mat),
    "\n"
  )

  cat(
    "Final-MBR ROIs          :",
    ncol(final_qc$expr_mat),
    "\n"
  )

  cat(
    "Consensus FAIL excluded :",
    sum(
      qc_compare$Cross_Report_Status ==
        "FAIL_CONSENSUS",
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "MBR-rescued REVIEW      :",
    sum(
      qc_compare$Cross_Report_Status ==
        "REVIEW_MBR_RESCUED",
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Final-only FAIL REVIEW  :",
    sum(
      qc_compare$Cross_Report_Status ==
        "REVIEW_FINAL_FAIL",
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Other warning REVIEW    :",
    sum(
      qc_compare$Cross_Report_Status ==
        "REVIEW_WARNING",
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Final retained ROIs     :",
    ncol(clean_expr_mat),
    "\n"
  )

  cat(
    "70% complete proteins   :",
    nrow(zonation_expr_mat),
    "\n"
  )

  cat("\nNo imputation was performed.\n")
  cat("=====================================================\n")

  settings <- data.frame(
    Dataset_ID = dataset_id, Condition = if (startsWith(dataset_id, "SPF")) "SPF" else "GF",
    Dataset_Number = sub("^(SPF|GF)", "", dataset_id),
    Fail_NMAD = FAIL_NMAD, Warning_NMAD = WARNING_NMAD,
    Secondary_NMAD = SECONDARY_NMAD, Protein_Completeness = PROTEIN_COMPLETENESS,
    Manual_Exclude = paste(MANUAL_EXCLUDE, collapse = ";"),
    PNG_DPI = png_dpi, stringsAsFactors = FALSE
  )
  utils::write.csv(settings, paste0(OUTPUT_PREFIX, "_QC_Settings.csv"), row.names = FALSE)
  data.frame(
    Dataset_ID = dataset_id,
    Condition = if (startsWith(dataset_id, "SPF")) "SPF" else "GF",
    Dataset_Number = sub("^(SPF|GF)", "", dataset_id),
    FirstPass_ROIs = ncol(first_qc$expr_mat),
    FinalMBR_ROIs = ncol(final_qc$expr_mat),
    Consensus_FAIL = sum(qc_compare$Cross_Report_Status == "FAIL_CONSENSUS", na.rm = TRUE),
    Manual_Excluded = sum(qc_compare$Final_QC_Action == "EXCLUDE_MANUAL", na.rm = TRUE),
    Retained_ROIs = ncol(clean_expr_mat),
    Proteins_70pct = nrow(zonation_expr_mat),
    Manual_Exclude = paste(MANUAL_EXCLUDE, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

