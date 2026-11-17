# Spatial Proteomics dashboard
# UI THEME: Bertis company blue v2 (2026-11-17)
# Install once: install.packages(c("shiny", "shinydashboard", "DT", "DBI", "duckdb", "callr", "base64enc"))
# Run: shiny::runApp()
#
# BRANDING: edit APP_CONFIG below.
# LOGO: put your company logo at www/company_logo.png beside this app.R.
# A text brand is shown automatically when that file is absent.

APP_CONFIG <- list(
  company_name = "Bertis",
  app_title = "Spatial Proteomics for Hepatic zonation",
  app_subtitle = "DIA-NN Based Analysis Workspace",
  logo_file = "company_logo.png",     # Relative to www/
  # Company color palette
  light_color = "#74BEEB",            # RGB 116 / 190 / 235
  accent_color = "#1191D0",           # RGB 17 / 145 / 208
  primary_color = "#0072BB",          # RGB 0 / 114 / 187
  sidebar_color = "#254D83",          # RGB 43 / 91 / 153
  deep_color = "#0153A5",             # RGB 1 / 83 / 165
  background_color = "#F4F8FC",
  upload_limit_gib = 10,
  db_memory_limit = "10GB",            # Buffer manager, not total R process RAM
  preview_rows = 20L,
  max_explorer_rows = 1000L,
  max_run_summary_rows = 5000L,
  qc_script = "codes/QC_Code.R",            # Absolute path or relative to app root
  comparison_script = "codes/Comparison_Pipeline.R",
  roi_metadata_script = "codes/ROI_Metadata_Generation.R",
  continuous_zonation_script = "codes/Continuous_Zonation_Total_Separation.R",
  zonation_coefficient_script = "codes/Zonation_Coefficient_Beta.R",
  trajectory_selected_script = "codes/Trajectory_Selected.R",
  spatial_expression_script = "codes/Spatial_Expression_Shiny.R",
  trajectory_cache_dir = "Trajectory_Cache",
  qc_results_dir = "Analysis_results",      # Persistent output base folder
  # Example defaults: qc_script = "D:/Spatial_Proteomics/codes/QC_Code.R"
  #                   qc_results_dir = "D:/Spatial_Proteomics/Analysis_results"
  qc_table_preview_rows = 1000L,
  qc_png_dpi = 600L
)

required_packages <- c("shiny", "shinydashboard", "DT", "DBI", "duckdb", "callr", "base64enc")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_packages)) {
  stop(paste0(
    "Install the required packages first: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "), "))"
  ), call. = FALSE)
}
library(shiny)
library(shinydashboard)

APP_ROOT <- normalizePath(".", winslash = "/", mustWork = TRUE)

options(shiny.maxRequestSize = APP_CONFIG$upload_limit_gib * 1024^10)

# Optional: set a temporary work directory on a disk with sufficient free space.
# options(spatial.temp_parent = "D:/ShinyTemp")
# No uploaded datasets are put under www/. Work files are session-specific.


# Resolve paths entered in the app. Relative paths are based on the app folder.
# Explorer's "Copy as path" quotes and either Windows slash direction are accepted.
resolve_user_path <- function(value, label = "path", base_dir = APP_ROOT) {
  if (is.null(value) || length(value) != 1L || is.na(value)) {
    stop(paste(label, " is required."), call. = FALSE)
  }
  path <- trimws(as.character(value))
  if (nchar(path) >= 2L) {
    first <- substr(path, 1L, 1L)
    last <- substr(path, nchar(path), nchar(path))
    if (first %in% c('"', "'") && identical(first, last)) {
      path <- substr(path, 2L, nchar(path) - 1L)
    }
  }
  path <- trimws(path)
  if (!nzchar(path)) stop(paste(label, " is required."), call. = FALSE)
  path <- gsub("\\", "/", path, fixed = TRUE)
  path <- path.expand(path)
  if (grepl("^[A-Za-z][A-Za-z0-9+.-]*://", path)) {
    stop("Enter a file-system path instead of a web address.", call. = FALSE)
  }
  windows_drive <- grepl("^[A-Za-z]:", path)
  if (windows_drive && !grepl("^[A-Za-z]:/", path)) {
    stop("Enter a Windows absolute path in the form D:/folder.", call. = FALSE)
  }
  if (.Platform$OS.type != "windows" && (windows_drive || startsWith(path, "//"))) {
    stop("R is not currently running on Windows. Enter a path on the server where R is running.",
         call. = FALSE)
  }
  if (!grepl("^(/|[A-Za-z]:/)", path)) path <- file.path(base_dir, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

inspect_qc_paths <- function(script_value, results_value) {
  script <- resolve_user_path(script_value, "QC script file path")
  results <- resolve_user_path(results_value, "Results directory")
  if (tolower(tools::file_ext(script)) != "r") {
    stop("The QC script must be an .R file.", call. = FALSE)
  }
  if (!file.exists(script) || dir.exists(script)) {
    stop(paste("Cannot find the QC script:", script), call. = FALSE)
  }
  if (file.access(script, 4L) != 0L) {
    stop(paste("Cannot read the QC script:", script), call. = FALSE)
  }
  if (file.exists(results) && !dir.exists(results)) {
    stop("Specify a directory, not a file, for the results path.", call. = FALSE)
  }
  list(script_path = normalizePath(script, winslash = "/", mustWork = TRUE),
       results_root = results)
}

inspect_comparison_script <- function(script_value) {
  script <- resolve_user_path(script_value, "Comparison analysis script file path")
  if (tolower(tools::file_ext(script)) != "r") {
    stop("The comparison analysis script must be an .R file.", call. = FALSE)
  }
  if (!file.exists(script) || dir.exists(script)) {
    stop(paste("Cannot find the comparison analysis script:", script), call. = FALSE)
  }
  if (file.access(script, 4L) != 0L) {
    stop(paste("Cannot read the comparison analysis script:", script), call. = FALSE)
  }
  normalizePath(script, winslash = "/", mustWork = TRUE)
}

inspect_analysis_script <- function(script_value, label) {
  script <- resolve_user_path(script_value, paste(label, "script file path"))
  if (tolower(tools::file_ext(script)) != "r") {
    stop(paste(label, " script must be an .R file."), call. = FALSE)
  }
  if (!file.exists(script) || dir.exists(script)) {
    stop(paste(label, " script cannot be found:", script), call. = FALSE)
  }
  if (file.access(script, 4L) != 0L) {
    stop(paste(label, " script cannot be read:", script), call. = FALSE)
  }
  normalizePath(script, winslash = "/", mustWork = TRUE)
}

path_inside <- function(path, directory) {
  path <- sub("/+$", "", normalizePath(path, winslash = "/", mustWork = TRUE))
  directory <- sub("/+$", "", normalizePath(directory, winslash = "/", mustWork = TRUE))
  if (.Platform$OS.type == "windows") {
    path <- tolower(path); directory <- tolower(directory)
  }
  identical(path, directory) || startsWith(path, paste0(directory, "/"))
}

find_analysis_session_root <- function(path) {
  current <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (grepl("^Analysis_[0-9]{8}_[0-9]{6}", basename(current))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Cannot find the analysis session directory: ", path, call. = FALSE)
    }
    current <- parent
  }
}

analysis_stage_dir <- function(path, stage_name) {
  session_root <- find_analysis_session_root(path)
  stage_dir <- file.path(session_root, stage_name)
  if (!dir.exists(stage_dir) &&
      !dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create the analysis-stage directory: ", stage_dir, call. = FALSE)
  }
  normalizePath(stage_dir, winslash = "/", mustWork = TRUE)
}

create_timestamped_run_dir <- function(root, prefix, code_length = 4L,
                                       max_attempts = 100L) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S", tz = "Asia/Seoul")
  alphabet <- c(as.character(0:9), letters[1:6])
  for (attempt in seq_len(max_attempts)) {
    short_code <- paste(sample(alphabet, code_length, replace = TRUE),
                        collapse = "")
    run_dir <- file.path(root, paste(prefix, stamp, short_code, sep = "_"))
    if (dir.create(run_dir, showWarnings = FALSE)) {
      return(normalizePath(run_dir, winslash = "/", mustWork = TRUE))
    }
  }
  stop(prefix, " results directory cannot be created.", call. = FALSE)
}

create_qc_run_dir <- function(results_root) {
  # The first QC execution creates one persistent analysis-session hierarchy.
  if (!dir.exists(results_root)) {
    suppressWarnings(dir.create(results_root, recursive = TRUE, showWarnings = FALSE))
  }
  if (!dir.exists(results_root)) {
    stop(paste("Cannot create the results directory:", results_root), call. = FALSE)
  }
  root <- normalizePath(results_root, winslash = "/", mustWork = TRUE)
  # R removes its session temp directory when it exits; results must survive that.
  if (path_inside(root, tempdir())) {
    stop("Specify the results directory outside the R temporary directory.", call. = FALSE)
  }
  probe <- tempfile(".qc_write_test_", tmpdir = root)
  if (!isTRUE(suppressWarnings(file.create(probe)))) {
    stop(paste("Cannot write to the results directory:", root), call. = FALSE)
  }
  unlink(probe)
  session_dir <- create_timestamped_run_dir(root, "Analysis")
  stage_names <- c("01_QC_Results", "02_Comparison",
                   "03_Continuous", "04_Trajectory")
  stage_dirs <- file.path(session_dir, stage_names)
  created <- vapply(stage_dirs, dir.create, logical(1), showWarnings = FALSE)
  if (!all(created)) {
    stop("Cannot create the analysis-session stage directories.", call. = FALSE)
  }
  normalizePath(stage_dirs[[1]], winslash = "/", mustWork = TRUE)
}

create_comparison_run_dir <- function(qc_run_dir) {
  root <- analysis_stage_dir(qc_run_dir, "02_Comparison")
  create_timestamped_run_dir(root, "Comparison")
}

create_continuous_run_dir <- function(comparison_run_dir) {
  root <- analysis_stage_dir(comparison_run_dir, "03_Continuous")
  run_dir <- create_timestamped_run_dir(root, "Continuous")
  metadata_dir <- file.path(run_dir, "01_ROI_Metadata")
  zonation_dir <- file.path(run_dir, "02_Continuous_Zonation")
  if (!dir.create(metadata_dir) || !dir.create(zonation_dir)) {
    stop("Cannot create the Continuous analysis stage directories.", call. = FALSE)
  }
  list(
    root = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
    metadata_dir = normalizePath(metadata_dir, winslash = "/", mustWork = TRUE),
    zonation_dir = normalizePath(zonation_dir, winslash = "/", mustWork = TRUE)
  )
}

create_trajectory_run_dir <- function(comparison_run_dir) {
  stage_root <- analysis_stage_dir(comparison_run_dir, "04_Trajectory")
  trajectory_parent <- file.path(stage_root, "01_Trajectory")
  spatial_map_dir <- file.path(stage_root, "02_Spatial_Map")
  for (directory in c(trajectory_parent, spatial_map_dir)) {
    if (!dir.exists(directory) &&
        !dir.create(directory, recursive = TRUE, showWarnings = FALSE)) {
      stop("Cannot create the Trajectory analysis stage directories.",
           call. = FALSE)
    }
  }
  run_dir <- create_timestamped_run_dir(trajectory_parent, "Trajectory")
  beta_dir <- file.path(run_dir, "Beta")
  curve_dir <- file.path(run_dir, "Curve")
  if (!dir.create(beta_dir) || !dir.create(curve_dir)) {
    stop("Cannot create the Trajectory analysis stage directories.", call. = FALSE)
  }
  list(
    root = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
    trajectory_dir = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
    beta_dir = normalizePath(beta_dir, winslash = "/", mustWork = TRUE),
    curve_dir = normalizePath(curve_dir, winslash = "/", mustWork = TRUE),
    spatial_map_dir = normalizePath(spatial_map_dir, winslash = "/", mustWork = TRUE)
  )
}

create_spatial_map_run_dir <- function(analysis_path) {
  stage_root <- analysis_stage_dir(analysis_path, "04_Trajectory")
  spatial_parent <- file.path(stage_root, "02_Spatial_Map")
  if (!dir.exists(spatial_parent) &&
      !dir.create(spatial_parent, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create the parent directory for Spatial Map results.", call. = FALSE)
  }
  create_timestamped_run_dir(spatial_parent, "SpatialMap")
}

# Parse user-entered manual ROI exclusions.
# One dataset per line, for example:
#   GF1-4: ROI5, ROI29-ROI31, ROI63
#   SPF1-2: 3, 7, 10-12
parse_manual_exclusions <- function(value) {
  if (is.null(value) || !length(value) || all(is.na(value))) return(list())
  lines <- unlist(strsplit(paste(as.character(value), collapse = "\n"),
                           "\r?\n", perl = TRUE), use.names = FALSE)
  lines <- trimws(sub("#.*$", "", lines))
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return(list())

  numbers_by_dataset <- list()
  for (line_number in seq_along(lines)) {
    matched <- regmatches(
      lines[[line_number]],
      regexec("^\\s*([^:=]+?)\\s*[:=]\\s*(.*?)\\s*$", lines[[line_number]], perl = TRUE)
    )[[1]]
    if (length(matched) != 3L) {
      stop(sprintf(
        "The format of Manual exclude line %d is invalid. Example: GF1-4: ROI5, ROI29-ROI31",
        line_number
      ), call. = FALSE)
    }
    dataset_id <- toupper(gsub("_", "-", trimws(matched[[2]]), fixed = TRUE))
    if (!grepl("^(GF|SPF)[0-9]+(?:-[0-9]+)?$", dataset_id, perl = TRUE)) {
      stop(sprintf("The sample name on Manual exclude line %d is invalid: %s",
                   line_number, matched[[2]]), call. = FALSE)
    }
    roi_text <- trimws(matched[[3]])
    if (!nzchar(roi_text)) next
    tokens <- unlist(strsplit(gsub("[,;]", " ", roi_text), "\\s+", perl = TRUE),
                     use.names = FALSE)
    tokens <- tokens[nzchar(tokens)]
    roi_numbers <- integer()
    for (token in tokens) {
      cleaned <- gsub("ROI", "", toupper(token), fixed = TRUE)
      if (grepl("^[0-9]+$", cleaned)) {
        roi_numbers <- c(roi_numbers, as.integer(cleaned))
      } else if (grepl("^[0-9]+-[0-9]+$", cleaned)) {
        bounds <- as.integer(strsplit(cleaned, "-", fixed = TRUE)[[1]])
        if (bounds[[1]] > bounds[[2]]) {
          stop(sprintf("The start of the Manual exclude range is greater than the end: %s", token),
               call. = FALSE)
        }
        roi_numbers <- c(roi_numbers, seq.int(bounds[[1]], bounds[[2]]))
      } else {
        stop(sprintf("Invalid Manual exclude ROI format: %s", token),
             call. = FALSE)
      }
    }
    if (anyNA(roi_numbers) || any(roi_numbers < 1L | roi_numbers > 96L)) {
      stop(sprintf("ROI numbers for %s must be between 1 and 96.", dataset_id), call. = FALSE)
    }
    numbers_by_dataset[[dataset_id]] <- unique(c(
      numbers_by_dataset[[dataset_id]], roi_numbers
    ))
  }
  lapply(numbers_by_dataset, function(x) paste0("ROI", sort(unique(x))))
}

format_manual_exclusions <- function(value) {
  if (!length(value)) return("None")
  paste(vapply(names(value), function(id) {
    paste0(id, ": ", paste(value[[id]], collapse = ", "))
  }, character(1)), collapse = " / ")
}

format_count <- function(x) {
  if (is.null(x) || length(x) == 0 || is.na(x[[1]])) return("\u2014")
  format(as.numeric(x[[1]]), big.mark = ",", scientific = FALSE, trim = TRUE)
}

register_diann_source <- function(con, path, view_name, required = character()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  sql_path <- as.character(DBI::dbQuoteString(con, path))
  sql_view <- as.character(DBI::dbQuoteIdentifier(con, view_name))
  DBI::dbExecute(con, paste("CREATE OR REPLACE VIEW", sql_view,
                            "AS SELECT * FROM read_parquet(", sql_path, ")"))
  schema <- DBI::dbGetQuery(con, paste("DESCRIBE", sql_view))
  missing <- setdiff(required, schema$column_name)
  if (length(missing)) {
    stop(paste("Required columns are missing from the final report:", paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  schema$column_name
}

query_diann_summary <- function(con) {
  # Only one result row is returned to R, even for a multi-gigabyte input.
  DBI::dbGetQuery(con, paste0(
    "SELECT COUNT(*) AS n_rows, ",
    "COUNT(DISTINCT NULLIF(\"Run\", '')) AS n_runs, ",
    "COUNT(DISTINCT NULLIF(\"Protein.Group\", '')) AS n_groups ",
    "FROM diann_input"
  ))
}

query_diann_runs <- function(con) {
  DBI::dbGetQuery(con, paste0(
    "SELECT COALESCE(NULLIF(CAST(\"Run\" AS VARCHAR), ''), '(missing)') AS run_name, ",
    "COUNT(*) AS n_rows, ",
    "COUNT(DISTINCT NULLIF(\"Protein.Group\", '')) AS n_groups ",
    "FROM diann_input GROUP BY 1 ORDER BY 1 LIMIT ",
    as.integer(APP_CONFIG$max_run_summary_rows + 1L)
  ))
}

query_comparison_groups <- function(con) {
  run_names <- DBI::dbGetQuery(
    con,
    'SELECT DISTINCT CAST("Run" AS VARCHAR) AS run_name FROM diann_input WHERE "Run" IS NOT NULL'
  )$run_name
  run_names <- unique(as.character(run_names[!is.na(run_names)]))
  if (!length(run_names)) return(list(GF = character(), SPF = character()))
  base <- basename(gsub("\\", "/", run_names, fixed = TRUE))
  upper <- toupper(base)
  pattern <- "(?<![[:alnum:]])(?:SPF|GF)[0-9]+(?:[-_][0-9]+)?(?![[:alnum:]-])"
  tokens <- regmatches(upper, gregexpr(pattern, upper, perl = TRUE))
  ambiguous <- vapply(tokens, function(x) {
    length(unique(gsub("_", "-", x, fixed = TRUE))) > 1L
  }, logical(1))
  if (any(ambiguous)) {
    stop("Multiple GF/SPF groups were detected in a single Run name: ",
         paste(head(run_names[ambiguous], 5L), collapse = "; "), call. = FALSE)
  }
  ids <- sort(unique(unlist(lapply(tokens, function(x) {
    gsub("_", "-", x, fixed = TRUE)
  }), use.names = FALSE)))
  list(GF = ids[startsWith(ids, "GF")], SPF = ids[startsWith(ids, "SPF")])
}

build_explorer_sql <- function(con, columns, run_text = "", protein_text = "", limit = 100L) {
  # Identifiers and values are quoted separately; text filters are literal substrings.
  select_sql <- paste(as.character(DBI::dbQuoteIdentifier(con, columns)), collapse = ", ")
  clauses <- character()
  if (nzchar(trimws(run_text))) {
    clauses <- c(clauses, paste0(
      "contains(CAST(\"Run\" AS VARCHAR), ",
      as.character(DBI::dbQuoteString(con, trimws(run_text))), ")"
    ))
  }
  if (nzchar(trimws(protein_text))) {
    clauses <- c(clauses, paste0(
      "contains(CAST(\"Protein.Group\" AS VARCHAR), ",
      as.character(DBI::dbQuoteString(con, trimws(protein_text))), ")"
    ))
  }
  where_sql <- if (length(clauses)) paste("WHERE", paste(clauses, collapse = " AND ")) else ""
  sprintf("SELECT %s FROM diann_input %s LIMIT %d", select_sql, where_sql, as.integer(limit))
}


# QC_Code.R integration:
# Preferred interface:
#   run_qc <- function(first_pass_path, report_path, output_dir) {
#     # Use the two file paths and save outputs under output_dir.
#     # Return NULL when saving files yourself, or:
#     list(tables = list(summary = data.frame(...)),
#          plots = list(figure = ggplot_object),
#          files = c("relative/or/absolute/file.pdf"))
#   }
# Legacy top-level scripts receive first_pass_path, report_path, output_dir
# and qc_inputs in their source environment. They must write results to output_dir.
# QC_Code.R now supplies run_qc() with automatic GF/SPF IDs and dataset-scoped QC.
qc_worker <- function(first_pass_path, report_path, output_dir, script_path,
                       manual_exclude_by_dataset = list(),
                       db_memory_limit = "1GB", png_dpi = 600L) {
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  script_path <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
  first_pass_path <- normalizePath(first_pass_path, winslash = "/", mustWork = TRUE)
  report_path <- normalizePath(report_path, winslash = "/", mustWork = TRUE)
  on.exit(grDevices::graphics.off(), add = TRUE)
  env <- new.env(parent = globalenv())
  env$first_pass_path <- first_pass_path
  env$report_path <- report_path
  env$output_dir <- output_dir
  env$manual_exclude_by_dataset <- manual_exclude_by_dataset
  env$qc_inputs <- list(first_pass_path = first_pass_path,
                       report_path = report_path, output_dir = output_dir,
                       manual_exclude_by_dataset = manual_exclude_by_dataset)
  # Catch base graphics produced directly by the script. ggplot objects can also
  # be returned in plots; no large data frames are sent back to the Shiny process.
  grDevices::png(file.path(output_dir, "QC_base_%03d.png"),
                 width = 1600, height = 1000, res = 150)
  cat("QC started:", format(Sys.time()), "\n")
  cat("Input 1:", basename(first_pass_path), "\nInput 2:", basename(report_path), "\n")
  sourced <- source(script_path, local = env, chdir = FALSE, encoding = "UTF-8")
  if (exists("qc_results", envir = env, inherits = FALSE) && !is.null(env$qc_results)) {
    result <- env$qc_results
  } else if (exists("run_qc", envir = env, inherits = FALSE) &&
             is.function(env$run_qc)) {
    result <- env$run_qc(first_pass_path = first_pass_path,
                         report_path = report_path, output_dir = output_dir,
                         manual_exclude_by_dataset = manual_exclude_by_dataset,
                         db_memory_limit = db_memory_limit, png_dpi = png_dpi)
  } else {
    result <- sourced$value
  }
  if (inherits(result, "ggplot")) result <- list(plots = list(QC_plot = result))
  if (is.data.frame(result)) result <- list(tables = list(QC_table = result))
  if (!is.list(result)) result <- list()

  empty_excluded_rois <- data.frame(
    Dataset_ID = character(), ROI = character(), ROI_Number = integer(),
    Final_QC_Action = character(), Exclusion_Source = character(),
    Original_Run = character(), stringsAsFactors = FALSE
  )
  excluded_rois <- result$excluded_rois
  excluded_file <- file.path(output_dir, "QC_Excluded_ROIs.csv")
  if (!is.data.frame(excluded_rois) && file.exists(excluded_file)) {
    excluded_rois <- tryCatch(
      utils::read.csv(excluded_file, stringsAsFactors = FALSE,
                      check.names = FALSE),
      error = function(e) NULL
    )
  }
  if (!is.data.frame(excluded_rois)) excluded_rois <- empty_excluded_rois
  if (!nrow(excluded_rois) && length(manual_exclude_by_dataset)) {
    manual_rows <- lapply(names(manual_exclude_by_dataset), function(dataset_id) {
      rois <- unique(as.character(manual_exclude_by_dataset[[dataset_id]]))
      if (!length(rois)) return(NULL)
      data.frame(
        Dataset_ID = dataset_id,
        ROI = rois,
        ROI_Number = suppressWarnings(as.integer(sub(
          "^ROI0*([0-9]+)$", "\\1", toupper(rois), perl = TRUE
        ))),
        Final_QC_Action = "EXCLUDE_MANUAL",
        Exclusion_Source = "Manual",
        Original_Run = NA_character_,
        stringsAsFactors = FALSE
      )
    })
    manual_rows <- manual_rows[!vapply(manual_rows, is.null, logical(1))]
    if (length(manual_rows)) excluded_rois <- do.call(rbind, manual_rows)
  }
  if (nrow(excluded_rois)) {
    if (!"Dataset_ID" %in% names(excluded_rois)) {
      stop("The QC excluded-ROI result does not contain a Dataset_ID column.", call. = FALSE)
    }
    if (!"ROI_Number" %in% names(excluded_rois)) {
      if (!"ROI" %in% names(excluded_rois)) {
        stop("The QC excluded-ROI result does not contain an ROI or ROI_Number column.",
             call. = FALSE)
      }
      excluded_rois$ROI_Number <- suppressWarnings(as.integer(sub(
        "^ROI0*([0-9]+)$", "\\1",
        toupper(trimws(as.character(excluded_rois$ROI))), perl = TRUE
      )))
    }
    if (anyNA(excluded_rois$ROI_Number)) {
      stop("The QC excluded-ROI result contains ROI numbers that cannot be parsed.",
           call. = FALSE)
    }
    excluded_rois$Dataset_ID <- toupper(gsub(
      "_", "-", trimws(as.character(excluded_rois$Dataset_ID)), fixed = TRUE
    ))
    excluded_rois$ROI_Number <- as.integer(excluded_rois$ROI_Number)
    excluded_rois$ROI <- paste0("ROI", excluded_rois$ROI_Number)
    for (column in c("Final_QC_Action", "Exclusion_Source", "Original_Run")) {
      if (!column %in% names(excluded_rois)) excluded_rois[[column]] <- NA_character_
      excluded_rois[[column]] <- as.character(excluded_rois[[column]])
    }
    excluded_rois <- unique(excluded_rois[, c(
      "Dataset_ID", "ROI", "ROI_Number", "Final_QC_Action",
      "Exclusion_Source", "Original_Run"
    ), drop = FALSE])
  } else {
    excluded_rois <- empty_excluded_rois
  }
  utils::write.csv(excluded_rois, excluded_file, row.names = FALSE, na = "")

  safe_label <- function(x, fallback) {
    if (is.null(x) || is.na(x) || !nzchar(x)) return(fallback)
    x <- gsub("[^A-Za-z0-9_-]+", "_", x)
    if (!nzchar(gsub("_", "", x))) fallback else x
  }
  if (is.data.frame(result$tables)) result$tables <- list(QC_table = result$tables)
  if (is.list(result$tables)) {
    for (i in seq_along(result$tables)) {
      table <- result$tables[[i]]
      if (!is.data.frame(table)) next
      nm <- if (!is.null(names(result$tables))) names(result$tables)[[i]] else NULL
      filename <- sprintf("QC_table_%03d_%s.csv", i, safe_label(nm, "table"))
      utils::write.csv(table, file.path(output_dir, filename),
                       row.names = FALSE, na = "", fileEncoding = "UTF-8")
    }
  }
  if (inherits(result$plots, "ggplot") || inherits(result$plots, "recordedplot")) {
    result$plots <- list(QC_plot = result$plots)
  }
  if (is.list(result$plots)) {
    for (i in seq_along(result$plots)) {
      plot <- result$plots[[i]]
      nm <- if (!is.null(names(result$plots))) names(result$plots)[[i]] else NULL
      filename <- file.path(output_dir, sprintf(
        "QC_plot_%03d_%s.png", i, safe_label(nm, "plot")))
      if (inherits(plot, "ggplot")) {
        if (!requireNamespace("ggplot2", quietly = TRUE)) stop("The ggplot2 package is required.")
        ggplot2::ggsave(filename, plot = plot, width = 10, height = 6, dpi = 150)
      } else if (inherits(plot, "recordedplot")) {
        grDevices::png(filename, width = 1500, height = 900, res = 150)
        replayPlot(plot)
        grDevices::dev.off()
      } else {
        message("Unsupported plot object: ", i, ". Save it to output_dir explicitly.")
      }
    }
  }
  if (!is.null(result$summary)) {
    writeLines(capture.output(print(result$summary)),
               file.path(output_dir, "QC_summary.txt"), useBytes = TRUE)
  }
  if (is.character(result$files)) {
    export_dir <- file.path(output_dir, "exported")
    for (i in seq_along(result$files)) {
      f <- result$files[[i]]
      candidate <- if (file.exists(file.path(output_dir, f))) file.path(output_dir, f) else f
      if (!file.exists(candidate) || dir.exists(candidate)) next
      candidate <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
      if (!startsWith(candidate, paste0(output_dir, "/"))) {
        dir.create(export_dir, showWarnings = FALSE)
        destination <- file.path(export_dir, paste0(sprintf("%03d_", i), basename(candidate)))
        if (!file.copy(candidate, destination, overwrite = TRUE)) {
          stop("Failed to copy the QC result file: ", basename(candidate))
        }
      }
    }
  }
  grDevices::graphics.off()
  paths <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  if (length(paths)) paths <- paths[!dir.exists(paths)]
  if (length(paths)) paths <- paths[!is.na(file.size(paths)) & file.size(paths) > 0]
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  ext <- tolower(tools::file_ext(paths))
  kind <- ifelse(ext %in% c("png", "jpg", "jpeg"), "plot",
                 ifelse(ext %in% c("csv", "tsv"), "table", "file"))
  kind[basename(paths) == "QC_execution.log"] <- "log"
  if (!length(paths) || !any(kind != "log")) {
    stop(paste("The QC script finished, but no result files are available to display.",
               "Save PNG/CSV/PDF files in output_dir, or",
               "return tables/plots/files lists from run_qc()."))
  }
  cat("QC completed:", format(Sys.time()), "\n")
  list(files = data.frame(
    name = substring(paths, nchar(output_dir) + 2L), path = paths,
    kind = kind, size_bytes = file.size(paths), stringsAsFactors = FALSE
  ), dataset_summary = if (is.data.frame(result$summary)) result$summary else NULL,
     excluded_rois = excluded_rois)
}

comparison_worker <- function(report_path, output_dir, script_path,
                              gf_sets, spf_sets, excluded_rois = NULL,
                              db_memory_limit = "1GB", png_dpi = 600L) {
  report_path <- normalizePath(report_path, winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  script_path <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
  gf_sets <- unique(as.character(gf_sets))
  spf_sets <- unique(as.character(spf_sets))
  if (length(gf_sets) < 3L || length(spf_sets) < 3L) {
    stop("Comparison analysis requires at least three sets each for GF and SPF.", call. = FALSE)
  }
  on.exit(grDevices::graphics.off(), add = TRUE)
  env <- new.env(parent = globalenv())
  sourced <- source(script_path, local = env, chdir = FALSE, encoding = "UTF-8")
  if (!exists("run_comparison", envir = env, inherits = FALSE) ||
      !is.function(env$run_comparison)) {
    stop("Comparison_Pipeline.R does not contain run_comparison().", call. = FALSE)
  }
  if (!"excluded_rois" %in% names(formals(env$run_comparison))) {
    stop(paste0(
      "Comparison_Pipeline.R does not support QC excluded-ROI integration. ",
      "Replace it with the latest file included in the app package."
    ), call. = FALSE)
  }
  cat("Comparison started:", format(Sys.time()), "\n")
  cat("GF sets:", paste(gf_sets, collapse = ", "), "\n")
  cat("SPF sets:", paste(spf_sets, collapse = ", "), "\n")
  result <- env$run_comparison(
    report_path = report_path,
    output_dir = output_dir,
    gf_set_patterns = stats::setNames(gf_sets, gf_sets),
    spf_set_patterns = stats::setNames(spf_sets, spf_sets),
    excluded_rois = excluded_rois,
    db_memory_limit = db_memory_limit,
    png_dpi = png_dpi
  )
  grDevices::graphics.off()
  paths <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  if (length(paths)) paths <- paths[!dir.exists(paths)]
  if (length(paths)) paths <- paths[!is.na(file.size(paths)) & file.size(paths) > 0]
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  ext <- tolower(tools::file_ext(paths))
  kind <- ifelse(ext %in% c("png", "jpg", "jpeg"), "plot",
                 ifelse(ext %in% c("csv", "tsv"), "table", "file"))
  kind[basename(paths) == "Comparison_execution.log"] <- "log"
  if (!length(paths) || !any(kind != "log")) {
    stop("Comparison analysis finished, but no result files are available to display.", call. = FALSE)
  }
  cat("Comparison completed:", format(Sys.time()), "\n")
  list(
    files = data.frame(
      name = substring(paths, nchar(output_dir) + 2L),
      path = paths, kind = kind, size_bytes = file.size(paths),
      stringsAsFactors = FALSE
    ),
    summary = if (is.list(result) && is.data.frame(result$summary)) result$summary else NULL,
    gf_sets = gf_sets,
    spf_sets = spf_sets
  )
}

collect_analysis_files <- function(output_dir, log_name) {
  paths <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
  if (length(paths)) paths <- paths[!dir.exists(paths)]
  if (length(paths)) paths <- paths[!is.na(file.size(paths)) & file.size(paths) > 0]
  if (!length(paths)) {
    return(data.frame(name = character(), path = character(), kind = character(),
                      size_bytes = numeric(), stringsAsFactors = FALSE))
  }
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  ext <- tolower(tools::file_ext(paths))
  kind <- ifelse(ext %in% c("png", "jpg", "jpeg"), "plot",
                 ifelse(ext %in% c("csv", "tsv"), "table", "file"))
  kind[basename(paths) == log_name] <- "log"
  data.frame(
    name = substring(paths, nchar(output_dir) + 2L), path = paths,
    kind = kind, size_bytes = file.size(paths), stringsAsFactors = FALSE
  )
}

roi_metadata_worker <- function(input_file, output_dir, script_path) {
  collect_worker_files <- function(output_dir, log_name) {
    paths <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
    if (length(paths)) paths <- paths[!dir.exists(paths)]
    if (length(paths)) paths <- paths[!is.na(file.size(paths)) & file.size(paths) > 0]
    if (!length(paths)) {
      return(data.frame(name = character(), path = character(), kind = character(),
                        size_bytes = numeric(), stringsAsFactors = FALSE))
    }
    paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
    ext <- tolower(tools::file_ext(paths))
    kind <- ifelse(ext %in% c("png", "jpg", "jpeg"), "plot",
                   ifelse(ext %in% c("csv", "tsv"), "table", "file"))
    kind[basename(paths) == log_name] <- "log"
    data.frame(
      name = substring(paths, nchar(output_dir) + 2L), path = paths,
      kind = kind, size_bytes = file.size(paths), stringsAsFactors = FALSE
    )
  }
  input_file <- normalizePath(input_file, winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  script_path <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
  on.exit(grDevices::graphics.off(), add = TRUE)
  env <- new.env(parent = globalenv())
  source(script_path, local = env, chdir = FALSE, encoding = "UTF-8")
  if (!exists("run_roi_metadata", envir = env, inherits = FALSE) ||
      !is.function(env$run_roi_metadata)) {
    stop("ROI_Metadata_Generation.R does not contain run_roi_metadata().",
         call. = FALSE)
  }
  cat("ROI metadata generation started:", format(Sys.time()), "\n")
  cat("Coordinate input:", input_file, "\n")
  result <- env$run_roi_metadata(input_file = input_file, output_dir = output_dir)
  metadata_file <- if (is.list(result)) result$metadata_file else NULL
  if (is.null(metadata_file)) metadata_file <- file.path(output_dir, "ROI_metadata.csv")
  if (!file.exists(metadata_file)) {
    stop("ROI metadata generation finished, but ROI_metadata.csv cannot be found.",
         call. = FALSE)
  }
  metadata_file <- normalizePath(metadata_file, winslash = "/", mustWork = TRUE)
  grDevices::graphics.off()
  files <- collect_worker_files(output_dir, "ROI_Metadata_execution.log")
  cat("ROI metadata generation completed:", format(Sys.time()), "\n")
  list(
    files = files,
    summary = if (is.list(result) && is.data.frame(result$summary)) result$summary else NULL,
    metadata_file = metadata_file
  )
}

continuous_zonation_worker <- function(report_path, metadata_path, output_dir,
                                       script_path) {
  collect_worker_files <- function(output_dir, log_name) {
    paths <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
    if (length(paths)) paths <- paths[!dir.exists(paths)]
    if (length(paths)) paths <- paths[!is.na(file.size(paths)) & file.size(paths) > 0]
    if (!length(paths)) {
      return(data.frame(name = character(), path = character(), kind = character(),
                        size_bytes = numeric(), stringsAsFactors = FALSE))
    }
    paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
    ext <- tolower(tools::file_ext(paths))
    kind <- ifelse(ext %in% c("png", "jpg", "jpeg"), "plot",
                   ifelse(ext %in% c("csv", "tsv"), "table", "file"))
    kind[basename(paths) == log_name] <- "log"
    data.frame(
      name = substring(paths, nchar(output_dir) + 2L), path = paths,
      kind = kind, size_bytes = file.size(paths), stringsAsFactors = FALSE
    )
  }
  report_path <- normalizePath(report_path, winslash = "/", mustWork = TRUE)
  metadata_path <- normalizePath(metadata_path, winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  script_path <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
  on.exit(grDevices::graphics.off(), add = TRUE)
  env <- new.env(parent = globalenv())
  source(script_path, local = env, chdir = FALSE, encoding = "UTF-8")
  if (!exists("run_continuous_zonation", envir = env, inherits = FALSE) ||
      !is.function(env$run_continuous_zonation)) {
    stop(paste0(basename(script_path),
                " does not contain run_continuous_zonation()."), call. = FALSE)
  }
  cat("Continuous zonation analysis started:", format(Sys.time()), "\n")
  cat("report.parquet:", report_path, "\n")
  cat("ROI metadata:", metadata_path, "\n")
  result <- env$run_continuous_zonation(
    report_path = report_path,
    metadata_path = metadata_path,
    output_dir = output_dir
  )
  grDevices::graphics.off()
  files <- collect_worker_files(output_dir, "Continuous_Zonation_execution.log")
  if (!nrow(files) || !any(files$kind != "log")) {
    stop("Continuous zonation finished, but no result files are available to display.",
         call. = FALSE)
  }
  cat("Continuous zonation analysis completed:", format(Sys.time()), "\n")
  list(
    files = files,
    summary = if (is.list(result) && is.data.frame(result$summary)) result$summary else NULL,
    metadata_file = metadata_path
  )
}

parse_trajectory_targets <- function(value, max_targets = 50L) {
  if (is.null(value) || !length(value) || all(is.na(value))) return(character())
  targets <- unlist(strsplit(paste(as.character(value), collapse = "\n"),
                             "[,;[:space:]]+", perl = TRUE), use.names = FALSE)
  targets <- trimws(targets)
  targets <- targets[!is.na(targets) & nzchar(targets)]
  targets <- targets[!duplicated(toupper(targets))]
  if (length(targets) > max_targets) {
    stop("The maximum number of targets per run is ", max_targets, ".",
         call. = FALSE)
  }
  targets
}

trajectory_analysis_worker <- function(input_dir, beta_output_dir,
                                       curve_output_dir, beta_script_path,
                                       curve_script_path, targets,
                                       cache_dir = NULL) {
  collect_output_files <- function(root) {
    paths <- list.files(root, recursive = TRUE, full.names = TRUE)
    if (length(paths)) paths <- paths[!dir.exists(paths)]
    if (length(paths)) paths <- paths[!is.na(file.size(paths)) & file.size(paths) > 0]
    if (!length(paths)) {
      return(data.frame(name = character(), path = character(), kind = character(),
                        size_bytes = numeric(), stringsAsFactors = FALSE))
    }
    paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
    ext <- tolower(tools::file_ext(paths))
    data.frame(
      name = substring(paths, nchar(root) + 2L),
      path = paths,
      kind = ifelse(ext == "png", "plot",
                    ifelse(ext %in% c("csv", "tsv"), "table", "file")),
      size_bytes = file.size(paths), stringsAsFactors = FALSE
    )
  }
  targets <- trimws(as.character(targets))
  targets <- targets[!is.na(targets) & nzchar(targets)]
  targets <- targets[!duplicated(toupper(targets))]
  if (!length(targets)) {
    stop("Enter at least one UniProt accession or gene name.", call. = FALSE)
  }
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  beta_output_dir <- normalizePath(beta_output_dir, winslash = "/", mustWork = TRUE)
  curve_output_dir <- normalizePath(curve_output_dir, winslash = "/", mustWork = TRUE)
  beta_script_path <- normalizePath(beta_script_path, winslash = "/", mustWork = TRUE)
  curve_script_path <- normalizePath(curve_script_path, winslash = "/", mustWork = TRUE)
  cache_dir <- if (is.null(cache_dir) || !nzchar(trimws(as.character(cache_dir)))) {
    NULL
  } else {
    normalizePath(cache_dir, winslash = "/", mustWork = FALSE)
  }
  on.exit(grDevices::graphics.off(), add = TRUE)

  cache_index <- NULL
  if (!is.null(cache_dir)) {
    cache_index_file <- file.path(cache_dir, "Cache_Index.csv")
    if (file.exists(cache_index_file)) {
      cache_index <- tryCatch(
        utils::read.csv(cache_index_file, check.names = FALSE,
                        stringsAsFactors = FALSE, fileEncoding = "UTF-8"),
        error = function(e) NULL
      )
      required_cache <- c("Protein.Group", "Protein", "Beta_Status",
                          "Trajectory_Status", "Beta_PNG", "Trajectory_PNG")
      if (!is.null(cache_index) &&
          length(setdiff(required_cache, names(cache_index)))) cache_index <- NULL
    }
  }
  cache_path <- function(value) {
    value <- trimws(as.character(value))
    if (is.na(value) || !nzchar(value)) return(NA_character_)
    value <- gsub("\\", "/", value, fixed = TRUE)
    if (grepl("^(/|[A-Za-z]:/)", value)) value else file.path(cache_dir, value)
  }
  cached_target <- function(target) {
    if (is.null(cache_index) || !nrow(cache_index)) return(NULL)
    wanted <- toupper(trimws(target))
    id_values <- ifelse(is.na(cache_index$Protein.Group), "",
                        as.character(cache_index$Protein.Group))
    id_hit <- toupper(trimws(id_values)) == wanted |
      vapply(strsplit(id_values, "[;,|]", perl = TRUE), function(values) {
        any(toupper(trimws(values)) == wanted)
      }, logical(1))
    gene_hit <- vapply(strsplit(ifelse(is.na(cache_index$Protein), "",
                                      as.character(cache_index$Protein)),
                                "[;,|]", perl = TRUE), function(values) {
      any(toupper(trimws(values)) == wanted)
    }, logical(1))
    rows <- cache_index[id_hit | gene_hit, , drop = FALSE]
    if (nrow(rows) != 1L) return(NULL)
    rows
  }

  beta_env <- new.env(parent = globalenv())
  source(beta_script_path, local = beta_env, chdir = FALSE, encoding = "UTF-8")
  if (!exists("run_zonation_coefficient_beta", envir = beta_env, inherits = FALSE) ||
      !is.function(beta_env$run_zonation_coefficient_beta)) {
    stop(basename(beta_script_path),
         " does not contain run_zonation_coefficient_beta().", call. = FALSE)
  }
  curve_env <- new.env(parent = globalenv())
  source(curve_script_path, local = curve_env, chdir = FALSE, encoding = "UTF-8")
  if (!exists("run_selected_trajectory", envir = curve_env, inherits = FALSE) ||
      !is.function(curve_env$run_selected_trajectory)) {
    stop(basename(curve_script_path),
         " does not contain run_selected_trajectory().", call. = FALSE)
  }

  cat("Trajectory analysis started:", format(Sys.time()), "\n")
  cat("Targets (", length(targets), "): ", paste(targets, collapse = ", "), "\n",
      sep = "")
  cat("Continuous output:", input_dir, "\n")
  beta_plots <- list()
  curve_plots <- list()
  summary_rows <- vector("list", length(targets))
  for (i in seq_along(targets)) {
    target <- targets[[i]]
    cat("\n[", i, "/", length(targets), "] Target: ", target, "\n", sep = "")
    cache_row <- cached_target(target)
    beta_cache_files <- if (!is.null(cache_row)) {
      cache_path(cache_row$Beta_PNG[[1]])
    } else character()
    curve_cache_files <- if (!is.null(cache_row)) {
      cache_path(cache_row$Trajectory_PNG[[1]])
    } else character()
    beta_cached <- length(beta_cache_files) == 1L &&
      identical(cache_row$Beta_Status[[1]], "Created") &&
      all(file.exists(beta_cache_files))
    curve_cached <- length(curve_cache_files) == 1L &&
      identical(cache_row$Trajectory_Status[[1]], "Created") &&
      all(file.exists(curve_cache_files))

    beta_result <- if (beta_cached) {
      cat("  Beta: cache hit\n")
      list(png_file = beta_cache_files[[1]], pdf_file = NA_character_,
           protein_group = cache_row$Protein.Group[[1]])
    } else tryCatch(
      beta_env$run_zonation_coefficient_beta(
        input_dir = input_dir, output_dir = beta_output_dir, target = target
      ), error = function(e) e
    )
    curve_result <- if (curve_cached) {
      cat("  Trajectory: cache hit\n")
      list(png_file = curve_cache_files[[1]], pdf_file = NA_character_,
           protein_group = cache_row$Protein.Group[[1]],
           protein = cache_row$Protein[[1]])
    } else tryCatch(
      curve_env$run_selected_trajectory(
        input_dir = input_dir, output_dir = curve_output_dir, target = target
      ), error = function(e) e
    )

    beta_ok <- !inherits(beta_result, "error") &&
      file.exists(beta_result$png_file) && file.size(beta_result$png_file) > 0
    curve_ok <- !inherits(curve_result, "error") &&
      file.exists(curve_result$png_file) && file.size(curve_result$png_file) > 0
    if (beta_ok) {
      beta_plots[[length(beta_plots) + 1L]] <- data.frame(
        Target = target, PNG = beta_result$png_file, PDF = beta_result$pdf_file,
        stringsAsFactors = FALSE
      )
    }
    if (curve_ok) {
      curve_plots[[length(curve_plots) + 1L]] <- data.frame(
        Target = target, PNG = curve_result$png_file, PDF = curve_result$pdf_file,
        stringsAsFactors = FALSE
      )
    }
    messages <- character()
    if (inherits(beta_result, "error")) {
      messages <- c(messages, paste0("Beta: ", conditionMessage(beta_result)))
    }
    if (inherits(curve_result, "error")) {
      messages <- c(messages, paste0("Trajectory: ", conditionMessage(curve_result)))
    }
    summary_rows[[i]] <- data.frame(
      Target = target,
      Beta_Status = if (beta_ok) if (beta_cached) "Cached" else "Created" else "Failed",
      Trajectory_Status = if (curve_ok) if (curve_cached) "Cached" else "Created" else "Failed",
      Protein = if (curve_ok) curve_result$protein else NA_character_,
      Protein.Group = if (curve_ok) curve_result$protein_group else
        if (beta_ok) beta_result$protein_group else NA_character_,
      Message = paste(c(
        if (beta_cached) "Beta cache reused" else character(),
        if (curve_cached) "Trajectory cache reused" else character(),
        messages
      ), collapse = " | "),
      stringsAsFactors = FALSE
    )
  }
  grDevices::graphics.off()

  beta_plots <- if (length(beta_plots)) do.call(rbind, beta_plots) else
    data.frame(Target = character(), PNG = character(), PDF = character())
  curve_plots <- if (length(curve_plots)) do.call(rbind, curve_plots) else
    data.frame(Target = character(), PNG = character(), PDF = character())
  summary <- do.call(rbind, summary_rows)
  if (!nrow(beta_plots) && !nrow(curve_plots)) {
    stop("No plots were generated for the entered targets. Check the gene names or UniProt accessions.\n",
         paste(summary$Message, collapse = "\n"), call. = FALSE)
  }
  run_root <- dirname(beta_output_dir)
  files <- collect_output_files(run_root)
  cat("Trajectory analysis completed:", format(Sys.time()), "\n")
  list(
    files = files,
    targets = targets,
    beta_plots = beta_plots,
    curve_plots = curve_plots,
    summary = summary
  )
}

spatial_map_worker <- function(report_path, output_root, script_path, targets,
                               coordinate_path, roi_metadata_path,
                               analysis_ready_path, histology_root,
                               require_all_sets = TRUE) {
  required_packages <- c("tidyverse", "arrow", "magick")
  missing_packages <- required_packages[!vapply(
    required_packages, requireNamespace, logical(1), quietly = TRUE
  )]
  if (length(missing_packages)) {
    stop("Required R packages for Spatial Map are missing: ",
         paste(missing_packages, collapse = ", "), call. = FALSE)
  }
  report_path <- normalizePath(report_path, winslash = "/", mustWork = TRUE)
  output_root <- normalizePath(output_root, winslash = "/", mustWork = TRUE)
  script_path <- normalizePath(script_path, winslash = "/", mustWork = TRUE)
  coordinate_path <- normalizePath(coordinate_path, winslash = "/", mustWork = TRUE)
  roi_metadata_path <- normalizePath(roi_metadata_path, winslash = "/",
                                     mustWork = TRUE)
  analysis_ready_path <- if (!is.null(analysis_ready_path) &&
                             length(analysis_ready_path) == 1L &&
                             !is.na(analysis_ready_path) &&
                             nzchar(analysis_ready_path) &&
                             file.exists(analysis_ready_path)) {
    normalizePath(analysis_ready_path, winslash = "/", mustWork = TRUE)
  } else ""
  histology_root <- normalizePath(histology_root, winslash = "/", mustWork = TRUE)
  script_text <- readLines(script_path, warn = FALSE, encoding = "UTF-8")
  if (!any(grepl("SPATIAL_EXPRESSION_SHINY_API", script_text, fixed = TRUE)) ||
      !any(grepl("SPATIAL_EXPRESSION_MULTI_TARGET_API", script_text, fixed = TRUE)) ||
      !any(grepl("SPATIAL_MAP_TRAJECTORY_PANEL_API", script_text, fixed = TRUE)) ||
      !any(grepl("SPATIAL_MAP_SAMPLE_LABEL_API", script_text, fixed = TRUE)) ||
      !any(grepl("2026-11-17-v7-trajectory-consistency", script_text, fixed = TRUE)) ||
      !any(grepl("SPATIAL_MAP_COMBINED_SAMPLE_HEADER_API", script_text, fixed = TRUE)) ||
      !any(grepl("SPATIAL_MAP_LEGEND_LAYOUT_API", script_text, fixed = TRUE)) ||
      !any(grepl("SPATIAL_MAP_SHINY_SPLIT_PANEL_API", script_text, fixed = TRUE)) ||
      !any(grepl("SPATIAL_MAP_TRAJECTORY_DATA_API", script_text, fixed = TRUE)) ||
      !any(grepl("SPATIAL_MAP_COORDINATE_FILE", script_text, fixed = TRUE))) {
    stop(
      "The current Spatial Expression script is not the multi-target Shiny-compatible version. ",
      "Use the latest codes/Spatial_Expression_Shiny.R file.",
      call. = FALSE
    )
  }
  targets <- trimws(as.character(targets))
  targets <- targets[!is.na(targets) & nzchar(targets)]
  targets <- targets[!duplicated(toupper(targets))]
  if (!length(targets)) {
    stop("Enter at least one UniProt accession or gene name.", call. = FALSE)
  }
  if (length(targets) > 10L) {
    stop("Spatial Map can process up to 10 targets at a time.",
         call. = FALSE)
  }
  safe_targets <- gsub("[^A-Za-z0-9_-]", "_", targets)
  if (anyDuplicated(tolower(safe_targets))) {
    stop("Some targets resolve to the same result-directory name: ",
         paste(targets[duplicated(tolower(safe_targets))], collapse = ", "),
         call. = FALSE)
  }
  env_values <- c(
    SPATIAL_MAP_INPUT_FILE = report_path,
    SPATIAL_MAP_OUTPUT_ROOT = output_root,
    SPATIAL_MAP_TARGET = targets[[1]],
    SPATIAL_MAP_COORDINATE_FILE = coordinate_path,
    SPATIAL_MAP_ROI_METADATA_FILE = roi_metadata_path,
    SPATIAL_MAP_HISTOLOGY_ROOT = histology_root,
    SPATIAL_MAP_REQUIRE_ALL_SETS = if (isTRUE(require_all_sets)) "true" else "false"
  )
  if (nzchar(analysis_ready_path)) {
    env_values <- c(
      env_values,
      SPATIAL_MAP_ANALYSIS_READY_FILE = analysis_ready_path
    )
  }
  old_values <- Sys.getenv(names(env_values), unset = NA_character_)
  on.exit({
    for (name in names(env_values)) {
      old <- old_values[[name]]
      if (is.na(old)) {
        Sys.unsetenv(name)
      } else {
        do.call(Sys.setenv, stats::setNames(list(old), name))
      }
    }
    grDevices::graphics.off()
    if (exists(".spatial_map_report_cache", envir = globalenv(), inherits = FALSE)) {
      rm(".spatial_map_report_cache", envir = globalenv())
    }
  }, add = TRUE)
  do.call(Sys.setenv, as.list(env_values))
  cat("Spatial Map started:", format(Sys.time()), "\n")
  cat("Targets:", paste(targets, collapse = ", "), "\n")
  cat("report.parquet:", report_path, "\n")
  cat("Coordinate CSV:", coordinate_path, "\n")
  cat("ROI metadata:", roi_metadata_path, "\n")
  cat("Trajectory analysis-ready data:",
      if (nzchar(analysis_ready_path)) analysis_ready_path else
        "not available; rebuild from ROI metadata and report.parquet",
      "\n")
  cat("Histology folder:", histology_root, "\n")

  target_summary <- lapply(targets, function(target) {
    safe_target <- gsub("[^A-Za-z0-9_-]", "_", target)
    target_dir <- file.path(output_root, safe_target)
    cat("\n--- Spatial Map target:", target, "---\n")
    Sys.setenv(SPATIAL_MAP_TARGET = target)
    started <- Sys.time()
    error_message <- NULL
    tryCatch(
      source(script_path, local = new.env(parent = globalenv()),
             chdir = FALSE, encoding = "UTF-8"),
      error = function(e) error_message <<- conditionMessage(e),
      finally = grDevices::graphics.off()
    )
    overlay_count <- if (dir.exists(target_dir)) {
      length(list.files(target_dir, pattern = "_Spatial_MAP[.]png$",
                        full.names = TRUE))
    } else 0L
    status <- if (is.null(error_message) && overlay_count > 0L) {
      "Success"
    } else {
      "Failed"
    }
    if (is.null(error_message) && overlay_count == 0L) {
      error_message <- "The final Spatial MAP PNG was not generated."
    }
    data.frame(
      Target = target,
      Safe_Target = safe_target,
      Status = status,
      Overlay_Images = overlay_count,
      Elapsed_Seconds = round(as.numeric(difftime(Sys.time(), started,
                                                   units = "secs")), 2),
      Message = if (is.null(error_message)) "" else error_message,
      stringsAsFactors = FALSE
    )
  })
  target_summary <- do.call(rbind, target_summary)
  utils::write.csv(target_summary,
                   file.path(output_root, "Spatial_Map_Target_Summary.csv"),
                   row.names = FALSE)
  successful_targets <- target_summary$Target[target_summary$Status == "Success"]
  if (!length(successful_targets)) {
    stop("Failed to generate Spatial Maps for the entered targets.\n",
         paste(target_summary$Message, collapse = "\n"), call. = FALSE)
  }
  paths <- list.files(output_root, recursive = TRUE, full.names = TRUE)
  if (length(paths)) paths <- paths[!dir.exists(paths)]
  if (length(paths)) paths <- paths[!is.na(file.size(paths)) & file.size(paths) > 0]
  if (!length(paths)) {
    stop("The Spatial Map script finished, but no result files were generated.", call. = FALSE)
  }
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  ext <- tolower(tools::file_ext(paths))
  files <- data.frame(
    name = substring(paths, nchar(output_root) + 2L),
    path = paths,
    kind = ifelse(ext %in% c("png", "jpg", "jpeg"), "plot",
                  ifelse(ext %in% c("csv", "tsv"), "table", "file")),
    size_bytes = file.size(paths), stringsAsFactors = FALSE
  )
  first_path_part <- sub("/.*$", "", files$name)
  target_index <- match(first_path_part, target_summary$Safe_Target)
  files$target <- target_summary$Target[target_index]
  image_basename <- basename(files$name)
  combined_image <- grepl(
    "_(GF|SPF)_Trajectory_SpatialMap_Combined[.]png$",
    image_basename,
    ignore.case = TRUE
  )
  combined_group <- sub(
    ".*_(GF|SPF)_Trajectory_SpatialMap_Combined[.]png$",
    "\\1-Combined",
    image_basename,
    ignore.case = TRUE
  )
  trajectory_panel_image <- grepl(
    "_(GF|SPF)_Trajectory_Combined[.]png$",
    image_basename,
    ignore.case = TRUE
  )
  trajectory_panel_group <- sub(
    ".*_(GF|SPF)_Trajectory_Combined[.]png$",
    "\\1-Combined",
    image_basename,
    ignore.case = TRUE
  )
  files$view <- ifelse(
    trajectory_panel_image,
    "trajectory_combined",
    ifelse(
      combined_image,
      "legacy_combined",
      ifelse(grepl("_Spatial_MAP[.]png$", image_basename),
             "individual", "supporting")
    )
  )
  files$sample <- ifelse(
    files$kind == "plot",
    ifelse(
      trajectory_panel_image,
      trajectory_panel_group,
      ifelse(combined_image, combined_group,
             sub("_.*$", "", image_basename))
    ),
    NA_character_
  )
  images <- files[files$kind == "plot", , drop = FALSE]
  if (!nrow(images)) {
    stop("No Spatial Map PNG results were generated.", call. = FALSE)
  }
  cat("Spatial Map completed:", format(Sys.time()), "\n")
  list(files = files, images = images, targets = targets,
       successful_targets = successful_targets, summary = target_summary)
}

page_heading <- function(kicker, title, description) {
  div(class = "page-heading",
      div(class = "eyebrow", kicker),
      h2(title), p(description))
}

section_box <- function(title, ..., width = 12) {
  shinydashboard::box(title = title, ..., width = width,
                      solidHeader = FALSE, status = "primary")
}

metric_cards <- function(info, counts) {
  fluidRow(
    valueBox(
      if (is.null(info)) "\u2014" else sprintf("%.2f GiB", info$size / 1024^3),
      "Input file size", icon = icon("file"), color = "teal", width = 3
    ),
    valueBox(format_count(if (!is.null(counts)) counts$n_runs), "Runs",
             icon = icon("layer-group"), color = "aqua", width = 3),
    valueBox(format_count(if (!is.null(counts)) counts$n_groups), "Protein groups",
             icon = icon("cubes"), color = "navy", width = 3),
    valueBox(format_count(if (!is.null(counts)) counts$n_rows), "Report rows",
             icon = icon("table"), color = "light-blue", width = 3)
  )
}

logo_path <- file.path(APP_ROOT, "logo", APP_CONFIG$logo_file)
logo_exists <- file.exists(logo_path) && !dir.exists(logo_path)
brand <- if (logo_exists) {
  # Embed the small branding image directly in the page. This avoids a broken
  # image when app.R is sourced from a different working directory or when a
  # reverse proxy changes the app's URL prefix.
  logo_src <- base64enc::dataURI(file = logo_path, mime = "image/png")
  tags$img(src = logo_src, alt = APP_CONFIG$company_name, class = "company-logo")
} else {
  div(class = "brand-fallback", span(class = "brand-mark", "SP"),
      span(class = "brand-name", APP_CONFIG$company_name))
}

css <- "
body, .main-header .logo, h1, h2, h3, h4, .small-box {
  font-family: 'Segoe UI', 'Noto Sans KR', Arial, sans-serif;
}
html body.skin-blue .main-header .logo,
html body.skin-blue .main-header .logo:hover,
html body.skin-blue .main-header .navbar {
  background-color: @DEEP@ !important;
}
html body.skin-blue .main-header .navbar .sidebar-toggle {
  color: #FFFFFF !important;
}
html body.skin-blue .main-header .navbar .sidebar-toggle:hover {
  background-color: @PRIMARY@ !important;
}
.main-header .logo { padding: 0 18px; text-align: center; }
.company-logo { max-width: 210px; max-height: 34px; width: auto; height: auto;
  object-fit: contain; display: inline-block; vertical-align: middle; }
.brand-fallback { display: flex; align-items: center; gap: 11px; height: 50px; }
.brand-mark { border: 1px solid rgba(255,255,255,.65); border-radius: 7px;
  width: 31px; height: 31px; line-height: 29px; text-align: center; font-size: 12px;
  font-weight: 700; letter-spacing: .5px; flex: 0 0 auto; }
.brand-name { font-size: 15px; font-weight: 700; overflow: hidden; text-overflow: ellipsis; }
.header-caption { position: absolute !important; left: 55px; top: 0; color: #fff;
  line-height: 50px; font-size: 18px; font-weight: 500; pointer-events: none; }
.header-tag { color: #D9F0FC; padding: 17px 21px; font-size: 10px;
  line-height: 16px; letter-spacing: 1.4px; font-weight: 700; }
html body.skin-blue .wrapper,
html body.skin-blue .main-sidebar,
html body.skin-blue .left-side {
  background-color: @SIDEBAR@ !important;
}
html body.skin-blue .sidebar-menu > li.header { color: #D9F0FC; background: transparent;
  font-size: 10px; letter-spacing: 1.8px; padding: 27px 20px 12px; }
html body.skin-blue .sidebar-menu > li > a { color: #E8F5FC !important; padding: 15px 18px;
  border-left: 3px solid transparent; font-weight: 500;
  display: flex; align-items: center; }
html body.skin-blue .sidebar-menu > li > a > i {
  display: inline-flex !important;
  align-items: center;
  justify-content: center;
  width: 26px !important;
  min-width: 26px;
  margin-right: 10px !important;
  font-size: 15px;
  color: @LIGHT@ !important;
  text-align: center;
  box-sizing: border-box;
}
html body.skin-blue .sidebar-menu > li.active > a,
html body.skin-blue .sidebar-menu > li:hover > a {
  color: #FFFFFF !important;
  background-color: @DEEP@ !important;
  border-left-color: @LIGHT@ !important;
}
html body.skin-blue .sidebar-menu > li.active > a {
  background-color: @PRIMARY@ !important;
}
.sidebar-note { color: #D9F0FC; padding: 27px 22px; font-size: 11px; line-height: 1.8; }
.sidebar-note strong { color: #FFFFFF; display: block; font-size: 12px; }
.sidebar-file { border-top: 1px solid rgba(255,255,255,.08); margin: 24px 18px 0;
  padding: 19px 4px 0; color: #D9F0FC; font-size: 11px; }
.sidebar-file-name { display: block; margin-top: 8px; overflow: hidden;
  text-overflow: ellipsis; white-space: nowrap; color: #f0f5f8; }
.content-wrapper, .right-side { background-color: @BG@ !important; }
.content { padding: 25px 30px 35px; max-width: 1900px; margin: auto; }
.page-heading { margin: 3px 0 23px; }
.eyebrow { color: @PRIMARY@; font-size: 10px; font-weight: 700; letter-spacing: 2px; }
.page-heading h2 { color: #243740; font-size: 29px; font-weight: 600; margin: 8px 0; }
.page-heading p { color: #74848d; font-size: 13px; line-height: 1.7; margin: 0; }
.small-box { border-radius: 8px; box-shadow: 0 3px 10px rgba(32,46,51,.06); }
.small-box .inner { padding: 19px 20px; }
.small-box h3 { font-size: 29px; font-weight: 650; margin: 0 0 7px; }
.small-box p { font-size: 12px; opacity: .9; }
.small-box .icon { top: 11px; right: 17px; font-size: 56px; opacity: .65; }
.bg-teal { background-color: @DEEP@ !important; }
.bg-aqua { background-color: @ACCENT@ !important; }
.bg-navy { background-color: @SIDEBAR@ !important; }
.bg-light-blue { background-color: @LIGHT@ !important; color: @DEEP@ !important; }
.bg-light-blue a, .bg-light-blue .icon { color: @DEEP@ !important; }
.box { border-radius: 8px; border-top-width: 0; box-shadow: 0 2px 9px rgba(32,46,51,.055);
  margin-bottom: 22px; }
.box.box-primary { border-top-color: @PRIMARY@ !important; }
.box-header { padding: 20px 21px 15px; border-bottom: 1px solid #eef2f5; }
.box-header .box-title { font-size: 15px; font-weight: 600; color: #2b414c; }
.box-body { padding: 21px; }
.btn { border-radius: 5px; font-size: 12px; padding: 9px 14px; }
/* Display all primary run/query buttons in the same style as the Run QC button */
.btn-primary {
  background-color: @PRIMARY@ !important;
  border-color: @PRIMARY@ !important;
  color: #FFFFFF !important;
}

/* On hover or selection */
.btn-primary:hover,
.btn-primary:focus {
  background-color: @DEEP@ !important;
  border-color: @DEEP@ !important;
  color: #FFFFFF !important;
}

/* While the button is pressed */
.btn-primary:active,
.btn-primary.active,
.open > .dropdown-toggle.btn-primary {
  background-color: @SIDEBAR@ !important;
  border-color: @SIDEBAR@ !important;
  color: #FFFFFF !important;
}


.btn-default { border-color: #d8e1e6; color: #3d5662; }
.form-control { border-radius: 5px; border-color: #dce4e9; box-shadow: none; }
.form-control:focus {
  border-color: @ACCENT@ !important;
  box-shadow: 0 0 0 1px rgba(17,145,208,.16) !important;
}
label { font-size: 12px; color: #415964; font-weight: 600; }
.upload-zone { border: 1.5px dashed @LIGHT@; background: #F3FAFE; padding: 23px 18px 14px;
  border-radius: 7px; margin-bottom: 17px; }
.upload-glyph { color: @PRIMARY@; font-size: 28px; margin-bottom: 12px; }
.upload-zone .shiny-input-container { width: 100% !important; margin-bottom: 8px; }
.upload-zone .progress { margin-bottom: 0; height: 6px; border-radius: 10px; }
.upload-zone .progress-bar { background-color: @ACCENT@ !important; }
.muted { font-size: 12px; color: #7b8f9a; line-height: 1.8; }
.micro-label { color: #8c9ca5; text-transform: uppercase; font-size: 10px;
  letter-spacing: 1.4px; font-weight: 700; }
.status-banner { display: flex; gap: 13px; align-items: flex-start; border-radius: 6px;
  padding: 16px; background: #EAF6FD; color: @DEEP@; margin-bottom: 18px; }
.status-banner .fa, .status-banner .fas { margin-top: 2px; font-size: 17px; }
.status-ready { background: #EAF6FD; color: @DEEP@; }
.status-error { background: #fff0ef; color: #a64540; }
.status-banner strong { display: block; margin-bottom: 4px; font-size: 13px; }
.file-meta { display: grid; grid-template-columns: 110px 1fr; gap: 12px 15px;
  font-size: 12px; margin: 17px 0; }
.file-meta dt { color: #8a9aa3; font-weight: 400; }
.file-meta dd { margin: 0; color: #334d58; overflow-wrap: anywhere; }
.checklist { display: flex; flex-wrap: wrap; gap: 8px; margin-top: 18px; }
.check-pill { background: #EAF6FD; color: @DEEP@; padding: 6px 10px;
  border-radius: 4px; font-size: 11px; }
.empty-state { text-align: center; padding: 42px 20px; color: #8ba0aa; }
.empty-state i { display: block; font-size: 32px; color: #b9cbd2; margin-bottom: 15px; }
.empty-state strong { display: block; font-size: 14px; color: #526d7a; margin-bottom: 8px; }
.spatial-map-selector-actions { display: flex; flex-wrap: wrap; gap: 7px;
  margin: -2px 0 18px; }
.spatial-map-selector-actions .btn { margin: 0; }
.spatial-map-grid { display: grid;
  grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
  gap: 18px; margin: 4px 0 20px; }
.spatial-map-grid.combined-view { grid-template-columns: 1fr; }
.spatial-map-card { min-width: 0; overflow: hidden; background: #fff;
  border: 1px solid #e3eaee; border-radius: 7px;
  box-shadow: 0 2px 7px rgba(32,46,51,.055); }
.spatial-map-card-title { padding: 11px 13px; border-bottom: 1px solid #edf2f5;
  color: #2e4d59; font-size: 12px; font-weight: 700; }
.spatial-map-card .shiny-image-output { width: 100%; min-height: 180px; }
.spatial-map-card img { display: block; width: 100%; height: auto; }
.spatial-combined-stack { display: flex; flex-direction: column; gap: 22px;
  margin: 4px 0 20px; }
.spatial-combined-section { background: #fff; border: 1px solid #e3eaee;
  border-radius: 7px; overflow: hidden;
  box-shadow: 0 2px 7px rgba(32,46,51,.055); }
.spatial-combined-title { padding: 12px 16px; border-bottom: 1px solid #edf2f5;
  color: #2e4d59; font-size: 13px; font-weight: 700; }
.spatial-combined-row { display: grid;
  grid-template-columns: minmax(0,1fr) minmax(0,1fr);
  gap: 14px; padding: 14px; align-items: start; }
.spatial-combined-trajectory { min-width: 0; overflow: hidden; }
.spatial-combined-trajectory img { display: block; width: 100%; height: auto; }
.spatial-combined-histology { display: grid;
  grid-template-columns: repeat(2, minmax(0,1fr)); gap: 12px; min-width: 0; }
.spatial-combined-histology .spatial-map-card-title { font-size: 13px;
  font-weight: 600; background: #fff; }
.preview-note { margin: 0 0 15px; }
.table { font-size: 12px; }
table.dataTable thead th { background: #EAF6FD; color: @DEEP@; font-weight: 600;
  border-bottom: 1px solid #dfe8ed !important; }
table.dataTable tbody td { padding: 10px 9px; border-top-color: #edf2f6; }
.dataTables_wrapper { font-size: 12px; }
.guide-step { display: flex; gap: 15px; padding: 12px 0 19px; }
.guide-number { flex: 0 0 31px; width: 31px; height: 31px; line-height: 31px;
  border-radius: 50%; background: #D9F0FC; color: @DEEP@; text-align: center; font-weight: 700; }
.content-wrapper a:not(.btn), .right-side a:not(.btn) { color: @PRIMARY@; }
.nav-tabs-custom > .nav-tabs > li.active { border-top-color: @PRIMARY@ !important; }
.pagination > .active > a, .pagination > .active > a:hover,
.pagination > .active > a:focus { background-color: @PRIMARY@; border-color: @PRIMARY@; }
.guide-step strong { color: #2e4754; font-size: 13px; }
.guide-step p { color: #7b8e98; font-size: 12px; margin: 6px 0 0; line-height: 1.8; }
.app-footer { color: #9aaab3; font-size: 10px; letter-spacing: .5px;
  text-align: right; margin-top: 10px; padding-bottom: 8px; }
@media (max-width: 991px) {
  .header-tag { display: none !important; }
  .content { padding: 21px 18px; }
  .spatial-combined-row { grid-template-columns: 1fr; }
}
@media (max-width: 767px) {
  .main-header .logo { text-align: center; }
  .brand-fallback { justify-content: center; }
  .header-caption { font-size: 15px; }
  .page-heading h2 { font-size: 24px; }
  .small-box h3 { font-size: 23px; }
}
"
for (pair in list(c("@PRIMARY@", APP_CONFIG$primary_color),
                  c("@ACCENT@", APP_CONFIG$accent_color),
                  c("@SIDEBAR@", APP_CONFIG$sidebar_color),
                  c("@DEEP@", APP_CONFIG$deep_color),
                  c("@LIGHT@", APP_CONFIG$light_color),
                  c("@BG@", APP_CONFIG$background_color))) {
  css <- gsub(pair[[1]], pair[[2]], css, fixed = TRUE)
}


ui <- dashboardPage(
  dashboardHeader(
    title = brand, titleWidth = 250,
    tags$li(class = "dropdown header-caption", APP_CONFIG$app_title),
    tags$li(class = "dropdown header-tag", "SPATIAL PROTEOMICS WORKSPACE")
  ),
  dashboardSidebar(
    width = 250,
    sidebarMenu(
      id = "main_tabs",
      tags$li(class = "header", "ANALYSIS WORKSPACE"),
      menuItem("Data Input", tabName = "data_input", icon = icon("file-import"), selected = TRUE),
      menuItem("QC Results", tabName = "qc_results", icon = icon("clipboard-check")),
      menuItem("Comparison analysis", tabName = "comparison", icon = icon("scale-balanced")),
      menuItem("Continuous analysis", tabName = "continuous", icon = icon("wave-square")),
      menuItem("Trajectory analysis", tabName = "trajectory", icon = icon("chart-line")),
      menuItem("Report Overview", tabName = "dashboard", icon = icon("chart-bar")),
      menuItem("Data Explorer", tabName = "explorer", icon = icon("table")),
      menuItem("Guide", tabName = "guide", icon = icon("circle-question"))
    ),
    uiOutput("sidebar_file"),
    div(class = "sidebar-note", strong(APP_CONFIG$app_title),
        APP_CONFIG$app_subtitle, tags$br(), "Paired Parquet / QC")
  ),
  dashboardBody(
    tags$head(
      tags$style(HTML(paste0(css, "
#qc_image img { max-width:100%; height:auto; }
#qc_log { max-height:320px; overflow:auto; font-size:11px; }
#comparison_image img { max-width:100%; height:auto; }
#comparison_log { max-height:320px; overflow:auto; font-size:11px; }
#roi_metadata_image img, #zonation_image img { max-width:100%; height:auto; }
#roi_metadata_log, #zonation_log { max-height:320px; overflow:auto; font-size:11px; }
#trajectory_beta_image img, #trajectory_curve_image img { max-width:100%; height:auto; }
#trajectory_log { max-height:320px; overflow:auto; font-size:11px; }

body.fixed .main-sidebar {
  overflow-y: auto;
  height: 100vh;
}
"))),
      tags$script(HTML("
document.addEventListener('DOMContentLoaded', function () {
  document.body.classList.add('fixed');
});
"))
    ),
    tabItems(
      tabItem(
        tabName = "data_input",
        page_heading("PAIRED REPORT INPUT", "DIA-NN QC Input",
                     "Automatically detects the sample number after GF or SPF in Run names and performs QC separately for each sample."),
        uiOutput("input_metrics"),
        fluidRow(
          section_box("01  /  First-pass report", width = 6,
            div(class = "upload-zone",
                div(class = "upload-glyph", icon("file-arrow-up")),
                fileInput("first_pass_file", "report-first-pass.parquet",
                          accept = ".parquet", multiple = FALSE,
                          buttonLabel = "Browse", placeholder = "report-first-pass.parquet"),
                p(class = "muted", paste0("First input · maximum ", APP_CONFIG$upload_limit_gib, " GiB"))),
            uiOutput("first_pass_status")),
          section_box("02  /  Final report", width = 6,
            div(class = "upload-zone",
                div(class = "upload-glyph", icon("file-arrow-up")),
                fileInput("report_file", "report.parquet",
                          accept = ".parquet", multiple = FALSE,
                          buttonLabel = "Browse", placeholder = "report.parquet"),
                p(class = "muted", paste0("Second input · maximum ", APP_CONFIG$upload_limit_gib, " GiB"))),
            uiOutput("report_status"))
        ),
        fluidRow(section_box("03  /  QC paths and execution",
          fluidRow(
            column(6,
              textInput("qc_script_path", "QC script file path (.R)",
                        value = APP_CONFIG$qc_script,
                        placeholder = "D:/Spatial_Proteomics/codes/QC_Code.R", width = "100%")),
            column(6,
              textInput("qc_results_path", "Results directory",
                        value = APP_CONFIG$qc_results_dir,
                        placeholder = "D:/Spatial_Proteomics/Analysis_results", width = "100%"))
          ),
          p(class = "muted",
            "Enter paths on the PC or server running R. Relative paths are resolved from the app.R directory."),
          uiOutput("qc_path_status"),
          textAreaInput(
            "manual_exclusions", "Manual exclude ROIs (optional)", value = "", rows = 4,
            width = "100%",
            placeholder = paste(
              "GF1-4: ROI5, ROI29-ROI31, ROI63",
              "SPF1-2: 3, 7, 10-12",
              sep = "\n"
            )
          ),
          p(class = "muted",
            "Enter one sample per line. Separate ROI labels or numbers with commas; ranges such as ROI29-ROI31 are supported. Leave blank to disable manual exclusions."),
          uiOutput("manual_exclusion_status"),
          uiOutput("qc_readiness"),
          actionButton("run_qc", "Run QC", icon = icon("play"), class = "btn-primary"),
          actionButton("cancel_qc", "Stop", icon = icon("stop")),
          actionButton("calculate_summary", "Summarize final report",
                       icon = icon("calculator"), class = "btn-primary"),
          p(class = "muted", style = "margin-top:14px;",
            "Automatically distinguishes sample names such as GF1-4, GF2-2, and SPF1-2. Multiple samples are analyzed independently. Changing an input file resets existing QC results."),
          uiOutput("summary_status")
        )),
        fluidRow(section_box("04  /  Input preview",
          radioButtons("preview_source", "File to preview",
                       choices = c("First-pass report" = "first_pass", "Final report" = "report"),
                       selected = "first_pass", inline = TRUE),
          p(class = "muted", paste0("Displays up to ", APP_CONFIG$preview_rows, " rows from the selected file.")),
          uiOutput("preview_empty"), DT::DTOutput("preview")
        ))
      ),
      tabItem(
        tabName = "qc_results",
        page_heading("QUALITY CONTROL", "QC Results",
                     "Review the plots, tables, result files, and execution log generated by QC_Code.R."),
        uiOutput("qc_result_status"),
        uiOutput("qc_run_location"),
        fluidRow(section_box("Auto-detected datasets and QC summary",
          p(class = "muted", "Sample names and result prefixes are detected automatically from Runs in both reports. Different samples are not merged even when their ROI numbers match."),
          DT::DTOutput("qc_dataset_summary")
        )),
        fluidRow(
          section_box("QC figures", width = 7,
            selectInput("qc_plot_choice", "Select plot", choices = character()),
            imageOutput("qc_image", height = "auto"),
            uiOutput("qc_plot_empty")),
          section_box("Result files", width = 5,
            DT::DTOutput("qc_file_table"),
            selectInput("qc_file_choice", "Result to download", choices = character()),
            downloadButton("download_qc_file", "Download selected file"))
        ),
        fluidRow(section_box("QC tables",
          selectInput("qc_table_choice", "Select table", choices = character()),
          p(class = "muted", paste0("The screen displays up to ", APP_CONFIG$qc_table_preview_rows,
                                    " rows. Download the complete result from the file list above.")),
          DT::DTOutput("qc_table"), uiOutput("qc_table_empty")
        )),
        fluidRow(section_box("Execution log",
          p(class = "muted", "Displays the first 2,000 log lines. Download the complete log using the button below."),
          downloadButton("download_qc_log", "Download execution log"),
          verbatimTextOutput("qc_log")
        ))
      ),
      tabItem(
        tabName = "comparison",
        page_heading("POST-QC COMPARISON", "Comparison analysis",
                     "After QC, compares the selected GF and SPF sets from report.parquet."),
        uiOutput("comparison_status"),
        fluidRow(section_box("01  /  Comparison settings",
          textInput(
            "comparison_script_path", "Comparison analysis script file path (.R)",
            value = APP_CONFIG$comparison_script,
            placeholder = "D:/Spatial_Proteomics/codes/Comparison_Pipeline.R", width = "100%"
          ),
          uiOutput("comparison_script_status"),
          fluidRow(
            column(6,
              checkboxGroupInput("comparison_gf_sets", "GF groups",
                                 choices = character(), selected = character())),
            column(6,
              checkboxGroupInput("comparison_spf_sets", "SPF groups",
                                 choices = character(), selected = character()))
          ),
          p(class = "muted",
            "Groups are detected automatically from Run names in report.parquet. Select at least three GF and three SPF sets for set-level statistics."),
          uiOutput("comparison_readiness"),
          actionButton("run_comparison", "Run comparison", icon = icon("play"),
                       class = "btn-primary"),
          actionButton("cancel_comparison", "Stop", icon = icon("stop"))
        )),
        uiOutput("comparison_run_location"),
        fluidRow(section_box("02  /  Comparison summary",
          DT::DTOutput("comparison_summary"), uiOutput("comparison_summary_empty")
        )),
        fluidRow(
          section_box("Comparison figures", width = 7,
            selectInput("comparison_plot_choice", "Select plot", choices = character()),
            imageOutput("comparison_image", height = "auto"),
            uiOutput("comparison_plot_empty")),
          section_box("Result files", width = 5,
            DT::DTOutput("comparison_file_table"),
            selectInput("comparison_file_choice", "Result to download", choices = character()),
            downloadButton("download_comparison_file", "Download selected file"))
        ),
        fluidRow(section_box("Comparison tables",
          selectInput("comparison_table_choice", "Select table", choices = character()),
          p(class = "muted", paste0("The screen displays up to ", APP_CONFIG$qc_table_preview_rows,
                                    " rows.")),
          DT::DTOutput("comparison_table"), uiOutput("comparison_table_empty")
        )),
        fluidRow(section_box("Execution log",
          downloadButton("download_comparison_log", "Download execution log"),
          verbatimTextOutput("comparison_log")
        ))
      ),
      tabItem(
        tabName = "continuous",
        page_heading("POST-COMPARISON SPATIAL MODELING", "Continuous analysis",
          "After comparison, generates ROI metadata from the reference-coordinate CSV and applies the resulting Spatial_S to continuous CV-to-PV zonation analysis."),
        uiOutput("continuous_overall_status"),
        tabBox(id = "continuous_analysis_step", width = 12,
          tabPanel("01 Generation ROI metadata",
            fluidRow(section_box("ROI metadata input and run",
              fluidRow(
                column(6,
                  fileInput(
                    "roi_coordinate_file", "Image_Reference_Coordinates.csv",
                    accept = c(".csv", "text/csv"), multiple = FALSE,
                    buttonLabel = "Browse", placeholder = "Select reference-coordinate CSV"
                  )),
                column(6,
                  textInput(
                    "roi_metadata_script_path", "ROI metadata script file path (.R)",
                    value = APP_CONFIG$roi_metadata_script,
                    placeholder = "D:/Spatial_Proteomics/codes/ROI_Metadata_Generation.R",
                    width = "100%"
                  ))
              ),
              p(class = "muted",
                "The CSV must contain Image, Group, Mouse, Set, and X/Y coordinates for ROI01, ROI08, ROI89, CV, and PV."),
              uiOutput("roi_metadata_readiness"),
              actionButton("run_roi_metadata", "Generate ROI metadata",
                           icon = icon("location-dot"), class = "btn-primary"),
              actionButton("cancel_roi_metadata", "Stop", icon = icon("stop")),
              uiOutput("roi_metadata_run_location")
            )),
            fluidRow(section_box("ROI metadata results",
              fluidRow(
                column(6,
                  selectInput("roi_metadata_plot_choice", "Select QC plot",
                              choices = character()),
                  imageOutput("roi_metadata_image", height = "auto"),
                  uiOutput("roi_metadata_plot_empty")),
                column(6,
                  DT::DTOutput("roi_metadata_file_table"),
                  selectInput("roi_metadata_file_choice", "Result to download",
                              choices = character()),
                  downloadButton("download_roi_metadata_file", "Download selected file"))
              ),
              selectInput("roi_metadata_table_choice", "Select table", choices = character()),
              DT::DTOutput("roi_metadata_table"),
              uiOutput("roi_metadata_table_empty"),
              downloadButton("download_roi_metadata_log", "Download execution log"),
              verbatimTextOutput("roi_metadata_log")
            ))
          ),
          tabPanel("02 Continuous CV→PV zonation",
            fluidRow(section_box("Continuous zonation input and run",
              textInput(
                "continuous_zonation_script_path", "Continuous zonation script file path (.R)",
                value = APP_CONFIG$continuous_zonation_script,
                placeholder = "D:/Spatial_Proteomics/codes/Continuous_Zonation_Total_Separation.R",
                width = "100%"
              ),
              p(class = "muted",
                "Automatically uses ROI_metadata.csv generated in Step 1 and the currently connected report.parquet."),
              uiOutput("zonation_readiness"),
              actionButton("run_zonation", "Run Continuous analysis",
                           icon = icon("play"), class = "btn-primary"),
              actionButton("cancel_zonation", "Stop", icon = icon("stop")),
              uiOutput("zonation_run_location")
            )),
            fluidRow(section_box("Continuous zonation results",
              DT::DTOutput("zonation_summary"), uiOutput("zonation_summary_empty"),
              fluidRow(
                column(7,
                  selectInput("zonation_plot_choice", "Select plot",
                              choices = character()),
                  imageOutput("zonation_image", height = "auto"),
                  uiOutput("zonation_plot_empty")),
                column(5,
                  DT::DTOutput("zonation_file_table"),
                  selectInput("zonation_file_choice", "Result to download",
                              choices = character()),
                  downloadButton("download_zonation_file", "Download selected file"))
              ),
              selectInput("zonation_table_choice", "Select table", choices = character()),
              DT::DTOutput("zonation_table"), uiOutput("zonation_table_empty"),
              downloadButton("download_zonation_log", "Download execution log"),
              verbatimTextOutput("zonation_log")
            ))
          )
        )
      ),
      tabItem(
        tabName = "trajectory",
        page_heading("TARGETED SPATIAL VISUALIZATION", "Trajectory analysis",
          "01 Trajectory generates CV-to-PV trajectories, and 02 Spatial Map generates ROI abundance maps over histology images."),
        tabBox(id = "trajectory_analysis_step", width = 12,
          tabPanel("01 Trajectory",
        uiOutput("trajectory_status"),
        fluidRow(section_box("01  /  Select target and run",
          textAreaInput(
            "trajectory_target", "UniProt accession or gene name (single / multiple)",
            value = "", rows = 3,
            placeholder = "Single: Q05421\nMultiple: Q05421, Ehhadh, Hmgcs1",
            width = "100%"
          ),
          fluidRow(
            column(6,
              textInput(
                "zonation_coefficient_script_path", "Zonation coefficient script (.R)",
                value = APP_CONFIG$zonation_coefficient_script,
                placeholder = "D:/Spatial_Proteomics/codes/Zonation_Coefficient_Beta.R",
                width = "100%"
              )),
            column(6,
              textInput(
                "trajectory_selected_script_path", "Trajectory script (.R)",
                value = APP_CONFIG$trajectory_selected_script,
                placeholder = "D:/Spatial_Proteomics/codes/Trajectory_Selected.R",
                width = "100%"
              ))
          ),
          textInput(
            "trajectory_cache_path", "Single-protein cache directory",
            value = APP_CONFIG$trajectory_cache_dir,
            placeholder = "D:/Spatial_Proteomics/Trajectory_Cache",
            width = "100%"
          ),
          p(class = "muted",
            "Separate multiple targets with commas, semicolons, spaces, or line breaks. Targets listed in Cache_Index.csv are loaded immediately without reanalysis; only uncached targets are calculated."),
          uiOutput("trajectory_readiness"),
          actionButton("run_trajectory", "Generate both plots", icon = icon("play"),
                       class = "btn-primary"),
          actionButton("cancel_trajectory", "Stop", icon = icon("stop")),
          uiOutput("trajectory_run_location")
        )),
        fluidRow(
          section_box("Zonation coefficient (Beta)", width = 6,
            selectInput("trajectory_beta_choice", "Target to display", choices = character()),
            imageOutput("trajectory_beta_image", height = "auto"),
            uiOutput("trajectory_beta_empty"),
            downloadButton("download_trajectory_beta_png", "Download PNG")
          ),
          section_box("GF/SPF CV→PV trajectory", width = 6,
            selectInput("trajectory_curve_choice", "Target to display (linked to Beta selection)",
                        choices = character()),
            imageOutput("trajectory_curve_image", height = "auto"),
            uiOutput("trajectory_curve_empty"),
            downloadButton("download_trajectory_curve_png", "Download PNG")
          )
        ),
        fluidRow(section_box("Analysis summary and files",
          DT::DTOutput("trajectory_summary"),
          DT::DTOutput("trajectory_file_table"),
          selectInput("trajectory_file_choice", "Result to download", choices = character()),
          downloadButton("download_trajectory_file", "Download selected file")
        )),
        fluidRow(section_box("Execution log",
          downloadButton("download_trajectory_log", "Download execution log"),
          verbatimTextOutput("trajectory_log")
        ))
          ),
          tabPanel("02 Spatial Map",
            uiOutput("spatial_map_status"),
            fluidRow(section_box("02  /  Spatial Map input and run",
              selectizeInput(
                "spatial_map_target", "UniProt accession or gene name (up to 10)",
                choices = character(), selected = character(), multiple = TRUE,
                width = "100%",
                options = list(
                  create = TRUE,
                  createOnBlur = TRUE,
                  persist = FALSE,
                  maxItems = 10,
                  delimiter = ",",
                  plugins = list("remove_button"),
                  placeholder = "Enter a target and press Enter: Q05421, Ehhadh, Cyp2e1"
                )
              ),
              p(class = "muted",
                "Enter each target and press Enter to add it as a tag. Up to 10 targets can be analyzed at once."),
              textInput(
                "spatial_expression_script_path", "Spatial Expression script (.R)",
                value = APP_CONFIG$spatial_expression_script,
                placeholder = "D:/Spatial_Proteomics/codes/Spatial_Expression_Shiny.R",
                width = "100%"
              ),
              textInput(
                "spatial_histology_root", "Histology image directory",
                value = "",
                placeholder = "D:/Spatial_Proteomics/Histology_Images",
                width = "100%"
              ),
              checkboxInput(
                "spatial_require_all_sets",
                "Require all eight GF/SPF sets", value = TRUE
              ),
              p(class = "muted",
                paste0(
                  "Automatically reuses coordinates from 01 Generation ROI metadata. ",
                  "When Continuous analysis-ready values are available, the same values and normalization as 01 Trajectory are used, ",
                  "and Spatial Map runs independently of whether 01 Trajectory has been run. ",
                  "The Histology directory must contain images named by set, such as GF1-2.png."
                )),
              uiOutput("spatial_map_readiness"),
              actionButton("run_spatial_map", "Generate Spatial Map",
                           icon = icon("map"), class = "btn-primary"),
              actionButton("cancel_spatial_map", "Stop", icon = icon("stop")),
              uiOutput("spatial_map_run_location")
            )),
            fluidRow(section_box("Spatial Map images",
              checkboxInput(
                "spatial_map_show_combined",
                "Show combined Trajectory + Spatial Map panel",
                value = TRUE
              ),
              p(class = "muted",
                "When selected, trajectory panels appear on the left and sample histology Spatial Maps appear as a 2x2 card grid on the right. When cleared, only sample Spatial Maps are shown."),
              checkboxGroupInput(
                "spatial_map_target_choice", "Target to display",
                choices = character(), inline = TRUE
              ),
              div(class = "spatial-map-selector-actions",
                actionButton("spatial_map_target_select_all", "All targets",
                             icon = icon("check-double"), class = "btn-primary btn-sm"),
                actionButton("spatial_map_target_select_first", "First target only",
                             icon = icon("check"), class = "btn-primary btn-sm"),
                actionButton("spatial_map_target_clear", "Clear target selection",
                             icon = icon("xmark"), class = "btn-primary btn-sm")
              ),
              checkboxGroupInput(
                "spatial_map_image_choice", "Samples to display together",
                choices = character(), inline = TRUE
              ),
              div(class = "spatial-map-selector-actions",
                actionButton("spatial_map_select_all", "Select all",
                             icon = icon("check-double"), class = "btn-primary btn-sm"),
                actionButton("spatial_map_select_gf", "GF only",
                             icon = icon("check"), class = "btn-primary btn-sm"),
                actionButton("spatial_map_select_spf", "SPF only",
                             icon = icon("check"), class = "btn-primary btn-sm"),
                actionButton("spatial_map_clear_selection", "Clear selection",
                             icon = icon("xmark"), class = "btn-primary btn-sm")
              ),
              uiOutput("spatial_map_image_grid"),
              downloadButton("download_spatial_map_image", "Download selected images")
            )),
            fluidRow(section_box("Spatial Map result files",
              h4("Target execution summary"),
              DT::DTOutput("spatial_map_target_summary"),
              DT::DTOutput("spatial_map_file_table"),
              selectInput("spatial_map_file_choice", "Result to download",
                          choices = character()),
              downloadButton("download_spatial_map_file", "Download selected file")
            )),
            fluidRow(section_box("Spatial Map execution log",
              downloadButton("download_spatial_map_log", "Download execution log"),
              verbatimTextOutput("spatial_map_log")
            ))
          )
        )
      ),
      tabItem(
        tabName = "dashboard",
        page_heading("DATASET OVERVIEW", "Final Report Overview",
                     "Basic observation summary for the final report.parquet. Review QC_Code.R results under QC Results."),
        uiOutput("dashboard_metrics"),
        fluidRow(section_box("Dataset summary",
          uiOutput("dashboard_status"),
          actionButton("refresh_summary", "Recalculate full summary", icon = icon("arrows-rotate"),
                       class = "btn-primary"),
          p(class = "muted", style = "display:inline; margin-left:12px;",
            "Summary based on the input file · no additional FDR filtering")
        )),
        fluidRow(
          section_box("Protein groups by Run", width = 6,
                      plotOutput("protein_plot", height = "330px")),
          section_box("Report rows by Run", width = 6,
                      plotOutput("rows_plot", height = "330px"))
        ),
        fluidRow(section_box("Run summary",
          p(class = "muted",
            "Report rows is the number of rows in the report. Protein groups is the number of unique group identifiers and may differ from the number of individual proteins."),
          uiOutput("run_limit_note"), DT::DTOutput("run_table")
        ))
      ),
      tabItem(
        tabName = "explorer",
        page_heading("INSPECT YOUR DATA", "Data Explorer",
                     "Select columns and search criteria to inspect data from the final report.parquet."),
        fluidRow(
          section_box("Query options", width = 3,
            selectizeInput("explorer_columns", "Columns to display", choices = NULL, multiple = TRUE),
            textInput("run_filter", "Run name contains", placeholder = "Example: GF1-3"),
            textInput("protein_filter", "Protein group contains", placeholder = "Example: P12345"),
            numericInput("explorer_limit", "Maximum rows to display", value = 100,
                         min = 20, max = APP_CONFIG$max_explorer_rows, step = 20),
            actionButton("query_preview", "Query data", icon = icon("magnifying-glass"),
                         class = "btn-primary btn-block"),
            p(class = "muted", style = "margin-top:12px;",
              "Text searches are case-sensitive. Query results contain only the rows within the display limit.")
          ),
          section_box("Query results", width = 9,
            textOutput("explorer_note"),
            DT::DTOutput("explorer_table"),
            tags$br(),
            downloadButton("download_preview", "Current query results CSV")
          )
        )
      ),

      tabItem(
        tabName = "guide",
        page_heading("WORKSPACE GUIDE", "Getting Started", "Run the complete QC-to-spatial-analysis workflow"),
        fluidRow(
          section_box("Recommended workflow", width = 7,
            div(class = "guide-step", span(class = "guide-number", "1"),
                div(strong("Upload two Parquet files"),
                    p("Select report-first-pass.parquet in the first field and report.parquet in the second, then specify the QC script path and results directory."))),
            div(class = "guide-step", span(class = "guide-number", "2"),
                div(strong("Run QC"),
                    p("If needed, enter Manual exclude ROIs for each sample and click Run QC. Execution status and logs are available under QC Results."))),
            div(class = "guide-step", span(class = "guide-number", "3"),
                div(strong("Review and download results"),
                    p("All results are stored in separate QC, Comparison, Continuous, and Trajectory stage directories under date/time sessions in Analysis_results."))),
            div(class = "guide-step", span(class = "guide-number", "4"),
                div(strong("Comparison analysis"),
                    p("After QC, select the GF/SPF groups detected automatically from report.parquet on the Comparison tab and run the comparison analysis."))),
            div(class = "guide-step", span(class = "guide-number", "5"),
                div(strong("Continuous CV-to-PV analysis"),
                    p("After comparison, open Continuous analysis. In 01 Generation ROI metadata, upload Image_Reference_Coordinates.csv and generate ROI_metadata.csv. Then run 02 Continuous CV-to-PV zonation using the generated Spatial_S coordinates."))),
            div(class = "guide-step", span(class = "guide-number", "6"),
                div(strong("Trajectory and Spatial Map"),
                    p("In 01 Trajectory, enter one or more UniProt accessions or gene names to generate zonation-coefficient and GF/SPF trajectory plots. In 02 Spatial Map, enter targets and select the histology image directory to generate sample-level histology overlays. Spatial Map can run independently of 01 Trajectory after ROI metadata is available.")))
          ),
          section_box("Input and result notes", width = 5,
            tags$dl(class = "file-meta",
              tags$dt("Input 1"), tags$dd("report-first-pass.parquet"),
              tags$dt("Input 2"), tags$dd("report.parquet"),
              tags$dt("File limit"), tags$dd(paste("per file", APP_CONFIG$upload_limit_gib, "GiB")),
              tags$dt("QC analysis"), tags$dd("QC_Code.R / paired Parquet inputs / optional excluded ROIs"),
              tags$dt("Comparison analysis"), tags$dd("Comparison_Pipeline.R / QC-filtered report / selected GF and SPF sets"),
              tags$dt("ROI metadata"), tags$dd("Image_Reference_Coordinates.csv -> ROI_metadata.csv"),
              tags$dt("Continuous analysis"), tags$dd("ROI metadata + report.parquet -> CV-to-PV zonation results"),
              tags$dt("Trajectory analysis"), tags$dd("UniProt accession or gene name / continuous-analysis results / optional protein cache"),
              tags$dt("Spatial Map"), tags$dd("ROI metadata + report.parquet + sample-matched histology images"),
              tags$dt("Result location"), tags$dd("Analysis_results / analysis session / stage-specific folders")),
            p(class = "muted",
              "QC calculations and decision criteria are defined in QC_Code.R. Manually excluded ROIs are removed from the QC-filtered data used by downstream comparison analysis. The app summary and Data Explorer continue to show the connected final report.parquet and are separate from QC-filtered results."),
            p(class = "muted",
              "Upload Image_Reference_Coordinates.csv once in Continuous analysis. The saved ROI metadata is reused by Spatial Map, so it does not need to be uploaded again. Spatial Map can run independently of 01 Trajectory after ROI metadata is available."),
            p(class = "muted",
              "Connected inputs are used throughout the current analysis session. Results are preserved in stage-specific folders under the corresponding session in Analysis_results and are not deleted when the app closes.")
          )
        )
      )
    ),
    div(class = "app-footer", APP_CONFIG$company_name, " / ", APP_CONFIG$app_title)
  ),
  skin = "blue", title = paste(APP_CONFIG$company_name, APP_CONFIG$app_title, sep = " | ")
)

server <- function(input, output, session) {
  temp_parent <- getOption("spatial.temp_parent", tempdir())
  dir.create(temp_parent, recursive = TRUE, showWarnings = FALSE)
  work_dir <- tempfile("spatial_diann_", tmpdir = temp_parent)
  if (!dir.create(work_dir)) stop("Cannot create the temporary working directory.")
  spill_dir <- file.path(work_dir, "spill"); dir.create(spill_dir)
  con <- tryCatch(DBI::dbConnect(
    duckdb::duckdb(), dbdir = file.path(work_dir, "session.duckdb"), bigint = "numeric",
    config = list(memory_limit = APP_CONFIG$db_memory_limit, threads = "2",
                  preserve_insertion_order = "false",
                  temp_directory = normalizePath(spill_dir, winslash = "/", mustWork = TRUE))
  ), error = function(e) { unlink(work_dir, recursive = TRUE); stop(e) })

  source_pair <- reactiveVal(list(first_pass = NULL, report = NULL))
  source_info <- reactive(source_pair()$report) # Existing overview/explorer: final report.
  upload_errors <- reactiveVal(list(first_pass = NULL, report = NULL))
  summary_result <- reactiveVal(NULL); run_result <- reactiveVal(NULL)
  summary_error <- reactiveVal(NULL); explorer_result <- reactiveVal(NULL)
  explorer_message <- reactiveVal("Connect the final report.parquet before querying.")
  qc_job <- reactiveVal(NULL); qc_state <- reactiveVal("idle")
  qc_message <- reactiveVal("Upload both files, then run QC.")
  qc_manifest <- reactiveVal(NULL); qc_output_dir <- reactiveVal(NULL)
  qc_run_paths <- reactiveVal(NULL)
  qc_log_tick <- reactiveVal(0L)
  comparison_groups <- reactiveVal(list(GF = character(), SPF = character()))
  comparison_job <- reactiveVal(NULL); comparison_state <- reactiveVal("idle")
  comparison_message <- reactiveVal("Select report.parquet groups and complete QC.")
  comparison_manifest <- reactiveVal(NULL); comparison_output_dir <- reactiveVal(NULL)
  comparison_run_paths <- reactiveVal(NULL); comparison_log_tick <- reactiveVal(0L)
  continuous_run_paths <- reactiveVal(NULL)
  roi_metadata_job <- reactiveVal(NULL); roi_metadata_state <- reactiveVal("idle")
  roi_metadata_message <- reactiveVal("After comparison is complete, connect the reference-coordinate CSV.")
  roi_metadata_manifest <- reactiveVal(NULL); roi_metadata_file <- reactiveVal(NULL)
  roi_metadata_log_tick <- reactiveVal(0L)
  zonation_job <- reactiveVal(NULL); zonation_state <- reactiveVal("idle")
  zonation_message <- reactiveVal("Generate ROI metadata first.")
  zonation_manifest <- reactiveVal(NULL); zonation_log_tick <- reactiveVal(0L)
  trajectory_job <- reactiveVal(NULL); trajectory_state <- reactiveVal("idle")
  trajectory_message <- reactiveVal("Complete Continuous analysis first.")
  trajectory_manifest <- reactiveVal(NULL); trajectory_run_paths <- reactiveVal(NULL)
  trajectory_log_tick <- reactiveVal(0L)
  spatial_map_job <- reactiveVal(NULL); spatial_map_state <- reactiveVal("idle")
  spatial_map_message <- reactiveVal("Check report.parquet and the Spatial Map inputs.")
  spatial_map_manifest <- reactiveVal(NULL); spatial_map_run_paths <- reactiveVal(NULL)
  spatial_map_log_tick <- reactiveVal(0L)

  kill_process <- function(job) {
    if (is.null(job)) return(invisible(NULL))
    if (tryCatch(job$is_alive(), error = function(e) FALSE)) {
      tryCatch(job$kill_tree(), error = function(e) try(job$kill(), silent = TRUE))
      try(job$wait(timeout = 1000), silent = TRUE)
    }
    invisible(NULL)
  }
  session$onSessionEnded(function() {
    kill_process(isolate(qc_job()))
    kill_process(isolate(comparison_job()))
    kill_process(isolate(roi_metadata_job()))
    kill_process(isolate(zonation_job()))
    kill_process(isolate(trajectory_job()))
    kill_process(isolate(spatial_map_job()))
    try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
    # Only session database/work files are deleted; persistent QC results stay intact.
    unlink(work_dir, recursive = TRUE)
  })
  reset_trajectory <- function(message = "Complete Continuous analysis first.") {
    kill_process(isolate(trajectory_job()))
    trajectory_job(NULL); trajectory_state("idle"); trajectory_manifest(NULL)
    trajectory_run_paths(NULL); trajectory_message(message)
    updateSelectInput(session, "trajectory_file_choice", choices = character())
    updateSelectInput(session, "trajectory_beta_choice", choices = character())
    updateSelectInput(session, "trajectory_curve_choice", choices = character())
  }
  reset_spatial_map <- function(message = "Check report.parquet and the Spatial Map inputs.") {
    kill_process(isolate(spatial_map_job()))
    spatial_map_job(NULL); spatial_map_state("idle"); spatial_map_manifest(NULL)
    spatial_map_run_paths(NULL); spatial_map_message(message)
    updateCheckboxGroupInput(session, "spatial_map_target_choice",
                             choices = character(), selected = character())
    updateCheckboxGroupInput(session, "spatial_map_image_choice",
                             choices = character(), selected = character())
    updateSelectInput(session, "spatial_map_file_choice", choices = character())
  }
  reset_continuous <- function(message = "After comparison is complete, connect the reference-coordinate CSV.") {
    kill_process(isolate(roi_metadata_job()))
    kill_process(isolate(zonation_job()))
    roi_metadata_job(NULL); roi_metadata_state("idle"); roi_metadata_manifest(NULL)
    roi_metadata_file(NULL); roi_metadata_message(message)
    zonation_job(NULL); zonation_state("idle"); zonation_manifest(NULL)
    zonation_message("Generate ROI metadata first.")
    reset_trajectory("Continuous analysis results have changed. Run the analysis again.")
    reset_spatial_map("Continuous analysis results have changed. Check the inputs again.")
    continuous_run_paths(NULL)
    for (id in c("roi_metadata_plot_choice", "roi_metadata_table_choice",
                 "roi_metadata_file_choice", "zonation_plot_choice",
                 "zonation_table_choice", "zonation_file_choice")) {
      updateSelectInput(session, id, choices = character())
    }
  }
  reset_qc <- function() {
    kill_process(isolate(qc_job()))
    kill_process(isolate(comparison_job()))
    qc_job(NULL); qc_state("idle"); qc_manifest(NULL); qc_output_dir(NULL); qc_run_paths(NULL)
    qc_message("The input files have changed. Check both files, then run QC.")
    comparison_job(NULL); comparison_state("idle"); comparison_manifest(NULL)
    comparison_output_dir(NULL); comparison_run_paths(NULL)
    comparison_message("The input files have changed. Complete QC, then run comparison analysis.")
    reset_continuous("The input files have changed. Run comparison again.")
    updateSelectInput(session, "qc_plot_choice", choices = character())
    updateSelectInput(session, "qc_table_choice", choices = character())
    updateSelectInput(session, "qc_file_choice", choices = character())
    updateSelectInput(session, "comparison_plot_choice", choices = character())
    updateSelectInput(session, "comparison_table_choice", choices = character())
    updateSelectInput(session, "comparison_file_choice", choices = character())
  }
  qc_path_settings <- reactive({
    invalidateLater(2000, session)
    tryCatch(inspect_qc_paths(input$qc_script_path, input$qc_results_path),
             error = function(e) list(error = conditionMessage(e)))
  })
  manual_exclusion_settings <- reactive({
    tryCatch(list(value = parse_manual_exclusions(input$manual_exclusions)),
             error = function(e) list(error = conditionMessage(e)))
  })
  comparison_script_settings <- reactive({
    invalidateLater(2000, session)
    tryCatch(list(path = inspect_comparison_script(input$comparison_script_path)),
             error = function(e) list(error = conditionMessage(e)))
  })
  roi_metadata_script_settings <- reactive({
    invalidateLater(2000, session)
    tryCatch(list(path = inspect_analysis_script(
      input$roi_metadata_script_path, "ROI metadata")),
      error = function(e) list(error = conditionMessage(e)))
  })
  zonation_script_settings <- reactive({
    invalidateLater(2000, session)
    tryCatch(list(path = inspect_analysis_script(
      input$continuous_zonation_script_path, "Continuous zonation")),
      error = function(e) list(error = conditionMessage(e)))
  })
  trajectory_script_settings <- reactive({
    invalidateLater(2000, session)
    tryCatch(list(
      beta_path = inspect_analysis_script(
        input$zonation_coefficient_script_path, "Zonation coefficient"),
      curve_path = inspect_analysis_script(
        input$trajectory_selected_script_path, "Trajectory")
    ), error = function(e) list(error = conditionMessage(e)))
  })
  spatial_expression_script_settings <- reactive({
    invalidateLater(2000, session)
    tryCatch({
      path <- inspect_analysis_script(
        input$spatial_expression_script_path, "Spatial Expression")
      script_text <- readLines(path, warn = FALSE, encoding = "UTF-8")
      if (!any(grepl("SPATIAL_EXPRESSION_SHINY_API", script_text, fixed = TRUE)) ||
          !any(grepl("SPATIAL_EXPRESSION_MULTI_TARGET_API", script_text,
                     fixed = TRUE)) ||
          !any(grepl("SPATIAL_MAP_TRAJECTORY_PANEL_API", script_text,
                     fixed = TRUE)) ||
          !any(grepl("SPATIAL_MAP_SAMPLE_LABEL_API", script_text,
                     fixed = TRUE)) ||
          !any(grepl("2026-11-17-v7-trajectory-consistency", script_text,
                     fixed = TRUE)) ||
          !any(grepl("SPATIAL_MAP_COMBINED_SAMPLE_HEADER_API", script_text,
                     fixed = TRUE)) ||
          !any(grepl("SPATIAL_MAP_LEGEND_LAYOUT_API", script_text,
                     fixed = TRUE)) ||
          !any(grepl("SPATIAL_MAP_SHINY_SPLIT_PANEL_API", script_text,
                     fixed = TRUE)) ||
          !any(grepl("SPATIAL_MAP_TRAJECTORY_DATA_API", script_text,
                     fixed = TRUE))) {
        stop("Select the latest Spatial_Expression_Shiny.R that supports ROI-metadata-based trajectory matching.",
             call. = FALSE)
      }
      list(path = path)
    }, error = function(e) list(error = conditionMessage(e)))
  })

  handle_upload <- function(file, key) {
    reset_qc()
    pair <- isolate(source_pair()); pair[[key]] <- NULL; source_pair(pair)
    errors <- isolate(upload_errors()); errors[[key]] <- NULL; upload_errors(errors)
    view <- if (key == "first_pass") "diann_first_pass" else "diann_input"
    expected <- if (key == "first_pass") "report-first-pass.parquet" else "report.parquet"
    if (key == "report") {
      comparison_groups(list(GF = character(), SPF = character()))
      updateCheckboxGroupInput(session, "comparison_gf_sets",
                               choices = character(), selected = character())
      updateCheckboxGroupInput(session, "comparison_spf_sets",
                               choices = character(), selected = character())
      summary_result(NULL); run_result(NULL); summary_error(NULL); explorer_result(NULL)
      explorer_message("Checking the final report.")
      updateSelectizeInput(session, "explorer_columns", choices = character(),
                           selected = character(), server = TRUE)
      updateTextInput(session, "run_filter", value = "")
      updateTextInput(session, "protein_filter", value = "")
    }
    result <- tryCatch(withProgress(message = paste("Checking file:", expected), value = 0, {
      DBI::dbExecute(con, paste("DROP VIEW IF EXISTS", view))
      if (tolower(file$name) != tolower(expected)) {
        stop(paste("Select", expected, "in this input field."), call. = FALSE)
      }
      if (file$size <= 0) stop("The file is empty.", call. = FALSE)
      if (file$size > APP_CONFIG$upload_limit_gib * 1024^3) stop("The per-file upload limit was exceeded.", call. = FALSE)
      required <- c("Run", "Protein.Group", "Q.Value", "PG.Q.Value", "PG.MaxLFQ")
      columns <- register_diann_source(con, file$datapath, view, required)
      incProgress(0.5)
      preview <- DBI::dbGetQuery(con, sprintf(
        "SELECT * FROM %s LIMIT %d", view, APP_CONFIG$preview_rows))
      if (!nrow(preview)) stop("No data rows were found.", call. = FALSE)
      groups <- if (key == "report") query_comparison_groups(con) else NULL
      incProgress(0.5)
      list(name = file$name, size = file$size, columns = columns, preview = preview,
           comparison_groups = groups, extension = "parquet",
           path = normalizePath(file$datapath, winslash = "/", mustWork = TRUE))
    }), error = function(e) e)
    if (inherits(result, "error")) {
      try(DBI::dbExecute(con, paste("DROP VIEW IF EXISTS", view)), silent = TRUE)
      errors <- isolate(upload_errors()); errors[[key]] <- conditionMessage(result)
      upload_errors(errors)
      if (key == "report") explorer_message("Resolve the final report file-validation error.")
    } else {
      pair <- isolate(source_pair()); pair[[key]] <- result; source_pair(pair)
      updateRadioButtons(session, "preview_source", selected = key)
      if (key == "report") {
        comparison_groups(result$comparison_groups)
        updateCheckboxGroupInput(
          session, "comparison_gf_sets",
          choices = result$comparison_groups$GF, selected = result$comparison_groups$GF
        )
        updateCheckboxGroupInput(
          session, "comparison_spf_sets",
          choices = result$comparison_groups$SPF, selected = result$comparison_groups$SPF
        )
        preferred <- intersect(c("Run", "Protein.Group", "Precursor.Id", "Genes",
                                 "PG.MaxLFQ", "Q.Value", "Global.PG.Q.Value"), result$columns)
        updateSelectizeInput(session, "explorer_columns", choices = result$columns,
                             selected = preferred, server = TRUE)
        explorer_message("Select columns and search criteria, then click Query data.")
      }
    }
  }
  observeEvent(input$first_pass_file, handle_upload(input$first_pass_file, "first_pass"))
  observeEvent(input$report_file, handle_upload(input$report_file, "report"))
  observeEvent(input$manual_exclusions, {
    if (identical(qc_state(), "idle") && is.null(qc_manifest())) return()
    reset_qc()
    qc_message("Manual exclude settings have changed. Run QC again.")
    comparison_message("Complete QC again to apply the changed excluded ROIs.")
    showNotification("Existing QC/Comparison result links were reset because Manual exclude settings changed.",
                     type = "message", duration = 7)
  }, ignoreInit = TRUE)

  observeEvent(input$run_qc, {
    if (!is.null(comparison_job())) {
      showNotification("Comparison analysis is running.", type = "message"); return()
    }
    if (!is.null(qc_job())) {
      showNotification("QC is running.", type = "message"); return()
    }
    pair <- source_pair()
    if (is.null(pair$first_pass) || is.null(pair$report)) {
      showNotification("Upload and validate both Parquet files.", type = "error"); return()
    }
    manual_settings <- manual_exclusion_settings()
    if (!is.null(manual_settings$error)) {
      showNotification(manual_settings$error, type = "error", duration = 10); return()
    }
    manual_exclusions <- manual_settings$value
    paths <- tryCatch({
      value <- inspect_qc_paths(input$qc_script_path, input$qc_results_path)
      value$run_dir <- create_qc_run_dir(value$results_root)
      value$results_root <- normalizePath(value$results_root, winslash = "/", mustWork = TRUE)
      value$manual_exclusions <- format_manual_exclusions(manual_exclusions)
      value
    }, error = function(e) e)
    if (inherits(paths, "error")) {
      showNotification(conditionMessage(paths), type = "error", duration = 10)
      return()
    }
    # Snapshot paths for this run. Editing the fields affects only the next run.
    qc_manifest(NULL)
    updateSelectInput(session, "qc_plot_choice", choices = character())
    updateSelectInput(session, "qc_table_choice", choices = character())
    updateSelectInput(session, "qc_file_choice", choices = character())
    out <- paths$run_dir
    qc_run_paths(paths); qc_output_dir(out)
    qc_state("running"); qc_message(paste(basename(paths$script_path), "is running."))
    result <- tryCatch(callr::r_bg(
      func = qc_worker,
      args = list(first_pass_path = pair$first_pass$path, report_path = pair$report$path,
                  output_dir = out, script_path = paths$script_path,
                  manual_exclude_by_dataset = manual_exclusions,
                  db_memory_limit = APP_CONFIG$db_memory_limit, png_dpi = APP_CONFIG$qc_png_dpi),
      wd = out, stdout = file.path(out, "QC_execution.log"), stderr = "2>&1",
      supervise = TRUE, user_profile = FALSE, system_profile = FALSE
    ), error = function(e) e)
    if (inherits(result, "error")) {
      qc_state("error"); qc_message(conditionMessage(result))
    } else {
      qc_job(result)
      updateTabItems(session, "main_tabs", selected = "qc_results")
    }
  })
  observeEvent(input$cancel_qc, {
    if (is.null(qc_job())) return()
    kill_process(qc_job()); qc_job(NULL)
    qc_state("cancelled"); qc_message("QC was stopped by the user.")
  })
  observe({
    job <- qc_job()
    if (is.null(job)) return()
    invalidateLater(750, session)
    qc_log_tick(isolate(qc_log_tick()) + 1L)
    if (tryCatch(job$is_alive(), error = function(e) FALSE)) return()
    qc_job(NULL)
    result <- tryCatch(job$get_result(), error = function(e) e)
    if (inherits(result, "error")) {
      qc_state("error"); qc_message(conditionMessage(result))
    } else {
      # A compact manifest returns from the worker; figures/tables stay on disk.
      qc_manifest(result); qc_state("success")
      ids <- if (!is.null(result$dataset_summary)) result$dataset_summary$Dataset_ID else character()
      excluded_count <- if (is.data.frame(result$excluded_rois)) {
        nrow(result$excluded_rois)
      } else 0L
      qc_message(paste(
        "QC complete.",
        if (length(ids)) paste("Detected samples:", paste(ids, collapse = ", ")) else "",
        paste0("ROIs excluded from comparison: ", excluded_count, "."),
        "Review the summary, plots, and tables for each sample."
      ))
      files <- result$files
      ids <- as.character(seq_len(nrow(files)))
      updateSelectInput(session, "qc_file_choice",
                         choices = stats::setNames(ids, files$name))
      plot_rows <- which(files$kind == "plot")
      table_rows <- which(files$kind == "table")
      updateSelectInput(session, "qc_plot_choice",
                         choices = stats::setNames(ids[plot_rows], files$name[plot_rows]))
      updateSelectInput(session, "qc_table_choice",
                         choices = stats::setNames(ids[table_rows], files$name[table_rows]))
    }
  })

  observeEvent(input$run_comparison, {
    if (!is.null(comparison_job())) {
      showNotification("Comparison analysis is running.", type = "message"); return()
    }
    if (!is.null(qc_job()) || qc_state() == "running") {
      showNotification("Run comparison analysis after QC finishes.", type = "error"); return()
    }
    if (!identical(qc_state(), "success") || is.null(qc_output_dir())) {
      showNotification("Complete QC successfully before running comparison analysis.",
                       type = "error", duration = 8); return()
    }
    report <- source_info()
    if (is.null(report)) {
      showNotification("Connect report.parquet.", type = "error"); return()
    }
    gf_sets <- unique(as.character(input$comparison_gf_sets))
    spf_sets <- unique(as.character(input$comparison_spf_sets))
    qc_result <- qc_manifest()
    excluded_rois <- if (!is.null(qc_result) &&
                            is.data.frame(qc_result$excluded_rois)) {
      qc_result$excluded_rois
    } else {
      data.frame()
    }
    detected <- comparison_groups()
    if (length(gf_sets) < 3L || length(spf_sets) < 3L) {
      showNotification("Select at least three GF and three SPF groups.",
                       type = "error", duration = 8); return()
    }
    if (!all(gf_sets %in% detected$GF) || !all(spf_sets %in% detected$SPF)) {
      showNotification("Only groups detected in report.parquet can be selected.",
                       type = "error", duration = 8); return()
    }
    paths <- tryCatch({
      script <- inspect_comparison_script(input$comparison_script_path)
      out <- create_comparison_run_dir(qc_output_dir())
      list(
        script_path = script, run_dir = out,
        gf_sets = paste(gf_sets, collapse = ", "),
        spf_sets = paste(spf_sets, collapse = ", "),
        excluded_count = nrow(excluded_rois)
      )
    }, error = function(e) e)
    if (inherits(paths, "error")) {
      showNotification(conditionMessage(paths), type = "error", duration = 10); return()
    }
    reset_continuous("Comparison was rerun. Generate ROI metadata after it completes.")
    comparison_manifest(NULL)
    updateSelectInput(session, "comparison_plot_choice", choices = character())
    updateSelectInput(session, "comparison_table_choice", choices = character())
    updateSelectInput(session, "comparison_file_choice", choices = character())
    comparison_run_paths(paths); comparison_output_dir(paths$run_dir)
    comparison_state("running")
    comparison_message(paste(basename(paths$script_path), "is running."))
    result <- tryCatch(callr::r_bg(
      func = comparison_worker,
      args = list(
        report_path = report$path, output_dir = paths$run_dir,
        script_path = paths$script_path, gf_sets = gf_sets, spf_sets = spf_sets,
        excluded_rois = excluded_rois,
        db_memory_limit = APP_CONFIG$db_memory_limit,
        png_dpi = APP_CONFIG$qc_png_dpi
      ),
      wd = paths$run_dir,
      stdout = file.path(paths$run_dir, "Comparison_execution.log"),
      stderr = "2>&1", supervise = TRUE,
      user_profile = FALSE, system_profile = FALSE
    ), error = function(e) e)
    if (inherits(result, "error")) {
      comparison_state("error"); comparison_message(conditionMessage(result))
    } else {
      comparison_job(result)
    }
  })

  observeEvent(input$cancel_comparison, {
    if (is.null(comparison_job())) return()
    kill_process(comparison_job()); comparison_job(NULL)
    comparison_state("cancelled")
    comparison_message("Comparison analysis was stopped by the user.")
  })

  observe({
    job <- comparison_job()
    if (is.null(job)) return()
    invalidateLater(750, session)
    comparison_log_tick(isolate(comparison_log_tick()) + 1L)
    if (tryCatch(job$is_alive(), error = function(e) FALSE)) return()
    comparison_job(NULL)
    result <- tryCatch(job$get_result(), error = function(e) e)
    if (inherits(result, "error")) {
      comparison_state("error"); comparison_message(conditionMessage(result))
    } else {
      comparison_manifest(result); comparison_state("success")
      comparison_message(paste0(
        "Comparison complete. GF: ", paste(result$gf_sets, collapse = ", "),
        " / SPF: ", paste(result$spf_sets, collapse = ", ")
      ))
      files <- result$files
      ids <- as.character(seq_len(nrow(files)))
      plot_rows <- which(files$kind == "plot")
      table_rows <- which(files$kind == "table")
      updateSelectInput(session, "comparison_file_choice",
                        choices = stats::setNames(ids, files$name))
      updateSelectInput(session, "comparison_plot_choice",
                        choices = stats::setNames(ids[plot_rows], files$name[plot_rows]))
      updateSelectInput(session, "comparison_table_choice",
                        choices = stats::setNames(ids[table_rows], files$name[table_rows]))
    }
  })

  observeEvent(input$roi_coordinate_file, {
    reset_continuous("A new reference-coordinate CSV was selected. Generate ROI metadata.")
  }, ignoreInit = TRUE)

  observeEvent(input$run_roi_metadata, {
    if (!is.null(roi_metadata_job()) || !is.null(zonation_job())) {
      showNotification("A Continuous analysis job is running.", type = "message")
      return()
    }
    if (!identical(comparison_state(), "success") || is.null(comparison_output_dir())) {
      showNotification("Complete comparison analysis successfully before running this step.",
                       type = "error", duration = 8)
      return()
    }
    coordinate <- input$roi_coordinate_file
    if (is.null(coordinate) || !file.exists(coordinate$datapath)) {
      showNotification("Select Image_Reference_Coordinates.csv.", type = "error")
      return()
    }
    if (tolower(tools::file_ext(coordinate$name)) != "csv") {
      showNotification("The reference-coordinate file must be CSV.", type = "error")
      return()
    }
    paths <- tryCatch({
      script <- inspect_analysis_script(input$roi_metadata_script_path, "ROI metadata")
      dirs <- create_continuous_run_dir(comparison_output_dir())
      stored_coordinate <- file.path(
        dirs$metadata_dir, "Image_Reference_coordinates.csv"
      )
      if (!file.copy(coordinate$datapath, stored_coordinate, overwrite = TRUE)) {
        stop("Failed to preserve the reference-coordinate CSV in the Analysis results directory.",
             call. = FALSE)
      }
      c(dirs, list(metadata_script = script,
                   coordinate_name = coordinate$name,
                   coordinate_path = normalizePath(stored_coordinate, winslash = "/",
                                                   mustWork = TRUE)))
    }, error = function(e) e)
    if (inherits(paths, "error")) {
      showNotification(conditionMessage(paths), type = "error", duration = 10)
      return()
    }
    reset_continuous("Generating ROI metadata.")
    continuous_run_paths(paths)
    roi_metadata_state("running")
    roi_metadata_message(paste(basename(paths$metadata_script), "is running."))
    result <- tryCatch(callr::r_bg(
      func = roi_metadata_worker,
      args = list(
        input_file = paths$coordinate_path,
        output_dir = paths$metadata_dir,
        script_path = paths$metadata_script
      ),
      wd = paths$metadata_dir,
      stdout = file.path(paths$metadata_dir, "ROI_Metadata_execution.log"),
      stderr = "2>&1", supervise = TRUE,
      user_profile = FALSE, system_profile = FALSE
    ), error = function(e) e)
    if (inherits(result, "error")) {
      roi_metadata_state("error"); roi_metadata_message(conditionMessage(result))
    } else {
      roi_metadata_job(result)
    }
  })

  observeEvent(input$cancel_roi_metadata, {
    if (is.null(roi_metadata_job())) return()
    kill_process(roi_metadata_job()); roi_metadata_job(NULL)
    roi_metadata_state("cancelled")
    roi_metadata_message("ROI metadata generation was stopped by the user.")
  })

  observe({
    job <- roi_metadata_job()
    if (is.null(job)) return()
    invalidateLater(750, session)
    roi_metadata_log_tick(isolate(roi_metadata_log_tick()) + 1L)
    if (tryCatch(job$is_alive(), error = function(e) FALSE)) return()
    roi_metadata_job(NULL)
    result <- tryCatch(job$get_result(), error = function(e) e)
    if (inherits(result, "error")) {
      roi_metadata_state("error"); roi_metadata_message(conditionMessage(result))
    } else {
      roi_metadata_manifest(result); roi_metadata_file(result$metadata_file)
      roi_metadata_state("success")
      roi_metadata_message("ROI metadata generation complete. Continuous zonation can now be run.")
      files <- result$files; ids <- as.character(seq_len(nrow(files)))
      updateSelectInput(session, "roi_metadata_file_choice",
                        choices = stats::setNames(ids, files$name))
      plot_rows <- which(files$kind == "plot")
      table_rows <- which(files$kind == "table")
      updateSelectInput(session, "roi_metadata_plot_choice",
                        choices = stats::setNames(ids[plot_rows], files$name[plot_rows]))
      updateSelectInput(session, "roi_metadata_table_choice",
                        choices = stats::setNames(ids[table_rows], files$name[table_rows]))
    }
  })

  observeEvent(input$run_zonation, {
    if (!is.null(zonation_job()) || !is.null(roi_metadata_job())) {
      showNotification("A Continuous analysis job is running.", type = "message")
      return()
    }
    if (!identical(comparison_state(), "success")) {
      showNotification("Complete comparison analysis first.", type = "error")
      return()
    }
    if (!identical(roi_metadata_state(), "success") ||
        is.null(roi_metadata_file()) || !file.exists(roi_metadata_file())) {
      showNotification("Generate ROI metadata successfully first.", type = "error")
      return()
    }
    report <- source_info(); paths <- continuous_run_paths()
    if (is.null(report) || is.null(paths)) {
      showNotification("report.parquet or Continuous run information is unavailable.", type = "error")
      return()
    }
    prepared <- tryCatch({
      script <- inspect_analysis_script(input$continuous_zonation_script_path,
                                        "Continuous zonation")
      out <- paths$zonation_dir
      if (length(list.files(out, all.files = FALSE))) {
        out <- create_timestamped_run_dir(paths$root,
                                          "02_Continuous_Zonation")
      }
      list(script = script, output_dir = out)
    }, error = function(e) e)
    if (inherits(prepared, "error")) {
      showNotification(conditionMessage(prepared), type = "error", duration = 10)
      return()
    }
    paths$zonation_dir <- prepared$output_dir
    paths$zonation_script <- prepared$script
    continuous_run_paths(paths)
    reset_trajectory("Rerunning Continuous zonation.")
    zonation_manifest(NULL); zonation_state("running")
    zonation_message(paste(basename(prepared$script), "is running."))
    for (id in c("zonation_plot_choice", "zonation_table_choice",
                 "zonation_file_choice")) updateSelectInput(session, id, choices = character())
    result <- tryCatch(callr::r_bg(
      func = continuous_zonation_worker,
      args = list(
        report_path = report$path,
        metadata_path = roi_metadata_file(),
        output_dir = prepared$output_dir,
        script_path = prepared$script
      ),
      wd = prepared$output_dir,
      stdout = file.path(prepared$output_dir, "Continuous_Zonation_execution.log"),
      stderr = "2>&1", supervise = TRUE,
      user_profile = FALSE, system_profile = FALSE
    ), error = function(e) e)
    if (inherits(result, "error")) {
      zonation_state("error"); zonation_message(conditionMessage(result))
    } else {
      zonation_job(result)
    }
  })

  observeEvent(input$cancel_zonation, {
    if (is.null(zonation_job())) return()
    kill_process(zonation_job()); zonation_job(NULL)
    zonation_state("cancelled")
    zonation_message("Continuous zonation was stopped by the user.")
  })

  observe({
    job <- zonation_job()
    if (is.null(job)) return()
    invalidateLater(750, session)
    zonation_log_tick(isolate(zonation_log_tick()) + 1L)
    if (tryCatch(job$is_alive(), error = function(e) FALSE)) return()
    zonation_job(NULL)
    result <- tryCatch(job$get_result(), error = function(e) e)
    if (inherits(result, "error")) {
      zonation_state("error"); zonation_message(conditionMessage(result))
    } else {
      zonation_manifest(result); zonation_state("success")
      zonation_message("Continuous CV-to-PV zonation analysis is complete.")
      trajectory_message("Enter a UniProt accession or gene name.")
      files <- result$files; ids <- as.character(seq_len(nrow(files)))
      updateSelectInput(session, "zonation_file_choice",
                        choices = stats::setNames(ids, files$name))
      plot_rows <- which(files$kind == "plot")
      table_rows <- which(files$kind == "table")
      updateSelectInput(session, "zonation_plot_choice",
                        choices = stats::setNames(ids[plot_rows], files$name[plot_rows]))
      updateSelectInput(session, "zonation_table_choice",
                        choices = stats::setNames(ids[table_rows], files$name[table_rows]))
    }
  })

  observeEvent(input$run_trajectory, {
    if (!is.null(trajectory_job())) {
      showNotification("Trajectory analysis is running.", type = "message")
      return()
    }
    if (!identical(zonation_state(), "success")) {
      showNotification("Complete Continuous analysis successfully first.",
                       type = "error", duration = 8)
      return()
    }
    targets <- tryCatch(parse_trajectory_targets(input$trajectory_target),
                        error = function(e) e)
    if (inherits(targets, "error")) {
      showNotification(conditionMessage(targets), type = "error", duration = 8)
      return()
    }
    if (!length(targets)) {
      showNotification("Enter at least one UniProt accession or gene name.",
                       type = "error")
      return()
    }
    paths <- continuous_run_paths()
    if (is.null(paths) || is.null(paths$zonation_dir) ||
        !dir.exists(paths$zonation_dir)) {
      showNotification("Cannot find the Continuous analysis results directory.",
                       type = "error")
      return()
    }
    required_inputs <- file.path(paths$zonation_dir, c(
      "GF_vs_SPF_Zonation_Interaction.csv",
      "Spatial_Proteomics_AnalysisReady.csv",
      "AllSamples_Zonation_Coefficients.csv"
    ))
    if (any(!file.exists(required_inputs))) {
      showNotification(
        paste("Required Continuous result files are missing:",
              paste(basename(required_inputs[!file.exists(required_inputs)]),
                    collapse = ", ")),
        type = "error", duration = 12
      )
      return()
    }
    prepared <- tryCatch({
      settings <- trajectory_script_settings()
      if (!is.null(settings$error)) stop(settings$error, call. = FALSE)
      # Keep targeted outputs directly below the Comparison run.  Nesting them
      # below Continuous can exceed the Windows MAX_PATH limit once the target
      # name and figure suffix are appended.
      dirs <- isolate(trajectory_run_paths())
      reusable <- !is.null(dirs) && !is.null(dirs$root) &&
        dir.exists(dirs$root) && !is.null(dirs$beta_dir) &&
        !is.null(dirs$curve_dir) && dir.exists(dirs$beta_dir) &&
        dir.exists(dirs$curve_dir) &&
        !length(list.files(dirs$beta_dir, all.files = FALSE)) &&
        !length(list.files(dirs$curve_dir, all.files = FALSE))
      if (!reusable) dirs <- create_trajectory_run_dir(dirname(paths$root))
      cache_dir <- resolve_user_path(input$trajectory_cache_path,
                                     "Trajectory cache directory")
      c(dirs, list(
        input_dir = normalizePath(paths$zonation_dir, winslash = "/", mustWork = TRUE),
        beta_script = settings$beta_path,
        curve_script = settings$curve_path,
        targets = targets,
        cache_dir = cache_dir
      ))
    }, error = function(e) e)
    if (inherits(prepared, "error")) {
      showNotification(conditionMessage(prepared), type = "error", duration = 10)
      return()
    }

    reset_trajectory("Running Trajectory analysis.")
    trajectory_run_paths(prepared)
    trajectory_state("running")
    trajectory_message(paste0(length(targets), " target plots are being generated."))
    result <- tryCatch(callr::r_bg(
      func = trajectory_analysis_worker,
      args = list(
        input_dir = prepared$input_dir,
        beta_output_dir = prepared$beta_dir,
        curve_output_dir = prepared$curve_dir,
        beta_script_path = prepared$beta_script,
        curve_script_path = prepared$curve_script,
        targets = targets,
        cache_dir = prepared$cache_dir
      ),
      wd = prepared$root,
      stdout = file.path(prepared$root, "Trajectory_execution.log"),
      stderr = "2>&1", supervise = TRUE,
      user_profile = FALSE, system_profile = FALSE
    ), error = function(e) e)
    if (inherits(result, "error")) {
      trajectory_state("error"); trajectory_message(conditionMessage(result))
    } else {
      trajectory_job(result)
    }
  })

  observeEvent(input$cancel_trajectory, {
    if (is.null(trajectory_job())) return()
    kill_process(trajectory_job()); trajectory_job(NULL)
    trajectory_state("cancelled")
    trajectory_message("Trajectory analysis was stopped by the user.")
  })

  observe({
    job <- trajectory_job()
    if (is.null(job)) return()
    invalidateLater(750, session)
    trajectory_log_tick(isolate(trajectory_log_tick()) + 1L)
    if (tryCatch(job$is_alive(), error = function(e) FALSE)) return()
    trajectory_job(NULL)
    result <- tryCatch(job$get_result(), error = function(e) e)
    if (inherits(result, "error")) {
      trajectory_state("error"); trajectory_message(conditionMessage(result))
    } else {
      trajectory_manifest(result); trajectory_state("success")
      both_created <- sum(result$summary$Beta_Status != "Failed" &
                            result$summary$Trajectory_Status != "Failed")
      cached <- sum(result$summary$Beta_Status == "Cached" &
                      result$summary$Trajectory_Status == "Cached")
      trajectory_message(paste0(
        length(result$targets), " targets processed · both plots generated for ",
        both_created, " · fully cached ", cached, " targets"
      ))
      files <- result$files
      ids <- as.character(seq_len(nrow(files)))
      updateSelectInput(session, "trajectory_file_choice",
                        choices = stats::setNames(ids, files$name))
      beta_ids <- as.character(seq_len(nrow(result$beta_plots)))
      curve_ids <- as.character(seq_len(nrow(result$curve_plots)))
      updateSelectInput(session, "trajectory_beta_choice",
                        choices = stats::setNames(beta_ids, result$beta_plots$Target),
                        selected = if (length(beta_ids)) beta_ids[[1]] else character())
      updateSelectInput(session, "trajectory_curve_choice",
                        choices = stats::setNames(curve_ids, result$curve_plots$Target),
                        selected = if (length(curve_ids)) curve_ids[[1]] else character())
    }
  })

  observeEvent(input$trajectory_beta_choice, {
    manifest <- trajectory_manifest()
    beta_id <- suppressWarnings(as.integer(input$trajectory_beta_choice))
    if (is.null(manifest) || length(beta_id) != 1L || is.na(beta_id) ||
        beta_id < 1L || beta_id > nrow(manifest$beta_plots)) return()
    selected_target <- manifest$beta_plots$Target[[beta_id]]
    curve_id <- match(selected_target, manifest$curve_plots$Target)
    if (!is.na(curve_id)) {
      updateSelectInput(session, "trajectory_curve_choice",
                        selected = as.character(curve_id))
    }
  }, ignoreInit = TRUE)

  observeEvent(input$run_spatial_map, {
    if (!is.null(spatial_map_job())) {
      showNotification("Spatial Map analysis is running.", type = "message")
      return()
    }
    report <- source_info()
    if (is.null(report) || is.null(report$path) || !file.exists(report$path)) {
      showNotification("Connect report.parquet under Data Input.",
                       type = "error", duration = 8)
      return()
    }
    targets <- tryCatch(parse_trajectory_targets(input$spatial_map_target,
                                                   max_targets = 10L),
                        error = function(e) e)
    if (inherits(targets, "error") || !length(targets)) {
      showNotification(
        if (inherits(targets, "error")) conditionMessage(targets) else
          "Enter at least one UniProt accession or gene name.",
        type = "error", duration = 8
      )
      return()
    }
    prepared <- tryCatch({
      settings <- spatial_expression_script_settings()
      if (!is.null(settings$error)) stop(settings$error, call. = FALSE)
      histology_root <- resolve_user_path(input$spatial_histology_root,
                                          "Histology image directory")
      if (!dir.exists(histology_root)) {
        stop("Cannot find the Histology image directory: ", histology_root,
             call. = FALSE)
      }
      continuous_paths <- continuous_run_paths()
      if (is.null(continuous_paths) || !dir.exists(continuous_paths$root)) {
        stop("Cannot find the ROI metadata run path.", call. = FALSE)
      }
      coordinate_path <- continuous_paths$coordinate_path
      if (is.null(coordinate_path) || !file.exists(coordinate_path)) {
        stop(paste0(
          "Cannot find Image_Reference_Coordinates.csv used by 01 Generation ROI metadata. ",
          "Regenerate ROI metadata."
        ), call. = FALSE)
      }
      metadata_path <- roi_metadata_file()
      if (is.null(metadata_path) || !file.exists(metadata_path)) {
        stop(
          "Cannot find ROI_metadata.csv generated by 01 Generation ROI metadata. Regenerate ROI metadata.",
          call. = FALSE
        )
      }
      # Spatial Map has its own output hierarchy and does not create or reuse
      # a 01_Trajectory run. ROI metadata provides the authoritative CV-PV
      # Spatial_S. If Continuous zonation output exists, its analysis-ready
      # values are reused so the curves exactly match 01 Trajectory.
      run_dir <- create_spatial_map_run_dir(continuous_paths$root)
      staged_coordinate <- file.path(run_dir, "Image_Reference_coordinates.csv")
      if (!file.copy(coordinate_path, staged_coordinate, overwrite = TRUE)) {
        stop("Failed to copy the coordinate CSV to the Spatial Map run directory.",
             call. = FALSE)
      }
      staged_metadata <- file.path(run_dir, "ROI_metadata.csv")
      if (!file.copy(metadata_path, staged_metadata, overwrite = TRUE)) {
        stop("Failed to copy ROI_metadata.csv to the Spatial Map run directory.",
             call. = FALSE)
      }
      analysis_ready_source <- file.path(
        continuous_paths$zonation_dir,
        "Spatial_Proteomics_AnalysisReady.csv"
      )
      staged_analysis_ready <- NULL
      if (file.exists(analysis_ready_source)) {
        staged_analysis_ready <- file.path(
          run_dir, "Spatial_Proteomics_AnalysisReady.csv"
        )
        if (!file.copy(analysis_ready_source, staged_analysis_ready,
                       overwrite = TRUE)) {
          stop("Failed to copy Continuous analysis-ready data to the Spatial Map run directory.",
               call. = FALSE)
        }
      }
      list(
        root = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
        report_path = normalizePath(report$path, winslash = "/", mustWork = TRUE),
        script_path = settings$path,
        targets = targets,
        coordinate_path = normalizePath(staged_coordinate, winslash = "/",
                                        mustWork = TRUE),
        coordinate_name = basename(coordinate_path),
        roi_metadata_path = normalizePath(staged_metadata, winslash = "/",
                                          mustWork = TRUE),
        analysis_ready_path = if (is.null(staged_analysis_ready)) NULL else
          normalizePath(staged_analysis_ready, winslash = "/", mustWork = TRUE),
        histology_root = normalizePath(histology_root, winslash = "/",
                                       mustWork = TRUE),
        require_all_sets = isTRUE(input$spatial_require_all_sets)
      )
    }, error = function(e) e)
    if (inherits(prepared, "error")) {
      showNotification(conditionMessage(prepared), type = "error", duration = 12)
      return()
    }

    reset_spatial_map("Generating Spatial Maps.")
    spatial_map_run_paths(prepared)
    spatial_map_state("running")
    spatial_map_message(paste0(
      length(prepared$targets), " target Spatial Maps are being generated: ",
      paste(prepared$targets, collapse = ", ")
    ))
    result <- tryCatch(callr::r_bg(
      func = spatial_map_worker,
      args = list(
        report_path = prepared$report_path,
        output_root = prepared$root,
        script_path = prepared$script_path,
        targets = prepared$targets,
        coordinate_path = prepared$coordinate_path,
        roi_metadata_path = prepared$roi_metadata_path,
        analysis_ready_path = prepared$analysis_ready_path,
        histology_root = prepared$histology_root,
        require_all_sets = prepared$require_all_sets
      ),
      wd = prepared$root,
      stdout = file.path(prepared$root, "Spatial_Map_execution.log"),
      stderr = "2>&1", supervise = TRUE,
      user_profile = FALSE, system_profile = FALSE
    ), error = function(e) e)
    if (inherits(result, "error")) {
      spatial_map_state("error"); spatial_map_message(conditionMessage(result))
    } else {
      spatial_map_job(result)
    }
  })

  observeEvent(input$cancel_spatial_map, {
    if (is.null(spatial_map_job())) return()
    kill_process(spatial_map_job()); spatial_map_job(NULL)
    spatial_map_state("cancelled")
    spatial_map_message("Spatial Map was stopped by the user.")
  })

  observe({
    job <- spatial_map_job()
    if (is.null(job)) return()
    invalidateLater(750, session)
    spatial_map_log_tick(isolate(spatial_map_log_tick()) + 1L)
    if (tryCatch(job$is_alive(), error = function(e) FALSE)) return()
    spatial_map_job(NULL)
    result <- tryCatch(job$get_result(), error = function(e) e)
    if (inherits(result, "error")) {
      spatial_map_state("error"); spatial_map_message(conditionMessage(result))
    } else {
      spatial_map_manifest(result); spatial_map_state("success")
      spatial_map_message(paste0(
        length(result$successful_targets), "/", length(result$targets),
        " targets complete · images ", nrow(result$images),
        " · total files ", nrow(result$files), ""
      ))
      file_ids <- as.character(seq_len(nrow(result$files)))
      updateSelectInput(session, "spatial_map_file_choice",
                        choices = stats::setNames(file_ids, result$files$name))
      image_priority <- ifelse(
        grepl("_(GF|SPF)_Trajectory_Combined[.]png$", result$images$name), 1L,
        ifelse(grepl("_Spatial_MAP[.]png$", result$images$name), 2L,
               ifelse(grepl("_Trajectory_SpatialMap_Combined[.]png$",
                            result$images$name), 3L,
               ifelse(grepl("_Spatial_Map_Publication[.]png$",
                            result$images$name), 4L, 5L)))
      )
      image_target_order <- match(result$images$target,
                                  result$successful_targets)
      image_order <- order(image_target_order, image_priority,
                           result$images$name, na.last = TRUE)
      result$images <- result$images[image_order, , drop = FALSE]
      spatial_map_manifest(result)
      sample_choices <- spatial_map_sample_choices(result)
      updateCheckboxGroupInput(
        session, "spatial_map_target_choice",
        choices = stats::setNames(result$successful_targets,
                                  result$successful_targets),
        selected = result$successful_targets, inline = TRUE
      )
      updateCheckboxGroupInput(
        session, "spatial_map_image_choice",
        choices = stats::setNames(sample_choices, sample_choices),
        selected = sample_choices,
        inline = TRUE
      )
      register_spatial_map_image_outputs(
        as.character(seq_len(nrow(result$images)))
      )
    }
  })

  update_spatial_map_image_selection <- function(group = c("all", "gf", "spf", "none")) {
    group <- match.arg(group)
    manifest <- spatial_map_manifest()
    samples <- spatial_map_sample_choices(manifest)
    if (!length(samples)) return(invisible(NULL))
    selected <- switch(
      group,
      all = samples,
      gf = samples[grepl("^GF", samples, ignore.case = TRUE)],
      spf = samples[grepl("^SPF", samples, ignore.case = TRUE)],
      none = character()
    )
    updateCheckboxGroupInput(
      session, "spatial_map_image_choice", selected = selected,
      inline = TRUE
    )
  }
  observeEvent(input$spatial_map_show_combined, {
    manifest <- spatial_map_manifest()
    if (is.null(manifest) || is.null(manifest$images) || !nrow(manifest$images)) {
      return()
    }
    samples <- spatial_map_sample_choices(manifest)
    updateCheckboxGroupInput(
      session,
      "spatial_map_image_choice",
      choices = stats::setNames(samples, samples),
      selected = samples,
      inline = TRUE
    )
  }, ignoreInit = TRUE)
  observeEvent(input$spatial_map_select_all,
               update_spatial_map_image_selection("all"))
  observeEvent(input$spatial_map_select_gf,
               update_spatial_map_image_selection("gf"))
  observeEvent(input$spatial_map_select_spf,
               update_spatial_map_image_selection("spf"))
  observeEvent(input$spatial_map_clear_selection,
               update_spatial_map_image_selection("none"))

  observeEvent(input$spatial_map_target_select_all, {
    manifest <- spatial_map_manifest()
    if (is.null(manifest) || !length(manifest$successful_targets)) return()
    updateCheckboxGroupInput(session, "spatial_map_target_choice",
                             selected = manifest$successful_targets,
                             inline = TRUE)
  })
  observeEvent(input$spatial_map_target_select_first, {
    manifest <- spatial_map_manifest()
    if (is.null(manifest) || !length(manifest$successful_targets)) return()
    updateCheckboxGroupInput(session, "spatial_map_target_choice",
                             selected = manifest$successful_targets[[1]],
                             inline = TRUE)
  })
  observeEvent(input$spatial_map_target_clear, {
    updateCheckboxGroupInput(session, "spatial_map_target_choice",
                             selected = character(), inline = TRUE)
  })

  observeEvent(list(input$calculate_summary, input$refresh_summary), {
    req(source_info())
    if (qc_state() == "running") { showNotification("Calculate the summary after QC is complete."); return() }
    summary_result(NULL); run_result(NULL); summary_error(NULL)
    result <- tryCatch(withProgress(message = "Calculating full-file summary...", value = 0, {
      totals <- query_diann_summary(con)
      incProgress(0.4, detail = "Calculating statistics by Run")
      runs <- query_diann_runs(con)
      incProgress(0.6)
      list(totals = totals, runs = runs)
    }), error = function(e) e)
    if (inherits(result, "error")) {
      summary_error(conditionMessage(result))
      showNotification("Summary calculation failed. Check the error message.",
                       type = "error", duration = 8)
    } else {
      summary_result(result$totals); run_result(result$runs)
      updateTabItems(session, "main_tabs", selected = "dashboard")
    }
  }, ignoreInit = TRUE)

  observeEvent(input$query_preview, {
    info <- source_info()
    req(info)
    explorer_result(NULL)
    columns <- intersect(input$explorer_columns, info$columns)
    if (!length(columns)) {
      explorer_message("Select at least one column to display.")
      return()
    }
    limit <- suppressWarnings(as.integer(input$explorer_limit))
    if (length(limit) != 1L || is.na(limit)) limit <- 100L
    limit <- max(1L, min(limit, APP_CONFIG$max_explorer_rows))
    sql <- build_explorer_sql(con, columns,
                               input$run_filter, input$protein_filter, limit)
    result <- tryCatch(withProgress(message = "Querying data that match the criteria...", value = 0, {
      value <- DBI::dbGetQuery(con, sql)
      incProgress(1)
      value
    }), error = function(e) e)
    if (inherits(result, "error")) {
      explorer_message(paste("Query failed:", conditionMessage(result)))
    } else {
      explorer_result(result)
      explorer_message(sprintf("Query result: %s rows / display limit: %s rows. This is not the total number of matching rows.",
                               format_count(nrow(result)), format_count(limit)))
    }
  })


  output$input_metrics <- renderUI({
    pair <- source_pair()
    connected <- sum(c(!is.null(pair$first_pass), !is.null(pair$report)))
    total <- sum(vapply(pair, function(x) if (is.null(x)) 0 else x$size, numeric(1)))
    state_labels <- c(idle = "Ready", running = "Running", success = "Completed",
                      error = "Error", cancelled = "Stopped")
    fluidRow(
      valueBox(paste0(connected, " / 2"), "Input files", icon("file"), color = "teal", width = 3),
      valueBox(sprintf("%.2f GiB", total / 1024^3), "Total input size",
               icon("database"), color = "aqua", width = 3),
      valueBox(unname(state_labels[[qc_state()]]), "QC status",
               icon("clipboard-check"), color = "navy", width = 3),
      valueBox(if (is.null(qc_manifest())) "\u2014" else nrow(qc_manifest()$files),
               "Result files", icon("folder"), color = "light-blue", width = 3)
    )
  })
  output$dashboard_metrics <- renderUI(metric_cards(source_info(), summary_result()))
  output$sidebar_file <- renderUI({
    pair <- source_pair()
    continuous_paths <- continuous_run_paths()
    coordinate_connected <- !is.null(continuous_paths) &&
      !is.null(continuous_paths$coordinate_path) &&
      file.exists(continuous_paths$coordinate_path)
    div(class = "sidebar-file", span(class = "micro-label", "PAIRED INPUTS"),
        span(class = "sidebar-file-name",
             if (is.null(pair$first_pass)) "1. First-pass: not connected" else "1. First-pass: connected"),
        span(class = "sidebar-file-name",
             if (is.null(pair$report)) "2. Report: not connected" else "2. Report: connected"),
        span(class = "sidebar-file-name",
             if (coordinate_connected) {
               "3. ROI coordinates: connected"
             } else {
               "3. ROI coordinates: not connected"
             }))
  })
  upload_status <- function(key) {
    info <- source_pair()[[key]]; error <- upload_errors()[[key]]
    if (!is.null(error)) {
      return(div(class = "status-banner status-error", icon("circle-exclamation"),
                 div(strong("File validation failed"), span(error))))
    }
    if (is.null(info)) {
      return(div(class = "status-banner", icon("circle-info"),
                 div(strong("No file selected"), span("Select the Parquet file with the required name."))))
    }
    div(class = "status-banner status-ready", icon("circle-check"),
        div(strong("File connected"),
            span(sprintf("%.3f GiB / %d columns", info$size / 1024^3, length(info$columns)))))
  }
  output$first_pass_status <- renderUI(upload_status("first_pass"))
  output$report_status <- renderUI(upload_status("report"))


  output$qc_path_status <- renderUI({
    settings <- qc_path_settings()
    if (!is.null(settings$error)) {
      return(div(class = "status-banner status-error", icon("circle-exclamation"),
                 div(strong("Path settings require attention"), span(settings$error))))
    }
    tagList(
      tags$dl(class = "file-meta",
        tags$dt("QC script"), tags$dd(settings$script_path),
        tags$dt("Results directory"), tags$dd(settings$results_root)),
      p(class = "muted",
        "If the storage directory does not exist, it is created when QC runs. The first QC run creates an Analysis_date_time session and four stage directories. Path changes made during a run apply to the next run.")
    )
  })
  output$qc_run_location <- renderUI({
    paths <- qc_run_paths()
    req(paths)
    tags$dl(class = "file-meta",
      tags$dt("Script for this run"), tags$dd(paths$script_path),
      tags$dt("Results directory for this run"), tags$dd(paths$run_dir),
      tags$dt("Manual exclude"), tags$dd(paths$manual_exclusions))
  })

  output$manual_exclusion_status <- renderUI({
    settings <- manual_exclusion_settings()
    if (!is.null(settings$error)) {
      return(div(class = "status-banner status-error", icon("circle-exclamation"),
                 div(strong("Manual exclude input requires attention"), span(settings$error))))
    }
    exclusions <- settings$value
    if (!length(exclusions)) {
      return(p(class = "muted", "Manual exclude: None"))
    }
    div(class = "status-banner status-ready", icon("circle-check"),
        div(strong("Check Manual exclude input"),
            span(format_manual_exclusions(exclusions))))
  })

  output$qc_readiness <- renderUI({
    pair <- source_pair()
    both <- !is.null(pair$first_pass) && !is.null(pair$report)
    settings <- qc_path_settings()
    script_ok <- is.null(settings$error)
    manual_ok <- is.null(manual_exclusion_settings()$error)
    text <- if (!both) "Connect both Parquet files." else if (!script_ok) {
      "Both input files are ready. Check the QC path settings above."
    } else if (!manual_ok) {
      "Both input files and QC paths are ready. Check the Manual exclude input format."
    } else {
      "Both input files, QC paths, and Manual exclude settings are ready. Write permission for the results directory is checked when QC runs."
    }
    ready <- both && script_ok && manual_ok
    div(class = paste("status-banner", if (ready) "status-ready" else ""),
        icon(if (ready) "circle-check" else "circle-info"),
        div(strong("QC readiness"), span(text)))
  })
  output$qc_result_status <- renderUI({
    state <- qc_state()
    css_state <- if (state == "success") "status-ready" else if (state == "error") "status-error" else ""
    div(class = paste("status-banner", css_state), icon("clipboard-check"),
        div(strong(paste("QC status:", state)), span(qc_message())))
  })
  output$comparison_script_status <- renderUI({
    settings <- comparison_script_settings()
    if (!is.null(settings$error)) {
      return(div(class = "status-banner status-error", icon("circle-exclamation"),
                 div(strong("Comparison analysis script path requires attention"), span(settings$error))))
    }
    tags$dl(class = "file-meta",
            tags$dt("Comparison analysis script"), tags$dd(settings$path))
  })
  output$comparison_status <- renderUI({
    state <- comparison_state()
    css_state <- if (state == "success") "status-ready" else if (state == "error") "status-error" else ""
    div(class = paste("status-banner", css_state), icon("scale-balanced"),
        div(strong(paste("Comparison status:", state)), span(comparison_message())))
  })
  output$comparison_readiness <- renderUI({
    groups <- comparison_groups()
    gf_selected <- unique(as.character(input$comparison_gf_sets))
    spf_selected <- unique(as.character(input$comparison_spf_sets))
    script_ok <- is.null(comparison_script_settings()$error)
    report_ok <- !is.null(source_info())
    qc_ok <- identical(qc_state(), "success")
    qc_result <- qc_manifest()
    excluded_count <- if (!is.null(qc_result) &&
                             is.data.frame(qc_result$excluded_rois)) {
      nrow(qc_result$excluded_rois)
    } else 0L
    selection_ok <- length(gf_selected) >= 3L && length(spf_selected) >= 3L
    ready <- report_ok && qc_ok && script_ok && selection_ok
    text <- if (!report_ok) {
      "Connect report.parquet to display the GF/SPF group checklist."
    } else if (!length(groups$GF) || !length(groups$SPF)) {
      "Both GF and SPF groups could not be detected from report.parquet Run names."
    } else if (!qc_ok) {
      "Complete QC successfully first. Groups can be selected now."
    } else if (!script_ok) {
      "Check the comparison analysis script path."
    } else if (!selection_ok) {
      "Select at least three GF and three SPF groups."
    } else {
      paste0("Ready to run · GF ", length(gf_selected), " / SPF ",
             length(spf_selected), " selected · QC excluded ROIs ",
             excluded_count, " applied")
    }
    div(class = paste("status-banner", if (ready) "status-ready" else ""),
        icon(if (ready) "circle-check" else "circle-info"),
        div(strong("Comparison readiness"), span(text)))
  })
  output$comparison_run_location <- renderUI({
    paths <- comparison_run_paths(); req(paths)
    tags$dl(class = "file-meta",
      tags$dt("Script for this run"), tags$dd(paths$script_path),
      tags$dt("Results directory for this run"), tags$dd(paths$run_dir),
      tags$dt("Selected GF"), tags$dd(paths$gf_sets),
      tags$dt("Selected SPF"), tags$dd(paths$spf_sets),
      tags$dt("QC excluded ROIs"), tags$dd(paste0(paths$excluded_count, " ROIs")))
  })
  output$continuous_overall_status <- renderUI({
    metadata_state <- roi_metadata_state()
    zone_state <- zonation_state()
    css_state <- if (zone_state == "success") "status-ready" else if (
      metadata_state == "error" || zone_state == "error") "status-error" else ""
    div(class = paste("status-banner", css_state), icon("wave-square"),
      div(
        strong(paste0("ROI metadata: ", metadata_state,
                      " · Continuous zonation: ", zone_state)),
        span(if (zone_state != "idle") zonation_message() else roi_metadata_message())
      ))
  })
  output$roi_metadata_readiness <- renderUI({
    coordinate <- input$roi_coordinate_file
    coordinate_ok <- !is.null(coordinate) && file.exists(coordinate$datapath) &&
      identical(tolower(tools::file_ext(coordinate$name)), "csv")
    comparison_ok <- identical(comparison_state(), "success") &&
      !is.null(comparison_output_dir())
    script_ok <- is.null(roi_metadata_script_settings()$error)
    ready <- coordinate_ok && comparison_ok && script_ok &&
      is.null(roi_metadata_job()) && is.null(zonation_job())
    text <- if (!comparison_ok) {
      "Complete Comparison analysis successfully first."
    } else if (!coordinate_ok) {
      "Select Image_Reference_Coordinates.csv."
    } else if (!script_ok) {
      roi_metadata_script_settings()$error
    } else {
      "Ready to generate ROI metadata."
    }
    div(class = paste("status-banner", if (ready) "status-ready" else ""),
        icon(if (ready) "circle-check" else "circle-info"),
        div(strong("Step 1 readiness"), span(text)))
  })
  output$roi_metadata_run_location <- renderUI({
    paths <- continuous_run_paths(); req(paths)
    tags$dl(class = "file-meta",
      tags$dt("Reference-coordinate file"), tags$dd(paths$coordinate_name),
      tags$dt("Metadata script"), tags$dd(paths$metadata_script),
      tags$dt("Metadata results"), tags$dd(paths$metadata_dir),
      if (!is.null(roi_metadata_file())) tagList(
        tags$dt("Generated metadata"), tags$dd(roi_metadata_file()))
    )
  })
  output$zonation_readiness <- renderUI({
    comparison_ok <- identical(comparison_state(), "success")
    metadata_ok <- identical(roi_metadata_state(), "success") &&
      !is.null(roi_metadata_file()) && file.exists(roi_metadata_file())
    report_ok <- !is.null(source_info())
    script_ok <- is.null(zonation_script_settings()$error)
    ready <- comparison_ok && metadata_ok && report_ok && script_ok &&
      is.null(roi_metadata_job()) && is.null(zonation_job())
    text <- if (!comparison_ok) {
      "Complete comparison analysis first."
    } else if (!metadata_ok) {
      "Generate ROI metadata in Step 1 first."
    } else if (!report_ok) {
      "Connect report.parquet."
    } else if (!script_ok) {
      zonation_script_settings()$error
    } else {
      paste("Ready to run · metadata:", basename(roi_metadata_file()))
    }
    div(class = paste("status-banner", if (ready) "status-ready" else ""),
        icon(if (ready) "circle-check" else "circle-info"),
        div(strong("Step 2 readiness"), span(text)))
  })
  output$zonation_run_location <- renderUI({
    paths <- continuous_run_paths(); req(paths)
    tags$dl(class = "file-meta",
      tags$dt("Zonation script"), tags$dd(if (is.null(paths$zonation_script))
        input$continuous_zonation_script_path else paths$zonation_script),
      tags$dt("Metadata used"), tags$dd(if (is.null(roi_metadata_file()))
        "Not generated yet" else roi_metadata_file()),
      tags$dt("Zonation results"), tags$dd(paths$zonation_dir))
  })
  output$trajectory_status <- renderUI({
    state <- trajectory_state()
    css_state <- if (state == "success") "status-ready" else if (state == "error") {
      "status-error"
    } else ""
    div(class = paste("status-banner", css_state), icon("chart-line"),
        div(strong(paste("Trajectory analysis:", state)),
            span(trajectory_message())))
  })
  output$trajectory_readiness <- renderUI({
    targets <- tryCatch(parse_trajectory_targets(input$trajectory_target),
                        error = function(e) e)
    target_ok <- !inherits(targets, "error") && length(targets) > 0L
    continuous_ok <- identical(zonation_state(), "success")
    settings <- trajectory_script_settings()
    script_ok <- is.null(settings$error)
    paths <- continuous_run_paths()
    required <- if (!is.null(paths) && !is.null(paths$zonation_dir)) {
      file.path(paths$zonation_dir, c(
        "GF_vs_SPF_Zonation_Interaction.csv",
        "Spatial_Proteomics_AnalysisReady.csv",
        "AllSamples_Zonation_Coefficients.csv"
      ))
    } else character()
    input_ok <- length(required) == 3L && all(file.exists(required))
    ready <- continuous_ok && input_ok && script_ok && target_ok &&
      is.null(trajectory_job())
    text <- if (!continuous_ok) {
      "Complete Continuous analysis successfully first."
    } else if (!input_ok) {
      paste0("The three required CSV files cannot be found in the Continuous results: ",
             "GF_vs_SPF_Zonation_Interaction.csv, Spatial_Proteomics_AnalysisReady.csv, ",
             "AllSamples_Zonation_Coefficients.csv")
    } else if (!script_ok) {
      settings$error
    } else if (!target_ok) {
      if (inherits(targets, "error")) conditionMessage(targets) else
        "Enter at least one UniProt accession or gene name."
    } else {
      paste0("Ready to run · targets ", length(targets), ": ",
             paste(targets, collapse = ", "))
    }
    div(class = paste("status-banner", if (ready) "status-ready" else ""),
        icon(if (ready) "circle-check" else "circle-info"),
        div(strong("Trajectory readiness"), span(text)))
  })
  output$trajectory_run_location <- renderUI({
    paths <- trajectory_run_paths(); req(paths)
    tags$dl(class = "file-meta",
      tags$dt("Input Continuous results"), tags$dd(paths$input_dir),
      tags$dt("Targets"), tags$dd(paste(paths$targets, collapse = ", ")),
      tags$dt("Cache directory"), tags$dd(paths$cache_dir),
      tags$dt("Beta script"), tags$dd(paths$beta_script),
      tags$dt("Trajectory script"), tags$dd(paths$curve_script),
      tags$dt("Results directory"), tags$dd(paths$root))
  })
  output$spatial_map_status <- renderUI({
    state <- spatial_map_state()
    css_state <- if (state == "success") "status-ready" else if (state == "error") {
      "status-error"
    } else ""
    div(class = paste("status-banner", css_state), icon("map"),
        div(strong(paste("Spatial Map:", state)),
            span(spatial_map_message())))
  })
  output$spatial_map_readiness <- renderUI({
    targets <- tryCatch(parse_trajectory_targets(input$spatial_map_target,
                                                   max_targets = 10L),
                        error = function(e) e)
    target_ok <- !inherits(targets, "error") && length(targets) > 0L
    report_ok <- !is.null(source_info()) && !is.null(source_info()$path) &&
      file.exists(source_info()$path)
    continuous_paths <- continuous_run_paths()
    coordinate_ok <- !is.null(continuous_paths) &&
      !is.null(continuous_paths$coordinate_path) &&
      file.exists(continuous_paths$coordinate_path)
    settings <- spatial_expression_script_settings()
    script_ok <- is.null(settings$error)
    histology <- tryCatch(
      resolve_user_path(input$spatial_histology_root, "Histology image directory"),
      error = function(e) e
    )
    histology_ok <- !inherits(histology, "error") && dir.exists(histology)
    ready <- target_ok && report_ok && coordinate_ok && script_ok &&
      histology_ok && is.null(spatial_map_job())
    text <- if (!report_ok) {
      "Connect report.parquet under Data Input."
    } else if (!script_ok) {
      settings$error
    } else if (!target_ok) {
      if (inherits(targets, "error")) conditionMessage(targets) else
        "Enter at least one UniProt accession or gene name."
    } else if (!coordinate_ok) {
      "Cannot find the ROI coordinate CSV. Run 01 Generation ROI metadata first."
    } else if (!histology_ok) {
      if (inherits(histology, "error")) conditionMessage(histology) else
        "Cannot find the Histology image directory."
    } else {
      paste0("Ready to run · targets ", length(targets), ": ",
             paste(targets, collapse = ", "),
             " · histology: ", basename(histology))
    }
    div(class = paste("status-banner", if (ready) "status-ready" else ""),
        icon(if (ready) "circle-check" else "circle-info"),
        div(strong("Spatial Map readiness"), span(text)))
  })
  output$spatial_map_run_location <- renderUI({
    paths <- spatial_map_run_paths(); req(paths)
    tags$dl(class = "file-meta",
      tags$dt("Targets"), tags$dd(paste(paths$targets, collapse = ", ")),
      tags$dt("report.parquet"), tags$dd(paths$report_path),
      tags$dt("Coordinate CSV"), tags$dd(paths$coordinate_name),
      tags$dt("Histology directory"), tags$dd(paths$histology_root),
      tags$dt("Spatial script"), tags$dd(paths$script_path),
      tags$dt("Results directory"), tags$dd(paths$root))
  })
  output$summary_status <- renderUI({
    if (!is.null(summary_error())) {
      return(div(class = "status-banner status-error",
                 div(strong("Final report summary failed"), span(summary_error()))))
    }
    if (!is.null(summary_result())) p(class = "muted", "The final report summary has been calculated.")
  })
  output$dashboard_status <- renderUI({
    if (!is.null(summary_error())) return(p(class = "muted", summary_error()))
    if (is.null(source_info())) return(p(class = "muted", "Connect report.parquet under Data Input."))
    if (is.null(summary_result())) return(p(class = "muted", "Calculate the final report summary."))
    p(class = "muted", icon("circle-check"), " report.parquet summary complete")
  })
  selected_input <- reactive({
    pair <- source_pair()
    if (identical(input$preview_source, "report")) pair$report else pair$first_pass
  })
  output$preview_empty <- renderUI({
    if (is.null(selected_input())) div(class = "empty-state", icon("table"),
      strong("Input preview"), span("A preview appears after the selected input file is connected."))
  })
  output$preview <- DT::renderDT({
    info <- selected_input(); req(info)
    DT::datatable(info$preview, rownames = FALSE, escape = TRUE,
                  options = list(scrollX = TRUE, pageLength = 10, lengthChange = FALSE))
  }, server = TRUE)

  selected_qc_file <- function(id, kind = NULL) {
    manifest <- qc_manifest(); req(manifest, id)
    index <- suppressWarnings(as.integer(id))
    req(length(index) == 1L, !is.na(index), index >= 1L, index <= nrow(manifest$files))
    row <- manifest$files[index, , drop = FALSE]
    if (!is.null(kind)) req(row$kind[[1]] == kind)
    req(file.exists(row$path[[1]]))
    path <- normalizePath(row$path[[1]], winslash = "/", mustWork = TRUE)
    root <- normalizePath(qc_output_dir(), winslash = "/", mustWork = TRUE)
    req(startsWith(path, paste0(root, "/")))
    row
  }
  output$qc_dataset_summary <- DT::renderDT({
    manifest <- qc_manifest()
    req(manifest, !is.null(manifest$dataset_summary))
    DT::datatable(manifest$dataset_summary, rownames = FALSE, escape = TRUE,
                  options = list(scrollX = TRUE, pageLength = 10))
  }, server = TRUE)
  output$qc_image <- renderImage({
    row <- selected_qc_file(input$qc_plot_choice, "plot")
    extension <- tolower(tools::file_ext(row$path[[1]]))
    list(src = row$path[[1]],
         contentType = if (extension == "png") "image/png" else "image/jpeg",
         alt = row$name[[1]])
  }, deleteFile = FALSE)
  output$qc_plot_empty <- renderUI({
    manifest <- qc_manifest()
    if (is.null(manifest) || !any(manifest$files$kind == "plot")) {
      div(class = "empty-state", icon("chart-bar"), strong("QC figures"),
          span("PNG/JPG results appear here after QC is complete."))
    }
  })
  output$qc_table_empty <- renderUI({
    manifest <- qc_manifest()
    if (is.null(manifest) || !any(manifest$files$kind == "table")) {
      p(class = "muted", "CSV/TSV results appear here after QC is complete.")
    }
  })
  output$qc_file_table <- DT::renderDT({
    manifest <- qc_manifest(); req(manifest)
    files <- manifest$files
    visible <- data.frame(File = files$name, Type = files$kind,
                           KiB = round(files$size_bytes / 1024, 1))
    DT::datatable(visible, rownames = FALSE, escape = TRUE,
                  options = list(scrollX = TRUE, pageLength = 5, lengthChange = FALSE))
  }, server = TRUE)
  output$qc_table <- DT::renderDT({
    row <- selected_qc_file(input$qc_table_choice, "table")
    path <- row$path[[1]]
    separator <- if (tolower(tools::file_ext(path)) == "tsv") "\t" else ","
    sql <- sprintf(paste0(
      "SELECT * FROM read_csv(%s, delim = %s, header = true, all_varchar = true, ",
      "sample_size = 20480) LIMIT %d"
    ), as.character(DBI::dbQuoteString(con, path)),
       as.character(DBI::dbQuoteString(con, separator)), APP_CONFIG$qc_table_preview_rows)
    value <- tryCatch(DBI::dbGetQuery(con, sql), error = function(e) e)
    if (inherits(value, "error")) validate(need(FALSE, conditionMessage(value)))
    DT::datatable(value, rownames = FALSE, escape = TRUE,
                  options = list(scrollX = TRUE, pageLength = 10))
  }, server = TRUE)
  output$download_qc_file <- downloadHandler(
    filename = function() basename(selected_qc_file(input$qc_file_choice)$path[[1]]),
    content = function(file) {
      path <- selected_qc_file(input$qc_file_choice)$path[[1]]
      if (!file.copy(path, file, overwrite = TRUE)) stop("Failed to download the result file.")
    }
  )
  output$qc_log <- renderText({
    qc_log_tick(); qc_state()
    out <- qc_output_dir()
    if (is.null(out)) return("QC has not been run yet.")
    log <- file.path(out, "QC_execution.log")
    if (!file.exists(log)) return("Preparing the log file.")
    paste(readLines(log, n = 2000L, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  })
  output$download_qc_log <- downloadHandler(
    filename = function() "QC_execution.log",
    content = function(file) {
      req(qc_output_dir())
      log <- file.path(qc_output_dir(), "QC_execution.log")
      req(file.exists(log))
      if (!file.copy(log, file, overwrite = TRUE)) stop("Failed to download the execution log.")
    }
  )

  selected_comparison_file <- function(id, kind = NULL) {
    manifest <- comparison_manifest(); req(manifest, id)
    index <- suppressWarnings(as.integer(id))
    req(length(index) == 1L, !is.na(index), index >= 1L,
        index <= nrow(manifest$files))
    row <- manifest$files[index, , drop = FALSE]
    if (!is.null(kind)) req(row$kind[[1]] == kind)
    req(file.exists(row$path[[1]]))
    path <- normalizePath(row$path[[1]], winslash = "/", mustWork = TRUE)
    root <- normalizePath(comparison_output_dir(), winslash = "/", mustWork = TRUE)
    req(startsWith(path, paste0(root, "/")))
    row
  }
  output$comparison_summary <- DT::renderDT({
    manifest <- comparison_manifest()
    req(manifest, !is.null(manifest$summary))
    DT::datatable(manifest$summary, rownames = FALSE, escape = TRUE,
                  options = list(scrollX = TRUE, pageLength = 5, lengthChange = FALSE))
  }, server = TRUE)
  output$comparison_summary_empty <- renderUI({
    manifest <- comparison_manifest()
    if (is.null(manifest) || is.null(manifest$summary)) {
      p(class = "muted", "Selected groups and a summary of key results appear after comparison analysis is complete.")
    }
  })
  output$comparison_image <- renderImage({
    row <- selected_comparison_file(input$comparison_plot_choice, "plot")
    extension <- tolower(tools::file_ext(row$path[[1]]))
    list(src = row$path[[1]],
         contentType = if (extension == "png") "image/png" else "image/jpeg",
         alt = row$name[[1]])
  }, deleteFile = FALSE)
  output$comparison_plot_empty <- renderUI({
    manifest <- comparison_manifest()
    if (is.null(manifest) || !any(manifest$files$kind == "plot")) {
      div(class = "empty-state", icon("chart-line"), strong("Comparison figures"),
          span("PNG/JPG results appear here after comparison analysis is complete."))
    }
  })
  output$comparison_table_empty <- renderUI({
    manifest <- comparison_manifest()
    if (is.null(manifest) || !any(manifest$files$kind == "table")) {
      p(class = "muted", "CSV/TSV results appear here after comparison analysis is complete.")
    }
  })
  output$comparison_file_table <- DT::renderDT({
    manifest <- comparison_manifest(); req(manifest)
    files <- manifest$files
    visible <- data.frame(File = files$name, Type = files$kind,
                          KiB = round(files$size_bytes / 1024, 1))
    DT::datatable(visible, rownames = FALSE, escape = TRUE,
                  options = list(scrollX = TRUE, pageLength = 5, lengthChange = FALSE))
  }, server = TRUE)
  output$comparison_table <- DT::renderDT({
    row <- selected_comparison_file(input$comparison_table_choice, "table")
    path <- row$path[[1]]
    separator <- if (tolower(tools::file_ext(path)) == "tsv") "\t" else ","
    sql <- sprintf(paste0(
      "SELECT * FROM read_csv(%s, delim = %s, header = true, all_varchar = true, ",
      "sample_size = 20480) LIMIT %d"
    ), as.character(DBI::dbQuoteString(con, path)),
       as.character(DBI::dbQuoteString(con, separator)), APP_CONFIG$qc_table_preview_rows)
    value <- tryCatch(DBI::dbGetQuery(con, sql), error = function(e) e)
    if (inherits(value, "error")) validate(need(FALSE, conditionMessage(value)))
    DT::datatable(value, rownames = FALSE, escape = TRUE,
                  options = list(scrollX = TRUE, pageLength = 10))
  }, server = TRUE)
  output$download_comparison_file <- downloadHandler(
    filename = function() basename(selected_comparison_file(
      input$comparison_file_choice)$path[[1]]),
    content = function(file) {
      path <- selected_comparison_file(input$comparison_file_choice)$path[[1]]
      if (!file.copy(path, file, overwrite = TRUE)) {
        stop("Failed to download the comparison result file.")
      }
    }
  )
  output$comparison_log <- renderText({
    comparison_log_tick(); comparison_state()
    out <- comparison_output_dir()
    if (is.null(out)) return("Comparison analysis has not been run yet.")
    log <- file.path(out, "Comparison_execution.log")
    if (!file.exists(log)) return("Preparing the log file.")
    paste(readLines(log, n = 2000L, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  })
  output$download_comparison_log <- downloadHandler(
    filename = function() "Comparison_execution.log",
    content = function(file) {
      req(comparison_output_dir())
      log <- file.path(comparison_output_dir(), "Comparison_execution.log")
      req(file.exists(log))
      if (!file.copy(log, file, overwrite = TRUE)) {
        stop("Failed to download the comparison execution log.")
      }
    }
  )

  selected_analysis_file <- function(manifest, id, root, kind = NULL) {
    req(manifest, id, root)
    index <- suppressWarnings(as.integer(id))
    req(length(index) == 1L, !is.na(index), index >= 1L,
        index <= nrow(manifest$files))
    row <- manifest$files[index, , drop = FALSE]
    if (!is.null(kind)) req(row$kind[[1]] == kind)
    req(file.exists(row$path[[1]]))
    path <- normalizePath(row$path[[1]], winslash = "/", mustWork = TRUE)
    root <- normalizePath(root, winslash = "/", mustWork = TRUE)
    req(identical(path, root) || startsWith(path, paste0(root, "/")))
    row
  }
  selected_roi_metadata_file <- function(id, kind = NULL) {
    paths <- continuous_run_paths()
    req(paths)
    selected_analysis_file(roi_metadata_manifest(), id, paths$metadata_dir, kind)
  }
  selected_zonation_file <- function(id, kind = NULL) {
    paths <- continuous_run_paths()
    req(paths)
    selected_analysis_file(zonation_manifest(), id, paths$zonation_dir, kind)
  }
  selected_trajectory_file <- function(id, kind = NULL) {
    paths <- trajectory_run_paths()
    req(paths)
    selected_analysis_file(trajectory_manifest(), id, paths$root, kind)
  }
  selected_spatial_map_file <- function(id, kind = NULL) {
    paths <- spatial_map_run_paths()
    req(paths)
    selected_analysis_file(spatial_map_manifest(), id, paths$root, kind)
  }
  spatial_map_display_image_ids <- function(manifest = spatial_map_manifest()) {
    if (is.null(manifest) || is.null(manifest$images) || !nrow(manifest$images)) {
      return(character())
    }
    trajectory_ids <- which(grepl(
      "_(GF|SPF)_Trajectory_Combined[.]png$",
      basename(manifest$images$name)
    ))
    individual_ids <- which(grepl(
      "_Spatial_MAP[.]png$",
      basename(manifest$images$name)
    ))
    ids <- if (isTRUE(input$spatial_map_show_combined)) {
      c(trajectory_ids, individual_ids)
    } else individual_ids
    if (!length(ids)) ids <- seq_len(nrow(manifest$images))
    as.character(ids)
  }
  spatial_map_sample_choices <- function(manifest = spatial_map_manifest()) {
    if (is.null(manifest) || is.null(manifest$images) || !nrow(manifest$images)) {
      return(character())
    }
    ids <- which(grepl(
      "_Spatial_MAP[.]png$", basename(manifest$images$name)
    ))
    if (!length(ids)) return(character())
    samples <- unique(manifest$images$sample[ids])
    samples[!is.na(samples) & nzchar(samples)]
  }
  spatial_map_visible_image_ids <- function(
      manifest = spatial_map_manifest(),
      targets = input$spatial_map_target_choice,
      samples = input$spatial_map_image_choice) {
    ids <- spatial_map_display_image_ids(manifest)
    if (!length(ids) || !length(targets) || !length(samples)) {
      return(character())
    }
    rows <- manifest$images[as.integer(ids), , drop = FALSE]
    selected_groups <- unique(ifelse(
      grepl("^GF", samples, ignore.case = TRUE), "GF",
      ifelse(grepl("^SPF", samples, ignore.case = TRUE), "SPF", NA_character_)
    ))
    selected_groups <- selected_groups[!is.na(selected_groups)]
    is_individual <- rows$view == "individual"
    is_trajectory <- rows$view == "trajectory_combined"
    trajectory_group <- sub("-Combined$", "", rows$sample)
    keep <- !is.na(rows$target) & rows$target %in% targets & (
      (is_individual & !is.na(rows$sample) & rows$sample %in% samples) |
      (isTRUE(input$spatial_map_show_combined) & is_trajectory &
         trajectory_group %in% selected_groups)
    )
    ids[keep]
  }
  spatial_map_image_label <- function(name) {
    filename <- basename(name)
    sample_name <- sub("_.*$", "", filename)
    if (grepl("^(GF|SPF)", sample_name, ignore.case = TRUE)) {
      sample_name
    } else {
      filename
    }
  }
  selected_spatial_map_images <- function(ids) {
    manifest <- spatial_map_manifest(); paths <- spatial_map_run_paths()
    req(manifest, paths, length(ids) > 0L, !is.null(manifest$images))
    index <- unique(suppressWarnings(as.integer(ids)))
    req(length(index) > 0L, all(!is.na(index)), all(index >= 1L),
        all(index <= nrow(manifest$images)))
    rows <- manifest$images[index, , drop = FALSE]
    req(all(file.exists(rows$path)))
    normalized_paths <- vapply(
      rows$path, normalizePath, character(1),
      winslash = "/", mustWork = TRUE
    )
    root <- normalizePath(paths$root, winslash = "/", mustWork = TRUE)
    req(all(normalized_paths == root |
            startsWith(normalized_paths, paste0(root, "/"))))
    rows$path <- normalized_paths
    rows
  }
  selected_spatial_map_image <- function(id) {
    rows <- selected_spatial_map_images(id)
    req(nrow(rows) == 1L)
    rows
  }
  register_spatial_map_image_outputs <- function(ids) {
    for (id in ids) {
      local({
        image_id <- as.character(id)
        output[[paste0("spatial_map_image_", image_id)]] <- renderImage({
          row <- selected_spatial_map_image(image_id)
          extension <- tolower(tools::file_ext(row$path[[1]]))
          list(
            src = row$path[[1]],
            contentType = if (extension == "png") "image/png" else "image/jpeg",
            alt = row$name[[1]]
          )
        }, deleteFile = FALSE)
      })
    }
  }
  render_csv_table <- function(row) {
    path <- row$path[[1]]
    separator <- if (tolower(tools::file_ext(path)) == "tsv") "\t" else ","
    sql <- sprintf(paste0(
      "SELECT * FROM read_csv(%s, delim = %s, header = true, all_varchar = true, ",
      "sample_size = 20480) LIMIT %d"
    ), as.character(DBI::dbQuoteString(con, path)),
       as.character(DBI::dbQuoteString(con, separator)), APP_CONFIG$qc_table_preview_rows)
    value <- tryCatch(DBI::dbGetQuery(con, sql), error = function(e) e)
    if (inherits(value, "error")) validate(need(FALSE, conditionMessage(value)))
    DT::datatable(value, rownames = FALSE, escape = TRUE,
                  options = list(scrollX = TRUE, pageLength = 10))
  }

  output$roi_metadata_file_table <- DT::renderDT({
    manifest <- roi_metadata_manifest(); req(manifest)
    files <- manifest$files
    DT::datatable(
      data.frame(File = files$name, Type = files$kind,
                 KiB = round(files$size_bytes / 1024, 1)),
      rownames = FALSE, escape = TRUE,
      options = list(scrollX = TRUE, pageLength = 6, lengthChange = FALSE)
    )
  }, server = TRUE)
  output$roi_metadata_image <- renderImage({
    row <- selected_roi_metadata_file(input$roi_metadata_plot_choice, "plot")
    extension <- tolower(tools::file_ext(row$path[[1]]))
    list(src = row$path[[1]],
         contentType = if (extension == "png") "image/png" else "image/jpeg",
         alt = row$name[[1]])
  }, deleteFile = FALSE)
  output$roi_metadata_plot_empty <- renderUI({
    manifest <- roi_metadata_manifest()
    if (is.null(manifest) || !any(manifest$files$kind == "plot")) {
      div(class = "empty-state", icon("location-crosshairs"),
          strong("ROI coordinate QC"),
          span("Coordinate QC plots appear after ROI metadata is generated."))
    }
  })
  output$roi_metadata_table <- DT::renderDT({
    render_csv_table(selected_roi_metadata_file(input$roi_metadata_table_choice, "table"))
  }, server = TRUE)
  output$roi_metadata_table_empty <- renderUI({
    manifest <- roi_metadata_manifest()
    if (is.null(manifest) || !any(manifest$files$kind == "table")) {
      p(class = "muted", "CSV/TSV results appear after ROI metadata is generated.")
    }
  })
  output$download_roi_metadata_file <- downloadHandler(
    filename = function() basename(selected_roi_metadata_file(
      input$roi_metadata_file_choice)$path[[1]]),
    content = function(file) {
      path <- selected_roi_metadata_file(input$roi_metadata_file_choice)$path[[1]]
      if (!file.copy(path, file, overwrite = TRUE)) stop("Failed to download the ROI metadata result.")
    }
  )
  output$roi_metadata_log <- renderText({
    roi_metadata_log_tick(); roi_metadata_state()
    paths <- continuous_run_paths()
    if (is.null(paths)) return("ROI metadata has not been run yet.")
    log <- file.path(paths$metadata_dir, "ROI_Metadata_execution.log")
    if (!file.exists(log)) return("Preparing the log file.")
    paste(readLines(log, n = 2000L, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  })
  output$download_roi_metadata_log <- downloadHandler(
    filename = function() "ROI_Metadata_execution.log",
    content = function(file) {
      paths <- continuous_run_paths(); req(paths)
      log <- file.path(paths$metadata_dir, "ROI_Metadata_execution.log")
      req(file.exists(log))
      if (!file.copy(log, file, overwrite = TRUE)) stop("Failed to download the ROI metadata log.")
    }
  )

  output$zonation_summary <- DT::renderDT({
    manifest <- zonation_manifest(); req(manifest, !is.null(manifest$summary))
    DT::datatable(manifest$summary, rownames = FALSE, escape = TRUE,
                  options = list(scrollX = TRUE, pageLength = 10, lengthChange = FALSE))
  }, server = TRUE)
  output$zonation_summary_empty <- renderUI({
    manifest <- zonation_manifest()
    if (is.null(manifest) || is.null(manifest$summary)) {
      p(class = "muted", "An analysis summary appears after Continuous zonation is complete.")
    }
  })
  output$zonation_file_table <- DT::renderDT({
    manifest <- zonation_manifest(); req(manifest)
    files <- manifest$files
    DT::datatable(
      data.frame(File = files$name, Type = files$kind,
                 KiB = round(files$size_bytes / 1024, 1)),
      rownames = FALSE, escape = TRUE,
      options = list(scrollX = TRUE, pageLength = 8, lengthChange = FALSE)
    )
  }, server = TRUE)
  output$zonation_image <- renderImage({
    row <- selected_zonation_file(input$zonation_plot_choice, "plot")
    extension <- tolower(tools::file_ext(row$path[[1]]))
    list(src = row$path[[1]],
         contentType = if (extension == "png") "image/png" else "image/jpeg",
         alt = row$name[[1]])
  }, deleteFile = FALSE)
  output$zonation_plot_empty <- renderUI({
    manifest <- zonation_manifest()
    if (is.null(manifest) || !any(manifest$files$kind == "plot")) {
      div(class = "empty-state", icon("chart-line"),
          strong("Continuous zonation figures"),
          span("PNG/JPG results appear after analysis is complete."))
    }
  })
  output$zonation_table <- DT::renderDT({
    render_csv_table(selected_zonation_file(input$zonation_table_choice, "table"))
  }, server = TRUE)
  output$zonation_table_empty <- renderUI({
    manifest <- zonation_manifest()
    if (is.null(manifest) || !any(manifest$files$kind == "table")) {
      p(class = "muted", "CSV/TSV results appear after analysis is complete.")
    }
  })
  output$download_zonation_file <- downloadHandler(
    filename = function() basename(selected_zonation_file(
      input$zonation_file_choice)$path[[1]]),
    content = function(file) {
      path <- selected_zonation_file(input$zonation_file_choice)$path[[1]]
      if (!file.copy(path, file, overwrite = TRUE)) stop("Failed to download the Zonation result.")
    }
  )
  output$zonation_log <- renderText({
    zonation_log_tick(); zonation_state()
    paths <- continuous_run_paths()
    if (is.null(paths)) return("Continuous zonation has not been run yet.")
    log <- file.path(paths$zonation_dir, "Continuous_Zonation_execution.log")
    if (!file.exists(log)) return("Preparing the log file.")
    paste(readLines(log, n = 2000L, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  })
  output$download_zonation_log <- downloadHandler(
    filename = function() "Continuous_Zonation_execution.log",
    content = function(file) {
      paths <- continuous_run_paths(); req(paths)
      log <- file.path(paths$zonation_dir, "Continuous_Zonation_execution.log")
      req(file.exists(log))
      if (!file.copy(log, file, overwrite = TRUE)) stop("Failed to download the Zonation log.")
    }
  )

  output$trajectory_summary <- DT::renderDT({
    manifest <- trajectory_manifest(); req(manifest, !is.null(manifest$summary))
    DT::datatable(manifest$summary, rownames = FALSE, escape = TRUE,
                  options = list(scrollX = TRUE, pageLength = 5,
                                 lengthChange = FALSE, searching = FALSE))
  }, server = TRUE)
  output$trajectory_file_table <- DT::renderDT({
    manifest <- trajectory_manifest(); req(manifest)
    files <- manifest$files
    DT::datatable(
      data.frame(File = files$name, Type = files$kind,
                 KiB = round(files$size_bytes / 1024, 1)),
      rownames = FALSE, escape = TRUE,
      options = list(scrollX = TRUE, pageLength = 8, lengthChange = FALSE)
    )
  }, server = TRUE)
  selected_trajectory_plot <- function(type, id) {
    manifest <- trajectory_manifest(); req(manifest, id)
    plots <- if (identical(type, "beta")) manifest$beta_plots else manifest$curve_plots
    index <- suppressWarnings(as.integer(id))
    req(length(index) == 1L, !is.na(index), index >= 1L, index <= nrow(plots))
    row <- plots[index, , drop = FALSE]
    req(file.exists(row$PNG[[1]]))
    row
  }
  output$trajectory_beta_image <- renderImage({
    row <- selected_trajectory_plot("beta", input$trajectory_beta_choice)
    list(src = row$PNG[[1]], contentType = "image/png",
         alt = paste(row$Target[[1]], "zonation coefficient"))
  }, deleteFile = FALSE)
  output$trajectory_curve_image <- renderImage({
    row <- selected_trajectory_plot("curve", input$trajectory_curve_choice)
    list(src = row$PNG[[1]], contentType = "image/png",
         alt = paste(row$Target[[1]], "GF SPF trajectory"))
  }, deleteFile = FALSE)
  output$trajectory_beta_empty <- renderUI({
    manifest <- trajectory_manifest()
    if (is.null(manifest) || !nrow(manifest$beta_plots)) {
      div(class = "empty-state", icon("circle-nodes"),
          strong("Zonation coefficient"),
          span("No Beta plots were generated. Check the status of each target in the summary."))
    }
  })
  output$trajectory_curve_empty <- renderUI({
    manifest <- trajectory_manifest()
    if (is.null(manifest) || !nrow(manifest$curve_plots)) {
      div(class = "empty-state", icon("chart-line"),
          strong("GF/SPF trajectory"),
          span("No trajectory plots were generated. Check the status of each target in the summary."))
    }
  })
  trajectory_plot_download <- function(type, extension, input_id, fallback) {
    downloadHandler(
      filename = function() {
        row <- selected_trajectory_plot(type, input[[input_id]])
        basename(row[[extension]][[1]])
      },
      content = function(file) {
        row <- selected_trajectory_plot(type, input[[input_id]])
        path <- row[[extension]][[1]]; req(file.exists(path))
        if (!file.copy(path, file, overwrite = TRUE)) stop(fallback)
      }
    )
  }
  output$download_trajectory_beta_png <- trajectory_plot_download(
    "beta", "PNG", "trajectory_beta_choice",
    "Failed to download the Zonation coefficient PNG.")
  output$download_trajectory_curve_png <- trajectory_plot_download(
    "curve", "PNG", "trajectory_curve_choice",
    "Failed to download the Trajectory PNG.")
  output$download_trajectory_file <- downloadHandler(
    filename = function() basename(selected_trajectory_file(
      input$trajectory_file_choice)$path[[1]]),
    content = function(file) {
      path <- selected_trajectory_file(input$trajectory_file_choice)$path[[1]]
      if (!file.copy(path, file, overwrite = TRUE)) {
        stop("Failed to download the Trajectory result.")
      }
    }
  )
  output$trajectory_log <- renderText({
    trajectory_log_tick(); trajectory_state()
    paths <- trajectory_run_paths()
    if (is.null(paths)) return("Trajectory analysis has not been run yet.")
    log <- file.path(paths$root, "Trajectory_execution.log")
    if (!file.exists(log)) return("Preparing the log file.")
    paste(readLines(log, n = 2000L, warn = FALSE, encoding = "UTF-8"),
          collapse = "\n")
  })
  output$download_trajectory_log <- downloadHandler(
    filename = function() "Trajectory_execution.log",
    content = function(file) {
      paths <- trajectory_run_paths(); req(paths)
      log <- file.path(paths$root, "Trajectory_execution.log")
      req(file.exists(log))
      if (!file.copy(log, file, overwrite = TRUE)) {
        stop("Failed to download the Trajectory execution log.")
      }
    }
  )

  output$spatial_map_target_summary <- DT::renderDT({
    manifest <- spatial_map_manifest(); req(manifest, manifest$summary)
    summary <- manifest$summary
    DT::datatable(
      data.frame(
        Target = summary$Target,
        Status = summary$Status,
        Overlay_Images = summary$Overlay_Images,
        Elapsed_Seconds = summary$Elapsed_Seconds,
        Message = summary$Message,
        check.names = FALSE
      ),
      rownames = FALSE, escape = TRUE,
      options = list(scrollX = TRUE, paging = FALSE, searching = FALSE,
                     info = FALSE)
    )
  }, server = FALSE)

  output$spatial_map_file_table <- DT::renderDT({
    manifest <- spatial_map_manifest(); req(manifest)
    files <- manifest$files
    DT::datatable(
      data.frame(Target = ifelse(is.na(files$target), "—", files$target),
                 File = files$name, Type = files$kind,
                 KiB = round(files$size_bytes / 1024, 1)),
      rownames = FALSE, escape = TRUE,
      options = list(scrollX = TRUE, pageLength = 10, lengthChange = FALSE)
    )
  }, server = TRUE)
  output$spatial_map_image_grid <- renderUI({
    manifest <- spatial_map_manifest()
    if (is.null(manifest) || is.null(manifest$images) || !nrow(manifest$images)) {
      return(div(class = "empty-state", icon("map"),
                 strong("Spatial Map"),
                 span("After analysis runs, histology-overlay images for eight samples appear here.")))
    }
    selected_targets <- input$spatial_map_target_choice
    selected_samples <- input$spatial_map_image_choice
    if (!length(selected_targets)) {
      return(div(class = "empty-state", icon("bullseye"),
                 strong("Select at least one target to display."),
                 span("Select one or more targets using the checkboxes above.")))
    }
    if (!length(selected_samples)) {
      return(div(class = "empty-state", icon("images"),
                 strong("Select samples to display."),
                 span("Use the checkboxes above or the GF only/SPF only buttons.")))
    }
    ids <- spatial_map_visible_image_ids(
      manifest, selected_targets, selected_samples
    )
    if (!length(ids)) {
      return(div(class = "empty-state", icon("triangle-exclamation"),
                 strong("No images match the selected combination."),
                 span("Check the target and sample selections.")))
    }
    rows <- selected_spatial_map_images(ids)
    rows$image_id <- ids

    if (isTRUE(input$spatial_map_show_combined)) {
      selected_groups <- unique(ifelse(
        grepl("^GF", selected_samples, ignore.case = TRUE), "GF",
        ifelse(grepl("^SPF", selected_samples, ignore.case = TRUE),
               "SPF", NA_character_)
      ))
      selected_groups <- selected_groups[!is.na(selected_groups)]
      sections <- list()
      section_index <- 0L

      for (target_name in selected_targets) {
        for (group_name in c("GF", "SPF")) {
          if (!group_name %in% selected_groups) next
          group_samples <- selected_samples[
            grepl(paste0("^", group_name), selected_samples,
                  ignore.case = TRUE)
          ]
          map_rows <- rows[
            rows$target == target_name & rows$view == "individual" &
              rows$sample %in% group_samples,
            , drop = FALSE
          ]
          trajectory_rows <- rows[
            rows$target == target_name &
              rows$view == "trajectory_combined" &
              rows$sample == paste0(group_name, "-Combined"),
            , drop = FALSE
          ]
          if (!nrow(map_rows) && !nrow(trajectory_rows)) next
          if (nrow(map_rows)) {
            map_rows <- map_rows[order(match(map_rows$sample, group_samples)),
                                 , drop = FALSE]
          }

          trajectory_ui <- if (nrow(trajectory_rows)) {
            imageOutput(
              paste0("spatial_map_image_", trajectory_rows$image_id[[1]]),
              width = "100%", height = "auto"
            )
          } else {
            div(class = "empty-state", icon("chart-line"),
                strong("Cannot find the Trajectory panel."),
                span("Regenerate it using the latest Spatial Expression script."))
          }

          map_cards <- lapply(seq_len(nrow(map_rows)), function(i) {
            div(
              class = "spatial-map-card",
              div(class = "spatial-map-card-title", map_rows$sample[[i]]),
              imageOutput(
                paste0("spatial_map_image_", map_rows$image_id[[i]]),
                width = "100%", height = "auto"
              )
            )
          })

          section_index <- section_index + 1L
          sections[[section_index]] <- div(
            class = "spatial-combined-section",
            div(
              class = "spatial-combined-title",
              paste0(target_name, " · ", group_name, "-Combined")
            ),
            div(
              class = "spatial-combined-row",
              div(class = "spatial-combined-trajectory", trajectory_ui),
              div(class = "spatial-combined-histology", tagList(map_cards))
            )
          )
        }
      }
      if (!length(sections)) {
        return(div(class = "empty-state", icon("triangle-exclamation"),
                   strong("No combined panels are available to display."),
                   span("Check the target and GF/SPF sample selections.")))
      }
      return(div(class = "spatial-combined-stack", tagList(sections)))
    }

    cards <- lapply(seq_len(nrow(rows)), function(i) {
      id <- rows$image_id[[i]]
      div(
        class = "spatial-map-card",
        div(class = "spatial-map-card-title",
            paste0(rows$target[[i]], " · ", rows$sample[[i]])),
        imageOutput(paste0("spatial_map_image_", id),
                    width = "100%", height = "auto")
      )
    })
    div(
      class = "spatial-map-grid",
      tagList(cards)
    )
  })
  output$download_spatial_map_image <- downloadHandler(
    filename = function() {
      rows <- selected_spatial_map_images(spatial_map_visible_image_ids())
      if (nrow(rows) == 1L) return(basename(rows$path[[1]]))
      targets <- unique(rows$target)
      target_label <- if (length(targets) == 1L) targets[[1]] else "MultiTarget"
      paste0(gsub("[^A-Za-z0-9_-]", "_", target_label),
             "_Selected_Spatial_Maps.zip")
    },
    content = function(file) {
      rows <- selected_spatial_map_images(spatial_map_visible_image_ids())
      if (nrow(rows) == 1L) {
        if (!file.copy(rows$path[[1]], file, overwrite = TRUE)) {
          stop("Failed to download the Spatial Map PNG.")
        }
        return(invisible(NULL))
      }
      zip_dir <- tempfile("spatial-map-download-")
      if (!dir.create(zip_dir, recursive = TRUE, showWarnings = FALSE)) {
        stop("Cannot create a temporary directory for the selected images.")
      }
      on.exit(unlink(zip_dir, recursive = TRUE, force = TRUE), add = TRUE)
      zip_names <- make.unique(basename(rows$path))
      copied <- file.copy(rows$path, file.path(zip_dir, zip_names),
                          overwrite = TRUE)
      if (!all(copied)) stop("Failed to copy the selected images.")
      old_wd <- getwd()
      on.exit(setwd(old_wd), add = TRUE)
      setwd(zip_dir)
      zip_path <- normalizePath(file, winslash = "/", mustWork = FALSE)
      suppressWarnings(utils::zip(zipfile = zip_path, files = zip_names))
      if (!file.exists(file) || is.na(file.size(file)) || file.size(file) <= 0) {
        stop("Failed to create the ZIP for selected images.")
      }
    }
  )
  output$download_spatial_map_file <- downloadHandler(
    filename = function() basename(selected_spatial_map_file(
      input$spatial_map_file_choice)$path[[1]]),
    content = function(file) {
      path <- selected_spatial_map_file(input$spatial_map_file_choice)$path[[1]]
      if (!file.copy(path, file, overwrite = TRUE)) {
        stop("Failed to download the Spatial Map result.")
      }
    }
  )
  output$spatial_map_log <- renderText({
    spatial_map_log_tick(); spatial_map_state()
    paths <- spatial_map_run_paths()
    if (is.null(paths)) return("Spatial Map has not been run yet.")
    log <- file.path(paths$root, "Spatial_Map_execution.log")
    if (!file.exists(log)) return("Preparing the log file.")
    paste(readLines(log, n = 2000L, warn = FALSE, encoding = "UTF-8"),
          collapse = "\n")
  })
  output$download_spatial_map_log <- downloadHandler(
    filename = function() "Spatial_Map_execution.log",
    content = function(file) {
      paths <- spatial_map_run_paths(); req(paths)
      log <- file.path(paths$root, "Spatial_Map_execution.log")
      req(file.exists(log))
      if (!file.copy(log, file, overwrite = TRUE)) {
        stop("Failed to download the Spatial Map execution log.")
      }
    }
  )

  draw_run_plot <- function(metric, color) {
    runs <- run_result()
    validate(need(!is.null(runs), "Plots appear after the full summary is calculated."))
    runs <- head(runs, APP_CONFIG$max_run_summary_rows)
    n_show <- min(20L, nrow(runs), APP_CONFIG$max_run_summary_rows)
    runs <- runs[order(runs[[metric]], decreasing = TRUE), , drop = FALSE]
    runs <- runs[seq_len(n_show), , drop = FALSE]
    runs <- runs[rev(seq_len(nrow(runs))), , drop = FALSE]
    labels <- ifelse(nchar(runs$run_name) > 29L,
                     paste0(substr(runs$run_name, 1, 26), "..."), runs$run_name)
    old <- par(mar = c(4, 11, 2, 1), bg = "white")
    on.exit(par(old))
    barplot(runs[[metric]], names.arg = labels, horiz = TRUE, las = 1,
            col = color, border = NA, cex.names = 0.75,
            xlab = if (metric == "n_groups") "Protein groups" else "Report rows")
    title(main = paste("Top", n_show, "Run (based on Runs included in the table)"), cex.main = 0.85)
  }
  output$protein_plot <- renderPlot({ draw_run_plot("n_groups", APP_CONFIG$primary_color) }, res = 110)
  output$rows_plot <- renderPlot({
    draw_run_plot("n_rows", APP_CONFIG$light_color)
  }, res = 110)

  output$run_limit_note <- renderUI({
    runs <- run_result()
    if (!is.null(runs) && nrow(runs) > APP_CONFIG$max_run_summary_rows) {
      p(class = "muted", paste("The table and plots are limited to the first",
                              APP_CONFIG$max_run_summary_rows,
                              "Runs in Run-name order. The metrics above are based on the entire file."))
    }
  })
  output$run_table <- DT::renderDT({
    runs <- run_result(); req(runs)
    DT::datatable(head(runs, APP_CONFIG$max_run_summary_rows), rownames = FALSE,
      colnames = c("Run", "Report rows", "Protein groups"), escape = TRUE,
      options = list(scrollX = TRUE, pageLength = 10))
  }, server = TRUE)

  output$explorer_note <- renderText(explorer_message())
  output$explorer_table <- DT::renderDT({
    df <- explorer_result(); req(!is.null(df))
    DT::datatable(df, rownames = FALSE, escape = TRUE,
      options = list(scrollX = TRUE, pageLength = 20))
  }, server = TRUE)
  output$download_preview <- downloadHandler(
    filename = function() "diann_query_preview.csv",
    content = function(file) {
      df <- explorer_result(); req(!is.null(df))
      utils::write.csv(df, file, row.names = FALSE, na = "", fileEncoding = "UTF-8")
    }
  )
}

shinyApp(ui = ui, server = server)
