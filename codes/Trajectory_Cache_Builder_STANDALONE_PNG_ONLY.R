# Standalone fast trajectory-cache builder.
# Contains the cache-compatible Beta and trajectory implementations.

# Targeted GF vs SPF zonation-coefficient plot
#
# Sourcing this file only defines functions. The Shiny app calls
# run_zonation_coefficient_beta() after Continuous analysis has finished.

.require_zonation_packages <- function() {
  packages <- c("ggplot2", "ggrepel")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Required package(s) missing: ", paste(missing, collapse = ", "),
         ". Install them before running this analysis.", call. = FALSE)
  }
}

.clean_target <- function(target) {
  if (is.null(target) || length(target) != 1L || is.na(target)) {
    stop("Enter one UniProt accession or gene name.", call. = FALSE)
  }
  target <- trimws(as.character(target))
  if (!nzchar(target)) stop("Enter one UniProt accession or gene name.", call. = FALSE)
  target
}

.target_tokens <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  strsplit(x, "[;,|]", perl = TRUE)
}

.matches_target <- function(x, target) {
  wanted <- toupper(trimws(target))
  full_match <- toupper(trimws(ifelse(is.na(x), "", as.character(x)))) == wanted
  full_match | vapply(.target_tokens(x), function(tokens) {
    any(toupper(trimws(tokens)) == wanted)
  }, logical(1))
}

.safe_target_name <- function(x) {
  value <- gsub("[^A-Za-z0-9._-]", "_", x)
  value <- substr(value, 1L, 80L)
  if (nzchar(value)) value else "Protein"
}

prepare_zonation_coefficient_data <- function(input_dir) {
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  input_file <- file.path(input_dir, "GF_vs_SPF_Zonation_Interaction.csv")
  if (!file.exists(input_file)) {
    stop("Continuous analysis output is missing: ", input_file, call. = FALSE)
  }
  dat <- utils::read.csv(input_file, check.names = FALSE,
                         stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  required <- c("Protein.Group", "Protein", "Beta_GF", "Beta_SPF",
                "Delta_Beta", "P_interaction", "Q_interaction")
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop(basename(input_file), " is missing column(s): ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  dat$Protein.Group <- trimws(as.character(dat$Protein.Group))
  dat$Protein <- trimws(as.character(dat$Protein))
  dat$Beta_GF <- suppressWarnings(as.numeric(dat$Beta_GF))
  dat$Beta_SPF <- suppressWarnings(as.numeric(dat$Beta_SPF))
  dat$Q_interaction <- suppressWarnings(as.numeric(dat$Q_interaction))
  dat <- dat[is.finite(dat$Beta_GF) & is.finite(dat$Beta_SPF), , drop = FALSE]
  if (!nrow(dat)) stop("No finite Beta_GF/Beta_SPF values were found.", call. = FALSE)
  dat
}

run_zonation_coefficient_beta <- function(input_dir, output_dir, target,
                                           q_cutoff = 0.05,
                                           png_dpi = 400L,
                                           prepared_data = NULL) {
  .require_zonation_packages()
  target <- .clean_target(target)
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) {
    stop("Could not create output directory: ", output_dir, call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

  dat <- if (is.null(prepared_data)) {
    prepare_zonation_coefficient_data(input_dir)
  } else {
    prepared_data
  }

  hit <- .matches_target(dat$Protein, target) |
    .matches_target(dat$Protein.Group, target)
  highlight <- dat[hit, , drop = FALSE]
  if (!nrow(highlight)) {
    stop("No exact match for '", target,
         "' in Protein or Protein.Group. Check the gene symbol or UniProt accession.",
         call. = FALSE)
  }
  matched_ids <- unique(trimws(as.character(highlight$Protein.Group)))
  matched_ids <- matched_ids[nzchar(matched_ids)]
  if (length(matched_ids) > 1L) {
    stop("Gene name '", target, "' matches multiple Protein.Group IDs: ",
         paste(matched_ids, collapse = ", "),
         ". Enter one UniProt accession to select a single protein group.",
         call. = FALSE)
  }

  dat$Interaction <- ifelse(!is.na(dat$Q_interaction) &
                              dat$Q_interaction < q_cutoff,
                            "Significant", "NS")
  dat$Highlight <- hit
  highlight <- dat[dat$Highlight, , drop = FALSE]
  highlight$Direction <- ifelse(
    highlight$Beta_GF < 0, "CV-enriched in GF",
    ifelse(highlight$Beta_GF > 0, "PV-enriched in GF", "No GF direction")
  )
  label <- trimws(as.character(highlight$Protein))
  label[is.na(label) | !nzchar(label)] <- highlight$Protein.Group[
    is.na(label) | !nzchar(label)
  ]
  highlight$Plot_Label <- label

  background_ns <- dat[!dat$Highlight & dat$Interaction == "NS", , drop = FALSE]
  background_sig <- dat[!dat$Highlight & dat$Interaction == "Significant", , drop = FALSE]
  axis_limit <- max(abs(c(dat$Beta_GF, dat$Beta_SPF)), na.rm = TRUE)
  if (!is.finite(axis_limit) || axis_limit == 0) axis_limit <- 1
  axis_limit <- ceiling(axis_limit * 10) / 10 * 1.08
  colors <- c("CV-enriched in GF" = "#6A3D9A",
              "PV-enriched in GF" = "#E67E22",
              "No GF direction" = "#168AAD")

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        linewidth = 0.45, color = "grey45") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        linewidth = 0.45, color = "grey45") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                         linewidth = 0.65, color = "grey25") +
    ggplot2::geom_point(data = background_ns,
                        ggplot2::aes(x = Beta_GF, y = Beta_SPF),
                        shape = 16, size = 1.7, alpha = 0.35,
                        color = "#BDBDBD") +
    ggplot2::geom_point(data = background_sig,
                        ggplot2::aes(x = Beta_GF, y = Beta_SPF),
                        shape = 17, size = 2.0, alpha = 0.65,
                        color = "#4D4D4D") +
    ggplot2::geom_point(data = highlight,
                        ggplot2::aes(x = Beta_GF, y = Beta_SPF,
                                     color = Direction, shape = Interaction),
                        size = 4.2, alpha = 1, stroke = 0.8) +
    ggrepel::geom_text_repel(
      data = highlight,
      ggplot2::aes(x = Beta_GF, y = Beta_SPF,
                   label = Plot_Label, color = Direction),
      size = 3.8, fontface = "italic", box.padding = 0.5,
      point.padding = 0.4, min.segment.length = 0,
      segment.color = "grey55", segment.linewidth = 0.35,
      max.overlaps = Inf, seed = 123, show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(name = "Selected protein", values = colors) +
    ggplot2::scale_shape_manual(name = "GF x SPF interaction",
                                values = c("NS" = 16, "Significant" = 17)) +
    ggplot2::coord_equal(xlim = c(-axis_limit, axis_limit),
                         ylim = c(-axis_limit, axis_limit), clip = "off") +
    ggplot2::labs(
      title = "GF vs SPF liver zonation coefficients",
      subtitle = paste0("Selected target: ", target),
      x = "GF zonation coefficient (Beta)",
      y = "SPF zonation coefficient (Beta)"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(size = 11),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(12, 70, 12, 12)
    )

  stem <- paste0(.safe_target_name(target), "_Beta")
  png_file <- file.path(output_dir, paste0(stem, ".png"))
  table_file <- file.path(output_dir, paste0(stem, "_Matched.csv"))
  ggplot2::ggsave(png_file, p, width = 10, height = 8, dpi = png_dpi,
                  bg = "white")
  utils::write.csv(highlight, table_file, row.names = FALSE, na = "",
                   fileEncoding = "UTF-8")

  created <- c(png_file, table_file)
  if (any(!file.exists(created))) {
    missing_created <- created[!file.exists(created)]
    stop("Output file creation failed: ", paste(missing_created, collapse = ", "),
         ". On Windows, use a shorter QC results root if the full path is still too long.",
         call. = FALSE)
  }

  cat("Zonation coefficient plot completed for:", target, "\n")
  cat("Matched Protein.Group:", paste(matched_ids, collapse = ", "), "\n")
  list(
    plot = p,
    png_file = normalizePath(png_file, winslash = "/", mustWork = TRUE),
    pdf_file = NA_character_,
    table_file = normalizePath(table_file, winslash = "/", mustWork = TRUE),
    target = target,
    protein_group = paste(matched_ids, collapse = ", ")
  )
}

# Targeted GF/SPF CV-to-PV trajectory plot
#
# Sourcing this file only defines functions. The Shiny app calls
# run_selected_trajectory() with a UniProt accession or gene name.

.require_trajectory_packages <- function() {
  packages <- c("ggplot2")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Required package(s) missing: ", paste(missing, collapse = ", "),
         ". Install them before running this analysis.", call. = FALSE)
  }
}

.trajectory_clean_target <- function(target) {
  if (is.null(target) || length(target) != 1L || is.na(target)) {
    stop("Enter one UniProt accession or gene name.", call. = FALSE)
  }
  target <- trimws(as.character(target))
  if (!nzchar(target)) stop("Enter one UniProt accession or gene name.", call. = FALSE)
  target
}

.trajectory_tokens <- function(x) {
  x <- ifelse(is.na(x), "", as.character(x))
  strsplit(x, "[;,|]", perl = TRUE)
}

.trajectory_matches <- function(x, target) {
  wanted <- toupper(trimws(target))
  full_match <- toupper(trimws(ifelse(is.na(x), "", as.character(x)))) == wanted
  full_match | vapply(.trajectory_tokens(x), function(tokens) {
    any(toupper(trimws(tokens)) == wanted)
  }, logical(1))
}

.trajectory_safe_name <- function(x) {
  value <- gsub("[^A-Za-z0-9._-]", "_", x)
  value <- substr(value, 1L, 80L)
  if (nzchar(value)) value else "Protein"
}

.fit_trajectory_loess <- function(df, span, n_predict,
                                  min_points, min_unique_s) {
  df <- df[is.finite(df$Spatial_S) & is.finite(df$RobustScaled) &
             df$Spatial_S >= 0 & df$Spatial_S <= 1, , drop = FALSE]
  n_points <- nrow(df)
  n_unique <- length(unique(df$Spatial_S))
  empty <- data.frame(Spatial_S = numeric(), Predicted_Z = numeric())
  if (n_points < min_points || n_unique < min_unique_s) {
    return(list(curve = empty, status = "Insufficient data", notes = "",
                n_points = n_points, n_unique = n_unique))
  }
  notes <- character()
  result <- tryCatch({
    fit <- withCallingHandlers(
      stats::loess(RobustScaled ~ Spatial_S, data = df, span = span,
                   degree = 2, na.action = stats::na.exclude,
                   control = stats::loess.control(surface = "direct")),
      warning = function(w) {
        notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning")
      }
    )
    grid <- data.frame(Spatial_S = seq(min(df$Spatial_S), max(df$Spatial_S),
                                        length.out = n_predict))
    prediction <- withCallingHandlers(
      stats::predict(fit, newdata = grid),
      warning = function(w) {
        notes <<- c(notes, conditionMessage(w)); invokeRestart("muffleWarning")
      }
    )
    curve <- data.frame(Spatial_S = grid$Spatial_S,
                        Predicted_Z = as.numeric(prediction))
    curve <- curve[is.finite(curve$Predicted_Z), , drop = FALSE]
    status <- if (nrow(curve) < 2L) "No usable curve" else if (length(notes)) {
      "Fitted with warning"
    } else "Fitted"
    if (nrow(curve) < 2L) curve <- empty
    list(curve = curve, status = status,
         notes = paste(unique(notes), collapse = " | "),
         n_points = n_points, n_unique = n_unique)
  }, error = function(e) {
    list(curve = empty, status = "Fit failed", notes = conditionMessage(e),
         n_points = n_points, n_unique = n_unique)
  })
  result
}

prepare_selected_trajectory_data <- function(input_dir) {
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  analysis_file <- file.path(input_dir, "Spatial_Proteomics_AnalysisReady.csv")
  coefficient_file <- file.path(input_dir, "AllSamples_Zonation_Coefficients.csv")
  missing_files <- c(analysis_file, coefficient_file)[
    !file.exists(c(analysis_file, coefficient_file))
  ]
  if (length(missing_files)) {
    stop("Continuous analysis output file(s) missing: ",
         paste(missing_files, collapse = ", "), call. = FALSE)
  }

  dat <- utils::read.csv(analysis_file, check.names = FALSE,
                         stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  coefficients <- utils::read.csv(coefficient_file, check.names = FALSE,
                                  stringsAsFactors = FALSE, fileEncoding = "UTF-8")
  required_data <- c("Group", "Protein.Group", "Protein", "Spatial_S", "RobustScaled")
  required_coefficients <- c("Protein.Group", "Beta", "Q_value")
  missing_data <- setdiff(required_data, names(dat))
  missing_coefficients <- setdiff(required_coefficients, names(coefficients))
  if (length(missing_data)) {
    stop(basename(analysis_file), " is missing column(s): ",
         paste(missing_data, collapse = ", "), call. = FALSE)
  }
  if (length(missing_coefficients)) {
    stop(basename(coefficient_file), " is missing column(s): ",
         paste(missing_coefficients, collapse = ", "), call. = FALSE)
  }

  dat$Group <- toupper(trimws(as.character(dat$Group)))
  dat$Protein.Group <- trimws(as.character(dat$Protein.Group))
  dat$Protein <- trimws(as.character(dat$Protein))
  dat$Spatial_S <- suppressWarnings(as.numeric(dat$Spatial_S))
  dat$RobustScaled <- suppressWarnings(as.numeric(dat$RobustScaled))
  dat <- dat[dat$Group %in% c("GF", "SPF") &
               !is.na(dat$Protein.Group) & nzchar(dat$Protein.Group), , drop = FALSE]
  if (!all(c("GF", "SPF") %in% unique(dat$Group))) {
    stop("Both GF and SPF must be present in the Group column.", call. = FALSE)
  }
  coefficients$Protein.Group <- trimws(as.character(coefficients$Protein.Group))
  coefficients$Beta <- suppressWarnings(as.numeric(coefficients$Beta))
  coefficients$Q_value <- suppressWarnings(as.numeric(coefficients$Q_value))
  list(data = dat, coefficients = coefficients)
}

run_selected_trajectory <- function(input_dir, output_dir, target,
                                    loess_span = 0.6,
                                    n_predict = 200L,
                                    min_points = 10L,
                                    min_unique_s = 5L,
                                    png_dpi = 400L,
                                    prepared = NULL) {
  .require_trajectory_packages()
  target <- .trajectory_clean_target(target)
  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) {
    stop("Could not create output directory: ", output_dir, call. = FALSE)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

  prepared <- if (is.null(prepared)) prepare_selected_trajectory_data(input_dir) else prepared
  dat <- prepared$data
  coefficients <- prepared$coefficients
  required_data <- c("Group", "Protein.Group", "Protein", "Spatial_S", "RobustScaled")
  required_coefficients <- c("Protein.Group", "Beta", "Q_value")

  hit <- .trajectory_matches(dat$Protein, target) |
    .trajectory_matches(dat$Protein.Group, target)
  matched <- dat[hit, , drop = FALSE]
  if (!nrow(matched)) {
    stop("No exact match for '", target,
         "' in Protein or Protein.Group. Check the gene symbol or UniProt accession.",
         call. = FALSE)
  }
  selected_ids <- unique(matched$Protein.Group)
  if (length(selected_ids) > 1L) {
    stop("Gene name '", target, "' matches multiple Protein.Group IDs: ",
         paste(selected_ids, collapse = ", "),
         ". Enter one UniProt accession to select a single protein group.",
         call. = FALSE)
  }
  protein_id <- selected_ids[[1]]
  plot_data <- dat[dat$Protein.Group == protein_id &
                     is.finite(dat$Spatial_S) & is.finite(dat$RobustScaled) &
                     dat$Spatial_S >= 0 & dat$Spatial_S <= 1, required_data,
                   drop = FALSE]
  if (!nrow(plot_data)) stop("No usable observations for ", target, ".", call. = FALSE)

  protein_labels <- unique(plot_data$Protein[!is.na(plot_data$Protein) &
                                               nzchar(plot_data$Protein)])
  plot_label <- if (length(protein_labels)) protein_labels[[1]] else protein_id
  group_fits <- lapply(c("GF", "SPF"), function(group) {
    .fit_trajectory_loess(plot_data[plot_data$Group == group, , drop = FALSE],
                          loess_span, n_predict, min_points, min_unique_s)
  })
  names(group_fits) <- c("GF", "SPF")
  curves <- do.call(rbind, lapply(names(group_fits), function(group) {
    curve <- group_fits[[group]]$curve
    if (!nrow(curve)) return(NULL)
    curve$Group <- group
    curve$Protein.Group <- protein_id
    curve
  }))
  if (is.null(curves)) {
    curves <- data.frame(Spatial_S = numeric(), Predicted_Z = numeric(),
                         Group = character(), Protein.Group = character())
  }
  audit <- data.frame(
    Group = c("GF", "SPF"),
    Protein.Group = protein_id,
    Protein = plot_label,
    N_Points = vapply(group_fits, `[[`, integer(1), "n_points"),
    N_Unique_S = vapply(group_fits, `[[`, integer(1), "n_unique"),
    Fit_Status = vapply(group_fits, `[[`, character(1), "status"),
    Fit_Notes = vapply(group_fits, `[[`, character(1), "notes"),
    stringsAsFactors = FALSE
  )

  coefficient_row <- coefficients[
    trimws(as.character(coefficients$Protein.Group)) == protein_id,
    required_coefficients, drop = FALSE
  ]
  if (nrow(coefficient_row) > 1L) coefficient_row <- coefficient_row[1L, , drop = FALSE]
  beta_note <- if (nrow(coefficient_row)) {
    paste0("All samples: Beta = ", format(as.numeric(coefficient_row$Beta[[1]]), digits = 3),
           ", Q = ", format(as.numeric(coefficient_row$Q_value[[1]]), digits = 3))
  } else "All-sample Beta/Q metadata unavailable"

  absent <- setdiff(c("GF", "SPF"), unique(plot_data$Group))
  no_curve <- setdiff(unique(plot_data$Group), unique(curves$Group))
  caption <- character()
  if (length(absent)) caption <- c(caption, paste("No observations:", paste(absent, collapse = ", ")))
  if (length(no_curve)) caption <- c(caption, paste("No LOESS curve:", paste(no_curve, collapse = ", ")))

  values <- c(plot_data$RobustScaled, curves$Predicted_Z)
  y_range <- range(values[is.finite(values)], na.rm = TRUE)
  padding <- max(diff(y_range) * 0.08, 0.1)
  y_limits <- y_range + c(-padding, padding)
  colors <- c(GF = "#71468E", SPF = "#D96A1C")
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(x = Spatial_S, y = RobustScaled, color = Group, shape = Group),
      size = 1.3, alpha = 0.10, stroke = 0
    ) +
    ggplot2::geom_line(
      data = curves,
      ggplot2::aes(x = Spatial_S, y = Predicted_Z,
                   color = Group, linetype = Group, group = Group),
      linewidth = 1.1, alpha = 0.95
    ) +
    ggplot2::scale_color_manual(name = "Group", values = colors,
                                limits = c("GF", "SPF"), drop = FALSE) +
    ggplot2::scale_linetype_manual(name = "Group",
                                   values = c(GF = "solid", SPF = "dashed"),
                                   limits = c("GF", "SPF"), drop = FALSE) +
    ggplot2::scale_shape_manual(name = "Group", values = c(GF = 16, SPF = 17),
                                limits = c("GF", "SPF"), drop = FALSE) +
    ggplot2::scale_x_continuous(
      limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("CV (0)", "0.25", "0.50", "0.75", "PV (1)"),
      expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = "dotted",
                        linewidth = 0.3, alpha = 0.4) +
    ggplot2::coord_cartesian(ylim = y_limits) +
    ggplot2::labs(
      title = plot_label,
      subtitle = paste("GF and SPF along the CV-PV axis |", beta_note),
      x = "Spatial ratio S", y = "Robust-scaled protein abundance",
      caption = if (length(caption)) paste(caption, collapse = " | ") else NULL
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(alpha = 1))) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "right",
      plot.margin = ggplot2::margin(10, 15, 10, 10),
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      axis.text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 13)
    )

  stem <- paste0(.trajectory_safe_name(target), "_Trajectory")
  png_file <- file.path(output_dir, paste0(stem, ".png"))
  observations_file <- file.path(output_dir, paste0(stem, "_Observations.csv"))
  curves_file <- file.path(output_dir, paste0(stem, "_LOESS_Curves.csv"))
  audit_file <- file.path(output_dir, paste0(stem, "_Fit_Audit.csv"))
  ggplot2::ggsave(png_file, p, width = 8, height = 6, dpi = png_dpi, bg = "white")
  utils::write.csv(plot_data, observations_file, row.names = FALSE, na = "",
                   fileEncoding = "UTF-8")
  utils::write.csv(curves, curves_file, row.names = FALSE, na = "",
                   fileEncoding = "UTF-8")
  utils::write.csv(audit, audit_file, row.names = FALSE, na = "",
                   fileEncoding = "UTF-8")

  created <- c(png_file, observations_file, curves_file, audit_file)
  if (any(!file.exists(created))) {
    missing_created <- created[!file.exists(created)]
    stop("Output file creation failed: ", paste(missing_created, collapse = ", "),
         ". On Windows, use a shorter QC results root if the full path is still too long.",
         call. = FALSE)
  }

  cat("GF/SPF trajectory completed for:", target, "\n")
  cat("Matched Protein.Group:", protein_id, "\n")
  list(
    plot = p,
    png_file = normalizePath(png_file, winslash = "/", mustWork = TRUE),
    pdf_file = NA_character_,
    observations_file = normalizePath(observations_file, winslash = "/", mustWork = TRUE),
    curves_file = normalizePath(curves_file, winslash = "/", mustWork = TRUE),
    audit_file = normalizePath(audit_file, winslash = "/", mustWork = TRUE),
    target = target,
    protein_group = protein_id,
    protein = plot_label
  )
}

# Build reusable single-protein Beta and GF/SPF trajectory figures.
#
# Place this single standalone file in the codes folder.
# The compatible Beta and trajectory functions are embedded below.
# Keep app(C)_F.R one level above codes.
#
# Run this script independently from Shiny.  It creates:
#   Trajectory_Cache/Cache_Index.csv
#   Trajectory_Cache/Beta/*.png
#   Trajectory_Cache/Curve/*.png
#
# The Shiny app checks this cache before running a new analysis.

CACHE_BUILDER_VERSION <- "2026-11-13-v6-standalone-png-only"

############################################################
# USER SETTINGS
############################################################

# Folder produced by Continuous analysis. It must contain:
#   GF_vs_SPF_Zonation_Interaction.csv
#   Spatial_Proteomics_AnalysisReady.csv
#   AllSamples_Zonation_Coefficients.csv
# Leave empty on Windows to select the folder interactively.
CONTINUOUS_OUTPUT_DIR <- ""

# Empty means all Protein.Group entries.  To prebuild only a panel, enter IDs:
# TARGET_PROTEIN_GROUPS <- c("Q05421", "P12345")
TARGET_PROTEIN_GROUPS <- character(0)

# Existing complete cache entries are skipped by default, so the script can be
# stopped and resumed safely.
OVERWRITE_EXISTING <- FALSE

# Auxiliary per-protein CSV files are not needed by Shiny and consume space.
KEEP_AUXILIARY_CSV <- FALSE

# Lower this to 200 for smaller preview files; 300 is publication-friendly.
PNG_DPI <- 300L

############################################################
# LOCATE APP FILES
############################################################

find_this_script <- function() {
  frames <- sys.frames()
  candidates <- unlist(lapply(frames, function(frame) {
    if (exists("ofile", envir = frame, inherits = FALSE)) frame$ofile else NULL
  }), use.names = FALSE)
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", args[grepl("^--file=", args)])
  candidates <- c(candidates, file_arg)
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  if (!length(candidates)) return(NA_character_)
  normalizePath(tail(candidates, 1L), winslash = "/", mustWork = TRUE)
}

script_file <- find_this_script()
CODES_DIR <- if (is.na(script_file)) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  dirname(script_file)
}
cat("Trajectory cache builder version:", CACHE_BUILDER_VERSION, "\n")
cat("Running builder file:",
    if (is.na(script_file)) "<RStudio working-directory mode>" else script_file,
    "\n")
cat("Analysis code folder:", CODES_DIR, "\n\n")
APP_DIR <- if (tolower(basename(CODES_DIR)) == "codes") dirname(CODES_DIR) else CODES_DIR
CACHE_DIR <- file.path(APP_DIR, "Trajectory_Cache")
BETA_DIR <- file.path(CACHE_DIR, "Beta")
CURVE_DIR <- file.path(CACHE_DIR, "Curve")

if (!nzchar(trimws(CONTINUOUS_OUTPUT_DIR))) {
  if (.Platform$OS.type == "windows") {
    CONTINUOUS_OUTPUT_DIR <- utils::choose.dir(
      default = APP_DIR,
      caption = "Continuous analysis 결과 폴더를 선택하세요"
    )
  } else {
    stop("Set CONTINUOUS_OUTPUT_DIR at the top of this script.", call. = FALSE)
  }
}
if (is.na(CONTINUOUS_OUTPUT_DIR) || !nzchar(CONTINUOUS_OUTPUT_DIR)) {
  stop("Continuous analysis folder selection was cancelled.", call. = FALSE)
}
CONTINUOUS_OUTPUT_DIR <- normalizePath(CONTINUOUS_OUTPUT_DIR,
                                       winslash = "/", mustWork = TRUE)

required_inputs <- file.path(CONTINUOUS_OUTPUT_DIR, c(
  "GF_vs_SPF_Zonation_Interaction.csv",
  "Spatial_Proteomics_AnalysisReady.csv",
  "AllSamples_Zonation_Coefficients.csv"
))
if (any(!file.exists(required_inputs))) {
  stop("Required Continuous analysis file(s) missing: ",
       paste(basename(required_inputs[!file.exists(required_inputs)]),
             collapse = ", "), call. = FALSE)
}

dir.create(BETA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CURVE_DIR, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(BETA_DIR) || !dir.exists(CURVE_DIR)) {
  stop("Could not create cache folders: ", CACHE_DIR, call. = FALSE)
}
CACHE_DIR <- normalizePath(CACHE_DIR, winslash = "/", mustWork = TRUE)
BETA_DIR <- normalizePath(BETA_DIR, winslash = "/", mustWork = TRUE)
CURVE_DIR <- normalizePath(CURVE_DIR, winslash = "/", mustWork = TRUE)

############################################################
# LOAD ANALYSIS DATA ONCE
############################################################

beta_runner <- run_zonation_coefficient_beta
curve_runner <- run_selected_trajectory
if (!is.function(beta_runner) || !is.function(curve_runner)) {
  missing_runners <- c(
    if (!is.function(beta_runner)) "run_zonation_coefficient_beta",
    if (!is.function(curve_runner)) "run_selected_trajectory"
  )
  stop(
    "Required analysis function(s) are missing: ",
    paste(missing_runners, collapse = ", "),
    ". Re-download the standalone cache builder.",
    call. = FALSE
  )
}

# Newer analysis scripts expose preparation functions so the large Continuous
# tables are read only once.  Older scripts can still be used: the builder
# falls back to their original per-protein readers instead of stopping with a
# cache-version error.
beta_prepare <- prepare_zonation_coefficient_data
curve_prepare <- prepare_selected_trajectory_data
beta_supports_prepared <- is.function(beta_prepare) &&
  "prepared_data" %in% names(formals(beta_runner))
curve_supports_prepared <- is.function(curve_prepare) &&
  "prepared" %in% names(formals(curve_runner))

if (beta_supports_prepared && curve_supports_prepared) {
  cat("Reading Continuous analysis tables once (fast cache mode)...\n")
} else {
  warning(
    "Legacy analysis script detected. Cache creation will continue in ",
    "compatibility mode, but may be slower. Replacing both analysis scripts ",
    "with the versions supplied in the app package enables fast cache mode.",
    call. = FALSE
  )
}

beta_data <- if (beta_supports_prepared) {
  beta_prepare(CONTINUOUS_OUTPUT_DIR)
} else NULL

trajectory_prepared <- if (curve_supports_prepared) {
  curve_prepare(CONTINUOUS_OUTPUT_DIR)
} else NULL

trajectory_data <- if (!is.null(trajectory_prepared)) {
  trajectory_prepared$data
} else {
  fallback_file <- file.path(
    CONTINUOUS_OUTPUT_DIR, "Spatial_Proteomics_AnalysisReady.csv"
  )
  fallback_data <- utils::read.csv(
    fallback_file, check.names = FALSE, stringsAsFactors = FALSE,
    fileEncoding = "UTF-8"
  )
  fallback_required <- c("Protein.Group", "Protein")
  fallback_missing <- setdiff(fallback_required, names(fallback_data))
  if (length(fallback_missing)) {
    stop(
      basename(fallback_file), " is missing column(s): ",
      paste(fallback_missing, collapse = ", "), call. = FALSE
    )
  }
  fallback_data$Protein.Group <- trimws(as.character(fallback_data$Protein.Group))
  fallback_data$Protein <- trimws(as.character(fallback_data$Protein))
  fallback_data
}

protein_ids <- unique(trimws(as.character(trajectory_data$Protein.Group)))
protein_ids <- protein_ids[!is.na(protein_ids) & nzchar(protein_ids)]
if (length(TARGET_PROTEIN_GROUPS)) {
  requested <- unique(trimws(as.character(TARGET_PROTEIN_GROUPS)))
  requested <- requested[!is.na(requested) & nzchar(requested)]
  missing_requested <- setdiff(toupper(requested), toupper(protein_ids))
  if (length(missing_requested)) {
    warning("Requested Protein.Group value(s) not found: ",
            paste(missing_requested, collapse = ", "), call. = FALSE)
  }
  protein_ids <- protein_ids[toupper(protein_ids) %in% toupper(requested)]
}
protein_ids <- sort(protein_ids)
if (!length(protein_ids)) stop("No Protein.Group entries to cache.", call. = FALSE)

first_label <- function(id) {
  values <- unique(trajectory_data$Protein[trajectory_data$Protein.Group == id])
  values <- trimws(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values)) values[[1]] else id
}
catalog <- data.frame(
  Protein.Group = protein_ids,
  Protein = vapply(protein_ids, first_label, character(1)),
  stringsAsFactors = FALSE
)

############################################################
# RESUMABLE CACHE INDEX
############################################################

index_file <- file.path(CACHE_DIR, "Cache_Index.csv")
index_columns <- c(
  "Protein.Group", "Protein", "Beta_Status", "Trajectory_Status",
  "Beta_PNG", "Beta_PDF", "Trajectory_PNG", "Trajectory_PDF",
  "Message", "Updated"
)
index <- data.frame(
  Protein.Group = catalog$Protein.Group,
  Protein = catalog$Protein,
  Beta_Status = "Pending",
  Trajectory_Status = "Pending",
  Beta_PNG = "", Beta_PDF = "", Trajectory_PNG = "", Trajectory_PDF = "",
  Message = "", Updated = "",
  stringsAsFactors = FALSE
)
if (file.exists(index_file)) {
  old <- tryCatch(
    utils::read.csv(index_file, check.names = FALSE,
                    stringsAsFactors = FALSE, fileEncoding = "UTF-8"),
    error = function(e) NULL
  )
  if (!is.null(old) && !length(setdiff(index_columns, names(old)))) {
    matched <- match(index$Protein.Group, old$Protein.Group)
    available <- which(!is.na(matched))
    copy_columns <- setdiff(index_columns, c("Protein.Group", "Protein"))
    index[available, copy_columns] <- old[matched[available], copy_columns]
  }
}

relative_to_cache <- function(path) {
  if (is.null(path) || !length(path) || is.na(path) || !nzchar(path)) return("")
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  prefix <- paste0(CACHE_DIR, "/")
  if (startsWith(path, prefix)) substring(path, nchar(prefix) + 1L) else path
}
absolute_from_cache <- function(path) {
  if (is.null(path) || !length(path) || is.na(path) || !nzchar(path)) return(NA_character_)
  path <- gsub("\\", "/", path, fixed = TRUE)
  if (grepl("^(/|[A-Za-z]:/)", path)) path else file.path(CACHE_DIR, path)
}
complete_files <- function(row, prefix) {
  fields <- if (prefix == "Beta") "Beta_PNG" else "Trajectory_PNG"
  paths <- vapply(row[fields], absolute_from_cache, character(1))
  all(!is.na(paths) & file.exists(paths) & file.size(paths) > 0)
}
write_index <- function() {
  utils::write.csv(index, index_file, row.names = FALSE, na = "",
                   fileEncoding = "UTF-8")
}

source_info <- data.frame(
  File = basename(required_inputs),
  Size_Bytes = file.size(required_inputs),
  Modified = format(file.info(required_inputs)$mtime, "%Y-%m-%d %H:%M:%S"),
  Continuous_Output_Dir = CONTINUOUS_OUTPUT_DIR,
  stringsAsFactors = FALSE
)
utils::write.csv(source_info, file.path(CACHE_DIR, "Cache_Source.csv"),
                 row.names = FALSE, na = "", fileEncoding = "UTF-8")
write_index()

############################################################
# BUILD ALL SINGLE-PROTEIN FIGURES
############################################################

cat("Cache folder:", CACHE_DIR, "\n")
cat("Protein.Group entries:", nrow(index), "\n")
cat("Existing complete entries will be skipped:", !OVERWRITE_EXISTING, "\n\n")

for (i in seq_len(nrow(index))) {
  protein_id <- index$Protein.Group[[i]]
  beta_done <- identical(index$Beta_Status[[i]], "Created") &&
    complete_files(index[i, , drop = FALSE], "Beta")
  curve_done <- identical(index$Trajectory_Status[[i]], "Created") &&
    complete_files(index[i, , drop = FALSE], "Trajectory")
  if (!OVERWRITE_EXISTING && beta_done && curve_done) {
    cat(sprintf("[%d/%d] %s - cached, skipped\n", i, nrow(index), protein_id))
    next
  }

  cat(sprintf("[%d/%d] %s (%s)\n", i, nrow(index), protein_id,
              index$Protein[[i]]))
  messages <- character()

  if (OVERWRITE_EXISTING || !beta_done) {
    beta_args <- list(
      input_dir = CONTINUOUS_OUTPUT_DIR,
      output_dir = BETA_DIR,
      target = protein_id,
      png_dpi = PNG_DPI
    )
    if (beta_supports_prepared) beta_args$prepared_data <- beta_data
    beta_result <- tryCatch(
      do.call(beta_runner, beta_args), error = function(e) e
    )
    if (inherits(beta_result, "error")) {
      index$Beta_Status[[i]] <- "Failed"
      messages <- c(messages, paste0("Beta: ", conditionMessage(beta_result)))
    } else {
      index$Beta_Status[[i]] <- "Created"
      index$Beta_PNG[[i]] <- relative_to_cache(beta_result$png_file)
      index$Beta_PDF[[i]] <- relative_to_cache(beta_result$pdf_file)
      if (!KEEP_AUXILIARY_CSV && file.exists(beta_result$table_file)) {
        unlink(beta_result$table_file)
      }
    }
  }

  if (OVERWRITE_EXISTING || !curve_done) {
    curve_args <- list(
      input_dir = CONTINUOUS_OUTPUT_DIR,
      output_dir = CURVE_DIR,
      target = protein_id,
      png_dpi = PNG_DPI
    )
    if (curve_supports_prepared) curve_args$prepared <- trajectory_prepared
    curve_result <- tryCatch(
      do.call(curve_runner, curve_args), error = function(e) e
    )
    if (inherits(curve_result, "error")) {
      index$Trajectory_Status[[i]] <- "Failed"
      messages <- c(messages, paste0("Trajectory: ", conditionMessage(curve_result)))
    } else {
      index$Trajectory_Status[[i]] <- "Created"
      index$Trajectory_PNG[[i]] <- relative_to_cache(curve_result$png_file)
      index$Trajectory_PDF[[i]] <- relative_to_cache(curve_result$pdf_file)
      if (!KEEP_AUXILIARY_CSV) {
        auxiliary <- c(curve_result$observations_file, curve_result$curves_file,
                       curve_result$audit_file)
        unlink(auxiliary[file.exists(auxiliary)])
      }
    }
  }

  index$Message[[i]] <- paste(messages, collapse = " | ")
  index$Updated[[i]] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  write_index()
  grDevices::graphics.off()
}

complete_both <- index$Beta_Status == "Created" &
  index$Trajectory_Status == "Created"
cat("\nCache build finished.\n")
cat("Both figures created:", sum(complete_both), "/", nrow(index), "\n")
cat("Failed or incomplete:", sum(!complete_both), "\n")
cat("Index:", index_file, "\n")
