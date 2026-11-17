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
  pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
  observations_file <- file.path(output_dir, paste0(stem, "_Observations.csv"))
  curves_file <- file.path(output_dir, paste0(stem, "_LOESS_Curves.csv"))
  audit_file <- file.path(output_dir, paste0(stem, "_Fit_Audit.csv"))
  ggplot2::ggsave(png_file, p, width = 8, height = 6, dpi = png_dpi, bg = "white")
  ggplot2::ggsave(pdf_file, p, width = 8, height = 6, bg = "white")
  utils::write.csv(plot_data, observations_file, row.names = FALSE, na = "",
                   fileEncoding = "UTF-8")
  utils::write.csv(curves, curves_file, row.names = FALSE, na = "",
                   fileEncoding = "UTF-8")
  utils::write.csv(audit, audit_file, row.names = FALSE, na = "",
                   fileEncoding = "UTF-8")

  created <- c(png_file, pdf_file, observations_file, curves_file, audit_file)
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
    pdf_file = normalizePath(pdf_file, winslash = "/", mustWork = TRUE),
    observations_file = normalizePath(observations_file, winslash = "/", mustWork = TRUE),
    curves_file = normalizePath(curves_file, winslash = "/", mustWork = TRUE),
    audit_file = normalizePath(audit_file, winslash = "/", mustWork = TRUE),
    target = target,
    protein_group = protein_id,
    protein = plot_label
  )
}
