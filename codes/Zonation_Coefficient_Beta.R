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
  pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
  table_file <- file.path(output_dir, paste0(stem, "_Matched.csv"))
  ggplot2::ggsave(png_file, p, width = 10, height = 8, dpi = png_dpi,
                  bg = "white")
  ggplot2::ggsave(pdf_file, p, width = 10, height = 8, bg = "white")
  utils::write.csv(highlight, table_file, row.names = FALSE, na = "",
                   fileEncoding = "UTF-8")

  created <- c(png_file, pdf_file, table_file)
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
    pdf_file = normalizePath(pdf_file, winslash = "/", mustWork = TRUE),
    table_file = normalizePath(table_file, winslash = "/", mustWork = TRUE),
    target = target,
    protein_group = paste(matched_ids, collapse = ", ")
  )
}
