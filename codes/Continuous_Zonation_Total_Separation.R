############################################################
# ==========================================================
# Mouse Liver Spatial Proteomics
#
# Continuous CV-PV Zonation Analysis
# Combined Total + GF/SPF Separate Figure 3-style workflow
#
# INPUT
#   1. report.parquet
#      - 8-set joint DIA-NN
#      - MBR ON
#      - low-quality ROIs already excluded before DIA-NN
#
#   2. ROI_metadata.csv
#      Required information:
#      Run
#      Image
#      Group
#      Mouse
#      Set
#      ROI
#      Spatial_S
#
#
# ANALYSIS
#
# report.parquet
#       ↓
# cRAP / CON__ identification and information export
#       ↓
# remove cRAP / CON__ proteins
#       ↓
# robust Run name normalization
#       ↓
# DIA-NN ↔ ROI metadata matching
#       ↓
# 70% protein completeness
#       ↓
# log2 abundance
#       ↓
# no imputation
#       ↓
# protein-wise RobustScaler (median / IQR)
#       ↓
# LMM:
#
# RobustScaled ~ Spatial_S + (1 | Mouse)
#
#       ↓
# Beta / P / BH-Q
#       ↓
# F3 style volcano
# continuous gradients
# canonical marker validation
#
#
# GF vs SPF:
#
# RobustScaled ~ Spatial_S * Group + (1 | Mouse)
#
############################################################


# Shiny/command-line entry point. Sourcing this file only defines the
# function; the analysis starts when run_continuous_zonation() is called.
run_continuous_zonation <- function(
  report_path,
  metadata_path,
  output_dir,
  completeness_cutoff = 0.70,
  q_cutoff = 0.05,
  beta_cutoff = 1,
  top_label_n = 15L,
  heatmap_protein_set = "all",
  heatmap_top_n = 100L,
  heatmap_n_bins = 20L,
  heatmap_min_mice_per_bin = 1L,
  heatmap_create_multipage = TRUE,
  heatmap_proteins_per_page = 100L
) {



############################################################
# 0. PACKAGES
############################################################

library(tidyverse)
library(arrow)

library(lme4)
library(lmerTest)

library(ggrepel)
library(patchwork)



############################################################
# 1. USER SETTINGS
############################################################


# ----------------------------------------------------------
# Input files
# ----------------------------------------------------------

PARQUET_FILE <- normalizePath(report_path, winslash = "/", mustWork = TRUE)

METADATA_FILE <- normalizePath(metadata_path, winslash = "/", mustWork = TRUE)


# ----------------------------------------------------------
# Output directory
# ----------------------------------------------------------

OUTPUT_DIR <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

dir.create(
  OUTPUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)


# ----------------------------------------------------------
# Analysis thresholds
# ----------------------------------------------------------

COMPLETENESS_CUTOFF <- as.numeric(completeness_cutoff)

Q_CUTOFF <- as.numeric(q_cutoff)

BETA_CUTOFF <- as.numeric(beta_cutoff)


# ----------------------------------------------------------
# Volcano labels
# ----------------------------------------------------------

TOP_LABEL_N <- as.integer(top_label_n)

if (any(!is.finite(c(COMPLETENESS_CUTOFF, Q_CUTOFF, BETA_CUTOFF))) ||
    COMPLETENESS_CUTOFF <= 0 || COMPLETENESS_CUTOFF > 1 ||
    Q_CUTOFF <= 0 || Q_CUTOFF > 1 || BETA_CUTOFF < 0 ||
    is.na(TOP_LABEL_N) || TOP_LABEL_N < 0L) {
  stop("Invalid continuous-analysis thresholds.", call. = FALSE)
}


# ----------------------------------------------------------
# Scaling method
# ----------------------------------------------------------
# RobustScaler equivalent used in this script:
#   (log2 abundance - protein median) / protein IQR
#
# Scaling is performed across all retained ROIs for each
# protein, before LMM fitting.




############################################################
# 2. LOAD DIA-NN report.parquet
############################################################

message(
  "\nLoading DIA-NN report.parquet ..."
)


diann <- read_parquet(
  PARQUET_FILE
)


cat(
  "\nDIA-NN rows:",
  nrow(diann),
  "\n"
)

cat(
  "DIA-NN columns:",
  ncol(diann),
  "\n"
)



############################################################
# 3. CHECK REQUIRED DIA-NN COLUMNS
############################################################

required_diann <- c(
  "Run",
  "Protein.Group"
)


missing_diann <- setdiff(
  required_diann,
  names(diann)
)


if (
  length(missing_diann) > 0
) {
  
  stop(
    paste0(
      "\nMissing DIA-NN columns: ",
      paste(
        missing_diann,
        collapse = ", "
      )
    )
  )
  
}



############################################################
# 4. SELECT PROTEIN ABUNDANCE COLUMN
############################################################

candidate_abundance <- c(
  "PG.MaxLFQ",
  "PG.Quantity",
  "Protein.Quantity"
)


available_abundance <- candidate_abundance[
  candidate_abundance %in%
    names(diann)
]


if (
  length(available_abundance) == 0
) {
  
  stop(
    paste0(
      "\nNo suitable protein abundance column was found.\n",
      "Please check names(diann).\n"
    )
  )
  
}


ABUNDANCE_COL <- available_abundance[1]


cat(
  "\nUsing abundance column:",
  ABUNDANCE_COL,
  "\n"
)



############################################################
# 5. SELECT PROTEIN LABEL
############################################################

if (
  "Genes" %in% names(diann)
) {
  
  PROTEIN_LABEL_COL <- "Genes"
  
} else if (
  "Protein.Names" %in% names(diann)
) {
  
  PROTEIN_LABEL_COL <- "Protein.Names"
  
} else {
  
  PROTEIN_LABEL_COL <- "Protein.Group"
  
}


cat(
  "Protein label column:",
  PROTEIN_LABEL_COL,
  "\n"
)



############################################################
# 5A. IDENTIFY, REPORT, AND REMOVE cRAP CONTAMINANTS
############################################################

# cRAP contaminants are removed before protein_long is
# created. Therefore, they cannot enter completeness,
# scaling, LMM, GF/SPF comparison, or heatmap calculations.
#
# Available protein annotation columns are checked for an
# explicit cRAP prefix or the commonly used DIA-NN
# contaminant prefix CON__.

crap_search_columns <- intersect(
  c(
    "Protein.Group",
    "Protein.Ids",
    "Protein.Names",
    "Genes"
  ),
  names(diann)
)


if (
  length(crap_search_columns) == 0
) {
  
  stop(
    "\nNo protein annotation column was available for cRAP detection."
  )
  
}


crap_pattern <- regex(
  paste0(
    "(^|[;|[:space:]])",
    "(?:REV__)?",
    "(?:cRAP(?:$|[_:\\-|0-9])|CON__)"
  ),
  ignore_case = TRUE
)


is_crap_row <- rep(
  FALSE,
  nrow(diann)
)


for (
  column_name in crap_search_columns
) {
  
  annotation_text <- replace_na(
    as.character(
      diann[[column_name]]
    ),
    ""
  )
  
  is_crap_row <-
    is_crap_row |
    str_detect(
      annotation_text,
      crap_pattern
    )
  
}


diann <- diann %>%
  mutate(
    Is_cRAP = is_crap_row
  )



############################################################
# 5B. SAVE cRAP INFORMATION BEFORE REMOVAL
############################################################

crap_information <- diann %>%
  filter(Is_cRAP) %>%
  transmute(
    Run = as.character(Run),
    Protein.Group = as.character(Protein.Group),
    Protein = as.character(
      .data[[PROTEIN_LABEL_COL]]
    ),
    Abundance = as.numeric(
      .data[[ABUNDANCE_COL]]
    )
  )


crap_summary <- tibble(
  Metric = c(
    "DIA-NN rows before cRAP removal",
    "Rows identified as cRAP or CON__",
    "Unique cRAP or CON__ protein groups",
    "Runs containing cRAP or CON__",
    "DIA-NN rows after cRAP removal"
  ),
  Value = c(
    nrow(diann),
    nrow(crap_information),
    n_distinct(crap_information$Protein.Group),
    n_distinct(crap_information$Run),
    nrow(diann) - nrow(crap_information)
  )
)


crap_protein_list <- crap_information %>%
  group_by(
    Protein.Group,
    Protein
  ) %>%
  summarise(
    N_Run = n_distinct(Run),
    N_Detected = sum(
      !is.na(Abundance) &
        Abundance > 0
    ),
    Median_Abundance = if (
      all(is.na(Abundance))
    ) {
      NA_real_
    } else {
      median(
        Abundance,
        na.rm = TRUE
      )
    },
    .groups = "drop"
  ) %>%
  arrange(
    desc(N_Run),
    Protein.Group
  )


crap_run_summary <- crap_information %>%
  group_by(Run) %>%
  summarise(
    N_cRAP_Protein_Group = n_distinct(Protein.Group),
    N_cRAP_Rows = n(),
    .groups = "drop"
  ) %>%
  arrange(Run)


write_csv(
  crap_summary,
  file.path(
    OUTPUT_DIR,
    "cRAP_Removal_Summary.csv"
  )
)


write_csv(
  crap_protein_list,
  file.path(
    OUTPUT_DIR,
    "cRAP_Protein_List_Before_Removal.csv"
  )
)


write_csv(
  crap_run_summary,
  file.path(
    OUTPUT_DIR,
    "cRAP_Run_Summary_Before_Removal.csv"
  )
)


cat(
  "\n=============================================\n"
)

cat(
  "cRAP / CONTAMINANT CHECK BEFORE ANALYSIS\n"
)

cat(
  "=============================================\n"
)

print(crap_summary)



############################################################
# 5C. REMOVE cRAP BEFORE ALL DOWNSTREAM ANALYSES
############################################################

diann <- diann %>%
  filter(!Is_cRAP) %>%
  select(-Is_cRAP)


cat(
  "\ncRAP removal completed. Remaining DIA-NN rows:",
  nrow(diann),
  "\n"
)



############################################################
# 6. EXTRACT PROTEIN-LEVEL DATA
############################################################

protein_long <- diann %>%
  
  transmute(
    
    Run =
      as.character(Run),
    
    Protein.Group =
      as.character(Protein.Group),
    
    Protein =
      as.character(
        .data[[PROTEIN_LABEL_COL]]
      ),
    
    Abundance =
      as.numeric(
        .data[[ABUNDANCE_COL]]
      )
    
  ) %>%
  
  filter(
    
    !is.na(Run),
    
    !is.na(Protein.Group)
    
  )



############################################################
# 7. COLLAPSE DUPLICATE Protein × Run
############################################################

protein_long <- protein_long %>%
  
  group_by(
    Run,
    Protein.Group,
    Protein
  ) %>%
  
  summarise(
    
    Abundance = {
      
      x <- Abundance
      
      if (
        all(is.na(x))
      ) {
        
        NA_real_
        
      } else {
        
        max(
          x,
          na.rm = TRUE
        )
        
      }
      
    },
    
    .groups = "drop"
    
  )



############################################################
# 8. LOAD ROI METADATA
############################################################

metadata <- read_csv(
  METADATA_FILE,
  show_col_types = FALSE
)



############################################################
# 9. CHECK METADATA COLUMNS
############################################################

required_metadata <- c(
  
  "Run",
  "Image",
  "Group",
  "Mouse",
  "Set",
  "ROI",
  "Spatial_S"
  
)


missing_metadata <- setdiff(
  required_metadata,
  names(metadata)
)


if (
  length(missing_metadata) > 0
) {
  
  stop(
    paste0(
      "\nMissing ROI metadata columns: ",
      paste(
        missing_metadata,
        collapse = ", "
      )
    )
  )
  
}



############################################################
# 10. FORMAT ROI METADATA
############################################################

metadata <- metadata %>%
  
  mutate(
    
    Run =
      as.character(Run),
    
    Image =
      as.character(Image),
    
    Group =
      factor(
        Group,
        levels = c(
          "GF",
          "SPF"
        )
      ),
    
    Mouse =
      factor(Mouse),
    
    Set =
      factor(Set),
    
    ROI =
      as.character(ROI),
    
    Spatial_S =
      as.numeric(Spatial_S)
    
  )



############################################################
# 11. CHECK Spatial_S RANGE
############################################################

if (
  any(
    metadata$Spatial_S < 0 |
    metadata$Spatial_S > 1,
    na.rm = TRUE
  )
) {
  
  stop(
    "\nSpatial_S contains values outside 0-1."
  )
  
}



############################################################
# ==========================================================
#
# ROBUST DIA-NN / ROI METADATA RUN MATCHING
#
# ==========================================================
############################################################



############################################################
# 12. FUNCTION TO NORMALIZE RUN NAME
#
# Examples:
#
# GF1-2_ROI1
# GF1-2_ROI01
# GF1-2-ROI1
# GF1-2_ROI_01
# 20260904_Mouse_Liver_GF1-2_ROI01.raw
# X:/folder/GF1-2_ROI01.raw
#
# are converted to:
#
# GF1-2_ROI01
#
############################################################

normalize_run_key <- function(x) {
  
  x <- as.character(x)
  
  
  # --------------------------------------------------------
  # Convert Windows path to generic path
  # --------------------------------------------------------
  
  x <- gsub(
    "\\\\",
    "/",
    x
  )
  
  
  # --------------------------------------------------------
  # Keep file name only
  # --------------------------------------------------------
  
  x <- basename(x)
  
  
  # --------------------------------------------------------
  # Remove common extensions
  # --------------------------------------------------------
  
  x <- str_remove(
    x,
    regex(
      "\\.(raw|mzml|d)$",
      ignore_case = TRUE
    )
  )
  
  
  # --------------------------------------------------------
  # Extract Set
  #
  # GF1-2
  # SPF1-4
  # --------------------------------------------------------
  
  set_part <- str_extract(
    
    x,
    
    "(GF|SPF)[0-9]+-[0-9]+"
    
  )
  
  
  # --------------------------------------------------------
  # Extract ROI number
  #
  # ROI1
  # ROI01
  # ROI_01
  # ROI-01
  # --------------------------------------------------------
  
  roi_match <- str_match(
    
    x,
    
    regex(
      "ROI[_-]?([0-9]+)",
      ignore_case = TRUE
    )
    
  )
  
  
  roi_num <- suppressWarnings(
    as.integer(
      roi_match[, 2]
    )
  )
  
  
  # --------------------------------------------------------
  # Return normalized key
  # --------------------------------------------------------
  
  ifelse(
    
    !is.na(set_part) &
      !is.na(roi_num),
    
    paste0(
      
      set_part,
      
      "_ROI",
      
      sprintf(
        "%02d",
        roi_num
      )
      
    ),
    
    NA_character_
    
  )
  
}



############################################################
# 13. CREATE DIA-NN RUN KEY
############################################################

protein_long <- protein_long %>%
  
  mutate(
    
    Run_original =
      Run,
    
    Run_Key =
      normalize_run_key(
        Run
      )
    
  )



############################################################
# 14. CREATE METADATA RUN KEY
############################################################

metadata <- metadata %>%
  
  mutate(
    
    Metadata_Run =
      Run,
    
    Run_Key =
      normalize_run_key(
        Run
      )
    
  )



############################################################
# 15. CHECK PARSING
############################################################

cat(
  "\n=============================================\n"
)

cat(
  "RUN NAME PARSING CHECK\n"
)

cat(
  "=============================================\n"
)


cat(
  "\nDIA-NN rows with unparsed Run:",
  sum(
    is.na(
      protein_long$Run_Key
    )
  ),
  "\n"
)


cat(
  "Metadata rows with unparsed Run:",
  sum(
    is.na(
      metadata$Run_Key
    )
  ),
  "\n"
)



############################################################
# 16. UNIQUE RUN KEYS
############################################################

diann_keys <- unique(
  na.omit(
    protein_long$Run_Key
  )
)


metadata_keys <- unique(
  na.omit(
    metadata$Run_Key
  )
)


matched_keys <- intersect(
  diann_keys,
  metadata_keys
)


cat(
  "\nUnique DIA-NN keys:",
  length(
    diann_keys
  ),
  "\n"
)


cat(
  "Unique metadata keys:",
  length(
    metadata_keys
  ),
  "\n"
)


cat(
  "Matched keys:",
  length(
    matched_keys
  ),
  "\n"
)



############################################################
# 17. SHOW EXAMPLE RUN NAMES
############################################################

cat(
  "\n--- DIA-NN Run examples ---\n"
)

print(
  
  head(
    
    unique(
      protein_long %>%
        select(
          Run_original,
          Run_Key
        )
    ),
    
    20
    
  )
  
)


cat(
  "\n--- Metadata Run examples ---\n"
)

print(
  
  head(
    
    unique(
      metadata %>%
        select(
          Metadata_Run,
          Run_Key
        )
    ),
    
    20
    
  )
  
)



############################################################
# 18. UNMATCHED DIA-NN RUNS
############################################################

unmatched_diann <- setdiff(
  diann_keys,
  metadata_keys
)


if (
  length(
    unmatched_diann
  ) > 0
) {
  
  cat(
    "\nDIA-NN Run keys absent from metadata:\n"
  )
  
  print(
    unmatched_diann
  )
  
}



############################################################
# 19. UNMATCHED METADATA RUNS
#
# This can be NORMAL because ROI_metadata.csv contains
# all 96 ROI coordinates, while low-quality ROIs were
# removed before joint DIA-NN.
############################################################

unmatched_metadata <- setdiff(
  metadata_keys,
  diann_keys
)


if (
  length(
    unmatched_metadata
  ) > 0
) {
  
  cat(
    "\nMetadata Run keys absent from DIA-NN:\n"
  )
  
  print(
    unmatched_metadata
  )
  
}



############################################################
# 20. SAVE RUN MATCH REPORT
############################################################

run_match_report <- tibble(
  
  Run_Key =
    union(
      diann_keys,
      metadata_keys
    )
  
) %>%
  
  mutate(
    
    In_DIANN =
      Run_Key %in%
      diann_keys,
    
    In_Metadata =
      Run_Key %in%
      metadata_keys,
    
    Matched =
      In_DIANN &
      In_Metadata
    
  )


write_csv(
  
  run_match_report,
  
  file.path(
    OUTPUT_DIR,
    "Run_Matching_Report.csv"
  )
  
)



############################################################
# 21. STOP IF ZERO MATCHES
############################################################

if (
  length(
    matched_keys
  ) == 0
) {
  
  stop(
    paste0(
      "\nNo normalized DIA-NN Runs matched ROI_metadata.csv.\n\n",
      "Please inspect:\n",
      "1. unique(diann$Run)\n",
      "2. metadata$Run\n",
      "3. Run_Matching_Report.csv\n"
    )
  )
  
}



############################################################
# 22. PREPARE METADATA FOR JOIN
############################################################

metadata_join <- metadata %>%
  
  select(
    
    Run_Key,
    
    Metadata_Run,
    
    Image,
    Group,
    Mouse,
    Set,
    ROI,
    Spatial_S,
    
    everything(),
    
    -Run
    
  )



############################################################
# 23. JOIN DIA-NN WITH ROI METADATA
############################################################

protein_long <- protein_long %>%
  
  inner_join(
    
    metadata_join,
    
    by =
      "Run_Key"
    
  )



############################################################
# 24. CHECK MATCHED QC-RETAINED ROIs
############################################################

cat(
  "\nSuccessfully matched QC-retained ROIs:",
  n_distinct(
    protein_long$Run_Key
  ),
  "\n"
)


if (
  n_distinct(
    protein_long$Run_Key
  ) == 0
) {
  
  stop(
    "\nNo ROI remained after metadata join."
  )
  
}



############################################################
# 25. SAVE ACTUAL QC-RETAINED ROI LIST
############################################################

retained_roi_list <- protein_long %>%
  
  distinct(
    
    Run_Key,
    
    Run_original,
    
    Metadata_Run,
    
    Image,
    Group,
    Mouse,
    Set,
    ROI,
    Spatial_S
    
  ) %>%
  
  arrange(
    Group,
    Set,
    ROI
  )


write_csv(
  
  retained_roi_list,
  
  file.path(
    OUTPUT_DIR,
    "QC_Retained_ROI_List.csv"
  )
  
)



############################################################
# 26. ROI COUNT PER SET
############################################################

roi_summary <- retained_roi_list %>%
  
  count(
    
    Group,
    Set,
    
    name =
      "N_ROI"
    
  )


print(
  roi_summary
)


write_csv(
  
  roi_summary,
  
  file.path(
    OUTPUT_DIR,
    "ROI_Count_Per_Set.csv"
  )
  
)



############################################################
# ==========================================================
#
# PROTEIN COMPLETENESS
#
# ==========================================================
############################################################



############################################################
# 27. TOTAL QC-RETAINED ROI NUMBER
############################################################

n_roi_total <- n_distinct(
  protein_long$Run_Key
)


cat(
  "\nTotal QC-retained ROIs:",
  n_roi_total,
  "\n"
)



############################################################
# 28. OVERALL PROTEIN COMPLETENESS
############################################################

completeness_table <- protein_long %>%
  
  group_by(
    Protein.Group,
    Protein
  ) %>%
  
  summarise(
    
    N_detected =
      sum(
        
        !is.na(
          Abundance
        ) &
          
          Abundance > 0
        
      ),
    
    .groups =
      "drop"
    
  ) %>%
  
  mutate(
    
    Total_ROI =
      n_roi_total,
    
    Completeness =
      N_detected /
      Total_ROI
    
  )



############################################################
# 29. APPLY 70% COMPLETENESS
############################################################

retained_proteins <- completeness_table %>%
  
  filter(
    
    Completeness >=
      COMPLETENESS_CUTOFF
    
  )


cat(
  "\nProteins before completeness filtering:",
  nrow(
    completeness_table
  ),
  "\n"
)


cat(
  "Proteins after >=70% completeness:",
  nrow(
    retained_proteins
  ),
  "\n"
)



############################################################
# 30. SAVE COMPLETENESS TABLES
############################################################

write_csv(
  
  completeness_table,
  
  file.path(
    OUTPUT_DIR,
    "Protein_Completeness_All.csv"
  )
  
)


write_csv(
  
  retained_proteins,
  
  file.path(
    OUTPUT_DIR,
    "Protein_Completeness_Retained70.csv"
  )
  
)



############################################################
# 31. KEEP ≥70% COMPLETE PROTEINS
############################################################

dat <- protein_long %>%
  
  semi_join(
    
    retained_proteins,
    
    by = c(
      "Protein.Group",
      "Protein"
    )
    
  )



############################################################
# 32. LOG2 TRANSFORMATION
#
# NO IMPUTATION
############################################################

dat <- dat %>%
  
  mutate(
    
    Log2Abundance =
      case_when(
        
        !is.na(
          Abundance
        ) &
          Abundance > 0
        ~ log2(
          Abundance
        ),
        
        TRUE
        ~ NA_real_
        
      )
    
  )



############################################################
# 33. PROTEIN-WISE ROBUST SCALING
#
# Paper-matched RobustScaler concept:
#   RobustScaled = (Log2Abundance - median) / IQR
#
# This replaces the previous mean/SD Z-score normalization.
# Scaling is performed protein-wise across all retained ROIs.
# No imputation is performed.
############################################################

dat <- dat %>%
  
  group_by(
    Protein.Group,
    Protein
  ) %>%
  
  mutate(
    
    Protein_Median =
      median(
        Log2Abundance,
        na.rm = TRUE
      ),
    
    Protein_IQR =
      IQR(
        Log2Abundance,
        na.rm = TRUE,
        type = 7
      ),
    
    RobustScaled =
      if_else(
        
        is.finite(
          Protein_IQR
        ) &
          Protein_IQR > 0,
        
        (
          Log2Abundance -
            Protein_Median
        ) /
          Protein_IQR,
        
        NA_real_
        
      )
    
  ) %>%
  
  ungroup()



############################################################
# 34. REMOVE UNDEFINED MODEL VALUES
############################################################

dat <- dat %>%
  
  filter(
    
    is.finite(
      RobustScaled
    ),
    
    !is.na(
      Spatial_S
    ),
    
    !is.na(
      Mouse
    )
    
  )



############################################################
# 35. SAVE ANALYSIS-READY DATA
############################################################

write_csv(
  
  dat,
  
  file.path(
    OUTPUT_DIR,
    "Spatial_Proteomics_AnalysisReady.csv"
  )
  
)



############################################################
# ==========================================================
#
# METADATA / DESIGN CHECK
#
# ==========================================================
############################################################



############################################################
# 36. MOUSE STRUCTURE CHECK
############################################################

mouse_summary <- dat %>%
  
  distinct(
    
    Run_Key,
    Group,
    Mouse,
    Set
    
  ) %>%
  
  count(
    
    Group,
    Mouse,
    
    name =
      "N_ROI"
    
  )


cat(
  "\nMouse structure:\n"
)


print(
  mouse_summary
)


write_csv(
  
  mouse_summary,
  
  file.path(
    OUTPUT_DIR,
    "Mouse_Structure_Check.csv"
  )
  
)



############################################################
# 37. SPATIAL SCORE SUMMARY
############################################################

spatial_summary <- retained_roi_list %>%
  
  group_by(
    
    Image,
    Group,
    Set
    
  ) %>%
  
  summarise(
    
    N_ROI =
      n(),
    
    Min_S =
      min(
        Spatial_S,
        na.rm = TRUE
      ),
    
    Q1_S =
      quantile(
        Spatial_S,
        0.25,
        na.rm = TRUE
      ),
    
    Median_S =
      median(
        Spatial_S,
        na.rm = TRUE
      ),
    
    Q3_S =
      quantile(
        Spatial_S,
        0.75,
        na.rm = TRUE
      ),
    
    Max_S =
      max(
        Spatial_S,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
    
  )


write_csv(
  
  spatial_summary,
  
  file.path(
    OUTPUT_DIR,
    "Spatial_S_Summary.csv"
  )
  
)



############################################################
# 38. SPATIAL SCORE DISTRIBUTION
############################################################

p_spatial_distribution <-
  
  retained_roi_list %>%
  
  ggplot(
    
    aes(
      x = Spatial_S
    )
    
  ) +
  
  geom_histogram(
    bins = 20
  ) +
  
  facet_wrap(
    ~ Image,
    scales = "free_y"
  ) +
  
  scale_x_continuous(
    
    limits = c(
      0,
      1
    )
    
  ) +
  
  labs(
    
    title =
      "Distribution of CV-PV spatial coordinates",
    
    x =
      "Spatial_S (CV = 0 → PV = 1)",
    
    y =
      "Number of QC-retained ROIs"
    
  ) +
  
  theme_classic()


ggsave(
  
  file.path(
    OUTPUT_DIR,
    "Spatial_S_Distribution.pdf"
  ),
  
  p_spatial_distribution,
  
  width = 12,
  height = 8
  
)



############################################################
# ==========================================================
#
# PART I
#
# CONTINUOUS ZONATION LMM
#
# RobustScaled ~ Spatial_S + (1 | Mouse)
#
# ==========================================================
############################################################



############################################################
# 39. ROBUST ZONATION LMM FUNCTION
############################################################

fit_zonation_lmm <- function(df) {
  
  tryCatch({
    
    ########################################################
    # Clean observations
    ########################################################
    
    df_model <- df %>%
      
      filter(
        
        !is.na(
          RobustScaled
        ),
        
        !is.na(
          Spatial_S
        ),
        
        !is.na(
          Mouse
        )
        
      )
    
    
    ########################################################
    # Minimum observation requirement
    ########################################################
    
    if (
      nrow(
        df_model
      ) < 5
    ) {
      
      return(
        
        tibble(
          
          Beta =
            NA_real_,
          
          SE =
            NA_real_,
          
          DF =
            NA_real_,
          
          Statistic =
            NA_real_,
          
          P_value =
            NA_real_,
          
          N =
            nrow(
              df_model
            ),
          
          N_Mouse =
            n_distinct(
              df_model$Mouse
            ),
          
          Singular =
            NA
          
        )
        
      )
      
    }
    
    
    ########################################################
    # Need spatial variation
    ########################################################
    
    if (
      n_distinct(
        df_model$Spatial_S
      ) < 3
    ) {
      
      return(
        
        tibble(
          
          Beta =
            NA_real_,
          
          SE =
            NA_real_,
          
          DF =
            NA_real_,
          
          Statistic =
            NA_real_,
          
          P_value =
            NA_real_,
          
          N =
            nrow(
              df_model
            ),
          
          N_Mouse =
            n_distinct(
              df_model$Mouse
            ),
          
          Singular =
            NA
          
        )
        
      )
      
    }
    
    
    ########################################################
    # Need ≥2 Mouse levels
    ########################################################
    
    if (
      n_distinct(
        df_model$Mouse
      ) < 2
    ) {
      
      return(
        
        tibble(
          
          Beta =
            NA_real_,
          
          SE =
            NA_real_,
          
          DF =
            NA_real_,
          
          Statistic =
            NA_real_,
          
          P_value =
            NA_real_,
          
          N =
            nrow(
              df_model
            ),
          
          N_Mouse =
            n_distinct(
              df_model$Mouse
            ),
          
          Singular =
            NA
          
        )
        
      )
      
    }
    
    
    ########################################################
    # Fit model
    ########################################################
    
    model <- lmerTest::lmer(
      
      RobustScaled ~
        Spatial_S +
        (1 | Mouse),
      
      data =
        df_model,
      
      REML =
        FALSE
      
    )
    
    
    ########################################################
    # Coefficient table
    ########################################################
    
    coef_table <-
      summary(
        model
      )$coefficients
    
    
    ########################################################
    # Check coefficient
    ########################################################
    
    if (
      !"Spatial_S" %in%
      rownames(
        coef_table
      )
    ) {
      
      stop(
        "Spatial_S coefficient not found."
      )
      
    }
    
    
    ########################################################
    # P-value column
    ########################################################
    
    p_col <- grep(
      
      "^Pr",
      
      colnames(
        coef_table
      ),
      
      value =
        TRUE
      
    )
    
    
    ########################################################
    # DF column
    ########################################################
    
    df_col <- grep(
      
      "^df$",
      
      colnames(
        coef_table
      ),
      
      value =
        TRUE,
      
      ignore.case =
        TRUE
      
    )
    
    
    ########################################################
    # Extract coefficient
    ########################################################
    
    beta <-
      coef_table[
        "Spatial_S",
        "Estimate"
      ]
    
    
    se <-
      coef_table[
        "Spatial_S",
        "Std. Error"
      ]
    
    
    statistic <-
      coef_table[
        "Spatial_S",
        "t value"
      ]
    
    
    if (
      length(
        p_col
      ) >= 1
    ) {
      
      pvalue <-
        coef_table[
          "Spatial_S",
          p_col[1]
        ]
      
    } else {
      
      pvalue <-
        NA_real_
      
    }
    
    
    if (
      length(
        df_col
      ) >= 1
    ) {
      
      df_value <-
        coef_table[
          "Spatial_S",
          df_col[1]
        ]
      
    } else {
      
      df_value <-
        NA_real_
      
    }
    
    
    ########################################################
    # Singular model
    ########################################################
    
    singular_flag <-
      lme4::isSingular(
        model,
        tol = 1e-4
      )
    
    
    ########################################################
    # Return
    ########################################################
    
    tibble(
      
      Beta =
        as.numeric(
          beta
        ),
      
      SE =
        as.numeric(
          se
        ),
      
      DF =
        as.numeric(
          df_value
        ),
      
      Statistic =
        as.numeric(
          statistic
        ),
      
      P_value =
        as.numeric(
          pvalue
        ),
      
      N =
        nrow(
          df_model
        ),
      
      N_Mouse =
        n_distinct(
          df_model$Mouse
        ),
      
      Singular =
        singular_flag
      
    )
    
  },
  
  
  error = function(e) {
    
    tibble(
      
      Beta =
        NA_real_,
      
      SE =
        NA_real_,
      
      DF =
        NA_real_,
      
      Statistic =
        NA_real_,
      
      P_value =
        NA_real_,
      
      N =
        nrow(
          df
        ),
      
      N_Mouse =
        n_distinct(
          df$Mouse
        ),
      
      Singular =
        NA
      
    )
    
  })
  
}



############################################################
# 40. RUN LMM FOR EVERY PROTEIN
############################################################

message(
  "\nRunning continuous CV-PV zonation models..."
)


zonation_results <- dat %>%
  
  group_by(
    Protein.Group,
    Protein
  ) %>%
  
  nest() %>%
  
  mutate(
    
    Model_Result =
      map(
        data,
        fit_zonation_lmm
      )
    
  ) %>%
  
  select(
    -data
  ) %>%
  
  unnest(
    Model_Result
  ) %>%
  
  ungroup()



############################################################
# 41. VERIFY P_value COLUMN
############################################################

cat(
  "\nLMM output columns:\n"
)


print(
  names(
    zonation_results
  )
)


if (
  !"P_value" %in%
  names(
    zonation_results
  )
) {
  
  stop(
    "\nP_value column was not generated."
  )
  
}



############################################################
# 42. MODEL SUCCESS SUMMARY
############################################################

model_summary <- zonation_results %>%
  
  summarise(
    
    Total_Proteins =
      n(),
    
    Successful_Models =
      sum(
        !is.na(
          P_value
        )
      ),
    
    Failed_Models =
      sum(
        is.na(
          P_value
        )
      ),
    
    Singular_Models =
      sum(
        Singular == TRUE,
        na.rm = TRUE
      )
    
  )


cat(
  "\nLMM model summary:\n"
)


print(
  model_summary
)


write_csv(
  
  model_summary,
  
  file.path(
    OUTPUT_DIR,
    "LMM_Model_Summary.csv"
  )
  
)



############################################################
# 43. BH MULTIPLE TESTING CORRECTION
############################################################

zonation_results <- zonation_results %>%
  
  mutate(
    
    Q_value =
      p.adjust(
        P_value,
        method = "BH"
      ),
    
    MinusLog10Q =
      case_when(
        
        !is.na(
          Q_value
        )
        ~ -log10(
          pmax(
            Q_value,
            1e-300
          )
        ),
        
        TRUE
        ~ NA_real_
        
      ),
    
    Zonation =
      case_when(
        
        !is.na(
          Q_value
        ) &
          
          Q_value <
          Q_CUTOFF &
          
          Beta <=
          -BETA_CUTOFF
        ~ "Central",
        
        !is.na(
          Q_value
        ) &
          
          Q_value <
          Q_CUTOFF &
          
          Beta >=
          BETA_CUTOFF
        ~ "Portal",
        
        TRUE
        ~ "NS"
        
      )
    
  )



############################################################
# 44. SAVE ZONATION RESULTS
############################################################

write_csv(
  
  zonation_results,
  
  file.path(
    OUTPUT_DIR,
    "AllSamples_Zonation_Coefficients.csv"
  )
  
)



############################################################
# 45. ZONATION SUMMARY
############################################################

zonation_summary <- zonation_results %>%
  
  count(
    Zonation
  )


print(
  zonation_summary
)


write_csv(
  
  zonation_summary,
  
  file.path(
    OUTPUT_DIR,
    "Zonation_Summary.csv"
  )
  
)



############################################################
# ==========================================================
#
# F3 STYLE
# ZONATION VOLCANO
#
# ==========================================================
############################################################



############################################################
# 46. LABEL STRONGLY ZONATED PROTEINS
############################################################

top_central <- zonation_results %>%
  
  filter(
    Zonation ==
      "Central"
  ) %>%
  
  arrange(
    Q_value,
    Beta
  ) %>%
  
  slice_head(
    n =
      TOP_LABEL_N
  )


top_portal <- zonation_results %>%
  
  filter(
    Zonation ==
      "Portal"
  ) %>%
  
  arrange(
    Q_value,
    desc(
      Beta
    )
  ) %>%
  
  slice_head(
    n =
      TOP_LABEL_N
  )


top_labels <- bind_rows(
  
  top_central,
  
  top_portal
  
)



############################################################
# 47. VOLCANO PLOT
############################################################

p_volcano <- ggplot(
  
  zonation_results,
  
  aes(
    
    x =
      Beta,
    
    y =
      MinusLog10Q,
    
    color =
      Zonation
    
  )
  
) +
  
  geom_point(
    
    size =
      1.7,
    
    alpha =
      0.65,
    
    na.rm =
      TRUE
    
  ) +
  
  geom_vline(
    
    xintercept = c(
      -BETA_CUTOFF,
      BETA_CUTOFF
    ),
    
    linetype =
      "dashed"
    
  ) +
  
  geom_hline(
    
    yintercept =
      -log10(
        Q_CUTOFF
      ),
    
    linetype =
      "dashed"
    
  ) +
  
  geom_text_repel(
    
    data =
      top_labels,
    
    aes(
      label =
        Protein
    ),
    
    size =
      3,
    
    max.overlaps =
      Inf
    
  ) +
  
  labs(
    
    title =
      "Continuous liver protein zonation",
    
    subtitle =
      "QC-retained ROIs | 8-set joint DIA-NN MBR",
    
    x =
      "Zonation coefficient (Beta)\nCentral ←          → Portal",
    
    y =
      expression(
        -log[10](Q)
      ),
    
    color =
      "Zonation"
    
  ) +
  
  theme_classic(
    base_size = 13
  ) +
  
  theme(
    
    legend.position =
      "top"
    
  )



############################################################
# 48. SAVE VOLCANO
############################################################

ggsave(
  
  file.path(
    OUTPUT_DIR,
    "Zonation_Volcano.pdf"
  ),
  
  p_volcano,
  
  width = 9,
  height = 7
  
)


ggsave(
  
  file.path(
    OUTPUT_DIR,
    "Zonation_Volcano.png"
  ),
  
  p_volcano,
  
  width = 9,
  height = 7,
  
  dpi = 300
  
)



############################################################
# ==========================================================
#
# F3 STYLE
# CONTINUOUS PROTEIN GRADIENTS
#
# ==========================================================
############################################################



############################################################
# 49. SELECT TOP ZONATED PROTEINS
############################################################

top_gradient_proteins <- zonation_results %>%
  
  filter(
    
    !is.na(
      Q_value
    ),
    
    Q_value <
      Q_CUTOFF,
    
    abs(
      Beta
    ) >=
      BETA_CUTOFF
    
  ) %>%
  
  arrange(
    Q_value
  ) %>%
  
  slice_head(
    n = 20
  ) %>%
  
  pull(
    Protein.Group
  )



############################################################
# 50. GRADIENT DATA
############################################################

gradient_data <- dat %>%
  
  filter(
    
    Protein.Group %in%
      top_gradient_proteins
    
  )



############################################################
# 51. ORDER BY BETA
############################################################

protein_order <- zonation_results %>%
  
  filter(
    
    Protein.Group %in%
      top_gradient_proteins
    
  ) %>%
  
  arrange(
    Beta
  ) %>%
  
  pull(
    Protein
  )


gradient_data <- gradient_data %>%
  
  mutate(
    
    Protein =
      factor(
        
        Protein,
        
        levels =
          unique(
            protein_order
          )
        
      )
    
  )



############################################################
# 52. GRADIENT PLOT
############################################################

p_gradient <- ggplot(
  
  gradient_data,
  
  aes(
    
    x =
      Spatial_S,
    
    y =
      RobustScaled
    
  )
  
) +
  
  geom_point(
    
    size =
      0.6,
    
    alpha =
      0.25
    
  ) +
  
  geom_smooth(
    
    method =
      "loess",
    
    span =
      0.4,
    
    se =
      FALSE,
    
    linewidth =
      0.9
    
  ) +
  
  facet_wrap(
    
    ~ Protein,
    
    scales =
      "free_y",
    
    ncol =
      5
    
  ) +
  
  scale_x_continuous(
    
    limits =
      c(
        0,
        1
      ),
    
    breaks =
      c(
        0,
        0.5,
        1
      ),
    
    labels =
      c(
        "CV",
        "0.5",
        "PV"
      )
    
  ) +
  
  labs(
    
    title =
      "Top continuous zonated proteins",
    
    x =
      "CV → PV spatial position",
    
    y =
      "Robust-scaled protein abundance"
    
  ) +
  
  theme_classic(
    base_size = 10
  )



############################################################
# 53. SAVE GRADIENT PLOT
############################################################

ggsave(
  
  file.path(
    OUTPUT_DIR,
    "Top20_Continuous_Gradient.pdf"
  ),
  
  p_gradient,
  
  width = 13,
  height = 9
  
)


ggsave(
  
  file.path(
    OUTPUT_DIR,
    "Top20_Continuous_Gradient.png"
  ),
  
  p_gradient,
  
  width = 13,
  height = 9,
  
  dpi = 300
  
)



############################################################
# ==========================================================
#
# CANONICAL LIVER ZONATION VALIDATION
#
# ==========================================================
############################################################



############################################################
# 54. CANONICAL MARKERS
############################################################

marker_genes <- c(
  
  # -----------------------------------------
  # Pericentral / CV
  # -----------------------------------------
  
  "Glul",
  "Cyp2e1",
  "Cyp1a2",
  "Oat",
  
  # -----------------------------------------
  # Periportal / PV
  # -----------------------------------------
  
  "Cps1",
  "Ass1",
  "Pck1",
  "Hal",
  "Arg1"
  
)



############################################################
# 55. EXTRACT CANONICAL MARKER RESULTS
############################################################

marker_results <- zonation_results %>%
  
  filter(
    
    str_to_lower(
      Protein
    ) %in%
      
      str_to_lower(
        marker_genes
      )
    
  ) %>%
  
  arrange(
    Beta
  )


write_csv(
  
  marker_results,
  
  file.path(
    OUTPUT_DIR,
    "Known_Zonation_Marker_Results.csv"
  )
  
)


cat(
  "\nCanonical marker results:\n"
)


print(
  marker_results
)



############################################################
# 56. CANONICAL MARKER DATA
############################################################

marker_groups <- marker_results$Protein.Group


marker_data <- dat %>%
  
  filter(
    
    Protein.Group %in%
      marker_groups
    
  )



############################################################
# 57. CANONICAL MARKER GRADIENT PLOT
############################################################

p_marker <- ggplot(
  
  marker_data,
  
  aes(
    
    x =
      Spatial_S,
    
    y =
      RobustScaled,
    
    color =
      Group
    
  )
  
) +
  
  geom_point(
    
    size =
      0.6,
    
    alpha =
      0.20
    
  ) +
  
  geom_smooth(
    
    method =
      "loess",
    
    span =
      0.4,
    
    se =
      FALSE,
    
    linewidth =
      1
    
  ) +
  
  facet_wrap(
    
    ~ Protein,
    
    scales =
      "free_y",
    
    ncol =
      3
    
  ) +
  
  scale_x_continuous(
    
    limits =
      c(
        0,
        1
      ),
    
    breaks =
      c(
        0,
        0.5,
        1
      ),
    
    labels =
      c(
        "CV",
        "0.5",
        "PV"
      )
    
  ) +
  
  labs(
    
    title =
      "Canonical liver zonation markers",
    
    x =
      "CV → PV spatial position",
    
    y =
      "Robust-scaled protein abundance",
    
    color =
      "Group"
    
  ) +
  
  theme_classic(
    base_size = 11
  )



############################################################
# 58. SAVE CANONICAL MARKER PLOT
############################################################

ggsave(
  
  file.path(
    OUTPUT_DIR,
    "Canonical_Zonation_Markers_GF_SPF.pdf"
  ),
  
  p_marker,
  
  width = 11,
  height = 9
  
)


ggsave(
  
  file.path(
    OUTPUT_DIR,
    "Canonical_Zonation_Markers_GF_SPF.png"
  ),
  
  p_marker,
  
  width = 11,
  height = 9,
  dpi = 300
  
)



############################################################
# ==========================================================
#
# PART II
#
# GF vs SPF ZONATION COMPARISON
#
# ==========================================================
############################################################



############################################################
# 59. GROUP-SPECIFIC TOTAL ROI COUNTS
############################################################

group_roi_counts <- protein_long %>%
  
  distinct(
    Group,
    Run_Key
  ) %>%
  
  count(
    
    Group,
    
    name =
      "N_Total_ROI"
    
  )



############################################################
# 60. GROUP-SPECIFIC COMPLETENESS
############################################################

group_completeness <- protein_long %>%
  
  group_by(
    
    Group,
    Protein.Group,
    Protein
    
  ) %>%
  
  summarise(
    
    N_detected =
      sum(
        
        !is.na(
          Abundance
        ) &
          
          Abundance > 0
        
      ),
    
    .groups =
      "drop"
    
  ) %>%
  
  left_join(
    
    group_roi_counts,
    
    by =
      "Group"
    
  ) %>%
  
  mutate(
    
    Completeness =
      N_detected /
      N_Total_ROI
    
  )



############################################################
# 61. SAVE GROUP COMPLETENESS
############################################################

write_csv(
  
  group_completeness,
  
  file.path(
    OUTPUT_DIR,
    "GF_SPF_Protein_Completeness.csv"
  )
  
)



############################################################
# 62. REQUIRE ≥70% COMPLETENESS IN BOTH GROUPS
############################################################

proteins_both70 <- group_completeness %>%
  
  filter(
    
    Completeness >=
      COMPLETENESS_CUTOFF
    
  ) %>%
  
  distinct(
    
    Group,
    Protein.Group
    
  ) %>%
  
  count(
    
    Protein.Group,
    
    name =
      "N_Groups"
    
  ) %>%
  
  filter(
    
    N_Groups ==
      2
    
  ) %>%
  
  pull(
    Protein.Group
  )


cat(
  "\nProteins >=70% complete in BOTH GF and SPF:",
  length(
    proteins_both70
  ),
  "\n"
)



############################################################
# 63. GF/SPF ANALYSIS DATA
############################################################

dat_GFSPF <- dat %>%
  
  filter(
    
    Protein.Group %in%
      proteins_both70
    
  )



############################################################
# 64. ROBUST GF × SPF INTERACTION MODEL
#
# GF is reference group
#
# Beta_GF =
#   coefficient for Spatial_S
#
# Delta_Beta =
#   Spatial_S:GroupSPF
#
# Beta_SPF =
#   Beta_GF + Delta_Beta
#
############################################################

fit_interaction_lmm <- function(df) {
  
  tryCatch({
    
    ########################################################
    # Clean model data
    ########################################################
    
    df_model <- df %>%
      
      filter(
        
        !is.na(
          RobustScaled
        ),
        
        !is.na(
          Spatial_S
        ),
        
        !is.na(
          Group
        ),
        
        !is.na(
          Mouse
        )
        
      )
    
    
    ########################################################
    # Minimum conditions
    ########################################################
    
    if (
      
      nrow(
        df_model
      ) < 10 ||
      
      n_distinct(
        df_model$Spatial_S
      ) < 3 ||
      
      n_distinct(
        df_model$Group
      ) < 2 ||
      
      n_distinct(
        df_model$Mouse
      ) < 2
      
    ) {
      
      return(
        
        tibble(
          
          Beta_GF =
            NA_real_,
          
          Beta_SPF =
            NA_real_,
          
          Delta_Beta =
            NA_real_,
          
          SE =
            NA_real_,
          
          DF =
            NA_real_,
          
          Statistic =
            NA_real_,
          
          P_interaction =
            NA_real_,
          
          N =
            nrow(
              df_model
            ),
          
          N_Mouse =
            n_distinct(
              df_model$Mouse
            ),
          
          Singular =
            NA
          
        )
        
      )
      
    }
    
    
    ########################################################
    # Fit interaction LMM
    ########################################################
    
    model <- lmerTest::lmer(
      
      RobustScaled ~
        Spatial_S *
        Group +
        (1 | Mouse),
      
      data =
        df_model,
      
      REML =
        FALSE
      
    )
    
    
    ########################################################
    # Coefficient table
    ########################################################
    
    coef_table <-
      summary(
        model
      )$coefficients
    
    
    main_term <-
      "Spatial_S"
    
    
    interaction_term <-
      "Spatial_S:GroupSPF"
    
    
    ########################################################
    # Check required coefficients
    ########################################################
    
    if (
      
      !main_term %in%
      rownames(
        coef_table
      ) ||
      
      !interaction_term %in%
      rownames(
        coef_table
      )
      
    ) {
      
      stop(
        "Required interaction coefficient missing."
      )
      
    }
    
    
    ########################################################
    # P-value column
    ########################################################
    
    p_col <- grep(
      
      "^Pr",
      
      colnames(
        coef_table
      ),
      
      value =
        TRUE
      
    )
    
    
    ########################################################
    # DF column
    ########################################################
    
    df_col <- grep(
      
      "^df$",
      
      colnames(
        coef_table
      ),
      
      value =
        TRUE,
      
      ignore.case =
        TRUE
      
    )
    
    
    ########################################################
    # Beta GF
    ########################################################
    
    beta_GF <-
      coef_table[
        main_term,
        "Estimate"
      ]
    
    
    ########################################################
    # GF-SPF difference
    ########################################################
    
    delta_beta <-
      coef_table[
        interaction_term,
        "Estimate"
      ]
    
    
    ########################################################
    # SPF Beta
    ########################################################
    
    beta_SPF <-
      beta_GF +
      delta_beta
    
    
    ########################################################
    # Interaction SE
    ########################################################
    
    se <-
      coef_table[
        interaction_term,
        "Std. Error"
      ]
    
    
    ########################################################
    # Interaction statistic
    ########################################################
    
    statistic <-
      coef_table[
        interaction_term,
        "t value"
      ]
    
    
    ########################################################
    # Interaction P
    ########################################################
    
    if (
      length(
        p_col
      ) >= 1
    ) {
      
      pvalue <-
        coef_table[
          interaction_term,
          p_col[1]
        ]
      
    } else {
      
      pvalue <-
        NA_real_
      
    }
    
    
    ########################################################
    # DF
    ########################################################
    
    if (
      length(
        df_col
      ) >= 1
    ) {
      
      df_value <-
        coef_table[
          interaction_term,
          df_col[1]
        ]
      
    } else {
      
      df_value <-
        NA_real_
      
    }
    
    
    ########################################################
    # Singular
    ########################################################
    
    singular_flag <-
      lme4::isSingular(
        model,
        tol = 1e-4
      )
    
    
    ########################################################
    # Return
    ########################################################
    
    tibble(
      
      Beta_GF =
        as.numeric(
          beta_GF
        ),
      
      Beta_SPF =
        as.numeric(
          beta_SPF
        ),
      
      Delta_Beta =
        as.numeric(
          delta_beta
        ),
      
      SE =
        as.numeric(
          se
        ),
      
      DF =
        as.numeric(
          df_value
        ),
      
      Statistic =
        as.numeric(
          statistic
        ),
      
      P_interaction =
        as.numeric(
          pvalue
        ),
      
      N =
        nrow(
          df_model
        ),
      
      N_Mouse =
        n_distinct(
          df_model$Mouse
        ),
      
      Singular =
        singular_flag
      
    )
    
  },
  
  
  error = function(e) {
    
    tibble(
      
      Beta_GF =
        NA_real_,
      
      Beta_SPF =
        NA_real_,
      
      Delta_Beta =
        NA_real_,
      
      SE =
        NA_real_,
      
      DF =
        NA_real_,
      
      Statistic =
        NA_real_,
      
      P_interaction =
        NA_real_,
      
      N =
        nrow(
          df
        ),
      
      N_Mouse =
        n_distinct(
          df$Mouse
        ),
      
      Singular =
        NA
      
    )
    
  })
  
}



############################################################
# 65. RUN GF vs SPF INTERACTION MODELS
############################################################

message(
  "\nRunning GF vs SPF zonation interaction models..."
)


GFSPF_results <- dat_GFSPF %>%
  
  group_by(
    
    Protein.Group,
    Protein
    
  ) %>%
  
  nest() %>%
  
  mutate(
    
    Model_Result =
      map(
        
        data,
        
        fit_interaction_lmm
        
      )
    
  ) %>%
  
  select(
    -data
  ) %>%
  
  unnest(
    Model_Result
  ) %>%
  
  ungroup()



############################################################
# 66. BH CORRECTION
############################################################

GFSPF_results <- GFSPF_results %>%
  
  mutate(
    
    Q_interaction =
      p.adjust(
        
        P_interaction,
        
        method =
          "BH"
        
      ),
    
    Significant_Interaction =
      case_when(
        
        !is.na(
          Q_interaction
        ) &
          
          Q_interaction <
          0.05
        ~ "Significant",
        
        TRUE
        ~ "NS"
        
      )
    
  )



############################################################
# 67. SAVE GF-SPF RESULTS
############################################################

write_csv(
  
  GFSPF_results,
  
  file.path(
    OUTPUT_DIR,
    "GF_vs_SPF_Zonation_Interaction.csv"
  )
  
)



############################################################
# 68. GF BETA vs SPF BETA SCATTER
############################################################

p_GFSPF <- ggplot(
  
  GFSPF_results,
  
  aes(
    
    x =
      Beta_GF,
    
    y =
      Beta_SPF
    
  )
  
) +
  
  geom_hline(
    
    yintercept =
      0,
    
    linetype =
      "dashed"
    
  ) +
  
  geom_vline(
    
    xintercept =
      0,
    
    linetype =
      "dashed"
    
  ) +
  
  geom_abline(
    
    slope =
      1,
    
    intercept =
      0,
    
    linetype =
      "dotted"
    
  ) +
  
  geom_point(
    
    aes(
      
      shape =
        Significant_Interaction
      
    ),
    
    size =
      2,
    
    alpha =
      0.65,
    
    na.rm =
      TRUE
    
  ) +
  
  labs(
    
    title =
      "GF vs SPF liver zonation coefficients",
    
    subtitle =
      "Distance from diagonal reflects altered spatial gradient",
    
    x =
      "GF zonation coefficient (Beta)",
    
    y =
      "SPF zonation coefficient (Beta)",
    
    shape =
      "Interaction"
    
  ) +
  
  theme_classic(
    base_size = 13
  )



############################################################
# 69. SAVE GF/SPF BETA PLOT
############################################################

ggsave(
  
  file.path(
    OUTPUT_DIR,
    "GF_vs_SPF_Beta_Scatter.pdf"
  ),
  
  p_GFSPF,
  
  width = 8,
  height = 7
  
)



############################################################
# 70. TOP 30 ALTERED ZONATION PROTEINS
############################################################

top_changed <- GFSPF_results %>%
  
  filter(
    
    !is.na(
      Q_interaction
    )
    
  ) %>%
  
  arrange(
    Q_interaction
  ) %>%
  
  slice_head(
    n = 30
  )



############################################################
# 71. DELTA BETA PLOT
############################################################

p_delta <- ggplot(
  
  top_changed,
  
  aes(
    
    x =
      reorder(
        Protein,
        Delta_Beta
      ),
    
    y =
      Delta_Beta
    
  )
  
) +
  
  geom_col() +
  
  coord_flip() +
  
  geom_hline(
    
    yintercept =
      0,
    
    linetype =
      "dashed"
    
  ) +
  
  labs(
    
    title =
      "Top GF-SPF changes in spatial zonation",
    
    subtitle =
      expression(
        Delta * beta ==
          beta[SPF] -
          beta[GF]
      ),
    
    x =
      NULL,
    
    y =
      expression(
        Delta * beta
      )
    
  ) +
  
  theme_classic(
    base_size = 12
  )



############################################################
# 72. SAVE DELTA BETA PLOT
############################################################

ggsave(
  
  file.path(
    OUTPUT_DIR,
    "GF_vs_SPF_Top30_DeltaBeta.pdf"
  ),
  
  p_delta,
  
  width = 8,
  height = 9
  
)



############################################################
# 73. SIGNIFICANT GF-SPF INTERACTIONS
############################################################

significant_GFSPF <- GFSPF_results %>%
  
  filter(
    
    !is.na(
      Q_interaction
    ),
    
    Q_interaction <
      0.05
    
  ) %>%
  
  arrange(
    Q_interaction
  )


write_csv(
  
  significant_GFSPF,
  
  file.path(
    OUTPUT_DIR,
    "GF_vs_SPF_Significant_Zonation_Changes.csv"
  )
  
)



############################################################
# ==========================================================
#
# PART III
#
# GF AND SPF SPATIAL_S HEATMAPS
#
# X-axis: Spatial_S (CV = 0 -> PV = 1)
# Y-axis: Protein
# Fill: median robust-scaled protein abundance
#
# Each mouse is summarized first within each Spatial_S bin.
# Group heatmap values are then calculated as the median of
# the mouse-level medians, so mice with more ROIs do not
# receive disproportionate weight.
#
# ==========================================================
############################################################



############################################################
# 74. HEATMAP SETTINGS
############################################################

# Protein selection mode:
#   "zonated"   = all proteins passing Q and Beta cutoffs
#   "top"       = top HEATMAP_TOP_N zonated proteins by Q
#   "canonical" = marker_genes defined in Section 54
#   "all"       = all proteins retained for the analysis

HEATMAP_PROTEIN_SET <- as.character(heatmap_protein_set)

HEATMAP_TOP_N <- as.integer(heatmap_top_n)

# Number of equal-width Spatial_S bins between 0 and 1
HEATMAP_N_BINS <- as.integer(heatmap_n_bins)

# A heatmap cell is shown when at least this many mice
# contribute to the corresponding group/protein/bin.
HEATMAP_MIN_MICE_PER_BIN <- as.integer(heatmap_min_mice_per_bin)

# When all proteins are selected, readable protein labels are
# provided in multi-page PDFs.
HEATMAP_CREATE_MULTIPAGE <- isTRUE(heatmap_create_multipage)

HEATMAP_PROTEINS_PER_PAGE <- as.integer(heatmap_proteins_per_page)



############################################################
# 75. SELECT PROTEINS FOR THE HEATMAP
############################################################

if (
  HEATMAP_PROTEIN_SET == "zonated"
) {
  
  heatmap_selection <- zonation_results %>%
    filter(
      !is.na(Q_value),
      Q_value < Q_CUTOFF,
      abs(Beta) >= BETA_CUTOFF
    )
  
} else if (
  HEATMAP_PROTEIN_SET == "top"
) {
  
  heatmap_selection <- zonation_results %>%
    filter(
      !is.na(Q_value),
      Q_value < Q_CUTOFF,
      abs(Beta) >= BETA_CUTOFF
    ) %>%
    arrange(Q_value, desc(abs(Beta))) %>%
    slice_head(n = HEATMAP_TOP_N)
  
} else if (
  HEATMAP_PROTEIN_SET == "canonical"
) {
  
  heatmap_selection <- zonation_results %>%
    filter(
      str_to_lower(Protein) %in%
        str_to_lower(marker_genes)
    )
  
} else if (
  HEATMAP_PROTEIN_SET == "all"
) {
  
  # Use every cRAP-free protein retained in dat, regardless
  # of LMM significance or whether Beta was estimated.
  heatmap_selection <- dat %>%
    distinct(
      Protein.Group,
      Protein
    ) %>%
    left_join(
      zonation_results %>%
        select(
          Protein.Group,
          Protein,
          Beta,
          Q_value,
          Zonation
        ),
      by = c(
        "Protein.Group",
        "Protein"
      )
    )
  
} else {
  
  stop(
    paste0(
      "\nUnknown HEATMAP_PROTEIN_SET: ",
      HEATMAP_PROTEIN_SET,
      "\nUse zonated, top, canonical, or all."
    )
  )
  
}


heatmap_selection <- heatmap_selection %>%
  distinct(
    Protein.Group,
    .keep_all = TRUE
  ) %>%
  mutate(
    Protein_Base = case_when(
      is.na(Protein) | Protein == "" ~ Protein.Group,
      TRUE ~ Protein
    )
  ) %>%
  group_by(Protein_Base) %>%
  mutate(
    Protein_Label = if (
      n() > 1
    ) {
      paste0(
        Protein_Base,
        " [",
        Protein.Group,
        "]"
      )
    } else {
      Protein_Base
    }
  ) %>%
  ungroup()


if (
  nrow(heatmap_selection) == 0
) {
  
  stop(
    paste0(
      "\nNo proteins were selected for the heatmap using ",
      "HEATMAP_PROTEIN_SET = '",
      HEATMAP_PROTEIN_SET,
      "'."
    )
  )
  
}


cat(
  "\nProteins selected for GF/SPF heatmaps:",
  nrow(heatmap_selection),
  "\n"
)


cat(
  "cRAP-free proteins retained in dat:",
  n_distinct(dat$Protein.Group),
  "\n"
)


if (
  HEATMAP_PROTEIN_SET == "all" &&
  n_distinct(heatmap_selection$Protein.Group) !=
    n_distinct(dat$Protein.Group)
) {
  
  stop(
    "\nThe all-protein heatmap selection does not match dat."
  )
  
}


all_retained_heatmap_ids <- dat %>%
  distinct(Protein.Group) %>%
  pull(Protein.Group)


zonated_heatmap_ids <- zonation_results %>%
  filter(
    !is.na(Q_value),
    Q_value < Q_CUTOFF,
    !is.na(Beta),
    abs(Beta) >= BETA_CUTOFF
  ) %>%
  distinct(Protein.Group) %>%
  pull(Protein.Group)


actual_heatmap_ids <- heatmap_selection %>%
  distinct(Protein.Group) %>%
  pull(Protein.Group)


heatmap_selection_check <- tibble(
  HEATMAP_PROTEIN_SET = HEATMAP_PROTEIN_SET,
  All_Retained_Protein_N = length(
    all_retained_heatmap_ids
  ),
  Zonated_Protein_N = length(
    zonated_heatmap_ids
  ),
  Actual_Heatmap_Protein_N = length(
    actual_heatmap_ids
  ),
  Actual_Equals_All_Retained = setequal(
    actual_heatmap_ids,
    all_retained_heatmap_ids
  ),
  Actual_Equals_Zonated = setequal(
    actual_heatmap_ids,
    zonated_heatmap_ids
  )
)


print(heatmap_selection_check)


write_csv(
  heatmap_selection_check,
  file.path(
    OUTPUT_DIR,
    "Heatmap_Protein_Selection_Check.csv"
  )
)



############################################################
# 76. ASSIGN Spatial_S BINS
############################################################

heatmap_mouse_bins <- dat %>%
  semi_join(
    heatmap_selection,
    by = "Protein.Group"
  ) %>%
  filter(
    Group %in% c("GF", "SPF"),
    !is.na(Spatial_S),
    Spatial_S >= 0,
    Spatial_S <= 1,
    is.finite(RobustScaled)
  ) %>%
  mutate(
    Spatial_Bin = pmin(
      HEATMAP_N_BINS,
      floor(Spatial_S * HEATMAP_N_BINS) + 1L
    ),
    Spatial_S_Mid =
      (Spatial_Bin - 0.5) / HEATMAP_N_BINS
  ) %>%
  group_by(
    Group,
    Mouse,
    Protein.Group,
    Spatial_Bin,
    Spatial_S_Mid
  ) %>%
  summarise(
    Mouse_Median = median(
      RobustScaled,
      na.rm = TRUE
    ),
    N_ROI_Mouse = n(),
    .groups = "drop"
  )



############################################################
# 77. CALCULATE GROUP-LEVEL HEATMAP VALUES
############################################################

heatmap_summary <- heatmap_mouse_bins %>%
  group_by(
    Group,
    Protein.Group,
    Spatial_Bin,
    Spatial_S_Mid
  ) %>%
  summarise(
    Heatmap_Value = median(
      Mouse_Median,
      na.rm = TRUE
    ),
    N_Mouse = n_distinct(Mouse),
    N_ROI = sum(N_ROI_Mouse),
    .groups = "drop"
  ) %>%
  mutate(
    Heatmap_Value = if_else(
      N_Mouse >= HEATMAP_MIN_MICE_PER_BIN,
      Heatmap_Value,
      NA_real_
    )
  )


# Create the complete GF/SPF x protein x Spatial_S grid so
# bins without observations appear explicitly as missing.

heatmap_grid <- tidyr::expand_grid(
  Group = c("GF", "SPF"),
  Protein.Group = heatmap_selection$Protein.Group,
  Spatial_Bin = seq_len(HEATMAP_N_BINS)
) %>%
  mutate(
    Spatial_S_Mid =
      (Spatial_Bin - 0.5) / HEATMAP_N_BINS
  ) %>%
  left_join(
    heatmap_summary,
    by = c(
      "Group",
      "Protein.Group",
      "Spatial_Bin",
      "Spatial_S_Mid"
    )
  ) %>%
  left_join(
    heatmap_selection %>%
      select(
        Protein.Group,
        Protein_Label,
        Beta,
        Q_value,
        Zonation
      ),
    by = "Protein.Group"
  )



############################################################
# 78. SET A COMMON PROTEIN ORDER AND COLOR SCALE
############################################################

# Negative Beta (CV-enriched) proteins are placed at the top,
# followed by positive Beta (PV-enriched) proteins.

heatmap_protein_order <- heatmap_selection %>%
  arrange(
    is.na(Beta),
    Beta,
    Q_value
  ) %>%
  pull(Protein_Label) %>%
  unique()


heatmap_grid <- heatmap_grid %>%
  mutate(
    Group = factor(
      Group,
      levels = c("GF", "SPF")
    ),
    Protein_Label = factor(
      Protein_Label,
      levels = rev(heatmap_protein_order)
    )
  )


# The 98th percentile prevents a few extreme cells from
# compressing the color contrast. The same limits are used
# for GF and SPF, allowing direct visual comparison.

heatmap_fill_limit <- quantile(
  abs(heatmap_grid$Heatmap_Value),
  probs = 0.98,
  na.rm = TRUE
)


if (
  !is.finite(heatmap_fill_limit) ||
  heatmap_fill_limit <= 0
) {
  heatmap_fill_limit <- 1
}



############################################################
# 79. HEATMAP PLOTTING FUNCTION
############################################################

# A single overview containing thousands of proteins cannot
# display readable Y-axis labels. Labels are shown here only
# for small selections; the multi-page PDFs created below
# always show every protein label.

show_overview_protein_labels <-
  nrow(heatmap_selection) <= 150

plot_spatial_heatmap <- function(group_name) {
  
  ggplot(
    heatmap_grid %>%
      filter(Group == group_name),
    aes(
      x = Spatial_S_Mid,
      y = Protein_Label,
      fill = Heatmap_Value
    )
  ) +
    geom_tile(
      width = 1 / HEATMAP_N_BINS,
      height = 0.95
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      labels = c("CV", "0.5", "PV"),
      expand = c(0, 0)
    ) +
    scale_y_discrete(
      drop = FALSE
    ) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(
        -heatmap_fill_limit,
        heatmap_fill_limit
      ),
      oob = scales::squish,
      na.value = "grey90"
    ) +
    labs(
      title = paste0(
        group_name,
        " CV-PV protein zonation"
      ),
      subtitle = paste0(
        "Spatial_S bins = ",
        HEATMAP_N_BINS,
        " | protein set = ",
        HEATMAP_PROTEIN_SET
      ),
      x = "Spatial_S (CV = 0 -> PV = 1)",
      y = NULL,
      fill = "Median\nrobust-scaled\nabundance"
    ) +
    theme_minimal(
      base_size = 10
    ) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = if (
        show_overview_protein_labels
      ) {
        element_text(size = 6)
      } else {
        element_blank()
      },
      axis.ticks.y = if (
        show_overview_protein_labels
      ) {
        element_line()
      } else {
        element_blank()
      },
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
}



############################################################
# 80. CREATE GF AND SPF HEATMAPS
############################################################

p_heatmap_GF <- plot_spatial_heatmap("GF")

p_heatmap_SPF <- plot_spatial_heatmap("SPF")

p_heatmap_combined <-
  p_heatmap_GF +
  p_heatmap_SPF +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "right"
  )



############################################################
# 81. SAVE HEATMAP DATA AND FIGURES
############################################################

heatmap_file_tag <- paste0(
  HEATMAP_PROTEIN_SET,
  "_",
  n_distinct(heatmap_selection$Protein.Group),
  "proteins"
)


write_csv(
  heatmap_grid %>%
    mutate(
      Protein_Label = as.character(Protein_Label),
      Group = as.character(Group)
    ) %>%
    arrange(
      Group,
      Beta,
      Spatial_Bin
    ),
  file.path(
    OUTPUT_DIR,
    paste0(
      "GF_SPF_Spatial_S_Heatmap_Data_",
      heatmap_file_tag,
      ".csv"
    )
  )
)


# Increase the overview height with the number of proteins,
# while keeping it within a practical PDF range. When more
# than 150 proteins are selected, labels are supplied by the
# multi-page PDFs instead of the overview.

heatmap_height <- max(
  6,
  min(
    40,
    2.5 + 0.18 * nrow(heatmap_selection)
  )
)


ggsave(
  file.path(
    OUTPUT_DIR,
    paste0(
      "GF_Spatial_S_Heatmap_Overview_",
      heatmap_file_tag,
      ".pdf"
    )
  ),
  p_heatmap_GF,
  width = 10,
  height = heatmap_height,
  limitsize = FALSE
)


ggsave(
  file.path(
    OUTPUT_DIR,
    paste0(
      "GF_Spatial_S_Heatmap_Overview_",
      heatmap_file_tag,
      ".png"
    )
  ),
  p_heatmap_GF,
  width = 10,
  height = heatmap_height,
  dpi = 300,
  limitsize = FALSE
)


ggsave(
  file.path(
    OUTPUT_DIR,
    paste0(
      "SPF_Spatial_S_Heatmap_Overview_",
      heatmap_file_tag,
      ".pdf"
    )
  ),
  p_heatmap_SPF,
  width = 10,
  height = heatmap_height,
  limitsize = FALSE
)


ggsave(
  file.path(
    OUTPUT_DIR,
    paste0(
      "SPF_Spatial_S_Heatmap_Overview_",
      heatmap_file_tag,
      ".png"
    )
  ),
  p_heatmap_SPF,
  width = 10,
  height = heatmap_height,
  dpi = 300,
  limitsize = FALSE
)


ggsave(
  file.path(
    OUTPUT_DIR,
    paste0(
      "GF_SPF_Spatial_S_Heatmap_Combined_",
      heatmap_file_tag,
      ".pdf"
    )
  ),
  p_heatmap_combined,
  width = 18,
  height = heatmap_height,
  limitsize = FALSE
)



############################################################
# 82. MULTI-PAGE HEATMAPS WITH ALL PROTEIN LABELS
############################################################

if (HEATMAP_CREATE_MULTIPAGE) {
  
  protein_page_table <- tibble(
    Protein_Label = heatmap_protein_order
  ) %>%
    mutate(
      Heatmap_Page = ceiling(
        row_number() /
          HEATMAP_PROTEINS_PER_PAGE
      )
    )
  
  
  heatmap_page_data <- heatmap_grid %>%
    mutate(
      Protein_Label = as.character(Protein_Label),
      Group = as.character(Group)
    ) %>%
    left_join(
      protein_page_table,
      by = "Protein_Label"
    )
  
  
  n_heatmap_pages <- max(
    protein_page_table$Heatmap_Page
  )
  
  
  plot_spatial_heatmap_page <- function(
      group_name,
      page_number
  ) {
    
    page_proteins <- protein_page_table %>%
      filter(
        Heatmap_Page == page_number
      ) %>%
      pull(Protein_Label)
    
    
    page_plot_data <- heatmap_page_data %>%
      filter(
        Group == group_name,
        Heatmap_Page == page_number
      ) %>%
      mutate(
        Protein_Label = factor(
          Protein_Label,
          levels = rev(page_proteins)
        )
      )
    
    
    ggplot(
      page_plot_data,
      aes(
        x = Spatial_S_Mid,
        y = Protein_Label,
        fill = Heatmap_Value
      )
    ) +
      geom_tile(
        width = 1 / HEATMAP_N_BINS,
        height = 0.95
      ) +
      scale_x_continuous(
        limits = c(0, 1),
        breaks = c(0, 0.5, 1),
        labels = c("CV", "0.5", "PV"),
        expand = c(0, 0)
      ) +
      scale_y_discrete(
        drop = TRUE
      ) +
      scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        limits = c(
          -heatmap_fill_limit,
          heatmap_fill_limit
        ),
        oob = scales::squish,
        na.value = "grey90"
      ) +
      labs(
        title = paste0(
          group_name,
          " CV-PV protein zonation"
        ),
        subtitle = paste0(
          "Protein set = ",
          HEATMAP_PROTEIN_SET,
          " | Page ",
          page_number,
          " of ",
          n_heatmap_pages
        ),
        x = "Spatial_S (CV = 0 -> PV = 1)",
        y = NULL,
        fill = "Median\nrobust-scaled\nabundance"
      ) +
      theme_minimal(
        base_size = 10
      ) +
      theme(
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 6),
        plot.title = element_text(face = "bold"),
        legend.position = "right"
      )
    
  }
  
  
  write_csv(
    protein_page_table,
    file.path(
      OUTPUT_DIR,
      paste0(
        "Spatial_S_Heatmap_Protein_Page_Index_",
        heatmap_file_tag,
        ".csv"
      )
    )
  )
  
  
  for (
    group_name in c("GF", "SPF")
  ) {
    
    pdf_file <- file.path(
      OUTPUT_DIR,
      paste0(
        group_name,
        "_Spatial_S_Heatmap_Multipage_",
        heatmap_file_tag,
        ".pdf"
      )
    )
    
    grDevices::pdf(
      pdf_file,
      width = 10,
      height = 20,
      onefile = TRUE
    )
    
    for (
      page_number in seq_len(n_heatmap_pages)
    ) {
      
      print(
        plot_spatial_heatmap_page(
          group_name = group_name,
          page_number = page_number
        )
      )
      
    }
    
    grDevices::dev.off()
    
  }
  
}



############################################################
# ==========================================================
#
# GF / SPF INDEPENDENT ZONATION ANALYSIS
#
# The Total analysis above uses all retained GF and SPF ROIs.
# This additional branch starts again from matched abundance data and applies
# completeness, robust scaling and LMM independently within each group.
# Mouse is defined as Set in this branch, matching Separation.R.
#
# ==========================================================
############################################################

SEPARATE_OUTPUT_DIR <- file.path(OUTPUT_DIR, "GF_SPF_Separate")
dir.create(SEPARATE_OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

write_csv(
  retained_roi_list,
  file.path(SEPARATE_OUTPUT_DIR, "QC_Retained_ROI_List_All.csv")
)

run_group_zonation <- function(group_name, protein_long_all, retained_roi_all) {
  message("\n============================================================")
  message("STARTING INDEPENDENT ZONATION ANALYSIS: ", group_name)
  message("============================================================")

  group_dir <- file.path(SEPARATE_OUTPUT_DIR, group_name)
  dir.create(group_dir, showWarnings = FALSE, recursive = TRUE)

  # Separation.R defines each Set as one biological mouse. Keep the original
  # Mouse value for audit while using Set for the random intercept.
  protein_group <- protein_long_all %>%
    filter(as.character(Group) == group_name) %>%
    mutate(
      Mouse_original = as.character(Mouse),
      Group = factor(as.character(Group), levels = group_name),
      Set = droplevels(factor(Set)),
      Mouse = droplevels(factor(Set))
    )

  retained_roi_group <- retained_roi_all %>%
    filter(as.character(Group) == group_name) %>%
    mutate(
      Mouse_original = as.character(Mouse),
      Group = factor(as.character(Group), levels = group_name),
      Set = droplevels(factor(Set)),
      Mouse = droplevels(factor(Set))
    )

  if (!nrow(protein_group)) {
    stop(paste0("No protein data found for group: ", group_name), call. = FALSE)
  }

  mouse_mapping <- protein_group %>%
    distinct(Group, Mouse_original, Set, Mouse) %>%
    arrange(Set)

  if (n_distinct(protein_group$Mouse) < 2L) {
    stop(paste0(group_name, " has fewer than two Set/Mouse levels."), call. = FALSE)
  }

  write_csv(
    mouse_mapping,
    file.path(group_dir, paste0(group_name, "_Set_Mouse_Mapping.csv"))
  )
  write_csv(
    retained_roi_group,
    file.path(group_dir, paste0(group_name, "_QC_Retained_ROI_List.csv"))
  )

  roi_summary_group <- retained_roi_group %>% count(Set, name = "N_ROI")
  write_csv(
    roi_summary_group,
    file.path(group_dir, paste0(group_name, "_ROI_Count_Per_Set.csv"))
  )

  n_roi_group <- n_distinct(protein_group$Run_Key)
  completeness_table_group <- protein_group %>%
    group_by(Protein.Group, Protein) %>%
    summarise(
      N_detected = sum(!is.na(Abundance) & Abundance > 0),
      .groups = "drop"
    ) %>%
    mutate(
      Total_ROI = n_roi_group,
      Completeness = N_detected / Total_ROI
    )

  retained_proteins_group <- completeness_table_group %>%
    filter(Completeness >= COMPLETENESS_CUTOFF)

  if (!nrow(retained_proteins_group)) {
    stop(paste0(group_name, ": no proteins passed the completeness cutoff."),
         call. = FALSE)
  }

  write_csv(
    completeness_table_group,
    file.path(group_dir, paste0(group_name, "_Protein_Completeness_All.csv"))
  )
  write_csv(
    retained_proteins_group,
    file.path(group_dir, paste0(group_name, "_Protein_Completeness_Retained.csv"))
  )

  # Group-specific completeness and group-specific robust scaling.
  dat_group <- protein_group %>%
    semi_join(retained_proteins_group, by = c("Protein.Group", "Protein")) %>%
    mutate(
      Log2Abundance = if_else(
        !is.na(Abundance) & Abundance > 0,
        log2(Abundance),
        NA_real_
      )
    ) %>%
    group_by(Protein.Group, Protein) %>%
    mutate(
      Protein_Median = median(Log2Abundance, na.rm = TRUE),
      Protein_IQR = IQR(Log2Abundance, na.rm = TRUE, type = 7),
      RobustScaled = if_else(
        is.finite(Protein_IQR) & Protein_IQR > 0,
        (Log2Abundance - Protein_Median) / Protein_IQR,
        NA_real_
      )
    ) %>%
    ungroup() %>%
    filter(
      is.finite(RobustScaled),
      !is.na(Spatial_S),
      is.finite(Spatial_S),
      !is.na(Mouse)
    )

  write_csv(
    dat_group,
    file.path(group_dir, paste0(group_name, "_Spatial_Proteomics_AnalysisReady.csv"))
  )

  p_spatial_group <- retained_roi_group %>%
    ggplot(aes(x = Spatial_S)) +
    geom_histogram(bins = 20) +
    facet_wrap(~ Image, scales = "free_y") +
    scale_x_continuous(limits = c(0, 1)) +
    labs(
      title = paste0(group_name, " distribution of CV-PV spatial coordinates"),
      x = "Spatial_S (CV = 0 → PV = 1)",
      y = "Number of QC-retained ROIs"
    ) +
    theme_classic()

  ggsave(
    file.path(group_dir, paste0(group_name, "_Spatial_S_Distribution.pdf")),
    p_spatial_group, width = 12, height = 8
  )
  ggsave(
    file.path(group_dir, paste0(group_name, "_Spatial_S_Distribution.png")),
    p_spatial_group, width = 12, height = 8, dpi = 300
  )

  message("Running ", group_name, " group-specific zonation LMMs...")
  zonation_results_group <- dat_group %>%
    group_by(Protein.Group, Protein) %>%
    nest() %>%
    mutate(Model_Result = map(data, fit_zonation_lmm)) %>%
    select(-data) %>%
    unnest(Model_Result) %>%
    ungroup() %>%
    mutate(
      Q_value = p.adjust(P_value, method = "BH"),
      MinusLog10Q = if_else(
        !is.na(Q_value),
        -log10(pmax(Q_value, 1e-300)),
        NA_real_
      ),
      Zonation = case_when(
        !is.na(Q_value) & Q_value < Q_CUTOFF & Beta <= -BETA_CUTOFF ~ "Central",
        !is.na(Q_value) & Q_value < Q_CUTOFF & Beta >= BETA_CUTOFF ~ "Portal",
        TRUE ~ "NS"
      )
    )

  write_csv(
    zonation_results_group,
    file.path(group_dir, paste0(group_name, "_Zonation_Coefficients.csv"))
  )

  zonation_summary_group <- zonation_results_group %>%
    count(Zonation, name = "N_Proteins") %>%
    mutate(Group = group_name, .before = 1)
  write_csv(
    zonation_summary_group,
    file.path(group_dir, paste0(group_name, "_Zonation_Summary.csv"))
  )

  # Independently fitted group-specific volcano.
  top_labels_group <- bind_rows(
    zonation_results_group %>%
      filter(Zonation == "Central") %>%
      arrange(Q_value, Beta) %>%
      slice_head(n = TOP_LABEL_N),
    zonation_results_group %>%
      filter(Zonation == "Portal") %>%
      arrange(Q_value, desc(Beta)) %>%
      slice_head(n = TOP_LABEL_N)
  )

  p_volcano_group <- ggplot(
    zonation_results_group,
    aes(x = Beta, y = MinusLog10Q, color = Zonation)
  ) +
    geom_point(size = 1.7, alpha = 0.65, na.rm = TRUE) +
    geom_vline(xintercept = c(-BETA_CUTOFF, BETA_CUTOFF), linetype = "dashed") +
    geom_hline(yintercept = -log10(Q_CUTOFF), linetype = "dashed") +
    geom_text_repel(
      data = top_labels_group,
      aes(label = Protein),
      size = 3,
      max.overlaps = Inf
    ) +
    labs(
      title = paste0(group_name, " continuous liver protein zonation"),
      subtitle = paste0(
        group_name,
        " only | group-specific completeness and RobustScaler"
      ),
      x = "Zonation coefficient (Beta)\nCentral ←          → Portal",
      y = expression(-log[10](Q)),
      color = "Zonation"
    ) +
    theme_classic(base_size = 13) +
    theme(legend.position = "top")

  ggsave(
    file.path(group_dir, paste0(group_name, "_Zonation_Volcano.pdf")),
    p_volcano_group, width = 9, height = 7
  )
  ggsave(
    file.path(group_dir, paste0(group_name, "_Zonation_Volcano.png")),
    p_volcano_group, width = 9, height = 7, dpi = 300
  )

  # Group-specific Top20 gradients.
  top_gradient_ids <- zonation_results_group %>%
    filter(
      !is.na(Q_value),
      Q_value < Q_CUTOFF,
      abs(Beta) >= BETA_CUTOFF
    ) %>%
    arrange(Q_value) %>%
    slice_head(n = 20) %>%
    pull(Protein.Group)

  p_gradient_group <- NULL
  if (length(top_gradient_ids)) {
    protein_order_group <- zonation_results_group %>%
      filter(Protein.Group %in% top_gradient_ids) %>%
      arrange(Beta) %>%
      pull(Protein)

    gradient_data_group <- dat_group %>%
      filter(Protein.Group %in% top_gradient_ids) %>%
      mutate(Protein = factor(Protein, levels = unique(protein_order_group)))

    p_gradient_group <- ggplot(
      gradient_data_group,
      aes(x = Spatial_S, y = RobustScaled)
    ) +
      geom_point(size = 0.6, alpha = 0.25) +
      geom_smooth(
        method = "loess", span = 0.4, se = FALSE, linewidth = 0.9
      ) +
      facet_wrap(~ Protein, scales = "free_y", ncol = 5) +
      scale_x_continuous(
        limits = c(0, 1),
        breaks = c(0, 0.5, 1),
        labels = c("CV", "0.5", "PV")
      ) +
      labs(
        title = paste0(group_name, " top continuous zonated proteins"),
        subtitle = "Group-specific Top20",
        x = "CV → PV spatial position",
        y = "Group-specific robust-scaled protein abundance"
      ) +
      theme_classic(base_size = 10)

    ggsave(
      file.path(group_dir, paste0(group_name, "_Top20_Continuous_Gradient.pdf")),
      p_gradient_group, width = 13, height = 9
    )
    ggsave(
      file.path(group_dir, paste0(group_name, "_Top20_Continuous_Gradient.png")),
      p_gradient_group, width = 13, height = 9, dpi = 300
    )
  } else {
    warning(paste0(group_name, ": no proteins passed the gradient thresholds."))
  }

  marker_results_group <- zonation_results_group %>%
    filter(str_to_lower(Protein) %in% str_to_lower(marker_genes)) %>%
    arrange(Beta)
  write_csv(
    marker_results_group,
    file.path(group_dir, paste0(group_name, "_Known_Zonation_Marker_Results.csv"))
  )

  p_marker_group <- NULL
  if (nrow(marker_results_group)) {
    marker_order_group <- unique(marker_results_group$Protein)
    marker_data_group <- dat_group %>%
      filter(Protein.Group %in% marker_results_group$Protein.Group) %>%
      mutate(Protein = factor(Protein, levels = marker_order_group))

    p_marker_group <- ggplot(
      marker_data_group,
      aes(x = Spatial_S, y = RobustScaled)
    ) +
      geom_point(size = 0.6, alpha = 0.20) +
      geom_smooth(method = "loess", span = 0.4, se = FALSE, linewidth = 1) +
      facet_wrap(~ Protein, scales = "free_y", ncol = 3) +
      scale_x_continuous(
        limits = c(0, 1),
        breaks = c(0, 0.5, 1),
        labels = c("CV", "0.5", "PV")
      ) +
      labs(
        title = paste0(group_name, " canonical liver zonation markers"),
        x = "CV → PV spatial position",
        y = "Group-specific robust-scaled protein abundance"
      ) +
      theme_classic(base_size = 11)

    ggsave(
      file.path(group_dir, paste0(group_name, "_Canonical_Zonation_Markers.pdf")),
      p_marker_group, width = 11, height = 9
    )
    ggsave(
      file.path(group_dir, paste0(group_name, "_Canonical_Zonation_Markers.png")),
      p_marker_group, width = 11, height = 9, dpi = 300
    )
  }

  group_summary <- tibble(
    Group = group_name,
    Total_QC_Retained_ROI = n_roi_group,
    N_Mouse = n_distinct(dat_group$Mouse),
    Proteins_Before_Completeness = nrow(completeness_table_group),
    Proteins_Retained = nrow(retained_proteins_group),
    Successful_Zonation_LMM = sum(!is.na(zonation_results_group$P_value)),
    Failed_Zonation_LMM = sum(is.na(zonation_results_group$P_value)),
    Singular_Zonation_LMM = sum(zonation_results_group$Singular == TRUE, na.rm = TRUE),
    Strong_Central = sum(zonation_results_group$Zonation == "Central", na.rm = TRUE),
    Strong_Portal = sum(zonation_results_group$Zonation == "Portal", na.rm = TRUE)
  )

  write_csv(
    group_summary,
    file.path(group_dir, paste0(group_name, "_Analysis_Summary.csv"))
  )

  list(
    Group = group_name,
    data = dat_group,
    zonation_results = zonation_results_group,
    summary = group_summary,
    p_volcano = p_volcano_group,
    p_gradient = p_gradient_group,
    p_marker = p_marker_group
  )
}

GF_result <- run_group_zonation("GF", protein_long, retained_roi_list)
SPF_result <- run_group_zonation("SPF", protein_long, retained_roi_list)

GF_SPF_separate_summary <- bind_rows(GF_result$summary, SPF_result$summary)
write_csv(
  GF_SPF_separate_summary,
  file.path(SEPARATE_OUTPUT_DIR, "GF_SPF_Separate_Analysis_Summary.csv")
)

p_separate_volcano <-
  GF_result$p_volcano +
  SPF_result$p_volcano +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

ggsave(
  file.path(SEPARATE_OUTPUT_DIR, "GF_SPF_Separate_Volcano_SideBySide.pdf"),
  p_separate_volcano, width = 18, height = 7
)
ggsave(
  file.path(SEPARATE_OUTPUT_DIR, "GF_SPF_Separate_Volcano_SideBySide.png"),
  p_separate_volcano, width = 18, height = 7, dpi = 300
)

if (!is.null(GF_result$p_gradient) && !is.null(SPF_result$p_gradient)) {
  p_separate_gradient <-
    GF_result$p_gradient +
    SPF_result$p_gradient +
    plot_layout(ncol = 1)

  ggsave(
    file.path(SEPARATE_OUTPUT_DIR, "GF_SPF_Separate_Gradient_Combined.pdf"),
    p_separate_gradient, width = 13, height = 18, limitsize = FALSE
  )
  ggsave(
    file.path(SEPARATE_OUTPUT_DIR, "GF_SPF_Separate_Gradient_Combined.png"),
    p_separate_gradient, width = 13, height = 18, dpi = 300, limitsize = FALSE
  )
}

GF_SPF_beta_descriptive <- GF_result$zonation_results %>%
  select(
    Protein.Group, Protein,
    Beta_GF = Beta, SE_GF = SE, Q_GF = Q_value, Zonation_GF = Zonation
  ) %>%
  full_join(
    SPF_result$zonation_results %>%
      select(
        Protein.Group, Protein,
        Beta_SPF = Beta, SE_SPF = SE, Q_SPF = Q_value, Zonation_SPF = Zonation
      ),
    by = c("Protein.Group", "Protein")
  ) %>%
  mutate(Delta_Beta_Descriptive = Beta_SPF - Beta_GF)

write_csv(
  GF_SPF_beta_descriptive,
  file.path(SEPARATE_OUTPUT_DIR, "GF_SPF_Independent_Beta_Descriptive_Comparison.csv")
)


############################################################
# ==========================================================
#
# FINAL SUMMARY
#
# ==========================================================
############################################################



############################################################
# 83. ANALYSIS SUMMARY
############################################################

analysis_summary <- tibble(
  
  Metric = c(
    
    "DIA-NN unique normalized Runs",
    
    "Metadata unique normalized Runs",
    
    "Matched QC-retained Runs",

    "DIA-NN rows identified as cRAP or CON__",

    "Unique cRAP or CON__ protein groups removed",
    
    "Proteins before completeness filtering",
    
    "Proteins >=70% overall",
    
    "Successful zonation LMM",
    
    "Failed zonation LMM",
    
    "Singular zonation LMM",
    
    "Strong central zonated proteins",
    
    "Strong portal zonated proteins",
    
    "Proteins >=70% in both GF and SPF",
    
    "Significant GF-SPF zonation interactions",

    "Proteins in Spatial_S heatmap",

    "Multi-page heatmap pages per group"
    
  ),
  
  
  Value = c(
    
    length(
      diann_keys
    ),
    
    length(
      metadata_keys
    ),
    
    length(
      matched_keys
    ),

    nrow(
      crap_information
    ),

    n_distinct(
      crap_information$Protein.Group
    ),
    
    nrow(
      completeness_table
    ),
    
    nrow(
      retained_proteins
    ),
    
    sum(
      !is.na(
        zonation_results$P_value
      )
    ),
    
    sum(
      is.na(
        zonation_results$P_value
      )
    ),
    
    sum(
      zonation_results$Singular ==
        TRUE,
      na.rm = TRUE
    ),
    
    sum(
      zonation_results$Zonation ==
        "Central",
      na.rm = TRUE
    ),
    
    sum(
      zonation_results$Zonation ==
        "Portal",
      na.rm = TRUE
    ),
    
    length(
      proteins_both70
    ),
    
    sum(
      GFSPF_results$Q_interaction <
        0.05,
      na.rm = TRUE
    ),

    n_distinct(
      heatmap_selection$Protein.Group
    ),

    if (
      HEATMAP_CREATE_MULTIPAGE
    ) {
      n_heatmap_pages
    } else {
      0
    }
    
  )
  
)



############################################################
# 84. SAVE SUMMARY
############################################################

write_csv(
  
  analysis_summary,
  
  file.path(
    OUTPUT_DIR,
    "Analysis_Summary.csv"
  )
  
)



############################################################
# 85. FINAL OUTPUT
############################################################

cat(
  "\n==============================================\n"
)

cat(
  "ANALYSIS COMPLETED\n"
)

cat(
  "==============================================\n\n"
)


print(
  analysis_summary
)

cat(
  "\nGF/SPF independent analysis summary:\n"
)

print(
  GF_SPF_separate_summary
)


cat(
  "\nOutput folder:\n"
)

cat(
  OUTPUT_DIR,
  "\n"
)

invisible(list(
  summary = analysis_summary,
  separate_summary = GF_SPF_separate_summary,
  output_dir = normalizePath(OUTPUT_DIR, winslash = "/", mustWork = TRUE)
))
}
