############################################################
# Mouse Liver Spatial Proteomics
# Automatic ROI Metadata Generator
#
# Grid:
#   10 x 10
#   Four corner positions excluded
#   Total = 96 ROIs
#
# Required Fiji coordinates per image:
#   ROI01
#   ROI08
#   ROI89
#   CV
#   PV
#
# INPUT:
#   Image_Reference_Coordinates.csv
#
# OUTPUT:
#   ROI_metadata.csv
#   image-specific coordinate files
#   QC plots
############################################################


# Shiny/command-line entry point. Sourcing this file only defines the
# function; the analysis starts when run_roi_metadata() is called.
run_roi_metadata <- function(input_file, output_dir, n_row = 10L, n_col = 10L) {


############################################################
# 0. Packages
############################################################

library(tidyverse)


############################################################
# 1. USER SETTINGS
############################################################

# ----------------------------------------------------------
# Input file
# ----------------------------------------------------------

INPUT_FILE <- normalizePath(input_file, winslash = "/", mustWork = TRUE)


# ----------------------------------------------------------
# Output folder
# ----------------------------------------------------------

OUTPUT_DIR <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)

dir.create(
  OUTPUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)


# ----------------------------------------------------------
# Grid structure
# ----------------------------------------------------------

N_ROW <- as.integer(n_row)
N_COL <- as.integer(n_col)

if (is.na(N_ROW) || is.na(N_COL) || N_ROW < 2L || N_COL < 2L) {
  stop("n_row and n_col must be integers greater than 1.", call. = FALSE)
}


############################################################
# 2. READ IMAGE REFERENCE COORDINATES
############################################################

ref <- read_csv(
  INPUT_FILE,
  show_col_types = FALSE
)


############################################################
# 3. CHECK REQUIRED COLUMNS
############################################################

required_columns <- c(
  
  "Image",
  "Group",
  "Mouse",
  "Set",
  
  "ROI01_X",
  "ROI01_Y",
  
  "ROI08_X",
  "ROI08_Y",
  
  "ROI89_X",
  "ROI89_Y",
  
  "CV_X",
  "CV_Y",
  
  "PV_X",
  "PV_Y"
  
)


missing_columns <- setdiff(
  required_columns,
  names(ref)
)


if (length(missing_columns) > 0) {
  
  stop(
    
    paste0(
      "\nMissing columns:\n",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
    
  )
  
}


############################################################
# 4. BASIC CHECK
############################################################

cat(
  "\nNumber of images:",
  nrow(ref),
  "\n"
)

print(
  ref %>%
    select(
      Image,
      Group,
      Mouse,
      Set
    )
)


############################################################
# 5. CREATE MASTER 10 x 10 GRID
############################################################

master_grid <- expand_grid(
  
  Row = 1:N_ROW,
  
  Column = 1:N_COL
  
) %>%
  
  # --------------------------------------------------------
# Remove four corners
# --------------------------------------------------------

filter(
  
  !(Row == 1  & Column == 1),
  
  !(Row == 1  & Column == 10),
  
  !(Row == 10 & Column == 1),
  
  !(Row == 10 & Column == 10)
  
) %>%
  
  arrange(
    Row,
    Column
  ) %>%
  
  mutate(
    
    ROI_Number = row_number(),
    
    ROI = sprintf(
      "ROI%02d",
      ROI_Number
    )
    
  )


cat(
  "\nNumber of generated ROI positions:",
  nrow(master_grid),
  "\n"
)


if (nrow(master_grid) != 96) {
  
  stop(
    "Grid generation error: ROI number is not 96."
  )
  
}


############################################################
# 6. FUNCTION:
# Generate all 96 ROI coordinates for one image
############################################################

generate_roi_coordinates <- function(image_row) {
  
  
  ##########################################################
  # 6-1. Extract image-specific coordinates
  ##########################################################
  
  ROI01_X <- image_row$ROI01_X
  ROI01_Y <- image_row$ROI01_Y
  
  ROI08_X <- image_row$ROI08_X
  ROI08_Y <- image_row$ROI08_Y
  
  ROI89_X <- image_row$ROI89_X
  ROI89_Y <- image_row$ROI89_Y
  
  CV_X <- image_row$CV_X
  CV_Y <- image_row$CV_Y
  
  PV_X <- image_row$PV_X
  PV_Y <- image_row$PV_Y
  
  
  ##########################################################
  # 6-2. Horizontal vector
  #
  # ROI01 = Row1 Col2
  # ROI08 = Row1 Col9
  #
  # 7 intervals
  ##########################################################
  
  H_X <- (
    ROI08_X -
      ROI01_X
  ) / 7
  
  
  H_Y <- (
    ROI08_Y -
      ROI01_Y
  ) / 7
  
  
  ##########################################################
  # 6-3. Vertical vector
  #
  # ROI01 = Row1 Col2
  # ROI89 = Row10 Col2
  #
  # 9 intervals
  ##########################################################
  
  V_X <- (
    ROI89_X -
      ROI01_X
  ) / 9
  
  
  V_Y <- (
    ROI89_Y -
      ROI01_Y
  ) / 9
  
  
  ##########################################################
  # 6-4. Generate all coordinates
  ##########################################################
  
  out <- master_grid %>%
    
    mutate(
      
      Image = image_row$Image,
      
      Group = image_row$Group,
      
      Mouse = image_row$Mouse,
      
      Set = image_row$Set,
      
      
      # ----------------------------------------------------
      # Coordinate calculation
      #
      # ROI01 = Row1 Col2
      # ----------------------------------------------------
      
      X =
        ROI01_X +
        (Column - 2) * H_X +
        (Row - 1) * V_X,
      
      Y =
        ROI01_Y +
        (Column - 2) * H_Y +
        (Row - 1) * V_Y,
      
      
      # ----------------------------------------------------
      # CV/PV coordinates
      # ----------------------------------------------------
      
      CV_X =
        CV_X,
      
      CV_Y =
        CV_Y,
      
      PV_X =
        PV_X,
      
      PV_Y =
        PV_Y,
      
      
      # ----------------------------------------------------
      # Distance to CV
      # ----------------------------------------------------
      
      d_CV = sqrt(
        
        (X - CV_X)^2 +
          (Y - CV_Y)^2
        
      ),
      
      
      # ----------------------------------------------------
      # Distance to PV
      # ----------------------------------------------------
      
      d_PV = sqrt(
        
        (X - PV_X)^2 +
          (Y - PV_Y)^2
        
      ),
      
      
      # ----------------------------------------------------
      # CV → PV normalized spatial score
      #
      # CV = 0
      # PV = 1
      # ----------------------------------------------------
      
      Spatial_S =
        d_CV /
        (
          d_CV +
            d_PV
        ),
      
      
      # ----------------------------------------------------
      # Save vector information
      # useful for QC
      # ----------------------------------------------------
      
      Horizontal_Step_X =
        H_X,
      
      Horizontal_Step_Y =
        H_Y,
      
      Vertical_Step_X =
        V_X,
      
      Vertical_Step_Y =
        V_Y
      
    )
  
  
  return(
    out
  )
  
}


############################################################
# 7. GENERATE ALL 8 IMAGE COORDINATES
############################################################

all_roi_metadata <- map_dfr(
  
  seq_len(
    nrow(ref)
  ),
  
  function(i) {
    
    generate_roi_coordinates(
      ref[i, ]
    )
    
  }
  
)


############################################################
# 8. CHECK TOTAL ROI NUMBER
############################################################

expected_n <- nrow(ref) * 96


cat(
  "\nExpected total ROI:",
  expected_n,
  "\n"
)


cat(
  "Generated total ROI:",
  nrow(all_roi_metadata),
  "\n"
)


if (
  nrow(all_roi_metadata) !=
  expected_n
) {
  
  stop(
    "Unexpected number of ROI rows."
  )
  
}


############################################################
# 9. CHECK SPATIAL SCORE RANGE
############################################################

if (
  any(
    all_roi_metadata$Spatial_S < 0 |
    all_roi_metadata$Spatial_S > 1,
    na.rm = TRUE
  )
) {
  
  stop(
    "Spatial_S outside 0-1 range detected."
  )
  
}


############################################################
# 10. CREATE RUN COLUMN
#
# IMPORTANT:
# Modify this section according to your DIA-NN Run naming.
#
# Example:
# GF1-2_ROI01
############################################################

all_roi_metadata <- all_roi_metadata %>%
  
  mutate(
    
    Run = paste0(
      Set,
      "_",
      ROI
    )
    
  )


############################################################
# 11. FINAL COLUMN ORDER
############################################################

ROI_metadata <- all_roi_metadata %>%
  
  select(
    
    Run,
    
    Image,
    Group,
    Mouse,
    Set,
    
    ROI,
    ROI_Number,
    
    Row,
    Column,
    
    X,
    Y,
    
    CV_X,
    CV_Y,
    
    PV_X,
    PV_Y,
    
    d_CV,
    d_PV,
    
    Spatial_S,
    
    Horizontal_Step_X,
    Horizontal_Step_Y,
    
    Vertical_Step_X,
    Vertical_Step_Y
    
  )


############################################################
# 12. SAVE FINAL ROI_metadata.csv
############################################################

write_csv(
  
  ROI_metadata,
  
  file.path(
    OUTPUT_DIR,
    "ROI_metadata.csv"
  )
  
)


############################################################
# 13. SAVE IMAGE-SPECIFIC FILES
############################################################

split_metadata <- split(
  
  ROI_metadata,
  
  ROI_metadata$Image
  
)


walk2(
  
  split_metadata,
  
  names(split_metadata),
  
  function(df, image_name) {
    
    write_csv(
      
      df,
      
      file.path(
        
        OUTPUT_DIR,
        
        paste0(
          image_name,
          "_ROI_coordinates.csv"
        )
        
      )
      
    )
    
  }
  
)


############################################################
# 14. IMAGE SUMMARY TABLE
############################################################

image_summary <- ROI_metadata %>%
  
  group_by(
    Image,
    Group,
    Mouse,
    Set
  ) %>%
  
  summarise(
    
    N_ROI =
      n(),
    
    Spatial_S_Min =
      min(
        Spatial_S,
        na.rm = TRUE
      ),
    
    Spatial_S_Max =
      max(
        Spatial_S,
        na.rm = TRUE
      ),
    
    Spatial_S_Median =
      median(
        Spatial_S,
        na.rm = TRUE
      ),
    
    Mean_dCV =
      mean(
        d_CV,
        na.rm = TRUE
      ),
    
    Mean_dPV =
      mean(
        d_PV,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
    
  )


write_csv(
  
  image_summary,
  
  file.path(
    OUTPUT_DIR,
    "Image_Spatial_Summary.csv"
  )
  
)


############################################################
# 15. CREATE QC PLOTS AUTOMATICALLY
############################################################

QC_DIR <- file.path(
  OUTPUT_DIR,
  "QC_Plots"
)


dir.create(
  QC_DIR,
  showWarnings = FALSE
)


image_names <- unique(
  ROI_metadata$Image
)


for (
  image_name in image_names
) {
  
  
  ##########################################################
  # Select one image
  ##########################################################
  
  df_plot <- ROI_metadata %>%
    
    filter(
      Image == image_name
    )
  
  
  ##########################################################
  # Plot grid + CV + PV
  ##########################################################
  
  p <- ggplot(
    
    df_plot,
    
    aes(
      x = X,
      y = Y
    )
    
  ) +
    
    geom_point(
      
      aes(
        color = Spatial_S
      ),
      
      size = 4
      
    ) +
    
    geom_text(
      
      aes(
        label = ROI_Number
      ),
      
      size = 2.5
      
    ) +
    
    geom_point(
      
      aes(
        x = CV_X,
        y = CV_Y
      ),
      
      inherit.aes = FALSE,
      
      size = 6,
      
      shape = 21,
      
      fill = "white"
      
    ) +
    
    geom_text(
      
      aes(
        x = CV_X,
        y = CV_Y,
        label = "CV"
      ),
      
      inherit.aes = FALSE,
      
      size = 5,
      
      fontface = "bold",
      
      vjust = -1
      
    ) +
    
    geom_point(
      
      aes(
        x = PV_X,
        y = PV_Y
      ),
      
      inherit.aes = FALSE,
      
      size = 6,
      
      shape = 21,
      
      fill = "white"
      
    ) +
    
    geom_text(
      
      aes(
        x = PV_X,
        y = PV_Y,
        label = "PV"
      ),
      
      inherit.aes = FALSE,
      
      size = 5,
      
      fontface = "bold",
      
      vjust = -1
      
    ) +
    
    coord_fixed() +
    
    # Fiji coordinate system
    scale_y_reverse() +
    
    labs(
      
      title =
        paste0(
          image_name,
          " | ROI spatial coordinates"
        ),
      
      subtitle =
        "Spatial_S: CV = 0 → PV = 1",
      
      x =
        "X",
      
      y =
        "Y",
      
      color =
        "Spatial_S"
      
    ) +
    
    theme_classic(
      base_size = 12
    )
  
  
  ##########################################################
  # Save QC plot
  ##########################################################
  
  ggsave(
    
    filename =
      file.path(
        QC_DIR,
        paste0(
          image_name,
          "_ROI_Spatial_QC.png"
        )
      ),
    
    plot =
      p,
    
    width =
      8,
    
    height =
      8,
    
    dpi =
      300
    
  )
  
}


############################################################
# 16. COMBINED FACET QC
############################################################

p_all <- ggplot(
  
  ROI_metadata,
  
  aes(
    x = X,
    y = Y,
    color = Spatial_S
  )
  
) +
  
  geom_point(
    size = 2
  ) +
  
  geom_text(
    
    aes(
      label = ROI_Number
    ),
    
    size = 1.8
    
  ) +
  
  geom_point(
    
    aes(
      x = CV_X,
      y = CV_Y
    ),
    
    inherit.aes = FALSE,
    
    size = 4,
    
    shape = 21,
    
    fill = "white"
    
  ) +
  
  geom_text(
    
    aes(
      x = CV_X,
      y = CV_Y,
      label = "CV"
    ),
    
    inherit.aes = FALSE,
    
    size = 3.5,
    
    fontface = "bold"
    
  ) +
  
  geom_point(
    
    aes(
      x = PV_X,
      y = PV_Y
    ),
    
    inherit.aes = FALSE,
    
    size = 4,
    
    shape = 21,
    
    fill = "white"
    
  ) +
  
  geom_text(
    
    aes(
      x = PV_X,
      y = PV_Y,
      label = "PV"
    ),
    
    inherit.aes = FALSE,
    
    size = 3.5,
    
    fontface = "bold"
    
  ) +
  
  facet_wrap(
    ~ Image,
    scales = "fixed"
  ) +
  
  scale_y_reverse() +
  
  labs(
    
    title =
      "All liver spatial proteomics ROI grids",
    
    subtitle =
      "CV = 0 → PV = 1",
    
    x =
      "X",
    
    y =
      "Y",
    
    color =
      "Spatial_S"
    
  ) +
  
  theme_classic(
    base_size = 10
  )


ggsave(
  
  filename =
    file.path(
      OUTPUT_DIR,
      "All_Images_ROI_Spatial_QC.png"
    ),
  
  plot =
    p_all,
  
  width =
    14,
  
  height =
    10,
  
  dpi =
    300
  
)


############################################################
# 17. SPATIAL SCORE DISTRIBUTION
############################################################

p_distribution <- ggplot(
  
  ROI_metadata,
  
  aes(
    x = Spatial_S,
    fill = Group
  )
  
) +
  
  geom_histogram(
    
    bins = 20,
    
    alpha = 0.5,
    
    position = "identity"
    
  ) +
  
  facet_wrap(
    ~ Image
  ) +
  
  scale_x_continuous(
    
    limits =
      c(
        0,
        1
      )
    
  ) +
  
  labs(
    
    title =
      "Distribution of CV-PV spatial scores",
    
    x =
      "Spatial_S (CV = 0 → PV = 1)",
    
    y =
      "Number of ROIs"
    
  ) +
  
  theme_classic()


ggsave(
  
  filename =
    file.path(
      OUTPUT_DIR,
      "Spatial_S_Distribution.png"
    ),
  
  plot =
    p_distribution,
  
  width =
    12,
  
  height =
    8,
  
  dpi =
    300
  
)


############################################################
# 18. FINAL SUMMARY
############################################################

cat(
  "\n==============================================\n"
)

cat(
  "ROI metadata generation completed.\n"
)

cat(
  "==============================================\n\n"
)


cat(
  "Images:",
  n_distinct(
    ROI_metadata$Image
  ),
  "\n"
)


cat(
  "Total ROIs:",
  nrow(
    ROI_metadata
  ),
  "\n"
)


cat(
  "ROI per image:",
  nrow(
    ROI_metadata
  ) /
    n_distinct(
      ROI_metadata$Image
    ),
  "\n"
)


cat(
  "\nOutput directory:\n",
  OUTPUT_DIR,
  "\n"
)


cat(
  "\nMain output:\n",
  file.path(
    OUTPUT_DIR,
    "ROI_metadata.csv"
  ),
  "\n"
)


print(
  image_summary
)

invisible(list(
  summary = image_summary,
  metadata_file = normalizePath(
    file.path(OUTPUT_DIR, "ROI_metadata.csv"),
    winslash = "/", mustWork = TRUE
  ),
  output_dir = normalizePath(OUTPUT_DIR, winslash = "/", mustWork = TRUE)
))
}
