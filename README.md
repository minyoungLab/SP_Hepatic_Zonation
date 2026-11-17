# Interactive Spatial Proteomics Analysis of Hepatic CV–PV Zonation

This R Shiny-based spatial bioimage analysis platform integrates quality control (QC), group comparison, hepatic lobule central vein-to-portal vein (CV→PV) zonation, protein trajectory analysis, and histology-based spatial mapping into a unified workflow for spatial proteomics data from germ-free (GF) and specific pathogen-free (SPF) mouse liver tissues. By automatically registering region of interest (ROI) coordinates on histological images with mass spectrometry-based protein quantification data, the platform visualizes spatial protein expression patterns in image form.

The analytical framework applies a linear mixed-effects model incorporating sample-level random intercepts to control for baseline expression differences across samples, enabling precise analysis across diverse spatial proteomics specimens. At its core, the method quantifies both the direction and magnitude of protein expression gradients along the zonation trajectory. By detecting subtle zonation patterns that are difficult to observe with conventional zonation isolation methods, the platform provides deeper insights into the spatially dependent functions of liver tissue.

Although currently configured for CV→PV spatial zonation analysis in GF and SPF mouse livers, the application can be readily adapted to other tissues by adjusting the ROI parsing rules and analysis settings to investigate how spatial organization drives changes in biological function.

# Demo
If you would like to view the live application, please visit the following webpage:

(https://sp-hepatic-zonation.shinyapps.io/SP-liver/)

## License

MIT License. See `LICENSE` for details.

## Citation

Please cite the Zenodo:

> C. [https://10.5281/zenodo.24304897].
