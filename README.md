# trajLagged

GAM-based Preprocessing and Bootstrap Analysis for Single-Cell Trajectory Data

## Overview

The `trajLagged` package provides functions for GAM-based smoothing and preprocessing of single-cell trajectory data, along with bootstrap analysis methods for statistical inference. The package includes functions for data transformation, GAM smoothing, weighted sum calculations, and bootstrap confidence interval estimation with pseudotime density weighting.

## Installation

You can install the development version of trajLagged from GitHub:

```r
# Install devtools if you haven't already
install.packages("devtools")

# Install trajLagged
devtools::install_github("your-username/trajLagged")
```

## Key Features

- **GAM-based smoothing**: Smooth trajectory data using Generalized Additive Models
- **Flexible preprocessing**: Handle both binned and unbinned data
- **Bootstrap analysis**: Perform statistical inference with confidence intervals
- **Pseudotime weighting**: Use density-based weighting for more robust bootstrap sampling
- **GAM flatness testing**: Statistical tests to determine if trajectories are flat
- **Efficient implementation**: Optimized for performance with large datasets

## Main Functions

### Preprocessing Functions

- `preprocess_gam()`: GAM-based preprocessing from single-cell experiment objects
- `preprocess_gam_from_df()`: GAM-based preprocessing from data frames
- `transform_data()`: Transform single-cell experiment data (needs customization)

### Analysis Functions

- `sumWeight()`: Calculate weighted sums with normalization
- `get_value_gam()`: Get analysis values from GAM results
- `calculate_pseudotime_weights()`: Calculate pseudotime density weights

### GAM Flatness Testing Functions

- `test_gam_flatness()`: Comprehensive statistical tests for GAM flatness
- `test_trajectory_flatness()`: Test flatness for trajectory GAMs
- `plot_gam_flatness()`: Diagnostic plots for GAM flatness testing

### Pipeline Functions

- `pipeline_gam()`: Complete GAM analysis pipeline
- `pipeline_from_df()`: Pipeline starting from data frames
- `bootstrap_pipeline_gam()`: Bootstrap analysis with statistical inference

## Usage Example

```r
library(trajLagged)

# Basic GAM preprocessing
result <- preprocess_gam(sce1, sce2, gene = "GENE1", binned = FALSE)

# Complete pipeline analysis
analysis_result <- pipeline_gam(sce1, sce2, gene = "GENE1", binned = FALSE)

# Bootstrap analysis with confidence intervals
bootstrap_result <- bootstrap_pipeline_gam(
  sce1, sce2,
  gene = "GENE1",
  binned = FALSE,
  n_bootstrap = 100,
  m_prop = 0.8
)

# Bootstrap analysis with GAM flatness testing
bootstrap_with_flatness <- bootstrap_pipeline_gam(
  sce1, sce2,
  gene = "GENE1",
  n_bootstrap = 100,
  test_flatness = TRUE,
  flatness_method = "all"
)

# Whole gene analysis with flatness testing
whole_gene_results <- run_whole_gene(
  sce1, sce2,
  geneList = c("GENE1", "GENE2", "GENE3"),
  bootstrap = TRUE,
  test_flatness = TRUE,
  flatness_method = "summary"
)

# Access results
print(bootstrap_result$original_value)
print(bootstrap_result$ci_lower)
print(bootstrap_result$ci_upper)
print(bootstrap_result$excludes_zero)

# Access flatness test results
if (!is.null(bootstrap_with_flatness$flatness_tests)) {
  print(bootstrap_with_flatness$flatness_tests)
}
```

## Data Requirements

The package expects single-cell experiment objects with:
- Pseudotime information stored in `$pseudotime`
- Gene expression data accessible through standard methods
- Compatible data structure for the `transform_data()` function

**Note**: You'll need to customize the `transform_data()` function based on your specific data structure.

## Dependencies

- R (>= 3.5.0)
- mgcv (for GAM fitting)
- dplyr (for data manipulation)
- stats (for statistical functions)

## License

MIT License

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## Citation

If you use this package in your research, please cite:

```
[Your citation information here]
``` 