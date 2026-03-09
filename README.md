# Dose Response Pro v18.1

**A browser-based dose-response meta-analysis tool implementing the Greenland & Longnecker two-stage generalized least squares (GLS) method.**

No installation required. Open `index.html` in any modern browser to begin.

**Live demo:** [https://mahmood726-cyber.github.io/dose-response-pro/](https://mahmood726-cyber.github.io/dose-response-pro/)

## Quick Start

1. Open `index.html` in Chrome, Firefox, Safari, or Edge
2. Click **Load Demo Data** or import your own CSV
3. Click **Run Analysis**
4. Explore Results, Plots, Sensitivity, and Bias tabs

Or serve locally:
```bash
python -m http.server 8000
# Visit http://localhost:8000
```

## Features

| Feature | Description |
|---------|-------------|
| **GLS Method** | Greenland & Longnecker two-stage dose-response meta-analysis |
| **Linear Model** | Simple linear dose-response trend |
| **Quadratic Model** | Non-linear dose-response with curvature |
| **Spline Model** | Restricted cubic spline subgroup analysis (exploratory) |
| **Sensitivity Analysis** | Leave-one-out with change metrics and influence diagnostics |
| **Subgroup Analysis** | Stratified heterogeneity exploration |
| **Publication Bias** | Funnel plots and asymmetry diagnostics |
| **CSV Import** | Flexible column detection with smart parsing |
| **R Code Export** | `metafor`-compatible reproducibility export |
| **Interactive Plots** | Dose-response curves, forest summaries, funnel/bias plots |

## Validation

Dose Response Pro has been validated against R packages (`dosresmeta` 2.2.0, `metafor` 4.8.0, `mvmeta` 1.0.3) on R 4.5.2:

- **3/3 validation scenarios passed** (linear, quadratic, multi-study)
- Maximum coefficient difference: 7.56e-06
- Maximum SE difference: 4.38e-06
- Maximum tau-squared difference: 6.65e-09

Validation script: [`tests/validate_dose_response_pro.R`](tests/validate_dose_response_pro.R)
Validation results: [`tests/r_validation_results.json`](tests/r_validation_results.json)

## Sample Datasets

Teaching datasets are included in `sample_data/`:

| Dataset | Demonstrates |
|---------|--------------|
| `linear_trend_teaching.csv` | Clear linear dose-response |
| `u_shaped_curve_teaching.csv` | Non-linear U-shaped curve |
| `high_heterogeneity_teaching.csv` | High I-squared, subgroup differences |
| `saturation_effect_teaching.csv` | Exponential saturation |
| `edge_case_zero_cases.csv` | Handling zero events |

## Comparison with R Packages

| Feature | Dose Response Pro | dosresmeta | metafor | mvmeta |
|---------|-------------------|------------|---------|--------|
| GLS Method | Yes | Yes | Manual | Yes |
| Interactive UI | Yes | No | No | No |
| Real-time plots | Yes | No | No | No |
| No installation | Yes | No | No | No |
| Sensitivity analysis | GUI | Manual | Manual | Manual |
| R code export | Yes | N/A | N/A | N/A |
| Cost | Free | Free | Free | Free |

## Repository Structure

```
dose-response-pro/
  index.html                          # Main application
  sample_dose_response_data.csv       # Default demo dataset
  sample_data/                        # Teaching datasets
  tests/
    validate_dose_response_pro.R      # R validation script
    r_validation_results.json         # Validation output
  docs/
    Getting_Started_Guide.md
    GLS_Method_Documentation.md
    Validation_Results_v18.1_Corrected.md
    Complete_Documentation.md
    Computational_Complexity.md
```

## Documentation

| For | Read |
|-----|------|
| First-time users | [Getting Started Guide](docs/Getting_Started_Guide.md) |
| Methodology | [GLS Method Documentation](docs/GLS_Method_Documentation.md) |
| Validation | [Validation Results](docs/Validation_Results_v18.1_Corrected.md) |
| Full API | [Complete Documentation](docs/Complete_Documentation.md) |
| Performance | [Computational Complexity](docs/Computational_Complexity.md) |

## Citation

```bibtex
@software{dose_response_pro_v18,
  title = {Dose Response Pro v18.1: A Browser-Based Dose-Response Meta-Analysis Tool},
  author = {Ahmad, Mahmood and Kumar, Niraj and Dar, Bilaal and Khan, Laiba and Woo, Andrew},
  year = {2026},
  url = {https://github.com/mahmood726-cyber/dose-response-pro},
  version = {18.1},
  license = {MIT}
}
```

## License

MIT License. See [LICENSE](LICENSE) for details.
