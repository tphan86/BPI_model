# Bacteriophage–Bacteria Interaction (BPI) Model

MATLAB code for fitting a generalized Lotka-Volterra competition model to growth
experimental data from phage-challenged bacterial populations. This repository
accompanies the manuscript:

> **Phage pressure modulates competitive dynamics between susceptible and resistant
> bacterial populations: a mathematical modeling study**
> Phan, Shrestha, et al. *(under revision, Bulletin of Mathematical Biology)*

---

## Overview

When bacteria are exposed to bacteriophage, some cells may acquire resistance through
mutation or phenotypic switching. The resulting population contains two competing
subpopulations — phage-susceptible (S) and phage-resistant (R) cells — whose dynamics
depend on the multiplicity of infection (MOI). This code fits an ODE-based competition
model to growth data from four bacterial strain–phage combinations and provides tools
for formal model comparison and profile likelihood analysis.

---

## Repository structure

```
BPI_model/
├── PA103_Rd3_Cocktail/     # P. aeruginosa PA103 challenged with phage cocktail
├── PAK_Rd3_Cocktail/       # P. aeruginosa PAK challenged with phage cocktail
├── Y3650_Fern/             # P. larvae Y-3650 challenged with phage Fern IDv1
├── 25747_Fern/             # P. larvae 25747 challenged with phage Fern IDv1
├── README.md
└── LICENSE
```

Each subfolder is self-contained and follows the same internal structure:

```
<strain_phage>/
├── data/                   # Experimental OD600 data (loaded from Excel)
├── control_fit/            # Stage 1: logistic growth fit to phage-free controls
├── model_fit/              # Stage 2: full model fit to phage treatment data
├── profile_likelihood/     # Profile likelihood analysis (P. larvae strains only)
└── figures/                # Output figures
```

---

## The mathematical model

The model describes the coupled dynamics of phage-susceptible bacteria S(t) (CFU/mL),
phage-resistant bacteria R(t) (CFU/mL), and phage P(t) (PFU/mL) using a system of
ordinary differential equations:

```
dS/dt =  r·S·(1 − (a_ss·S + a_rs·R)/K)  −  δ·S  −  φ·P·S

dR/dt =  γ·r·R·(1 − (a_sr·S + a_rr·R)/K)  +  δ·S  −  ε·φ·P·R

dP/dt =  β·φ·(S + ε·R)·P  −  m·P
```

### Parameters

| Symbol   | Description                                        | Units              |
|----------|----------------------------------------------------|--------------------|
| r        | Intrinsic growth rate of susceptible bacteria      | h⁻¹                |
| γ        | Fitness cost of resistance (γ ∈ (0,1])             | —                  |
| K        | Carrying capacity (fixed at 10⁹ CFU/mL)            | CFU/mL             |
| a_ss     | Intraspecific competition coefficient of S         | —                  |
| a_rr     | Intraspecific competition coefficient of R         | —                  |
| a_rs     | Interspecific competition: effect of R on S        | —                  |
| a_sr     | Interspecific competition: effect of S on R        | —                  |
| δ        | Rate of resistance acquisition (S → R)             | h⁻¹                |
| φ        | Phage adsorption rate to S                         | mL·PFU⁻¹·h⁻¹      |
| ε        | Partial resistance efficiency (ε = 0: complete)    | —                  |
| β        | Phage burst size                                   | PFU·CFU⁻¹          |
| m        | Phage decay rate                                   | h⁻¹                |
| S₀, R₀   | Initial densities of S and R (fitted)              | CFU/mL             |
| P₀       | Initial phage density (set from experiment)        | PFU/mL             |

OD₆₀₀ measurements are converted to CFU/mL using the calibration factor
1 OD₆₀₀ ≈ 10⁹ CFU/mL before fitting.

---

## Parameterization: M1 vs. M2

Two parameterizations of the model are compared:

- **M1 (shared)** — the interspecific competition coefficients a_rs and a_sr are
  constrained to be identical across all MOI treatments. This parameterization uses
  5 globally shared parameters (r, a_ss, a_rs, a_sr, γ) plus treatment-specific
  initial conditions.

- **M2 (treatment-specific)** — a_rs and a_sr are fitted independently for each MOI
  treatment, reducing the shared parameter set to 3 (r, a_ss, γ). This
  parameterization allows the effective competitive landscape to vary with phage
  pressure.

Model selection between M1 and M2 is performed using AIC and BIC (see below).

---

## Fitting procedure

Parameter estimation is carried out in two stages.

**Stage 1 — Control fitting.**  
A logistic growth model is fitted to the phage-free control data using `lsqnonlin`
(MATLAB's trust-region-reflective nonlinear least squares solver). This step estimates
the bacterial growth parameters r, a_ss, and the initial susceptible cell density S₀,
with K fixed at 10⁹ CFU/mL.

**Stage 2 — Full model fitting.**  
The control-estimated parameters are fixed, and the full ODE model is fitted globally
to all four MOI treatment groups (Low, Medium, Medium-High, High) simultaneously by
minimizing a single combined sum of squared residuals. The ODE system is integrated
using MATLAB's stiff solver `ode15s` (RelTol = 10⁻⁸, AbsTol = 10⁻¹⁰) with
non-negativity constraints on all state variables. Multiple random initial conditions
are tested to reduce sensitivity to local minima.

### MOI treatment groupings

Replicates are grouped by MOI level and observed growth phenotype:

- **P. aeruginosa (PA103, PAK):** grouped by phage dilution (Low: 10⁻⁹–10⁻⁸;
  Medium: two-phase growth dilutions; Medium-High: delayed growth; High: suppression).
- **P. larvae (Y-3650, 25747):** grouped by observed growth pattern, with some
  individual replicates from the same nominal dilution assigned to different MOI
  categories based on their phenotype (delayed growth, two-phase growth, complete
  suppression).

### Goodness of fit

Fit quality is reported using:

- **R²** (coefficient of determination), computed per treatment and overall
- **RMSE** (root mean squared error)
- **AIC and BIC** for formal model comparison between M1 and M2 (P. larvae strains)

---

## Profile likelihood analysis

For the two *P. larvae* strains (Y-3650 and 25747), the code includes a profile
likelihood analysis to assess whether a_rs and a_sr can be shared across MOI
treatments. Each coefficient is fixed at 40 log-spaced values spanning the range of
treatment-specific estimates, and all other parameters are re-optimized at each grid
point using `lsqnonlin`. Profile R² curves are produced per treatment to identify
conflicting optima that would preclude a shared parameterization.

This analysis directly supports the manuscript's use of M2 (treatment-specific a_rs
and a_sr) for the *P. larvae* datasets.

---

## Requirements

| Software | Version |
|----------|---------|
| MATLAB   | R2023b or later |
| Global Optimization Toolbox | Required |
| Optimization Toolbox | Required (`lsqnonlin`) |

No additional third-party MATLAB toolboxes are required. Input data files are provided
in Excel format (`.xlsx`) in each strain subfolder.

---

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/tphan86/BPI_model.git
   ```

2. Open MATLAB and navigate to the subfolder for the strain of interest, e.g.:
   ```matlab
   cd BPI_model/Y3650_Fern
   ```

3. Run the control fit first:
   ```matlab
   run_control_fit.m
   ```
   This fits a logistic growth model to the phage-free control data and saves the
   estimated parameters for use in Stage 2.

4. Run the full model fit:
   ```matlab
   run_model_fit.m
   ```
   This fits M1 and M2 to all four MOI treatments simultaneously and outputs
   per-treatment R², RMSE, AIC, BIC, fitted parameter values, and figures.

5. (P. larvae only) Run the profile likelihood analysis:
   ```matlab
   run_profile_likelihood.m
   ```
   This generates profile R² curves for a_rs and a_sr across all four treatments.

Each script produces figures that are saved to the `figures/` subdirectory.

---

## Output figures

Each strain folder produces the following figures:

| Figure | Description |
|--------|-------------|
| Control fit | Logistic growth fit to phage-free data (OD₆₀₀ and CFU/mL scales) |
| M1 vs. M2 fit comparison | Side-by-side model fits for all four MOI treatments |
| Residual comparison | Residual plots for M1 and M2 across treatments |
| Optimal parameter values | Per-treatment optimal a_rs and a_sr under M1 and M2 |
| Profile likelihood curves | Per-treatment R² as a function of fixed a_rs or a_sr (P. larvae only) |

---

## Data

Experimental data are provided in Excel files within each strain subfolder. Time is
recorded in hours for *P. aeruginosa* and in minutes for *P. larvae*; the code converts
*P. larvae* time to hours automatically. OD₆₀₀ values are converted to CFU/mL using
the factor 1 OD₆₀₀ = 10⁹ CFU/mL prior to fitting.

---

## Citation

If you use this code, please cite the accompanying manuscript (citation will be updated
upon acceptance).

---

## License

MIT License. See `LICENSE` for details.

---

## Contact

Tuan Phan — tphan86@github  
For questions about the experimental data or biological interpretation, please refer
to the manuscript.
