# GERT-Azeotrope-Part2
Scripts relacionados ao paper Azeotrope- Part2
# GERT Azeotrope Benchmark — Paper 2

Reproducible code, figures, and data for *The GERT Cohesive-Fraction Criterion
Predicts Azeotropes Across the Water and Organic–Organic Domains: A Benchmark
Without System-Specific Fitting Against Quantum Solvation Models* (V. P. Dutra).

## Overview

A single dimensionless criterion — a system sits at a phase boundary when the
cohesive fraction **φ = f_M / (f_M + f_L) = ½** — is applied to two domains,
with the cohesive and dispersive forces identified per domain. All calibration
constants are held fixed from the foundational study (Paper 1); **nothing is
tuned to the azeotrope data** (no system-specific fitting).

- **Water domain (asymmetric host network).** φ classifies 89 of 91 unambiguous
  aqueous systems without new parameters (existence + sign + composition via the
  Margules identity W = −RT ln γ∞; the asymmetric aqueous amines via a
  dominant-channel form of f_M).
- **Organic–organic domain (symmetric).** The deviation reduces to the exact
  order parameter **ΔC = (α₁−α₂)(β₂−β₁)**, gated for existence by a Trouton
  volatility term |ln r| = (88/R)|ΔT_b|/T_m.

## Key results (organic–organic, n = 2,141)

| Metric                  | GERT (1 law, 0 fitted params) | SM12 (quantum solvation) |
| ----------------------- | ----------------------------- | ------------------------ |
| Raw accuracy            | 73.1%                         | 73.9%                    |
| Balanced accuracy       | 73.8%                         | 48.9%                    |
| Recall none / pos / neg | 73 / 74 / 75                  | 84 / 61 / 1              |

GERT ties SM12 on raw accuracy, leads by 24.9 points on balanced accuracy, and
recovers 75% of the maximum-boiling (negative) azeotropes that SM12 misses (1%).

## Files

### Scripts (`scripts/`)

| File                             | Purpose                                                     |
| -------------------------------- | ----------------------------------------------------------- |
| `gert_orgorg_unified.py`         | Canonical organic–organic law; prints the headline metrics  |
| `gert_orgorg_negatives.py`       | Standalone maximum-boiling (negative) classifier            |
| `gert_si_orgorg_table.py`        | Writes the full 2,141-system SI table (CSV)                 |
| `gert_water_level1.py`           | Water existence-level classification + Fig. 3               |
| `gert_water_level2.py`           | Water composition profiles (Fig. 2)                         |
| `gert_water_homologous.py`       | Homologous-series φ crossing (Fig. 1)                       |
| `gert_water_type12_retab.py`     | Type 1+2 retabulation (amines, dominant channel)            |
| `gert_si_water_table.py`         | Writes the curated water SI table (CSV)                     |
| `_kt.py`, `_hansen.py`, `_tb.py` | Tabulated descriptors (Kamlet–Taft, Hansen, boiling points) |

### Figures (`figures/`)

`fig1_homologous_phi.png`, `fig2_phi_profiles.png`, `fig3_classification_plane.png`
(water); `fig_orgorg_type_plane.png`, `fig_orgorg_negative_plane.png`,
`fig_orgorg_gert_vs_smx.png` (organic–organic). `fig_unified_criterion.png` is an
optional cross-domain figure (not used in the main text).

### Data / SI (`data/`)

`GERT_Paper2_SI_orgorg_2141.csv` — all 2,141 organic–organic systems with
φ_type, ΔC, |ln r|, GERT and SM12 predictions, observed behaviour.
`GERT_Paper2_SI_water_curated.csv` — curated water existence set.

### Paper (`paper/`)

`GERT_Paper2_COMPLETE.tex` (single-file LaTeX), `GERT_Paper2_COMPLETE.pdf`,
`GERT_Paper2_refs.bib`, `GERT_Paper2_SI.pdf`.

## Reproduction

```bash
pip install pandas numpy matplotlib openpyxl
# place the benchmark spreadsheet (Roese et al. 2019 SI) in scripts/:
#   ie9b03694_si_002.xlsx
cd scripts
python3 gert_orgorg_unified.py     # -> GERT raw 73.1% / SM12 73.9%
python3 gert_si_orgorg_table.py    # -> ../data/GERT_Paper2_SI_orgorg_2141.csv
python3 gert_water_level1.py       # -> water Fig. 3 + existence verdicts
```

## Data sources

- **Benchmark / SMx predictions:** Roese, Margulis, Schmidt, Uzat, Heintz,
  Paluch, *Ind. Eng. Chem. Res.* **58** (2019) 22626–22632, DOI
  10.1021/acs.iecr.9b03694 (SI: `ie9b03694_si_002.xlsx`).
- **Kamlet–Taft α, β:** Marcus, *The Properties of Solvents* (Wiley, 1998).
- **Hansen δ_d, δ_p, V_m:** Hansen, *Hansen Solubility Parameters* (CRC, 2007).
- **Water γ∞:** Brouwer & Schuur, *Sep. Purif. Technol.* **272** (2021) 118727,
  built on the Kojima database (*Fluid Phase Equilib.* **131**, 1997).
- **Observed azeotropy:** Gmehling et al., *Azeotropic Data* (Wiley-VCH).

## Note on the water table

The curated water table (`GERT_Paper2_SI_water_curated.csv`) reproduces the
existence set behind Fig. 3. The **complete 91-system enumeration** additionally
requires the per-solute water γ∞ values from the Brouwer 2021 SI (not bundled
here for licensing); with that file, `gert_water_level1.py` extends to the full
domain.

## License

Code and data: CC BY 4.0.
