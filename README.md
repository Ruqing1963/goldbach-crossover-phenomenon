# The Crossover Phenomenon in Hardy-Littlewood Goldbach Formula

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18123132.svg)](https://doi.org/10.5281/zenodo.18123132)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Paper](https://img.shields.io/badge/Paper-PDF-red.svg)](https://github.com/Ruqing1963/goldbach-crossover-phenomenon/blob/main/paper/Chen_Goldbach_FINAL.pdf)

## Overview

This repository contains the computational code, analysis scripts, and manuscript for the research paper:

**"The Crossover Phenomenon in Hardy-Littlewood Goldbach Formula: Computational Evidence of Scale-Dependent Performance and Asymptotic Dominance"**

*by Ruqing Chen (2026)*

GUT Geoservice Inc., Montreal, Quebec, Canada

## Key Findings

- **Crossover phenomenon** at N ≈ 10⁵ where optimal method transitions from series expansion to logarithmic integral formulation
- **Three distinct computational regimes** with different error characteristics
- **Peak accuracy advantage** of 843-fold for integral method at N ≈ 5.5×10⁷
- **Non-monotonic error structure** revealing complex arithmetic resonances (advantage ratios vary from 2× to 843× depending on N's number-theoretic properties)
- **Factor-of-2 correction** in counting methodology resolving prior literature ambiguity

## Dataset

Complete dataset (21,511 verification points spanning N = 10³ to N = 10⁸) available at:

**Zenodo:** [https://doi.org/10.5281/zenodo.18123132](https://doi.org/10.5281/zenodo.18123132)

## Repository Structure

```
goldbach-crossover-phenomenon/
├── README.md                      # This file
├── LICENSE                        # MIT License
├── paper/                         # Research paper
│   ├── Chen_Goldbach_FINAL.pdf   # Final manuscript (PDF)
│   ├── manuscript.tex            # LaTeX source
│   └── figures/                  # Figures for paper
│       └── crossover_figure.png
├── code/                          # Analysis and computation code
│   ├── goldbach_counter.py       # Goldbach pair counting
│   ├── hardy_littlewood.py       # Hardy-Littlewood formula implementations
│   ├── data_collection.py        # Main data collection script
│   └── visualization.py          # Figure generation
├── data/                          # Sample data (full dataset on Zenodo)
│   └── sample_data.csv           # Representative sample (1000 points)
└── docs/                          # Documentation
    └── methodology.md            # Detailed methodology
```

## Quick Start

### Prerequisites

```bash
Python 3.10 or higher
NumPy >= 1.24.0
SciPy >= 1.10.0
matplotlib >= 3.7.0
pandas >= 2.0.0
```

### Installation

```bash
# Clone the repository
git clone https://github.com/Ruqing1963/goldbach-crossover-phenomenon.git
cd goldbach-crossover-phenomenon

# Install dependencies
pip install numpy scipy matplotlib pandas
```

### Basic Usage

```python
from code.hardy_littlewood import hardy_littlewood_series, hardy_littlewood_integral
from code.goldbach_counter import count_goldbach_ordered

# Example: Compare methods for N = 100,000
N = 100000
actual = count_goldbach_ordered(N)
series_pred = hardy_littlewood_series(N)
integral_pred = hardy_littlewood_integral(N)

print(f"Actual G(N): {actual}")
print(f"Series prediction: {series_pred:.2f}")
print(f"Integral prediction: {integral_pred:.2f}")
print(f"Series bias: {(series_pred - actual) / actual * 100:.2f}%")
print(f"Integral bias: {(integral_pred - actual) / actual * 100:.2f}%")
```

### Running Full Analysis

```bash
# Generate sample data (WARNING: takes several hours for full range)
python code/data_collection.py --max-n 1000000 --output sample_output.csv

# Generate visualizations
python code/visualization.py --input sample_output.csv --output figures/
```

## Methodology

### Counting Method

We employ **ordered-pair** counting consistent with Hardy-Littlewood's circle method derivation:

- For N = 10: ordered pairs (3,7), (5,5), (7,3) give G(10) = 3
- This resolves a factor-of-2 definitional ambiguity in prior literature

### Hardy-Littlewood Implementations

**Series Expansion (4th order):**
```
r(N) ≈ 2C₂S(N) × N/log²(N) × [1 + a₁/log(N) + a₂/log²(N) + ...]
```
where C₂ ≈ 0.6601618158 (twin prime constant)

**Logarithmic Integral:**
```
r(N) ≈ 2C₂S(N) × ∫₂^(N-2) dt / [log(t)log(N-t)]
```
Evaluated using adaptive Gaussian quadrature

### Singular Series S(N)

Computed using exact formula accounting for prime factorization of N/2:
```
S(N) = ∏ (1 - 1/(p-1)²) / (1 - 2/(p-1)²)
       p|N, p odd
```

Values typically range from 1.0 to 3.0 depending on N's arithmetic properties.

## Results Summary

### Performance by Scale

| N Range | Series \|Bias\| | Integral \|Bias\| | Winner | Adv. Ratio |
|---------|----------------|------------------|--------|-----------|
| 10³-10⁴ | 4.72% | 5.38% | Series | 0.88× |
| 10⁴-10⁵ | 1.75% | 1.89% | Series | 0.93× |
| 10⁵-10⁶ | 0.86% | 0.85% | Parity | 1.01× |
| 10⁶-10⁷ | 0.72% | 0.32% | Integral | 2.25× |
| 10⁷-10⁸ | 0.65% | 0.12% | Integral | 5.67× |

### Key Milestones

| N | G(N) | Series Bias | Integral Bias | Advantage Ratio |
|---|------|-------------|---------------|----------------|
| 10³ | 56 | +7.25% | +4.74% | 1.5× |
| 10⁵ | 1,620 | +0.06% | -0.04% | 1.5× |
| 10⁷ | 77,614 | -0.65% | +0.04% | 18.4× |
| 5.5×10⁷ | 451,658 | -0.61% | +0.10% | **843×** |
| 10⁸ | 582,800 | -0.60% | -0.06% | 9.8× |

### Non-Monotonic Behavior

The peak advantage occurs at N ≈ 5.5×10⁷ rather than at the maximum N tested (10⁸), revealing that asymptotic error structure is governed by complex arithmetic resonances (S(N) and ω(N) modulation) rather than simple power-law decay.

## Citation

If you use this code or dataset in your research, please cite:

```bibtex
@article{chen2026goldbach,
  title={The Crossover Phenomenon in Hardy-Littlewood Goldbach Formula: 
         Computational Evidence of Scale-Dependent Performance and 
         Asymptotic Dominance},
  author={Chen, Ruqing},
  journal={[Submitted to Mathematics of Computation]},
  year={2026},
  note={Code: \url{https://github.com/Ruqing1963/goldbach-crossover-phenomenon}, 
        Dataset: \url{https://doi.org/10.5281/zenodo.18123132}}
}
```

## License

- **Code:** MIT License (see LICENSE file)
- **Data:** CC BY 4.0 (available at Zenodo)
- **Paper:** All rights reserved (submitted manuscript)

## Contact

**Ruqing Chen**  
GUT Geoservice Inc.  
Montreal, Quebec, Canada  
Email: ruqing@hotmail.com

## Acknowledgments

This research received no specific funding. All computations performed on personal hardware. Special thanks to the OEIS community for maintaining sequence A006307 used for validation.

## Related Resources

- **OEIS A006307:** Number of Goldbach partitions of 2n  
  [https://oeis.org/A006307](https://oeis.org/A006307)

- **Hardy & Littlewood (1923):** Original asymptotic formula  
  Acta Mathematica, 44(1), 1-70

- **Oliveira e Silva et al. (2014):** Large-scale Goldbach verification  
  Mathematics of Computation, 83(288), 2033-2060
