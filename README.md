# The PARNI Scheme

A collection of algorithms implementing the **Point-wise Adaptive Random Neighbourhood Informed (PARNI)** Markov chain Monte Carlo framework for high-dimensional Bayesian variable selection and structure learning.

## Papers

This repository contains code accompanying the following publications:

1. **Adaptive random neighbourhood informed Markov chain Monte Carlo for high-dimensional Bayesian variable selection**
   Xitong Liang, Samuel Livingstone, Jim Griffin.
   *Statistics and Computing*, 32(5), Article 84, 2022.
   [Paper](https://link.springer.com/article/10.1007/s11222-022-10137-8)

2. **Adaptive MCMC for Bayesian Variable Selection in Generalised Linear Models and Survival Models**
   Xitong Liang, Samuel Livingstone, Jim Griffin.
   *Entropy*, 25(9), 1310, 2023.
   [Paper](https://www.mdpi.com/1099-4300/25/9/1310)

3. **Structure Learning with Adaptive Random Neighborhood Informed MCMC**
   Alberto Caron\*, Xitong Liang\*, Samuel Livingstone, Jim Griffin.
   *NeurIPS 2023*.
   [Paper](https://proceedings.neurips.cc/paper_files/paper/2023/file/8027ace571384361920665f1d1b69758-Paper-Conference.pdf)

   \*Equal contribution.

## Folder Guide

| Folder | Description |
|--------|-------------|
| [`R/`](R/) | PARNI, ASI, and ADS samplers for Bayesian variable selection in linear regression (Paper 1). Includes simulation studies and real data experiments. |
| [`GLM_PARNI/`](GLM_PARNI/) | Extension to GLMs and survival models (Paper 2) — logistic, Cox, and Weibull regression with approximate Laplace approximation and Polya-Gamma data augmentation. |
| [`Structure-Learning/`](Structure-Learning/) | PARNI adapted for DAG structure learning (Paper 3) — explores adjacency matrices over DAG space with PC-algorithm skeleton constraints. Includes comparisons with order and partition MCMC. |
| [`python/`](python/) | Python implementation of the linear regression BVS samplers as a modular package (`bvs/`), with PARNI, ASI, and ADS proposals, parallel tempering, and random g-prior support. |

See the README in each folder for detailed usage instructions and file descriptions.

## Getting Started

### Requirements

- R (>= 4.0)
- R packages: `MASS`, `survival` (for Cox models), `pgdraw` (for Polya-Gamma sampling)

### Quick Start (Linear Regression Variable Selection)

```r
# Set working directory to R/
setwd("R")

# Load all sampler functions
source("main functions.R")

# Set up hyperparameters
source("make_hyper_par.R")

# Run a simulation example
source("simulation/simulation_example.R")
```

### Quick Start (GLM Variable Selection)

```r
# Set working directory to GLM_PARNI/
setwd("GLM_PARNI")

# Load all GLM functions
source("main.R")
```

### Quick Start (Structure Learning)

```r
# Set working directory to Structure-Learning/
setwd("Structure-Learning")

# Run examples
source("0_examples.R")
```

## Citation

If you use this code, please cite the relevant paper(s):

```bibtex
@article{liang2022adaptive,
  title={Adaptive random neighbourhood informed {M}arkov chain {M}onte {C}arlo for high-dimensional {B}ayesian variable selection},
  author={Liang, Xitong and Livingstone, Samuel and Griffin, Jim},
  journal={Statistics and Computing},
  volume={32},
  number={5},
  pages={84},
  year={2022},
  publisher={Springer}
}

@article{liang2023adaptive,
  title={Adaptive {MCMC} for {B}ayesian Variable Selection in Generalised Linear Models and Survival Models},
  author={Liang, Xitong and Livingstone, Samuel and Griffin, Jim},
  journal={Entropy},
  volume={25},
  number={9},
  pages={1310},
  year={2023},
  publisher={MDPI}
}

@inproceedings{caron2023structure,
  title={Structure Learning with Adaptive Random Neighborhood Informed {MCMC}},
  author={Caron, Alberto and Liang, Xitong and Livingstone, Samuel and Griffin, Jim},
  booktitle={Advances in Neural Information Processing Systems},
  volume={36},
  year={2023}
}
```

## License

Please contact the authors for licensing information.
