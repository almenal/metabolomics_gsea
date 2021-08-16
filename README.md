# Evaluating Gene Set Analysis methods on metabolomics data
## GSEA and Globaltest

---

This repository contains scripts and ulitity functions to perform an evaluation of GSEA and Globaltest for metabolomics data.

## Usage 

### 1. Metabolite misidentification

To assess the effects of metabolite misidentification on each method, a simulation study was performed in which compounds were replaced by compounds with a similar mass (within 20ppm of the original mass). For example, the KEGG compound C00526 (Deoxyuridine) with mass 228.0746 can be replaced by C09317, with mass 228.0786.

```bash 
Rscript misIDreplicates.R
```

### 2. Semi-synthetic data simulations

The lack of a ground truth limits the comparison of these methods to mere degree of agreement. 
In order to overcome this, we established an approach based on simulations using semi-synthetic data. 
The first step is to produce a null dataset, which is achieved by permuting the phenotype labels.
Next, a pathway to be enriched is chosen according to a pre-defined criterion (see section `Metropolis coefficient`).
The entries of the data matrix corresponding to the metabolites in the selected pathway are enriched for the case group by adding random numbers drawn from a normal distribution (mean and standard deviation are parameters than can be altered).
GSEA and Globaltest are applied to this semi-synthetic dataset, and the whole process is repeated N times (250 by default).
Our approach is not fully synthetic because we make use of the original values of the data matrices, however we do use random numbers to enrich the selected compounds with the decoy signal.


The following code runs 250 simulations using GSEA and Globaltest on a dataset enriched with a signal magnitude of 1 and standard deviation of 0.5. 

```bash 
Rscript misIDreplicates.R --globaltest --gsea \
  --niter 250 --npaths 1 --seed 0  \
  --signal_mode "add_distrib" --signal_mean 1 --signal_sd 0.5 \
  --outputDir "Robjs/semisynth_"  --output "sims_gt_gsea" 
```

#### 2.1 Metropolis coefficient

This study uses real pathways, instead of chosing a random set of metabolites, to select the variables to enrich.
However, it is often the case that only some compounds of a pathway are present in the dataset.
Infering pathway activity from a small number of components can be challenging, thus our is based on two desired charactersitics:

1. Pathways should have ideally as many compounds in the dataset as possible
2. The compounds in the pathway should ideally not be part of other pathways, to avoid secondary pathways from being enriched as a side effect


## Datasets

|Dataset|# Case/Control|# Variables|
|:---:|:---:|:---:|
|[Stevens _et al._ (2018)](https://doi.org/10.1007/s11306-018-1393-1)| 332/667 | 352 |
|[Su _et al._ (2020)](https://doi.org/10.1016/j.cell.2020.10.037)| 254/133 | 223 |


