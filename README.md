# Evaluating Gene Set Analysis methods on metabolomics data: GSEA and Globaltest

This repository contains scripts and ulitity functions to perform an evaluation of GSEA and Globaltest for metabolomics data.

## Usage 

### 0. Data pre-processing, cleaning, normalisation, imputation

The pre-processing pipeline is fully implemented in `scripts/pre-processing.R`, the details of which can be found in the Methods section of the thesis.
The raw data (as provided in the sources, i.e. supplementary information of publications or data repositories like MetaboLights) are provided under `./resources/`, and are the starting point of the pipeline. The steps include:
  
  - Removal of variables with > 50% NAs
- SVD-based imputation
- Log-transform
- Centering, Pareto scaling

The output files are serialized `.rds` objects saved under `./Robjs/` for easy import into later scripts.

> **CAUTION:** Parts of this script involve the download of a large scRNA-Seq dataset that is later reduced to pseudobulk data. This can take time and be computationally intensive, which is why all `.rds` are already provided under `./Robjs/`

```bash 
Rscript scripts/pre-processing.R 
```


### 1. Metabolite misidentification

To assess the effects of metabolite misidentification on each method, a simulation study was performed in which compounds were replaced by compounds with a similar mass (within 20ppm of the original mass). For example, the KEGG compound C00526 (Deoxyuridine) with mass 228.0746 can be replaced by C09317, with mass 228.0786.

```bash 
Rscript misIDreplicates.R --globaltest --gsea --niter 250 --outputDir "./Robjs/misID" --seed 42
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
Rscript simulations_fgsea_gt.R --globaltest --gsea \
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
  |[Su _et al._ (2020)](https://doi.org/10.1016/j.cell.2020.10.037) metabolomics| 254/133 | 223 |
  |[Su _et al._ (2020)](https://doi.org/10.1016/j.cell.2020.10.037) transcriptomics| 254/16 | 10,988 |
  
  
  
