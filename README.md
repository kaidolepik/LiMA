# (I-)LiMA: robust inference of molecular mediation from summary statistics

This repository contains the source code for **Integrated Likelihood based Mediation Analysis (I-LiMA)**, along with:
* Reproducible workflows and scripts for comparing I-LiMA with alternative mediation analysis approaches, including a Mendelian randomization (MR)-based baseline (univariable MR to estimate total effects, multivariable MR to estimate direct effects) that is referred to as MR+MVMR
* Code for reproducing all **figures** and **tables** presented in the associated paper

## Citation

If you use the code or workflows in this repository, please consider citing the associated article:

Lepik K, Auwerx C, Sadler MC, van der Graaf A, Ojavee SE, Kutalik Z. (I-)LiMA: robust inference of molecular mediation from summary statistics linking classical epidemiological risk factors to cardiovascular outcomes.

## Usage

This repository uses [Nextflow](https://www.nextflow.io/) (version 23.04) to run workflows. To facilitate environment setup, all additional dependencies can be managed via [Docker](https://www.docker.com/).

Clone the respository:

```
git clone https://github.com/kaidolepik/LiMA.git
```

### Running simulations

The workflow relies on the following external data sets not included in the repository. Please refer to the **Methods** and **Data Availability** sections of the paper for details on data access:
* [cis-eQTL summary data](https://eqtlgen.org/cis-eqtls.html) from the eQTLGen Consortium
* Individual-level gene expression data (RPKM) from the [CoLaus](https://www.colaus-psycolaus.ch/professionals/colaus) cohort

With the data in place, simulations can be run as follows, where the `--type` argument refers to a configuration file in `simulations/data/raw/conf/` that defines the simulation parameters:

```
cd LiMA/simulations
docker build -t simulations .
nextflow run simulations.nf --type=test -profile local_docker
```

### Running real data-based analyses

The workflow relies on the following external data sets not included in the repository. Please refer to the **Methods** and **Data Availability** sections of the paper for details on data access:
* UK Biobank GWAS summary statistics from the [Neale Lab](https://www.nealelab.is/uk-biobank)
* Genome reference data from [UK10K](https://www.uk10k.org/data_access.html)
* Metabolite data from the [Lotta et al. 2021](https://www.nature.com/articles/s41588-020-00751-5) study
* Protein data from the [INTERVAL study](https://www.nature.com/articles/s41586-018-0175-2)

With the data in place (with UKBB traits configured in `real_data/data/raw/conf/`), the mediation analysis can be run as follows, where the `--mediators` argument specifies which mediator data set to use:

```
cd LiMA/real_data
docker build -t real_data .
nextflow run mediation_application.nf --mediators=Lotta_et_al_2021 -profile local_docker
```

### Reproducing figures and tables

All figures and tables (excluding graphical abstracts) are generated directly from the output of the Nextflow workflows. Refer to the following scripts:
* `simulations/scripts/figure_prep.R` for simulation-based figures and tables
* `real_data/scripts/figure_real_prep.R` for real data-based analyses
