# Multi-omics EEN Analysis Presented in [Häcker and Siebert et al.](https://www.medrxiv.org/content/10.1101/2023.12.21.23300351v1)

## Prerequisites

The code in this repository was tested with R version 4.2.1 and python version
3.8.11.

For required `R` packages and the used versions please see `seesionInfo.txt`

For required `python` packages and the used versions please see
`requirements.txt`

## Feature Selection, Community Extraction and Mode Evaluation
The main file of the analysis is `multi_omics.R`, essentially it's just a
wrapper for finding the optimal number of features to select, performing
evaluations and subsequent correlation network analysis.

The file uses fixed metabolome data and works with variable microbiome data that
has to be passed via a command line interface together with a result path.

For example
```shell
Rscript multi_omics.R Data/16S_data.RData Results
```

The microbiome data should be an RData file containing the processed zOTU data as
well as the taxonomic annotation. Since result files produced always have the
same name, different analyses should be saved to different folders.

Other options to pass regard the parameter to use for selecting the optimal
number of metabolomic and microbial features. To see them invoke

```shell
Rscript multi_omics.R --help
```

To recreate all experiments we provide a bash script including all parameters
(and optionally a slurm-adapted version). It can be run (linux and macOS only)
by

```shell
bash run_experiments.sh && bash run_downstream_analysis.sh
```

## Network Plotting

For plotting the network, a combination of `R` and `python` are used.
The required scripts are wrapped in `run_network_plotting.sh` and can be run via

```shell
bash run_newtork_plotting.sh
```

For parameter information please use

```shell
Rscript community_network_plotting/community_network_data.R --help
```

and

```shell
python community_network_plotting/community_network_plotting.py --help
```
