# ComTreeOpt

The R package **ComTreeOpt** is a Community Tree Optimizer using Open Tree of Life. It enhances the phylogenetic resolution of community trees for the analysis of phylogenetic diversity and dispersion measures.

## Installation
You can install the latest version with the commands:
```r
install.packages("devtools")
devtools::install_github("iantrotter/ComTreeOpt")
```


## Usage

```r
library(ComTreeOpt)
?ComTreeOpt
data(FSN_species)
tree <- ComTreeOpt(FSN_species)
plot(tree, show.node.label=T)

data(EifelDataset_species)
tree <- ComTreeOpt(EifelDataset_species)
plot(tree, show.node.label=T)
```

## Citation
Please cite the package as:
Gastauer, Markus, Cecílio Frois Caldeira, Ian Trotter, Silvio Junio Ramos, and João Augusto Alves Meira Neto (2018). Optimizing Community Trees Using the Open Tree of Life Increases the Reliability of Phylogenetic Diversity and Dispersion Indices. Ecological Informatics 46, 192–98. https://doi.org/10.1016/j.ecoinf.2018.06.008.

## References
Gastauer, M., & Meira-Neto, J. A. A. (2013a). Avoiding inaccuracies in tree calibration and phylogenetic community analysis using Phylocom 4.2. Ecological Informatics, 15, 85–90.

Gastauer, M., & Meira-Neto, J. A. A. (2013b). Interactions, Environmental Sorting and Chance: Phylostructure of a Tropical Forest Assembly. Folia Geobotanica, 1–17.

Gastauer, M., & Meira Neto, J. A. A. (2017). Updated angiosperm family tree for analyzing phylogenetic diversity and community structure . Acta Botanica Brasilica . scielo

## License

This package is free and open source software, licensed under GPL-2.

