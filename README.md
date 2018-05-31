# ComTreeOpt

The R package **ComTreeOpt** is a 

## Installation


You can install the latest version with the commands:
```r
install.packages("devtools")
devtools::install_github("username/packagename")
```


## Usage

```r
library(ComTreeOpt)
?ComTreeOpt
data(FSN_species)
tree <- ComTreeOpt(family, genus, species)
plot(tree, show.node.label=T)
```

## References
Gastauer, M., & Meira-Neto, J. A. A. (2013a). Avoiding inaccuracies in tree calibration and phylogenetic community analysis using Phylocom 4.2. Ecological Informatics, 15, 85–90.

Gastauer, M., & Meira-Neto, J. A. A. (2013b). Interactions, Environmental Sorting and Chance: Phylostructure of a Tropical Forest Assembly. Folia Geobotanica, 1–17.

Gastauer, M., & Meira Neto, J. A. A. (2017). Updated angiosperm family tree for analyzing phylogenetic diversity and community structure . Acta Botanica Brasilica . scielo

## License

This package is free and open source software, licensed under GPL-2.

