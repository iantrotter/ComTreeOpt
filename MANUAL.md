
# Community Tree Optimizer (ComTreeOpt) Manual

Dear user,

In this manual you find comands and comments about the function "ComTreeOpt" to facilitate its application. This is a preliminary version, please test it and let us know any bugs, problems or general comments. 
For that, please get in contact per Email (markus.gastauer@itv.org, mgasti@hotmail.com), per phone (+55 91 3213 5554) or per skype (mgasti).

Best wishes,

Markus Gastauer


# Instalation of the package 'ComTreeOpt'

Currently, the package 'ComTreeOpt' is hosted at the Ian Trotter's github: 
https://github.com/iantrotter/ComTreeOpt

For package installation from github, you need the package 'devtools':

```r
install.packages("devtools") # if 'devtools' is already installed, you may pull this step
library("devtools") 
install_github("iantrotter/ComTreeOpt") # don't scare, this command will install many packages on your maschine
library("ComTreeOpt") 
```

On some machines, the installation or application of 'ComTreeOpt' requires additional packages. Please install them manually and let us know them, so that we can update this in the original code. Thanks!

For some steps in the last section of this script, you need the 'picante' package as well:

```r
library(picante)
```

# Running example data

We provide two example datasets in the package. Running these examples is not very time-consuming, but assures that all components are installed correctly.

```r
data(FSN_species)
tree <- ComTreeOpt(FSN_species)
plot(tree, show.node.label=T)

data(EifelDataset_species)
tree <- ComTreeOpt(EifelDataset_species)
plot(tree, show.node.label=T)

# both examples show a lot of warnings, that you can visualize using 
warnings()
```

Herein, all unresolved species are listed. Species remain unresolved due to different reasons:
1. Incomplete identification: Morpho-Species identified to genus or family level only.
2. Invalid species names or valid species names not available in the Open Tree of Life.
3. Species from families that contain only one or two species in the list (i.e., no polytomies to be resolved in descending branches).

If we look at the Eifel dataset, you can see that Betonica officinalis is currently unresolvable. To enhance resolution of your community tree, you may investigate, why B. officinalis is not found in the Time Tree of Life database or look for synonyms of the name. Furthermore, a lot of species are from families that have only one or two species in the dataset, but there is nothing we can do to enhance phylogenetic resolution.


# Running user-specific lists

Before you upload your user-specific file, let's have a look at the example 'species' files:

```r
class(EifelDataset_species) # it is a dataframe
str(EifelDataset_species) # with a single column called 'V1', 
	# and all entries in this column are in the format that you already know from the phylocom package:
	# family/genus/species_name

# if your data are in the same format, you may upload them after setting your working directory

species <- read.csv("species.txt", header=F) # don't forget to change to the correct file name!
class(species) # ok, it is a dataframe!
str(species) # stop, something is different! It is not the same format as the example file, so we run the 
# following command:
species$V1 <- as.character(species$V1) 
str(species) # Now it is the same format as 'EifelDataset_species', guaranteeing that the following comand will run:
tree <- ComTreeOpt(species) 
```

Your community tree is optimized as suggested by Open Tree of Life!


# Working with the optimized community tree

Following some examples how to plot, export or use the optimized community tree to compute phylogenetic diversity or phylogenetic community structure; for the later, please keep in mind that ultrametric, calibrated trees are required. 

Plotting:
```r
plot(tree, show.node.label=T)
```

Exporting the optimized community tree as newick file to your working directory:
```r
write.tree(tree, "tree.txt")
```

In your phylocom cmd window, you may now calibrate the exported tree with the following command (when an 'ages' is provided), please copy the following code:

```r
phylocom bladj -f tree.txt > tree_dt.txt
```

Alternatively, you may calibrate your tree in R using 'brranching' package.

If you decide to continue the computation of phylogenetic diversity and community structure in your cmd window, we will end here.

But you may re-import the optimized, calibrated community tree from your working directory to R and run your analysis here:

```r
tree_dt <- read.tree("tree_dt.txt")
tree_dt
plot(tree_dt)

# import 'sample' file
sample <- readsample("sample.txt")

# compute NRI
phy <- cophenetic(tree_dt)
NRI <- ses.mpd(sample, phy, null.model = "phylogeny.pool")
```
