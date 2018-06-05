#' @title Enhance user-specific community trees
#'
#' @description Based on an interface with Open Tree of Life, the Community Tree Optimizer (ComTreeOpt) classifies families from user-specific lists into resolvable and non-resolvable families, builds resolved subtrees for resolvable families and inserts them in the megatree R20160415.new (Gastauer & Meira-Neto, 2017) previously pruned to species from non-resolvable families.
#'
#' @param input_data (1) A data.frame with family, genus and species columns specified by parameters family.colname, genus.colname and species.colname. Or: (2) The name of a file containing the species data in the format 'family/genus/species_name' (one entry per line) and no headers. Or (3) A character or factor vector with entries of the form 'family/genus/species_name'.
#' @param family.colname The name of the column that contains the family names.
#' @param genus.colname The name of the column that contains the genus names.
#' @param species.colname The name of the column that contains the species names.
#' @param megatree.uri URI to a megatree. The default is the R20160415.new megatree (Gastauer & Meira-Neto, 2017).
#' @return A phylogenetic tree containing the user-specified species, based on the Open Tree of Life.
#' @examples
#' data(FSN_species);
#' tree <- ComTreeOpt(FSN_species);
#' plot(tree, show.node.label=T)
#'
#' data(EifelDataset_species)
#' tree <- ComTreeOpt(EifelDataset_species)
#' plot(tree, show.node.label=T)
#' @references Gastauer, M. &  Meira-Neto, J. A. A. (2017) Updated angiosperm family tree for analyzing phylogenetic diversity and community structure. \emph{Acta Botanica Brasília}, vol.31, n.2, pp.191-198. DOI: \url{http://dx.doi.org/10.1590/0102-33062016abb0306}.
#ComTreeOpt <- function(input_data, family.colname = "family", genus.colname = "genus", species.colname = "species", megatree.uri = "https://www.dropbox.com/s/rc8jtmb6b2hdweo/R20160415.new?dl=0") {
ComTreeOpt <- function(input_data, family.colname = "family", genus.colname = "genus", species.colname = "species", megatree.uri = "https://raw.githubusercontent.com/iantrotter/ComTreeOpt/master/data/R20160415.new.db") {
  #### Checking the input to the function ####
  if (class(input_data) == "data.frame" && ncol(input_data) >= 3) {
    # If input_data is a data.frame, get the relevant columns, as specified by the function arguments
    message(paste("Input data is a 'data.frame': attempting to use columns '", family.colname, "' for family, '", genus.colname, "' for genus, and '", species.colname, "' for species.", sep=""));
    family <- input_data[,family.colname];
    genus  <- input_data[,genus.colname];
    species <- input_data[,species.colname];
  } else if (class(input_data) == "data.frame" && ncol(input_data) == 1) {
    # If input_data is a data.frame with a single column, it is assumed to be in
    # family/genus/species_name format
    tmp <- matrix(unlist(strsplit(input_data[,1], split="/")), ncol = 3, byrow = TRUE);
    family <-  tmp[,1];
    genus  <-  tmp[,2];
    species <- tmp[,3];
    rm(tmp);
  } else if (class(input_data) %in% c("character", "factor") && length(input_data) == 1) {
    # If input_data is a single string/factor, assume it is a filename and read the file contents
    # File format is assumed to be family/genus/species_name with no header
    message(paste("Input data is of type 'character' or 'factor', and has length == 1: interpreted as the name of a file in the format 'family/genus/species_name', with no file header.", sep=""))
    tmp <- read.table(as.character(input_data), header=F, sep="/", stringsAsFactors = F);
    family  <- as.character(tmp$V1);
    genus   <- as.character(tmp$V2);
    species <- as.character(tmp$V3);
    rm(tmp);
  } else if (class(input_data) %in% c("character", "factor") && length(input_data) > 1) {
    # If input_data is an array of several strings/factors, assume it is a list of species in a particular format
    message("Input data is of type 'character' or 'factor', and has length > 1: interpreted as an array of strings in the format 'family/genus/species_name'.")
    tmp <- matrix(unlist(strsplit(input_data, split="/")), ncol = 3, byrow = TRUE);
    family <-  tmp[,1];
    genus  <-  tmp[,2];
    species <- tmp[,3];
    rm(tmp);
  } else {
    stop("Invalid input.");
  }

  # Check that all the input is of the same length, if not then throw an error
  if (length(family) != length(genus) || length(genus) != length(species) || length(family) != length(species)) {
    stop("Input is not the same length.");
  }
  # Make sure all the input are vectors of characters
  family <- as.character(family);
  genus  <- as.character(genus);
  species <- as.character(species);

  #### Querying the taxonomy of the species and processing the reply ####
  # Substitute underscore by space in the species name, just in case
  species <- gsub("_", " ", species);

  # Maximum number of species that can be sent in a single request
  max.q <- 250;

  # Iterate through the species and query the Open Tree of Life, using
  # the tnrs_match_names function from the  rotl package
  n <- 0
  df.data <- data.frame();
  while (n < length(species)) {
    nn <- min(n+max.q,length(species));
    df.data <-rbind(df.data, rotl::tnrs_match_names(names = as.vector(species[(n+1):nn])));
    n <- nn;
  }
  rm(max.q, n, nn); # Clean up temporary variables

  # Create a data frame with all the variables
  df.data <- cbind(data.frame(family=family, genus=genus, species=species, stringsAsFactors=F), df.data);
  rm(family, genus, species); # Clean up workspace

  # All species from all families that tnrs_match is unable to match are
  # marked with TRUE in a column name selSp
  err.ind1 <- (df.data$approximate_match==TRUE
               | is.na(df.data$approximate_match)
               | df.data$flags=="INCERTAE_SEDIS_INHERITED");
  for (i in 1:nrow(df.data)) {
    if (err.ind1[i]) {
      warning(paste("Unable to resolve species:", df.data[i,"species"]));
    }
  }
  err.families = unique(df.data$family[err.ind1]);
  df.data$selSp <- F;
  df.data$selSp[df.data$family %in% err.families] <- T;
  rm(err.ind1, err.families);

  # Families with only one or two species should also be marked
  selNb <- data.frame(t(t(table(df.data$family))))
  df.data$selNb <- selNb$Freq[match(df.data$family, selNb$Var1)]
  df.data$selSp <- ifelse(df.data$selNb<=2, T, df.data$selSp)
  for (i in 1:nrow(df.data)) {
    if (df.data$selNb[i] <= 2) {
      warning(paste("Only one or two species within the determined family of species:", df.data[i,"species"]));
    }
  }

  rm(selNb);

  # Separate between dirty and clean data
  df.clean <- df.data[!df.data$selSp,names(df.data)[1:10]];
  df.dirty <- df.data[df.data$selSp,names(df.data)[1:10]];

  # Create a species File, including species from df.dirty and families from df.clean
  dirty.species <- paste(df.dirty$family, df.dirty$genus, gsub(" ", "_", df.dirty$species), sep="/");
  clean.families <- unique(df.clean$family);
  species <- c(dirty.species, clean.families)

  # Create tree
  #  megatree.uri <- paste(find.package("ComTreeOpt", lib.loc=NULL, quiet = TRUE),"/data/R20160415.new",sep="");
  tree <- brranching::phylomatic(species, get="POST", taxnames=F, treeuri=megatree.uri);
  rm(species, dirty.species);

  # Loop through families and paste onto the tree
  for (f in clean.families) {
    fam.ind <- f==df.clean$family;
    fam.tree <- rotl::tol_induced_subtree(ott_ids = df.clean$ott_id[fam.ind], label_format = "id");
    fam.tree$tip.label <- gsub("ott", "", fam.tree$tip.label);
    fam.tree$tip.label <- as.character(df.clean$species[match(fam.tree$tip.label, df.clean[["ott_id"]])])
    nodes <- matrix(NA,nrow=fam.tree$Nnode,ncol=1);
    df.nodes<-data.frame(nodes);
    df.nodes[1,1] <- f;
    fam.tree$node.label <- df.nodes$nodes;
    fam.tree$root.edge <- 0;
    tree$tip.label <- gsub(f, "NA", tree$tip.label);
    tree <- phytools::paste.tree(tree, fam.tree);
    tree$tip.label <- gsub("NA", f, tree$tip.label);
  }
  rm(f, fam.ind, fam.tree, nodes, df.nodes);

  # Return the tree
  return(tree);

}
