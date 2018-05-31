#' @title Enhance user-specific community trees
#'
#' @description Based on an interface with Open Tree of Life, the Community Tree Optimizer (ComTreeOpt) classifies families from user-specific lists into resolvable and non-resolvable families, builds resolved subtrees for resolvable families and inserts them in the megatree R20160415.new (Gastauer & Meira-Neto, 2017) previously pruned to species from non-resolvable families.
#'
#' @param family A list of families.
#' @param genus A list of genuses.
#' @param species A list of species. The list of families, genus and species must have the same number of entries.
#' @param megatree.uri URI to a megatree. The default value is the R20160415.new megatree (Gastauer & Meira-Neto, 2017).
#' @return A phylogenetic tree containing the user-specified species, based on the Open Tree of Life.
#' @examples
#' data(FSN_species);
#' tree <- ComTreeOpt(family, genus, species);
#' plot(tree, show.node.label=T)
#' @references Gastauer, M. &  Meira-Neto, J. A. A. (2017) Updated angiosperm family tree for analyzing phylogenetic diversity and community structure. \emph{Acta Botanica Brasília}, vol.31, n.2, pp.191-198. DOI: \url{http://dx.doi.org/10.1590/0102-33062016abb0306}.
ComTreeOpt <- function(family, genus, species, megatree.uri = "https://www.dropbox.com/s/rc8jtmb6b2hdweo/R20160415.new?dl=0") {
  #### Checking the input to the function ####
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
  max.q <- 2500;

  # Iterate through the species and query the Open Tree of Life, using
  # the tnrs_match_names function from the  rotl package
  n <- 0
  df.data <- data.frame();
  while (n < length(species)) {
    nn <- min(n+max.q,length(species))
    df.data <-rbind(df.data, tnrs_match_names(names = as.vector(species[n+1:nn])));
    n <- nn;
  }
  rm(max.q, n, nn); # Clean up temporary variables

  # Create a data frame with all the variables
  df.data <- cbind(data.frame(family=family, genus=genus, species=species, stringsAsFactors=F), df.data);
  rm(family, genus, species); # Clean up workspace

  #### Divide the results into a good part (clean) and a bad part (dirty) ####
  #err.ind <- (is.na(df.data$approximate_match) # Is missing data in the approximate_match field
  #          | df.data$approximate_match==T     # approximate_match was returned as TRUE
  #          | is.na(df.data$ott_id)            # Is missing the ott_id
  #          | df.data$ott_id==T                # ott_id was returned as TRUE
  #          | is.na(df.data$is_synonym)        # Is missing data to tell whether it is a synonym
  #          | df.data$is_synonym==T            # The species name is a synonym
  #          | df.data$flags=="INCERTAE_SEDIS_INHERITED"); # Inherits from incertae cedis

  # All species from all families that tnrs_match is unable to match are
  # marked with TRUE in a column name selSp
  err.ind1 <- (df.data$approximate_match==TRUE
               | is.na(df.data$approximate_match)
               | df.data$flags=="INCERTAE_SEDIS_INHERITED");
  err.families = unique(df.data$family[err.ind1]);
  df.data$selSp <- F;
  df.data$selSp[df.data$family %in% err.families] <- T;
  rm(err.ind1, err.families);

  # Families with only one or two species should also be marked
  selNb <- data.frame(t(t(table(df.data$family))))
  df.data$selNb <- selNb$Freq[match(df.data$family, selNb$Var1)]
  df.data$selSp <- ifelse(df.data$selNb<=2, T, df.data$selSp)
  rm(selNb);

  # Separate between dirty and clean data
  df.clean <- df.data[!df.data$selSp,names(df.data)[1:10]];
  df.dirty <- df.data[df.data$selSp,names(df.data)[1:10]];

  # Create a species File, including species from df.dirty and families from df.clean
  dirty.species <- paste(df.dirty$family, df.dirty$genus, gsub(" ", "_", df.dirty$species), sep="/");
  for (ds in dirty.species) {
    warning(paste("Unable to resolve species:", ds));
  }
  rm(ds);
  clean.families <- unique(df.clean$family);
  species <- c(dirty.species, clean.families)

  # Create tree
  tree <- phylomatic(species, get="POST", taxnames=F, treeuri=megatree.uri);
  rm(species, dirty.species);

  # Loop through families and paste onto the tree
  for (f in clean.families) {
    fam.ind <- f==df.clean$family;
    fam.tree <- tol_induced_subtree(ott_ids = df.clean$ott_id[fam.ind], label_format = "id");
    fam.tree$tip.label <- gsub("ott", "", fam.tree$tip.label);
    fam.tree$tip.label <- as.character(df.clean$species[match(fam.tree$tip.label, df.clean[["ott_id"]])])
    nodes <- matrix(NA,nrow=fam.tree$Nnode,ncol=1);
    df.nodes<-data.frame(nodes);
    df.nodes[1,1] <- f;
    fam.tree$node.label <- df.nodes$nodes;
    fam.tree$root.edge <- 0;
    tree$tip.label <- gsub(f, "NA", tree$tip.label);
    tree <- paste.tree(tree, fam.tree);
    tree$tip.label <- gsub("NA", f, tree$tip.label);
  }
  rm(f, fam.ind, fam.tree, nodes, df.nodes);

  # Return the tree
  return(tree);

}
