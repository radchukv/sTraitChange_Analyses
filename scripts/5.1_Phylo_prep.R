## A script to prepare vertebrate mega-tree
## preparing 100 vertebrate phylogenies

library(ape)
library(tidyverse)
library(magrittr)
library(RRphylo)
library(rtrees)
library(geiger)

# set seed for reproducibility of the results
set.seed(100)
# 1. Prepare trees for each class -----------------------------------------

# read in the data
# 1. read in the full dataset (to be able to get unique species names)
dat <- readRDS(file = './output_forSEM_temp/all_SEM.RDS')  # check whether this is the right directory to use....
sp_only <- dat %>%
  mutate(Sp1 = case_when(
    Species == 'Cyanistes caeruleus' ~ 'Parus caeruleus',
    Species == 'Thalasseus sandvicensis' ~ 'Sterna sandvicensis',
    Species == 'Setophaga caerulescens' ~ 'Dendroica caerulescens',
    Species == 'Thalassarche melanophris' ~ 'Thalassarche melanophrys',
    Species == 'Ichthyaetus audouinii' ~ 'Larus audouinii',
    Species == 'Stercorarius maccormicki' ~ 'Catharacta maccormicki',
    TRUE ~ Species)) %>%
  distinct(Sp1, .keep_all = TRUE)


sp_only$Species_und <- unlist(lapply(1:nrow(sp_only), FUN = function(x){
  binary <- strsplit(as.character(sp_only$Sp1[x]), " ")
  Underscore <- paste(binary[[1]][1], binary[[1]][2], sep = "_")
}))


# generating class-specific trees using rtrees
# 1. birds
Bird_multi <- get_tree(sp_only$Species_und[sp_only$Taxon %in% c('Bird')],
                       taxon = 'bird', show_grafted = TRUE)  ## this is a multiphylo object
bird_1 <- Bird_multi[[1]]
plot(bird_1, label.offset = 0.3, cex = 0.6)


# 2. mammals
Mammal_multi <- get_tree(sp_only$Species_und[sp_only$Taxon %in% c('Mammal')],
                         taxon = 'mammal',
                         mammal_tree = "vertlife", show_grafted = TRUE)

mamm_1 <- Mammal_multi[[1]]
plot(mamm_1)

# 3. reptiles
Rept_multi <- get_tree(sp_only$Species_und[sp_only$Taxon %in% c('Reptilia')],
                       taxon = 'reptile',
                       show_grafted = TRUE)
rept_1 <- Rept_multi[[1]]
plot(rept_1)


# 4. fish
Fish <- get_tree(sp_only$Species_und[sp_only$Taxon %in% c('Fish')],
                       taxon = 'fish', fish_tree = 'timetree',
                       show_grafted = TRUE)
plot(Fish)

# 5. turtles
turtles_Thomson <- read.tree('./data/phylogenies/Thomson_TurtlePhylogeny_sampleFromPosterior.100OnlyT.tre')  # 100 random trees
class(turtles_Thomson[[1]])
turtles_Thomson[[1]]$tip.label

# I will have to read in the labels separately and replace them here
spLab <- read.csv('./data/phylogenies/turtle_spLabels.csv')
spLab$Species_und <- unlist(lapply(1:nrow(spLab), FUN = function(x){
  text <- strsplit(spLab$Species[x], '_')
  Underscore <- paste0(text[[1]][1], '_', text[[1]][2])
}))

# and now extract the two turtles that I have in the dataset and still then need to add them to the order tree
# i.e. the time of convergence will have to be checked and added to the backbone tree
two_turtl <- keep.tip(turtles_Thomson, as.character(spLab$Code[spLab$Species_und %in% c("Chrysemys_picta", "Chelonia_mydas")]))
two_turtl

rename.tips.phylo <- function(tree, names) {
  tree$tip.label <- names
  return(tree)
}

trtl_trees_renamed <- lapply(two_turtl, rename.tips.phylo,
                               names = c(spLab$Species_und[spLab$Code == 50],
                               spLab$Species_und[spLab$Code == 32]))
class(trtl_trees_renamed) <- "multiPhylo"

turtles_1 <- trtl_trees_renamed[[1]]

# get all tress into multiPhylo object
tree_list <- list(fish = Fish, squam = rept_1,
                  turtl = turtles_1,
                  birds = bird_1, mamm = mamm_1)
class(tree_list) <- "multiPhylo"


# 2. prepare the class-level tree -----------------------------------------
tip.labels <- c("fish", "squam", "turtl", "birds", "mamm") # same as the names in the multiPhylo object

## Make a tree with just classes:
edge <- matrix(c(7, 5,
                 9, 4,
                 9, 3,
                 8, 9,
                 8, 2,
                 7, 8,
                 6, 7,
                 6, 1), byrow=TRUE, ncol=2)
## Dates from Timetree of life (timetree.org)
edge.length <- c(319, 261, 261, 280-261, 280, 319-280, 429-319, 429)
Nnode <- 4
ordertree <- list(edge=edge, Nnode=Nnode, tip.label=tip.labels, edge.length=edge.length)
class(ordertree) <- 'phylo'


plot(ordertree)
tiplabels()
nodelabels()
edgelabels(ordertree$edge.length, bg="white", col="black", font=0.8)


# 3. put the class-level trees on the order-tree --------------------------

otax <- data.frame("Class" = ordertree$tip.label, "Superclass"=c("", rep("Tetrapoda", 4)))
rownames(otax) <- ordertree$tip.label
classtree <- nodelabel.phylo(ordertree, otax, ncores=1)

vert_tree <- glomogram.phylo(classtree, tree_list)

plot(vert_tree, cex = 0.8)
vert_tree # 72 species (reindeer 2 subspecies collapsed into 1

# save this tree
ape::write.tree(vert_tree, file = "./data/phylogenies_100/vert1.tre")


# write 100 trees created by simply taking one tree from each multi-tree, whenever those are available (I do not
# Draw them randomly, simply take in order)
for(i in 2:100){
  # pick respective tree from the multi-trees
  bird <- Bird_multi[[i]]
  mamm <- Mammal_multi[[i]]
  rept <- Rept_multi[[i]]
  trtl <- trtl_trees_renamed[[i]]
  tree_list <- list(fish = Fish, squam = rept,
                    turtl = trtl, birds = bird, mamm = mamm)
  class(tree_list) <- "multiPhylo"

  vert_tree <- glomogram.phylo(classtree, tree_list)
  ape::write.tree(vert_tree, file = paste0("./data/phylogenies_100/vert", i, ".tre"))
}
