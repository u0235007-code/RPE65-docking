# writes ChimeraX script to superimpose molecular structure of lutein molecules pairwise
setwd("~/Docking/dietary_lutein")
# read libraries
library(readr)
library(tidyr)
library(tidyverse)
library(tibble)
library(ggplot2)
library(geometry)
# read the source filenames
lutein.ligand <- read.csv("root_names_luteinligand.csv", stringsAsFactors = FALSE, header = FALSE, col.names = "lutein")
# build the ChimeraX script to superimpose beta from lutein to replace epsilon of dietary lutein
ChimeraX.script <- tibble("command" = character())
indexX <- as.integer(1)
ChimeraX.script[indexX, "command"] <- "close"
indexX <- 1 + indexX
ChimeraX.script[indexX, "command"] <- "cd ~/Docking/dietary_lutein"
indexX <- 1 + indexX
ChimeraX.script[indexX, "command"] <- "open *_luteinligand.pdb"
indexX <- 1 + indexX
ChimeraX.script[indexX, "command"] <- str_c("align #", 1, " toAtoms #", 1, " move residues")
indexX <- 1 + indexX
#
#build list of lutein names
lutein.pair <- tibble("lutein1" = character(), "lutein2" = character(), "rmsd" = character(), "identical" = integer())
indexP <- as.integer(1)
indexM <- as.integer(1)
set.seed(1103)
while (indexM <= nrow(lutein.ligand)) {
  lutein.pair[indexP:(indexP+nrow(lutein.ligand)), "lutein1"] <- lutein.ligand[indexM, "lutein"]
  indexN <- as.integer(1)
  while (indexN <= nrow(lutein.ligand)) {
    lutein.pair[indexP, "lutein2"] <- lutein.ligand[indexN, "lutein"]
    ChimeraX.script[indexX, "command"] <- str_c("align #", indexM, " toAtoms #", indexN, " move residues")
    indexX <- 1 + indexX
    indexN <- 1 + indexN
    indexP <- 1 + indexP
  } #for each lutein molecule
  indexM <- 1 + indexM
} # for each lutein molecule compare to each lutein
write.table(ChimeraX.script, "compare.luteins.cxc", row.names = FALSE, quote = FALSE, col.names = FALSE)
#
# read the rmsd values and join with the lutein pairs
lutein.rmsd.output <- read.csv("compare.luteins.rmsd.csv", stringsAsFactors = FALSE, header = FALSE)
lutein.rmsd.position <- str_locate(lutein.rmsd.output[ , ], "angstrom")
lutein.rmsd.position[, 1] <- lutein.rmsd.position[, 1] - 6
lutein.rmsd.position[, 2] <- lutein.rmsd.position[, 2] - str_length("angstrom") -  1
lutein.pair[1:nrow(lutein.rmsd.position), "rmsd"]<- str_sub(lutein.rmsd.output[ , ], lutein.rmsd.position[, 1:2])
# treat data as vectors to avoid while loop; must match vector length precisely
indexM <- as.integer(1)
while (indexM <= nrow(lutein.rmsd.position)) {
  lutein.pair[indexM, "identical"] <- str_count(lutein.pair[indexM, "lutein1"], as.character(lutein.pair[indexM, "lutein2"]))
  indexM <- 1 + indexM
} # for each lutein1 compare to lutein2 and check if identical
write.table(lutein.pair[1:nrow(lutein.rmsd.position), ], "lutein.pair.csv", sep = ", ", row.names = FALSE, quote = FALSE)
#
# read the source filenames for zeaxanthin ligands
zeaxanthin.ligand <- read.csv("root_names_zeaxanthinligand.csv", stringsAsFactors = FALSE, header = FALSE, col.names = "zeaxanthin")
#build list of zeaxanthin pairs
zeaxanthin.pair <- tibble("zeaxanthin1" = character(), "zeaxanthin2" = character(), "rmsd" = character(), "identical" = integer())
indexP <- as.integer(1)
indexM <- as.integer(1)
set.seed(1103)
while (indexM <= nrow(zeaxanthin.ligand)) {
  zeaxanthin.pair[indexP:(indexP+nrow(zeaxanthin.ligand)), "zeaxanthin1"] <- zeaxanthin.ligand[indexM, "zeaxanthin"]
  indexN <- as.integer(1)
  while (indexN <= nrow(zeaxanthin.ligand)) {
    zeaxanthin.pair[indexP, "zeaxanthin2"] <- zeaxanthin.ligand[indexN, "zeaxanthin"]
    indexN <- 1 + indexN
    indexP <- 1 + indexP
  } #for each zeaxanthin molecule
  indexM <- 1 + indexM
} # for each zeaxanthin molecule compare to each zeaxanthin
# read the rmsd values and join with the lutein pairs
zeaxanthin.rmsd.output <- read.csv("compare.zeaxanthins.rmsd.csv", stringsAsFactors = FALSE, header = FALSE)
zeaxanthin.rmsd.position <- str_locate(zeaxanthin.rmsd.output[ , ], "angstrom")
zeaxanthin.rmsd.position[, 1] <- zeaxanthin.rmsd.position[, 1] - 6
zeaxanthin.rmsd.position[, 2] <- zeaxanthin.rmsd.position[, 2] - str_length("angstrom") -  1
zeaxanthin.pair[1:nrow(zeaxanthin.rmsd.position), "rmsd"]<- str_sub(zeaxanthin.rmsd.output[ , ], zeaxanthin.rmsd.position[, 1:2])
# treat data as vectors to avoid while loop; must match vector length precisely
indexM <- as.integer(1)
while (indexM <= nrow(zeaxanthin.rmsd.position)) {
  zeaxanthin.pair[indexM, "identical"] <- str_count(zeaxanthin.pair[indexM, "zeaxanthin1"], as.character(zeaxanthin.pair[indexM, "zeaxanthin2"]))
  indexM <- 1 + indexM
} # for each lutein1 compare to lutein2 and check if identical
write.table(zeaxanthin.pair[1:nrow(zeaxanthin.rmsd.position), ], "zeaxanthin.pair.csv", sep = ", ", row.names = FALSE, quote = FALSE)
#
# rebuild ChimeraX script for 121 molecules in place of 122
ChimeraX.script <- tibble("command" = character())
indexX <- as.integer(1)
ChimeraX.script[indexX, "command"] <- "close"
indexX <- 1 + indexX
ChimeraX.script[indexX, "command"] <- "cd ~/Docking/dietary_lutein"
indexX <- 1 + indexX
ChimeraX.script[indexX, "command"] <- "open *_mesozeaxanthinligand.pdb"
indexX <- 1 + indexX
ChimeraX.script[indexX, "command"] <- str_c("align #", 1, " toAtoms #", 1, " move residues")
indexX <- 1 + indexX
#
# read the source filenames for meso-zeaxanthin ligands
zeaxanthin.ligand <- read.csv("root_names_mesozeaxanthinligand.csv", stringsAsFactors = FALSE, header = FALSE, col.names = "zeaxanthin")
#build list of meso-zeaxanthin pairs
zeaxanthin.pair <- tibble("zeaxanthin1" = character(), "zeaxanthin2" = character(), "rmsd" = character(), "identical" = integer())
indexP <- as.integer(1)
indexM <- as.integer(1)
set.seed(1103)
while (indexM <= nrow(zeaxanthin.ligand)) {
  zeaxanthin.pair[indexP:(indexP+nrow(zeaxanthin.ligand)), "zeaxanthin1"] <- zeaxanthin.ligand[indexM, "zeaxanthin"]
  indexN <- as.integer(1)
  while (indexN <= nrow(zeaxanthin.ligand)) {
    zeaxanthin.pair[indexP, "zeaxanthin2"] <- zeaxanthin.ligand[indexN, "zeaxanthin"]
    ChimeraX.script[indexX, "command"] <- str_c("align #", indexM, " toAtoms #", indexN, " move residues")
    indexX <- 1 + indexX
    indexN <- 1 + indexN
    indexP <- 1 + indexP
  } #for each meso-zeaxanthin molecule
  indexM <- 1 + indexM
} # for each meso-zeaxanthin molecule compare to each meso-zeaxanthin
write.table(ChimeraX.script, "compare.meso.zeaxanthins.cxc", row.names = FALSE, quote = FALSE, col.names = FALSE)
#
# read the rmsd values and join with the meso-zeaxanthin pairs
zeaxanthin.rmsd.output <- read.csv("compare.meso.zeaxanthin.rmsd.csv", stringsAsFactors = FALSE, header = FALSE)
zeaxanthin.rmsd.position <- str_locate(zeaxanthin.rmsd.output[ , ], "angstrom")
zeaxanthin.rmsd.position[, 1] <- zeaxanthin.rmsd.position[, 1] - 6
zeaxanthin.rmsd.position[, 2] <- zeaxanthin.rmsd.position[, 2] - str_length("angstrom") -  1
zeaxanthin.pair[1:nrow(zeaxanthin.rmsd.position), "rmsd"]<- str_sub(zeaxanthin.rmsd.output[ , ], zeaxanthin.rmsd.position[, 1:2])
# treat data as vectors to avoid while loop; must match vector length precisely
indexM <- as.integer(1)
while (indexM <= nrow(zeaxanthin.rmsd.position)) {
  zeaxanthin.pair[indexM, "identical"] <- str_count(zeaxanthin.pair[indexM, "zeaxanthin1"], as.character(zeaxanthin.pair[indexM, "zeaxanthin2"]))
  indexM <- 1 + indexM
} # for each lutein1 compare to lutein2 and check if identical
write.table(zeaxanthin.pair[1:nrow(zeaxanthin.rmsd.position), ], "meso.zeaxanthin.pair.csv", sep = ", ", row.names = FALSE, quote = FALSE)
