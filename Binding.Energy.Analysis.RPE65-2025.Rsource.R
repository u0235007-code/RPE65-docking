# R script for extracting energy terms from vina --score_only log files
#
# define working directory, here and also in ChimeraX script
work.dir <- "~/research/docking/dietary_lutein_2025/RPE65_dietary_lutein"
ligand.dir <- "~/research/docking/AnnabelLee_2024/CHPC-2024/LIGANDS"
gromacs.dir <- "~/research/gromacs/RPE65_dietary_lutein"
ligand.code <- c("emixustat", "STQ", "VFQ", "8UH", "BCR")
carotenoid.code <- c("STQ", "VFQ", "8UH", "BCR")
xanthophyll.code <- c("STQ", "VFQ", "8UH")
ligand.name <- c("retinoid", "lutein", "meso-zeaxanthin", "zeaxanthin", "beta-carotene")
ligand.mode <- c(2, 5)
#ligand.name.csv.file <- c("stq2_names_pdbqt.csv", "stq5_names_pdbqt.csv", 
#                        "vfq2_names_pdbqt.csv", "vfq5_names_pdbqt.csv", 
#                        "zea2_names_pdbqt.csv", "zea5_names_pdbqt.csv")
ligand.name.csv.file <- c("emixustat_names_pdbqt.csv",
                          "stq2_names_pdbqt.csv", 
                          "vfq2_names_pdbqt.csv", 
                          "zea2_names_pdbqt.csv",
                          "bcr2_names_pdbqt.csv")
setwd(work.dir)
# clear the R environment use this to clear the environment
# rm(list = ls())
# read libraries, only necessary at start of new Rstudio session, resolve conflicts by reading libraries a second time
library(readr)
library(tidyr)
library(tidyverse)
library(tibble)
library(ggplot2)
library(geometry)
#
receptor <- "RPE65"
input.directory <- str_c("INPUT", receptor, sep = "_")
#
# find the landmark for bottom of chunnel
file.name.receptor.pdbqt <- str_c(input.directory, "/cow_7l0e_chainA_rigid.pdbqt")
receptor.pdbqt <- read.delim(file.name.receptor.pdbqt, header = FALSE, sep = "", stringsAsFactors = FALSE)
receptor.chunnel.landmark <- numeric(3)
receptor.chunnel.landmark[1] <- receptor.pdbqt[5007, "V7"]
receptor.chunnel.landmark[2] <- receptor.pdbqt[5007, "V8"]
receptor.chunnel.landmark[3] <- receptor.pdbqt[5007, "V9"]
receptor.pdbqt[5007,]
rm(receptor.pdbqt)
#
#
# function for cross product
vector.cross <- function(a, b) {
  if ((length(a) == 3) & (length(b) == 3)) {
    i1 <- c(2, 3, 1)
    i2 <- c(3, 1, 2)
    return(a[i1] * b[i2] - a[i2] * b[i1])
  } #if 3D vectors
  else
    stop("Warning: vector.cross passed non 3D vectors")
} #return cross-product of two 3D vectors
#
# function for vector normalization
normalize.vector <- function(a) {
  if (abs(sum(a * a)) > 1E-18) {
    return(a / sqrt(sum(a * a)))
  }
  else 
    stop("Warning: normalize.vector passed vector length close to zero")
} # return vector normalized by length, meaning of unit length
#
# function for dihedral angle
dihedral.angle <- function(a, b, c, d) {
  p <- normalize.vector(b - a)
  q <- normalize.vector(c - b)
  r <- normalize.vector(d - c)
  n1 <- vector.cross(p, q)
  n2 <- vector.cross(q, r)
  m1 <- vector.cross(n1, q)
  x <- sum(n1 * n2)
  y <- sum(m1 * n2)
  return(-180.0 * atan2(y, x) / pi)
} # return the dihedral angle defined by 4 atoms
#
# function to find coordinates for atoms of interest
find.atom.interest <- function(pdbqt) {
  if ((ncol(pdbqt) != 13) & (ncol(pdbqt) != 12))
    stop("warning: unexpected number of columns in .pdbqt file")
  atom.interest <- data.frame("O3" = numeric(3), "O23" = numeric(3), "C1" = numeric(3), "C3" = numeric(3), "C6" = numeric(3), "C7" = numeric(3), "C8" = numeric(3), "C12" = numeric(3), "C21" = numeric(3), "C23" = numeric(3), "C24" = numeric(3), "C26" = numeric(3), "C27" = numeric(3), "C28" = numeric(3), "N" = numeric(3))
  atom.interest.name <- c("O3", "O23", "C1", "C3", "C6", "C7", "C8", "C12", "C21", "C23", "C24", "C26", "C27", "C28", "N")
  ix <- (ncol(pdbqt) - 6)
  iy <- 1 + ix
  iz <- 1 + iy
  for (indexS in 1:nrow(pdbqt)) {
    indexN <- as.integer(1)
    found.atom <- FALSE
    while ((indexN <= length(atom.interest.name)) & (!found.atom)) {
      if(pdbqt[indexS, "V3"] == atom.interest.name[indexN]) {
        atom.interest[1, pdbqt[indexS, "V3"]] <- pdbqt[indexS, ix]
        atom.interest[2, pdbqt[indexS, "V3"]] <- pdbqt[indexS, iy]
        atom.interest[3, pdbqt[indexS, "V3"]] <- pdbqt[indexS, iz]
        found.atom <- TRUE
      } #if atom name matches
      indexN <- 1 + indexN
    } #for each atom of interest check if name matches
  } # while searching for atoms
  return(atom.interest)
} # locate all atoms of interest and assign x, y, z coordinates from pdbqt file
#
# function to find coordinates for atoms of interest
find.atom.interest.Y338 <- function(pdbqt) {
  if ((ncol(pdbqt) != 13) & (ncol(pdbqt) != 12))
    stop("warning: unexpected number of columns in .pdbqt file")
  atom.interest <- data.frame("OH" = numeric(3), "HH" = numeric(3))
  atom.interest.name <- c("OH", "HH")
  ix <- (ncol(pdbqt) - 6)
  iy <- 1 + ix
  iz <- 1 + iy
  for (indexS in c(1,2)) {
    indexN <- as.integer(1)
    found.atom <- FALSE
    while ((indexN <= length(atom.interest.name)) & (!found.atom)) {
      if(pdbqt[indexS, "V3"] == atom.interest.name[indexN]) {
        atom.interest[1, pdbqt[indexS, "V3"]] <- pdbqt[indexS, ix]
        atom.interest[2, pdbqt[indexS, "V3"]] <- pdbqt[indexS, iy]
        atom.interest[3, pdbqt[indexS, "V3"]] <- pdbqt[indexS, iz]
        found.atom <- TRUE
      } #if atom name matches
      indexN <- 1 + indexN
    } #for each atom of interest check if name matches
  } # while searching for atoms
  return(atom.interest)
} # locate all atoms of interest and assign x, y, z coordinates from pdbqt file
#
# Collect energy terms for all ligands and poses
max_mode <- as.integer(5)
dock.id <- as.integer(1)
dietary.xanthophyll <- tibble("dock.id" = integer(), "file.name.ligand.pdbqt" = character(), 
                              "file.name.flex.pdbqt" = character(), "winner" = logical(), 
                              "receptor" = character(), "ligand.name" = character(), 
                              "ligand.name.index" = integer(), "ligand.type" = character(), 
                              "rotation" = integer(), "binding.mode" = integer(), 
                              "buried" = character(), "measure.type" = character(), 
                              "measure.value" = numeric())
index.X <- as.integer(1)
index.ligand.csv.file.name <- as.integer(1)
# ligand <- ligand.code[5]
# rotation.mode <- as.integer(2)
for (ligand in ligand.code) {
  for (rotation.mode in c(2)) {
    if (ligand %in% carotenoid.code) {
      ligand.type <- str_c(ligand, rotation.mode)
    } else { #for carotenoid ligands distinguish between rotating (rotation.mode 5) and rigid (rotation.mode 2)
      ligand.type <- ligand
    } #exception for emixustat related ligands; note that else goes on same line as closing bracket of if {}
    split.directory <- str_c("SPLIT", receptor, ligand.type, sep = "_")
    score.only.directory <- str_c("SCOREONLY", receptor, ligand.type, sep = "_")
    vina.directory <- str_c("OUTPUT", receptor, ligand.type, sep = "_")
    ligand.names.df <- read.csv(str_c(ligand.dir, "/", ligand.type, "/", ligand.name.csv.file[index.ligand.csv.file.name]), header = FALSE, col.names = c("ligand_name"), stringsAsFactors = FALSE)
    index.ligand.csv.file.name <- 1 + index.ligand.csv.file.name
    index.M <- as.integer(1)
    for (ligand.name in ligand.names.df[, "ligand_name"]) {
      # Affinity from vina OUTPUT
      file.name.vina.log <- str_c(vina.directory, "/vina_", receptor, "_", ligand.name, ".log")
      if (file.exists(file.name.vina.log) ) {
        vina.log <- read.delim(file.name.vina.log, header = FALSE, sep = "", stringsAsFactors = FALSE)
      } # check if file exists before opening
      for (binding.mode in c(1:max_mode)) {
        #
        #define filename patterns here
        file.name.scoreonly.log <- str_c(score.only.directory, "/split_", receptor, "_", ligand.name, "_mode_", binding.mode, ".log")
        file.name.split.ligand.pdbqt <- str_c(split.directory, "/split_", receptor, "_", ligand.name, "_ligand_mode_", binding.mode, ".pdbqt")
        file.name.split.flex.pdbqt <- str_c(split.directory, "/split_", receptor, "_", ligand.name, "_flex_mode_", binding.mode, ".pdbqt")
        if (file.exists(file.name.scoreonly.log) ) {
          score.only.log <- read.delim(file.name.scoreonly.log, header = FALSE, sep = ":", stringsAsFactors = FALSE)
          if (nrow(score.only.log) == 39) {
            #
            # Affinity is Binding Energy as reported by Vina Score Only
            dietary.xanthophyll[index.X, "dock.id"] <- dock.id
            dietary.xanthophyll[index.X, "file.name.ligand.pdbqt"] <- str_c(split.directory, "/split_", receptor, "_", ligand.name, "_ligand_mode_", binding.mode, ".pdbqt") 
            dietary.xanthophyll[index.X, "file.name.flex.pdbqt"] <- str_c(split.directory, "/split_", receptor, "_", ligand.name, "_flex_mode_", binding.mode, ".pdbqt") 
            dietary.xanthophyll[index.X, "receptor"] <- receptor
            dietary.xanthophyll[index.X, "ligand.name"] <- ligand.name
            dietary.xanthophyll[index.X, "ligand.name.index"] <- index.M
            dietary.xanthophyll[index.X, "ligand.type"] <- ligand.type
            dietary.xanthophyll[index.X, "rotation"] <- as.integer(rotation.mode)
            dietary.xanthophyll[index.X, "binding.mode"] <- as.integer(binding.mode)
            #
            # find OH and HH atoms belonging to flex residue Y338
            split.flex.pdbqt <- read.delim(file.name.split.flex.pdbqt, header = FALSE, sep = "", stringsAsFactors = FALSE, skip = 52)
            atom.interest.Y338 <- find.atom.interest.Y338(split.flex.pdbqt)
            #
            # Determine which end is buried by distance and also get atom coordinates for torsion angles
            split.ligand.pdbqt <- read.delim(file.name.split.ligand.pdbqt, header = FALSE, sep = "", stringsAsFactors = FALSE, skip = 20)
            atom.interest <- find.atom.interest(split.ligand.pdbqt)
            if (ligand %in% xanthophyll.code) {
              distance.O3 <- sum((receptor.chunnel.landmark - atom.interest[, "O3"]) * (receptor.chunnel.landmark - atom.interest[, "O3"]))
              distance.O3 <- sqrt(distance.O3)
              distance.O23 <- sum((receptor.chunnel.landmark - atom.interest[, "O23"]) * (receptor.chunnel.landmark - atom.interest[, "O23"]))
              distance.O23 <- sqrt(distance.O23)
              if (distance.O3 < distance.O23) {
                dietary.xanthophyll[index.X, "buried"] <- "beta.O3"
              } # if beta ionone ring with O3 is buried
              if (distance.O3 > distance.O23) {
                if (ligand == "STQ") {
                  dietary.xanthophyll[index.X, "buried"] <- "epsilon.O23"
                } # if the other ring is an epsilon ionone ring
                if (ligand != "STQ") {
                  dietary.xanthophyll[index.X, "buried"] <- "beta.O23"
                } # if the other ring is a beta ring
              } # if ionone ring with O23 is buried
            } #if this is a xanthophyll
            if (ligand == "BCR") {
              distance.C3 <- sum((receptor.chunnel.landmark - atom.interest[, "C3"]) * (receptor.chunnel.landmark - atom.interest[, "C3"]))
              distance.C3 <- sqrt(distance.C3)
              distance.C23 <- sum((receptor.chunnel.landmark - atom.interest[, "C23"]) * (receptor.chunnel.landmark - atom.interest[, "C23"]))
              distance.C23 <- sqrt(distance.C23)
              if (distance.C3 < distance.C23) {
                dietary.xanthophyll[index.X, "buried"] <- "beta.C3"
              } # if beta ionone ring with O3 is buried
              if (distance.C3 > distance.C23) {
                dietary.xanthophyll[index.X, "buried"] <- "beta.C23"
              } # if ionone ring with O23 is buried
            } #if beta carotene use C3 and C23 in place of oxygen atoms
            if (ligand == "emixustat") {
              distance.N <- sum((receptor.chunnel.landmark - atom.interest[, "N"]) * (receptor.chunnel.landmark - atom.interest[, "N"]))
              distance.N <- sqrt(distance.N)
              distance.C12 <- sum((receptor.chunnel.landmark - atom.interest[, "C12"]) * (receptor.chunnel.landmark - atom.interest[, "C12"]))
              distance.C12 <- sqrt(distance.C12)
              if (distance.N < distance.C12) {
                dietary.xanthophyll[index.X, "buried"] <- "amine.N"
              } # if primary amine is closer to landmark
              if (distance.N > distance.C12) {
                dietary.xanthophyll[index.X, "buried"] <- "difluoro.C12"
              } # if difluoro group is closer to landmark
            } #if difluoro-emixustat use N and C12 with difluoro group for orientation determination
            dietary.xanthophyll[index.X, "measure.type"] <- "Binding.Energy.Score.Only"
            row.of.interest <- which(str_detect(score.only.log[ ,1], "Estimated"))
            position.kcal <- str_locate(score.only.log[row.of.interest, 2], "kcal")
            dietary.xanthophyll[index.X, "measure.value"] <- as.double(str_sub(score.only.log[row.of.interest, 2], 1, (position.kcal[1, "start"] -2)))
            indexbase <- index.X
            index.X <- 1 + index.X
            #
            # ligand internal energy
            dietary.xanthophyll[index.X,] <- dietary.xanthophyll[indexbase,]
            dietary.xanthophyll[index.X, "measure.type"] <- "Ligand.Internal.Energy"
            row.of.interest <- which(str_detect(score.only.log[ ,1], "   Ligand                         "))
            position.kcal <- str_locate(score.only.log[row.of.interest, 2], "kcal")
            dietary.xanthophyll[index.X, "measure.value"] <- as.double(str_sub(score.only.log[row.of.interest, 2], 1, (position.kcal[1, "start"] -2)))
            index.X <- 1 + index.X
            #
            # Affinity from Vina Output file
            if (nrow(vina.log) > 36) {
              dietary.xanthophyll[index.X,] <- dietary.xanthophyll[indexbase,]
              dietary.xanthophyll[index.X, "measure.type"] <- "Binding.Energy.Vina"
              row.of.interest <- which(str_detect(vina.log[ , 1], as.character(binding.mode)))
              dietary.xanthophyll[index.X, "measure.value"] <- as.double(vina.log[row.of.interest, 2])
              index.X <- 1 + index.X
              #
              # difference Vina and Score Only energy
              dietary.xanthophyll[index.X,] <- dietary.xanthophyll[indexbase,]
              dietary.xanthophyll[index.X, "measure.type"] <- "Binding.Energy.Difference"
              dietary.xanthophyll[index.X, "measure.value"] <- dietary.xanthophyll[dietary.xanthophyll$dock.id == dock.id & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", "measure.value"]
              dietary.xanthophyll[index.X, "measure.value"] <- dietary.xanthophyll[index.X, "measure.value"] - 
                dietary.xanthophyll[dietary.xanthophyll$dock.id == dock.id & dietary.xanthophyll$measure.type == "Binding.Energy.Vina", "measure.value"]
              index.X <- 1 + index.X
            } # check there are data in the vina log file
            #
            # torsion angles for ionone rings
            if (ligand %in% carotenoid.code) {
              dietary.xanthophyll[index.X,] <- dietary.xanthophyll[indexbase,]
              dietary.xanthophyll[index.X, "measure.type"] <- "torsion.c6c7"
              dietary.xanthophyll[index.X, "measure.value"] <- dihedral.angle(atom.interest[, "C1"], atom.interest[, "C6"], atom.interest[, "C7"], atom.interest[, "C8"])
              index.X <- 1 + index.X
              dietary.xanthophyll[index.X,] <- dietary.xanthophyll[indexbase,]
              dietary.xanthophyll[index.X, "measure.type"] <- "torsion.c26c27"
              dietary.xanthophyll[index.X, "measure.value"] <- dihedral.angle(atom.interest[, "C21"], atom.interest[, "C26"], atom.interest[, "C27"], atom.interest[, "C28"])
              index.X <- 1 + index.X
            } #calculate torsion angles for carotenoids but not emixustat
            #
            # distances to landmark
            dietary.xanthophyll[index.X,] <- dietary.xanthophyll[indexbase,]
            dietary.xanthophyll[index.X, "measure.type"] <- "distance.landmark"
            if (ligand %in% xanthophyll.code) {
              dietary.xanthophyll[index.X, "measure.value"] <- min(distance.O23, distance.O3)
            } # if ligand is a xanthophyll use O23 or O3 distances
            if (ligand == "BCR") {
              dietary.xanthophyll[index.X, "measure.value"] <- min(distance.C23, distance.C3)
            } # if ligand is beta carotene use C23 or C3 distances
            if (ligand == "emixustat") {
              dietary.xanthophyll[index.X, "measure.value"] <- min(distance.N, distance.C12)
            } # if ligand is emixustat use N or C12 distances
            index.X <- 1 + index.X
            
            if (ligand %in% carotenoid.code) {
              #
              # distance Y338@HH to C24
              dietary.xanthophyll[index.X,] <- dietary.xanthophyll[indexbase,]
              dietary.xanthophyll[index.X, "measure.type"] <- "distance.HH.C24"
              distance.HH.C24 <- sum((atom.interest.Y338[, "HH"] - atom.interest[, "C24"]) * (atom.interest.Y338[, "HH"] - atom.interest[, "C24"]))
              distance.HH.C24 <- sqrt(distance.HH.C24)
              dietary.xanthophyll[index.X, "measure.value"] <- distance.HH.C24
              index.X <- 1 + index.X
              #
              # distance Y338@HH to C26
              dietary.xanthophyll[index.X,] <- dietary.xanthophyll[indexbase,]
              dietary.xanthophyll[index.X, "measure.type"] <- "distance.HH.C26"
              distance.HH.C26 <- sum((atom.interest.Y338[, "HH"] - atom.interest[, "C26"]) * (atom.interest.Y338[, "HH"] - atom.interest[, "C26"]))
              distance.HH.C26 <- sqrt(distance.HH.C26)
              dietary.xanthophyll[index.X, "measure.value"] <- distance.HH.C26
              index.X <- 1 + index.X
            } #calculate distance to Y338@HH for carotenoids only
            
            #
            # all done, advance to the next docking outcome
            dock.id <- 1 + dock.id
          } #if the log file contains energy terms
        } #if the score only file exists, open and extract energy values
      } #for first max_mode binding modes, only first one is guaranteed to exist
      index.M <- 1 + index.M
      } # for each ligand in ensemble derived from light harvesting complexes
    } # for each rotation mode
} # for each ligand code: emixustat, STQ, VFQ, 8UH, BCR
write.csv(dietary.xanthophyll, 
          str_c(receptor, "rigid.ligands.emixustat.xanthophylls.beta.carotene.binding.energy", "csv", sep = "."),
          row.names = FALSE, quote = FALSE)
#
# flag.loser to set winner flag FALSE for each card belonging to named pdbqt files
flag.loser <- function(xanthophyll, loser) {
  index.F <- as.integer(1)
  while (index.F <= nrow(loser)) {
    xanthophyll[xanthophyll$dock.id == loser[[index.F, "dock.id"]], "winner"] <- FALSE
    index.F <- index.F + 1
  } # for each identified loser set winner flag to FALSE
  return(xanthophyll)
} # sets winner flag to FALSE if named in loser
#
# Build the list of winners with high priority; combine criteria by omitting the <- TRUE command
index.R <- as.integer(1)
ligand.rotation <- c("rigid", "rotate")
frequency.df <- data.frame("ligand.rotation" = character(), "ligand.type" = character(), "buried" = character(), "frequency.winner" = integer(), "frequency.loser" = integer(), "frequency.total" = integer())
index.F <- as.integer(1)
for (rotation.mode in c(2)) {
  dietary.xanthophyll[, "winner"] <- TRUE
  frequency.df[index.F, "ligand.rotation"] <- ligand.rotation[index.R]
  frequency.df[index.F, "ligand.type"] <- "all"
  frequency.df[index.F, "buried"] <- "either"
  # remove rigid or rotating ligands
  loser.identity <- dietary.xanthophyll[dietary.xanthophyll$rotation != rotation.mode, ]
  dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
  frequency.df[index.F, "frequency.total"] <- 
    nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
  binding.energy <- data.frame(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", "measure.value"])
  top.10p <- quantile(binding.energy[, "measure.value"], c(0, 0.1, 0.9))[2]
  frequency.df[index.F, "frequency.winner"] <- 
    nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only"
                             & dietary.xanthophyll$measure.value < top.10p, ])
  frequency.df[index.F, "frequency.loser"] <- 
    nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only"
                             & dietary.xanthophyll$measure.value >= top.10p, ])
  index.F <- 1 + index.F
  for (ligand in ligand.code) {
    if (ligand %in% carotenoid.code) {
      ligand.type <- str_c(ligand, rotation.mode)
    } else {
      ligand.type <- ligand
    } # if not a carotenoid, i.e. emixustat ligand
    dietary.xanthophyll[, "winner"] <- TRUE
    loser.identity <- dietary.xanthophyll[dietary.xanthophyll$ligand.type != ligand.type, ]
    dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
    frequency.df[index.F, "ligand.rotation"] <- ligand.rotation[index.R]
    frequency.df[index.F, "ligand.type"] <- ligand.type
    frequency.df[index.F, "buried"] <- "either"
    frequency.df[index.F, "frequency.total"] <- 
      nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
    frequency.df[index.F, "frequency.winner"] <- 
      nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only"
                               & dietary.xanthophyll$measure.value < top.10p, ])
    frequency.df[index.F, "frequency.loser"] <- 
      nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only"
                               & dietary.xanthophyll$measure.value >= top.10p, ])
    index.F <- 1 + index.F
    for (buried in c("O23", "O3", "C23", "C3", "N", "C12")) {
      frequency.df[index.F, "ligand.rotation"] <- ligand.rotation[index.R]
      frequency.df[index.F, "ligand.type"] <- ligand.type
      frequency.df[index.F, "buried"] <- buried
      frequency.df[index.F, "frequency.total"] <-
        nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only"
                                 & str_detect(dietary.xanthophyll$buried, buried), ])
      frequency.df[index.F, "frequency.winner"] <-
        nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only"
                                 & str_detect(dietary.xanthophyll$buried, buried)
                                 & dietary.xanthophyll$measure.value < top.10p, ])
      frequency.df[index.F, "frequency.loser"] <-
        nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only"
                                 & str_detect(dietary.xanthophyll$buried, buried)
                                 & dietary.xanthophyll$measure.value >= top.10p, ])
      index.F <- 1 + index.F
    } #separate on basis of being buried
  } # calculate frequency for each ligand type
  #
  # Box Plot for binding energy score only by ligand type for winner top 10%
  dietary.xanthophyll[, "winner"] <- TRUE
  # remove rigid or rotating ligands
  loser.identity <- dietary.xanthophyll[dietary.xanthophyll$rotation != rotation.mode, ]
  dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
  loser.identity <- dietary.xanthophyll[dietary.xanthophyll$buried %in% c("beta.O3", "beta.C3"), ]
  dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
  #
  # select the top 10% winners for each ligand
  for (ligand in ligand.code) {
    if (ligand %in% carotenoid.code) {
      ligand.type <- str_c(ligand, rotation.mode)
    } else {
      ligand.type <- ligand
    } # if ligand is not carotenoid; i.e. emixustat
    # ligand.type <- "emixustat"
    binding.energy <- data.frame(dietary.xanthophyll[dietary.xanthophyll$winner & 
                                                       dietary.xanthophyll$ligand.type == ligand.type &
                                                       dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", "measure.value"])
    top.10p <- quantile(binding.energy[, "measure.value"], c(0, 0.1, 0.5, 0.9))[2]
    loser.identity <- dietary.xanthophyll[dietary.xanthophyll$winner & 
                                            dietary.xanthophyll$ligand.type == ligand.type &
                                            dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" & 
                                            dietary.xanthophyll$measure.value >= top.10p, ]
    dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
  } # for each ligand select the top 10%
  
  myplot <- ggplot(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type %in% c("Binding.Energy.Score.Only"), ], aes(x=ligand.type, y=measure.value)) +
    coord_cartesian(ylim = c(-15, -10.5), expand = TRUE, default = TRUE) +
    labs(y="Binding Energy Score Only", x="Ligand Type", title=str_c(receptor, " docked with dietary xanthophyll")) +
    theme_linedraw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_boxplot(alpha = 0.5, aes(color = ligand.type, fill = ligand.type), show.legend = FALSE, )
  ggsave(str_c(receptor, "C23buried.rotation.mode", ligand.rotation[index.R], "ligands.top10p", "Binding.Energy.ScoreOnly.BoxPlot.pdf", sep = "."), height=4, width=5.5, dpi=300)
  ggsave(str_c(receptor, "C23buried.rotation.mode", ligand.rotation[index.R], "ligands.top10p", "Binding.Energy.ScoreOnly.BoxPlot.png", sep = "."), height=4, width=5.5, dpi=300)
  # index.R <- 1 + index.R
} # repeat for both rigid and rotating ligands
#
# Function for writing a cif file
write.cif <- function(file.name, structure.df) {
  cif.header <- read.delim(str_c(input.directory, "/cif.header"), header = FALSE, sep = "", stringsAsFactors = FALSE, skip = 0)
  cif.df <- tibble("textline" = character())
  cif.df[1:nrow(cif.header), "textline"] <- cif.header
  structure.df[structure.df$t == "A", "t"] <- "C"
  index.C <- 1 + nrow(cif.df)
  index.A <- as.integer(1)
# pdb.column.names <- c("CARD", "ATOM.INDEX", "ATOM.NAME", "RES.TYPE", "CHAIN", "RES.ID", "x", "y", "z", "occ", "b", "t")
  while (index.A <= nrow(structure.df)) {
    if (structure.df[[index.A, "ATOM.NAME"]] == "FE") {
      element.length <- as.integer(2)
    } else {
      element.length <- as.integer(1)
    }
    cif.df[index.C, "textline"] <- str_c(structure.df[[index.A, "CARD"]], 
                                         structure.df[[index.A, "ATOM.INDEX"]],
                                         str_sub(structure.df[[index.A, "ATOM.NAME"]], 1, element.length),
                                         structure.df[[index.A, "ATOM.NAME"]],
                                         ".",
                                         structure.df[[index.A, "RES.TYPE"]],
                                         structure.df[[index.A, "CHAIN"]],
                                         "1",
                                         structure.df[[index.A, "RES.ID"]],
                                         "   . . ",
                                         structure.df[[index.A, "x"]],
                                         structure.df[[index.A, "y"]],
                                         structure.df[[index.A, "z"]],
                                         structure.df[[index.A, "occ"]],
                                         structure.df[[index.A, "b"]],
                                         " ? ? ? ? ? ? ",
                                         structure.df[[index.A, "RES.ID"]],
                                         structure.df[[index.A, "RES.TYPE"]],
                                         structure.df[[index.A, "CHAIN"]],
                                         structure.df[[index.A, "ATOM.NAME"]],
                                         "1",
                                         sep = " ")
    index.C <- 1 + index.C
    index.A <- 1 + index.A
  } # for each atom in the structure
  cif.df[index.C, "textline"] <- "# "
  index.C <- 1 + index.C
  write.table(cif.df, str_c(file.name, ".cif"), row.names = FALSE, quote = FALSE, col.names = FALSE)
  rm(cif.df)
  return()
} # writes coordinates as a .cif file

#
# Build the winner table and sort by binding energy as determined by score only
for (numberofwinners in c(10)) {
  all.winner.xanthophyll <- NULL
  for (ligand in ligand.code) {
    for (rotation.mode in c("2")) {
      if (ligand %in% carotenoid.code) {
        ligand.type <- str_c(ligand, rotation.mode)
      } else {
        ligand.type <- ligand
      }
      dietary.xanthophyll[, "winner"] <- TRUE
      # only accept this ligand type
      loser.identity <- dietary.xanthophyll[dietary.xanthophyll$ligand.type != ligand.type, ]
      dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
      # remove the buried O3 examples
      loser.identity <- dietary.xanthophyll[dietary.xanthophyll$buried %in% c("beta.O3", "beta.C3"), ]
      dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
      # rank by binding energy as measured by score only mode
      winner.xanthophyll <- dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ]
      by.binding.energy <- order(winner.xanthophyll[["measure.value"]])
      # merge into one data frame for all ligand types
      all.winner.xanthophyll <- rbind(all.winner.xanthophyll, winner.xanthophyll[by.binding.energy, ][1:numberofwinners, ])
    } # for rigid and rotating ligands
  } # for each ligand type 
  #
  # extract other measurements of interest for the winners
  all.winner.dietary.xanthophyll <- NULL
  index.W <- as.integer(1)
  for (index.W in 1:nrow(all.winner.xanthophyll)) {
    dock.id <- all.winner.xanthophyll[[index.W, "dock.id"]]
    all.winner.dietary.xanthophyll <- rbind(all.winner.dietary.xanthophyll, dietary.xanthophyll[dietary.xanthophyll$dock.id == dock.id, ])
  } # for each winner get measurements of interest
  #
  # plot measurements of interest binding energy score only
#  myplot <- ggplot(all.winner.xanthophyll[all.winner.xanthophyll$measure.type == "Binding.Energy.Score.Only" &
#                                                    all.winner.xanthophyll$buried %in% c("epsilon.O23", "beta.O23", "beta.C23", "amine.N") & , ], aes(x=ligand.type, y=measure.value)) +
  myplot <- ggplot(all.winner.xanthophyll[all.winner.xanthophyll$measure.value < -11.3, ], aes(x=ligand.type, y=measure.value)) +
    coord_cartesian(ylim = c(-15, -11), expand = TRUE, default = TRUE) +
    labs(y="Binding Energy Score Only", x="Ligand Type", title=str_c(receptor, " rigid.C23.buried", as.character(numberofwinners), "winner.ligands")) +
    theme_linedraw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_boxplot(alpha = 0.5, aes(color = ligand.type, fill = ligand.type), show.legend = FALSE)
  ggsave(str_c(receptor, "rigid.C23.buried", as.character(numberofwinners), "winner.ligands", "Binding.Energy.Score.Only.BoxPlot.pdf", sep = "."), height=3, width=5.5, dpi=300)
  #
  # plot measurements of interest distance to landmark
  myplot <- ggplot(all.winner.dietary.xanthophyll[all.winner.dietary.xanthophyll$measure.type == "distance.landmark", ], aes(x=ligand.type, y=measure.value)) +
    coord_cartesian(ylim = c(2, 6), expand = TRUE, default = TRUE) +
    labs(y="Distance to landmark atom", x="Ligand Type", title=str_c(receptor, " docked with dietary xanthophyll")) +
    theme_linedraw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_boxplot(alpha = 0.5, aes(color = ligand.type, fill = ligand.type))
  ggsave(str_c(receptor, "all", as.character(numberofwinners), "winner.ligands", "Distance.Landmark.BoxPlot.pdf", sep = "."), height=6, width=5.5, dpi=300)
  #
  # plot measurements of interest distance TYR338@HH to C24
  myplot <- ggplot(all.winner.dietary.xanthophyll[all.winner.dietary.xanthophyll$measure.type == "distance.HH.C24", ], aes(x=ligand.type, y=measure.value)) +
    coord_cartesian(ylim = c(0, 12), expand = TRUE, default = TRUE) +
    labs(y="Distance TYR338@HH to C24", x="Ligand Type", title=str_c(receptor, " docked with dietary xanthophyll")) +
    theme_linedraw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_boxplot(alpha = 0.5, aes(color = ligand.type, fill = ligand.type))
  ggsave(str_c(receptor, "all", as.character(numberofwinners), "winner.ligands", "Distance.TYR338HH.C24.BoxPlot.pdf", sep = "."), height=6, width=5.5, dpi=300)
  #
  # plot measurements of interest distance TYR338@HH to C26
  myplot <- ggplot(all.winner.dietary.xanthophyll[all.winner.dietary.xanthophyll$measure.type == "distance.HH.C26", ], aes(x=ligand.type, y=measure.value)) +
    coord_cartesian(ylim = c(0, 12), expand = TRUE, default = TRUE) +
    labs(y="Distance TYR338@HH to C26", x="Ligand Type", title=str_c(receptor, " docked with dietary xanthophyll")) +
    theme_linedraw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_boxplot(alpha = 0.5, aes(color = ligand.type, fill = ligand.type))
  ggsave(str_c(receptor, "all", as.character(numberofwinners), "winner.ligands", "Distance.TYR338HH.C26.BoxPlot.pdf", sep = "."), height=6, width=5.5, dpi=300)
  #
  # plot measurements of interest ligand internal energy
  myplot <- ggplot(all.winner.dietary.xanthophyll[all.winner.dietary.xanthophyll$measure.type == "Ligand.Internal.Energy", ], aes(x=ligand.type, y=measure.value)) +
    coord_cartesian(ylim = c(-2, 2), expand = TRUE, default = TRUE) +
    labs(y="Ligand internal energy", x="Ligand Type", title=str_c(receptor, " docked with dietary xanthophyll")) +
    theme_linedraw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_boxplot(alpha = 0.5, aes(color = ligand.type, fill = ligand.type))
  ggsave(str_c(receptor, "all", as.character(numberofwinners), "winner.ligands", "Ligand.Internal.Energy.BoxPlot.pdf", sep = "."), height=6, width=5.5, dpi=300)
  #
  # plot measurements of interest c26-c27 torsion angle
  myplot <- ggplot(all.winner.dietary.xanthophyll[all.winner.dietary.xanthophyll$measure.type == "torsion.c26c27", ], aes(x=ligand.type, y=measure.value)) +
    coord_cartesian(ylim = c(-180, 180), expand = TRUE, default = TRUE) +
    labs(y="C26-C27 torsion", x="Ligand Type", title=str_c(receptor, " docked with dietary xanthophyll")) +
    theme_linedraw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    geom_boxplot(alpha = 0.5, aes(color = ligand.type, fill = ligand.type))
  ggsave(str_c(receptor, "all", as.character(numberofwinners), "winner.ligands", "C26-C27.Torsion.Angle.BoxPlot.pdf", sep = "."), height=6, width=5.5, dpi=300)
  #
  # plot measurements of interest against each other requiring reorganization of the data
  dietary.xanthophyll.matrix <- data.frame(dock.id = integer(), buried = character(),
                                           ligand.type = character(),
                                           Binding.Energy.Score.Only = numeric(),
                                           Binding.Energy.Vina = numeric(), 
                                           torsion.c6c7 = numeric(),
                                           torsion.c26c27 = numeric(),
                                           Ligand.Internal.Energy = numeric())
  index.W <- as.integer(1)
  for (index.W in 1:nrow(all.winner.xanthophyll)) {
    dock.id <- all.winner.xanthophyll[[index.W, "dock.id"]]
    dietary.xanthophyll.matrix[index.W, "dock.id"] <- dock.id
    dietary.xanthophyll.matrix[index.W, "buried"] <- all.winner.xanthophyll[index.W, "buried"]
    dietary.xanthophyll.matrix[index.W, "ligand.type"] <- all.winner.xanthophyll[index.W, "ligand.type"]
    dietary.xanthophyll.matrix[index.W, "Binding.Energy.Score.Only"] <- all.winner.xanthophyll[index.W, "measure.value"]
    dietary.xanthophyll.matrix[index.W, "Binding.Energy.Vina"] <- all.winner.dietary.xanthophyll[all.winner.dietary.xanthophyll$dock.id == dock.id &
                                                                                                   all.winner.dietary.xanthophyll$measure.type == "Binding.Energy.Vina", 
                                                                                                 "measure.value"]
    dietary.xanthophyll.matrix[index.W, "torsion.c6c7"] <- all.winner.dietary.xanthophyll[all.winner.dietary.xanthophyll$dock.id == dock.id &
                                                                                            all.winner.dietary.xanthophyll$measure.type == "torsion.c6c7", 
                                                                                          "measure.value"]
    dietary.xanthophyll.matrix[index.W, "torsion.c26c27"] <- all.winner.dietary.xanthophyll[all.winner.dietary.xanthophyll$dock.id == dock.id &
                                                                                              all.winner.dietary.xanthophyll$measure.type == "torsion.c26c27", 
                                                                                            "measure.value"]
    dietary.xanthophyll.matrix[index.W, "Ligand.Internal.Energy"] <- all.winner.dietary.xanthophyll[all.winner.dietary.xanthophyll$dock.id == dock.id &
                                                                                                      all.winner.dietary.xanthophyll$measure.type == "Ligand.Internal.Energy", 
                                                                                                    "measure.value"]
    dietary.xanthophyll.matrix[index.W, "distance.HH.C24"] <- all.winner.dietary.xanthophyll[all.winner.dietary.xanthophyll$dock.id == dock.id &
                                                                                                      all.winner.dietary.xanthophyll$measure.type == "distance.HH.C24", 
                                                                                                    "measure.value"]
    dietary.xanthophyll.matrix[index.W, "distance.HH.C26"] <- all.winner.dietary.xanthophyll[all.winner.dietary.xanthophyll$dock.id == dock.id &
                                                                                                      all.winner.dietary.xanthophyll$measure.type == "distance.HH.C26", 
                                                                                                    "measure.value"]
  } # for each winner get measurements of interest into a matrix for easier plotting
  myplot <- ggplot(dietary.xanthophyll.matrix, aes(x=Binding.Energy.Vina, y=Binding.Energy.Score.Only)) +
    coord_cartesian(xlim = c(-15.0, -0.0), ylim = c(-15, -0), expand = TRUE, default = TRUE) +
    labs(y="Binding Energy Score Only", x="Binding Energy Vina Output", title=str_c(receptor, " docked with dietary xanthophyll ", ligand.type)) +
    theme_linedraw(base_size = 12) +
    geom_point(alpha = 0.50, aes(color = ligand.type))
  ggsave(str_c(receptor, "all.ligand", as.character(numberofwinners), "winners", "Binding.Energy.Linear.pdf", sep = "."), height=6, width=5.5, dpi=300)
  #
  # plot ligand internal energy as a function of torsion angles
  myplot <- ggplot(dietary.xanthophyll.matrix, aes(x=torsion.c26c27, y=Ligand.Internal.Energy)) +
    coord_cartesian(ylim = c(2, 8), xlim = c(2, 8), expand = TRUE, default = TRUE) +
    labs(x="C26-C27 torsion", y="Ligand Internal Energy", title=str_c(receptor, " docked with dietary xanthophyll ", ligand.type)) +
    theme_linedraw(base_size = 12) +
    geom_point(alpha = 0.5, aes(color = ligand.type))
  ggsave(str_c(receptor, "all.ligand", as.character(numberofwinners), "winners", "Internal.Energy.C26-C27.Torsion.Angle.pdf", sep = "."), height=6, width=5.5, dpi=300)
  #
  # plot TYR338@HH to C24 distance vs TYR338@HH to C26 distance
  myplot <- ggplot(dietary.xanthophyll.matrix, aes(x=distance.HH.C24, y=distance.HH.C26)) +
    coord_cartesian(ylim = c(2, 8), xlim = c(2, 8), expand = TRUE, default = TRUE) +
    labs(x="distance Y338@HH to C24", y="distance Y338@HH to C26", title=str_c(receptor, " docked with dietary xanthophyll ", ligand.type)) +
    theme_linedraw(base_size = 12) +
    geom_point(alpha = 0.5, aes(color = ligand.type))
  ggsave(str_c(receptor, "all.ligand", as.character(numberofwinners), "winners", "TYR338HH.C24vsC26.distance.pdf", sep = "."), height=6, width=5.5, dpi=300)
  #
  # repeat plots for individual ligand types
  for (ligand in carotenoid.code) {
    for (rotation.mode in c("2")) {
      ligand.type = str_c(ligand, rotation.mode)
      # make the plot comparing binding energy measured by score only and by vina output 
      myplot <- ggplot(dietary.xanthophyll.matrix[dietary.xanthophyll.matrix$ligand.type == ligand.type, ], aes(x=Binding.Energy.Vina, y=Binding.Energy.Score.Only)) +
        coord_cartesian(xlim = c(-15.0, -0.0), ylim = c(-15, -0), expand = TRUE, default = TRUE) +
        labs(y="Binding Energy Score Only", x="Binding Energy Vina Output", title=str_c(receptor, " docked with dietary xanthophyll ", ligand.type)) +
        theme_linedraw(base_size = 12) +
        geom_point(alpha = 0.50, aes(color = ligand.type))
      ggsave(str_c(receptor, ligand.type, as.character(numberofwinners), "winners", "Binding.Energy.Linear.pdf", sep = "."), height=6, width=5.5, dpi=300)
      #
      # plot ligand internal energy as a function of torsion angles
      myplot <- ggplot(dietary.xanthophyll.matrix[dietary.xanthophyll.matrix$ligand.type == ligand.type, ], aes(x=torsion.c26c27, y=Ligand.Internal.Energy)) +
        coord_cartesian(ylim = c(-2, 2), expand = TRUE, default = TRUE) +
        labs(x="C26-C27 torsion", y="Ligand Internal Energy", title=str_c(receptor, " docked with dietary xanthophyll ", ligand.type)) +
        theme_linedraw(base_size = 12) +
        geom_point(alpha = 0.5, aes(color = ligand.type))
      ggsave(str_c(receptor, ligand.type, as.character(numberofwinners), "winners", "Internal.Energy.C26-C27.Torsion.Angle.pdf", sep = "."), height=6, width=5.5, dpi=300)
      #
      # plot TYR338@HH to C24 distance vs TYR338@HH to C26 distance
      myplot <- ggplot(dietary.xanthophyll.matrix, aes(x=distance.HH.C24, y=distance.HH.C26)) +
        coord_cartesian(ylim = c(2, 8), xlim = c(2, 8), expand = TRUE, default = TRUE) +
        labs(x="distance Y338@HH to C24", y="distance Y338@HH to C26", title=str_c(receptor, " docked with dietary xanthophyll ", ligand.type)) +
        theme_linedraw(base_size = 12) +
        geom_point(alpha = 0.5, aes(color = ligand.type))
      ggsave(str_c(receptor, ligand.type, as.character(numberofwinners), "winners", "TYR338HH.C24vsC26.distance.pdf", sep = "."), height=6, width=5.5, dpi=300)
    } # for rigid and rotating ligands
  } # for each ligand STQ, VFQ, 8UH
} # create plots for each type of ligand for different categories top 50, top 5, top 10)
#
# build the ChimeraX script to open and view the winner ligands
ChimeraX.script <- tibble("command" = character())
index.X <- as.integer(1)
ChimeraX.script[index.X, "command"] <- "close"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- str_c("cd ", work.dir)
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- str_c("open ", file.name.receptor.pdbqt)
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- str_c("open ", file.name.receptor.pdbqt)
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "surface #1"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "color #1 cornflower blue target s"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "transparency #1 50 target s"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "lighting simple"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "set bgColor white"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "graphics silhouettes true"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "display #2/A:601,180,241,313,527"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "display #2/A:148,239,275"
index.X <- 1 + index.X
# ChimeraX.script[index.X, "command"] <- "display #2/A:TYR"
# index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "view #2/A:601 @< 6.5"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "hide #!*.1 models"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "distance #2/A:601 #2/A:180@NE2"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "distance #2/A:601 #2/A:241@NE2"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "distance #2/A:601 #2/A:313@NE2"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "distance #2/A:601 #2/A:527@NE2"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "~ribbon #1"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "hide #!*.1 models"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "clip model #3-83 off"
index.X <- 1 + index.X
#
#open winner xanthophylls, restricting to rigid lutein and meso-zeaxanthin ligands without rotation
all.winner.xanthophyll <- all.winner.xanthophyll[all.winner.xanthophyll$ligand.type %in% c("STQ2", "VFQ2", "8UH2", "BCR2"), ]
index.F <- as.integer(1)
while (index.F <= nrow(all.winner.xanthophyll)) {
  ChimeraX.script[index.X, "command"] <- str_c("open ", all.winner.xanthophyll[[index.F, "file.name.ligand.pdbqt"]])
  index.X <- 1 + index.X
  ChimeraX.script[index.X, "command"] <- str_c("open ", all.winner.xanthophyll[[index.F, "file.name.flex.pdbqt"]])
  index.X <- 1 + index.X
  index.F <- 1 + index.F
} # for each winner.xanthophyll open its ligand and its flexible residues
#
# count potential hydrogen bond partners for each winner xanthophyll
# hydrogen.bond.cutoff <- as.numeric(3.2)
# index.F <- as.integer(1)
# while (index.F <= nrow(all.winner.xanthophyll)) {
#   model.index <- 2.0* index.F + 1
#   ChimeraX.script[index.X, "command"] <- str_c("select #", model.index, "@O23 @<", hydrogen.bond.cutoff, " & #1,", (model.index+1), "@N*,O*")
#   index.X <- 1 + index.X
#   index.F <- 1 + index.F
# } # for each winner.xanthophyll test if there are any electronegative atoms within hydrogen bonding distance of O23 and O3
write.table(ChimeraX.script, str_c(receptor, "all.ligand.types", "open.winner.xanthophyll.cxc", sep = "."), row.names = FALSE, quote = FALSE, col.names = FALSE)
#
# Merge missing regions, flexible, and rigid portions of a receptor structure as starting point for MD
# Add the docked ligand and write the structures as .cif files with and without the ligand
#
# Also build the ChimeraX script to open the cif files for receptor + ligand, 
# add hydrogens and save as PDB files
ChimeraX.script <- tibble("command" = character())
index.X <- as.integer(1)
ChimeraX.script[index.X, "command"] <- "close all"
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- str_c("cd ", gromacs.dir)
index.X <- 1 + index.X
index.M <- as.integer(1)
#
missing.receptor.pdb <- read.delim(str_c(input.directory, "/AF-Q28175-F1-model_v4_aligned-coot-2_missing.pdb"), header = FALSE, sep = "", stringsAsFactors = FALSE, skip = 133)
rigid.receptor.notmissing.pdb <- read.delim(str_c(input.directory, "/7l0e.cif-coot-2_clean_rigid_notmissing.pdb"), header = FALSE, sep = "", stringsAsFactors = FALSE, skip = 46)
pdbqt.column.names <- c("CARD", "ATOM.INDEX", "ATOM.NAME", "RES.TYPE", "CHAIN", "RES.ID", "x", "y", "z", "occ", "b", "q", "t")
pdb.column.names <- c("CARD", "ATOM.INDEX", "ATOM.NAME", "RES.TYPE", "CHAIN", "RES.ID", "x", "y", "z", "occ", "b", "t")
colnames(missing.receptor.pdb) <- pdb.column.names
colnames(rigid.receptor.notmissing.pdb) <- pdb.column.names
#index.F <- as.integer(1)
for (index.F in c(1, 9, 13, 15)) {
  if (file.exists(all.winner.xanthophyll[[index.F, "file.name.flex.pdbqt"]])) {
    receptor.flex.pdbqt <- read.delim(all.winner.xanthophyll[[index.F, "file.name.flex.pdbqt"]], header = FALSE, sep = "", stringsAsFactors = FALSE, skip = 6)
    colnames(receptor.flex.pdbqt) <- pdbqt.column.names
    receptor.whole.pdbqt <- rbind(receptor.flex.pdbqt[receptor.flex.pdbqt$CARD == "ATOM", pdb.column.names],
                                  missing.receptor.pdb[missing.receptor.pdb$CARD == "ATOM", ], 
                                  rigid.receptor.notmissing.pdb[rigid.receptor.notmissing.pdb$CARD %in% c("ATOM", "HETATM"), ])
    by.resid <- order(as.integer(receptor.whole.pdbqt[ , "RES.ID"]))
    receptor.whole.pdbqt <- receptor.whole.pdbqt[by.resid, ]
    receptor.whole.pdbqt[, "ATOM.INDEX"] <- 1:nrow(receptor.whole.pdbqt)
    row.names(receptor.whole.pdbqt) <- 1:nrow(receptor.whole.pdbqt)
    receptor.whole.pdbqt[receptor.whole.pdbqt$t == "A", "t"] <- "C"
    file.name.whole.receptor <- str_c(gromacs.dir, "/", receptor, "_", all.winner.xanthophyll[[index.F, "ligand.name"]], "_mode_", all.winner.xanthophyll[[index.F, "binding.mode"]])
    write.cif(file.name.whole.receptor, receptor.whole.pdbqt)
    ligand.pdbqt <- read.delim(all.winner.xanthophyll[[index.F, "file.name.ligand.pdbqt"]], header = FALSE, sep = "", stringsAsFactors = FALSE, skip = 20)
    if (ncol(ligand.pdbqt) == 12) {
      ligand.pdbqt <- cbind(ligand.pdbqt[, 1:4], ligand.pdbqt[, 5], ligand.pdbqt[, 5:12])
      ligand.pdbqt[ , 5] <- str_sub(ligand.pdbqt[ , 5], 1, 1)
      ligand.pdbqt[ , 6] <- str_sub(ligand.pdbqt[ , 6], 2, -1)
    } #separate the chain from the resid if these are merged
    colnames(ligand.pdbqt) <- pdbqt.column.names
    file.name.ligand <- str_c(gromacs.dir, "/", all.winner.xanthophyll[[index.F, "ligand.name"]], "_mode_", all.winner.xanthophyll[[index.F, "binding.mode"]])
    file.name.whole.receptor.ligand <- str_c(gromacs.dir, "/", receptor, "_", all.winner.xanthophyll[[index.F, "ligand.name"]], "_mode_", all.winner.xanthophyll[[index.F, "binding.mode"]], ".with.ligand")
    receptor.ligand.complex <- rbind(receptor.whole.pdbqt, ligand.pdbqt[ligand.pdbqt$CARD %in% c("ATOM", "HETATM"), pdb.column.names])
    receptor.ligand.complex[, "ATOM.INDEX"] <- 1:nrow(receptor.ligand.complex)
    row.names(receptor.ligand.complex) <- 1:nrow(receptor.ligand.complex)
    receptor.ligand.complex[receptor.ligand.complex$t == "A", "t"] <- "C"
    write.cif(file.name.whole.receptor.ligand, receptor.ligand.complex)
    #
    ChimeraX.script[index.X, "command"] <- str_c("open ", file.name.whole.receptor.ligand, ".cif", " id ", index.M)
    index.X <- 1 + index.X
    ChimeraX.script[index.X, "command"] <- str_c("addh #", index.M, " & ligand hbond false")
    index.X <- 1 + index.X
    ChimeraX.script[index.X, "command"] <- str_c("save ", file.name.whole.receptor.ligand, ".pdb", " models #", index.M)
    index.X <- 1 + index.X
    ChimeraX.script[index.X, "command"] <- str_c("select #", index.M, " & ~ligand")
    index.X <- 1 + index.X
    ChimeraX.script[index.X, "command"] <- str_c("save ", file.name.whole.receptor, ".pdb", " models #", index.M, " select TRUE")
    index.X <- 1 + index.X
    ChimeraX.script[index.X, "command"] <- str_c("select #", index.M, " & ligand")
    index.X <- 1 + index.X
    ChimeraX.script[index.X, "command"] <- str_c("save ", file.name.ligand, ".pdb", " models #", index.M, " select TRUE")
    index.X <- 1 + index.X
    index.M <- 1 + index.M
  } #if the flexible residue file exists, merge with the rigid receptor and alphafold model
} # for select winner.xanthophyll receptor merge flexible, rigid, and missing receptor structures and save as cif without and with the ligand
write.table(ChimeraX.script, str_c(gromacs.dir, "/", receptor, ".cif2pdb.winner.xanthophyll.cxc"), row.names = FALSE, quote = FALSE, col.names = FALSE)
#
# plot measurements of interest binding energy score only
myplot <- ggplot(all.winner.dietary.xanthophyll[all.winner.dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" &
                                                  all.winner.dietary.xanthophyll$ligand.type %in% c("STQ2", "VFQ2", "8UH2", "BCR2") &
                                                  all.winner.dietary.xanthophyll$buried %in% c("epsilon.O23", "beta.O23", "beta.C23"), ], aes(x=ligand.type, y=measure.value)) +
  coord_cartesian(ylim = c(-15, -11), expand = TRUE, default = TRUE) +
  labs(y="Binding Energy Score Only", x="Ligand Type", title=str_c(receptor, " docked with dietary xanthophyll")) +
  theme_linedraw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_boxplot(alpha = 0.5, aes(color = ligand.type, fill = ligand.type), show.legend = FALSE)
ggsave(str_c(receptor, "rigid.C23.buried", as.character(numberofwinners), "winner.ligands", "Binding.Energy.Score.Only.BoxPlot.png", sep = "."), height=3, width=5.5, dpi=300)
ggsave(str_c(receptor, "rigid.C23.buried", as.character(numberofwinners), "winner.ligands", "Binding.Energy.Score.Only.BoxPlot.pdf", sep = "."), height=3, width=5.5, dpi=300)
#
#
# CODE DEVELOPMENT BELOW
#
# build the ChimeraX script to open pdbqt files for the ligand, add hydrogens, and save as PDB
ChimeraX.script <- tibble("command" = character())
index.X <- as.integer(1)
ChimeraX.script[index.X, "command"] <- str_c("cd ", work.dir)
index.X <- 1 + index.X
ChimeraX.script[index.X, "command"] <- "close #201"
index.X <- 1 + index.X
#
#open winner ligand pdbqt files, build missing hydrogens and save as pdb
index.F <- as.integer(1)
while (index.F <= nrow(winner.xanthophyll)) {
  if (file.exists(winner.xanthophyll[[index.F, "file.name.ligand.pdbqt"]])) {
    ChimeraX.script[index.X, "command"] <- str_c("open ", winner.xanthophyll[[index.F, "file.name.ligand.pdbqt"]], " id 201")
    index.X <- 1 + index.X
    ChimeraX.script[index.X, "command"] <- "addh #201 hbond false"
    index.X <- 1 + index.X
    file.name.ligand.pdb <- str_c(winner.directory, "/split_5flex", receptor, "_", winner.xanthophyll[[index.F, "ligand.name"]], "_ligand_", winner.xanthophyll[[index.F, "ligand.type"]], winner.xanthophyll[[index.F, "rotation"]], "_mode_", winner.xanthophyll[[index.F, "binding.mode"]], ".pdb")
    ChimeraX.script[index.X, "command"] <- str_c("save ", file.name.ligand.pdb, " models #201")
    index.X <- 1 + index.X
    ChimeraX.script[index.X, "command"] <- "close #201"
    index.X <- 1 + index.X
  } #if the ligand pdbqt file exists, add hydrogens, and save as pdb
  index.F <- 1 + index.F
} # for each winner.xanthophyll open its ligand pdbqt file, add hydrogens, and save as pdb
write.table(ChimeraX.script, str_c(receptor, ligand, "pdbqt2pdb.winner.xanthophyll.cxc", sep = "."), row.names = FALSE, quote = FALSE, col.names = FALSE)


dietary.xanthophyll[, "winner"] <- TRUE
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
# remove ligands with rotation for C6-C7 qand C26-C27 bonds
loser.identity <- dietary.xanthophyll[dietary.xanthophyll$rotation == 5, ]
dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
binding.energy <- data.frame(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" &
                                                   dietary.xanthophyll$rotation == 2, "measure.value"])
top.10p <- quantile(binding.energy[, "measure.value"], c(0, 0.1, 0.9))[2]
loser.identity <- dietary.xanthophyll[dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" &
                                        dietary.xanthophyll$measure.value > top.10p, ]
dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" 
                         & dietary.xanthophyll$ligand.type == "STQ2", ])
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" 
                         & dietary.xanthophyll$ligand.type == "VFQ2", ])
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" 
                         & dietary.xanthophyll$ligand.type == "8UH2", ])
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" 
                         & dietary.xanthophyll$ligand.type == "STQ2"
                         & dietary.xanthophyll$buried == "epsilon.O23", ])
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" 
                         & dietary.xanthophyll$ligand.type == "VFQ2"
                         & dietary.xanthophyll$buried == "beta.O23", ])
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" 
                         & dietary.xanthophyll$ligand.type == "8UH2"
                         & dietary.xanthophyll$buried == "beta.O23", ])
dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" 
                    & dietary.xanthophyll$ligand.type == "STQ2"
                    & dietary.xanthophyll$buried == "epsilon.O23", c("ligand.name", "binding.mode", "measure.value")]
dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" 
                    & dietary.xanthophyll$ligand.type == "VFQ2"
                    & dietary.xanthophyll$buried == "beta.O23", c("ligand.name", "binding.mode", "measure.value")]
dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" 
                    & dietary.xanthophyll$ligand.type == "VFQ2"
                    & dietary.xanthophyll$buried == "beta.O3", c("ligand.name", "binding.mode", "measure.value")]
dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" 
                    & dietary.xanthophyll$ligand.type == "8UH2"
                    & dietary.xanthophyll$buried == "beta.O23", c("ligand.name", "binding.mode", "measure.value")]
dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" 
                    & dietary.xanthophyll$ligand.type == "8UH2"
                    & dietary.xanthophyll$buried == "beta.O3", c("ligand.name", "binding.mode", "measure.value")]
#
# Box Plot for binding energy score only by ligand type for winner top 10%
myplot <- ggplot(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type %in% c("Binding.Energy.Score.Only"), ], aes(x=ligand.type, y=measure.value)) +
  coord_cartesian(ylim = c(-20, 2.5), expand = TRUE, default = TRUE) +
  labs(y="Binding Energy Score Only", x="Ligand Type", title=str_c(receptor, " docked with dietary xanthophyll")) +
  theme_linedraw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_boxplot(alpha = 0.5, aes(color = ligand.type, fill = ligand.type))
ggsave(str_c(receptor, "by.ligands.top10p", "Binding.Energy.ScoreOnly.BoxPlot.pdf", sep = "."), height=6, width=5.5, dpi=300)

#
#
# check binding energy differences
dietary.xanthophyll[, "winner"] <- TRUE
loser.identity <- dietary.xanthophyll[dietary.xanthophyll$measure.type == "Binding.Energy.Difference" & abs(dietary.xanthophyll$measure.value) < 5.0, "dock.id"]
dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
weird.energy.difference <- dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ]
#
# Exclude weird torsion angles
dietary.xanthophyll[, "winner"] <- TRUE
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
loser.identity <- dietary.xanthophyll[dietary.xanthophyll$measure.type == "torsion.c6c7" & dietary.xanthophyll$rotation == 5 & (dietary.xanthophyll$measure.value > -40.0 | dietary.xanthophyll$measure.value < -160), "dock.id"]
dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
loser.identity <- dietary.xanthophyll[dietary.xanthophyll$measure.type == "torsion.c26c27" & dietary.xanthophyll$rotation == 5 & (dietary.xanthophyll$measure.value > -40.0 | dietary.xanthophyll$measure.value < -160), "dock.id"]
dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
#
# Build the winner table and sort by binding energy as determined by score only
dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", c("ligand.name.index", "binding.mode", "buried", "measure.type", "measure.value")]
winner.xanthophyll <- dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ]
by.binding.energy <- order(winner.xanthophyll[["measure.value"]])
winner.xanthophyll <- winner.xanthophyll[by.binding.energy, ]
winner.xanthophyll[, c("ligand.name.index", "binding.mode", "buried", "measure.type", "measure.value")]


# Make the plots. May want ligand specific plots
binding.energy.tibble <- tibble(dietary.xanthophyll[dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", "measure.value"], 
                                dietary.xanthophyll[dietary.xanthophyll$measure.type == "Binding.Energy.Vina", "measure.value"], 
                                dietary.xanthophyll[dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", "rotation"], 
                                dietary.xanthophyll[dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", "binding.mode"],
                                dietary.xanthophyll[dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", "buried"],
                                dietary.xanthophyll[dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", "ligand.type"],
                                .name_repair = "unique")
y <- binding.energy.tibble[binding.energy.tibble$rotation == rotation.mode,]$measure.value...1
x <- binding.energy.tibble[binding.energy.tibble$rotation == rotation.mode,]$measure.value...2
binding.energy.rotation.lm <- lm(as.formula(y ~ 0 + x), binding.energy.tibble)
summary(binding.energy.rotation.lm)
#
# make the plot comparing binding energy measured by score only and by vina output
myplot <- ggplot(binding.energy.tibble[binding.energy.tibble$ligand.type == "STQ5", ], aes(x=measure.value...2, y=measure.value...1)) +
  coord_cartesian(xlim = c(-15.0, -0.0), ylim = c(-15, -0), expand = TRUE, default = TRUE) +
  labs(y="Binding Energy Score Only", x="Binding Energy Vina Output", title=str_c(receptor, " docked with dietary xanthophyll")) +
  theme_linedraw(base_size = 16) +
  #geom_line(data=binding.energy.rotation.lm, alpha = 0.5, aes(x=binding.energy.rotation.lm$model[,1], y=binding.energy.rotation.lm$model[,2]), color = "black") +
  geom_point(alpha = 0.20, aes(color = ligand.type))
ggsave(str_c(receptor, "all.ligands", "Binding.Energy.Linear.pdf", sep = "."), height=6, width=5.5, dpi=300)
#
# Box Plot for all measurements except torsion angles
myplot <- ggplot(dietary.xanthophyll[dietary.xanthophyll$measure.type %in% c("Binding.Energy.Score.Only", "Binding.Energy.Vina"), ], aes(x=measure.type, y=measure.value)) +
  coord_cartesian(ylim = c(-20, 2.5), expand = TRUE, default = TRUE) +
  labs(y="Value", x="Measurement", title=str_c(receptor, " docked with dietary xanthophyll")) +
  theme_linedraw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_boxplot(alpha = 0.5, aes(color = measure.type, fill = measure.type))
ggsave(str_c(receptor, "all.ligands", "Binding.Energy.BoxPlot.pdf", sep = "."), height=6, width=5.5, dpi=300)
#
# Box Plot for binding energy score only by ligand type
myplot <- ggplot(dietary.xanthophyll[dietary.xanthophyll$measure.type %in% c("Binding.Energy.Score.Only"), ], aes(x=ligand.type, y=measure.value)) +
  coord_cartesian(ylim = c(-20, 2.5), expand = TRUE, default = TRUE) +
  labs(y="Binding Energy Score Only", x="Ligand Type", title=str_c(receptor, " docked with dietary xanthophyll")) +
  theme_linedraw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_boxplot(alpha = 0.5, aes(color = ligand.type, fill = ligand.type))
ggsave(str_c(receptor, "by.ligands", "Binding.Energy.ScoreOnly.BoxPlot.pdf", sep = "."), height=6, width=5.5, dpi=300)
#
# Box Plot for binding energy vina by ligand type
myplot <- ggplot(dietary.xanthophyll[dietary.xanthophyll$measure.type %in% c("Binding.Energy.Vina"), ], aes(x=ligand.type, y=measure.value)) +
  coord_cartesian(ylim = c(-20, 2.5), expand = TRUE, default = TRUE) +
  labs(y="Binding Energy Vina", x="Ligand Type", title=str_c(receptor, " docked with dietary xanthophyll")) +
  theme_linedraw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_boxplot(alpha = 0.5, aes(color = ligand.type, fill = ligand.type))
ggsave(str_c(receptor, "by.ligands", "Binding.Energy.Vina.BoxPlot.pdf", sep = "."), height=6, width=5.5, dpi=300)
#
# Box Plot for torsion angle
myplot <- ggplot(dietary.xanthophyll[dietary.xanthophyll$measure.type %in% c("torsion.c6c7", "torsion.c26c27"), ], aes(x=measure.type, y=measure.value)) +
  coord_cartesian(ylim = c(-180, 180), expand = TRUE, default = TRUE) +
  labs(y="Value", x="Torsion Angle", title=str_c(receptor, " docked with dietary xanthophyll")) +
  theme_linedraw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_boxplot(alpha = 0.5, aes(color = measure.type, fill = measure.type))
ggsave(str_c(receptor, "all.ligands", "Torsion.Angle.BoxPlot.pdf", sep = "."), height=6, width=5.5, dpi=300)
#
# Box Plot for torsion angle c6-c7 by ligands
myplot <- ggplot(dietary.xanthophyll[dietary.xanthophyll$measure.type %in% c("torsion.c6c7"), ], aes(x=ligand.type, y=measure.value)) +
  coord_cartesian(ylim = c(-180, 180), expand = TRUE, default = TRUE) +
  labs(y="Torsion Angle C6-C7", x="Ligand Type", title=str_c(receptor, " docked with dietary xanthophyll")) +
  theme_linedraw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_boxplot(alpha = 0.5, aes(color = ligand.type, fill = ligand.type))
ggsave(str_c(receptor, "by.ligands", "Torsion.Angle.C6C7.BoxPlot.pdf", sep = "."), height=6, width=5.5, dpi=300)
#
# Box Plot for torsion angle c26-c27 by ligands
myplot <- ggplot(dietary.xanthophyll[dietary.xanthophyll$measure.type %in% c("torsion.c26c27"), ], aes(x=ligand.type, y=measure.value)) +
  coord_cartesian(ylim = c(-180, 180), expand = TRUE, default = TRUE) +
  labs(y="Torsion Angle C26-C27", x="Ligand Type", title=str_c(receptor, " docked with dietary xanthophyll")) +
  theme_linedraw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 )) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  geom_boxplot(alpha = 0.5, aes(color = ligand.type, fill = ligand.type))
ggsave(str_c(receptor, "by.ligands", "Torsion.Angle.C26C27.BoxPlot.pdf", sep = "."), height=6, width=5.5, dpi=300)
#
# Exclude Unfavroable Binding Energy > cutoff
dietary.xanthophyll[, "winner"] <- TRUE
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
cutoff.value <- as.numeric(-13.0)
loser.identity <- dietary.xanthophyll[dietary.xanthophyll$ligand.type == "stq" & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" & dietary.xanthophyll$measure.value > cutoff.value, "dock.id"]
dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
#
# Exclude if beta ring is buried
dietary.xanthophyll[, "winner"] <- TRUE
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
loser.identity <- dietary.xanthophyll[dietary.xanthophyll$buried == "beta.O3", "dock.id"]
dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
#
# Exclude if epsilon ring is buried
# loser.identity <- dietary.xanthophyll[dietary.xanthophyll$buried == "epsilon.O23", "dock.id"]
# dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
# nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
#
# Exclude weird torsion angles
dietary.xanthophyll[, "winner"] <- TRUE
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
loser.identity <- dietary.xanthophyll[dietary.xanthophyll$measure.type == "torsion.c6c7" & dietary.xanthophyll$rotation == 5 & (dietary.xanthophyll$measure.value > -40.0 | dietary.xanthophyll$measure.value < -160), "dock.id"]
dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
loser.identity <- dietary.xanthophyll[dietary.xanthophyll$measure.type == "torsion.c26c27" & dietary.xanthophyll$rotation == 5 & (dietary.xanthophyll$measure.value > -40.0 | dietary.xanthophyll$measure.value < -160), "dock.id"]
dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
#
# count docking outcomes for chi-squared tests
dietary.xanthophyll[, "winner"] <- TRUE
#
# only natural torsion angles allowed both ionone rings
loser.identity <- dietary.xanthophyll[dietary.xanthophyll$measure.type == "torsion.c6c7" & dietary.xanthophyll$rotation == 5 & (dietary.xanthophyll$measure.value > -40.0 | dietary.xanthophyll$measure.value < -160), "dock.id"]
dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
loser.identity <- dietary.xanthophyll[dietary.xanthophyll$measure.type == "torsion.c26c27" & dietary.xanthophyll$rotation == 5 & (dietary.xanthophyll$measure.value > -40.0 | dietary.xanthophyll$measure.value < -160), "dock.id"]
dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
#
# for STQ count which end is buried
loser.identity <- dietary.xanthophyll[dietary.xanthophyll$ligand.type != "stq" , "dock.id"]
dietary.xanthophyll <- flag.loser(dietary.xanthophyll, loser.identity)
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$buried == "epsilon.O23" & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$buried == "beta.O3" & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", ])
binding.energy <- data.frame(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only", "measure.value"])
top.10p <- quantile(binding.energy[, "measure.value"], c(0, 0.1, 0.9))[2]
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" & dietary.xanthophyll$measure.value < top.10p, ])
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" & dietary.xanthophyll$buried == "epsilon.O23" & dietary.xanthophyll$measure.value < top.10p, ])
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" & dietary.xanthophyll$buried == "beta.O3" & dietary.xanthophyll$measure.value < top.10p, ])
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" & dietary.xanthophyll$measure.value >= top.10p, ])
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" & dietary.xanthophyll$buried == "epsilon.O23" & dietary.xanthophyll$measure.value >= top.10p, ])
nrow(dietary.xanthophyll[dietary.xanthophyll$winner & dietary.xanthophyll$measure.type == "Binding.Energy.Score.Only" & dietary.xanthophyll$buried == "beta.O3" & dietary.xanthophyll$measure.value >= top.10p, ])
