#!/usr/bin/env Rscript
library(data.table)
library(dplyr)

# Collect arguments from the command line
args <- commandArgs(trailingOnly = TRUE)
intersectOut <- args[1]
fstOut <- args[2]
mafs.g1 <- args[3]
mafs.g2 <- args[4]
graph.me.fst <- args[5]

# Define headers
fst.headers <- c("A", "AB", "f", "FST", "Pvar")
intersect.headers <- c("Chr","bp")

# Read in data from .fst file and shared.pos file
intersect <- fread(intersectOut, sep = "\t", col.names = intersect.headers)
fst <- fread(fstOut, sep = "\t", col.names = fst.headers %>%
               select(c("Pvar")))
mafs1 <- fread(mafs.g1 %>%
                 select(c(chromo, position, major, minor, knownEM)))
mafs2 <- fread(mafs.g2 %>%
                 select(c(chromo, position, major, minor, knownEM)))

# Combine two datasets
FST <- cbind(intersect, fst, mafs1, mafs2)

# Save to .fst file
write.table(FST, graph.me.fst, col.names = F, row.names = F)