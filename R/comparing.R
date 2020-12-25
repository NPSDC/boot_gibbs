load("environment/seList.RData")
load("../../real_data/out_1/sim_counts_matrix.rda")

se100 <- seList[["100"]]
se100 <- computeSizeFactors(se100, counts_matrix[,1])
se100 <- computeConfInt(se100, sf = T)
covInds <- computeCoverage(counts_matrix[,1], se100, list(seq(nrow(counts_matrix))), prop = F)[[1]]
print(sapply(covInds, length))
print(nrow(counts_matrix))

commInds <- Reduce(intersect, covInds)

### Boot only
bootOnly <- setdiff(setdiff(covInds[[1]],covInds[[2]]), covInds[[3]])
print(length(bootOnly))

## Gibbs only
gibbsOnly <- setdiff(setdiff(covInds[[2]],covInds[[1]]), covInds[[3]])
print(length(gibbsOnly))

## Polee Only
poleeOnly <- setdiff(setdiff(covInds[[3]],covInds[[1]]), covInds[[2]])
print(length(poleeOnly))

library(VennDiagram)
g <- venn.diagram(covInds, category.names=names(covInds), filename=NULL, output=T)
grid.draw(g)