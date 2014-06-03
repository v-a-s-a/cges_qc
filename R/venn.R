## make some variant count venn diagrams

library(VennDiagram)
library(RColorBrewer)

"%&%" <- function(a,b) paste(a, b, sep="")

## turn a map file into a vector of variant IDs
varvec <- function(mapFile) {
  dat <- read.table(mapFile)
  return(dat$V1 %&% ':' %&% dat$V4)
}

args <- commandArgs(trailing=TRUE)

atlas.base <- args[1]
gatk.base <- args[2]
freebayes.base <- args[3]
mpileup.base <- args[4]
cges.base <- args[5]
pdf.file <- args[6]

## store variants in vectors

sets <- list('Atlas', 'GATK', 'Freebayes', 'Mpileup', 'CGES')
bases <- list(atlas.base, gatk.base, freebayes.base, mpileup.base, cges.base)
mapFiles <- lapply(bases, function(x) x %&% '.map')
var <- lapply(mapFiles, varvec)

print(lapply(var, head))

venn.plt <- draw.quad.venn(area1 = length(var[[1]]),
               area2 = length(var[[2]]),
               area3 = length(var[[3]]),
               area4 = length(var[[4]]),
               n12 = length(intersect(var[[1]], var[[2]])),
               n13 = length(intersect(var[[1]], var[[3]])),
               n14 = length(intersect(var[[1]], var[[4]])),
               n23 = length(intersect(var[[2]], var[[3]])),
               n24 = length(intersect(var[[2]], var[[4]])),
               n34 = length(intersect(var[[3]], var[[4]])),
               n123 = length(Reduce(intersect, var[c(1,2,3)])),
               n124 = length(Reduce(intersect, var[c(1,2,4)])),
               n134 = length(Reduce(intersect, var[c(1,3,4)])),
               n234 = length(Reduce(intersect, var[c(2,3,4)])),
               n1234 = length(var[[5]]),
               category = c('ATLAS', 'GATK', 'Freebayes', 'Mpileup'),
               fill = c( "#F8766D", "#C77CFF", "#00BFC4", "#FF0000"),
               scaled = TRUE)
pdf(pdf.file)
grid.draw(venn.plt)
dev.off()


