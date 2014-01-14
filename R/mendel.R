library(ggplot2)
library(RColorBrewer)
library(gridExtra)

calc_mendel <- function(x) {
  return( with(x, length(which(N>0))/length(N)) )
}

args <- commandArgs(trailing=TRUE)

atlas.file <- args[1]
gatk.file <- args[2]
freebayes.file <- args[3]
cges.file <- args[4]

## look at mendelian inconsistencies per trio
atlas.fmen <- read.table(atlas.file, header=T)
gatk.fmen <- read.table(gatk.file, header=T)
freebayes.fmen <- read.table(freebayes.file, header=T)
consensus.fmen <- read.table(cges.file, header=T)

## look at mendelian inconsistencies per locus
atlas.lmen <- read.table(gsub('.fmendel', '.lmendel', atlas.file), header=T)
gatk.lmen <- read.table(gsub('.fmendel','.lmendel', gatk.file), header=T)
freebayes.lmen <- read.table(gsub('.fmendel', '.lmendel', freebayes.file), header=T)
consensus.lmen <- read.table(gsub('.fmendel', '.lmendel', cges.file), header=T)

## count the total number of variants in the set 
freebayes.mvar <- length(freebayes.lmen$N)
consensus.mvar <- length(consensus.lmen$N)
atlas.mvar <- length(atlas.lmen$N)
gatk.mvar <- length(gatk.lmen$N)



trio_mvar <- list(atlas.mvar*sum(as.numeric(atlas.fmen$CHLD)),
             gatk.mvar*sum(as.numeric(gatk.fmen$CHLD)),
             freebayes.mvar*sum(as.numeric(freebayes.fmen$CHLD)),
             consensus.mvar*sum(as.numeric(consensus.fmen$CHLD)))

locus_mvar <- list(atlas.mvar, gatk.mvar, freebayes.mvar, consensus.mvar)

## format numerical data
atlas.lmen$N <- as.numeric(atlas.lmen$N)
gatk.lmen$N <- as.numeric(gatk.lmen$N)
freebayes.lmen$N <- as.numeric(freebayes.lmen$N)
consensus.lmen$N <- as.numeric(consensus.lmen$N)

atlas.fmen$N <- as.numeric(atlas.fmen$N)
gatk.fmen$N <- as.numeric(gatk.fmen$N)
freebayes.fmen$N <- as.numeric(freebayes.fmen$N)
consensus.fmen$N <- as.numeric(consensus.fmen$N)


## throw data into one large dataframe
callers <- list( 'Atlas', 'GATK', 'Freebayes', 'CGES')
fmendel.dat <- list(atlas.fmen$N, gatk.fmen$N, freebayes.fmen$N, consensus.fmen$N)
lmendel.dat <- list(which(atlas.lmen$N>0), which(gatk.lmen$N > 0), which(freebayes.lmen$N > 0), which(consensus.lmen$N > 0))
mendel <- data.frame( trio_error_rate = (unlist(lapply(fmendel.dat, sum)) / unlist(trio_mvar)) * 100,
                      locus_error_rate = unlist(lapply(lmendel.dat, length)) / unlist(locus_mvar) * 100,
                      Callers = unlist(callers)) 

mendel$locus_order_callers <- reorder(mendel$Callers, mendel$locus_error_rate)
mendel$trio_order_callers <- reorder(mendel$Callers, mendel$trio_error_rate)


## plot dat graph
locus_plt <- ggplot( mendel, aes(x=locus_order_callers, y=locus_error_rate, fill=Callers) ) + 
          geom_bar(stat="identity", show_guide = FALSE) +
          theme_bw() +
          xlab("Caller") +
          ylab("% of Sites with >1 Mendelian Errors") 

trio_plt <- ggplot( mendel, aes(x=locus_order_callers, y=trio_error_rate, fill=Callers) ) + 
          geom_bar(stat="identity", show_guide = FALSE) +
          theme_bw() +
          xlab("Caller") +
          ylab("% of Proband Genotypes with Mendelian errors") 

x11()
show(locus_plt)
x11()
show(trio_plt)
