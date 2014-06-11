library(ggplot2)
library(gridExtra)

calc_mendel <- function(x) {
  return( with(x, length(which(N>0))/length(N)) )
}

"%&%" <- function(a,b) paste(a, b, sep="")


args <- commandArgs(trailing=TRUE)

consensus.base <- args[1]
high.base <- args[2]
mid.base <- args[3]
union.base <- args[4]
pdf.file <- args[5]

## look at mendelian inconsistencies per trio
consensus.fmen <- read.table(consensus.base %&% ".fmendel", header=T)
high.fmen <- read.table(high.base %&% ".fmendel", header=T)
mid.fmen <- read.table(mid.base %&% ".fmendel", header=T)
union.fmen <- read.table(union.base %&% ".fmendel", header=T)

## look at mendelian inconsistencies per locus
consensus.lmen <- read.table(consensus.base %&% ".lmendel", header=T)
high.lmen <- read.table(high.base %&% ".lmendel", header=T)
mid.lmen <- read.table(mid.base %&% ".lmendel", header=T)
union.lmen <- read.table(union.base %&% ".lmendel", header=T)

## count the total number of variants in the set 
mid.mvar <- length(mid.lmen$N)
consensus.mvar <- length(consensus.lmen$N)
union.mvar <- length(union.lmen$N)
high.mvar <- length(high.lmen$N)



trio_mvar <- list(consensus.mvar*sum(as.numeric(consensus.fmen$CHLD)),
             high.mvar*sum(as.numeric(high.fmen$CHLD)),
             mid.mvar*sum(as.numeric(mid.fmen$CHLD)),
             union.mvar*sum(as.numeric(union.fmen$CHLD)))

locus_mvar <- list(consensus.mvar, high.mvar, mid.mvar, union.mvar)

## format numerical data
consensus.lmen$N <- as.numeric(consensus.lmen$N)
high.lmen$N <- as.numeric(high.lmen$N)
mid.lmen$N <- as.numeric(mid.lmen$N)
union.lmen$N <- as.numeric(union.lmen$N)

consensus.fmen$N <- as.numeric(consensus.fmen$N)
high.fmen$N <- as.numeric(high.fmen$N)
mid.fmen$N <- as.numeric(mid.fmen$N)
union.fmen$N <- as.numeric(union.fmen$N)


## throw data into one large dataframe
callers <- list( 'Consensus', '3of4', '2of4', 'Union')
fmendel.dat <- list(consensus.fmen$N,
                    high.fmen$N,
                    mid.fmen$N,
                    union.fmen$N)
lmendel.dat <- list(which(consensus.lmen$N>0),
                    which(high.lmen$N > 0),
                    which(mid.lmen$N > 0),
                    which(union.lmen$N > 0))
mendel <- data.frame( trio_error_rate = (unlist(lapply(fmendel.dat, sum)) / unlist(trio_mvar)) * 100,
                      locus_error_rate = (unlist(lapply(lmendel.dat, length)) / unlist(locus_mvar)) * 100,
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

pdf(pdf.file)
show(locus_plt)
show(trio_plt)
dev.off()
