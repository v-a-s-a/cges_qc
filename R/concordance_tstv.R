library(ggplot2)

"%&%" <- function(a,b) paste(a, b, sep="")

extract_tstv <- function(x) {
  return( x$COUNT[7] / x$COUNT[8] )
}

args <- commandArgs(trailing=TRUE)

consensus.base <- args[1]
high.base <- args[2]
mid.base <- args[3]
union.base <- args[4]
pdf.file <- args[5]

consensus.tstv <- read.table(consensus.base %&% ".TsTv.summary", header=T)
high.tstv <- read.table(high.base %&% ".TsTv.summary", header=T)
mid.tstv <- read.table(mid.base %&% ".TsTv.summary", header=T)
union.tstv <- read.table(union.base %&% ".TsTv.summary", header=T)

callers <- list('Consensus', '3of4', '2of4', 'Union')
tstv.dat <- list(consensus.tstv, high.tstv, mid.tstv, union.tstv)


tstv <- data.frame( tstv = unlist( lapply(tstv.dat, extract_tstv) ),
                    Callers = unlist(callers))
tstv$order_callers <- reorder(tstv$Callers, tstv$tstv)

plt <- ggplot( tstv, aes(x=order_callers, y=tstv, fill=Callers) ) +
        geom_bar(stat="identity", show_guide = FALSE) +
        labs(x="Callers", y="Ts/Tv") +
        theme_bw()
 
pdf(pdf.file)
show(plt)
dev.off()


