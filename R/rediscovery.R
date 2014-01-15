library(ggplot2)
library(gridExtra)

"%&%" <- function(a,b) paste(a, b, sep="")

#' Make a data frame of TsTV summary data
make_kg_fpath <- function(base) {
  fpath <- paste( base, '.kg', sep = '' )
  return(fpath)
}

make_evs_fpath <- function(base) {
  fpath <- paste( base, '.evs', sep = '' )
  return(fpath)
}


args <- commandArgs(trailing=TRUE)



atlas.base <- args[1]
gatk.base <- args[2]
freebayes.base <- args[3]
cges.base <- args[4]
pdf.file <- args[5]


set.names <- list( 'Atlas', 'GATK', 'Freebayes', 'CGES')
call.sets = list( atlas.base
                  gatk.base,
                  freebayes.base,
                  cges.base)


evs.fpaths <-  mapply(make_evs_fpath, call.sets)
evs.dat <- lapply( evs.fpaths, read.table )
evs <- data.frame( evs = unlist( lapply(evs.dat, subset, select=V2) ),
                  Callers = unlist(set.names))
evs$order_sets <- reorder(evs$Callers, evs$evs)


evs.plt <- ggplot( evs, aes(x=order_sets, y=evs, fill=Callers) ) +
        geom_bar(stat="identity", show_guide = FALSE) +
        labs(x="Caller", y="% of Variants Found in Exome Variant Server") +
        theme_bw() +
        colScale
kg.fpaths <-  mapply(make_kg_fpath, call.sets)
kg.dat <- lapply( kg.fpaths, read.table )

kg <- data.frame( kg = unlist( lapply(kg.dat, subset, select=V2) ),
                   Callers = unlist(set.names) )
kg$order_sets <- reorder(kg$Callers, kg$kg)


kg.plt <- ggplot( kg, aes(x=order_sets, y=kg, fill=Callers) ) +
        geom_bar(stat="identity", show_guide = FALSE) +
        labs(x="Caller", y="% of Variants Found in 1000 Genomes") +
        theme_bw()


pdf(pdf.file)
show(kg.plt)
show(evs.plt)
dev.off()
