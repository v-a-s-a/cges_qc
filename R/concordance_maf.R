library(ggplot2)

"%&%" <- function(a,b) paste(a, b, sep="")

parse_file <- function(x) {
  read.table(x, header=T)
}

args <- commandArgs(trailing=TRUE)

consensus.base <- args[1]
high.base <- args[2]
mid.base <- args[3]
low.base <- args[4]
pdf.file <- args[5]

callers <- list( 'Consensus', '3of4', '2of4', 'Union')
files <- list(consensus.base %&% '.frq',
              high.base %&% '.frq',
              mid.base %&% '.frq',
              low.base %&% '.frq')

## read in tables
dat <- lapply(files, parse_file)

## drop unnecessary columns
dat <- Map( function(dataf) dataf[ ,"MAF" ], dat )

## and ye shall be named
names(dat) <- callers

plt.dat <- stack(dat, callers)

plt <- ggplot(plt.dat, aes(x=values, fill=ind)) + facet_grid( ind ~ . ,scale= "free") +
      geom_histogram(show_guide = FALSE) +
      xlab('MAF') +
      ylab('Counts')

plt2 <- ggplot(plt.dat, aes(values, colour=ind)) + stat_ecdf() +
      labs(x='MAF', y='Proportion of Variants with MAF <= x') +
      theme_bw()

pdf(pdf.file)
show(plt)
show(plt2)
dev.off()

