library(ggplot2)

"%&%" <- function(a,b) paste(a, b, sep="")

parse_file <- function(x) {
  read.table(x, header=T)
}

args <- commandArgs(trailing=TRUE)

atlas.base <- args[1]
gatk.base <- args[2]
freebayes.base <- args[3]
cges.base <- args[4]
pdf.file <- args[5]

callers <- list( 'Atlas', 'GATK', 'Freebayes', 'CGES')
files <- list(atlas.base %&% '.frq',
              gatk.base %&% '.frq',
              freebayes.base %&% '.frq',
              cges.base %&% '.frq')

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

