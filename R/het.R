library(ggplot2)
library(plyr)

"%&%" <- function(a,b) paste(a, b, sep="")

extract_tstv <- function(x) {
  return( x$COUNT[7] / x$COUNT[8] )
}

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
files <- list(atlas.base %&% '.het',
              gatk.base %&% '.het',
              freebayes.base %&% '.het',
              cges.base %&% '.het')


## read in tables
dat <- lapply(files, parse_file)

## drop unnecessary columns
dat <- Map( function(dataf) dataf[ ,c("INDV", "F") ], dat )

## uniquify errors column names
dat <- Map( function(dataf, newname) rename(dataf, c("F"=newname)), dat, callers )

## Merge them into a single data frame properly ordered by sample
dat <- Reduce( function(...) merge(..., by.x=c('INDV'), by.y=c('INDV'), suffixes), dat )

## friendly form for plotting
plt.dat <- stack( dat, select = unlist(callers))


plt <- ggplot(plt.dat, aes(y=values, x=c(rep(1:length(dat$INDV), length(callers))), color = ind)) + facet_grid( ind ~ . ,scale= "free") +
      geom_point(show_guide = FALSE) +
      xlab('F statistics') +
      ylab('# of samples') +
      theme_bw()

pdf(pdf.file)
show(plt)
dev.off()

