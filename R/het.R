library(ggplot2)
library(plyr)

"%&%" <- function(a,b) paste(a, b, sep="")

extract_tstv <- function(x) {
  return( x$COUNT[7] / x$COUNT[8] )
}

parse_file <- function(x) {
  dat <- read.table(x, header=T)
  keepsamples <- read.table('/nas40t0/vasya/autism/resubmission/all_samples/original_het_samples.txt', header=T)
  subdat <- dat[dat$IID %in% keepsamples$INDV,]
}

args <- commandArgs(trailing=TRUE)

atlas.base <- args[1]
gatk.base <- args[2]
freebayes.base <- args[3]
mpileup.base <- args[4]
cges.base <- args[5]
pdf.file <- args[6]

callers <- list( 'Atlas', 'GATK', 'Freebayes', 'Mpileup', 'CGES')
files <- list(atlas.base %&% '.het',
              gatk.base %&% '.het',
              freebayes.base %&% '.het',
              mpileup.base %&% '.het',
              cges.base %&% '.het')


## read in tables
dat <- lapply(files, parse_file)

## drop unnecessary columns
dat <- Map( function(dataf) dataf[ ,c("IID", "F") ], dat )

## uniquify errors column names
dat <- Map( function(dataf, newname) rename(dataf, c("F"=newname)), dat, callers )

## Merge them into a single data frame properly ordered by sample
dat <- Reduce( function(...) merge(..., by.x=c('IID'), by.y=c('IID'), suffixes), dat )

## friendly form for plotting
plt.dat <- stack( dat, select = unlist(callers))


plt <- ggplot(plt.dat, aes(x=values, fill = ind)) + facet_grid( ind ~ . ,scale= "free") +
      geom_histogram(show_guide = FALSE) +
      xlab('F statistics') +
      ylab('# of samples') +
      theme_bw()

pdf(pdf.file)
show(plt)
dev.off()

