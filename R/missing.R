library(ggplot2)
library(plyr)

"%&%" <- function(a,b) paste(a, b, sep="")

parse_imen <- function(x) {
  read.table(x, header=T)
}

args <- commandArgs(trailing=TRUE)

atlas.base <- args[1]
gatk.base <- args[2]
freebayes.base <- args[3]
cges.base <- args[4]
pdf.file <- args[5]

sets <- list('ATLAS', 'GATK', 'Freebayes', 'CGES')
files <- list(atlas.base %&% '.imiss',
              gatk.base %&% '.imiss',
              freebayes.base %&% '.imiss',
              cges.base %&% '.imiss')

## read in tables
dat <- lapply(files, parse_imen)

## drop unnecessary columns
dat <- Map( function(x) x[,!(names(x) %in% c('MISS_PHENO','N_MISS','N_GENO'))], dat )

## uniquify errors column names
dat <- Map( function(dataf, newname) rename(dataf, c("F_MISS"=newname)), dat, sets )

## Merge them into a single data frame properly ordered by sample
dat <- Reduce( function(...) merge(..., by.x=c('FID','IID'), by.y=c('FID', 'IID'), suffixes), dat )

## friendly form for plotting
plt.dat <- stack( dat, select = unlist(sets) )


plt <- ggplot(plt.dat, aes(x=values, fill = ind)) +
  facet_grid( ind ~ . ,scale= "free") +
  geom_histogram(binwidth=0.003) +
  xlab('# of Samples') + ylab('Rate of Missing Genotypes')


pdf(pdf.file)
show(plt)
dev.off()

