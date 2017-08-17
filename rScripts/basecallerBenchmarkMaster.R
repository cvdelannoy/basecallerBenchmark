#'Minion basecaller benchmark V3

#' Load arguments:
#' 1. Pileup file
#' 2. Intermediate files path
#' 3. number of threads
args = commandArgs(TRUE)
pileupName <- args[1]
intPath <- args[2]
nThreads <- as.integer(args[3])

# arg check and prerequisites ------------------------------
.libPaths( c(.libPaths(), "rPackages") )
# Install packages if necessary
list.of.packages = c("ggplot2", "reshape2","data.table", "doMC", "foreach")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, lib="rPackages", repos="http://cran.xl-mirror.nl/")

library(ggplot2)

# Load functions
source("rScripts/basecallerBenchmarkFunctions.R")

# Overall identity, error rates ---------------------------------------

overallMetrics = overallStatistics(pileupName)
listNames = names(overallMetrics)
overallMetrics = data.frame(matrix(unlist(overallMetrics), nrow=4))
row.names(overallMetrics) = listNames
write.table(overallMetrics, file=paste0(intPath,"/basecallerMetrics/bcOverallMetrics_", 
                                        basename(pileupName),".txt"), col.names = F, row.names= T, quote = F)


# Overall rolling average identity ------------------------------------
raFilename = paste0(intPath,"/basecallerMetrics/bcRollingAvg_",basename(pileupName),".txt")
rollingAverageIdentity(pileupName, 10000, raFilename)

# Accuracy repeat regions ------------------------------
r = 3:15
u = c("A","C","G","T")
outDfHp = delRateInRepeatsParallel(r,pileupName,u, nThreads)

u2 = apply(combn(c("A","C","G","T"),2),2,paste, collapse='')
outDfDp = delRateInRepeatsParallel(r,pileupName,u2, nThreads)

outDf = base::merge(outDfHp[,c(1,3)], outDfDp[,c(1,3)], 
                    by="repeats", suffixes = c("_homopolymeric", "_dimeric"))
outDf = reshape2::melt(outDf, id.vars = "repeats")
names(outDf)[2:3] = c("RepeatType", "matchRatio")
outDf[,"RepeatType"] = sapply(as.character(outDf$RepeatType), function(x) strsplit(x, "matchRatio_")[[1]][2])
write.table(outDf, file = paste0(intPath,"/basecallerMetrics/repeatRegionsError_",basename(pileupName),".txt"))


# ggRepeatRatio = ggplot(data = outDf, aes(x=repeats)) + geom_point(aes(y=value, color = variable))

# ggRepeatRatio = ggplot(data = outDf, aes(x=repeats)) + geom_point(aes(y=deletionFraction, color = RepeatType))
# svg(filename= paste0(intPath,"/images/deletionRatio_", basename(pileupName),".svg"),
#     width = 10, height = 5, pointsize = 12)
# ggRepeatRatio
# dev.off()


