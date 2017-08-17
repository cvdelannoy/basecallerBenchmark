#!/usr/bin/Rscript
# Load arguments:
args = commandArgs(TRUE)
int = args[1]

.libPaths( c("rPackages",.libPaths()) )
trapOut = suppressMessages(library(data.table))
trapOut = suppressMessages(library(ggplot2))
source("rScripts/ggMultiplot.R")

# List all files ---------------------------------------------------------
bcFiles = list.files(paste0(int,"/basecallerMetrics"), full.names=T)
asFiles = list.files(paste0(int,"/assemblerMetrics"), full.names=T)
# Rolling mean graphs ----------------------------------------------------

rmFiles = grep("bcRollingAvg_", bcFiles, value=T)
rmdtList = list(); rmdtIndex = 1

for(rmf in rmFiles){
  # Calculate number of bases in pileup file
  n = system(paste0("wc -l ", rmf, collapse = ''), intern=T)
  n = as.numeric(stringr::str_match(n, pattern="^[0-9]+"))
  
  # If number of bases exceeds 100000, plotting may take unfeasibly long.
  # For now, skip if that is the case
  if (n>100000){
    print("Too many bases in pileup, plotting skipped")
    next
    # nRead = floor(n/100000)
    # 
    # # Load rows, average, create data point, repeat
    # dt = data.table(coverage = rep(NA_real_, 10000), identity = NA_real_)
    # con <- file(rmf)
    # open(con)
    # for (rl in 1:10000){
    #   curSet = strsplit(readLines(con = con, n = nRead), "\t")
    #   dtLine = list(coverage= mean(as.numeric(sapply(curSet, "[[", 2))), 
    #                 identity=mean(as.numeric(sapply(curSet, "[[", 3))))
    #   set(dt, rl, names(dt), dtLine)
  } else {
    dt = fread(rmf)
    setnames(dt, names(dt), c("position", "coverage", "identity"))
  }
  
  svg(filename = paste0(int,"/images/", basename(rmf),".svg"),
    width = 15, height=5, pointsize = 12)
  ggId = ggplot(data = dt, aes(x=position, y=identity)) + geom_line()
  ggCov = ggplot(data = dt, aes(x=position, y=coverage)) + geom_line()
  multiplot(ggId, ggCov, cols = 1)
  dev.off()
}

# Repeat graphs ----------------------------------------------------------
repFiles = grep("repeatRegionsError_", bcFiles, value=T)

# Create data table to contain repeat error data of all data sets
# (assuming same # of repeats have been assessed for all)
n = system(paste0("wc -l ", repFiles[1], collapse = ''), intern=T)
n = as.numeric(stringr::str_match(n, pattern="^[0-9]+"))-1
repDt = data.table(repeats = rep(NA_integer_, n*length(repFiles)), 
                   repeatType = NA_character_, 
                   matchRatio = NA_real_, 
                   dataset = NA_character_)
ri = 1
for (rf in repFiles){
  curRf = read.table(rf)
  set(repDt, ri:(ri+n-1), c("repeats", "repeatType", "matchRatio"), as.list(curRf))
  dsName = stringr::str_match( basename(rf),"(?<=repeatRegionsError_).+(?=.bam.pileup.txt)")
  set(repDt, ri:(ri+n-1), "dataset", dsName)
  ri = ri+n
}

# plot
svg(filename= paste0(int,"/images/deletionRatio.svg"),
    width = 10, height = 5, pointsize = 12)
ggplot(data = subset(repDt, !is.na(matchRatio)), aes(x=repeats)) + 
  geom_line(aes(y=matchRatio, color = repeatType, linetype = dataset))
dev.off()

# K-mer graphs -----------------------------------------------------------
kmFiles = grep("kmerCount_.+\\.txt", asFiles, value=T, perl=T)

if (length(kmFiles) != 0){
  for(kf in kmFiles){
    nameKf = as.character(stringr::str_match(kf,"(?<=kmerCount_).+(?=.txt)"))
    curKf = fread(kf, col.names = c("kmer", nameKf))
    setkey(curKf, kmer)
    set(curKf, ,2L, curKf[,2L]/sum(curKf[,2L]))
    if(!exists("kmDt")){kmDt = curKf}else{kmDt = kmDt[curKf]}
  }
  
  # Mark homoNT and diNT k-mers
  homoNt = grepl(pattern = "^(.)\\1*$", kmDt$kmer, perl=T)
  diNt = grepl(pattern = "(.)(?!\\1)(.)(\\1\\2)+(\\1)", kmDt$kmer, perl=T)
  kmDt[,"kmerType" := "other"]
  kmDt[homoNt,"kmerType" := "homopolymer repeat"]
  kmDt[diNt,"kmerType":="dimeric repeat"]
  
  # Plot
  kmNames = stringr::str_match(kmFiles,"(?<=kmerCount_).+(?=.txt)")
  graphCount = 1
  for(km1 in 1:length(kmNames)){
    for(km2 in (km1+1):length(kmNames)){
      kmgg = ggplot(kmDt, aes_string(x=kmNames[km1], y=kmNames[km2])) + 
        geom_point(aes(color=kmerType)) +
        geom_abline(intercept =0, slope = 1, color = "grey") +
        geom_text(data = subset(kmDt, kmerType!="other"), 
                  nudge_y = 0.0002, size = 3, check_overlap=F, 
                  aes_string(x = kmNames[km1], y=kmNames[km2], label="kmer")) +
        labs(color = "k-mer type", x=paste0("occurrence ", kmNames[km1]), y=paste0("occurrence ", kmNames[km2])) + 
        scale_color_manual(values=c("dimeric repeat" = "#F8766D","homopolymer repeat" =  "#00BFC4","other" = "#00BA38"))
  
      # save plot
      filename = paste0("int/images/kmerCountPlot_",kmNames[km1],"_",kmNames[km2],".png")
      png(filename= filename,
          width = 30, height = 15, units = "cm", res = 600, pointsize = 10)
      kmgg + theme_grey(base_size = 14) 
      dev.off()
      
      graphCount = graphCount + 1
    }
  }
}