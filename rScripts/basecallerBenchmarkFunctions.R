list.of.packages = c("stringr", "doMC", "foreach")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos="http://cran.xl-mirror.nl/", lib="~/rPackages/")

# Calculate deletion rate in repeats of lengths i, from pileup file, of repeat units u
delRateInRepeatsStreaming = function(i, pileupName, u){
    # Construct completed regex
    patFwd = ""
    for(n in 1:length(u)){
      subPat = paste0("((?<!",u[n],")",u[n],"{",i,"}(?!",u[n], "))",collapse='')
      if(n != length(u)){subPat = paste0(subPat,"|", collapse = '')}
      patFwd = paste0(patFwd,subPat)
    }
    # Construct reverse pattern regex
    patRev = chartr("ATGC","TACG",patFwd)
    
    # Construct data table
    dt = data.table::data.table(base =rep(NA_character_,1000), read = NA_character_)
    matches = 0; dels = 0; nRepeats = 0
    
    con <- file(pileupName)
    open(con)
    # Load first 1000 bases into data table 
    cl = strsplit(readLines(con, n=1000), "\t", fixed = T); cl = lapply(cl, "[", c(3,5))
    for(cn in c(1L,2L)){data.table::set(dt,,cn,sapply(cl,"[", cn))}
    #Define how much overlap must exist between batches to avoid missing hits
    overlapSize = as.integer(i*nchar(u[1]) + 2 - 1)
    while (length(cl)!=0){
      # Find positions of repeat regions
      curString = paste0(unlist(dt[,"base"]), collapse = '')
      foundFwd = stringr::str_locate_all(curString, patFwd)[[1]]
      foundRev = stringr::str_locate_all(curString, patRev)[[1]]
      # If found, retrieve indices for positions and count pileup markings
      allInd=c()
      if(nrow(foundFwd)!=0){
        ffInd = c()
        for (n in 1:nrow(foundFwd)){ffInd = c(ffInd, seq(foundFwd[n,1],foundFwd[n,2]))}
        matches = matches + sum(stringr::str_count(dt[ffInd,"read"], pattern = "\\."))
        allInd = c(allInd, ffInd)}
      if(nrow(foundRev)!=0){
        frInd = c()
        for (n in 1:nrow(foundRev)){frInd = c(frInd, seq(foundRev[n,1],foundRev[n,2]))}
        matches = matches + sum(stringr::str_count(dt[frInd,"read"], pattern = ","))
        allInd = c(allInd, frInd)}
      if(length(allInd)){
        nRepeats = nRepeats + nrow(foundFwd) + nrow(foundRev);
        dels = dels + sum(as.numeric(stringr::str_match_all(dt[allInd,"read"], pattern = "(?<=-)[0-9]+")[[1]]))}
      # Load new bases into data table
      cl = strsplit(readLines(con, n=1000-overlapSize), "\t", fixed = T)
      cl = lapply(cl, "[", c(3,5))
      dt[,names(dt) := data.table::shift(dt, n=(1000-overlapSize), type="lead")]
      if(!(length(cl)==0)){
        for(cn in c(1L,2L)){data.table::set(dt,(overlapSize+1):(overlapSize+length(cl)),cn,sapply(cl,"[", cn))}
        if(!length(cl)==1000){
          dt=dt[1:(length(cl)+overlapSize),]
        }
      }
    }
    close(con)
    coverage = (matches+dels)/(i*nchar(u[1])*nRepeats)
    outList = list(repeats = i, coverage = coverage, matchRatio = matches/(matches+dels))
    return(outList)
}

# Wrapper for delRatesInRepeatsStreaming, runs function in parallel
delRateInRepeatsParallel = function(r, pileupName, u, nCores){
  require(foreach)
  #TODO: check, all elements of u same nchar
  doMC::registerDoMC(cores = nCores)
  outDf = foreach::foreach(i = r, .packages = c("stringr","data.table")) %dopar% delRateInRepeatsStreaming(i, pileupName, u)
  outDf = data.frame(matrix(unlist(outDf), nrow=(max(r) - min(r)+1), byrow=T))
  colnames(outDf) = c("repeats", "coverage", "matchRatio")
  return(outDf)
}

overallStatistics = function(pileupName){
  con <- file(pileupName)
  open(con)
  coverage = identity = deletions = insertions = 0
  
  cl = unlist(strsplit(readLines(con, n=1), "\t", fixed = T))[c(4,5)]
  while(length(cl)!=0){
    coverage = coverage + as.numeric(cl[1])
    identity = identity + sum(stringr::str_count(cl[2], pattern = "\\.|\\,"))
    deletions = deletions + sum(as.numeric(unlist((stringr::str_match_all(cl[2], pattern = "(?<=-)[0-9]+")))))
    insertions = insertions + sum(as.numeric(unlist((stringr::str_match_all(cl[2], pattern = "(?<=\\+)[0-9]+")))))
    cl = unlist(strsplit(readLines(con, n=1), "\t", fixed = T))[c(4,5)]
  }
  close(con)
  outList = list(coverage = coverage, identity = identity, deletions = deletions, insertions = insertions)
  return(outList)
}

# Works, but way too slow: processes ~50000 bases per minute
rollingAverageIdentity = function(pileupName, w, outFile){
  if(file.exists(outFile)){file.remove(outFile)}
  file.create(outFile)
  # Create Table to store rolling average scores
  n = system(paste0("wc -l ", pileupName, collapse = ''), intern=T)
  n = as.numeric(stringr::str_match(n, pattern="^[0-9]+"))
  
  # If window bigger than n, make smaller rolling window:1/100 of n
  if(w > n){w = floor(n/100)}
  # outDt = data.table::data.table(position = rep(NA_integer_,n), coverage = NA_real_, meanIdentity= NA_real_)
  
  # create datatable to temporarily store results, occasionally flushed to file
  outDt = data.table::data.table(position = rep(NA_integer_,1000), tab1 = "\t", coverage = NA_real_, tab2="\t", identity = NA_real_, newLine = "\n")
  outDtRow = 1L
  
  # create datatable to allow streaming assessment of window size
  dtSize = as.integer(w)
  dt = data.table::data.table(position=rep(NA_integer_,dtSize), coverage = NA_integer_, read = NA_character_)
  
  con <- file(pileupName)
  open(con)
  for(i in 1:n){
    cl = as.list(unlist(strsplit(readLines(con, n=1), "\t", fixed = T))[c(2,4,5)])
    if(length(cl)==0){
      print("Error in rolling mean calculation\n:
              missing position value. Pileup file may be corrupted or stale file handle.");
      break
    }
    cl[[1]] = as.integer(cl[[1]])
    cl[[2]] = as.integer(cl[[2]])
    dt[,names(dt) := data.table::shift(dt, n=1, type="lead")]
    data.table::set(dt, dtSize, names(dt), cl)
    if(is.na(dt[1,"position"])){next}
    identity = sum(stringr::str_count(unlist(dt[,"read"]), pattern = "\\.|\\,"))
    coverage = sum(dt[,"coverage"])
    outList = as.list(c(cl[[1]], (coverage/w), (identity/coverage)))
    data.table::set(outDt, outDtRow, c("position", "coverage","identity"), outList)
    # outVec = c(cl[[1]], (coverage/w), (identity/coverage))
    if (outDtRow==1000){
      dtString = paste0("\"",unlist(t(outDt)),"\"", collapse = "")
      system(paste0("echo ", dtString," >> ", outFile), intern=T)
      outDtRow=0L
    }
    outDtRow = outDtRow+1L
  }
  # flush last rows
  dtString = paste0("\"",unlist(t(outDt[1:(outDtRow-1),])),"\"", collapse = "")
  system(paste0("echo ", dtString," >> ", outFile), intern=T)
  close(con)
  # return(outDt)
}
