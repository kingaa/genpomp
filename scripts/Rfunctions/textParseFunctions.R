## Function to parse genetic trees or infection trees all saved to one file
## returns a list of dataframes, each dataframe representing a tree
getTrees <- function(treeFile) {
  d <- readLines(treeFile)
  firstElement <- function(x) { unlist(strsplit(x,' '))[1] }
  f <- sapply(d,firstElement)
  headers <- as.numeric(which(f == 'node'))
  df <- split(d,cut(seq_along(d),breaks=c(headers-1,length(d)),labels=F))
  df <- lapply(df,function(x) read.table(text = x, header = T))
}


