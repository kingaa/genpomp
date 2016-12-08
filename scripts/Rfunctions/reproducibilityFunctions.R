## Functions to allow for any experiment to be easily reproduced

## Function to read the script being run into text that can be saved
getScript <- function(){
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  script <- readLines(paste0(getwd(),"/",script.name))
  return(script)
}

## Function to get current status of the git repository
getGitCommit <- function(){
  fileName <- paste0("lastCommit",runif(1),".txt")
  system(paste0("git cat-file commit HEAD > ",fileName))
  commit <- readLines(fileName)
  system(paste0("rm ", fileName))
  return(commit)
}

## Function to write parameters to a file that can be read by standalone executables
writeParams <- function(params,filename){
  df <- data.frame(params)
  write.table(df,file=filename,quote=FALSE,col.names=FALSE)
}

## Function to copy files for building the executable and make the executable
## Assumes that the working directory is 'standalone', sister to 'src'
## usermodelDirectory is the name of the directory that houses usermodel.h
makeExecutable <- function(mainfileStem, mainfiledir, srcdir, tmpdir, usermodelDirectory){
  includefiles <- c(
    "substmodel.h",
    "lmatrix.h",
    "basenode.h",
    "gnode.h",
    "node.h",
    "tree.h",
    "particleFilter.h",
    "mif.h",
    "abstractUsermodel.h",
    "abstractPerturb.h",
    "userUnif.h",
    "RngStream.h",
    "io.h",
    "reduceSeqs.h",
    "type_defs.h",
    "gsl_randist.h",
    "gsl_rng.h",
    "constants.h",
    "tmpdirMakefile"
  )

  cfiles <- c(
    "basenode.cc",
    "gnode.cc",
    "node.cc",
    "substmodel.cc",
    "io.cc",
    "reduceSeqs.cc",
    "userUnif.c",
    "RngStream.c",
    "peel.c",
    "constants.c",
    "beta.c",
    "exponential.c",
    "gamma.c",
    "gausszig.c",
    "gsl_unif.c"
  )

  for (file in c(includefiles, cfiles))
    file.copy(from = file.path(srcdir, file),to = tmpdir, overwrite = TRUE)
  file.copy(from = file.path(mainfiledir,paste0(mainfileStem, '.cc')), to = tmpdir, overwrite = TRUE)
  file.copy(from = file.path(usermodelDirectory, '/usermodel.h'), to = tmpdir, overwrite = TRUE)
  file.copy(from = file.path(usermodelDirectory, '/perturb.h'), to = tmpdir, overwrite = TRUE)
  wd <- setwd(tmpdir)
  on.exit(setwd(wd))
  file.rename("tmpdirMakefile","Makefile")
  system2("make",paste0(mainfileStem, '.exe'))
}

## Reads in a sequence file, strips the sequences, writes the file without
## sequences to tmpdir
## This allows for fitting models in genpomp to only diagnoses, for example
strip_sequences <- function(seq_file, tmpdir){
  data <- read.table(seq_file, header = F)
  names(data) <- c('time','sequence')
  data$sequence <- NA
  write.table(data, file = file.path(tmpdir,"no_seqs.txt"), row.names = F, col.names = F, quote = F)
  file.path(tmpdir,"no_seqs.txt")
}



