
##80############################################################################
## This is a script for running n particle filters over simulated data #########
## The usage is displayed when not passing the correct number of parameters ####
##80############################################################################

require(digest)

"This script runs iterated filtering
Usage: pfilter.R --data_file=<dat> --usermodel_dir=<usermodeldir> --params_file=<par> --result_dir=<dest> --n_reps=<nreps> --n_threads=<nthreads> --n_nested=<nnested> --n_branch_samples=<nbranchsamples> --n_particles=<nparticles> --environment=<environment> [--remove_invariant_sites --save_internals --reproducible] [PARAMS]...
" -> doc

library(docopt)
opt <- docopt::docopt(doc)

## Read in the parameters file
params_file <- opt$params_file
params_df <- read.table(file = params_file) 
params <- params_df$V2
names(params) <- params_df$V1  

## Load functions for reproducibility and calculating initial conditions
source('../Rfunctions/reproducibility_functions.R')
source('../Rfunctions/text_parse_functions.R')
source('../Rfunctions/initial_conditions.R')

## Build executable
tmp_dir <- tempdir()
usermodel_dir <- file.path("../../../src/models", opt$usermodel_dir)
if(opt$reproducible) {
  Sys.setenv(CMD_LINE_FLAGS='-DREPRODUCIBLE')
} else {
  Sys.unsetenv('CMD_LINE_FLAGS')  
}
makeExecutable('pfilter', mainfiledir = "../../standalone", tmpdir = tmp_dir, srcdir = "../../src",
               usermodelDirectory = usermodel_dir)

## Extract info on last git commit and extract this script to an object
pf_runscript <- getScript()
pf_gitcommit <- getGitCommit()

## Set algorithmic parameters
params['np'] <- as.numeric(opt$n_particles)
params['num_threads'] <- as.numeric(opt$n_threads)
params['num_nested'] <- as.numeric(opt$n_nested)
params['num_branch_samples'] <- as.numeric(opt$n_branch_samples)

## Set user specified parameters 
p <- opt$PARAMS
if(length(p) > 0){
  for(i in 1:length(p)){
    par <- strsplit(p[i], split = '=')
    params[par[[1]][1]] <- as.numeric(par[[1]][2])
  }
}

## Create a unique name for the results directory if requested
result_dir <- opt$result_dir
if(result_dir == "generate"){  
  result_dir <- paste(unlist(strsplit(date(), split = ' ')), collapse = '_')
} 
res_stem <- '/scratch/kingaa_flux/alxsmth/results/pfilter_results'
result_dir <- file.path(res_stem, result_dir)
if(!file.exists(result_dir)){
  print(paste0("Creating directory: ", result_dir))
  dir.create(result_dir)
}

## Get a unique id based on the parameter vector and some random deviates
unique_id <- digest(c(params, rnorm(10)))

## Run the particle filter over the data n_reps times
for(i in 1:opt$n_reps){
  ## Name a unique result directory
  rep_resdir <- file.path(result_dir, paste0(unique_id,'_replicate_', i, '/'))  
  ## Create a directory to write results to
  dir.create(rep_resdir)
  print(paste0("Writing results to: ", rep_resdir))
  
  ## Write parameters used for filtering to the results directory
  par_filename <- file.path(rep_resdir, 'params.txt')
  writeParams(params, par_filename)
  
  print(paste0(c(par_filename, opt$data_file, "tamuraNei", rep_resdir, opt$remove_invariant_sites, opt$save_internals), collapse = ' '))
  # Run the standalone executable
  start <- Sys.time()
  system2(file.path(tmp_dir, "pfilter.exe"),
          c(par_filename, opt$data_file, "tamuraNei", rep_resdir, opt$remove_invariant_sites, opt$save_internals))
  end <- Sys.time()
  
  ## Compute runtime
  runtime <- end - start
  ## Read in likelihoods to save
  condlogliks <- read.table( paste0(rep_resdir,"logliks.txt"), header = T)
  
  ## Obtain a unique id for this pfilter run
  pf_unique_id <- digest(condlogliks)
  
  ## Write unique ids to results file
  write(pf_unique_id, file = file.path(rep_resdir,'pf_unique_id.txt'))
  
  ## Save results, parameters, git commit info and this script to an Rdata file
  master_filename <- file.path(rep_resdir, "pfilter.Rdata")
  save(list = c('runtime', 'params', 'pf_runscript', 'pf_unique_id',
                'pf_gitcommit', 'environment', 'condlogliks'), file = master_filename)
}
