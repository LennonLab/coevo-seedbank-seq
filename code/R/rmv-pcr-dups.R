# Removing PCR duplicates using FastUniq v1.1
# https://sourceforge.net/projects/fastuniq/files/
# installed in my home directory at 
# `/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source`. 
fastuniq.path <- "/N/u/danschw/Carbonate/my_tools/FastUniq/FastUniq/source/fastuniq"

# phage or host?
  input.arg <- NULL
  
  input.arg <- commandArgs(trailingOnly = TRUE)
  # test if there is an argument: if not, return an error
  if ( length(input.arg)!=1 ) 
    stop("one argument must be supplied: phage OR host", call.=FALSE)
  
  if (! input.arg %in% c("phage", "host"))
    stop("one argument must be supplied: phage OR host", call.=FALSE)


# load libraries
  library(here, quietly = T, verbose = F)
  library(tidyverse, quietly = T, verbose = F)

######################################
# read and parse names of fastq files
######################################
  if (input.arg =="phage"){
    reads <-
      tibble(path = list.files(here("data/input/fastq", input.arg), 
                               pattern = "fastq", full.names = TRUE) ,
             f = list.files(here("data/input/fastq", input.arg), 
                            pattern = "fastq", full.names = FALSE) )%>%
      separate(1,into=c("run","num","trt", "line", "transfer", "bc", "read", "sufx"), sep=regex("-|_"), remove = F) %>% 
      mutate(unq.sample = interaction(trt, line, transfer, sep = '-'))
  }
  if (input.arg =="host"){
    reads <-
      tibble(path = list.files(here("data/input/fastq", input.arg), 
                            pattern = "fastq", full.names = TRUE) ,
             f = list.files(here("data/input/fastq", input.arg), 
                               pattern = "fastq", full.names = FALSE) )%>%
      separate(f,into=c("run","trt", "line", "transfer", "extract", "bc", "read", "sufx"),
               sep=regex("-|_"), remove = F) %>% 
      mutate(unq.sample = interaction(trt, line, transfer,extract, sep = '-'))
  }  


########################################
# make a folder for deduplicated data
########################################
#make directories
if (!dir.exists(here("data", "ddup-fastq")))
  dir.create(here("data", "ddup-fastq"))
if (!dir.exists(here("data", "ddup-fastq", input.arg)))
  dir.create(here("data", "ddup-fastq", input.arg))


############################
# write commands for batch job
############################
  
  # folder to save the batch scripts commands 
  if (!dir.exists(here("code/bash", "ddup-scripts")))
    dir.create(here("code/bash", "ddup-scripts"))
  if (!dir.exists(here("code/bash", "ddup-scripts", input.arg)))
    dir.create(here("code/bash", "ddup-scripts", input.arg))
  
for (current.reads in unique(reads$unq.sample)){
  
  # reads in current loop iteration
  paths <- reads %>% 
    filter(unq.sample == current.reads) %>% 
    select(path, f) %>% 
    mutate(f = str_remove(f, ".gz$")) %>% 
    mutate(out = paste0(here("data/ddup-fastq/"), input.arg, "/ddup-",f)) %>% 
    #fix extraction names to unify across sequencing runs
    mutate(out = str_replace(out, "SB", "rS")) %>% 
    mutate(out = str_replace(out, "veg", "rV"))
             

  # path for current batch script
    path.sh <- here("code/bash/ddup-scripts", input.arg,
                    paste0("deduplicate-", current.reads, ".sh"))
    
  #slrum header
    write_lines(c("#!/bin/bash
#SBATCH --mail-user=danschw@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL,BEGIN,END",
    paste0("#SBATCH --job-name=ddup-", current.reads)), 
                    path.sh)
    
    #copy and extract paired read data
    write_lines("\n#copy and extract paired read data",
                file = path.sh, append = TRUE)
    write_lines(paste("zcat", paths$path, ">", 
                      here("data", "ddup-fastq", input.arg, paths$f )),
                file = path.sh, append = TRUE)

    # make input list for FastUniq
    # Each list file should contain the file names of paired sequencing files
    input.list <- here("data/ddup-fastq", input.arg, paste0("input_list_",current.reads,".txt"))
    write_lines("\n# make input list for FastUniq",
                file = path.sh, append = TRUE)
    write_lines(paste("touch", input.list), 
                file = path.sh, append = TRUE)
    write_lines(paste("echo",paths$f, ">>", input.list),
                file = path.sh, append = TRUE)
    
    # write_lines(paths$f, file = input.list,append = FALSE)
    
    #FastUniq commands
    write_lines("\n#de-duplicate with FastUniq (local)",
                path.sh, append = TRUE)
    
    write_lines(paste("cd ", here("data", "ddup-fastq", input.arg)),
                path.sh, append = TRUE)

    sys.cmd <- paste(
      fastuniq.path,
      "-i", input.list, # input list
      "-t q", # Output sequence format: fastq
      "-o", paths$out[1], # first output file
      "-p", paths$out[2], # second output file
      "-c 1") # New serial numbers assigned by FastUniq
    
    write_lines(paste("\n", sys.cmd), path.sh, append = T)
    
    # clean up
    write_lines(paste("\n", "# delete copied inputs"),
                path.sh, append = T)
    write_lines(paste("rm ", here("data", "ddup-fastq", input.arg, paths$f)),
                      path.sh, append = T)
    write_lines(paste("\n", "# delete input list"),
                path.sh, append = T)
    write_lines(paste("rm ", input.list),
                path.sh, append = T)
    write_lines(paste("\n", "# compress outputs"),
                path.sh, append = T)
    write_lines(paste("gzip ", paths$out),
                path.sh, append = T)
    
     
    # submit job
    system(paste("sbatch", path.sh))
}
  
  
quit(save = "no")
